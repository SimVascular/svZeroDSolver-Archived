import copy
import numpy as np
import run_0d_solver_fast

class altered_bc_block:
    def __init__(self, name, type, segment_number, location, parameter_list):
        self.block_name = name
        self.type = type
        self.segment_number = segment_number
        self.parameter_list = parameter_list
        self.location = location

def compute_time_averaged_bc_value_for_single_cardiac_cycle(time, bc_values, cardiac_cycle_period):
    """
    Inputs:
        np.array time
            = defined over a single cardiac cycle (time[-1] = time[0] + cardiac_cycle_period)
        np.array bc_values
            = defined over a single cardiac cycle (bc_values[0] == bc_values[-1])
        float cardiac_cycle_period
    """
    time_averaged_value = (1.0/cardiac_cycle_period)*np.trapz(bc_values, time)
    return time_averaged_value

def use_steady_state_values_for_bcs(parameters):
    """
    Convert pulsatile BCs into equivalent steady BC by
        1. using mean value as the value of the steady BC
        2. removing capacitors from BCs (convert RCR and CORONARY into equivalent resistance-only BCs)

    Caveats:
        - removing the capacitors from the RCR and Coronary changes the 0d simulation results outputted by run_0d_solver_fast.run_network_util() because the internal variables originally present in the RCR and Coronary BCs are removed. These internal variables must be restored (and set to their steady state values) to use the steady state solutions as the initial conditions for the pulsatile simulations.
    """
    locations = ["inlet", "outlet"]
    num_variables_for_BCs = ({  "FLOW"          : 1, # [Q]
                                "PRESSURE"      : 1, # [P]
                                "RESISTANCE"    : 2, # [R, distal_pressure]
                                "RCR"           : 3  # [Rp, C, Rd]
                            })
    altered_bc_blocks = [] # a list of the modified BC blocks (modified meaning "converted from RCR or coronary to resistance"); these BC blocks lost their internal variable in the modification
    for location in locations:
        for segment_number in parameters["boundary_condition_types"][location]:
            datatable_name = parameters["boundary_condition_datatable_names"][location][segment_number]
            bc_type = parameters["boundary_condition_types"][location][segment_number]
            if bc_type in num_variables_for_BCs:
                num_variables = num_variables_for_BCs[bc_type]
                total_num_values = len(parameters["datatable_values"][datatable_name])
                if total_num_values == 0:
                    message = "Error. DATATABLE for segment #" + str(segment_number) + " is missing values."
                    raise RuntimeError(message)
                num_values_per_variable = int(total_num_values/num_variables)
                if num_values_per_variable >= 4: # num_values_per_variable == 2 means one data pt prescribed for BC, so BC is already steady; not possible for num_values_per_variable to be odd-valued
                    new_datatable_values = []
                    for k in range(num_variables):
                        start_index = k*num_values_per_variable
                        end_index = start_index + num_values_per_variable
                        time, bc_values = run_0d_solver_fast.extract_bc_time_and_values(start_index, end_index, parameters, segment_number, location)
                        if len(time) < 2:
                                message = "Error. len(time) < 2"
                                raise RuntimeError(message)
                        cardiac_cycle_period = time[-1] - time[0]
                        time_averaged_value = compute_time_averaged_bc_value_for_single_cardiac_cycle(time, bc_values, cardiac_cycle_period)
                        new_datatable_values = new_datatable_values + [time[0], time_averaged_value]
                    parameters["datatable_values"][datatable_name] = new_datatable_values
                if bc_type == "RCR":
                    t = parameters["datatable_values"][datatable_name][0]
                    Rp = parameters["datatable_values"][datatable_name][1]
                    Rd = parameters["datatable_values"][datatable_name][5]
                    equivalent_R = Rp + Rd # add resistances in series
                    parameters["datatable_values"][datatable_name] = [t, equivalent_R, t, 0] # setting distal pressure to zero b/c default RCR BC has a distal pressure of zero
                    parameters["boundary_condition_types"][location][segment_number] = "RESISTANCE"
                    altered_bc_blocks.append(altered_bc_block("BC" + str(segment_number) + "_" + location, bc_type, segment_number, location, [Rp]))
            elif bc_type == "CORONARY":
                # use mean intramyocardial pressure for the CORONARY BCs
                time_of_intramyocardial_pressure, bc_values_of_intramyocardial_pressure = run_0d_solver_fast.extract_bc_time_and_values(12, len(parameters["datatable_values"][datatable_name]), parameters, segment_number, location)
                if len(time_of_intramyocardial_pressure) >= 2:
                    cardiac_cycle_period = time_of_intramyocardial_pressure[-1] - time_of_intramyocardial_pressure[0]
                    time_averaged_value = compute_time_averaged_bc_value_for_single_cardiac_cycle(time_of_intramyocardial_pressure, bc_values_of_intramyocardial_pressure, cardiac_cycle_period)
                    parameters["datatable_values"][datatable_name] = parameters["datatable_values"][datatable_name][:12] + [time_of_intramyocardial_pressure[0], time_averaged_value, time_of_intramyocardial_pressure[-1], time_averaged_value]

                t = parameters["datatable_values"][datatable_name][12]
                Pv_distal_pressure = parameters["datatable_values"][datatable_name][11]
                Ra1 = parameters["datatable_values"][datatable_name][1]
                Ra2 = parameters["datatable_values"][datatable_name][3]
                Cim = parameters["datatable_values"][datatable_name][7]
                Rv1 =  parameters["datatable_values"][datatable_name][9]
                Pim =  parameters["datatable_values"][datatable_name][13]
                equivalent_R = Ra1 + Ra2 + Rv1
                parameters["datatable_values"][datatable_name] = [t, equivalent_R, t, Pv_distal_pressure]
                parameters["boundary_condition_types"][location][segment_number] = "RESISTANCE"
                altered_bc_blocks.append(altered_bc_block("BC" + str(segment_number) + "_" + location, bc_type, segment_number, location, [Ra1, Ra2, Cim, Pim]))
    return parameters, altered_bc_blocks

def restore_internal_variables_for_capacitance_based_bcs(y_f, ydot_f, var_name_list_f, altered_bc_blocks):
    """
    Restore the internal variables (and set them to their steady state values).
    """
    y0 = copy.deepcopy(y_f)
    ydot0 = copy.deepcopy(ydot_f)
    var_name_list = copy.deepcopy(var_name_list_f)
    for block in altered_bc_blocks:
        if block.type in ["RCR", "CORONARY"]:
            if block.location == "outlet":
                if block.type == "RCR":
                    # get BC parameters
                    Rp = block.parameter_list[0]
                    # get soltns
                    Pin = y0[var_name_list.index("P_V" + str(block.segment_number) + "_" + block.block_name)]
                    Qin = y0[var_name_list.index("Q_V" + str(block.segment_number) + "_" + block.block_name)]

                    P_internal = Pin - Rp * Qin
                    var_name_list.append("var_0_" + block.block_name)
                    y0 = np.append(y0, np.array(P_internal))
                if block.type == "CORONARY":
                    # get BC parameters
                    Ra1 = block.parameter_list[0]
                    Ra2 = block.parameter_list[1]
                    Cim = block.parameter_list[2]
                    Pim = block.parameter_list[3]
                    # get soltns
                    Pin = y0[var_name_list.index("P_V" + str(block.segment_number) + "_" + block.block_name)]
                    Qin = y0[var_name_list.index("Q_V" + str(block.segment_number) + "_" + block.block_name)]

                    Pd = Pin - Qin * (Ra1 + Ra2)
                    volume_internal = Cim * (Pd - Pim)
                    var_name_list.append("var_0_" + block.block_name)
                    y0 = np.append(y0, np.array(volume_internal))
                ydot0 = np.append(ydot0, np.zeros(1)) # the time derivative of the soltn is zero b/c at steady-state
            else:
                message = "Error. This function does not work for 'location' = " + block.location
                raise RuntimeError(message)
    return y0, ydot0, var_name_list
