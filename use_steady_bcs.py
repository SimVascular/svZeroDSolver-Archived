import numpy as np
import run_0d_solver_fast

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
    locations = ["inlet", "outlet"]
    num_variables_for_BCs = ({  "FLOW"          : 1, # [Q]
                                "PRESSURE"      : 1, # [P]
                                "RESISTANCE"    : 2, # [R, distal_pressure]
                                "RCR"           : 3  # [Rp, C, Rd]
                            })
    capacitor_indices = ({  "RCR"       : [3],    # DATATABLE format: [t, Rp, t, C, t, Rd]
                            "CORONARY"  : [5, 7]  # [t, Ra1, t, Ra2, t, Ca, t, Cc, t, Rv1, t, Pv_distal_pressure, P_im]
                        })
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
                        new_datatable_values = new_datatable_values + [time[0], time_averaged_value, time[-1], time_averaged_value]
                    parameters["datatable_values"][datatable_name] = new_datatable_values
                    if bc_type in capacitor_indices:
                        for ind in capacitor_indices[bc_type]:
                            parameters["datatable_values"][datatable_name][ind] = 0 # change capacitance to zero
                    # print("datatable_name = ", datatable_name)
                    # print('parameters["datatable_values"][datatable_name] = ', parameters["datatable_values"][datatable_name])
            elif bc_type == "CORONARY":
                # use mean intramyocardial pressure for the CORONARY BCs
                time_of_intramyocardial_pressure, bc_values_of_intramyocardial_pressure = run_0d_solver_fast.extract_bc_time_and_values(12, len(parameters["datatable_values"][datatable_name]), parameters, segment_number, location)
                if len(time_of_intramyocardial_pressure) >= 2:
                    cardiac_cycle_period = time_of_intramyocardial_pressure[-1] - time_of_intramyocardial_pressure[0]
                    time_averaged_value = compute_time_averaged_bc_value_for_single_cardiac_cycle(time_of_intramyocardial_pressure, bc_values_of_intramyocardial_pressure, cardiac_cycle_period)
                    parameters["datatable_values"][datatable_name] = parameters["datatable_values"][datatable_name][:12] + [time_of_intramyocardial_pressure[0], time_averaged_value, time_of_intramyocardial_pressure[-1], time_averaged_value]
                    for ind in capacitor_indices[bc_type]:
                        parameters["datatable_values"][datatable_name][ind] = 0 # change capacitance to zero
                    # print("datatable_name = ", datatable_name)
                    # print('parameters["datatable_values"][datatable_name] = ', parameters["datatable_values"][datatable_name])

    return parameters
