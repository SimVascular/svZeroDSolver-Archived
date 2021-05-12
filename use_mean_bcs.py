"""
All possible BCs:

Inlet BCs:
1. steady flow - 1 time pt prescribed
2. steady flow - 2 time pts prescribed
3. pulsatile flow - 3 or more time pts prescribed
4. steady pressure - 1 time pt prescribed
5. steady pressure - 2 time pts prescribed
6. pulsatile pressure - 3 or more time pts prescribed

Outlet BCs:
1. steady flow - 1 time pt prescribed
2. steady flow - 2 time pts prescribed
3. pulsatile flow - 3 or more time pts prescribed
4. steady pressure - 1 time pt prescribed
5. steady pressure - 2 time pts prescribed
6. pulsatile pressure - 3 or more time pts prescribed
7. steady resistance with steady distal pressure - 1 time pt prescribed for each
8. steady resistance with steady distal pressure - 2 time pts prescribed for each
9. pulsatile resistance with steady distal pressure - 3 or more time pts prescribed for each
10. steady rcr - 1 time pt prescribed for each variable (Rp, C, Rd)
11. steady rcr - 2 time pts prescribed for each variable (Rp, C, Rd)
12. pulsatile rcr - 3 or more time pts prescribed for each variable (Rp, C, Rd)
13. steady coronary - 1 time pt prescribed for intramyocadial pressure only
14. steady coronary - 2 time pts prescribed for intramyocadial pressure only
15. pulsatile coronary - 3 or more time pts prescribed for intramyocadial pressure only

todo: delete this entire comment section on the BCs
"""

def compute_time_averaged_bc_value_for_single_cardiac_cycle(time, bc_values, cardiac_cycle_period):
    time_averaged_value = (1.0/cardiac_cycle_period)*np.trapz(bc_values, time)
    return time_averaged_value

def use_mean_values_for_bcs(parameters):
    locations = ["inlet", "outlet"]
    num_variables_for_BCs = ({  "FLOW"          : 1, # [Q]
                                "PRESSURE"      : 1, # [P]
                                "RESISTANCE"    : 2, # [R, distal_pressure]
                                "RCR"           : 3  # [Rp, C, Rd]
                            })
    for location in locations:
        bc_datatable_types_segment_numbers = list(parameters["boundary_condition_types"][location].keys())
        for segment_number in bc_datatable_types_segment_numbers:
            datatable_name = parameters["boundary_condition_datatable_names"][location][segment_number]
            if parameters["boundary_condition_types"][location][segment_number] in num_variables_for_BCs:
                num_variables = num_variables_for_BCs[parameters["boundary_condition_types"][location][segment_number]]
                total_num_values = len(parameters["datatable_values"][datatable_name])
                if total_num_values == 0:
                    message = "Error. DATATABLE for segment #" + str(segment_number) + " is missing values."
                    raise RuntimeError(message)
                num_values_per_variable = int(total_num_values/num_variables)
                if num_values_per_variable >= 2: # num_values_per_variable == 1 means steady BC already
                    new_datatable_values = []
                    for k in range(num_variables):
                        start_index = k*num_values_per_variable
                        end_index = start_index + num_values_per_variable
                        time, bc_values = extract_bc_time_and_values(start_index, end_index, parameters, segment_number, location)
                        if len(time) < 2:
                                message = "Error. len(time) < 2"
                                raise RuntimeError(message)
                        cardiac_cycle_period = time[-1] - time[0]
                        time_averaged_value = compute_time_averaged_bc_value_for_single_cardiac_cycle(time, bc_values, cardiac_cycle_period)
                        new_datatable_values = new_datatable_values + [time[0], time_averaged_value, time[-1], time_averaged_value]
                    parameters["datatable_values"][datatable_name] = new_datatable_values
            elif parameters["boundary_condition_types"][location][segment_number] == "CORONARY":
                # use mean intramyocardial pressure for the CORONARY BCs
                time_of_intramyocardial_pressure, bc_values_of_intramyocardial_pressure = extract_bc_time_and_values(start_index = 12, end_index = len(parameters["datatable_values"][datatable_name]), parameters = parameters, segment_number = segment_number, location = location)
                if len(time) >= 2:
                    cardiac_cycle_period = time_of_intramyocardial_pressure[-1] - time_of_intramyocardial_pressure[0]
                    time_averaged_value = compute_time_averaged_bc_value_for_single_cardiac_cycle(time_of_intramyocardial_pressure, bc_values_of_intramyocardial_pressure, cardiac_cycle_period)
                    parameters["datatable_values"][datatable_name] = parameters["datatable_values"][datatable_name][:12] + [time_of_intramyocardial_pressure[0], time_averaged_value, time_of_intramyocardial_pressure[-1], time_averaged_value]
    return parameters

    last here - check that this code is correct
