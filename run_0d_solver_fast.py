"""
Author: Pham, Jonathan

This code simulates the 0D model described in the 0D solver input file by creating network_util_NR::LPNBlock objects for each 0D element and running the network_util_NR solver routines.

Available vessel modeling types:
    1) R (constant resistor)
    2) C (constant capacitor)
    3) L (constant inductor)
    4) RC (constant resistor-capacitor)
    5) RL (constant resistor-inductor)
    6) RCL (constant resistor-capacitor-inductor)
    7) custom, user-defined vessel model

Available junction modeling types:
    1) NORMAL_JUNCTION (mass conservation and pressure continuity only)
    2) custom, user-defined junction type

Available boundary conditions (BCs):
- inlet BCs:
        1) constant flow rate
        2) time-dependent flow rate
        3) constant pressure
        4) time-dependent pressure
        5) custom, user-defined inlet BC
- outlet BCs:
        1) constant resistor with user-prescribed constant distal pressure value
        2) time-dependent resistor with user-prescribed time-dependent distal pressure value
        3) constant RCR with a distal pressure of zero
        4) time-dependent RCR with a distal pressure of zero
        5) constant flow rate
        6) time-dependent flow rate
        7) constant pressure
        8) time-dependent pressure
        9) steady coronary with time-dependent intramyocadial pressure with user-prescribed constant distal pressure
        10) custom, user-defined outlet BC
"""

import sys
import os
import numpy as np
import scipy.interpolate
import matplotlib.pyplot as plt
from tqdm import tqdm
try:
    import oneD_to_zeroD_convertor
except ImportError:
    message = "Error. oneD_to_zeroD_convertor.py was not imported. This code is needed to extract relevant information from the 0D solver input file."
    raise ImportError(message)
try:
    import networkx as nx # only needed if you want to visualize the 0d model as a directed graph
except ImportError:
    print("\nnetworkx not found. networkx is needed only if you want to visualize your 0d model as a directed graph.")
try:
    import blocks as ntwku
    import connections
    import time_integration as time_int
except ImportError:
    message = "Error. 0D solver was not imported. This code is needed to create the network_util_NR::LPNBlock objects for the 0D elements and run the 0D simulations."
    raise ImportError(message)
try:
    from profilehooks import profile # only needed if you want to profile this script
except ImportError:
    print("\nprofilehooks.profile not found. profilehooks.profile is needed only if you want to profile this script.")
np.set_printoptions(threshold=sys.maxsize)
import importlib
import argparse

""" TODOs """
# currently here 8/26/20: should try runing code profiling on this python code when i have free time to see how to make this code more efficient, which will be useful for large 0d models, ie pulmonaries # https://stackoverflow.com/questions/582336/how-can-you-profile-a-python-script

def import_custom_0d_elements(custom_0d_elements_arguments_file_path):
    """
    Purpose:
        Import the user-defined custom 0d element file that specifies the arguments for the desired user-defined custom 0d elements.
    Inputs:
        string custom_0d_elements_arguments_file_path
            = path to user-defined custom 0d element file
    Returns:
        module custom_0d_elements_arguments
            = module to call custom 0d element arguments from
    """
    try:
        sys.path.append(os.path.dirname(custom_0d_elements_arguments_file_path))
        custom_0d_elements_arguments = importlib.import_module(os.path.basename(os.path.normpath(custom_0d_elements_arguments_file_path)))
        sys.path.pop()
        return custom_0d_elements_arguments
    except ImportError:
        message = "Error. Custom 0d elements desired, but 'custom_0d_elements_arguments_file_path' not provided or file not found."
        raise ImportError(message)
    except Exception:
        message = "Error. There is an issue with the custom 0d element file provided."
        raise Exception(message)

def create_custom_element(custom_element_class_name, custom_element_args):
    """
    Purpose:
        Create an instance of the custom 0d element class, custom_element_class_name
    Inputs:
        string custom_element_class_name
            = the name of the custom 0d element class
        dict custom_element_args
            = the arguments for the custom 0d element class's parameters
    Returns:
        custom_element_class_name instance_of_custom_class
            = an instance of the custom 0d element class, custom_element_class_name
    Reference: https://www.tutorialspoint.com/How-to-convert-a-string-to-a-Python-class-object
    """
    instance_of_custom_class = getattr(sys.modules[ntwku.__name__], custom_element_class_name)(**custom_element_args)
    return instance_of_custom_class

def create_unsteady_bc_value_function(time, bc_values):
    """
    Purpose:
        Create a time-dependent function for bc_values (boundary condition values) by interpolating bc_values over time with a cubic spline and periodic boundaries
    Caveats:
        1) bc_values must have the same start and end values, since periodic boundariess are used for the cubic spline
        2) bc_values and time must have the same length
    Inputs:
        list time
            = [list of time points]
        list bc_values
            = [list of bc_values]
    Returns:
        cubic spline for bc_values as a function of time
    Reference:
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.interpolate.CubicSpline.html
    """
    if len(bc_values) == 2 and bc_values[0] == bc_values[1]:
        return lambda x: bc_values[0]
    else:
        return scipy.interpolate.CubicSpline(np.array(time), np.array(bc_values), bc_type = 'periodic')

def create_bc_equations(parameters):
    """
    Purpose:
        Create equations for all prescribed, built-in boundary conditions as functions of time.
    Caveat:
        This assumes that if we are using a FLOW, PRESSURE, RESISTANCE, or RCR BC (for a DATATABLE of type LIST), then every parameter value for the BC must be specified over the same number of time points. For example, for a time-varying RCR BC, if Rp is prescribed for 5 different time points in the DATATABLE in the 0d solver input file, then C and Rd must be prescribed for 5 time points as well.
        Time-dependent boundary conditions must be prescribed over a single cardiac cycle only.
    Inputs:
        dict parameters
            -- created from function oneD_to_zeroD_convertor.extract_info_from_solver_input_file
    Returns:
        void but updates parameters to include:
            float cardiac_cycle_period
                = period of a cardiac cycle
            dict bc_equations
                =   {
                        "inlet" : {segment number : [list of equations for BC parameeter values as functions of time]}
                        "outlet" : {segment number : [list of equations for BC parameeter values as functions of time]}
                    }
    """
    bc_equations = {"inlet" : {}, "outlet" : {}}
    built_in_bc_types = {"inlet" : ["FLOW", "PRESSURE"], "outlet" : ["RESISTANCE", "RCR", "FLOW", "PRESSURE"]}
    temp = -10000000.0
    cardiac_cycle_period = temp # initialize as a large negative number for the below code to work correctly
    for location in list(built_in_bc_types.keys()):
        location_segments_of_model = location + "_segments_of_model"
        for i in range(0, len(parameters[location_segments_of_model])):
            segment_number = parameters[location_segments_of_model][i]
            if parameters["boundary_condition_types"][location][segment_number] in built_in_bc_types[location]:

                if parameters["boundary_condition_types"][location][segment_number] == "RESISTANCE":
                    num_variables = 2 # [R, distal_pressure]
                elif parameters["boundary_condition_types"][location][segment_number] == "RCR":
                    num_variables = 3 # [Rp, C, Rd]
                elif parameters["boundary_condition_types"][location][segment_number] in ["FLOW", "PRESSURE"]:
                    num_variables = 1 # [Q] or [P]
                total_num_values = len(parameters["datatable_values"][parameters["boundary_condition_datatable_names"][location][segment_number]])
                if total_num_values == 0:
                    message = "Error. DATATABLE for segment #" + str(segment_number) + " is missing values."
                    raise RuntimeError(message)
                num_values_per_variable = int(total_num_values/num_variables)
                list_of_unsteady_functions = []
                for k in range(num_variables):
                    start_index = k*num_values_per_variable
                    end_index = start_index + num_values_per_variable
                    time, bc_values = extract_bc_time_and_values(start_index, end_index, parameters, segment_number, location)
                    if len(time) > 1: # for time-dependent boundary conditions
                        if cardiac_cycle_period == temp: # this if statement should get called only once
                            if (time[-1] - time[0]) <= 0.0:
                                message = "Error. Boundary condition for segment #" + str(segment_number) + " has a prescribed negative or zero cardiac cycle period."
                                raise RuntimeError(message)
                            else:
                                cardiac_cycle_period = time[-1] - time[0]
                        else:
                            if cardiac_cycle_period != time[-1] - time[0]:
                                message = "Error. The cardiac cycle period is " + str(cardiac_cycle_period) + ", but the boundary condition for segment #" + str(segment_number)+ " is " + str(time[-1] - time[0]) + ". All boundary conditions, including the inlet and outlet boundary conditions, should have the same prescribed cardiac cycle period. Note that each boundary conditions must be prescribed over exactly one cardiac cycle."
                                raise RuntimeError(message)
                    elif len(time) == 1: # this allows for BCs with constant values (steady, time-independent BCs)
                        time.append(time[0] + 1.0)
                        bc_values.append(bc_values[0])
                    list_of_unsteady_functions.append(create_unsteady_bc_value_function(time, bc_values))

                bc_equations[location][segment_number] = list_of_unsteady_functions

    if cardiac_cycle_period == temp:
        cardiac_cycle_period = float(input("\nEnter the period of one cardiac cycle for your 0d simulation. All 0d simulation results will be periodic with a period of this value: "))
        if cardiac_cycle_period <= 0:
            message = "Error. The prescribed cardiac cycle period must be greater than 0."
            raise ValueError(message)
    parameters.update({"cardiac_cycle_period" : cardiac_cycle_period})
    parameters.update({"bc_equations" : bc_equations})

def extract_bc_time_and_values(start_index, end_index, parameters, segment_number, location):
    """
    Purpose:
        Extract the time points and boundary condition values prescribed for the boundary condition of the given segment number at the given location (inlet or outlet of model). Time points and BC values for each location and segment number extracted from the list, parameters["datatable_values"][datatable_name] (obtained from the 0d solver input file), where the time points and BC values alternate in the list.
        datatable_name is given by parameters["boundary_condition_datatable_names"][location][segment_number]
    Inputs:
        int start_index
            = index of parameters["datatable_values"][datatable_name] at which to start the extraction
        int end_index
            = index of parameters["datatable_values"][datatable_name] at which to terminate the extraction
        dict parameters
            -- created from function oneD_to_zeroD_convertor.extract_info_from_solver_input_file
        int segment_number
            = segment number of the 0d vessel element (in the ELEMENT card of the 0d solver input file) for which we want to extract the time points and BC values
        string location
            = "inlet" or "outlet"
            -- specifies whether the BC is an inlet BC or an outlet BC
    Returns:
        list time
            = list of time points prescribed for the BC parameter
        list bc_values
            = list of the BC parameter values
    """
    if location not in ["inlet", "outlet"]:
        message = "Error. 'location' can only be 'inlet' or 'outlet'."
        raise RuntimeError(message)
    time = []
    bc_values = []
    for i in range(start_index, end_index, 2):
        time.append(parameters["datatable_values"][parameters["boundary_condition_datatable_names"][location][segment_number]][i])
        bc_values.append(parameters["datatable_values"][parameters["boundary_condition_datatable_names"][location][segment_number]][i + 1])
    if time[0] != 0.0:
        message = "Error. The initial time for all boundary condition parameters be must 0."
        raise RuntimeError(message)
    return time, bc_values

def create_junction_blocks(parameters, custom_0d_elements_arguments):
    """
    Purpose:
        Create the junction blocks for the 0d model.
    Inputs:
        dict parameters
            -- created from function oneD_to_zeroD_convertor.extract_info_from_solver_input_file
        module custom_0d_elements_arguments
            = module to call custom 0d element arguments from
    Returns:
        void, but updates parameters["blocks"] to include the junction_blocks, where
            parameters["blocks"] = {block_name : block_object}
    """
    junction_blocks = {} # {block_name : block_object}
    for joint_name in list(parameters["junction_types"].keys()):
        if not joint_name.startswith("J") and not joint_name[1].isnumeric():
            message = "Error. Joint name, " + joint_name + ", is not 'J' followed by numeric values. The 0D solver assumes that all joint names are 'J' followed by numeric values in the 0d solver input file. Note that the joint names are the same as the junction names."
            raise RuntimeError(message)
        block_name = joint_name
        connecting_block_list = []
        flow_directions = []
        for segment_number in parameters["joint_name_to_segments_map"][joint_name]["inlet"]:
            connecting_block_list.append("V" + str(segment_number))
            flow_directions.append(-1)
        for segment_number in parameters["joint_name_to_segments_map"][joint_name]["outlet"]:
            connecting_block_list.append("V" + str(segment_number))
            flow_directions.append(+1)
        if (+1 in flow_directions) and (-1 in flow_directions):
            if parameters["junction_types"][joint_name] == "NORMAL_JUNCTION":
                junction_blocks[block_name] = ntwku.Junction(connecting_block_list = connecting_block_list, name = block_name, flow_directions = flow_directions)
            else: # this is a custom, user-defined junction block
                custom_0d_elements_arguments.junction_args[joint_name].update({"connecting_block_list" : connecting_block_list, "flow_directions" : flow_directions, "name" : block_name})
                junction_blocks[block_name] = create_custom_element(parameters["junction_types"][joint_name], custom_0d_elements_arguments.junction_args[joint_name])
        else:
            message = "Error. Junction block, " + block_name + ", must have at least 1 inlet connection and 1 outlet connection."
            raise RuntimeError(message)
    parameters["blocks"].update(junction_blocks)

def get_vessel_block_helpers(parameters):
    """
    Purpose:
        Create helper dictionaries to support the creation of the vessel blocks.
    Inputs:
        dict parameters
            -- created from function oneD_to_zeroD_convertor.extract_info_from_solver_input_file
    Returns:
        dict vessel_blocks_connecting_block_lists
            = {segment_number : connecting_block_list}
        dict vessel_blocks_flow_directions
            = {segment_number : flow_directions}
        dict vessel_blocks_names
            = {segment_number : block_name}
    """
    vessel_blocks_connecting_block_lists = {}
    vessel_blocks_flow_directions = {}
    vessel_blocks_names = {}
    for segment_number in parameters["segment_numbers_list"]:
        vessel_blocks_connecting_block_lists[segment_number] = []
        vessel_blocks_flow_directions[segment_number] = []
        vessel_blocks_names[segment_number] = "V" + str(segment_number)
    for location in ["inlet", "outlet"]:
        for segment_number in parameters[location + "_segments_of_model"]:
            if location == "inlet":
                bc_block_name = "BC" + str(segment_number) + "_inlet"
                vessel_blocks_flow_directions[segment_number].append(-1)
            else:
                bc_block_name = "BC" + str(segment_number) + "_outlet"
                vessel_blocks_flow_directions[segment_number].append(+1)
            vessel_blocks_connecting_block_lists[segment_number].append(bc_block_name)
    for joint_name in list(parameters["junction_types"].keys()):
        if not joint_name.startswith("J") and not joint_name[1].isnumeric():
            message = "Joint name, " + joint_name + ", is not 'J' followed by numeric values. The 0D solver assumes that all joint names are 'J' followed by numeric values in the 0d solver input file. Note that the joint names are the same as the junction names."
            raise RuntimeError(message)
        for location in ["inlet", "outlet"]:
            for segment_number in parameters["joint_name_to_segments_map"][joint_name][location]:
                vessel_blocks_connecting_block_lists[segment_number].append(joint_name)
                if location == "inlet":
                    vessel_blocks_flow_directions[segment_number].append(+1)
                else:
                    vessel_blocks_flow_directions[segment_number].append(-1)
    return vessel_blocks_connecting_block_lists, vessel_blocks_flow_directions, vessel_blocks_names

def create_vessel_blocks(parameters, custom_0d_elements_arguments):
    """
    Purpose:
        Create the vessel blocks for the 0d model.
    Inputs:
        dict parameters
            -- created from function oneD_to_zeroD_convertor.extract_info_from_solver_input_file
        module custom_0d_elements_arguments
            = module to call custom 0d element arguments from
    Returns:
        void, but updates parameters["blocks"] to include the vessel_blocks, where
            parameters["blocks"] = {block_name : block_object}
    """
    vessel_blocks = {} # {block_name : block_object}
    vessel_blocks_connecting_block_lists, vessel_blocks_flow_directions, vessel_blocks_names = get_vessel_block_helpers(parameters)
    for segment_number in parameters["segment_numbers_list"]:
        block_name = vessel_blocks_names[segment_number]
        connecting_block_list = vessel_blocks_connecting_block_lists[segment_number]
        flow_directions = vessel_blocks_flow_directions[segment_number]
        if parameters["segment_0d_types"][segment_number] == "R":
            R = parameters["segment_0d_values"][segment_number][0]
            vessel_blocks[block_name] = ntwku.Resistance(connecting_block_list = connecting_block_list, R = R, name = block_name, flow_directions = flow_directions)
        elif parameters["segment_0d_types"][segment_number] == "C":
            C = parameters["segment_0d_values"][segment_number][0]
            vessel_blocks[block_name] = ntwku.Capacitance(connecting_block_list = connecting_block_list, C = C, name = block_name, flow_directions = flow_directions)
        elif parameters["segment_0d_types"][segment_number] == "L":
            L = parameters["segment_0d_values"][segment_number][0]
            vessel_blocks[block_name] = ntwku.Inductance(connecting_block_list = connecting_block_list, L = L, name = block_name, flow_directions = flow_directions)
        elif parameters["segment_0d_types"][segment_number] == "RC":
            R = parameters["segment_0d_values"][segment_number][0]
            C = parameters["segment_0d_values"][segment_number][1]
            vessel_blocks[block_name] = ntwku.RCBlock(R = R, C = C, connecting_block_list = connecting_block_list, name = block_name, flow_directions = flow_directions)
        elif parameters["segment_0d_types"][segment_number] == "RL":
            R = parameters["segment_0d_values"][segment_number][0]
            L = parameters["segment_0d_values"][segment_number][1]
            vessel_blocks[block_name] = ntwku.RLBlock(R = R, L = L, connecting_block_list = connecting_block_list, name = block_name, flow_directions = flow_directions)
        elif parameters["segment_0d_types"][segment_number] == "RCL":
            R = parameters["segment_0d_values"][segment_number][0]
            C = parameters["segment_0d_values"][segment_number][1]
            L = parameters["segment_0d_values"][segment_number][2]
            vessel_blocks[block_name] = ntwku.RCLBlock(R = R, C = C, L = L, connecting_block_list = connecting_block_list, name = block_name, flow_directions = flow_directions)
        elif parameters["segment_0d_types"][segment_number] == "STENOSIS":
            R = parameters["segment_0d_values"][segment_number][0]
            stenosis_coefficient = parameters["segment_0d_values"][segment_number][1]
            vessel_blocks[block_name] = ntwku.StenosisBlock(R = R, stenosis_coefficient = stenosis_coefficient, connecting_block_list = connecting_block_list, name = block_name, flow_directions = flow_directions)
        else: # this is a custom, user-defined element block
            custom_0d_elements_arguments.vessel_args[segment_number].update({"connecting_block_list" : connecting_block_list, "flow_directions" : flow_directions, "name" : block_name})
            vessel_blocks[block_name] = create_custom_element(parameters["segment_0d_types"][segment_number], custom_0d_elements_arguments.vessel_args[segment_number])
    parameters["blocks"].update(vessel_blocks)

def create_outlet_bc_blocks(parameters, custom_0d_elements_arguments):
    """
    Purpose:
        Create the outlet bc (boundary condition) blocks for the 0d model.
    Inputs:
        dict parameters
            -- created from function oneD_to_zeroD_convertor.extract_info_from_solver_input_file
        module custom_0d_elements_arguments
            = module to call custom 0d element arguments from
    Returns:
        void, but updates parameters["blocks"] to include the outlet_bc_blocks, where
            parameters["blocks"] = {block_name : block_object}
    """
    outlet_bc_blocks = {} # {block_name : block_object}
    for segment_number in parameters["outlet_segments_of_model"]:
        if parameters["boundary_condition_types"]["outlet"][segment_number] != "NOBOUND":
            block_name = "BC" + str(segment_number) + "_outlet"
            connecting_block_list = ["V" + str(segment_number)]
            flow_directions = [-1]
            if parameters["boundary_condition_types"]["outlet"][segment_number] == "RESISTANCE":
                Rfunc = parameters["bc_equations"]["outlet"][segment_number][0]
                Pref_func = parameters["bc_equations"]["outlet"][segment_number][1]
                outlet_bc_blocks[block_name] = ntwku.UnsteadyResistanceWithDistalPressure(connecting_block_list = connecting_block_list, Rfunc = Rfunc, Pref_func = Pref_func, name = block_name, flow_directions = flow_directions)
            elif parameters["boundary_condition_types"]["outlet"][segment_number] == "RCR":
                Rp_func = parameters["bc_equations"]["outlet"][segment_number][0]
                C_func = parameters["bc_equations"]["outlet"][segment_number][1]
                Rd_func = parameters["bc_equations"]["outlet"][segment_number][2]
                Pref = 0
                Pref_func = create_unsteady_bc_value_function([0.0, 1.0], [Pref, Pref])
                outlet_bc_blocks[block_name] = ntwku.UnsteadyRCRBlockWithDistalPressure(Rp_func = Rp_func, C_func = C_func, Rd_func = Rd_func, Pref_func = Pref_func, connecting_block_list = connecting_block_list, name = block_name, flow_directions = flow_directions)
            elif parameters["boundary_condition_types"]["outlet"][segment_number] == "FLOW":
                Qfunc = parameters["bc_equations"]["outlet"][segment_number][0]
                outlet_bc_blocks[block_name] = ntwku.UnsteadyFlowRef(connecting_block_list = connecting_block_list, Qfunc = Qfunc, name = block_name, flow_directions = flow_directions)
            elif parameters["boundary_condition_types"]["outlet"][segment_number] == "PRESSURE":
                Pfunc = parameters["bc_equations"]["outlet"][segment_number][0]
                outlet_bc_blocks[block_name] = ntwku.UnsteadyPressureRef(connecting_block_list = connecting_block_list, Pfunc = Pfunc, name = block_name, flow_directions = flow_directions)
            elif parameters["boundary_condition_types"]["outlet"][segment_number] == "CORONARY":
                "Publication reference: Kim, H. J. et al. Patient-specific modeling of blood flow and pressure in human coronary arteries. Annals of Biomedical Engineering 38, 3195–3209 (2010)."
                Ra1 = parameters["datatable_values"][parameters["boundary_condition_datatable_names"]["outlet"][segment_number]][1]
                Ra2 = parameters["datatable_values"][parameters["boundary_condition_datatable_names"]["outlet"][segment_number]][3]
                Ca = parameters["datatable_values"][parameters["boundary_condition_datatable_names"]["outlet"][segment_number]][5]
                Cc = parameters["datatable_values"][parameters["boundary_condition_datatable_names"]["outlet"][segment_number]][7]
                Rv1 = parameters["datatable_values"][parameters["boundary_condition_datatable_names"]["outlet"][segment_number]][9]
                Pv_distal_pressure = parameters["datatable_values"][parameters["boundary_condition_datatable_names"]["outlet"][segment_number]][11]

                time_of_intramyocardial_pressure, bc_values_of_intramyocardial_pressure = extract_bc_time_and_values(start_index = 12, end_index = len(parameters["datatable_values"][parameters["boundary_condition_datatable_names"]["outlet"][segment_number]]), parameters = parameters, segment_number = segment_number, location = "outlet")

                if "cardiac_cycle_period" in parameters:
                    if time_of_intramyocardial_pressure[-1] - time_of_intramyocardial_pressure[0] != parameters["cardiac_cycle_period"]:
                        message = "Error. The time history of the intramyocadial pressure for the coronary boundary condition for segment #" + str(segment_number) + " does not have the same cardiac cycle period as the other boundary conditions.  All boundary conditions, including the inlet and outlet boundary conditions, should have the same prescribed cardiac cycle period. Note that each boundary conditions must be prescribed over exactly one cardiac cycle."
                        raise RuntimeError(message)
                else:
                    parameters.update({"cardiac_cycle_period": time_of_intramyocardial_pressure[-1] - time_of_intramyocardial_pressure[0]})

                # time_derivative_of_intramyocardial_pressure = np.zeros((len(time_of_intramyocardial_pressure), 2))
                # time_derivative_of_intramyocardial_pressure[:, 0] = time_of_intramyocardial_pressure
                # deltat_temp = time_of_intramyocardial_pressure[1] - time_of_intramyocardial_pressure[0]
                # extended_time = np.array([time_of_intramyocardial_pressure[0] - deltat_temp] + time_of_intramyocardial_pressure + [time_of_intramyocardial_pressure[-1] + deltat_temp]) # extend time beyond the first and last boundaries to get a periodic dPvdt
                # bc_values_of_intramyocardial_pressure_extended = np.array([bc_values_of_intramyocardial_pressure[-2]] + bc_values_of_intramyocardial_pressure + [bc_values_of_intramyocardial_pressure[1]])
                # intramyocardial_pressure_derivative = np.gradient(bc_values_of_intramyocardial_pressure_extended, extended_time)
                # time_derivative_of_intramyocardial_pressure[:, 1] = intramyocardial_pressure_derivative[1:-1]

                # plt.figure()
                # # plt.plot(time_derivative_of_intramyocardial_pressure[:, 0], time_derivative_of_intramyocardial_pressure[:, 1], 'b--', label = 'dPimdt')
                # ttt = np.zeros(len(time_of_intramyocardial_pressure)*2)
                # ttt[:len(time_of_intramyocardial_pressure)] = np.array(time_of_intramyocardial_pressure)
                # ttt[len(time_of_intramyocardial_pressure):] = np.array(time_of_intramyocardial_pressure) + time_of_intramyocardial_pressure[-1]
                # yyy = np.zeros(len(time_of_intramyocardial_pressure)*2)
                # yyy[:len(time_of_intramyocardial_pressure)] = np.array(bc_values_of_intramyocardial_pressure)
                # yyy[len(time_of_intramyocardial_pressure):] = np.array(bc_values_of_intramyocardial_pressure)
                # # plt.plot(time_of_intramyocardial_pressure, bc_values_of_intramyocardial_pressure, 'k-', label = 'Pim')
                # plt.plot(ttt, yyy)
                # plt.legend()
                # plt.show()

                # outlet_bc_blocks[block_name] = ntwku.OpenLoopCoronaryWithDistalPressureBlock(R1 = Ra1, C1 = Ca, R2 = Ra2, C2 = Cc, R3 = Rv1, dPvdt_f = time_derivative_of_intramyocardial_pressure, Pv = Pv_distal_pressure, cardiac_cycle_period = parameters["cardiac_cycle_period"], connecting_block_list = connecting_block_list, name = block_name, flow_directions = flow_directions)

                # Pim_func = lambda t: np.interp(t, time_of_intramyocardial_pressure, bc_values_of_intramyocardial_pressure)
                #
                # Pv_distal_pressure_func = lambda t: Pv_distal_pressure + (t - t)

                # Pv_distal_pressure_func= lambda t: np.interp(t, time_of_intramyocardial_pressure, bc_values_of_intramyocardial_pressure)

                # Pim_func = scipy.interpolate.CubicSpline(time_of_intramyocardial_pressure, bc_values_of_intramyocardial_pressure, bc_type = 'periodic')

                Pim_func = np.zeros((len(time_of_intramyocardial_pressure), 2))
                Pv_distal_pressure_func = np.zeros((len(time_of_intramyocardial_pressure), 2))

                Pim_func[:, 0] = time_of_intramyocardial_pressure
                Pv_distal_pressure_func[:, 0] = time_of_intramyocardial_pressure

                Pim_func[:, 1] = bc_values_of_intramyocardial_pressure
                Pv_distal_pressure_func[:, 1] = np.ones(len(time_of_intramyocardial_pressure))*Pv_distal_pressure
                # Pv_distal_pressure_func[:, 1] = bc_values_of_intramyocardial_pressure
                # plt.figure()
                # plt.plot(Pim_func[:, 0], Pim_func[:, 1], label = "Pim")
                # plt.plot(Pv_distal_pressure_func[:, 0], Pv_distal_pressure_func[:, 1], label = "Pv")
                # plt.legend()
                # plt.show()

                outlet_bc_blocks[block_name] = ntwku.OpenLoopCoronaryWithDistalPressureBlock_v2(Ra = Ra1, Ca = Ca, Ram = Ra2, Cim = Cc, Rv = Rv1, Pim = Pim_func, Pv = Pv_distal_pressure_func, cardiac_cycle_period = parameters["cardiac_cycle_period"], connecting_block_list = connecting_block_list, name = block_name, flow_directions = flow_directions)


            else: # this is a custom, user-defined outlet bc block
                custom_0d_elements_arguments.outlet_bc_args[segment_number].update({"connecting_block_list" : connecting_block_list, "flow_directions" : flow_directions, "name" : block_name})
                outlet_bc_blocks[block_name] = create_custom_element(parameters["boundary_condition_types"]["outlet"][segment_number], custom_0d_elements_arguments.outlet_bc_args[segment_number])
    parameters["blocks"].update(outlet_bc_blocks)

def create_inlet_bc_blocks(parameters, custom_0d_elements_arguments):
    """
    Purpose:
        Create the inlet bc (boundary condition) blocks for the 0d model.
    Inputs:
        dict parameters
            -- created from function oneD_to_zeroD_convertor.extract_info_from_solver_input_file
        module custom_0d_elements_arguments
            = module to call custom 0d element arguments from
    Returns:
        void, but updates parameters["blocks"] to include the inlet_bc_blocks, where
            parameters["blocks"] = {block_name : block_object}
    """
    inlet_bc_blocks = {} # {block_name : block_object}
    for segment_number in parameters["inlet_segments_of_model"]:
        block_name = "BC" + str(segment_number) + "_inlet"
        connecting_block_list = ["V" + str(segment_number)]
        flow_directions = [+1]
        if parameters["boundary_condition_types"]["inlet"][segment_number] == "FLOW":
            func = parameters["bc_equations"]["inlet"][segment_number][0]
            inlet_bc_blocks[block_name] = ntwku.UnsteadyFlowRef(Qfunc = func, connecting_block_list = connecting_block_list, name = block_name, flow_directions = flow_directions)
        elif parameters["boundary_condition_types"]["inlet"][segment_number] == "PRESSURE":
            func = parameters["bc_equations"]["inlet"][segment_number][0]
            inlet_bc_blocks[block_name] = ntwku.UnsteadyPressureRef(Pfunc = func, connecting_block_list = connecting_block_list, name = block_name, flow_directions = flow_directions)
        else: # this is a custom, user-defined inlet bc block
            custom_0d_elements_arguments.inlet_bc_args[segment_number].update({"connecting_block_list" : connecting_block_list, "flow_directions" : flow_directions, "name" : block_name})
            inlet_bc_blocks[block_name] = create_custom_element(parameters["boundary_condition_types"]["inlet"][segment_number], custom_0d_elements_arguments.inlet_bc_args[segment_number])
    parameters["blocks"].update(inlet_bc_blocks)

def create_LPN_blocks(parameters, custom_0d_elements_arguments):
    """
    Purpose:
        Create all LPNBlock objects for the 0d model.
    Inputs:
        dict parameters
            -- created from function oneD_to_zeroD_convertor.extract_info_from_solver_input_file
        module custom_0d_elements_arguments
            = module to call custom 0d element arguments from
    Returns:
        void, but updates parameters to include:
            dict blocks
                = {block_name : block_object}
                -- where the values are the junction, vessel, outlet BC, and inlet BC block objects
    """
    blocks = {}  # {block_name : block_object}
    parameters.update({"blocks" : blocks})
    create_junction_blocks(parameters, custom_0d_elements_arguments)
    create_vessel_blocks(parameters, custom_0d_elements_arguments)
    create_bc_equations(parameters)
    create_outlet_bc_blocks(parameters, custom_0d_elements_arguments)
    create_inlet_bc_blocks(parameters, custom_0d_elements_arguments)
    parameters.update({"block_names" : list(parameters["blocks"].keys())})

def set_solver_parameters(parameters):
    """
    Purpose:
        Set the 0d simulation time-stepping parameters
    Inputs:
        dict parameters
            -- created from function oneD_to_zeroD_convertor.extract_info_from_solver_input_file
    Returns:
        void, but updates parameters to include:
            int number_of_time_pts_per_cardiac_cycle
                = number of time steps to simulate for each cardiac cycle
            int number_of_cardiac_cycles
                = number of cardiac cycles to simulate
            float delta_t
                = constant time step size for the 0d simulation
            int total_number_of_simulated_time_steps
                = total number of time steps to simulate for the entire 0d simulation
    """
    delta_t = parameters["cardiac_cycle_period"]/(parameters["number_of_time_pts_per_cardiac_cycle"] - 1)
    parameters.update({"delta_t" : delta_t})
    total_number_of_simulated_time_steps = int((parameters["number_of_time_pts_per_cardiac_cycle"] - 1)*parameters["number_of_cardiac_cycles"] + 1)
    parameters.update({"total_number_of_simulated_time_steps" : total_number_of_simulated_time_steps})

def load_in_ics(var_name_list, ICs_dict):

    var_name_list_loaded = ICs_dict["var_name_list"]
    y_initial_loaded = ICs_dict["y"]
    ydot_initial_loaded = ICs_dict["ydot"]

    y_initial = np.zeros(len(var_name_list))
    ydot_initial = np.zeros(len(var_name_list))
    for i in range(len(var_name_list)):
        var_name = var_name_list[i]
        ind = var_name_list_loaded.index(var_name)
        y_initial[i] = y_initial_loaded[ind]
        ydot_initial[i] = ydot_initial_loaded[ind]

    return y_initial, ydot_initial

# @profile
def run_network_util(zero_d_solver_input_file_path, parameters, draw_directed_graph, use_ICs_from_npy_file, ICs_npy_file_path, save_y_ydot_to_npy, y_ydot_file_path, simulation_start_time):
    """
    Purpose:
        Run functions from network_util_NR to execute the 0d simulation and generate simulation results (pressure, flow rates, and wall shear stress).
    Inputs:
        string zero_d_solver_input_file_path
            = path to the 0d solver input file
        dict parameters
            -- created from function oneD_to_zeroD_convertor.extract_info_from_solver_input_file
        boolean draw_directed_graph
            = True to visualize the 0d model as a directed graph using networkx -- saves the graph to a .png file (hierarchical graph layout) and a networkx .dot file; False, otherwise. .dot file can be opened with neato from graphviz to visualize the directed in a different format.
    Returns:
        np.array zero_d_time
            = np.array of simulated time points
        np.array results_0d
            = np.array of the 0d simulation results, where the rows correspond to the each time point in zero_d_time and the column j corresponds to the solution for item j in var_name_list
        list var_name_list
            = list of the names of the 0d simulation results; most of the items in var_name_list are the QoIs + the names of the wires used in the 0d model (the wires connect the 0d blocks), where the wire names are usually comprised of the wire's inlet block name + "_" + the wire's outlet block name
                Example:
                    for var_name_list = ['P_V6_BC6_outlet', 'Q_V6_BC6_outlet'], then results_0d[:, i] holds the pressure (i = 0) or flow rate simulation result (i = 1) (both as np.arrays) for wire R6_BC6, which corresponds to cap segment #6
    """

    block_list = list(parameters["blocks"].values())
    connect_list, wire_dict = connections.connect_blocks_by_inblock_list(block_list)
    if draw_directed_graph == True:
        zero_d_input_file_name, zero_d_input_file_extension = os.path.splitext(zero_d_solver_input_file_path)
        directed_graph_file_path = zero_d_input_file_name + "_directed_graph"
        save_directed_graph(block_list, connect_list, directed_graph_file_path)
    neq = connections.compute_neq(block_list, wire_dict) # number of equations governing the 0d model
    for block in block_list: # run a consistency check
        connections.check_block_connection(block)
    var_name_list = connections.assign_global_ids(block_list, wire_dict) # assign solution variables with global ID

    # initialize solution structures
    if use_ICs_from_npy_file:
        ICs_dict = np.load(ICs_npy_file_path, allow_pickle = True).item()
        y_initial, ydot_initial = load_in_ics(var_name_list, ICs_dict)
    else:
        y_initial, ydot_initial = connections.initialize_solution_structures(neq) # initial conditions for all solutions are zero
    y_next = y_initial.copy()
    ydot_next = ydot_initial.copy()

    rho = 0.1
    args = {}
    args['Time step'] = parameters["delta_t"]
    args['rho'] = rho
    args['Wire dictionary'] = wire_dict
    args["check_jacobian"] = parameters["check_jacobian"]

    # y_next, ydot_next = min_ydot_least_sq_init(neq, 1e-8, y_initial, block_list, args, parameters["delta_t"], rho)

    print("starting simulation")

    ylist = [y_next.copy()]
    parameters["initial_time"] = simulation_start_time
    tlist = np.array([ parameters["initial_time"] + _*parameters["delta_t"] for _ in range(0, parameters["total_number_of_simulated_time_steps"])])

    # create time integration
    t_int = time_int.GenAlpha(rho, y_next)

    for t_current in tqdm(tlist[:-1]): # added this line on 10/14/20:
        args['Solution'] = y_next
        y_next, ydot_next = t_int.step(y_next, ydot_next, t_current, block_list, args, parameters["delta_t"])
        ylist.append(y_next)

    if save_y_ydot_to_npy:
        np.save(y_ydot_file_path, {"y" : y_next, "ydot" : ydot_next, "var_name_list" : var_name_list})
        print("var_name_list = ", var_name_list)

    results_0d = np.array(ylist)
    ylist, var_name_list = compute_wss(parameters, results_0d, ylist, var_name_list)
    results_0d = np.array(ylist)
    zero_d_time = tlist
    # print("var_name_list = ", var_name_list)
    return zero_d_time, results_0d, var_name_list

def compute_wss(parameters, results_0d, ylist, var_name_list): # currently here 8/13/20: need to recheck that this function correctly computes wss...
    """
    Purpose:
        Compute the wall shear stress (wss) using Poiseuille flow. Compute wss only for all blocks' inlet and outlet wires that have both flow and pressure results.
        wss = tau = 4.0*mu*Q/np.pi/R^3
    Inputs:
        dict parameters
            -- created from function oneD_to_zeroD_convertor.extract_info_from_solver_input_file
        np.array results_0d
            = np.array of the 0d simulation results, where the rows correspond to the each simulated time point and column j corresponds to the 0d solution for the solution variable name in var_name_list[j]
        list ylist
            = list version of results_0d, where ylist[i] is an np.array of the solutions for all solution variable names in var_name_list at time step #i
        list var_name_list
            = list of the 0d simulation results' solution variable names; most of the items in var_name_list are the QoIs + the names of the wires used in the 0d model (the wires connecting the 0d LPNBlock objects), where the wire names are usually comprised of the wire's inlet block name + "_" + the wire's outlet block name
                Example:
                    for var_name_list = ['P_V6_BC6_outlet', 'Q_V6_BC6_outlet'], then results_0d[:, i] holds the pressure (i = 0) or flow rate simulation result (i = 1) (both as np.arrays) for wire R6_BC6_outlet. This wire connects a resistance vessel block to an outlet BC block (specifically for vessel segment #6)
    Returns:
        list ylist, but updated to include the wss solutions at all simulated time points
        list var_name_list, but updated to include the wss solution variable names
            -- the wss solution variable names will be of the form: 'tau' + wire_name
    """
    length_of_original_var_name_list = len(var_name_list)
    for i in range(length_of_original_var_name_list):
        wire_name = var_name_list[i][2:]
        if "Q_" + wire_name in var_name_list and "P_" + wire_name in var_name_list and "tau_" + wire_name not in var_name_list and "var" not in wire_name:
            blocks_attached_to_wire = wire_name.split("_")
            if "BC" in wire_name:
                counter = 2
                if blocks_attached_to_wire[0].startswith("BC") and blocks_attached_to_wire[0][2:].isnumeric() and blocks_attached_to_wire[1] == "inlet": # at inlet segment
                    place = 0 # get the segment number of the inlet segment
                elif blocks_attached_to_wire[-2].startswith("BC") and blocks_attached_to_wire[-2][2:].isnumeric() and blocks_attached_to_wire[-1] == "outlet": # at outlet segment
                    place = -2 # get the segment number of the outlet segment
            elif "J" in wire_name:
                if blocks_attached_to_wire[0].startswith("J") and blocks_attached_to_wire[0][1:].isnumeric(): # at an outlet wire of a junction block
                    place = -1 # get the segment number of the daughter vessel
                elif blocks_attached_to_wire[-1].startswith("J") and blocks_attached_to_wire[-1][1:].isnumeric(): # at an inlet wire of a junction block
                    place = 0 # get the segment number of the parent vessel
                counter = 1
            segment_number = int(blocks_attached_to_wire[place][counter:])
            radius = parameters["radii"][segment_number]
            mu = parameters["mu"][segment_number]
            Q_index = var_name_list.index("Q_" + wire_name)
            tau = 4.0*mu*results_0d[:, Q_index]/(np.pi*radius**3)
            var_name_list.append("tau_" + wire_name)
            for j in range(len(ylist)):
                ylist[j] = list(ylist[j]) # results at time step #j # this line is needed because ylist[j] is an np.array, so we have to turn it into a list before we can append tau[j] to it
                ylist[j].append(tau[j])
    return ylist, var_name_list

def extract_0d_cap_results(zero_d_results_for_var_names): # currently here 8/18/20: need to run some simple test cases to make sure that this function works correctly
    """
    Purpose:
        Extract the 0d simulation results (pressure, flow rate, and wall shear stress) at only the model's caps.
    Inputs:
        dict zero_d_results_for_var_names
            =   {
                    "time" : np.array of simulated time points,

                    "flow" : {var_name : np.array of flow rate,

                    "pressure" : {var_name : np.array of pressure},

                    "wss" : {var_name : np.array of wall shear stress},

                    "internal" : {var_name : np.array of internal block solutions},

                        where var_name is an item in var_name_list (var_name_list generated from run_network_util)
                }
    Returns:
        dict zero_d_cap_results
            =   {
                    "inlet" :   {
                                    "flow" : {segment number : np.array of simulation results},

                                    "pressure" : {segment number : np.array of simulation results},

                                    "wss" : {segment number : np.array of simulation results}
                                }

                    "outlet" :   {
                                    "flow" : {segment number : np.array of simulation results},

                                    "pressure" : {segment number : np.array of simulation results},

                                    "wss" : {segment number : np.array of simulation results}
                                }

                    "time" : np.array of simulation time points
                }
    """
    zero_d_inlet_cap_results = {}
    zero_d_outlet_cap_results = {}
    for qoi in ["pressure", "flow", "wss"]:
        if qoi == "flow":
            solution_type = "Q_"
        elif qoi == "pressure":
            solution_type = "P_"
        elif qoi == "wss":
            solution_type = "tau_"
        zero_d_inlet_cap_results[qoi] = {}
        zero_d_outlet_cap_results[qoi] = {}
        var_names = list(zero_d_results_for_var_names[qoi].keys())
        for i in range(len(var_names)):
            if ("BC" in var_names[i]) and ("var" not in var_names[i]) and (var_names[i].startswith(solution_type)): # at a cap segment # currently here 12/3/20: should check to see if var_names[i] also includes "V" here
                var_names_split = var_names[i].split("_")
                index = [i for i, s in enumerate(var_names_split) if "BC" in s][0]
                segment_number = int(var_names_split[index][2:])
                if "inlet" in var_names[i]: # this is an inlet cap
                    zero_d_inlet_cap_results[qoi][segment_number] = zero_d_results_for_var_names[qoi][var_names[i]]
                else: # this is an outlet cap # currently here 9/24/20: should add a "elif "outlet" in var_names[i]
                    zero_d_outlet_cap_results[qoi][segment_number] = zero_d_results_for_var_names[qoi][var_names[i]]
    zero_d_cap_results = {"inlet" : zero_d_inlet_cap_results, "outlet" : zero_d_outlet_cap_results, "time" : zero_d_results_for_var_names["time"]}
    return zero_d_cap_results

def reformat_network_util_results(zero_d_time, results_0d, var_name_list):
    """
    Purpose:
        Reformat the 0d simulation results (results_0d) into a dictionary (zero_d_results_for_var_names)
    Inputs:
        np.array zero_d_time
            = np.array of simulated time points
        np.array results_0d
            = np.array of the 0d simulation results, where the rows correspond to the each simulated time point and column j corresponds to the 0d solution for the solution variable name in var_name_list[j]
        list var_name_list
            = list of the 0d simulation results' solution variable names; most of the items in var_name_list are the QoIs + the names of the wires used in the 0d model (the wires connecting the 0d LPNBlock objects), where the wire names are usually comprised of the wire's inlet block name + "_" + the wire's outlet block name
                Example:
                    for var_name_list = ['P_V6_BC6_outlet', 'Q_V6_BC6_outlet'], then results_0d[:, i] holds the pressure (i = 0) or flow rate simulation result (i = 1) (both as np.arrays) for wire R6_BC6_outlet. This wire connects a resistance vessel block to an outlet BC block (specifically for vessel segment #6)
    Returns:
        dict zero_d_results_for_var_names
            =   {
                    "time" : np.array of simulated time points,

                    "flow" : {var_name : np.array of flow rate,

                    "pressure" : {var_name : np.array of pressure},

                    "wss" : {var_name : np.array of wall shear stress},

                    "internal" : {var_name : np.array of internal block solution},

                        where var_name is an item in var_name_list (var_name_list generated from run_network_util)
                }
    """
    zero_d_results_for_var_names = {"flow" : {}, "pressure" : {}, "wss" : {}, "time" : zero_d_time, "internal" : {}}
    for i in range(len(var_name_list)):
        var_name = var_name_list[i]
        res = results_0d[:, i]
        if var_name.startswith("Q_"):
            zero_d_results_for_var_names["flow"][var_name] = res
        elif var_name.startswith("P_"):
            zero_d_results_for_var_names["pressure"][var_name] = res
        elif var_name.startswith("tau_"):
            zero_d_results_for_var_names["wss"][var_name] = res
        elif var_name.startswith("var_"):
            zero_d_results_for_var_names["internal"][var_name] = res
        else:
            message = "Error. There are unaccounted for solution variables here, for var_name = " + var_name
            raise RuntimeError(message)
    return zero_d_results_for_var_names

def extract_last_cardiac_cycle_simulation_results(time, results, number_of_time_pts_per_cardiac_cycle):
    """
    Purpose:
        Extract the simulation results for the last cardiac cycle for the given np.arrays, time and results
    Inputs:
        np.array time
            = array of time points
        np.array results
            = array of simulation results
        int number_of_time_pts_per_cardiac_cycle
            = number of simulated time points per cycle
    Returns:
        np.array time_for_last_cardiac_cycle
            = time, but condensed to contain just the values for the last cardiac cycle
        np.array results_for_last_cardiac_cycle
            = results, but condensed to contain just the values for the last cardiac cycle
    """
    time_for_last_cardiac_cycle = time[-1*number_of_time_pts_per_cardiac_cycle:]
    results_for_last_cardiac_cycle =  results[-1*number_of_time_pts_per_cardiac_cycle:]
    time_for_last_cardiac_cycle = time_for_last_cardiac_cycle - time_for_last_cardiac_cycle[0] + time[0]
    return time_for_last_cardiac_cycle, results_for_last_cardiac_cycle

def run_last_cycle_extraction_routines(cardiac_cycle_period, number_of_time_pts_per_cardiac_cycle, zero_d_results_for_var_names):
    """
    Purpose:
        Extract the last cardiac cycle waveform for all 0d simulation results
    Inputs:
        float cardiac_cycle_period
            = period of a cardiac cycle
        int number_of_time_pts_per_cardiac_cycle
            = number of simulated 0d time points per cycle
        dict zero_d_results_for_var_names
            =   {
                    "time" : np.array of simulated time points,

                    "flow" : {var_name : np.array of flow rate,

                    "pressure" : {var_name : np.array of pressure},

                    "wss" : {var_name : np.array of wall shear stress},

                    "internal" : {var_name : np.array of internal block solutions},

                        where var_name is an item in var_name_list (var_name_list generated from run_network_util)
                }
    Returns:
        dict zero_d_results_for_var_names_last_cycle
            =   zero_d_results_for_var_names, but for just the last simulated cardiac cycle
    """
    zero_d_results_for_var_names_last_cycle = {}
    for qoi in list(zero_d_results_for_var_names.keys()):
        if qoi != "time":
            zero_d_results_for_var_names_last_cycle[qoi] = {}
            for var_name in list(zero_d_results_for_var_names[qoi].keys()):
                time_for_last_cardiac_cycle, res = extract_last_cardiac_cycle_simulation_results(zero_d_results_for_var_names["time"], zero_d_results_for_var_names[qoi][var_name], number_of_time_pts_per_cardiac_cycle)
                zero_d_results_for_var_names_last_cycle[qoi][var_name] = res
    zero_d_results_for_var_names_last_cycle["time"] = time_for_last_cardiac_cycle

    # check that the cardiac cycle period is correctly given by zero_d_time
    errr = (cardiac_cycle_period - (zero_d_results_for_var_names_last_cycle["time"][-1] - zero_d_results_for_var_names_last_cycle["time"][0]))/cardiac_cycle_period*100.0
    if errr > 0.01:
        message = 'Error. cardiac_cycle_period != zero_d_results_for_var_names_last_cycle["time"][-1] - zero_d_results_for_var_names_last_cycle["time"][0]'
        raise RuntimeError(message)
    return zero_d_results_for_var_names_last_cycle

def save_directed_graph(block_list, connect_list, directed_graph_file_path):
    """
    Purpose:
        Visualize the 0d model as a directed graph -- save the graph in a hierarchical graph layout to a .png file; also save a networkx .dot file that can be opened with neato via graphviz to visualize the graph in a different layout.
    Inputs:
        list block_list
            = [list of all of the 0d LPNBlock objects]
        list connect_list
            = [list of (blockA_index, blockB_index)]
                where blockA and blockB are connected to each other, and blockA_index and blockB_index are the index locations at which the blockA and blockB objects are stored in block_list
        string directed_graph_file_path
            = name of the hierarchical graph .png file and networkx .dot file that will be saved
    Returns:
        void, but saves a .png file visualizing the 0d model as a directed graph, as well as a networkx .dot file that can be opened with neato via graphviz to visualize the graph in a different layout
    """
    plt.figure(figsize=(20, 11))
    G = nx.DiGraph()
    G.add_edges_from([(block_list[tpl[0]].name, block_list[tpl[1]].name) for tpl in connect_list])
    # nx.draw(G)
    # pos = nx.planar_layout(G)
    # pos = nx.spring_layout(G)
    pos = nx.nx_pydot.pydot_layout(G, prog='dot') # see other graph layouts: https://stackoverflow.com/questions/21978487/improving-python-networkx-graph-layout
    nx.draw_networkx_nodes(G,pos)
    nx.draw_networkx_labels(G,pos)
    nx.draw_networkx_edges(G,pos)

    plt.tight_layout()
    plt.savefig(directed_graph_file_path + ".png", format="PNG")
    # plt.show()
    nx.nx_pydot.write_dot(G, directed_graph_file_path + ".dot")
    plt.close("all")

    # To post-process the dot file:
    # You need graphviz installed
    # Then you could use one of its' native post-processors (neato.exe) for dot files from the command line:
    # neato.exe -Tpng test.dot -Gsplines=ortho -Gnodesep=1 -Goverlap=scale -o test.png
    # or maybe try: neato -Tpng -Gstart=7 -Gepsilon=.0000001 <name of dot file> -Goverlap=scale -Gsplines=true -Gdpi=300 -o <name of png file to save figure to> (the former command might yield overlaps
    # or: neato -Tpng -Gstart=7 0003_0001_lpn.dot -o 0003_0001_lpn.png
    # see the following references for other options:
    #       https://stackoverflow.com/questions/3967600/how-to-prevent-edges-in-graphviz-to-overlap-each-other
    #       https://www.graphviz.org/pdf/neatoguide.pdf

    # Gsplines=ortho ensures that connectors are horizontal or vertical

    # neato avoids overlaps, though there are a number of other post-processors to try

    # To extract coordinates and structure from the graphviz post-processor, use:
    # neato.exe test.dot -Gsplines=ortho -Gnodesep=1 -Goverlap=scale

def collapse_inlet_and_outlet_0d_results(zero_d_cap_results): # currently here 8/18/20: need to run some simple test cases to make sure that this function works correctly
    """
    Purpose:
        Collapse zero_d_cap_results into a simpler dictionary, zero_d_cap_results_collapsed. This collapse is performed only for 0d models with more than one vessel element (listed in the ELEMENT card of 0d solver input file).
    Caveat:
        This function effectively combines the inlet and outlet solutions into a single dictionary. As such, it assumes that inlet caps and outlet caps have unique segment numbers, which only occurs for 0d models with more than one vessel element.
    Inputs:
        dict zero_d_cap_results
            =   {
                    "inlet" :   {
                                    "flow" : {segment number : np.array of simulation results},

                                    "pressure" : {segment number : np.array of simulation results},

                                    "wss" : {segment number : np.array of simulation results}
                                }

                    "outlet" :   {
                                    "flow" : {segment number : np.array of simulation results},

                                    "pressure" : {segment number : np.array of simulation results},

                                    "wss" : {segment number : np.array of simulation results}
                                }

                    "time" : np.array of simulation time points
                }
    Returns:
        If the inlet and outlet caps have unique segment numbers, meaning the 0d model has more than one element in the 0d solver input file, then return:
            dict zero_d_cap_results_collapsed
                =   {
                        "flow" : {segment number : np.array of simulation results},

                        "pressure" : {segment number : np.array of simulation results},

                        "wss" : {segment number : np.array of simulation results},

                        "time" : np.array of simulation time points
                    }
        otherwise, return:
            dict zero_d_cap_results
    """
    all_segment_numbers = {"inlet" : [], "outlet" :[]}
    cap_types = list(all_segment_numbers.keys())
    for cap_type in cap_types:
        for qoi in list(zero_d_cap_results[cap_type].keys()):
            for segment_number in list(zero_d_cap_results[cap_type][qoi].keys()):
                if segment_number not in all_segment_numbers[cap_type]:
                    all_segment_numbers[cap_type].append(segment_number)
    if len(list(set(all_segment_numbers["inlet"]) & set(all_segment_numbers["outlet"]))) == 0: # collapse zero_d_cap_results only if the inlet and and outlet caps have unique segment numbers; non-unique segment numbers means that model has a single element/segment that has both the inlet and outlet cap
        zero_d_cap_results_collapsed = {"time" : zero_d_cap_results["time"]}
        for qoi in list(zero_d_cap_results["outlet"].keys()):
            zero_d_cap_results_collapsed[qoi] = {}
            for cap_type in ["inlet", "outlet"]:
                for segment_number in list(zero_d_cap_results[cap_type][qoi].keys()):
                    zero_d_cap_results_collapsed[qoi][segment_number] = zero_d_cap_results[cap_type][qoi][segment_number]
        return zero_d_cap_results_collapsed
    else:
        print("Warning. Model has only a single element/segment. Collapse is not possible. Returning the original input dictionary.\n")
        return zero_d_cap_results

def check_convergence_of_simulation_result_gradient(time, result, number_of_values_to_check = 20, tolerance = 0.01):
    """
    Purpose:
        Check if result has converged to a periodic state. result has converged to a periodic state if the last number_of_values_to_check values of the temporal gradient of result are within the given tolerance and if their signs do not fluctuate during the last number_of_values_to_check values.
    Caveat:
        result must be a time-averaged waveform
    Inputs:
        np.array time
            = np.array of time points
        np.array result
            = np.array of result values for which we want to check convergence
            -- must have the same length as time
        int number_of_values_to_check
            -- check the last number_of_values_to_check values in the temporal gradient of result for convergence
        float tolerance
            = convergence requires that the last number_of_values_to_check values of the temporal gradient of result are equal to or below the given tolerance
    Returns:
        boolean converged
            = True if result has converged; otherwise, False
            -- True requires fluctuating_gradient to be False
    """
    converged = False
    drdt = np.gradient(result, time) # temporal gradient of result
    counter = 0
    fluctuation_counter = 0
    fluctuating_gradient = False
    start_index = -1*number_of_values_to_check
    end_index = 0
    for i in range(start_index, end_index):
        if np.abs(drdt[i]) > tolerance:
            # print("drdt[i] = ", drdt[i])
            print("Warning: the gradient, " + str(np.abs(drdt[i])) + " of the result is not within the given tolerance, " + str(tolerance) + ". Results not converged yet.") # currently here 7/21/20: maybe i should put these warnings in a separate text file too, like an output file that logs all of the warnings and issues, so that the user can easily refer to all the errors in a single consistent coherent doc
        else:
            counter =  counter + 1
        previous_drdt = drdt[i]
        if i > start_index:
            if np.sign(drdt[i]) != np.sign(previous_drdt):
                fluctuation_counter = fluctuation_counter + 1
                print("Warning: the last " + number_of_values_to_check + " values of the gradient of the result exhibit fluctuating signs. This means that the results haven't converged.") # currently here 7/21/20: maybe i should put these warnings in a separate text file too, like an output file that logs all of the warnings and issues, so that the user can easily refer to all the errors in a single consistent coherent doc
    if fluctuation_counter > 0:
        fluctuating_gradient = True
    else:
        if counter == number_of_values_to_check:
            converged = True
    return converged

def check_convergence_of_simulation_result_dif(time, result, number_of_values_to_check = 20, tolerance = 0.01):
    """
    Purpose:
        Check if result has converged to a periodic state. result has converged to a periodic state if the percent differences between the last number_of_values_to_check values of result and the last value of result is equal to or smaller than the tolerance.
    Caveat:
        result must be a time-averaged waveform
    Inputs:
        np.array time
            = np.array of time points
        np.array result
            = np.array of result values for which we want to check convergence
            -- must have the same length as time
        float tolerance
            = convergence requires that the last 2 values of result is equal to or smaller than the tolerance
    Returns:
        boolean converged
            = True if result has converged; otherwise, False
    """
    start_index = -1*number_of_values_to_check
    for i in range(start_index, 0):
        if (100.0*np.abs((result[i] - result[-1])/result[-1]) > tolerance):
            return False
    return True

def check_convergence_of_simulation_result_max(time, result, number_of_values_to_check = 20, tolerance = 0.01):
    """
    Purpose:
        Check if result has converged to a periodic state. result has converged to a periodic state if the percent differences between the last number_of_values_to_check values of result and the last value of result is equal to or smaller than the tolerance.
    Caveat:
        result must be a time-averaged waveform
    Inputs:
        np.array time
            = np.array of time points
        np.array result
            = np.array of result values for which we want to check convergence
            -- must have the same length as time
        float tolerance
            = convergence requires that the last 2 values of result is equal to or smaller than the tolerance
    Returns:
        boolean converged
            = True if result has converged; otherwise, False
    """
    max_result = np.amax(result)
    start_index = -1*number_of_values_to_check
    for i in range(start_index, 0):
        if (100.0*np.abs((result[i] - result[-1])/max_result) > tolerance):
            return False
    return True

def compute_time_averaged_result(time, result, parameters):
    """
    Purpose:
        Compute a continuous time-average of result.
    Inputs:
        np.array time
            = np.array of time points
        np.array result
            = np.array of result values
        dict parameters
            -- created from function oneD_to_zeroD_convertor.extract_info_from_solver_input_file
    Returns:
        np.array time_of_time_averaged_result
            = an np.array of the time points for time_averaged_result
            -- this is an np.array of length = len(time) - parameters["number_of_time_pts_per_cardiac_cycle"] + 1
        np.array time_averaged_result
            = an np.array of the continuous time-average of result
            -- the time-averaged is computed via (1/T)*(integral of result over domain of length T), where T = period of a cardiac cycle
            -- this is an np.array of length = len(time) - parameters["number_of_time_pts_per_cardiac_cycle"] + 1
    """
    n = len(time) - parameters["number_of_time_pts_per_cardiac_cycle"]
    time_averaged_result = np.zeros(n + 1)
    for i in range(n, -1, -1):
        if (time[i + parameters["number_of_time_pts_per_cardiac_cycle"] - 1] - time[i] - parameters["cardiac_cycle_period"])/parameters["cardiac_cycle_period"]*100.0 > 0.001:
            print("bounds on numerical integral = ", time[i + parameters["number_of_time_pts_per_cardiac_cycle"] - 1] - time[i])
            print("true period = ", parameters["cardiac_cycle_period"])
            message = "Error. The bounds on the numerical integral do not amount to a single cardiac cycle period."
            raise RuntimeError(message)
        else:
            time_averaged_result[i] = (1.0/parameters["cardiac_cycle_period"])*np.trapz(result[i:i + parameters["number_of_time_pts_per_cardiac_cycle"]], time[i:i + parameters["number_of_time_pts_per_cardiac_cycle"]])
    time_of_time_averaged_result = time[parameters["number_of_time_pts_per_cardiac_cycle"] - 1:]
    return time_of_time_averaged_result, time_averaged_result

def run_convergence_check_of_0d_results(zero_d_results_for_var_names, parameters):# currently here 8/18/20: need to make sure that this function works correctly
    """
    Purpose:
        Check if all 0d simulation results in zero_d_results_for_var_names have converged to a periodic state.
    Inputs:
        dict zero_d_results_for_var_names
            =   {
                    "time" : np.array of simulated time points,

                    "flow" : {var_name : np.array of flow rate,

                    "pressure" : {var_name : np.array of pressure},

                    "wss" : {var_name : np.array of wall shear stress},

                    "internal" : {var_name : np.array of internal block solutions},

                        where var_name is an item in var_name_list (var_name_list generated from run_network_util)
                }
        dict parameters
            -- created from function oneD_to_zeroD_convertor.extract_info_from_solver_input_file
    Returns:
        boolean all_results_converged
            = True if all results in zero_d_results_for_var_names have converged to a periodic state; otherwise, False
    """
    all_results_converged = True
    time = zero_d_results_for_var_names["time"]
    for qoi in list(zero_d_results_for_var_names.keys()):
        if qoi != "time":
            for var_name in list(zero_d_results_for_var_names[qoi].keys()):
                result = zero_d_results_for_var_names[qoi][var_name]
                time_of_time_averaged_result, time_averaged_result = compute_time_averaged_result(time, result, parameters)
                if qoi == "flow" or "wss":

                    # Reference: Xiao N, Alastruey J, Figueroa C. A systematic comparison between 1-D and 3-D hemodynamics in compliant arterial models. International Journal of Numerical Methods in Biomedical Engineering 2014; 30:204–231

                    time_of_time_averaged_result, time_averaged_result = extract_last_cardiac_cycle_simulation_results(time_of_time_averaged_result, time_averaged_result, parameters["number_of_time_pts_per_cardiac_cycle"])

                    converged = check_convergence_of_simulation_result_max(time_of_time_averaged_result, time_averaged_result)

                    # converged = check_convergence_of_simulation_result_gradient(time_of_time_averaged_result, time_averaged_result)
                else:
                    converged = check_convergence_of_simulation_result_dif(time_of_time_averaged_result, time_averaged_result) # changed from using the gradient-based convergence check to the difference-based check on 10/21/20
                if not converged:
                    print("------> Warning: the 0d simulation result for " + var_name + " has not yet converged to a periodic state.\n")
                    all_results_converged = False
    return all_results_converged

def save_simulation_results(zero_d_solver_input_file_path, zero_d_results_for_var_names):
    """
    Purpose:
        Save the 0d simulation results to a .npy file.
    Usage:
        To open and load the .npy file to extract the 0d simulation results, use the following command:
            zero_d_results_for_var_names = np.load(/path/to/zero/d/simulation/results/npy/file, allow_pickle = True).item()
    Inputs:
        string zero_d_solver_input_file_path
            = path to the 0d solver input file
        dict zero_d_results_for_var_names
            =   {
                    "time" : np.array of simulated time points,

                    "flow" : {var_name : np.array of flow rate,

                    "pressure" : {var_name : np.array of pressure},

                    "wss" : {var_name : np.array of wall shear stress},

                    "internal" : {var_name : np.array of internal block solutions},

                        where var_name is an item in var_name_list (var_name_list generated from run_network_util)
                }
    Returns:
        void, but saves a .npy file storing the 0d simulation results (zero_d_results_for_var_names) as a dictionary.
    """
    zero_d_input_file_name, zero_d_input_file_extension = os.path.splitext(zero_d_solver_input_file_path)
    zero_d_simulation_results_file_path = zero_d_input_file_name + "_all_results"
    np.save(zero_d_simulation_results_file_path, zero_d_results_for_var_names)

def comparing_true_to_loaded(zero_d_results_for_var_names_true, zero_d_results_for_var_names_loaded, tolerance = 0.001): # currently here 8/18/20: need to make sure that this function works correctly

    # currently here 02/12/21: use the np.testing.assert_allclose( np_array_1 , np_array_2 ) function to optimize this function

    """
    Purpose:
        This function is only used to compare the 0d simulation results dictionary returned/outputted by set_up_and_run_0d_simulation and the dictionary saved in the .npy file by save_simulation_results. Theoretically, these two dictionaries should be exactly the same, but I use this function just to quickly check and make sure. Otherwise, I don't need this function and I certainly don't need to have it implemented in the final version that gets uploaded to SimVascular.
    Inputs:
        dict zero_d_results_for_var_names_true
            = zero_d_results_for_var_names, but obtained from the output of set_up_and_run_0d_simulation
        dict zero_d_results_for_var_names_loaded
            = zero_d_results_for_var_names, but obtained from the dict stored in the .npy file created by save_simulation_results
        recall that:
            dict zero_d_results_for_var_names
                =   {
                        "time" : np.array of simulated time points,

                        "flow" : {var_name : np.array of flow rate,

                        "pressure" : {var_name : np.array of pressure},

                        "wss" : {var_name : np.array of wall shear stress},

                        "internal" : {var_name : np.array of internal block solutions},

                            where var_name is an item in var_name_list (var_name_list generated from run_network_util)
                    }
        Returns:
            void, but if the two variations of zero_d_results_for_var_names match, then nothing happens. Otherwise; the program prints errors.
    """
    qois_true = list(zero_d_results_for_var_names_true.keys())
    qois_loaded = list(zero_d_results_for_var_names_loaded.keys())
    qois_true.sort()
    qois_loaded.sort()
    if qois_true.sort() == qois_loaded.sort():
        for qoi in qois_true:
            if qoi != "time":
                var_names_true = list(zero_d_results_for_var_names_true[qoi].keys())
                var_names_loaded = list(zero_d_results_for_var_names_loaded[qoi].keys())
                var_names_true.sort()
                var_names_loaded.sort()
                if var_names_true.sort() == var_names_loaded.sort():
                    for var_name in var_names_true:
                        if len(zero_d_results_for_var_names_true[qoi][var_name]) != len(zero_d_results_for_var_names_loaded[qoi][var_name]):
                            print("qoi = ", qoi)
                            print("var_name = ", var_name)
                            message = "Error. The 0d simulation results of the results output of set_up_and_run_0d_simulation do not have the same length as the data saved in the npy file."
                            raise RuntimeError(message)
                        else:
                            for i in range(len(zero_d_results_for_var_names_true[qoi][var_name])):
                                errr = 100.0*np.abs((zero_d_results_for_var_names_true[qoi][var_name][i] - zero_d_results_for_var_names_loaded[qoi][var_name][i])/zero_d_results_for_var_names_true[qoi][var_name][i])
                                if errr > tolerance:
                                    print("i = ", i)
                                    print("qoi = ", qoi)
                                    print("var_name = ", var_name)
                                    print("errr = ", errr)
                                    print("zero_d_results_for_var_names_true[qoi][var_name][i] = ", zero_d_results_for_var_names_true[qoi][var_name][i])
                                    print("zero_d_results_for_var_names_loaded[qoi][var_name][i] = ", zero_d_results_for_var_names_loaded[qoi][var_name][i])
                                    message = "Error. The 0d simulation results of the results output of set_up_and_run_0d_simulation do not match what was saved in the npy file."
                                    raise RuntimeError(message)
                else:
                    message = "Error. The var_names of the results output of set_up_and_run_0d_simulation do not match what was saved in the npy file."
                    raise RuntimeError(message)
            else:
                if len(zero_d_results_for_var_names_true[qoi]) != len(zero_d_results_for_var_names_loaded[qoi]):
                    message = "Error. The 0d simulation times of the results output of set_up_and_run_0d_simulation do not have the same length as the data saved in the npy file."
                    raise RuntimeError(message)
                else:
                    for i in range(len(zero_d_results_for_var_names_true[qoi])):
                        if zero_d_results_for_var_names_true[qoi][i] != zero_d_results_for_var_names_loaded[qoi][i]:
                            print("i = ", i)
                            message = "Error. The 0d simulation times of the results output of set_up_and_run_0d_simulation do not match what was saved in the npy file."
                            raise RuntimeError(message)
    else:
        message = "Error. The qois of the results output of set_up_and_run_0d_simulation do not match what was saved in the npy file."
        raise RuntimeError(message)

def set_up_and_run_0d_simulation(zero_d_solver_input_file_path, draw_directed_graph = False, last_cycle = True, save_0d_simulation_data = True, use_custom_0d_elements = False, custom_0d_elements_arguments_file_path = None, check_convergence = True, use_ICs_from_npy_file = False, ICs_npy_file_path = None, save_y_ydot_to_npy = False, y_ydot_file_path = None, check_jacobian = False, simulation_start_time = 0.0):
    """
    Purpose:
        Create all network_util_NR::LPNBlock objects for the 0d model and run the 0d simulation.
    Inputs:
        string zero_d_solver_input_file_path
            = path to the 0d solver input file
        boolean draw_directed_graph
            = True to visualize the 0d model as a directed graph using networkx -- saves the graph to a .png file (hierarchical graph layout) and a networkx .dot file; False, otherwise. .dot file can be opened with neato from graphviz to visualize the directed in a different format.
        boolean last_cycle
            = True to return 0d simulation results for only the last cycle; False to return the results for all simulated cycles
        boolean save_0d_simulation_data
            = True to save the 0d simulation results (zero_d_results_for_var_names) to a .npy file
        boolean use_custom_0d_elements
            = True to use user-defined, custom 0d elements in the 0d model; False, otherwire
        string custom_0d_elements_arguments_file_path
            = path to user-defined custom 0d element file
        boolean check_convergence
            = True to check if all 0d simulation results have converged to a periodic state; False, otherwise
        float simulation_start_time
            = time at which to begin the 0d simulation
            -- assumes the cardiac cycle begins at t = 0.0 and ends at t = cardiac_cycle_period
            -- 0.0 <= simulation_start_time <= cardiac_cycle_period
    Returns:
        dict zero_d_results_for_var_names
            =   {
                    "time" : np.array of simulated time points,

                    "flow" : {var_name : np.array of flow rate,

                    "pressure" : {var_name : np.array of pressure},

                    "wss" : {var_name : np.array of wall shear stress},

                    "internal" : {var_name : np.array of internal block solutions},

                        where var_name is an item in var_name_list (var_name_list generated from run_network_util)
                }
    """
    if use_custom_0d_elements:
        custom_0d_elements_arguments = import_custom_0d_elements(custom_0d_elements_arguments_file_path)
    else:
        custom_0d_elements_arguments = None
    parameters = oneD_to_zeroD_convertor.extract_info_from_solver_input_file(zero_d_solver_input_file_path)
    parameters["check_jacobian"] = check_jacobian
    create_LPN_blocks(parameters, custom_0d_elements_arguments)
    set_solver_parameters(parameters)
    zero_d_time, results_0d, var_name_list = run_network_util(zero_d_solver_input_file_path, parameters, draw_directed_graph, use_ICs_from_npy_file, ICs_npy_file_path, save_y_ydot_to_npy, y_ydot_file_path, simulation_start_time)
    print("0D simulation completed!\n")
    zero_d_results_for_var_names = reformat_network_util_results(zero_d_time, results_0d, var_name_list)
    if check_convergence:
        all_results_converged = run_convergence_check_of_0d_results(zero_d_results_for_var_names, parameters)
    if last_cycle == True:
        zero_d_results_for_var_names = run_last_cycle_extraction_routines(parameters["cardiac_cycle_period"], parameters["number_of_time_pts_per_cardiac_cycle"], zero_d_results_for_var_names)
    if save_0d_simulation_data:
        save_simulation_results(zero_d_solver_input_file_path, zero_d_results_for_var_names)
    return zero_d_results_for_var_names


# currently here 10/13/20:
# todos:
# 1. give users the option to
    # 1) specify the number of cycles to run for
    # 2) dont specify the number of cycles to run for, but instead run the 0d model until convergence
# 2. check convergence function to use difference in mean waveform value (mean flow and mean pressure) to check for periodic convergence rather than using the gradient method to check
#
# initialization options:
#     0. use zero initial conditions
#     1. simulate a resistor only model (including BCs) and use the converged solution as the initial condition for the original model
#     2. simulate the original model but with the mean flow value as the inflow BC and use a coarse time step size and then use the solution as the initial condition for the original model
#   3. simulate the original model but with a super coarse time step first and then use that solution as the initialization for a original model with a finer time step size
