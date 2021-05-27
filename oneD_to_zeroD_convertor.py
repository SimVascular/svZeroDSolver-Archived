# Copyright (c) Stanford University, The Regents of the University of
#               California, and others.
#
# All Rights Reserved.
#
# See Copyright-SimVascular.txt for additional details.
#
# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject
# to the following conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
# IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
# TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
# PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
# OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
# PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
# LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
# NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
# SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

"""
Author: Pham, Jonathan

This code converts a 1D model into a 0D model by reading in the 1D solver input file, extracting relevant geometry and material properties, and creating the resulting 0D solver input file. The 0D model can then be simulated with the 0D solver via the 0D solver linker code, run_0d_solver.py.

Available vessel models:
    1) R (steady resistor)
    2) C (steady capacitor)
    3) L (steady inductor)
    4) RC (steady resistor-capacitor)
    5) RL (steady resistor-inductor)
    6) RCL (steady resistor-capacitor-inductor)
    7) custom, user-defined vessel model

Available junction types:
    1) NORMAL_JUNCTION (mass conservation only)
    2) custom, user-defined junction type

Available material types:
    1) LINEAR
    2) OLUFSEN
"""

import sys
import os
import numpy as np

""" TODOs """
# need to test this code and make sure that every single variable output from each function yields the correct value

# currently here 7/21/20: need to create a single coherent document that explains how to use the entire 0d code, from the 1d to 0d convertor and the run_0d_solver and network_util_NR. and how to use the built in vessel model types (R, RC, RL, L, C, RCL), built junction types (NORMAL_JUNCTION), built in outlet bc types (RESISTANCE, RCR, FLOW, PRESSURE), and the built in inlet bc types (FLOW, PRESSURE). and how to use custom user-defined functions (requires the creation of the sample_special_user_defined_elements_v2.py file, which involves creation of a dictionary of arguments, and the creation of the special element classes in network_util_NR, and the specification of those special element types in the 0d solver input file -- note that for special vessel elements, element values should not be listed in the ELEMENT card like they are for the built in types, like the R, RL, RC, etc... and for the special BC (outlet or inlet), need to specify the name of the BC type (using the same name as the one used to define the special class in network_util) and the indication of "NONE" for the datatable name in the ELEMENT CARD and the INLET_BOUNDARY_CONDITION CARD.

# if i want to make this code better, i should make classes/objects for each CARD of the 1d/0d solver input file. Ie rather than having datatable_types and datatable_values are two separate entities, i make one class called datatable and then make an instance of that class, where the fields of that class will be the datatable_type and datatable_values

def create_Ehor_lambda_function(k1, k2, k3, material_type):
    """
    Purpose:
        Create a lambda function to compute E*h/r as a function of r. This value is used to compute vessel capacitance.
        Sources:    Alison Marsden's CME 285 Spring 2019 Course Reader, pg 68
                        http://simvascular.github.io/docs1DSimulation.html#format_material
                    http://simvascular.github.io/docs1DSimulation.html
    Available material models:
        LINEAR:
            Eh/r = k1
        OLUFSEN:
            Eh/r = k1*exp(k2*r) + k3
    Inputs:
        float k1
        float k2
        float k3
        string material_type
            = "LINEAR" or "OLUFSEN"
    Returns:
        lambda function
    """
    if material_type == "OLUFSEN":
        return lambda r: k1*np.exp(k2*r) + k3
    elif material_type == "LINEAR":
        return lambda r: k1
    else:
        message = "Error. Material type, " + material_type + ", has not been implemented. Only implemented materials are OLUSEN and LINEAR."
        raise ValueError(message)

def compute_average_segment_radius(segment_inlet_area, segment_outlet_area):
    """
    Purpose:
        Compute the average radius for a given vessel segment
    Inputs:
        float segment_inlet_area
            = cross-sectional area of a segment's inlet face
        float segment_outlet_area
            = cross-sectional area of a segment's outlet face
    Returns:
        float average_vessel_radius
            = (segment inlet radius + segment outlet radius)/2.0
    """
    return (np.sqrt(segment_inlet_area/np.pi) + np.sqrt(segment_outlet_area/np.pi))/2.0

def calc_segment_tangent_vector(segment_inlet_node_coords, segment_outlet_node_coords):
    """
    Purpose:
        Compute the (un-normalized) tanget vector for the vessel segment, where the vector points from the inlet node (placed at the origin) to the outlet node
    Inputs:
        list segment_inlet_node_coords
            = [x, y, z] coordinates of a segment's inlet node
        list segment_outlet_node_coords
            = [x, y, z] coordinates of a segment's outlet node
    Returns:
        list vessel_vector
            = tangent vector, in the form of [x, y, z] vector components, for the vessel segment
    """
    outlet_nodes_coords_shifted = np.asarray(segment_outlet_node_coords) - np.asarray(segment_inlet_node_coords)
    vessel_vector = outlet_nodes_coords_shifted.tolist()
    return vessel_vector

def extract_info_from_solver_input_file(solver_input_file_path, one_d_inlet_segment_number = 0):
    """
    Purpose:
        Extract geometric and material properties prescribed in the 1d or 0d solver input file
    Inputs:
        string solver_input_file_path
            = name of the 1d or 0d solver input file
        int one_d_inlet_segment_number
            = inlet segment number of the 1d model (this number is usually zero)
    Returns:
        dict parameters = # see comments below for description of each item in parameters
            {
                "model_name" : model_name,
                "node_coordinates" : node_coordinates,
                "nodes_of_segments" : nodes_of_segments,
                "segment_names" : segment_names,
                "segment_numbers_list" : segment_numbers_list,
                "radii" : radii,
                "lengths" : lengths,
                "segment_tangent_vectors" : segment_tangent_vectors,
                "material_names" : material_names,
                "material_types" : material_types,
                "rho" : rho,
                "mu" : mu,
                "Ehor_func" : Ehor_func,
                "Ehor_values" : Ehor_values,
                "segment_0d_types" : segment_0d_types,
                "segment_0d_values" : segment_0d_values,
                "joint_node_numbers" : joint_node_numbers,
                "joint_inlet_names" : joint_inlet_names,
                "joint_outlet_names" : joint_outlet_names,
                "joint_inlet_segments" : joint_inlet_segments,
                "joint_outlet_segments" : joint_outlet_segments,
                "joint_name_to_segments_map" : joint_name_to_segments_map,
                "junction_types" : junction_types,
                "boundary_condition_types" : boundary_condition_types,
                "boundary_condition_datatable_names" : boundary_condition_datatable_names,
                "datatable_types" : datatable_types,
                "datatable_values" : datatable_values,
                "outlet_segments_of_model" : outlet_segments_of_model,
                "inlet_segments_of_model" : inlet_segments_of_model,
                "number_of_time_pts_per_cardiac_cycle" : number_of_time_pts_per_cardiac_cycle,
                "number_of_cardiac_cycles" : number_of_cardiac_cycles
            }
    """

    model_name = "" # name prescribed in the MODEL card of the 1d/0d solver input file
    node_coordinates = {} # {node number : [x-coordinate, y-coord, z-coord]}
    nodes_of_segments = {} # {segment number : [inlet_node_number, outlet_node_number]}
    segment_names = {} # {segment number : segment name}
    radii = {} # {segment number : average segment radius}
    lengths = {} # {segment number : segment length}
    segment_tangent_vectors = {} # {segment number : segment tangent vector in the form of [x, y, z] vector components}
        # the segment tanget vector is unnormalized and has been translated so that its inlet node is located at the origin
    material_names = {} # {segment number : material card name}
    material_types = {} # {segment number : material type, i.e. OLUFSEN OR LINEAR}
    rho = {} # {segment number : density of blood}
    mu = {} # {segment number : dynamic viscosity of blood}
    Ehor_func = {} # {segment number : lambda function for Eh/r as a function of r}
    Ehor_values = {} # {segment number : Eh/r value}
    segment_0d_types = {} # {segment number : 0d element type, i.e. R, RC, RCL, RL, L, C}
    segment_0d_values = {} # {segment number : [list of 0d element values]}
        # example: if the 0d element type is R, then the 0d element value is [R_value]
        # example: if the 0d element type is RCL, then the 0d element values are [R_value, C_value, L_value]
    joint_node_numbers = {} # {joint name : joint node number}
    joint_inlet_names = {} # {joint name : joint inlet name}
    joint_outlet_names = {} # {joint name : joint outlet name}
    joint_inlet_segments = {} # {joint inlet name : [list of inlet segments at the joint]}
    joint_outlet_segments = {} # {joint outlet name : [list of outlet segments at the joint]}
    joint_name_to_segments_map = {} # { joint name : { "inlet" : [list of inlet segments at the joint], "outlet" : [list of outlet segments at the joint] } }
    junction_types = {} # {joint name : junction type, i.e. NORMAL_JUNCTION}
    boundary_condition_types = {"inlet" : {}, "outlet" : {}} # {"inlet" : {segment number : type of boundary condition, i.e. "FLOW", "PRESSURE"}, "outlet" : {segment number : type of boundary condition, i.e. "RESISTANCE", "RCR", "CORONARY", "FLOW", "PRESSURE"}}
    boundary_condition_datatable_names = {"inlet" : {}, "outlet" : {}} # {"inlet" : {segment number : datatable name for boundary condition}, "outlet" : {segment number : corresponding datatable name for boundary condition}}
    datatable_types = {} # {datatable name for boundary condition : type of datatable, i.e LIST}
    datatable_values = {} # {datatable name for boundary condition : [time1 value1 time2 value2 ...]}
        # value corresponds to the value of the BC, ie resistance value
    inlet_segments_of_model = [] # [model's inlet segments (these are attached to the inlet boundary conditions)]
    outlet_segments_of_model = [] # [model's outlet segments (these are attached to the outlet boundary conditions)]

    with open(solver_input_file_path) as f:
        while True:
            line = f.readline()
            if line == "": # signifies end of file
                break
            else:
                if line.startswith("MODEL"):
                    model_name = line.split()[1]
                elif line.startswith("NODE"):
                    line_list = line.split()
                    node_coordinates[int(line_list[1])] = [float(line_list[2]), float(line_list[3]), float(line_list[4])]
                elif line.startswith("JOINT "):
                    line_list = line.split()
                    joint_name = line_list[1]
                    if (not joint_name.startswith("J")) or (not joint_name[1].isnumeric()):
                        message = "Error in file, " + solver_input_file_path + ". Joint name, " + joint_name + ", is not 'J' followed by numeric values. The 0D solver assumes that all joint names are 'J' followed by numbers in the solver input file."
                        raise ValueError(message)
                    joint_node_numbers[joint_name] = int(line_list[2])
                    joint_inlet_names[joint_name] = line_list[3]
                    joint_outlet_names[joint_name] = line_list[4]
                elif line.startswith("JOINTINLET"):
                    line_list = line.split()
                    joint_in_name = line_list[1]
                    if (not joint_in_name.startswith("IN")) or (not joint_in_name[2].isnumeric()):
                        message = "Error in file, " + solver_input_file_path + ". Joint inlet name, " + joint_in_name + ", is not 'IN' followed by numbers. This code assumes that all joint inlet names is 'IN' followed by numbers in the solver input file."
                        raise ValueError(message)
                    joint_inlet_segments[joint_in_name] = [int(i) for i in line_list[3:]] # inlet segment #'s at joint
                    if len(joint_inlet_segments[joint_in_name]) != int(line_list[2]):
                        message = "Error in file, " + solver_input_file_path + ". Number of inlet segments for joint inlet, " + joint_in_name + ", does not match the number of prescribed inlet segments."
                        raise ValueError(message)
                elif line.startswith("JOINTOUTLET"):
                    line_list = line.split()
                    joint_out_name = line_list[1]
                    if (not joint_out_name.startswith("OUT")) or (not joint_out_name[3].isnumeric()):
                        message = "Error in file, " + solver_input_file_path + ". Joint outlet name, " + joint_out_name + ", is not 'OUT' followed by numbers. This code assumes that all joint outlet names is 'OUT' followed by numbers in the solver input file."
                        raise ValueError(message)
                    joint_outlet_segments[joint_out_name] = [int(i) for i in line_list[3:]] # outlet segment #'s at joint
                    if len(joint_outlet_segments[joint_out_name]) != int(line_list[2]):
                        message = "Error in file, " + solver_input_file_path + ". Number of outlet segments for joint outlet, " + joint_out_name + ", does not match the number of prescribed outlet segments."
                        raise ValueError(message)
                elif line.startswith("JUNCTION_MODEL"): # this section is only available in the 0d solver input file
                    line_list = line.split()
                    joint_name = line_list[1]
                    if (not joint_name.startswith("J")) or (not joint_name[1].isnumeric()):
                        message = "Error in file, " + solver_input_file_path + ". Joint name, " + joint_name + ", is not 'J' followed by numeric values. The 0D solver assumes that all joint names are 'J' followed by numbers in the 0d solver input file. Note that the joint names are the same as the junction names."
                        raise ValueError(message)
                    junction_types[joint_name] = line_list[2]
                elif line.startswith("SEGMENT"): # this section is only available in the 1d solver input file
                    line_list = line.split()
                    segment_number = int(line_list[2])
                    segment_names[segment_number] = line_list[1]
                    lengths[segment_number] = float(line_list[3])
                    nodes_of_segments[segment_number] = [int(line_list[5]), int(line_list[6])]
                    radii[segment_number] = compute_average_segment_radius(float(line_list[7]), float(line_list[8]))
                    boundary_condition_types["outlet"][segment_number] = line_list[15]
                    boundary_condition_datatable_names["outlet"][segment_number] = line_list[16]
                    if line_list[15] != "NOBOUND":
                        outlet_segments_of_model.append(segment_number)
                    material_names[segment_number] = line_list[10]
                    segment_tangent_vectors[segment_number] = calc_segment_tangent_vector(node_coordinates[int((line_list[5]))], node_coordinates[int((line_list[6]))])
                elif line.startswith("ELEMENT"): # this section is only available in the 0d solver input file
                    line_list = line.split()
                    segment_number = int(line_list[2])
                    segment_names[segment_number] = line_list[1]
                    radii[segment_number] = float(line_list[3])
                    lengths[segment_number] = float(line_list[4])
                    rho[segment_number] = float(line_list[5])
                    mu[segment_number] = float(line_list[6])
                    Ehor_values[segment_number] = float(line_list[7])
                    nodes_of_segments[segment_number] = [int(line_list[8]), int(line_list[9])]
                    segment_tangent_vectors[segment_number] = [float(line_list[10]), float(line_list[11]), float(line_list[12])]
                    boundary_condition_types["outlet"][segment_number] = line_list[13]
                    boundary_condition_datatable_names["outlet"][segment_number] = line_list[14]
                    if line_list[13] != "NOBOUND":
                        outlet_segments_of_model.append(segment_number)
                    segment_0d_types[segment_number] = line_list[15]
                    segment_0d_values[segment_number] = [float(i) for i in line_list[16:]]
                elif line.startswith("DATATABLE"):
                    line_list = line.split()
                    datatable_name = line_list[1]
                    datatable_types[datatable_name] = line_list[2]
                    if line_list[2] == "LIST":
                        datatable_values[datatable_name] = []
                        line = f.readline()
                        line_list = line.split()
                        while line_list[0] != "ENDDATATABLE":
                            datatable_values[datatable_name] = datatable_values[datatable_name] + [float(i) for i in line_list]
                            line = f.readline()
                            line_list = line.split()
                    else:
                        message = "Error in file, " + solver_input_file_path + ". DATATABLE type, " + line_list[2] + ", is not allowed. Only allowed type is 'LIST'."
                        raise ValueError(message)
                elif line.startswith("INLET_BOUNDARY_CONDITION"): # this section is only available in the 0d solver input file
                    line_list = line.split()
                    segment_number = int(line_list[1])
                    inlet_segments_of_model.append(segment_number)
                    boundary_condition_types["inlet"][segment_number] = line_list[2]
                    boundary_condition_datatable_names["inlet"][segment_number] = line_list[3]
                elif line.startswith("SOLVEROPTIONS"):
                    line_list = line.split()
                    if line_list[0] == "SOLVEROPTIONS": # 1D solver options
                        segment_number = one_d_inlet_segment_number # 1d models have only a single inlet segment number
                        inlet_segments_of_model.append(segment_number)
                        boundary_condition_types["inlet"][segment_number] = line_list[6]
                        boundary_condition_datatable_names["inlet"][segment_number] = line_list[5]
                    elif line_list[0] == "SOLVEROPTIONS_0D": # 0D solver options
                        number_of_time_pts_per_cardiac_cycle = int(line_list[1])
                        number_of_cardiac_cycles = int(line_list[2])
                    else:
                        message = "Error. Unidentified solver card, " + line_list[0] + "."
                        raise RuntimeError(message)
                elif line.startswith("MATERIAL"): # this section is only available in the 1d solver input file
                    line_list = line.split()
                    material_name = line_list[1]
                    segments_with_this_material_name = [k for k, v in material_names.items() if v == material_name]
                    rho.update({segment_number : float(line_list[3]) for segment_number in segments_with_this_material_name})
                    mu.update({segment_number : float(line_list[4]) for segment_number in segments_with_this_material_name})
                    material_type = line_list[2]
                    material_types.update({segment_number : material_type for segment_number in segments_with_this_material_name})
                    k1 = float(line_list[7])
                    k2 = float(line_list[8])
                    k3 = float(line_list[9])
                    Ehor_func.update({segment_number : create_Ehor_lambda_function(k1, k2, k3, material_type) for segment_number in segments_with_this_material_name})

    segment_numbers_list = list(segment_names.keys())
    segment_numbers_list.sort()
    if len(segment_numbers_list) == 1:
        if segment_numbers_list != inlet_segments_of_model:
            message = "Error in file, " + solver_input_file_path + ". This model has only 1 segment, segment #" + str(segment_numbers_list[0]) + ", but this segment isn't equal to inlet_segments_of_model[0], " + str(inlet_segments_of_model[0]) + "."
            raise RuntimeError(message)
        elif segment_numbers_list != outlet_segments_of_model:
            message = "Error in file, " + solver_input_file_path + ". This model has only 1 segment, segment #" + str(segment_numbers_list[0]) + ", but this segment isn't equal to outlet_segments_of_model[0], " + str(outlet_segments_of_model[0]) + "."
            raise RuntimeError(message)

    for joint_name in list(joint_node_numbers.keys()):
        joint_name_to_segments_map[joint_name] = {"inlet" : joint_inlet_segments[joint_inlet_names[joint_name]].copy(), "outlet" : joint_outlet_segments[joint_outlet_names[joint_name]].copy()}

    parameters = ({
        "model_name" : model_name,
        "node_coordinates" : node_coordinates,
        "nodes_of_segments" : nodes_of_segments,
        "segment_names" : segment_names,
        "segment_numbers_list" : segment_numbers_list,
        "radii" : radii,
        "lengths" : lengths,
        "segment_tangent_vectors" : segment_tangent_vectors,
        "material_names" : material_names,
        "material_types" : material_types,
        "rho" : rho,
        "mu" : mu,
        "Ehor_func" : Ehor_func,
        "Ehor_values" : Ehor_values,
        "segment_0d_types" : segment_0d_types,
        "segment_0d_values" : segment_0d_values,
        "joint_node_numbers" : joint_node_numbers,
        "joint_inlet_names" : joint_inlet_names,
        "joint_outlet_names" : joint_outlet_names,
        "joint_inlet_segments" : joint_inlet_segments,
        "joint_outlet_segments" : joint_outlet_segments,
        "joint_name_to_segments_map" : joint_name_to_segments_map,
        "junction_types" : junction_types,
        "boundary_condition_types" : boundary_condition_types,
        "boundary_condition_datatable_names" : boundary_condition_datatable_names,
        "datatable_types" : datatable_types,
        "datatable_values" : datatable_values,
        "outlet_segments_of_model" : outlet_segments_of_model,
        "inlet_segments_of_model" : inlet_segments_of_model,
        "number_of_time_pts_per_cardiac_cycle" : number_of_time_pts_per_cardiac_cycle,
        "number_of_cardiac_cycles" : number_of_cardiac_cycles
        })
    return parameters

def create_segment_0d_types(parameters, default_element_type, segment_0d_types_input_file_path):
    """
    Purpose:
        Create a dictionary of 0d element types for each segment and create/update the txt file, segment_0d_types_input_file_path, that tells the users these 0d element types.
    Inputs:
        dict parameters
            -- created from function extract_info_from_solver_input_file
        string default_element_type
            = default 0d element type, i.e. "R", "RC", "RCL", "RL", "L", "C", etc. to be used for the majority of the segments
        string segment_0d_types_input_file_path
            = path to the file that tells users the 0d element types for each segment
    Returns:
        void, but
            1) updates parameters to include/update:
                    dict segment_0d_types
                        = {segment number : 0d element type, i.e. R, RC, RCL}
            2) creates/updates the txt file, segment_0d_types_input_file_path, that tells the users these 0d element types
    """
    segment_0d_types = {}
    segment_0d_types_input_file = open(segment_0d_types_input_file_path, "w")
    segment_0d_types_input_file.write("# segment_name : segment_number : 0d_element_type" "\n")
    for segment_number in parameters["segment_numbers_list"]:
        segment_0d_types_input_file.write(parameters["segment_names"][segment_number] + " : " + str(segment_number) + " : " + default_element_type + "\n")
        segment_0d_types[segment_number] = default_element_type
    segment_0d_types_input_file.close()
    parameters.update({"segment_0d_types" : segment_0d_types})

def read_segment_0d_types_input_file(parameters, segment_0d_types_input_file_path):
    """
    Purpose:
        Extract the 0d element type used for each segment from the txt file, segment_0d_types_input_file_path
    Inputs:
        dict parameters
            -- created from function extract_info_from_solver_input_file
        string segment_0d_types_input_file_path
            = path to the file that tells users the 0d element types for each segment
    Returns:
        void, but updates parameters to include:
            dict segment_0d_types
                = {segment number : 0d element type, i.e. R, RC, RCL}
    """
    segment_0d_types = {}
    with open(segment_0d_types_input_file_path) as f:
        for line in f:
            if not line.startswith("#"):
                line_list = line.replace(":", " ").split()
                segment_0d_types[int(line_list[1])] = line_list[2]
    segment_numbers_list_temp = list(segment_0d_types.keys())
    segment_numbers_list_temp.sort()
    if parameters["segment_numbers_list"] != segment_numbers_list_temp:
        print("segments in 1d solver input file   : ", parameters["segment_numbers_list"])
        print("segments in segment_input.txt file : ", segment_numbers_list_temp)
        message = "Error in file, " + segment_0d_types_input_file_path + ". Segments listed in the file do not match the segments listed in the SEGMENT card of the parent 1D solver input file."
        raise RuntimeError(message)
    parameters.update({"segment_0d_types" : segment_0d_types})

def check_default_0d_element_type_with_user(parameters, default_element_type, segment_0d_types_input_file_path, check_default_element_type):
    """
    Purpose:
        Run function, create_segment_0d_types, and give the user the option to manually change the 0d element types used for each segment
    Inputs:
        dict parameters
            -- created from function extract_info_from_solver_input_file
        string default_element_type
            = default 0d element type, i.e. "R", "RC", "RCL", "RL", "L", "C", etc. to be used for the majority of the segments
        string segment_0d_types_input_file_path
            = path to the file that tells users the 0d element types for each segment
        boolen check_default_element_type
            = True to give user ability to manually change the 0d element types used for each segment; False otherwise
    Returns:
        void, but updates parameters to include/update:
            dict segment_0d_types
                = {segment number : 0d element type, i.e. R, RC, RCL}
    """
    create_segment_0d_types(parameters, default_element_type, segment_0d_types_input_file_path)
    if check_default_element_type == True:
        command_line_string = "All 1d segments/vessels will be modeled as " + default_element_type + " in the 0d surrogate model by default. Would you like to model your vessels with another type of element? The built-in 0d element options are R, C, L, RC, RL, and RCL. You are also allowed to specify a custom element type. If yes, first make the change in the txt file for the desired vessels/segments and then type 'y' in the command line; otherwise, type 'n' [y/n]: "
        change_segment_0d_types = input(command_line_string)
        while True:
            if change_segment_0d_types == "y" or change_segment_0d_types == "Y":
                read_segment_0d_types_input_file(parameters, segment_0d_types_input_file_path)
                break
            elif change_segment_0d_types == "n" or change_segment_0d_types == "N":
                break
            else:
                change_segment_0d_types = input("Error. Not a valid input. [y/n]: ")

def create_junction_types(parameters, default_junction_type, junction_types_input_file_path):
    """
    Purpose:
        Create a dictionary of junctions types for each joint and creates/updates the txt file, junction_types_input_file_path, that tells the users these junction types.
    Inputs:
        dict parameters
            -- created from function extract_info_from_solver_input_file
        string default_junction_type
            = default junction type, i.e. "NORMAL_JUNCTION", etc. to be used for the majority of the joints
        string junction_types_input_file_path
            = path to the file that tells users the junction types for each joint
    Returns:
        void, but
            1) updates parameters to include/update:
                dict junction_types
                    = {joint_name : junction type, i.e. NORMAL_JUNCTION}
            2)  creates/updates the txt file, junction_types_input_file_path, that tells the users these junction types
    """
    junction_types = {}
    junction_types_input_file = open(junction_types_input_file_path, "w")
    junction_types_input_file.write("# joint_name : junction_type" "\n")
    joint_names = list(parameters["joint_name_to_segments_map"].keys())
    joint_names.sort()
    for joint_name in joint_names:
        junction_types_input_file.write(joint_name + " : " + default_junction_type + "\n")
        junction_types[joint_name] = default_junction_type
    junction_types_input_file.close()
    parameters.update({"junction_types" : junction_types})

def read_junction_types_input_file(parameters, junction_types_input_file_path):
    """
    Purpose:
        Extract the junction type used for each joint from the txt file, junction_types_input_file_path
    Inputs:
        dict parameters
            -- created from function extract_info_from_solver_input_file
        string junction_types_input_file_path
            = path to the file that tells users the junction types for each joint
    Returns:
        void, but updates parameters to include:
            dict junction_types
                = {joint_name : junction type, i.e. NORMAL_JUNCTION}
    """
    junction_types = {}
    with open(junction_types_input_file_path) as f:
        for line in f:
            if not line.startswith("#"):
                line_list = line.replace(":", " ").split()
                joint_name = line_list[0]
                if (not joint_name.startswith("J")) or (not joint_name[1].isnumeric()):
                    message = "Error in file, " + junction_types_input_file_path + ". Joint name, " + joint_name + ", is not 'J' followed by numeric values. The 0D solver assumes that all joint names are 'J' followed by numbers. Please change the joint name to follow this criteria."
                    raise ValueError(message)
                junction_types[joint_name] = line_list[1]
    junction_names_temp1 = list(junction_types.keys())
    junction_names_temp1.sort()
    junction_names_temp2 = list(parameters["joint_name_to_segments_map"].keys())
    junction_names_temp2.sort()
    if junction_names_temp1 != junction_names_temp2:
        print("Joints in the 1d solver input file          : ", junction_names_temp2)
        print("Joints in the junction_types_input.txt file : ", junction_names_temp1)
        message = "Error in file, " + junction_types_input_file_path + ". Joints listed in the file do not match the joints listed in the JOINT card of the 1D solver input file."
        raise RuntimeError(message)
    parameters.update({"junction_types" : junction_types})

def check_default_junction_type_with_user(parameters, default_junction_type, junction_types_input_file_path, check_default_junction_type):
    """
    Purpose:
        Run create_junction_types and give the user the option to manually change the junction types used for each joint
    Inputs:
        dict parameters
            -- created from function extract_info_from_solver_input_file
        string default_junction_type
            = default junction type, i.e. "NORMAL_JUNCTION", etc. to be used for the majority of the joints
        string junction_types_input_file_path
            = path to the file that tells users the junction types for each joint
        boolean check_default_junction_type
            = True to give the user the ability to manually change the junction type used for each joint; False otherwise
    Returns:
        void, but updates parameters to include/update:
            dict junction_types
                = {joint_name : junction type, i.e. NORMAL_JUNCTION}
    """
    create_junction_types(parameters, default_junction_type, junction_types_input_file_path)
    if check_default_junction_type == True:
        command_line_string = "All junctions are modeled as " + default_junction_type + " by default. Would you like to model your junctions with another type? If yes, first make the change in the junction_types_input.txt file for the desired junctions and then type 'y' in the command line. [y/n]: "
        change_junction_types = input(command_line_string)
        while True:
            if change_junction_types == "y" or change_junction_types == "Y":
                read_junction_types_input_file(parameters, junction_types_input_file_path)
                break
            elif change_junction_types == "n" or change_junction_types == "N":
                break
            else:
                change_junction_types = input("Error. Not a valid input. [y/n]: ")

def compute_resistance(parameters, segment_number):
    """
    Purpose:
        Compute resistance for segment #segment_number using Poiseuille flow
    Inputs:
        dict parameters
            -- created from function extract_info_from_solver_input_file
        int segment_number
            = segment number for which we want to compute the resistance
    Returns:
        float R
            = resistance value = 8*mu*L/(pi*r^4)
    """
    return 8.0*parameters["mu"][segment_number]*parameters["lengths"][segment_number]/(np.pi*(parameters["radii"][segment_number]**4))

def compute_inductance(parameters, segment_number):
    """
    Purpose:
        Compute inductance for segment #segment_number using Poiseuille flow
    Inputs:
        dict parameters
            -- created from function extract_info_from_solver_input_file
        int segment_number
            = segment number for which we want to compute the inductance
    Returns:
        float L
            = inductance value = L*rho/(pi*r^2)
    """
    return parameters["lengths"][segment_number]*parameters["rho"][segment_number]/(np.pi*(parameters["radii"][segment_number]**2))

def compute_capacitance(parameters, segment_number):
    """
    Purpose:
        Compute capacitance for segment #segment_number using Poiseuille flow
    Inputs:
        dict parameters
            -- created from function extract_info_from_solver_input_file
        int segment_number
            = segment number for which we want to compute the capacitance
    Returns:
        float C
            = capacitance value = 3*L*pi*r^3/(2*E*h)
    """
    return 3.0*parameters["lengths"][segment_number]*np.pi*(parameters["radii"][segment_number]**2)/(2*parameters["Ehor_func"][segment_number](parameters["radii"][segment_number]))

def create_0d_elements(parameters):
    """
    Purpose:
        Compute the 0d element value (from Poiseuille flow) for the segments being modeled with one of the built-in 0d element types (R, C, L, RC, RL, RCL)
    Inputs:
        dict parameters
            -- created from function extract_info_from_solver_input_file
    Returns:
        void, but updates parameters to include:
            dict segment_0d_values
                = {segment_number : [list of 0d element values]}
    """
    segment_0d_values = {}
    for segment_number in parameters["segment_numbers_list"]:
        values_list = []
        temp_type = parameters["segment_0d_types"][segment_number]
        if temp_type == "R":
            values_list.append(compute_resistance(parameters, segment_number))
        elif temp_type == "RC":
            values_list.append(compute_resistance(parameters, segment_number))
            values_list.append(compute_capacitance(parameters, segment_number))
        elif temp_type == "RL":
            values_list.append(compute_resistance(parameters, segment_number))
            values_list.append(compute_inductance(parameters, segment_number))
        elif temp_type == "RCL":
            values_list.append(compute_resistance(parameters, segment_number))
            values_list.append(compute_capacitance(parameters, segment_number))
            values_list.append(compute_inductance(parameters, segment_number))
        elif temp_type == "L":
            values_list.append(compute_inductance(parameters, segment_number))
        elif temp_type == "C":
            values_list.append(compute_capacitance(parameters, segment_number))
        segment_0d_values[segment_number] = values_list
    parameters.update({"segment_0d_values" : segment_0d_values})

def get_zero_solver_input_file_path(one_d_solver_input_file_path, destination_folder, use_self_defined_name):
    """
    Purpose:
        Get the name of the 0d solver input file, its associated "segment_0d_types_input_file", and its associated "junction_types_input_file"
    Inputs:
        string one_d_solver_input_file_path
            = path to the 1d solver input file
        string destination_folder
            = path to directory in which you want to save the 0d solver input file and its auxiliary files
        boolean use_self_defined_name
            = True to allow the user to specify the name of the 0d solver input file; False to automatically make the 0d solver input file name the same as the 1d solver input file name (but with an extra '_0d' appended)
    Returns:
        string zero_d_solver_input_file_path
            = path to the 0d solver input file
        string segment_0d_types_input_file_path
            = path to the file that tells users the 0d element types for each segment
        string junction_types_input_file_path
            = path the the file that tells users the junction types for each joint
    """
    one_d_input_file_name = os.path.basename(one_d_solver_input_file_path)
    file_name, file_extension = os.path.splitext(one_d_input_file_name)
    if use_self_defined_name == True:
        zero_d_input_file_name = input("Enter a name for your 0d solver input file (do not include the file extension): ")
    else:
        zero_d_input_file_name = file_name + "_0d"
    segment_0d_types_input_file_name = zero_d_input_file_name + "_segment_0d_types_input" + ".txt"
    junction_types_input_file_name = zero_d_input_file_name + "_junction_types_input" + ".txt"
    zero_d_input_file_name = zero_d_input_file_name + file_extension
    zero_d_solver_input_file_path = os.path.join(destination_folder, zero_d_input_file_name)
    segment_0d_types_input_file_path = os.path.join(destination_folder, segment_0d_types_input_file_name)
    junction_types_input_file_path = os.path.join(destination_folder, junction_types_input_file_name)
    return zero_d_solver_input_file_path, segment_0d_types_input_file_path, junction_types_input_file_path

def create_0d_solver_input_file(parameters, zero_d_solver_input_file_path):
    """
    Purpose:
        Create the 0d solver input file
    Inputs:
        dict parameters
        string zero_d_solver_input_file_path
            = path to the 0d solver input file
    Returns:
        void but writes out the 0d solver input file
    """
    print("Creating 0d solver input file...")

    zero_d_input_file = open(zero_d_solver_input_file_path, "w")
    zero_d_input_file.write("# ================================\n# " + parameters["model_name"] + " MODEL" + "\n# ================================\n\n# ==========\n# MODEL CARD\n# ==========\n# - model_name\n\nMODEL " + parameters["model_name"] + "\n")

    # write nodes
    zero_d_input_file.write("\n### DO NOT CHANGE THIS SECTION - generated automatically\n#\n# ==========\n# NODE CARD\n# ==========\n# - node_number\n# - node_x_coord\n# - node_y_coord\n# - node_z_coord\n\n")
    for node_number in list(parameters["node_coordinates"].keys()):
        zero_d_input_file.write("NODE " + str(node_number) + ''.join(" " + str(i) for i in parameters["node_coordinates"][node_number]) + "\n")

    # write joints
    zero_d_input_file.write("\n### DO NOT CHANGE THIS SECTION - generated automatically\n#\n# ==========\n# JOINT CARD\n# ==========\n# - joint_name\n# - joint_node\n# - joint_inlet_name\n# - joint_outlet_name\n")
    zero_d_input_file.write("\n### DO NOT CHANGE THIS SECTION - generated automatically\n#\n# ================================\n# JOINTINLET AND JOINTOUTLET CARDS\n# ================================\n# - joint_inlet/outlet_name\n# - total_number_of_joint_inlet/outlet_segments\n# - list_of_joint_inlet/outlet_segments\n\n")
    list_of_joint_names = list(parameters["joint_node_numbers"].keys())
    list_of_joint_names.sort()
    for joint_name in list_of_joint_names:
        zero_d_input_file.write("JOINT " + joint_name + " " + str(parameters["joint_node_numbers"][joint_name]) + " " + parameters["joint_inlet_names"][joint_name] + " " + parameters["joint_outlet_names"][joint_name] + "\n")
        zero_d_input_file.write("JOINTINLET " + parameters["joint_inlet_names"][joint_name] + " " + str(len(parameters["joint_inlet_segments"][parameters["joint_inlet_names"][joint_name]])) + ''.join(" " + str(i) for i in parameters["joint_inlet_segments"][parameters["joint_inlet_names"][joint_name]]) + "\n")
        zero_d_input_file.write("JOINTOUTLET " + parameters["joint_outlet_names"][joint_name] + " " + str(len(parameters["joint_outlet_segments"][parameters["joint_outlet_names"][joint_name]])) + ''.join(" " + str(i) for i in parameters["joint_outlet_segments"][parameters["joint_outlet_names"][joint_name]]) + "\n\n")

    # write junction models
    zero_d_input_file.write("\n### DO NOT CHANGE THIS SECTION - generated automatically\n#\n# ==========\n# JUNCTION MODEL CARD\n# ==========\n# - joint_name\n# - junction_type, i.e. NORMAL_JUNCTION, etc.\n\n")
    junction_names = list(parameters["junction_types"].keys())
    junction_names.sort()
    for joint_name in junction_names:
        zero_d_input_file.write("JUNCTION_MODEL " + joint_name + " " + parameters["junction_types"][joint_name] + "\n")

    # write elements
    zero_d_input_file.write("\n# =============\n# ELEMENT CARD\n# =============\n# - 1d_segment_name\n# - 1d_segment_number / 0d_element_number\n# - average_radius\n# - length\n# - blood_density\n# - blood_dynamic_viscosity\n# - Eh/r\n# - inlet_node_number\n# - outlet_node_number\n# - vessel_vector_x_coor\n# - vessel_vector_y_coor\n# - vessel_vector_z_coor\n# - boundary_condition_type\n# - boundary_condition_datatable_name\n# - 0d_element_type, i.e. R, RC, RL, RCL, L, C, etc.\n# - 0d_element_values in format of: R or R C or R L or R C L or L or C or blank for custom 0d elements\n\n")
    for segment_number in parameters["segment_numbers_list"]:
        zero_d_input_file.write("ELEMENT " + parameters["segment_names"][segment_number] + " " + str(segment_number))
        zero_d_input_file.write(" " + str(parameters["radii"][segment_number]) + " ")
        zero_d_input_file.write(str(parameters["lengths"][segment_number]) + " ")
        zero_d_input_file.write(str(parameters["rho"][segment_number]) + " ")
        zero_d_input_file.write(str(parameters["mu"][segment_number]) + " ")
        zero_d_input_file.write(str(parameters["Ehor_func"][segment_number](parameters["radii"][segment_number])) + " ")
        zero_d_input_file.write(str(parameters["nodes_of_segments"][segment_number][0]) + " ") # inlet node of element
        zero_d_input_file.write(str(parameters["nodes_of_segments"][segment_number][1]) + " ") # outlet node of element
        zero_d_input_file.write(''.join(str(i) + " " for i in parameters["segment_tangent_vectors"][segment_number]))
        zero_d_input_file.write(parameters["boundary_condition_types"]["outlet"][segment_number] + " " + parameters["boundary_condition_datatable_names"]["outlet"][segment_number] + " ")
        zero_d_input_file.write(parameters["segment_0d_types"][segment_number])
        zero_d_input_file.write(''.join(" " + str(value) for value in parameters["segment_0d_values"][segment_number]) + "\n")

    # write datatables
    zero_d_input_file.write("\n# ===============\n# DATATABLE CARD\n# ===============\n# - datatable_name\n# - datatable_type\n# - datatable_values in format of:\n#\t\ttime0 value1\n#\t\ttime1 value1\n#\t\t...\n")
    list_of_bc_datatable_names = list(parameters["datatable_values"].keys())
    list_of_bc_datatable_names.sort()
    for dt_name in list_of_bc_datatable_names:
        zero_d_input_file.write("\nDATATABLE " + dt_name + " " + parameters["datatable_types"][dt_name] + "\n")
        for i in range(0, len(parameters["datatable_values"][dt_name]), 2):
            zero_d_input_file.write(str(parameters["datatable_values"][dt_name][i]) + " " + str(parameters["datatable_values"][dt_name][i + 1]) + "\n")
        zero_d_input_file.write("ENDDATATABLE\n")

    # write inflow bc
    zero_d_input_file.write("\n# ==============================\n# INLET BOUNDARY CONDITION CARD\n# ==============================\n# - inlet_segment_number\n# - inlet_boundary_condition_type\n# - inlet_boundary_condition_datatable_name\n\n")
    for i in range(0, len(parameters["inlet_segments_of_model"])):
        inlet_segment_number = parameters["inlet_segments_of_model"][i]
        zero_d_input_file.write("INLET_BOUNDARY_CONDITION " + str(inlet_segment_number) + " " + parameters["boundary_condition_types"]["inlet"][inlet_segment_number] + " " + parameters["boundary_condition_datatable_names"]["inlet"][inlet_segment_number] + "\n")

    zero_d_input_file.close()
    print("1d model has been successfully converted into a 0d model.")

def convert_1d_to_0d(one_d_solver_input_file_path, destination_folder, default_element_type, default_junction_type = "NORMAL_JUNCTION", one_d_inlet_segment_number = 0, check_default_element_type = False, check_default_junction_type = False, use_self_defined_name = False):
    """
    Purpose:
        Convert the 1d model (stored in the 1d solver input file) into its equivalent 0d model and create the 0d solver input file.
    Inputs:
        string one_d_solver_input_file_path
            = path to the 1d solver input file
        string destination_folder
            = path to directory in which you want to save the 0d input file and its auxiliary files
        string default_element_type
            = default 0d element type, i.e. "R", "RC", "RCL", "RL", "L", "C", etc. to be used for the majority of the segments
        string default_junction_type
            = default junction type, i.e. "NORMAL_JUNCTION", etc. to be used for the majority of the joints
        int one_d_inlet_segment_number
            = inlet segment number of the 1d model
        boolean check_default_element_type
            = True to give the user the option to change the types of some 0d elements, instead of using the default type for all segments
        boolean check_default_junction_type
            = True to give the user the option to change the junction type used for each joint; False otherwise
        boolean use_self_defined_name
            = True to allow the user to specify the name of the 0d solver input file; False to automatically make the 0d solver input file name the same as the 1d solver input file name (but with an extra '_0d' appended)
    Returns:
        string zero_d_solver_input_file_path
            = path to the 0d solver input file
    """
    parameters = extract_info_from_solver_input_file(one_d_solver_input_file_path, one_d_inlet_segment_number)
    zero_d_solver_input_file_path, segment_0d_types_input_file_path, junction_types_input_file_path = get_zero_solver_input_file_path(one_d_solver_input_file_path, destination_folder, use_self_defined_name)
    check_default_0d_element_type_with_user(parameters, default_element_type, segment_0d_types_input_file_path, check_default_element_type)
    check_default_junction_type_with_user(parameters, default_junction_type, junction_types_input_file_path, check_default_junction_type)
    create_0d_elements(parameters)
    create_0d_solver_input_file(parameters, zero_d_solver_input_file_path)
    return zero_d_solver_input_file_path

if __name__ == "__main__":
    """
    Inputs:
        string one_d_solver_input_file_path
            = path to the 1d solver input file
        string destination_folder
            = path to directory in which you want to save the 0d input file and its auxiliary files
        string default_element_type
            = default 0d element type, i.e. "R", "RC", "RCL", "RL", "L", "C", etc. to be used for the majority of the segments
        string default_junction_type
            = default junction type, i.e. "NORMAL_JUNCTION", etc. to be used for the majority of the joints
        int one_d_inlet_segment_number
            = inlet segment number of the 1d model
        boolean check_default_element_type
            = True to give the user the option to change the types of some 0d elements, instead of using the default type for all segments
        boolean check_default_junction_type
            = True to give the user the option to change the junction type used for each joint; False otherwise
        boolean use_self_defined_name
            = True to allow the user to specify the name of the 0d solver input file; False to automatically make the 0d solver input file name the same as the 1d solver input file name (but with an extra '_0d' appended)
    Outputs:
        string zero_d_solver_input_file_path
            = path to the 0d solver input file
    """

    destination_folder = "./zero_d_solver_test_cases/test_condense_0d/condense_0d_tests/coarse"
    one_d_input_file_name = "0003_0001.inp"
    one_d_solver_input_file_path = os.path.join(destination_folder, one_d_input_file_name)
    default_element_type = "R"
    default_junction_type = "NORMAL_JUNCTION"
    one_d_inlet_segment_number = 0
    check_default_element_type = False
    check_default_junction_type = False
    use_self_defined_name = False

    zero_d_solver_input_file_path = convert_1d_to_0d(one_d_solver_input_file_path, destination_folder, default_element_type, default_junction_type, one_d_inlet_segment_number, check_default_element_type, check_default_junction_type, use_self_defined_name)
