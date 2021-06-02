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

import sys
import os
import numpy as np

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
