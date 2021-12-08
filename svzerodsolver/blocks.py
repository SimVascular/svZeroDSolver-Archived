# coding=utf-8

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

import numpy as np
from collections import defaultdict


class Wire:
    """
    Wires connect LPNBlocks. Wires contain only a single pressure and flow value (which are part of the system variables). Wires must be attached to two LPNBlocks, one at each end.
    """
    def __init__(self, name = "NoNameWire"):
        self.name = name
        self.lpn_solution_ids = [None] * 2


class LPNBlock:
    def __init__(self, connecting_block_list = [], name = "NoName", flow_directions = []):
        self.name = name
        self.connecting_block_list = connecting_block_list
        self.n_eqn = 2
        self.n_block_var = 0
        self.connecting_wires_list = []

        # flow direction is -1 for flow into this block and +1 for flow out of this block
        self.flow_directions = flow_directions

        # solution IDs for the LPN block's internal solution variables
        self.lpn_solution_ids = []

        # block matrices
        self.mat = defaultdict(list)

        # row and column indices of block in global matrix
        self.global_col_id = []
        self.global_row_id = []

    def add_connecting_block(self, block, direction):
        # direction = +1 if flow sent from self to block
        #           = -1 if flow sent from block to self
        self.connecting_block_list.append(block)
        self.flow_directions.append(direction)

    def add_connecting_wire(self, new_wire):
        self.connecting_wires_list.append(new_wire)

    def update_constant(self):
        """
        Update solution- and time-independent blocks
        """
        pass

    def update_time(self, args):
        """
        Update time-dependent blocks
        """
        pass

    def update_solution(self, args):
        """
        Update solution-dependent blocks
        """
        pass

    def get_global_equation_ids(self, wire_dict, local_solution_id):
        """ 
        Return the index at which the local solution variable resides in the global vector of solution variables.
        """
        
        n_wire_var = len(self.connecting_block_list) * 2  # there are 2 soltns (P and Q) per wire
        if local_solution_id < n_wire_var: # wire solution variable
            var_type = local_solution_id % 2  # 0 --> P, 1 --> Q
            local_wire_id = int(local_solution_id / 2)

            # Example: Assume len(self.connecting_block_list) is 2. This can be a normal resistor block, which has 2 connections. Then this R block has 2 connecting wires. thus, This R block has 4 related solution variables (P_in, Q_in, P_out, Q_out).
            #     then for these are the var_types we get for each local_solution_id:
            #  local_solution_id  : var_type :  local_wire_id
            #         0           :     0    :       0        <---  var_type = pressure, local_wire_id = inlet wire
            #         1           :     1    :       0        <---  var_type = flow,     local_wire_id = inlet wire
            #         2           :     0    :       1        <---  var_type = pressure, local_wire_id = outlet wire
            #         3           :     1    :       1        <---  var_type = flow,     local_wire_id = outlet wire

            return wire_dict[self.connecting_wires_list[local_wire_id]].lpn_solution_ids[var_type]
        else: # block internal solution variable
            internal_var_local_id = local_solution_id - n_wire_var
            return self.lpn_solution_ids[internal_var_local_id]


class Junction(LPNBlock):
    """
    Junction points between LPN blocks with specified directions of flow
    """
    def __init__(self, connecting_block_list = None, name = "NoNameJunction", flow_directions = None):
        LPNBlock.__init__(self, connecting_block_list, name = name, flow_directions = flow_directions)
        self.n_eqn = len(self.connecting_block_list)  # number of equations = num of blocks that connect to this junction, where the equations are 1) mass conservation 2) inlet pressures = outlet pressures

    def add_connecting_block(self, block, direction):
        self.connecting_block_list.append(block)
        self.n_eqn = len(self.connecting_block_list)
        self.flow_directions.append(direction)

    def update_constant(self):
        # Number of variables per tuple = 2 * len(self.connecting_block_list)
        # Number of equations = len(self.connecting_block_list) - 1 Pressure equations, 1 flow equation
        # Format : P1,Q1,P2,Q2,P3,Q3, .., Pn,Qm
        self.mat['F'] = [(1.,) + (0,) * (2 * i + 1) + (-1,) + (0,) * (2 * len(self.connecting_block_list) - 2 * i - 3) for i in
                         range(len(self.connecting_block_list) - 1)]

        tmp = (0,)
        for d in self.flow_directions[:-1]:
            tmp += (d,)
            tmp += (0,)

        tmp += (self.flow_directions[-1],)
        self.mat['F'].append(tmp)


class BloodVessel(LPNBlock):
    """
    Stenosis:
        equation: delta_P = ( K_t * rho / ( 2 * (A_0)**2 ) ) * ( ( A_0 / A_s ) - 1 )**2 * Q * abs(Q) + R_poiseuille * Q
                          =               stenosis_coefficient                          * Q * abs(Q) + R_poiseuille * Q

        source: Mirramezani, M., Shadden, S.C. A distributed lumped parameter model of blood flow. Annals of Biomedical Engineering. 2020.
    """
    def __init__(self, R, C, L, stenosis_coefficient, connecting_block_list = None, name = "NoNameBloodVessel", flow_directions = None):
        LPNBlock.__init__(self, connecting_block_list, name = name, flow_directions = flow_directions)
        self.R = R  # poiseuille resistance value = 8 * mu * L / (pi * r**4)
        self.C = C
        self.L = L
        self.stenosis_coefficient = stenosis_coefficient

    # the ordering of the solution variables is : (P_in, Q_in, P_out, Q_out)

    def update_constant(self):
        self.mat['E'] = [(0, 0, 0, -self.L), (-self.C, self.C * self.R, 0, 0)]

    def update_solution(self, args):
        curr_y = args['Solution']  # the current solution for all unknowns in our 0D model
        wire_dict = args['Wire dictionary']
        Q_in = curr_y[wire_dict[self.connecting_wires_list[0]].lpn_solution_ids[1]]
        self.mat['F'] = [(1.0, -1.0 * self.stenosis_coefficient * np.abs(Q_in) - self.R, -1.0, 0), (0, 1.0, 0, -1.0)]
        self.mat['dF'] = [(0, -1.0 * self.stenosis_coefficient * np.abs(Q_in), 0, 0), (0,) * 4]


class UnsteadyResistanceWithDistalPressure(LPNBlock):
    def __init__(self, Rfunc, Pref_func, connecting_block_list = None, name = "NoNameUnsteadyResistanceWithDistalPressure", flow_directions=None):
        LPNBlock.__init__(self, connecting_block_list, name = name, flow_directions = flow_directions)
        self.n_eqn = 1
        self.Rfunc = Rfunc
        self.Pref_func = Pref_func

    def update_time(self, args):
        """
        the ordering is : (P_in,Q_in)
        """
        t = args['Time']
        self.mat['F'] = [(1., -1.0 * self.Rfunc(t))]
        self.mat['C'] = [-1.0 * self.Pref_func(t)]


class UnsteadyPressureRef(LPNBlock):
    """
    Unsteady P reference
    """
    def __init__(self, Pfunc, connecting_block_list = None, name = "NoNameUnsteadyPressureRef", flow_directions = None):
        LPNBlock.__init__(self, connecting_block_list, name = name, flow_directions = flow_directions)
        self.n_eqn = 1
        self.Pfunc = Pfunc

    def update_time(self, args):
        t = args['Time']
        self.mat['C'] = [-1.0 * self.Pfunc(t)]

    def update_constant(self):
        self.mat['F'] = [(1., 0.)]


class UnsteadyFlowRef(LPNBlock):
    """
    Flow reference
    """
    def __init__(self, Qfunc, connecting_block_list = None, name = "NoNameUnsteadyFlowRef", flow_directions = None):
        LPNBlock.__init__(self, connecting_block_list, name = name, flow_directions = flow_directions)
        self.n_eqn = 1
        self.Qfunc = Qfunc

    def update_time(self, args):
        t = args['Time']
        self.mat['C'] = [-1.0 * self.Qfunc(t)]

    def update_constant(self):
        self.mat['F'] = [(0, 1.)]


class UnsteadyRCRBlockWithDistalPressure(LPNBlock):
    """
    Unsteady RCR - time-varying RCR values
    Formulation includes additional variable : internal pressure proximal to capacitance.
    """
    def __init__(self, Rp_func, C_func, Rd_func, Pref_func, connecting_block_list = None, name = "NoNameUnsteadyRCRBlockWithDistalPressure", flow_directions = None):
        LPNBlock.__init__(self, connecting_block_list, name = name, flow_directions = flow_directions)
        self.n_eqn = 2
        self.n_block_var = 1
        self.Rp_func = Rp_func
        self.C_func = C_func
        self.Rd_func = Rd_func
        self.Pref_func = Pref_func

    def update_time(self, args):
        """
        unknowns = [P_in, Q_in, internal_var (Pressure at the intersection of the Rp, Rd, and C elements)]
        """
        t = args['Time']
        self.mat['E'] = [(0, 0, 0), (0, 0, -1.0 * self.Rd_func(t) * self.C_func(t))]
        self.mat['F'] = [(1., -self.Rp_func(t), -1.), (0.0, self.Rd_func(t), -1.0)]
        self.mat['C'] = [0, self.Pref_func(t)]


class OpenLoopCoronaryWithDistalPressureBlock(LPNBlock):
    """
    open-loop coronary BC = RCRCR BC
    Publication reference: Kim, H. J. et al. Patient-specific modeling of blood flow and pressure in human coronary arteries. Annals of Biomedical Engineering 38, 3195–3209 (2010)."
    """
    def __init__(self, Ra, Ca, Ram, Cim, Rv, Pim, Pv, cardiac_cycle_period, connecting_block_list = None, name = "NoNameCoronary", flow_directions = None):
        LPNBlock.__init__(self, connecting_block_list, name = name, flow_directions = flow_directions)
        self.n_eqn = 2
        self.n_block_var = 1
        self.Ra = Ra
        self.Ca = Ca
        self.Ram = Ram
        self.Cim = Cim
        self.Rv = Rv
        self.Pa = 0.0
        self.Pim = Pim
        self.Pv = Pv
        self.cardiac_cycle_period = cardiac_cycle_period

    def get_P_at_t(self, P, t):
        tt = P[:, 0]
        P_val = P[:, 1]
        ti, td = divmod(t, self.cardiac_cycle_period)
        P_tt = np.interp(td, tt, P_val)
        return P_tt

    def update_time(self, args):
        # For this open-loop coronary BC, the ordering of solution unknowns is : (P_in, Q_in, V_im)
        # where V_im is the volume of the second capacitor, Cim
        # Q_in is the flow through the first resistor
        # and P_in is the pressure at the inlet of the first resistor
        ttt = args['Time']
        Pim_value = self.get_P_at_t(self.Pim, ttt)
        Pv_value = self.get_P_at_t(self.Pv, ttt)
        self.mat['C'] = [-1.0 * self.Cim * Pim_value + self.Cim * Pv_value,
                         -1.0 * self.Cim * (self.Rv + self.Ram) * Pim_value + self.Ram * self.Cim * Pv_value]

    def update_constant(self):
        self.mat['E'] = [
            (-1.0 * self.Ca * self.Cim * self.Rv, self.Ra * self.Ca * self.Cim * self.Rv, -1.0 * self.Cim * self.Rv),
            (0.0, 0.0, -1.0 * self.Cim * self.Rv * self.Ram)]
        self.mat['F'] = [(0.0, self.Cim * self.Rv, -1.0),
                         (self.Cim * self.Rv, -1.0 * self.Cim * self.Rv * self.Ra, -1.0 * (self.Rv + self.Ram))]