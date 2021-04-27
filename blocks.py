# coding=utf-8

import numpy as np
from collections import defaultdict


class LPNVariable:
    def __init__(self, value, units, name="NoName", vtype='ArbitraryVariable'):
        self.type = vtype
        self.value = value
        # Two generic units accepted : SI, cgs. Conversion for special values applied
        self.units = units
        self.name = name


class PressureVariable(LPNVariable):

    def __init__(self, value, units='cgs', name='NoNamePressure'):
        LPNVariable.__init__(self, value=value, units=units, name=name, vtype='Pressure')

    def convert_to_cgs(self):
        if self.units == 'cgs':
            print("Variable: " + self.name + " already at cgs")
        elif self.units == 'SI':
            self.value *= 1.0E5
            self.units = 'cgs'
        elif self.units == 'mmHg':
            self.value *= 0.001333224
            self.units = 'cgs'
        else:
            raise Exception("Units " + self.units + " not recognized")

    def convert_to_mmHg(self):
        if self.units == 'cgs':
            self.value *= 750.06
            self.units = 'mmHg'
        elif self.units == 'SI':
            self.value = self.value * 7.50 * 1E-3
            self.units = 'mmHg'
        elif self.units == 'mmHg':
            print("Variable: " + self.name + " already at mmHg")
        else:
            raise Exception("Units " + self.units + " not recognized")


class FlowVariable(LPNVariable):

    def __init__(self, value, units='cgs', name='NoNameFlow'):
        LPNVariable.__init__(self, value=value, units=units, name=name, vtype='Flow')

    def convert_to_cgs(self):
        if self.units == 'cgs':
            print("Variable: " + self.name + " already at cgs")
        elif self.units == 'SI':
            self.value = self.value * 1.0E-6
            self.units = 'cgs'
        elif self.units == 'Lpm':  # litres per minute
            self.value = self.value * 16.6667
            self.units = 'cgs'
        else:
            raise Exception("Units " + self.units + " not recognized")

    def convert_to_Lpm(self):
        if self.units == 'cgs':
            self.value = self.value / 16.6667
            self.units = 'Lpm'
        elif self.units == 'SI':
            self.value = self.value / (16.6667 * 1.0E-6)
            self.units = 'Lpm'
        elif self.units == 'Lpm':
            print("Variable: " + self.name + " already at Lpm")
        else:
            raise Exception("Units " + self.units + " not recognized")


class wire:
    """
    Wires connect circuit elements and junctions
    They can only posses a single pressure and flow value (system variables)
    They can also only possess one element(or junction) at each end
    """
    def __init__(self, connecting_elements, Pval=0, Qval=0, name="NoNameWire", P_units='cgs', Q_units='cgs'):
        self.name = name
        self.type = 'Wire'
        self.P = PressureVariable(value=Pval, units=P_units, name=name + "_P")
        self.Q = FlowVariable(value=Qval, units=Q_units, name=name + "_Q")
        if len(connecting_elements) > 2:
            raise Exception('Wire cannot connect to more than two elements at a time. Use a junction LPN block')
        if type(connecting_elements) != tuple:
            raise Exception('Connecting elements to wire should be passed as a 2-tuple')
        self.connecting_elements = connecting_elements
        self.LPN_solution_ids = [None] * 2


class LPNBlock:
    def __init__(self, connecting_block_list=None, name="NoName", flow_directions=[]):
        if connecting_block_list == None:
            connecting_block_list = []
        self.connecting_block_list = connecting_block_list
        self.num_connections = len(connecting_block_list)
        self.name = name
        self.neq = 2
        self.n_connect = 2
        self.n_connect = None
        self.type = "ArbitraryBlock"
        self.num_block_vars = 0
        self.connecting_wires_list = []

        # -1 : Inflow to block, +1 outflow from block
        self.flow_directions = flow_directions

        # solution IDs for the LPN block's internal solution variables
        self.LPN_solution_ids = []

        # block matrices
        self.mat = defaultdict(list)

        # row and column indices of block in global matrix
        self.global_col_id = []
        self.global_row_id = []

    def check_block_consistency(self):
        if len(connecting_block_list) != self.n_connect:
            msg = self.name + " block can be connected only to " + str(self.n_connect) + " elements"
            raise Exception(msg)

    def add_connecting_block(self, block, direction):
        # Direction = +1 if flow sent to block
        #            = -1 if flow recvd from block
        self.connecting_block_list.append(block)
        self.num_connections = len(self.connecting_block_list)
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

    def eqids(self, wire_dict, local_eq):
        # EqID returns variable's location in solution vector

        nwirevars = self.num_connections * 2  # num_connections is multipled by 2 because each wire has 2 soltns (P and Q)
        if local_eq < nwirevars:
            vtype = local_eq % 2  # 0 --> P, 1 --> Q
            wnum = int(local_eq / 2)

            # example: assume num_connections is 2. this can be a normal resistor block, which has 2 connections. then this R block has 2 connecting wires. thus, this R block has 4 related solution variables/unknowns (P_in, Q_in, P_out, Q_out). note that local_eq = local ID.
            #     then for these are the vtypes we get for each local_eq:
            #         local_eq    :     vtype     :     wnum
            #         0            :     0        :    0        <---    vtype = pressure, wnum = inlet wire
            #         1            :    1        :    0        <---    vtype = flow, wnum = inlet wire
            #         2            :    0        :    1        <---    vtype = pressure, wnum = outlet wire
            #         3            :    1        :    1        <---    vtype = flow, wnum = outlet wire
            #    note that vtype represents whether the solution variable in local_eq (local ID) is a P or Q solution
            #        and wnum represents whether the solution variable in local_eq comes from the inlet wire or the outlet wire, for this LPNBlock with 2 connections (one inlet, one outlet)

            return wire_dict[self.connecting_wires_list[wnum]].LPN_solution_ids[vtype]
        else:  # this section will return the index at which the LPNBlock's  INTERNAL SOLUTION VARIABLES are stored in the global vector of solution unknowns/variables (i.e. I think RCR and OpenLoopCoronaryBlock have internal solution variables; these internal solution variables arent the P_in, Q_in, P_out, Q_out that correspond to the solutions on the attached wires, they are the solutions that are internal to the LPNBlock itself)
            vnum = local_eq - nwirevars
            return self.LPN_solution_ids[vnum]


class Junction(LPNBlock):
    """
    Junction points between LPN blocks with specified directions of flow
    """
    def __init__(self, connecting_block_list=None, name="NoNameJunction", flow_directions=None):
        LPNBlock.__init__(self, connecting_block_list, name=name, flow_directions=flow_directions)
        self.type = "Junction"
        self.neq = self.num_connections  # number of equations = num of blocks that connect to this junction, where the equations are 1) mass conservation 2) inlet pressures = outlet pressures

    def add_connecting_block(self, block, direction):
        self.connecting_block_list.append(block)
        self.num_connections = len(self.connecting_block_list)
        self.neq = self.num_connections
        self.flow_directions.append(direction)

    def update_constant(self):
        # Number of variables per tuple = 2*num_connections
        # Number of equations = num_connections-1 Pressure equations, 1 flow equation
        # Format : P1,Q1,P2,Q2,P3,Q3, .., Pn,Qm
        self.mat['F'] = [(1.,) + (0,) * (2 * i + 1) + (-1,) + (0,) * (2 * self.num_connections - 2 * i - 3) for i in
                         range(self.num_connections - 1)]

        tmp = (0,)
        for d in self.flow_directions[:-1]:
            tmp += (d,)
            tmp += (0,)

        tmp += (self.flow_directions[-1],)
        self.mat['F'].append(tmp)


class Resistance(LPNBlock):
    """
    Resistance element
    """
    def __init__(self, R, connecting_block_list=None, name="NoNameResistance", flow_directions=None):
        LPNBlock.__init__(self, connecting_block_list, name=name, flow_directions=flow_directions)
        self.type = "Resistance"
        self.R = R

    def update_constant(self):
        """
        For resistors, the ordering is : (P_in,Q_in,P_out,Q_out)
        """
        self.mat['F'] = [(1., -1. * self.R, -1., 0), (0, 1., 0, -1.)]


class StenosisBlock(LPNBlock):
    """
    Stenosis:
        equation:   delta_P = ( K_t * rho / ( 2 * (A_0)**2 ) ) * ( ( A_0 / A_s ) - 1 )**2 * Q * abs(Q)
                            =               stenosis_coefficient                          * Q * abs(Q)

        source: Mirramezani, M., Shadden, S.C. A distributed lumped parameter model of blood flow. Annals of Biomedical Engineering. 2020.
    """
    def __init__(self, R, stenosis_coefficient, connecting_block_list=None, name="NoNameRCL", flow_directions=None):
        LPNBlock.__init__(self, connecting_block_list, name=name, flow_directions=flow_directions)
        self.type = "Stenosis"
        self.R = R  # poiseuille resistance value = 8 * mu * L / (pi * r**4)
        self.stenosis_coefficient = stenosis_coefficient

    def update_solution(self, args):
        curr_y = args['Solution']  # the current solution for all unknowns in our 0D model
        wire_dict = args['Wire dictionary']
        Q_in = curr_y[wire_dict[self.connecting_wires_list[0]].LPN_solution_ids[1]]
        self.mat['F'] = [(1.0, -1.0 * self.stenosis_coefficient * np.abs(Q_in) - self.R, -1.0, 0), (0, 1.0, 0, -1.0)]
        self.mat['dF'] = [(0, -1.0 * self.stenosis_coefficient * np.abs(Q_in), 0, 0), (0,) * 4]


class UnsteadyResistance(LPNBlock):
    """
    Unsteady Resistance : delta_P = q*Rfunc(t)
    """
    def __init__(self, Rfunc, connecting_block_list=None, name="NoNameUnsteadyResistance", flow_directions=None):
        LPNBlock.__init__(self, connecting_block_list, name=name, flow_directions=flow_directions)
        self.type = "UnsteadyResistance"
        self.Rfunc = Rfunc

    def update_time(self, args):
        """
        For resistors, the ordering is : (P_in,Q_in,P_out,Q_out)
        """
        t = args['Time']
        self.mat['F'] = [(1., -1.0 * self.Rfunc(t), -1., 0), (0, 1., 0, -1.)]


class UnsteadyResistanceWithDistalPressure(LPNBlock):
    def __init__(self, Rfunc, Pref_func, connecting_block_list=None, name="NoNameUnsteadyResistanceWithDistalPressure",
                 flow_directions=None):
        LPNBlock.__init__(self, connecting_block_list, name=name, flow_directions=flow_directions)
        self.type = "UnsteadyResistanceWithDistalPressure"
        self.neq = 1
        self.Rfunc = Rfunc
        self.Pref_func = Pref_func

    def update_time(self, args):
        """
        the ordering is : (P_in,Q_in)
        """
        t = args['Time']
        self.mat['F'] = [(1., -1.0 * self.Rfunc(t))]
        self.mat['C'] = [-1.0 * self.Pref_func(t)]


class PressureRef(LPNBlock):
    """
    Pressure reference
    """
    def __init__(self, Pref, connecting_block_list=None, name="NoNamePressureRef", flow_directions=None):
        LPNBlock.__init__(self, connecting_block_list, name=name, flow_directions=flow_directions)
        self.type = "PressureRef"
        self.neq = 1
        self.Pref = Pref

    def update_constant(self):
        self.mat['F'] = [(1.,)]
        self.mat['C'] = [-1.0 * self.Pref]


class UnsteadyPressureRef(LPNBlock):
    """
    Unsteady P reference
    """
    def __init__(self, Pfunc, connecting_block_list=None, name="NoNameUnsteadyPressureRef", flow_directions=None):
        LPNBlock.__init__(self, connecting_block_list, name=name, flow_directions=flow_directions)
        self.type = "UnsteadyPressureRef"
        self.neq = 1
        self.n_connect = 1
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
    def __init__(self, Qfunc, connecting_block_list=None, name="NoNameUnsteadyFlowRef", flow_directions=None):
        LPNBlock.__init__(self, connecting_block_list, name=name, flow_directions=flow_directions)
        self.type = "UnsteadyFlowRef"
        self.neq = 1
        self.n_connect = 1
        self.Qfunc = Qfunc

    def update_time(self, args):
        t = args['Time']
        self.mat['C'] = [-1.0 * self.Qfunc(t)]

    def update_constant(self):
        self.mat['F'] = [(0, 1.)]


class Capacitance(LPNBlock):
    """
    Capacitance
    """
    def __init__(self, C, connecting_block_list=None, name="NoNameCapacitance", flow_directions=None):
        LPNBlock.__init__(self, connecting_block_list, name=name, flow_directions=flow_directions)
        self.type = "Capacitance"
        self.C = C

    def update_constant(self):
        self.mat['E'] = [(1.0 * self.C, 0, -1.0 * self.C, 0), (0, 0, 0, 0)]
        self.mat['F'] = [(0, -1.0, 0, 0), (0, 1., 0, -1.)]


class RCLBlock(LPNBlock):
    """
    RCL - constant resistor, capacitor, inductor - vessel representation
    Formulation includes additional variable : internal pressure proximal to capacitance.
    """
    def __init__(self, R, C, L, connecting_block_list=None, name="NoNameRCL", flow_directions=None):
        LPNBlock.__init__(self, connecting_block_list, name=name, flow_directions=flow_directions)
        self.type = "RCL"
        self.neq = 3
        self.num_block_vars = 1
        self.R = R
        self.C = C
        self.L = L

    def update_constant(self):
        self.mat['E'] = [(0, 0, 0, -self.L, 0), (0, 0, 0, 0, -self.C), (0, 0, 0, 0, 0)]
        self.mat['F'] = [(1., -self.R, -1., 0, 0), (0, 1., 0, -1., 0), (1., -self.R, 0, 0, -1.)]


class RCBlock(LPNBlock):
    """
    RC - constant resistor, capacitor - low inertia vessel
    """
    def __init__(self, R, C, connecting_block_list=None, name="NoNameRC", flow_directions=None):
        LPNBlock.__init__(self, connecting_block_list, name=name, flow_directions=flow_directions)
        self.type = "RC"
        self.R = R
        self.C = C

    def update_constant(self):
        self.mat['E'] = [(0, 0, 0, 0), (0, 0, -self.C, 0)]
        self.mat['F'] = [(1.0, -self.R, -1.0, 0), (0, 1., 0, -1.)]


class RLBlock(LPNBlock):
    """
    RL - constant resistor, inductor
    """
    def __init__(self, R, L, connecting_block_list=None, name="NoNameRL", flow_directions=None):
        LPNBlock.__init__(self, connecting_block_list, name=name, flow_directions=flow_directions)
        self.type = "RL"
        self.R = R
        self.L = L

    def update_constant(self):
        """
        the ordering of solution unknowns is : (P_in, Q_in, P_out, Q_out)
        """
        self.mat['E'] = [(0, 0, 0, -self.L), (0, 0, 0, 0)]
        self.mat['F'] = [(1.0, -self.R, -1.0, 0), (0, 1., 0, -1.)]


class RCRBlock(LPNBlock):
    """
    RCR - constant RCR - outflow representation
    Formulation includes additional variable : internal pressure proximal to capacitance.
    """
    def __init__(self, Rp, C, Rd, connecting_block_list=None, name="NoNameRCR", flow_directions=None):
        LPNBlock.__init__(self, connecting_block_list, name=name, flow_directions=flow_directions)
        self.type = "RCR"
        self.neq = 3
        self.num_block_vars = 1
        self.Rp = Rp
        self.C = C
        self.Rd = Rd

    def update_constant(self):
        self.mat['E'] = [(0, 0, 0, 0, 0), (0, 0, 0, 0, -self.C), (0, 0, 0, 0, 0)]
        self.mat['F'] = [(1.0, -self.Rp, -1.0, -self.Rd, 0), (0, 1., 0, -1., 0), (1., -self.Rp, 0, 0, -1.)]


class UnsteadyRCRBlock(LPNBlock):
    """
    Unsteady RCR - time-varying RCR values
    Formulation includes additional variable : internal pressure proximal to capacitance.
    """
    def __init__(self, Rp_func, C_func, Rd_func, connecting_block_list=None, name="NoNameUnsteadyRCR",
                 flow_directions=None):
        LPNBlock.__init__(self, connecting_block_list, name=name, flow_directions=flow_directions)
        self.type = "UnsteadyRCR"
        self.neq = 3
        self.num_block_vars = 1
        self.Rp_func = Rp_func
        self.C_func = C_func
        self.Rd_func = Rd_func

    def update_time(self, args):
        t = args['Time']
        self.mat['F'] = [(1.0, -self.Rp_func(t), -1.0, -self.Rd_func(t), 0), (0, 1., 0, -1., 0),
                         (1., -self.Rp_func(t), 0, 0, -1.)]

    def update_constant(self):
        self.mat['E'] = [(0, 0, 0, 0, 0), (0, 0, 0, 0, -self.C_func(t)), (0, 0, 0, 0, 0)]


class UnsteadyRCRBlockWithDistalPressure(LPNBlock):
    """
    Unsteady RCR - time-varying RCR values
    Formulation includes additional variable : internal pressure proximal to capacitance.
    """
    def __init__(self, Rp_func, C_func, Rd_func, Pref_func, connecting_block_list=None,
                 name="NoNameUnsteadyRCRBlockWithDistalPressure", flow_directions=None):
        LPNBlock.__init__(self, connecting_block_list, name=name, flow_directions=flow_directions)
        self.type = "UnsteadyRCRBlockWithDistalPressure"
        self.neq = 2
        self.n_connect = 1
        self.num_block_vars = 1
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
    def __init__(self, Ra, Ca, Ram, Cim, Rv, Pim, Pv, cardiac_cycle_period, connecting_block_list=None,
                 name="NoNameCoronary", flow_directions=None):
        LPNBlock.__init__(self, connecting_block_list, name=name, flow_directions=flow_directions)
        self.type = "OpenLoopCoronaryWithDistalPressureBlock_v2"
        self.neq = 2
        self.n_connect = 1
        self.num_block_vars = 1
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


class TimeDependentCapacitance(LPNBlock):
    """
    Time Varying Capacitance
    """
    def __init__(self, Cfunc, connecting_block_list=None, name="NoNameTimeDependentCapacitance", flow_directions=None):
        LPNBlock.__init__(self, connecting_block_list, name=name, flow_directions=flow_directions)
        self.type = "TimeDependentCapacitance"
        self.Cfunc = Cfunc

    def update_time(self, args):
        t = args['Time']
        self.mat['E'] = [(1.0 * self.Cfunc(t), 0, -1.0 * self.Cfunc(t), 0), (0, 0, 0, 0)]

    def update_constant(self):
        self.mat['F'] = [(0, -1.0, 0, 0), (0, 1., 0, -1.)]


class ChamberModel(LPNBlock):
    """
    Chamber Capacitance -- with direct prescription of pressure
    """
    def __init__(self, Activefunc, Passivefunc, Pref, connecting_block_list=None, name="NoNameChamber",
                 flow_directions=None):

        LPNBlock.__init__(self, connecting_block_list, name=name, flow_directions=flow_directions)
        self.neq = 3
        self.type = "Chamber"
        self.num_block_vars = 1  # Chamber volume
        self.Activefunc = Activefunc
        self.Passivefunc = Passivefunc
        self.Pref = Pref

    def linearize_func(self, t, v, func, eps_dv=1e-4):
        # Input  : func(t,v), typically t = t_n+1 ie t_n + dt
        # Return : d/dv (func) at (t,v) using central differences
        dv = eps_dv * v
        return (func(t, v + dv) - func(t, v - dv)) / (2. * dv)

    def linearize_func_t(self, t, v, func, dt=1e-3):
        return (func(t + dt, v) - func(t - dt, v)) / (2. * dt)

    def initialize_volume(self, sol_vec, Vu=1.0, scale_fact=1.5):
        sol_vec[self.LPN_solution_ids[0]] = scale_fact * Vu

    def update_solution(self, args):
        t = args['Time']
        curr_y = args['Solution']
        # dt = args['Time step']
        # rho = args['rho']
        # alpha_f = 1/(1+rho)

        wire_dict = args['Wire dictionary']
        qi = curr_y[wire_dict[self.connecting_wires_list[0]].LPN_solution_ids[1]]
        qo = curr_y[wire_dict[self.connecting_wires_list[1]].LPN_solution_ids[1]]

        v = curr_y[self.LPN_solution_ids[0]]

        # print self.name,' has volume: ',v
        if v <= 0:
            print("Time: ", t)
            raise Exception("Zero chamber volume detected in chamber: ", self.name)

        a_lin = self.linearize_func(t, v, self.Activefunc)
        p_lin = self.linearize_func(t, v, self.Passivefunc)

        # c_from_lin = -(self.Activefunc(t,v) + self.Passivefunc(t,v)) + v*(a_lin+p_lin)
        # print 'Volume ',v,'has pressure ',(self.Activefunc(t,v) + self.Passivefunc(t,v))/1333.33

        self.mat['E'] = [(0, 0, 0, 0, 1.), (0, 0, 0, 0, 0), (0, 0, 0, 0, 0)]
        self.mat['F'] = [(0, -1., 0, 1., 0), (-1., 0, 1., 0, 0), (1., 0, 0., 0., 0.)]
        self.mat['C'] = [0, 0, -(self.Activefunc(t, v) + self.Passivefunc(t, v)) - self.Pref]

        # self.mat['F'] = [(0,-1.,0,1.,0), (-1.,0,1.,0,0),  (1.,0,0.,0,-(a_lin+p_lin))]
        # self.mat['C'] = [0,0,c_from_lin+self.Pref]

        self.mat['dC'] = [(0,) * 5, (0,) * 5, (0, 0, 0, 0, -(a_lin + p_lin))]


class Inductance(LPNBlock):
    """
    Inductance
    """
    def __init__(self, L, connecting_block_list=None, name="NoNameInductance", flow_directions=None):
        LPNBlock.__init__(self, connecting_block_list, name=name, flow_directions=flow_directions)
        self.type = "Inductance"
        self.L = L

    def update_constant(self):
        self.mat['E'] = [(0.0, -1., 0.0, 0), (0, 0, 0, 0)]
        self.mat['F'] = [(1. / self.L, 0.0, -1. / self.L, 0), (0, 1., 0, -1.)]


class IdealDiode2(LPNBlock):
    """
    Ideal diode - state variable
    """
    def __init__(self, connecting_block_list=None, eps=1e-17, name="NoNameIdealDiode", flow_directions=None):

        LPNBlock.__init__(self, connecting_block_list, name=name, flow_directions=flow_directions)
        self.type = "IdealDiode"
        self.neq = 3
        self.num_block_vars = 1  # State : 0- shunt, 1- interrupt
        self.eps = eps
        # self.Peps = Peps

    def update_solution(self, args):
        # Needs to take current solution vector as argument
        t = args['Time']
        curr_y = args['Solution']
        wire_dict = args['Wire dictionary']
        rho = args['rho']
        dt = args['Time step']
        alpha_f = 1.0 / (1.0 + rho)
        Qi = curr_y[wire_dict[self.connecting_wires_list[0]].LPN_solution_ids[1]]
        Qo = curr_y[wire_dict[self.connecting_wires_list[1]].LPN_solution_ids[1]]
        Pi = curr_y[wire_dict[self.connecting_wires_list[0]].LPN_solution_ids[0]]
        Po = curr_y[wire_dict[self.connecting_wires_list[1]].LPN_solution_ids[0]]
        state = curr_y[self.LPN_solution_ids[0]]

        if state > 0:
            # Zero flow
            self.mat['F'] = [(0, 1, 0, 0, 0), (0, 0, 0, 1., 0)]
            self.mat['C'] = [0., 0.]

        else:
            # Perfect wire
            self.mat['F'] = [(1., 0, -1., 0, 0), (0, 1., 0, -1., 0)]
            self.mat['C'] = [0, 0]

        if Qi < -self.eps or Qo < -self.eps:
            self.mat['F'].append((0, 0, 0, 0, 1))
            self.mat['C'].append(-alpha_f - (1 - alpha_f) * state)
            # self.mat['C'].append(-1)
            # print 'Flow toggle'
            # print state

        elif Pi - Po > -self.eps:
            self.mat['F'].append((0, 0, 0, 0, 1))
            self.mat['C'].append(-(1 - alpha_f) * state)
            # self.mat['C'].append(0)
            # print 'P toggle'
            # print state
        else:
            self.mat['F'].append((0, 0, 0, 0, 1))
            self.mat['C'].append(-state)
            # print 'State maintain'
            # print state


class PressureSource(LPNBlock):
    """
    Pressure source (two terminal : analogous to a battery)
    """
    def __init__(self, Pfunction, Pref, connecting_block_list=None, name="NoNamePressureSource", flow_directions=None):
        LPNBlock.__init__(self, connecting_block_list, name=name, flow_directions=flow_directions)
        self.neq = 2
        self.type = "PressureSource"
        self.Pfunction = Pfunction
        self.Pref = Pref

    def update_time(self, args):
        t = args['Time']
        self.mat['C'] = [-1.0 * self.Pref, -1.0 * self.Pfunction(t)]

    def update_constant(self):
        self.mat['F'] = [(1., 0, 0, 0), (0, 0, 1., 0)]


class UnsteadyResistanceWithDistalPressure_special(LPNBlock):
    def __init__(self, Rfunc, Pref_func, connecting_block_list=None,
                 name="NoNameUnsteadyResistanceWithDistalPressure_special", flow_directions=None):
        LPNBlock.__init__(self, connecting_block_list, name=name, flow_directions=flow_directions)
        self.type = "UnsteadyResistanceWithDistalPressure_special"
        self.neq = 1
        self.n_connect = 1
        self.Rfunc = Rfunc
        self.Pref_func = Pref_func

    def update_time(self, args):
        # For resistors, the ordering is : (P_in,Q_in)
        t = args['Time']

        # currently here 7/19/20: since this resistor element with distal pressure is supposed to be a terminal element, so that i dont have to connect a pressureRef block downstream, how many unknowns and equations should this block have? I think its 2 unknowns (P_in,Q_in), but how many equations do i need?
        self.mat['F'] = [(1., -1.0 * self.Rfunc(t))]
        self.mat['C'] = [-1.0 * self.Pref_func(t)]


# -- Resistance
class UnsteadyResistance_special(LPNBlock):
    def __init__(self, Rfunc, connecting_block_list=None, name="NoNameUnsteadyResistance_special",
                 flow_directions=None):
        LPNBlock.__init__(self, connecting_block_list, name=name, flow_directions=flow_directions)
        self.type = "UnsteadyResistance_special"
        self.Rfunc = Rfunc

    def update_time(self, args):
        # For resistors, the ordering is : (P_in,Q_in,P_out,Q_out)
        t = args['Time']
        self.mat['F'] = [(1., -1.0 * self.Rfunc(t), -1., 0), (0, 1., 0, -1.)]


class UnsteadyFlowRef_special(LPNBlock):
    def __init__(self, Qfunc, connecting_block_list=None, name="NoNameUnsteadyFlowRef_special", flow_directions=None):
        LPNBlock.__init__(self, connecting_block_list, name=name, flow_directions=flow_directions)
        self.type = "UnsteadyFlowRef_special"
        self.neq = 1
        self.Qfunc = Qfunc

    def update_time(self, args):
        t = args['Time']
        self.mat['C'] = [-1.0 * self.Qfunc(t)]

    def update_constant(self):
        self.mat['F'] = [(0, 1.)]


class Junction_special(LPNBlock):
    def __init__(self, temp_parameter, connecting_block_list=None, name="NoNameJunction_special", flow_directions=None):
        LPNBlock.__init__(self, connecting_block_list, name=name, flow_directions=flow_directions)
        self.type = "Junction_special"
        self.neq = self.num_connections
        self.temp_parameter = temp_parameter

        print("------------------------------------------")
        print("Junction_special - click to continue")
        print("temp_parameter = ", self.temp_parameter)
        print("------------------------------------------")

    def add_connecting_block(self, block, direction):
        self.connecting_block_list.append(block)
        self.num_connections = len(self.connecting_block_list)
        self.neq = self.num_connections
        self.flow_directions.append(direction)
        # print self.name

    def update_constant(self):
        # Number of variables per tuple = 2*num_connections
        # Number of equations = num_connections-1 Pressure equations, 1 flow equation
        # Format : P1,Q1,P2,Q2,P3,Q3, .., Pn,Qm

        self.mat['E'] = [(0,) * (2 * self.num_connections)] * (self.num_connections)

        self.mat['F'] = [(1.,) + (0,) * (2 * i + 1) + (-1,) + (0,) * (2 * self.num_connections - 2 * i - 3) for i in
                         range(self.num_connections - 1)]

        tmp = (0,)
        for d in self.flow_directions[:-1]:
            tmp += (d,)
            tmp += (0,)

        tmp += (self.flow_directions[-1],)
        self.mat['F'].append(tmp)
