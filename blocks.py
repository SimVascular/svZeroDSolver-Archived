# coding=utf-8

class LPNVariable:
    def __init__(self,value,units,name="NoName",vtype='ArbitraryVariable'):
        self.type = vtype
        self.value = value
        # Two generic units accepted : SI, cgs. Conversion for special values applied
        self.units = units
        self.name = name


class PressureVariable(LPNVariable):

    def __init__(self,value,units='cgs',name='NoNamePressure'):
        LPNVariable.__init__(self,value=value,units=units,name=name,vtype='Pressure')

    def convert_to_cgs(self):
        if self.units == 'cgs':
            print("Variable: "+self.name+" already at cgs")
        elif self.units == 'SI' :
            self.value = self.value*1.0E5
            self.units = 'cgs'
        elif self.units == 'mmHg':
            self.value = self.value*0.001333224
            self.units = 'cgs'
        else :
            raise Exception("Units "+self.units+" not recognized")


    def convert_to_mmHg(self):
        if self.units == 'cgs':
            self.value = self.value*750.06
            self.units = 'mmHg'
        elif self.units == 'SI' :
            self.value = self.value*7.50*1E-3
            self.units = 'mmHg'
        elif self.units == 'mmHg':
            print("Variable: "+self.name+" already at mmHg")
        else :
            raise Exception("Units "+self.units+" not recognized")


class FlowVariable(LPNVariable):

    def __init__(self,value,units='cgs',name='NoNameFlow'):
        LPNVariable.__init__(self,value=value,units=units,name=name,vtype='Flow')

    def convert_to_cgs(self):
        if self.units == 'cgs':
            print("Variable: "+self.name+" already at cgs")
        elif self.units == 'SI' :
            self.value = self.value*1.0E-6
            self.units = 'cgs'
        elif self.units == 'Lpm' :  # litres per minute
            self.value = self.value*16.6667
            self.units = 'cgs'
        else :
            raise Exception("Units "+self.units+" not recognized")


    def convert_to_Lpm(self):
        if self.units == 'cgs':
            self.value = self.value/16.6667
            self.units = 'Lpm'
        elif self.units == 'SI' :
            self.value = self.value/(16.6667*1.0E-6)
            self.units = 'Lpm'
        elif self.units == 'Lpm':
            print("Variable: "+self.name+" already at Lpm")
        else :
            raise Exception("Units "+self.units+" not recognized")


#==========================================================
#                        LPN WIRES
#==========================================================

# Wires connect circuit elements and junctions
# They can only posses a single pressure and flow value (system variables)
# They can also only possess one element(or junction) at each end
class wire:
    def __init__(self,connecting_elements,Pval=0,Qval=0,name="NoNameWire",P_units='cgs',Q_units='cgs'):
        self.name=name
        self.type='Wire'
        self.P = PressureVariable(value=Pval,units=P_units,name=name+"_P")
        self.Q = FlowVariable(value=Qval,units=Q_units,name=name+"_Q")
        if len(connecting_elements) >2 :
            raise Exception('Wire cannot connect to more than two elements at a time. Use a junction LPN block')
        if type(connecting_elements)!=tuple :
            raise Exception('Connecting elements to wire should be passed as a 2-tuple')
        self.connecting_elements = connecting_elements

        self.LPN_solution_ids = [None]*2 # i think there are 2 items in the list here because each wire should have to solution variables (P and Q) # LPN_solution_ids is a list that contains the index at which this wire's P and Q solutions are stored in the global vector of solution variables/unknown

    def initialize_PQ(self,y,P,Q):
        y[self.LPN_solution_ids[0]] = P
        y[self.LPN_solution_ids[0]] = Q # I think this should be [1]] = Q, but also this initialize_PQ function is never used so it doesn't matter


#==========================================================
#                        LPN BLOCKS
#==========================================================


# -- GENERAL LPN BLOCK

class LPNBlock:
    def __init__(self,connecting_block_list=None,name = "NoName",flow_directions=None):
        if connecting_block_list == None:
            connecting_block_list = []
        self.connecting_block_list = connecting_block_list
        self.num_connections = len(connecting_block_list)
        self.name = name
        self.neq = 2
        self.type="ArbitraryBlock"
        self.num_block_vars = 0
        self.connecting_wires_list = []
        if flow_directions == None :
            self.flow_directions = [] # -1 : Inflow to block, +1 outflow from block
        else :
            self.flow_directions = flow_directions
        self.LPN_solution_ids = [] # LPN_solution_ids for LPNBlock contains the solution IDs for the LPN block's internal solution variables; basically, LPN_solution_ids stores the index at which the LPNBlock's internal solution variables are stored in the global vector of solution variables/unknowns
        self.emxcoe = []
        self.fmxcoe = []
        self.cveccoe = []

        # Tangent matrix coes
        self.demxcoe = [] # sum_k[ydot_k(d(E_ik)/dy_j)]
        self.dfmxcoe = [] # sum_k[y_k(d(F_ik)/dy_j)]
        self.dcmxcoe = [] # (d(C_i)/dy_j)

        self.global_col_id = [] # a list of the indices at which this LPNBlock's associated solution variables (Pin, Qin, Pout, Qout, and internal solutions) are stored in the global vector of solution variables/unknowns
        self.global_row_id = []

        self.is_nonlinear = False

    def check_block_consistency(self):
        return

    def add_connecting_block(self,block,direction):
        # Direction = +1 if flow sent to block
        #            = -1 if flow recvd from block
        self.connecting_block_list.append(block)
        self.num_connections = len(self.connecting_block_list)
        self.flow_directions.append(direction)
        # print self.name


    def add_connecting_wire(self,new_wire):
        self.connecting_wires_list.append(new_wire)
        # self.flow_directions.append(direction)

    def update_constant(self):
        pass

    def update_time(self, args):
        pass

    def update_solution(self, args):
        pass

    def eqids(self,wire_dict,local_eq):
        # EqID returns variable's location in solution vector

        nwirevars = self.num_connections*2 # num_connections is multipled by 2 because each wire has 2 soltns (P and Q)
        if local_eq < nwirevars :
            vtype = local_eq%2 # 0 --> P, 1 --> Q
            wnum = int(local_eq/2)

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
        else : # this section will return the index at which the LPNBlock's  INTERNAL SOLUTION VARIABLES are stored in the global vector of solution unknowns/variables (i.e. I think RCR and OpenLoopCoronaryBlock have internal solution variables; these internal solution variables arent the P_in, Q_in, P_out, Q_out that correspond to the solutions on the attached wires, they are the solutions that are internal to the LPNBlock itself)
            vnum = local_eq - nwirevars
            return self.LPN_solution_ids[vnum]


# -- JUNCTION
# Junction points between LPN blocks with specified directions of flow
class Junction(LPNBlock):
    def __init__(self,connecting_block_list=None,name="NoNameJunction",flow_directions=None):
        LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)
        self.type = "Junction"
        self.neq = self.num_connections # number of equations = num of blocks that connect to this junction, where the equations are 1) mass conservation 2) inlet pressures = outlet pressures


    def add_connecting_block(self,block,direction):
        self.connecting_block_list.append(block)
        self.num_connections = len(self.connecting_block_list)
        self.neq = self.num_connections
        self.flow_directions.append(direction)
        # print self.name

    def update_constant(self):

        # Number of variables per tuple = 2*num_connections
        # Number of equations = num_connections-1 Pressure equations, 1 flow equation
        # Format : P1,Q1,P2,Q2,P3,Q3, .., Pn,Qm

        # self.emxcoe = [(0,)*(2*self.num_connections)]*(self.num_connections)

        self.fmxcoe = [ (1.,)+(0,)*(2*i+1) + (-1,) + (0,)*(2*self.num_connections-2*i-3) for i in range(self.num_connections-1) ]

        tmp = (0,)
        for d in self.flow_directions[:-1]:
            tmp+=(d,)
            tmp+=(0,)

        tmp += (self.flow_directions[-1],)
        self.fmxcoe.append(tmp)
        # self.cveccoe = [0]*self.num_connections
        #
        # self.demxcoe = [(0,)*(2*self.num_connections)]*(self.num_connections)
        # self.dfmxcoe = [(0,)*(2*self.num_connections)]*(self.num_connections)
        # self.dcmxcoe = [(0,)*(2*self.num_connections)]*(self.num_connections)

# -- JUNCTION WITH PRESSURE LOSS
# Junction points between LPN blocks with specified directions of flow
class Junction_with_PressureLoss(LPNBlock):
    # Reference: Chnafa et al. Improved reduced-order modelling of cerebrovascular flow distribution by accounting for arterial bifurcation pressure drops. Journal of Biomechanics. 2017.
    # NOTES:    - this junction block with pressure losses currently works only for junctions with only a single inlet vessel, but it can have multiple outlet vessels
    #
    # TODOS:    - enable junctions with multiple inlets

    def __init__(self, fluid_density, segment_vectors, segment_radii, connecting_block_list = None, name = "NoNameJunction_with_PressureLoss", flow_directions = None):
        LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)
        self.type = "Junction_with_PressureLoss"
        self.neq = self.num_connections
        self.fluid_density = fluid_density # {segment number : fluid density for the vessel segment}
        self.segment_vectors = segment_vectors # {segment number : the (tangent) vector for the vessel segment (np array form: [x, y, z])}
        self.segment_radii = segment_radii # {segment number : segment's radius}

        # define the resistance corresponding to the pressure loss across the junction
        eps = 1.0e-7 # define a small value to add to the flow rate to avoid dividing by zero (when the flow rate is zero)

        # NEED TO CHECK THAT THE R_LOSS FUNCTION WAS DERIVED CORRECTLY IN MY NOTES - DONE, it is correct - 12/16/19

        # R_loss calculation is correct - I checked this - 12/17/19
        R_PressureLoss_func = lambda fluid_density, Q_dat, r_dat, Q_i, r_i, theta_i: (0.5*fluid_density*(Q_dat**2.0)/((np.pi**2.0)*(r_dat**4.0)))*((1.0/(Q_i + eps)) + (Q_i*((r_dat/r_i)**4.0)/((Q_dat + eps)**2.0)) - (2.0*((r_dat/r_i)**2.0)*np.cos(theta_i)/(Q_dat + eps)))
        self.R_PressureLoss_func = R_PressureLoss_func

    def add_connecting_block(self, block,direction):
        self.connecting_block_list.append(block)
        self.num_connections = len(self.connecting_block_list)
        self.neq = self.num_connections
        self.flow_directions.append(direction)

    def dRodQin(self, fluid_density, Q_dat, r_dat, Q_i, r_i, theta_i, R_PressureLoss_func): # this function calculates the derivative of the R_PressureLoss_func wrt Qin

        # Input  : R_PressureLoss_func(fluid_density, Q_dat, r_dat, Q_i, r_i, theta_i)
        # Return : d/dQ_inlet of (R_PressureLoss_func) using central difference

        eps_dq = 1.0e-3
        delta_Qin = eps_dq*Q_dat
        if abs(delta_Qin) < 1.0e-5:
            delta_Qin = 1.0e-5
        # derivative calculation is correct - I checked this - 12/17/19
        return (R_PressureLoss_func(fluid_density, Q_dat + delta_Qin, r_dat, Q_i, r_i, theta_i) - R_PressureLoss_func(fluid_density, Q_dat - delta_Qin, r_dat, Q_i, r_i, theta_i))/(2.0*delta_Qin)

    def dRodQi(self, fluid_density, Q_dat, r_dat, Q_i, r_i, theta_i, R_PressureLoss_func): # this function calculates the derivative of the R_PressureLoss_func wrt Qi

        # Input  : R_PressureLoss_func(fluid_density, Q_dat, r_dat, Q_i, r_i, theta_i)
        # Return : d/dQ_i of (R_PressureLoss_func) using central difference

        eps_dq = 1.0e-3
        delta_Qi = eps_dq*Q_i
        if abs(delta_Qi) < 1.0e-5:
            delta_Qi = 1.0e-5
        # derivative calculation is correct - I checked this - 12/17/19
        return (R_PressureLoss_func(fluid_density, Q_dat, r_dat, Q_i + delta_Qi, r_i, theta_i) - R_PressureLoss_func(fluid_density, Q_dat, r_dat, Q_i - delta_Qi, r_i, theta_i))/(2.0*delta_Qi)

    def calc_junction_angle(self, inlet_vessel_vector, outlet_vessel_vector): # calculate the angle between the parent vessel and the daughter vessel
        junction_angle = np.arccos(np.dot(inlet_vessel_vector, outlet_vessel_vector)/(np.linalg.norm(inlet_vessel_vector)*np.linalg.norm(outlet_vessel_vector))) # in radians
        # print "junction_angle = ", junction_angle
        return junction_angle

    def update_solution(self, args):

        # Number of variables per tuple = 2*num_connections
        # Number of equations = num_connections-1 Pressure equations, 1 flow equation
        # Format : P1,Q1,P2,Q2,P3,Q3, .., Pn,Qn

        curr_y = args['Solution'] # curr_y is the solution vector at the current time, t; curr_y is a 1D array, where each entry represents each wire, so entry i holds the solution for wire i
        wire_dict = args['Wire dictionary']
        # print "wire_dict[self.connecting_wires_list[0]].LPN_solution_ids[1] = ", wire_dict[self.connecting_wires_list[0]].LPN_solution_ids[1]
        # print "wire_dict[self.connecting_wires_list[0]] = ", wire_dict[self.connecting_wires_list[0]]
        # print "self.connecting_wires_list[0] = ", self.connecting_wires_list[0]
        # print "self.connecting_wires_list = ", self.connecting_wires_list
        # print "self.flow_directions = ", self.flow_directions # use flow_directions to find out if flow is an inlet or outlet flow to the junction
        # # print "LPN_solution_ids[1] = ", LPN_solution_ids[1]
        # raw_input("yup\n")
        # notes:     -    self.connecting_wires_list = a list that contains the names of the wires connected to this junction block
        #             -    self.connecting_wires_list[0] = the name of the wire in index 0, ie R0_J0
        #            -    wire_dict[self.connecting_wires_list[0]] houses the wire object with the name, self.connecting_wires_list[0]
        #             -     wire_dict[self.connecting_wires_list[0]].LPN_solution_ids[1] gives you the index in the solution array, curr_y, where that index pertains to the wire of interest, wire_dict[self.connecting_wires_list[0]]
        #                      - LPN_solution_ids[i] where i = 0 (pressure solution) or i = 1 (flow rate solution)
        # q = curr_y[wire_dict[self.connecting_wires_list[0]].LPN_solution_ids[1]]
        # print "q = ", q

        # obtain the inlet and outlet flows to the junction
        junction_inflows = [] # a list of the inlet flow rates where the length of junction_inflows is = number of connecting blocks; but obviously not all connecting blocks to this junction block are inlet blocks but some are outlet blocks instead; thus if junction_inflows[i] == "not", then this means that index i represents an outlet block instead of an inlet block, where we know if something is an inlet or outlet based on flow_directions
        junction_outflows = [] # a list of the outlet flow rates where the length of junction_outflows is = number of connecting blocks; but obviously not all connecting blocks to this junction block are outlet blocks but some are arelet blocks instead; thus if junction_outflows[i] == "not", then this means that index i represents an inlet block instead of an outlet block
        Q_dat = 0.0 # initialize Q_dat, the total inlet flow rate
        inlet_segment_numbers = [] # a list of the segment numbers corresponding to the inlet vessels of the junction; has the same structure as junction_inflows, where the "not" entry in index i of inlet_segment_numbers means that index i represents an outlet block instead of an inlet
        outlet_segment_numbers = [] # a list of the segment numbers corresponding to the outlet vessels of the junction; has the same structure as inlet_segment_numbers
        for i in range(len(self.flow_directions)): # recall that self.flow_directions and self.connecting_wires_list are both lists where the item in index i of both lists corresponds to the same wire (which is effectively the same block that connects to this junction block)
            if self.flow_directions[i] == -1: # this is an inflow to the junction block
                inflow_i = curr_y[wire_dict[self.connecting_wires_list[i]].LPN_solution_ids[1]] # flow rate coming into the junction
                junction_inflows.append(inflow_i)
                Q_dat += inflow_i
                junction_outflows.append("not")
                # determine the inlet segment number (for use in calculating the angle between parent/daughter segments to calculate R_loss)
                # NOTE: THIS METHOD to determine the segement number for the inlet CURRENTLY ASSUMES THAT THERE IS ONLY ONE INLET VESSEL TO THE JUNCTION (THIS WORKS ONLY FOR DIVERGING FLOWS AT JUNCTIONS (ARTERIAL SYSTEMS), NOT CONVERGING FLOWS (CONFLUENCES; VENUOUS SYSTEMS))
                inflow_wire_name = self.connecting_wires_list[i]
                inlet_segment_numbers.append(int(inflow_wire_name[-1])) # inflow_wire_name will be something like "R0_J0", so the last value in the string is the segment number of the parent (inlet) vessel
                outlet_segment_numbers.append("not")
            elif self.flow_directions[i] == 1: # this is an outflow from the junction block
                junction_inflows.append("not")
                outflow_i = curr_y[wire_dict[self.connecting_wires_list[i]].LPN_solution_ids[1]] # flow rate coming out of the junction
                junction_outflows.append(outflow_i)
                # determine the outlet segment number (for use in calculating the angle between parent/daughter segments to calculate R_loss)
                inlet_segment_numbers.append("not")
                outflow_wire_name = self.connecting_wires_list[i]
                outlet_segment_numbers.append(int(outflow_wire_name[-1])) # outflow_wire_name will be something like "J0_R1", so the last value in the string is the segment number of the daughter (outlet) vessel
            else:
                raw_input("Error. self.flow_directions contains an entry that is not 1 or -1. Quit the program and fix the error now.")

        # note:
        #         recall that for now, i am developing code to calculate pressure losses in junction blocks assuming that each junction has only 1 inlet, but can have more than 1 outlet
        #
        #         steps:
        #         1. from flow_direction, figure out the index that corresponds to the inlet (the inflow) - use python's built-in function, LIST.index(element), to find the index that holds the entry "element" in LIST
        #         2. now I know which index corresponds to the inflow; with this info, then in E, F, C dE, dF, dC, I know which index corresponds to the inflow and then I can place a "1" in this position, where the "1" is the coefficient for P_inlet in the pressure eqn for the junction, P_inlet - P_i - R_loss,i*Q_i = 0
        #         3. initialize F as an empty list
        #         4. loop through flow_directions and at the beginning of the loop, create a new tuple; then for each outlet, append zeros and the coefficient to this tuple; then at the end, append this tuple to F; note: the index in flow_direction that corresponds to the current outlet that you are dealing with tells you the index in F at which to place the "-1" and "R_loss" coefficients - before and after that index position, place zeros in F
        #         5. do the same for the derivative matrices, dE, dF, etc

        # print "\njunction_inflows = ", junction_inflows
        # print "\njunction_outflows = ", junction_outflows
        # print "\nQ_dat = ", Q_dat
        # print "\ninlet_segment_numbers = ", inlet_segment_numbers
        # print "\noutlet_segment_numbers = ", outlet_segment_numbers
        # raw_input("\ninside here 1 ---------------------------------------------------------")
         # CURRENTLY CHECKED UP TO HERE - 12/16/19

         # question: does aekaansh's code initialize all the flow rates to be zero initially? - look into this

        # construct the E, F, C, dE/dy, dF/dy, dC/dy matrices
        self.emxcoe = [(0.,)*(2*self.num_connections)]*(self.num_connections)
        self.demxcoe = [(0.,)*(2*self.num_connections)]*(self.num_connections)

        inlet_index = self.flow_directions.index(-1) # find the index in flow_directions that corresponds to the inlet
        # print "\ninlet_index = ", inlet_index
        r_dat = self.segment_radii[inlet_segment_numbers[inlet_index]] # inlet vessel radius
        # print "\nr_dat = ", r_dat
        inlet_vessel_vector = self.segment_vectors[inlet_segment_numbers[inlet_index]]
        # print "\ninlet_vessel_vector = ", inlet_vessel_vector
        self.fmxcoe = []
        self.dfmxcoe = []
        for i in range(len(self.flow_directions)): # note: len(self.flow_directions) = num_connections = number of blocks attached to this block = number of equations for this block
            # index i is like the index of the equation, except not really because there is also the mass conservation equation
            if i != inlet_index: # the coefficients in F and dF for the inlet vessel will be taken into account inherently in the code below
                # print "\noutlet_segment_numbers[i] = ", outlet_segment_numbers[i]
                eqn_temp_line_f = [0, 0]*self.num_connections
                # print "\neqn_temp_line_f = ", eqn_temp_line_f
                eqn_temp_line_f[inlet_index*2] = 1.0 # apply the correct coefficient for the inlet vessel's pressure to the F equation
                Q_i = junction_outflows[i]
                # print "\nQ_i = ", Q_i
                r_i = self.segment_radii[outlet_segment_numbers[i]]
                # print "\nr_i = ", r_i
                outlet_vessel_vector = self.segment_vectors[outlet_segment_numbers[i]]
                # print "\noutlet_vessel_vector = ", outlet_vessel_vector
                theta_i = self.calc_junction_angle(inlet_vessel_vector, outlet_vessel_vector)
                # print "\ntheta_i = ", theta_i
                rho_i = self.fluid_density[outlet_segment_numbers[i]]
                # print "\nrho_i = ", rho_i
                R = self.R_PressureLoss_func(rho_i, Q_dat, r_dat, Q_i, r_i, theta_i)
                # print "\nR =", R
                # print "\n-1.0*R = ", -1.0*R
                eqn_temp_line_f[i*2:i*2 + 2] = [-1.0, float(-1.0*R)] # apply the correct coefficient for the current outlet vessel to the F equation
                # print "\neqn_temp_line_f = ", eqn_temp_line_f
                self.fmxcoe.append(tuple(eqn_temp_line_f))

                eqn_temp_line_df = [0, 0]*self.num_connections
                # print "\neqn_temp_line_df = ", eqn_temp_line_df # just finished adding print statement here - 12/16/19
                dRodQin = self.dRodQin(rho_i, Q_dat, r_dat, Q_i, r_i, theta_i, self.R_PressureLoss_func)
                # print "\ndRodQin  =", dRodQin
                eqn_temp_line_df[inlet_index*2 + 1] = float(-1.0*dRodQin*Q_i)
                dRodQi = self.dRodQi(rho_i, Q_dat, r_dat, Q_i, r_i, theta_i, self.R_PressureLoss_func)
                # print "\ndRodQi  =", dRodQi
                eqn_temp_line_df[i*2 + 1] = float(-1.0*dRodQi*Q_i) # apply the correct coefficient for the current outlet vessel to the F equation
                # print "\n-1.0*dRodQin*Q_i = ", -1.0*dRodQin*Q_i
                # print "\n-1.0*dRodQi*Q_i = ", -1.0*dRodQi*Q_i
                # print "\neqn_temp_line_df = ", eqn_temp_line_df # just finished adding print statement here - 12/16/19
                self.dfmxcoe.append(tuple(eqn_temp_line_df))
                # raw_input("\ninside here 2 ---------------------------------------------------------")
        tmp = () # this is used to construct the mass/flow conservation equation
        for d in self.flow_directions:
            tmp += (0, d,)
        self.fmxcoe.append(tmp)
        self.dfmxcoe.append((0, 0, )*self.num_connections) # this line represents the coefficients for the mass conservation eqn in the junction # good
        # print "\nfmxcoe = ", self.fmxcoe
        # print "\ndfmxcoe = ", self.dfmxcoe
        # raw_input("\ninside here 3 ---------------------------------------------------------")
        self.cveccoe = [0.]*self.num_connections
        self.dcmxcoe = [(0.,)*(2*self.num_connections)]*(self.num_connections)

# -- Resistance
class Resistance(LPNBlock):
    def __init__(self,R,connecting_block_list=None,name="NoNameResistance",flow_directions=None):

        LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)
        self.type = "Resistance"
        self.R = R

    def check_block_consistency(self):
        if len(connecting_block_list) != 2:
            raise Exception("Resistance block can be connected only to two elements")

    def update_constant(self):

        # For resistors, the ordering is : (P_in,Q_in,P_out,Q_out)

        # self.emxcoe = [(0,)*4]*2
        self.fmxcoe = [(1.,-1.*self.R,-1.,0),(0,1.,0,-1.)]
        # self.cveccoe = [0]*2
        #
        # self.demxcoe = [(0,)*4]*2
        # self.dfmxcoe = [(0,)*4]*2
        # self.dcmxcoe = [(0,)*4]*2

# -- Flow dependent Resistance : delta_P = q*Rfunc(t,q)
class FlowDepResistance(LPNBlock):
    def __init__(self,Rfunc,connecting_block_list=None,name="NoNameFlowDepResistance",flow_directions=None):

        LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)
        self.type = "FlowDepResistance"
        self.Rfunc = Rfunc

    def check_block_consistency(self):
        if len(connecting_block_list) != 2:
            raise Exception("FlowDepResistance block can be connected only to two elements")

    def linearize_func(self,t,q,func,eps_dq=1e-3): # 3/7/21: this linearize function isnt really linearizing anything; it is just computing the derivative of "func" with respect to "q"

        # Input  : func(t,q), typically t = t_n+1 ie t_n + dt
        # Return : d/dq (func) at (t,q) using central differences

        dq = eps_dq*q # the step size, dq, for the central finite difference formula is arbitrary here, according to Aekaansh; here, he chooses to scale the step size by the flow rate, q
        if abs(dq) < 1e-5:
            dq = 1e-5
        return (func(t,q+dq)-func(t,q-dq))/(2.*dq)

    def update_solution(self, args):
        # delta_P = q*Rfunc(t,q) , so linearization yields:
        # delta_P(n+1) - q(n+1)*[Rfunc(t(n+1),q(n))+q(n)*dR] + q(n)^2*dR = 0

        t = args['Time']
        curr_y = args['Solution']
        # dt = args['Time step']
        # rho = args['rho']
        # alpha_f = 1/(1.0+rho)
        wire_dict = args['Wire dictionary']
        q = curr_y[wire_dict[self.connecting_wires_list[0]].LPN_solution_ids[1]]

        # pi = curr_y[wire_dict[self.connecting_wires_list[0]].LPN_solution_ids[0]]
        # po = curr_y[wire_dict[self.connecting_wires_list[1]].LPN_solution_ids[0]]

        dR = self.linearize_func(t,q,self.Rfunc)
        R = self.Rfunc(t,q)

        # c_lin = q*q*dR
        # self.emxcoe = [(0,)*4]*2
        # self.fmxcoe = [(1.,-1.*(self.Rfunc(t+dt*alpha_f,q)+q*dR),-1.,0),(0,1.,0,-1.)]
        # self.cveccoe = [c_lin,0]

        # self.fmxcoe = [(1./(R),-1.,-1./(R),0),(0,1.,0,-1.)]
        self.fmxcoe = [(1.,-1.*(R),-1.,0),(0,1.,0,-1.)]

        # self.dfmxcoe = [(0,(pi-po)*(-dR/(R*R)),0,0),(0,)*4]
        self.dfmxcoe = [(0,-dR*q,0,0),(0,)*4] # WHY IS dR MULTIPLED BY q HERE? SHOULDNT IT JUST BE dR? - ask aekaansh, 12/4/19: after I ask him, then I can create the dE, dF, dC matrices for the Junction_with_PressureLoss # 3/7/21: see derivation of the dfmxcoe matrix in page 1 of the physical composition notebook titled "Research Notebook #3"

    # def update_constant(self):
    #     self.emxcoe = [(0,)*4]*2 # "mxcoe" = "MatriX COEfficient"
    #     self.cveccoe = [0,0]
    #
    #     self.demxcoe = [(0,)*4]*2
    #     self.dcmxcoe = [(0,)*4]*2


"""
Stenosis:
    equation:   delta_P = ( K_t * rho / ( 2 * (A_0)**2 ) ) * ( ( A_0 / A_s ) - 1 )**2 * Q * abs(Q)
                        =               stenosis_coefficient                          * Q * abs(Q)

    source: Mirramezani, M., Shadden, S.C. A distributed lumped parameter model of blood flow. Annals of Biomedical Engineering. 2020.
"""
class StenosisBlock(LPNBlock):
    def __init__(self,R,stenosis_coefficient,connecting_block_list=None,name="NoNameRCL",flow_directions=None):

        LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)
        self.type = "Stenosis"
        self.R = R # poiseuille resistance value = 8 * mu * L / (pi * r**4)
        self.stenosis_coefficient = stenosis_coefficient
        self.is_nonlinear = True

    def check_block_consistency(self):
        if len(connecting_block_list) != 2:
            raise Exception("StenosisBlock can be connected only to two elements")

    def update_solution(self, args):
        curr_y = args['Solution'] # the current solution for all unknowns in our 0D model
        wire_dict = args['Wire dictionary']
        Q_in =  curr_y[wire_dict[self.connecting_wires_list[0]].LPN_solution_ids[1]]
        self.fmxcoe   = [ (1.0, -1.0 * self.stenosis_coefficient * np.abs(Q_in) - self.R, -1.0, 0), (0, 1.0, 0, -1.0)]
        self.dfmxcoe  = [ (0, -1.0 * self.stenosis_coefficient * np.abs(Q_in), 0, 0), (0,)*4] # 3/7/21: see derivation in GoodNotes on 3/7/21 in the notebook titled "Research"

    def update_constant(self):
        self.emxcoe   = [ (0,0,0,0), (0,0,0,0)]
        self.cveccoe  = [0, 0]
        self.demxcoe  = [(0,)*4]*2
        self.dcmxcoe  = [(0,)*4]*2


# -- Unsteady Resistance : delta_P = q*Rfunc(t)
class UnsteadyResistance(LPNBlock):
    def __init__(self,Rfunc,connecting_block_list=None,name="NoNameUnsteadyResistance",flow_directions=None):

        LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)
        self.type = "UnsteadyResistance"
        self.Rfunc = Rfunc

    def check_block_consistency(self):
        if len(connecting_block_list) != 2:
            raise Exception("UnsteadyResistance block can be connected only to two elements")

    def update_time(self, args):
        # For resistors, the ordering is : (P_in,Q_in,P_out,Q_out)
        t = args['Time']

        self.fmxcoe = [(1.,-1.0*self.Rfunc(t),-1.,0),(0,1.,0,-1.)]

    # def update_constant(self):
    #     self.emxcoe = [(0,)*4]*2
    #     self.cveccoe = [0]*2
    #
    #     self.demxcoe = [(0,)*4]*2
    #     self.dfmxcoe = [(0,)*4]*2
    #     self.dcmxcoe = [(0,)*4]*2

class UnsteadyResistanceWithDistalPressure(LPNBlock):
    def __init__(self,Rfunc,Pref_func,connecting_block_list=None,name="NoNameUnsteadyResistanceWithDistalPressure",flow_directions=None):

        LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)
        self.type = "UnsteadyResistanceWithDistalPressure"
        self.neq = 1
        self.Rfunc = Rfunc
        self.Pref_func = Pref_func

    def check_block_consistency(self):
        if len(connecting_block_list) != 1:
            raise Exception("UnsteadyResistanceWithDistalPressure block can be connected only to two elements")

    def update_time(self, args):
        # For resistors, the ordering is : (P_in,Q_in)
        t = args['Time']

        # currently here 7/19/20: since this resistor element with distal pressure is supposed to be a terminal element, so that i dont have to connect a pressureRef block downstream, how many unknowns and equations should this block have? I think its 2 unknowns (P_in,Q_in), but how many equations do i need?
        self.fmxcoe = [(1.,-1.0*self.Rfunc(t))]
        self.cveccoe = [-1.0*self.Pref_func(t)]

    def update_constant(self):
        self.emxcoe = [(0,)*2]

        self.demxcoe = [(0,)*2]
        self.dfmxcoe = [(0,)*2]
        self.dcmxcoe = [(0,)*2]

# -- Pressure reference
class PressureRef(LPNBlock):
    def __init__(self,Pref,connecting_block_list=None,name="NoNamePressureRef",flow_directions=None):

        LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)
        self.type = "PressureRef"
        self.neq = 1
        self.Pref=Pref

    def check_block_consistency(self):
        if len(connecting_block_list) != 1:
            raise Exception("PressureRef block can be connected only to one element")

    def update_constant(self):
        # self.emxcoe = [(0,)]
        self.fmxcoe = [(1.,)]
        self.cveccoe = [-1.0*self.Pref]

        # self.demxcoe = [(0,)]
        # self.dfmxcoe = [(0,)]
        # self.dcmxcoe = [(0,)]

# -- Unsteady P reference
class UnsteadyPressureRef(LPNBlock):
    def __init__(self,Pfunc,connecting_block_list=None,name="NoNameUnsteadyPressureRef",flow_directions=None):

        LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)
        self.type = "UnsteadyPressureRef"
        self.neq = 1
        self.Pfunc=Pfunc

    def check_block_consistency(self):
        if len(connecting_block_list) != 1:
            raise Exception("UnsteadyPressureRef block can be connected only to one element")

    def update_time(self, args):
        t = args['Time']
        self.cveccoe = [-1.0*self.Pfunc(t)]

    def update_constant(self):
        self.emxcoe = [(0,0)]
        self.fmxcoe = [(1.,0.)]

        self.demxcoe = [(0,)]
        self.dfmxcoe = [(0,)]
        self.dcmxcoe = [(0,)]




# -- Flow reference
class UnsteadyFlowRef(LPNBlock):
    def __init__(self,Qfunc,connecting_block_list=None,name="NoNameUnsteadyFlowRef",flow_directions=None):

        LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)
        self.type = "UnsteadyFlowRef"
        self.neq = 1
        self.Qfunc=Qfunc

    def check_block_consistency(self):
        if len(connecting_block_list) != 1:
            raise Exception("UnsteadyFlowRef block can be connected only to one element")

    def update_time(self, args):
        t = args['Time']
        self.cveccoe = [-1.0*self.Qfunc(t)]

    def update_constant(self):
        # self.emxcoe = [(0,0)]
        self.fmxcoe = [(0,1.)]

        # self.demxcoe = [(0,)]
        # self.dfmxcoe = [(0,)]
        # self.dcmxcoe = [(0,)]

# -- Capacitance
class Capacitance(LPNBlock):
    def __init__(self,C,connecting_block_list=None,name="NoNameCapacitance",flow_directions=None):

        LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)
        self.type = "Capacitance"
        self.C = C

    def check_block_consistency(self):
        if len(connecting_block_list) != 2:
            raise Exception("Capacitance block can be connected only to two elements")


    def update_constant(self):
        self.emxcoe = [(1.0*self.C,0,-1.0*self.C,0),(0,0,0,0)]
        self.fmxcoe = [(0,-1.0,0,0),(0,1.,0,-1.)]
        # self.cveccoe = [0,0]

        # self.demxcoe = [(0,)*4]*2
        # self.dfmxcoe = [(0,)*4]*2
        # self.dcmxcoe = [(0,)*4]*2



# -- RCL - constant resistor, capacitor, inductor - vessel representation
# Formulation includes additional variable : internal pressure proximal to capacitance.
# Trading off speed with enforced uniqueness
class RCLBlock(LPNBlock):
    def __init__(self,R,C,L,connecting_block_list=None,name="NoNameRCL",flow_directions=None):

        LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)
        self.type = "RCL"
        self.neq = 3
        self.num_block_vars = 1
        self.R = R
        self.C = C
        self.L = L

    def check_block_consistency(self):
        if len(connecting_block_list) != 2:
            raise Exception("RCL block can be connected only to two elements")

    def update_constant(self):
        self.emxcoe = [(0,0,0,-self.L,0),      (0,0,0,0,-self.C),   (0,0,0,0,0)]
        self.fmxcoe = [(1.,-self.R,-1.,0,0), (0,1.,0,-1.,0),      (1.,-self.R,0,0,-1.)]
        self.cveccoe = [0,0,0]

        self.demxcoe = [(0,)*5]*3
        self.dfmxcoe = [(0,)*5]*3
        self.dcmxcoe = [(0,)*5]*3


# -- RC - constant resistor, capacitor - low inertia vessel
class RCBlock(LPNBlock):
    def __init__(self,R,C,connecting_block_list=None,name="NoNameRC",flow_directions=None):

        LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)
        self.type = "RC"
        self.R = R
        self.C = C

    def check_block_consistency(self):
        if len(connecting_block_list) != 2:
            raise Exception("RC block can be connected only to two elements")

    def update_constant(self):
        self.emxcoe = [(0,0,0,0),      (0,0,-self.C,0)]
        self.fmxcoe = [(1.0,-self.R,-1.0,0), (0,1.,0,-1.)]
        # self.cveccoe = [0,0]
        #
        # self.demxcoe = [(0,)*4]*2
        # self.dfmxcoe = [(0,)*4]*2
        # self.dcmxcoe = [(0,)*4]*2

# -- RL - constant resistor, inductor
class RLBlock(LPNBlock):
    def __init__(self,R,L,connecting_block_list=None,name="NoNameRL",flow_directions=None):

        LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)
        self.type = "RL"
        self.R = R
        self.L = L

    def check_block_consistency(self):
        if len(connecting_block_list) != 2:
            raise Exception("RL block can be connected only to two elements")

    def update_constant(self):

        # For this RL block, the ordering of solution unknowns is : (P_in, Q_in, P_out, Q_out)

        self.emxcoe = [(0,0,0,-self.L),      (0,0,0,0)]
        self.fmxcoe = [(1.0,-self.R,-1.0,0), (0,1.,0,-1.)]
        self.cveccoe = [0,0]

        self.demxcoe = [(0,)*4]*2
        self.dfmxcoe = [(0,)*4]*2
        self.dcmxcoe = [(0,)*4]*2

# -- RCR - constant RCR - outflow representation
# Formulation includes additional variable : internal pressure proximal to capacitance.
# Trading off speed with enforced uniqueness
class RCRBlock(LPNBlock):
    def __init__(self,Rp,C,Rd,connecting_block_list=None,name="NoNameRCR",flow_directions=None):

        LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)
        self.type = "RCR"
        self.neq = 3
        self.num_block_vars = 1
        self.Rp = Rp
        self.C = C
        self.Rd = Rd

    def check_block_consistency(self):
        if len(connecting_block_list) != 2:
            raise Exception("RCR block can be connected only to two elements")

    def update_constant(self):
        self.emxcoe = [(0,0,0,0,0),      (0,0,0,0,-self.C),   (0,0,0,0,0)]
        self.fmxcoe = [(1.0,-self.Rp,-1.0,-self.Rd,0), (0,1.,0,-1.,0),      (1.,-self.Rp,0,0,-1.)]
        # self.cveccoe = [0,0,0]
        #
        # self.demxcoe = [(0,)*5]*3
        # self.dfmxcoe = [(0,)*5]*3
        # self.dcmxcoe = [(0,)*5]*3

# -- Unsteady RCR - time-varying RCR values
# Formulation includes additional variable : internal pressure proximal to capacitance.
class UnsteadyRCRBlock(LPNBlock):
    def __init__(self,Rp_func,C_func,Rd_func,connecting_block_list=None,name="NoNameUnsteadyRCR",flow_directions=None):

        LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)
        self.type = "UnsteadyRCR"
        self.neq = 3
        self.num_block_vars = 1
        self.Rp_func = Rp_func
        self.C_func = C_func
        self.Rd_func = Rd_func

    def check_block_consistency(self):
        if len(connecting_block_list) != 2:
            raise Exception("Unsteady RCR block can be connected only to two elements")

    def update_time(self, args):
        t = args['Time']
        self.fmxcoe = [(1.0,-self.Rp_func(t),-1.0,-self.Rd_func(t),0),  (0,1.,0,-1.,0),   (1.,-self.Rp_func(t),0,0,-1.)]

    def update_constant(self):
        self.emxcoe = [(0,0,0,0,0),      (0,0,0,0,-self.C_func(t)),   (0,0,0,0,0)]
        self.cveccoe = [0,0,0]

        self.demxcoe = [(0,)*5]*3
        self.dfmxcoe = [(0,)*5]*3
        self.dcmxcoe = [(0,)*5]*3


# -- Unsteady RCR - time-varying RCR values
# Formulation includes additional variable : internal pressure proximal to capacitance.
class UnsteadyRCRBlockWithDistalPressure(LPNBlock):
    def __init__(self,Rp_func,C_func,Rd_func,Pref_func,connecting_block_list=None,name="NoNameUnsteadyRCRBlockWithDistalPressure",flow_directions=None):

        LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)
        self.type = "UnsteadyRCRBlockWithDistalPressure"
        self.neq = 2
        self.num_block_vars = 1
        self.Rp_func = Rp_func
        self.C_func = C_func
        self.Rd_func = Rd_func
        self.Pref_func = Pref_func

    def check_block_consistency(self):
        if len(connecting_block_list) != 1:
            raise Exception("UnsteadyRCRBlockWithDistalPressure block can be connected only to two elements")

    def update_time(self, args):
        # unknowns = [P_in, Q_in, internal_var (Pressure at the intersection of the Rp, Rd, and C elements)]

        t = args['Time']
        self.emxcoe = [(0,0,0),      (0,0,-1.0*self.Rd_func(t)*self.C_func(t))]
        self.fmxcoe = [(1.,-self.Rp_func(t),-1.), (0.0, self.Rd_func(t),-1.0)]
        self.cveccoe = [0,self.Pref_func(t)]

    # def update_constant(self):
        # self.demxcoe = [(0,)*3]*2
        # self.dfmxcoe = [(0,)*3]*2
        # self.dcmxcoe = [(0,)*3]*2

# -- Open loop coronary block - RCRCR - pressure imposed on the second capacitor
class OpenLoopCoronaryBlock(LPNBlock):
    "Publication reference: Kim, H. J. et al. Patient-specific modeling of blood flow and pressure in human coronary arteries. Annals of Biomedical Engineering 38, 3195–3209 (2010)."

    "open-loop coronary BC = RCRCR BC"

    def __init__(self,R1,C1,R2,C2,R3,dPvdt_f,cardiac_cycle_period, connecting_block_list=None,name="NoNameCoronary",flow_directions=None):

        LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)
        self.type = "Coronary"
        self.neq = 3
        self.num_block_vars = 1
        self.R1 = R1
        self.C1 = C1
        self.R2 = R2
        self.C2 = C2
        self.R3 = R3
        self.dPvdt_f = dPvdt_f
        self.cardiac_cycle_period = cardiac_cycle_period
# Provide a [Nx2] array for the time-derivative of intramyocardial pressure. dPvdt_f[:,0] is the time info. dPvdt_f[:,1] is the pressure deriv. values.

    def check_block_consistency(self):
        if len(connecting_block_list) != 2:
            raise Exception("Coronary block can be connected only to two elements")
# The Open Loop Coronary block needs to attach a downstream pressure reference block.

    def dPvdt(self,dPvdt_f,t):
        tt=dPvdt_f[:,0]
        dPvdt_t = dPvdt_f[:,1]
        ti, td = divmod(t, self.cardiac_cycle_period)
        dPv = np.interp(td, tt, dPvdt_t)
        return dPv
# dPvdt should return a time derivative of intramyocardial pressurea at time t. Here this example assumes a heart cycle of 1.0. User may implement modified funtion for different heart cycle.

    def update_time(self, args):

        t = args['Time']
        dPv = self.dPvdt(self.dPvdt_f,t)

        # For this open-loop coronary BC, the ordering of solution unknowns is : (P_in, Q_in, P_out, Q_out, Q_internal)
        # where Q_internal is the flow through the second resistor in the RCRCR BC
        # and Q_out is the flow the third resistor, and Q_in is the flow through the first resistor
        # and P_in is the pressure at the inlet of the first resistor and P_out is the the pressure ath the outlet of the third resistor
        # Note that the first capacitor is attached to ground but the second capacitor is attached to the time-varying  intramyocardial pressure
        self.cveccoe = [0,0,self.C2*dPv]

    def update_constant(self):
        self.emxcoe = [(0,0,0,0,0), (-self.C1,self.C1*self.R1,0,0,0),   (-self.C2,self.C2*self.R1,0,0,self.C2*self.R2)]
        self.fmxcoe = [(1.0,-self.R1,-1.0,-self.R3,-self.R2), (0,1.0,0,0,-1.0),      (0,0,0,-1.0,1.0)]

        # self.demxcoe = [(0,)*5]*3
        # self.dfmxcoe = [(0,)*5]*3
        # self.dcmxcoe = [(0,)*5]*3

class OpenLoopCoronaryWithDistalPressureBlock(LPNBlock):
    "Publication reference: Kim, H. J. et al. Patient-specific modeling of blood flow and pressure in human coronary arteries. Annals of Biomedical Engineering 38, 3195–3209 (2010)."

    "open-loop coronary BC = RCRCR BC"

    def __init__(self,R1,C1,R2,C2,R3,dPvdt_f,Pv,cardiac_cycle_period,connecting_block_list=None,name="NoNameCoronary",flow_directions=None):

        LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)
        self.type = "OpenLoopCoronaryWithDistalPressureBlock"
        self.neq = 2
        self.num_block_vars = 1
        self.R1 = R1
        self.C1 = C1
        self.R2 = R2
        self.C2 = C2
        self.R3 = R3
        self.dPvdt_f = dPvdt_f
        self.Pv = Pv
        self.cardiac_cycle_period = cardiac_cycle_period
# Provide a [Nx2] array for the time-derivative of intramyocardial pressure. dPvdt_f[:,0] is the time info. dPvdt_f[:,1] is the pressure deriv. values.

    def check_block_consistency(self):
        if len(connecting_block_list) != 1:
            raise Exception("OpenLoopCoronaryWithDistalPressureBlock can be connected only to one elements")
# The Open Loop Coronary block needs to attach a downstream pressure reference block.

    def dPvdt(self,dPvdt_f,t):
        tt=dPvdt_f[:,0]
        dPvdt_t = dPvdt_f[:,1]
        ti, td = divmod(t, self.cardiac_cycle_period)
        dPv = np.interp(td, tt, dPvdt_t)
        return dPv
# dPvdt should return a time derivative of intramyocardial pressurea at time t. Here this example assumes a heart cycle of 1.0. User may implement modified funtion for different heart cycle.

    def update_time(self,args):

        # unknowns = [P_in, Q_in, flow rate through the second resistor in the RCRCR block]

        t = args['Time']
        dPv = self.dPvdt(self.dPvdt_f,t)

        # For this open-loop coronary BC, the ordering of solution unknowns is : (P_in, Q_in, P_out, Q_out, Q_internal)
        # where Q_internal is the flow through the second resistor in the RCRCR BC
        # and Q_out is the flow the third resistor, and Q_in is the flow through the first resistor
        # and P_in is the pressure at the inlet of the first resistor and P_out is the the pressure ath the outlet of the third resistor
        # Note that the first capacitor is attached to ground but the second capacitor is attached to the time-varying  intramyocardial pressure
        self.cveccoe = [0, -1.0*self.R3*self.C2*dPv - self.Pv]

    def update_constant(self):
        self.emxcoe = [(-1.0*self.C1, self.R1*self.C1, 0), (self.R3*self.C2, -1.0*self.R3*self.C2*self.R1, -1.0*self.R3*self.C2*self.R2)]
        self.fmxcoe = [(0, 1.0, -1.0), (1.0, -self.R1, -1.0*(self.R3 + self.R2))]

        # self.demxcoe = [(0,)*3]*self.neq
        # self.dfmxcoe = [(0,)*3]*self.neq
        # self.dcmxcoe = [(0,)*3]*self.neq

###################################################################
class OpenLoopCoronaryWithDistalPressureBlock_v2(LPNBlock):
    "Publication reference: Kim, H. J. et al. Patient-specific modeling of blood flow and pressure in human coronary arteries. Annals of Biomedical Engineering 38, 3195–3209 (2010)."

    "open-loop coronary BC = RCRCR BC"

    def __init__(self, Ra, Ca, Ram, Cim, Rv, Pim, Pv, cardiac_cycle_period, connecting_block_list = None, name = "NoNameCoronary", flow_directions = None):

        LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)
        self.type = "OpenLoopCoronaryWithDistalPressureBlock_v2"
        self.neq = 2
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
        # import scipy
        P_tt = np.interp(td, tt, P_val)
        # P_tt = scipy.interpolate.CubicSpline(tt, P_val, bc_type = 'periodic')(td)
        # import matplotlib.pyplot as plt
        # plt.figure()
        # plt.plot(P[:, 0], P[:, 1], label = "P")
        # plt.show()
        return P_tt

    def check_block_consistency(self):
        if len(connecting_block_list) != 1:
            raise Exception("OpenLoopCoronaryWithDistalPressureBlock_v2 can be connected only to one elements")

    def update_time(self,args):

        # For this open-loop coronary BC, the ordering of solution unknowns is : (P_in, Q_in, V_im)
        # where V_im is the volume of the second capacitor, Cim
        # Q_in is the flow through the first resistor
        # and P_in is the pressure at the inlet of the first resistor
        ttt = args['Time']

        Pim_value = self.get_P_at_t(self.Pim, ttt)
        Pv_value = self.get_P_at_t(self.Pv, ttt)

        # Pim_func = lambda t: np.interp(t, self.Pim[:, 0], self.Pim[:, 1])
        # Pv_func = lambda t: self.Pv[:, 1][0] + (t - t)

        # Pim_value = Pim_func(ttt)
        # print("Pim_value = ", Pim_value)
        # Pv_value = Pv_func(ttt)
        # print("Pv_value = ", Pv_value)
        # input("")
        # last here - there is something wrong here because i am doing the exact same thing as earlier, but now it is not working. the only difference that instead of making a lambda function in run_0d_solver and then sending it OpenLoopCoronaryWithDistalPressureBlock_v2, I sent a matrix to OpenLoopCoronaryWithDistalPressureBlock_v2 and use that matrix to create the lambda function from within OpenLoopCoronaryWithDistalPressureBlock_v2. These 2 methods should be equivalent, so whats going on here. did i actually do something different earlier??

        # print("Pim_value = ", Pim_value)
        self.cveccoe = [-1.0*self.Cim*Pim_value + self.Cim*Pv_value, -1.0*self.Cim*(self.Rv + self.Ram)*Pim_value + self.Ram*self.Cim*Pv_value]

    def update_constant(self):
        self.emxcoe = [(-1.0*self.Ca*self.Cim*self.Rv, self.Ra*self.Ca*self.Cim*self.Rv, -1.0*self.Cim*self.Rv), (0.0, 0.0, -1.0*self.Cim*self.Rv*self.Ram)]
        self.fmxcoe = [(0.0, self.Cim*self.Rv, -1.0), (self.Cim*self.Rv, -1.0*self.Cim*self.Rv*self.Ra, -1.0*(self.Rv + self.Ram))]

        # self.demxcoe = [(0,)*3]*self.neq
        # self.dfmxcoe = [(0,)*3]*self.neq
        # self.dcmxcoe = [(0,)*3]*self.neq
###############################################################

# -- Time Varying Capacitance
class TimeDependentCapacitance(LPNBlock):
    def __init__(self,Cfunc,connecting_block_list=None,name="NoNameTimeDependentCapacitance",flow_directions=None):

        LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)
        self.type = "TimeDependentCapacitance"
        self.Cfunc = Cfunc

    def check_block_consistency(self):
        if len(connecting_block_list) != 2:
            raise Exception("Capacitance block can be connected only to two elements")


    def update_time(self,args):
        t = args['Time']
        # dt = args['Time step']
        # rho = args['rho']
        # alpha_f = 1/(1.0+rho)
        self.emxcoe = [(1.0*self.Cfunc(t),0,-1.0*self.Cfunc(t),0),(0,0,0,0)]

    def update_constant(self):
        self.fmxcoe = [(0,-1.0,0,0),(0,1.,0,-1.)]
        # self.cveccoe = [0,0]
        #
        # self.demxcoe = [(0,)*4]*2
        # self.dfmxcoe = [(0,)*4]*2
        # self.dcmxcoe = [(0,)*4]*2


# -- Chamber Capacitance -- with direct prescription of pressure
class ChamberModel(LPNBlock):
    def __init__(self,Activefunc,Passivefunc,Pref,connecting_block_list=None,name="NoNameChamber",flow_directions=None):

        LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)

        self.neq = 3
        self.type = "Chamber"
        self.num_block_vars = 1 # Chamber volume
        self.Activefunc = Activefunc
        self.Passivefunc = Passivefunc
        self.Pref = Pref


    def linearize_func(self,t,v,func,eps_dv=1e-4):

        # Input  : func(t,v), typically t = t_n+1 ie t_n + dt
        # Return : d/dv (func) at (t,v) using central differences

        dv = eps_dv*v
        return (func(t,v+dv)-func(t,v-dv))/(2.*dv)

    def linearize_func_t(self,t,v,func,dt=1e-3):

        return (func(t+dt,v)-func(t-dt,v))/(2.*dt)

    def check_block_consistency(self):
        if len(connecting_block_list) != 2:
            raise Exception("Chamber block can be connected only to two elements")
        return

    def initialize_volume(self,sol_vec,Vu=1.0,scale_fact=1.5):

        sol_vec[self.LPN_solution_ids[0]] = scale_fact*Vu

    def update_solution(self,args):
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
        if v<=0:
            print("Time: ",t)
            raise Exception("Zero chamber volume detected in chamber: ",self.name)

        a_lin = self.linearize_func(t,v,self.Activefunc)
        p_lin = self.linearize_func(t,v,self.Passivefunc)

        # c_from_lin = -(self.Activefunc(t,v) + self.Passivefunc(t,v)) + v*(a_lin+p_lin)
        # print 'Volume ',v,'has pressure ',(self.Activefunc(t,v) + self.Passivefunc(t,v))/1333.33

        self.emxcoe = [(0,0,0,0,1.),     (0,0,0,0,0),     (0,0,0,0,0)]
        self.fmxcoe = [(0,-1.,0,1.,0), (-1.,0,1.,0,0),  (1.,0,0.,0.,0.)]
        self.cveccoe = [0,0,-(self.Activefunc(t,v)+self.Passivefunc(t,v))-self.Pref]

        # self.fmxcoe = [(0,-1.,0,1.,0), (-1.,0,1.,0,0),  (1.,0,0.,0,-(a_lin+p_lin))]
        # self.cveccoe = [0,0,c_from_lin+self.Pref]

        # self.demxcoe = [(0,)*5,(0,)*5,(0,)*5]
        # self.dfmxcoe = [(0,)*5,(0,)*5,(0,)*5]
        self.dcmxcoe = [(0,)*5,(0,)*5,(0,0,0,0,-(a_lin+p_lin))]



# -- Inductance
class Inductance(LPNBlock):
    def __init__(self,L,connecting_block_list=None,name="NoNameInductance",flow_directions=None):

        LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)
        self.type = "Inductance"
        self.L = L

    def check_block_consistency(self):
        if len(connecting_block_list) != 2:
            raise Exception("Inductance block can be connected only to two elements")

    def update_constant(self):
        self.emxcoe = [(0.0,-1.,0.0,0),(0,0,0,0)]
        self.fmxcoe = [(1./self.L,0.0,-1./self.L,0),(0,1.,0,-1.)]
        # self.cveccoe = [0,0]
        #
        # self.demxcoe = [(0,)*4]*2
        # self.dfmxcoe = [(0,)*4]*2
        # self.dcmxcoe = [(0,)*4]*2


# -- Ideal diode - state variable
class IdealDiode2(LPNBlock):
    def __init__(self,connecting_block_list=None,eps=1e-17,name="NoNameIdealDiode",flow_directions=None):

        LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)
        self.type = "IdealDiode"
        self.neq=3
        self.num_block_vars = 1 # State : 0- shunt, 1- interrupt
        self.eps =eps
        # self.Peps = Peps

    def check_block_consistency(self):
        if len(connecting_block_list) != 2:
            raise Exception("IdealDiode block can be connected only to two elements")

    def update_solution(self, args):
        # Needs to take current solution vector as argument
        t = args['Time']
        curr_y = args['Solution']
        wire_dict = args['Wire dictionary']
        rho = args['rho']
        dt = args['Time step']
        alpha_f = 1.0/(1.0+rho)
        Qi = curr_y[wire_dict[self.connecting_wires_list[0]].LPN_solution_ids[1]]
        Qo = curr_y[wire_dict[self.connecting_wires_list[1]].LPN_solution_ids[1]]
        Pi = curr_y[wire_dict[self.connecting_wires_list[0]].LPN_solution_ids[0]]
        Po = curr_y[wire_dict[self.connecting_wires_list[1]].LPN_solution_ids[0]]
        state = curr_y[self.LPN_solution_ids[0]]
        # self.emxcoe = [(0,)*5]*3
        # self.demxcoe = [(0,)*5]*3
        # self.dfmxcoe = [(0,)*5]*3
        # self.dcmxcoe = [(0,)*5]*3

        if state > 0 :
            # Zero flow
            self.fmxcoe = [(0,1,0,0,0),(0,0,0,1.,0)]
            self.cveccoe = [0.,0.]

        else :
            # Perfect wire
            self.fmxcoe = [(1.,0,-1.,0,0),(0,1.,0,-1.,0)]
            self.cveccoe = [0,0]

        if Qi < -self.eps or Qo < -self.eps :
            self.fmxcoe.append((0,0,0,0,1))
            self.cveccoe.append(-alpha_f-(1-alpha_f)*state)
            # self.cveccoe.append(-1)
            # print 'Flow toggle'
            # print state

        elif Pi - Po > -self.eps :
            self.fmxcoe.append((0,0,0,0,1))
            self.cveccoe.append(-(1-alpha_f)*state)
            # self.cveccoe.append(0)
            # print 'P toggle'
            # print state
        else :
            self.fmxcoe.append((0,0,0,0,1))
            self.cveccoe.append(-state)
            # print 'State maintain'
            # print state



# -- Pressure source (two terminal : analogous to a battery)
class PressureSource(LPNBlock):
    def __init__(self,Pfunction,Pref,connecting_block_list=None,name="NoNamePressureSource",flow_directions=None):

        LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)
        self.neq = 2
        self.type = "PressureSource"
        self.Pfunction=Pfunction
        self.Pref= Pref

    def check_block_consistency(self):
        if len(connecting_block_list) != 2:
            raise Exception("PressureSource block can be connected only to two elements")

    def update_time(self, args):
        t = args['Time']
        # rho = args['rho']
        # dt = args['Time step']
        # alpha_f = 1.0/(1.0+rho)
        self.cveccoe = [-1.0*self.Pref,-1.0*self.Pfunction(t)]
        # self.cveccoe = [-1.0*self.Pref,-1.0*self.Pfunction(t)]

    def update_constant(self):
        # self.emxcoe = [(0,0)]*2
        self.fmxcoe = [(1.,0,0,0),(0,0,1.,0)]

        # self.demxcoe = [(0,)*4]*2
        # self.dfmxcoe = [(0,)*4]*2
        # self.dcmxcoe = [(0,)*4]*2


# # -- Flow source (single terminal)
# class FlowSource(LPNBlock):
#     def __init__(self,Qfunction,connecting_block_list=None,name="NoNameFlowSource",flow_directions=None):

#         LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)
#         self.neq = 1
#         self.type = "FlowSource"
#         self.Qfunction=Qfunction

#     def check_block_consistency(self):
#         if len(connecting_block_list) != 1:
#             raise Exception("FlowSource block can be connected only to one elements")

#     def update_time(self,args):
#         t = args['Time']

#         # rho = args['rho']
#         # dt = args['Time step']
#         # alpha_f = 1.0/(1.0+rho)

#         self.emxcoe = [(0,0)]
#         self.fmxcoe = [(0,1.)]
#         # self.cveccoe = [-1.0*self.Qfunction(t+alpha_f*dt)]
#         self.cveccoe = [-1.0*self.Qfunction(t)]

#         self.demxcoe = [(0,)*4]*2
#         self.dfmxcoe = [(0,)*4]*2
#         self.dcmxcoe = [(0,)*4]*2


class UnsteadyResistanceWithDistalPressure_special(LPNBlock):
    def __init__(self,Rfunc,Pref_func,connecting_block_list=None,name="NoNameUnsteadyResistanceWithDistalPressure_special",flow_directions=None):

        LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)
        self.type = "UnsteadyResistanceWithDistalPressure_special"
        self.neq = 1
        self.Rfunc = Rfunc
        self.Pref_func = Pref_func


        # input("UnsteadyResistanceWithDistalPressure_special - click to continue")


    def check_block_consistency(self):
        if len(connecting_block_list) != 1:
            raise Exception("UnsteadyResistanceWithDistalPressure_special block can be connected only to two elements")

    def update_time(self, args):
        # For resistors, the ordering is : (P_in,Q_in)
        t = args['Time']

        # currently here 7/19/20: since this resistor element with distal pressure is supposed to be a terminal element, so that i dont have to connect a pressureRef block downstream, how many unknowns and equations should this block have? I think its 2 unknowns (P_in,Q_in), but how many equations do i need?
        self.fmxcoe = [(1.,-1.0*self.Rfunc(t))]
        self.cveccoe = [-1.0*self.Pref_func(t)]

    def update_constant(self):
        self.emxcoe = [(0,)*2]

        self.demxcoe = [(0,)*2]
        self.dfmxcoe = [(0,)*2]
        self.dcmxcoe = [(0,)*2]

# -- Resistance
class UnsteadyResistance_special(LPNBlock):
    def __init__(self,Rfunc,connecting_block_list=None,name="NoNameUnsteadyResistance_special",flow_directions=None):

        LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)
        self.type = "UnsteadyResistance_special"
        self.Rfunc = Rfunc

        # input("UnsteadyResistance_special - click to continue")


    def check_block_consistency(self):
        if len(connecting_block_list) != 2:
            raise Exception("UnsteadyResistance_special block can be connected only to two elements")

    def update_time(self, args):
        # For resistors, the ordering is : (P_in,Q_in,P_out,Q_out)
        t = args['Time']
        self.fmxcoe = [(1.,-1.0*self.Rfunc(t),-1.,0),(0,1.,0,-1.)]

    def update_constant(self):
        self.emxcoe = [(0,)*4]*2
        self.cveccoe = [0]*2

        self.demxcoe = [(0,)*4]*2
        self.dfmxcoe = [(0,)*4]*2
        self.dcmxcoe = [(0,)*4]*2

class UnsteadyFlowRef_special(LPNBlock):
    def __init__(self,Qfunc,connecting_block_list=None,name="NoNameUnsteadyFlowRef_special",flow_directions=None):

        LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)
        self.type = "UnsteadyFlowRef_special"
        self.neq = 1
        self.Qfunc=Qfunc

        # input("UnsteadyFlowRef_special - click to continue")


    def check_block_consistency(self):
        if len(connecting_block_list) != 1:
            raise Exception("UnsteadyFlowRef_special block can be connected only to one element")

    def update_time(self, args):
        t = args['Time']
        self.cveccoe = [-1.0*self.Qfunc(t)]

    def update_constant(self):
        self.emxcoe = [(0,0)]
        self.fmxcoe = [(0,1.)]

        self.demxcoe = [(0,)]
        self.dfmxcoe = [(0,)]
        self.dcmxcoe = [(0,)]

class Junction_special(LPNBlock):
    def __init__(self,temp_parameter,connecting_block_list=None,name="NoNameJunction_special",flow_directions=None):
        LPNBlock.__init__(self,connecting_block_list,name=name,flow_directions=flow_directions)
        self.type = "Junction_special"
        self.neq = self.num_connections
        self.temp_parameter = temp_parameter

        print("------------------------------------------")
        print("Junction_special - click to continue")
        print("temp_parameter = ", self.temp_parameter)
        print("------------------------------------------")


    def add_connecting_block(self,block,direction):
        self.connecting_block_list.append(block)
        self.num_connections = len(self.connecting_block_list)
        self.neq = self.num_connections
        self.flow_directions.append(direction)
        # print self.name

    def update_constant(self):

        # Number of variables per tuple = 2*num_connections
        # Number of equations = num_connections-1 Pressure equations, 1 flow equation
        # Format : P1,Q1,P2,Q2,P3,Q3, .., Pn,Qm

        self.emxcoe = [(0,)*(2*self.num_connections)]*(self.num_connections)

        self.fmxcoe = [ (1.,)+(0,)*(2*i+1) + (-1,) + (0,)*(2*self.num_connections-2*i-3) for i in range(self.num_connections-1) ]

        tmp = (0,)
        for d in self.flow_directions[:-1]:
            tmp+=(d,)
            tmp+=(0,)

        tmp += (self.flow_directions[-1],)
        self.fmxcoe.append(tmp)
        self.cveccoe = [0]*self.num_connections

        self.demxcoe = [(0,)*(2*self.num_connections)]*(self.num_connections)
        self.dfmxcoe = [(0,)*(2*self.num_connections)]*(self.num_connections)
        self.dcmxcoe = [(0,)*(2*self.num_connections)]*(self.num_connections)
