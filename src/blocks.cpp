// Copyright (c) Stanford University, The Regents of the University of
//               California, and others.
//
// All Rights Reserved.
//
// See Copyright-SimVascular.txt for additional details.
//
// Permission is hereby granted, free of charge, to any person obtaining
// a copy of this software and associated documentation files (the
// "Software"), to deal in the Software without restriction, including
// without limitation the rights to use, copy, modify, merge, publish,
// distribute, sublicense, and/or sell copies of the Software, and to
// permit persons to whom the Software is furnished to do so, subject
// to the following conditions:
//
// The above copyright notice and this permission notice shall be included
// in all copies or substantial portions of the Software.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
// IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
// PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER
// OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

//#include <string>
#include "blocks.hpp"

Wire::Wire(const std::string& name) : name(name) {} // this is a member initialization/initializer list (between the colon and the squigly bracket)

Wire::~Wire() {} // https://www.mycplus.com/tutorials/cplusplus-programming-tutorials/destructors/

void Wire::set_lpn_solution_ids(const std::array<int, 2>& solution_ids) 
{
    lpn_solution_ids = solution_ids;
}

std::string Wire::get_name() const 
{
    return name;
}

std::array<int, 2> Wire::get_lpn_solution_ids() const 
{
    return lpn_solution_ids;
}

Args::Args(double dt, double rho, bool check_jacobian, const std::unordered_map<std::string, Wire *>& wire_dict) : dt(dt), rho(rho), check_jacobian(check_jacobian), wire_dict(wire_dict) {}

Args::~Args() {}

double Args::get_dt() const
{
    return dt;
}

double Args::get_rho() const
{
    return rho;
}

bool Args::get_check_jacobian() const
{
    return check_jacobian;
}

 std::unordered_map<std::string, Wire *> Args::get_wire_dict() const
{
    return wire_dict;
}

double Args::get_time() const
{
    return time;
}

LPNBlock::LPNBlock() {}

LPNBlock::LPNBlock(const std::string& name, const std::vector<std::string>& connecting_block_list, const std::vector<int>& flow_directions) : name(name), connecting_block_list(connecting_block_list), flow_directions(flow_directions) {}

LPNBlock::~LPNBlock() {} // https://stackoverflow.com/questions/3065154/undefined-reference-to-vtable

void LPNBlock::set_n_eqn(int neqn) 
{
    n_eqn = neqn;
}

std::string LPNBlock::get_name() const 
{
    return name;
}

std::vector<std::string> LPNBlock::get_connecting_block_list() const 
{
    return connecting_block_list;
}

int LPNBlock::get_n_eqn() const 
{
    return n_eqn;
}

int LPNBlock::get_n_block_var() const 
{
    return n_block_var;
}

std::vector<std::string> LPNBlock::get_connecting_wires_list() const
{
    return connecting_wires_list;
}

std::vector<int> LPNBlock::get_flow_directions() const 
{
    return flow_directions;
}

std::vector<int> LPNBlock::get_lpn_solution_ids() const
{
    return lpn_solution_ids;
}

std::unordered_map<char, std::vector<std::vector<double>>> LPNBlock::get_mat() const
{
    return mat;
}

std::vector<int> LPNBlock::get_global_col_id() const
{
    return global_col_id;
}
std::vector<int> LPNBlock::get_global_row_id() const
{
    return global_row_id;
}

void LPNBlock::add_connecting_block(const std::string& block_name, int direction) 
{
    // direction = +1 if flow sent from this to block
    //           = -1 if flow sent from block to this
    connecting_block_list.push_back(block_name);
    flow_directions.push_back(direction);
}

void LPNBlock::add_connecting_wire(const std::string& wire_name) 
{
    connecting_wires_list.push_back(wire_name);
}

void LPNBlock::update_constant() {}

void LPNBlock::update_time(Args * args) {}

void LPNBlock::update_solution(Args * args) {}

int LPNBlock::get_global_equation_ids(std::unordered_map<std::string, Wire *> wire_dict, int local_solution_id)
{
    int n_wire_var = connecting_block_list.size() * 2;  // there are 2 soltns (P and Q) per wire
    if (local_solution_id < n_wire_var) { // wire solution variable
        int var_type = local_solution_id % 2; // 0 --> P, 1 --> Q
        int local_wire_id = local_solution_id / 2;

        // Example: Assume connecting_block_list.size() is 2. This can be a normal resistor block, which has 2 connections. Then this R block has 2 connecting wires. thus, This R block has 4 related solution variables (P_in, Q_in, P_out, Q_out).
        //     then for these are the var_types we get for each local_solution_id:
        //  local_solution_id  : var_type :  local_wire_id
        //         0           :     0    :       0        <---  var_type = pressure, local_wire_id = inlet wire
        //         1           :     1    :       0        <---  var_type = flow,     local_wire_id = inlet wire
        //         2           :     0    :       1        <---  var_type = pressure, local_wire_id = outlet wire
        //         3           :     1    :       1        <---  var_type = flow,     local_wire_id = outlet wire
        
        return wire_dict[connecting_wires_list[local_wire_id]]->get_lpn_solution_ids()[var_type];
    } else { // block internal solution variable
        int internal_var_local_id = local_solution_id - n_wire_var;
        return get_lpn_solution_ids()[internal_var_local_id];
    }
}

Junction::Junction(const std::string& name, const std::vector<std::string>& connecting_block_list, const std::vector<int>& flow_directions) : LPNBlock(name, connecting_block_list, flow_directions) 
{
    set_n_eqn(connecting_block_list.size()); // the equations are 1) mass conservation 2) inlet pressures = outlet pressures
}

Junction::~Junction() {}

void Junction::update_constant()
{
    // Number of variables per sub-vector = 2 * connecting_block_list.size()
    // Number of equations: connecting_block_list.size() - 1 pressure equations, 1 flow equation
    // Format : P1, Q1, P2, Q2, P3, Q3, .., Pn, Qm
    
    // std::unordered_map<char, std::vector<std::vector<double>>> mat
    
    for (int i = 0; i < connecting_block_list.size() - 1; i++) {
        std::vector<double> tmp {1};
        for (int j = 0; j < 2 * i + 1; j++) {
            tmp.push_back(0); 
        }
        tmp.push_back(-1); 
        for (int j = 0; j < 2 * connecting_block_list.size() - 2 * i - 3; j++) {
            tmp.push_back(0); 
        }
        mat['F'].push_back(tmp);
    }
    
    std::vector<double> tmp {0};
    for (int d : flow_directions) {
        tmp.push_back(d);
        tmp.push_back(0);
    }
    mat['F'].push_back(tmp); // std::vector<T>::push_back() creates a copy of the argument and stores it in the vector; https://stackoverflow.com/questions/2275076/is-stdvector-copying-the-objects-with-a-push-back
}