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

#ifndef BLOCKS_H
#define BLOCKS_H

#include <string>
#include <array>
#include <vector>
#include <unordered_map>

// todo: use static and const variables where necessary: https://www.learncpp.com/cpp-tutorial/const-class-objects-and-member-functions/ ; https://www.learncpp.com/cpp-tutorial/const-constexpr-and-symbolic-constants/

// todo: delete unused getters and setters -- maybe to make my code simpler and cleaner, I should delete them now and then add them back later when I actually them? according to https://stackoverflow.com/questions/1568091/why-use-getters-and-setters-accessors , maybe I should just use getters and setters for every variable

class LPNBlock; // forward declaration of LPNBlock; source: https://stackoverflow.com/questions/396084/headers-including-each-other-in-c

class Wire 
{
    std::string name;
    std::array<int, 2> lpn_solution_ids;

public:
    Wire(const std::string& name);
    ~Wire(); // https://www.learncpp.com/cpp-tutorial/destructors/
    
    // setters
    void set_lpn_solution_ids(const std::array<int, 2>& lpn_solution_ids);

    // getters
    std::string get_name() const;
    std::array<int, 2> get_lpn_solution_ids() const; // https://www.learncpp.com/cpp-tutorial/const-class-objects-and-member-functions/
};

class Args 
{
    double dt; // time step size
    double rho; // generalized-alpha rho parameter
    bool check_jacobian;
    std::unordered_map<std::string, Wire *> wire_dict;
    double time; // current time?? or time_af??
    
    // some_kind_of_eigen_vector yaf; // (the "Solution" key)
    
public:
    Args(double dt, double rho, bool check_jacobian, const std::unordered_map<std::string, Wire *>& wire_dict);
    ~Args();
    
    // getters
    double get_dt() const;
    double get_rho() const;
    bool get_check_jacobian() const;
    std::unordered_map<std::string, Wire *> get_wire_dict() const;
    double get_time() const;
};

class LPNBlock 
{
    int n_eqn;
    int n_block_var;
    std::vector<std::string> connecting_wires_list;
    
    // solution IDs for the LPN block's internal solution variables
    std::vector<int> lpn_solution_ids;

    // row and column indices of block in global matrices
    std::vector<int> global_col_id;
    std::vector<int> global_row_id;

public:
    LPNBlock(); // default constructor; https://stackoverflow.com/questions/31211319/no-matching-function-for-call-to-class-constructor
    LPNBlock(const std::string& name, const std::vector<std::string>& connecting_block_list, const std::vector<int>& flow_directions);
    virtual ~LPNBlock(); // https://stackoverflow.com/questions/3065154/undefined-reference-to-vtable ; https://www.geeksforgeeks.org/virtual-destructor/

    // getters
    std::string get_name() const;
    std::vector<std::string> get_connecting_block_list() const;
    int get_n_eqn() const; 
    int get_n_block_var() const;
    std::vector<std::string> get_connecting_wires_list() const;
    std::vector<int> get_flow_directions() const;
    std::vector<int> get_lpn_solution_ids() const; 
    std::unordered_map<char, std::vector<std::vector<double>>> get_mat() const;
    std::vector<int> get_global_col_id() const;
    std::vector<int> get_global_row_id() const;

    // misc
    
    void add_connecting_block(const std::string& block_name, int direction);
    void add_connecting_wire(const std::string& wire_name);
    
    /* Update solution- and time-independent blocks */
    virtual void update_constant();
    
    /* Update time-dependent blocks */
    virtual void update_time(Args * args);
    
    /* Update solution-dependent blocks */
    virtual void update_solution(Args * args);
    
    /* Return the index at which the local solution variable resides in the global vector of solution variables. */
    int get_global_equation_ids(std::unordered_map<std::string, Wire *> wire_dict, int local_solution_id);
    
protected:
    std::string name;
    std::vector<std::string> connecting_block_list;
    
    // flow direction is -1 for flow into this block and +1 for flow out of this block
    std::vector<int> flow_directions;
    
    // block matrices
    std::unordered_map<char, std::vector<std::vector<double>>> mat; // todo: need to create a function to initialize map to have 'E', 'F', 'C', 'dE', 'dF', 'dC' keys with default values of empty vectors for the std::vector<std::vector<double>> values
    
    // setters
    void set_n_eqn(int n_eqn);
};

class Junction : public LPNBlock 
{
public:
    Junction(const std::string& name, const std::vector<std::string>& connecting_block_list, const std::vector<int>& flow_directions);
    ~Junction();
    
    // misc
    void update_constant() override; // override keyword: https://www.programiz.com/cpp-programming/virtual-functions
};

// class BloodVessel : public LPNBlock 
// {
// 
// public:
//     // constructors
// 
//     // destructor
// 
//     // setters
// 
//     // getters
// };
// 
// class UnsteadyResistanceWithDistalPressure : public LPNBlock 
// {
// 
// public:
//     // constructors
//     UnsteadyResistanceWithDistalPressure(std::string name, std::vector<std::string> connecting_block_list, std::vector<int> flow_directions, Rfunc, Pref_func); 
//     11/10/21: last here - how to pass a function as an argument to another function, e.g., how to pass Rfunc (a function) as an argument to this UnsteadyResistanceWithDistalPressure constructor?
// 
//     // destructor
// 
//     // misc
//     void update_time(Args * args);
// };
// 
// class UnsteadyPressureRef : public LPNBlock 
// {
// 
// public:
//     // constructors
//     UnsteadyPressureRef(std::string name, std::vector<std::string> connecting_block_list, std::vector<int> flow_directions, Pfunc);
// 
//     // destructor
// 
//     // misc
//     void update_constant(Args * args);
//     void update_time(Args * args);
// };

class UnsteadyFlowRef : public LPNBlock 
{
    // last here - create this UnsteadyFlowRef function first, so that I can connect an instance of it to the junction instance and then test Junction and UnsteadyFlowRef blocks
    Qfunc - so I can make Qfunc a function pointer that points to an actual function somewhere else, but if I make to pass a spline in as the function for Qfunc, then how would I do that? I cant just write a brand new function to create a spline, can I? Maybe rather than using Eigen to create the spline, I could write a general function that solves for a general cubic spline and returns a function pointer to that spline function (which would take in as inputs, the coefficients defining the spline and a variable x, which dictates where to evaluate that spline).
public:
    // constructors
    UnsteadyFlowRef(std::string name, std::vector<std::string> connecting_block_list, std::vector<int> flow_directions, Qfunc);

    // destructor

    // misc
    void update_constant(Args * args);
    void update_time(Args * args);
};

// class UnsteadyRCRBlockWithDistalPressure : public LPNBlock 
// {
// 
// public:
//     // constructors
//     UnsteadyRCRBlockWithDistalPressure(std::string name, std::vector<std::string> connecting_block_list, std::vector<int> flow_directions, Rp_func, C_func, Rd_func, Pref_func);
// 
//     // destructor
// 
//     // misc
//     void update_time(Args * args);
// };
// 
// class OpenLoopCoronaryWithDistalPressureBlock : public LPNBlock 
// {
// 
// public:
//     // constructors
//     OpenLoopCoronaryWithDistalPressureBlock(std::string name, std::vector<std::string> connecting_block_list, std::vector<int> flow_directions, double Ra, double Ca, double Ram, double Cim, double Rv, Pim, Pv, double cardiac_cycle_period);
// 
//     // destructor
// 
//     // misc
//     double get_P_at_t(P, t);
//     void update_constant(Args * args);
//     void update_time(Args * args);
// };

#endif