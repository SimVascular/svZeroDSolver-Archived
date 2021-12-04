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

// todo: delete unused getters and setters

class LPNBlock; // forward declaration of LPNBlock; source: https://stackoverflow.com/questions/396084/headers-including-each-other-in-c

class Wire {
  // private
  // these fields are inherently private, since they are declared above the "public" space
  std::string name;
  std::array<int, 2> lpn_solution_ids;
  std::array<LPNBlock *, 2> connecting_block_list;

public:
  //constuctor
  Wire(const std::string& name, const std::array<LPNBlock *, 2>& connecting_block_list); // todo: LPNBlock also has a "connecting_block_list" field, but that one is a vector of strings instead. Therefore, consider changing Wire's connecting_block)list from an array of LPNBlock pointers to a vector of strings as well. This will create consistency between both classes' similarly-named fields. Do this after the entire implementation of svZeroDSolver is done
  
  // destructor
  ~Wire(); // https://www.learncpp.com/cpp-tutorial/destructors/
    
  // setters
  void set_lpn_solution_ids(const std::array<int, 2>& lpn_solution_ids);

  // getters
  std::string get_name() const;
  std::array<int, 2> get_lpn_solution_ids() const;
  std::array<LPNBlock *, 2> get_connecting_block_list() const;
};

class Args {
  // private
  // these fields are inherently private, since they are declared above the "public" space
  double dt; // time step size
  double rho; // generalized-alpha rho parameter
  bool check_jacobian;
  std::unordered_map<std::string, Wire *> wire_dict;
  double time; // current time?? or time_af??
  
  // some_kind_of_eigen_vector yaf; // (the "Solution" key)
  
public:
  // constructors
  Args(double dt, double rho, bool check_jacobian, const std::unordered_map<std::string, Wire *>& wire_dict);
  
  // destructor
  ~Args();
  
  // getters
  double get_dt() const;
  double get_rho() const;
  bool get_check_jacobian() const;
  std::unordered_map<std::string, Wire *> get_wire_dict();
  double get_time() const;
  
  // misc
  void update_time();
  void update_yaf();
  
  // last here - need to make unit tests and test cases for EVERYTHING INSIDE ARGS
};

class LPNBlock {
  // private
  // these fields are inherently private, since they are declared above the "public" space
  std::string name;
  std::string type;
  std::vector<std::string> connecting_block_list;
  std::vector<int> flow_directions;
  int num_connections;
  int neq;
  int num_block_vars;
  std::vector<std::string> connecting_wires_list;
  std::vector<int> lpn_solution_ids; // solution IDs for the LPN block's internal solution variables

  std::unordered_map<char, std::vector<std::vector<double>>> mat; // todo: need to create a function to initialize map to have 'E', 'F', 'C', 'dE', 'dF', 'dC' keys with default values of empty vectors for the std::vector<std::vector<double>> values

  // row and column indices of block in global matrices
  std::vector<int> global_col_id;
  std::vector<int> global_row_id;

public:
  // constructors
  LPNBlock(); // default constructor; https://stackoverflow.com/questions/31211319/no-matching-function-for-call-to-class-constructor
  LPNBlock(const std::string& name, const std::string& type, const std::vector<std::string>& connecting_block_list, const std::vector<int>& flow_directions);
  
  // destructor
  virtual ~LPNBlock(); // https://stackoverflow.com/questions/3065154/undefined-reference-to-vtable ; https://www.geeksforgeeks.org/virtual-destructor/

  // setters
  

  // getters
  std::string get_name() const;
  std::string get_type() const;
  std::vector<std::string> get_connecting_block_list() const;
  std::vector<int> get_flow_directions() const;
  int get_num_connections() const;
  int get_neq() const; 

  // misc
  void add_connecting_block(const std::string& block_name, int direction);
  void add_connecting_wire(const std::string& wire_name);
  // int get_equation_id(std::unordered_map<std::string, Wire *>, int); // todo: write function for get_equation_id
  virtual void update_constant(Args * args);
  virtual void update_time(Args * args);
  virtual void update_solution(Args * args); // todo: should this be a "const Args * args" (a pointer to a const Args)? (see "Pointers and const" in https://www.cplusplus.com/doc/tutorial/pointers/) I dont think so because it is possible that update_solution will change one of the fields in args
};

class Junction : public LPNBlock {
  // private
  // these fields are inherently private, since they are declared above the "public" space

public:
  // constructors

  // destructor

  // setters

  // getters
};

class BloodVessel : public LPNBlock {
  // private
  // these fields are inherently private, since they are declared above the "public" space

public:
  // constructors

  // destructor

  // setters

  // getters
};

class UnsteadyResistanceWithDistalPressure : public LPNBlock {
  // private
  // these fields are inherently private, since they are declared above the "public" space

public:
  // constructors
  UnsteadyResistanceWithDistalPressure(std::string name, std::vector<std::string> connecting_block_list, std::vector<int> flow_directions, Rfunc, Pref_func); 
    11/10/21: last here - how to pass a function as an argument to another function, e.g., how to pass Rfunc (a function) as an argument to this UnsteadyResistanceWithDistalPressure constructor?

  // destructor

  // misc
  void update_time(Args * args);
};

class UnsteadyPressureRef : public LPNBlock {
  // private
  // these fields are inherently private, since they are declared above the "public" space

public:
  // constructors
  UnsteadyPressureRef(std::string name, std::vector<std::string> connecting_block_list, std::vector<int> flow_directions, Pfunc);

  // destructor

  // misc
  void update_constant(Args * args);
  void update_time(Args * args);
};

class UnsteadyFlowRef : public LPNBlock {
  // private
  // these fields are inherently private, since they are declared above the "public" space

public:
  // constructors
  UnsteadyFlowRef(std::string name, std::vector<std::string> connecting_block_list, std::vector<int> flow_directions, Qfunc);

  // destructor

  // misc
  void update_constant(Args * args);
  void update_time(Args * args);
};

class UnsteadyRCRBlockWithDistalPressure : public LPNBlock {
  // private
  // these fields are inherently private, since they are declared above the "public" space

public:
  // constructors
  UnsteadyRCRBlockWithDistalPressure(std::string name, std::vector<std::string> connecting_block_list, std::vector<int> flow_directions, Rp_func, C_func, Rd_func, Pref_func);

  // destructor

  // misc
  void update_time(Args * args);
};

class OpenLoopCoronaryWithDistalPressureBlock : public LPNBlock {
  // private
  // these fields are inherently private, since they are declared above the "public" space

public:
  // constructors
  OpenLoopCoronaryWithDistalPressureBlock(std::string name, std::vector<std::string> connecting_block_list, std::vector<int> flow_directions, double Ra, double Ca, double Ram, double Cim, double Rv, Pim, Pv, double cardiac_cycle_period);

  // destructor

  // misc
  double get_P_at_t(P, t);
  void update_constant(Args * args);
  void update_time(Args * args);
};

todo: for unit testing in c++: how to do it best: throw exception or use assert or what?? -- last here - 12/3/21 -- do this first when I this code again

todo: review this code in its entirely to make sure that I recall and understand everything again, before continuing the below todo items

todo: do all todos first, so that I dont forget to do something

todo: do all last here's

#endif