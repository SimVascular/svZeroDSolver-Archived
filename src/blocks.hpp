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

class LPNVariable {
    // private
    // these fields are inherently private, since they are declared above the "public" space
    std::string name;
    std::string type;
    double value;

public:
        // constructor
    LPNVariable();
    LPNVariable(std::string name, std::string type, double value);

    // todo: add destructor: https://www.learncpp.com/cpp-tutorial/destructors/
    // todo: use static and const variables: https://www.learncpp.com/cpp-tutorial/const-class-objects-and-member-functions/ ; https://www.learncpp.com/cpp-tutorial/const-constexpr-and-symbolic-constants/

    // getter methods
    std::string GetName() const; // todo: is this function needed
    std::string GetType() const; // todo: is this function needed
    double GetValue() const; // todo: is this function needed
};

class PressureVariable : public LPNVariable {
public:
    PressureVariable(); // default constructor; https://stackoverflow.com/questions/31211319/no-matching-function-for-call-to-class-constructor
    PressureVariable(std::string name, std::string type, double value);
    // todo: need to add destructor
};

class FlowVariable : public LPNVariable {
public:
    FlowVariable();
    FlowVariable(std::string name, std::string type, double value);
    // todo: need to add destructor
};

class LPNBlock; // forward declaration of LPNBlock; source: https://stackoverflow.com/questions/396084/headers-including-each-other-in-c

class Wire {
    // private
    // these fields are inherently private, since they are declared above the "public" space
    std::string name;
    //PressureVariable P; // todo: do I need to set P via pointer/reference?
    //FlowVariable Q; // todo: do I need to set Q via pointer/reference?
    std::array<int, 2> lpn_solution_ids;
    std::array<LPNBlock *, 2> connecting_block_list;

public:
    //constuctor
    
    //todo: in the constructor of Wire, I think I should make sure to pass in the LPNBlocks by reference or by pointer or whatever. But do I pass in connecting_block_list as a pointer or do i pass in the LPNBlocks stored in that list by reference/pointer?
    Wire(std::string name, std::array<LPNBlock *, 2> connecting_block_list); // todo: LPNBlock also has a "connecting_block_list" field, but that one is a vector of strings instead. Therefore, consider changing Wire's connecting_block)list from an array of LPNBlock pointers to a vector of strings as well. This will create consistency between both classes' similarly-named fields.
    
    
    // todo: need to add destructor

    // setters
    void SetLPNSolutionIds(std::array<int, 2> lpn_solution_ids);
    //void SetP(PressureVariable P);
    //void SetQ(FlowVariable Q);

    // getters
    std::string GetName() const;
    //PressureVariable GetP() const;
    //FlowVariable GetQ() const;
    std::array<int, 2> GetLPNSolutionIds() const;
    std::array<LPNBlock *, 2> GetConnectingBlockList() const;
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
    Args(double dt, double rho, bool check_jacobian, std::unordered_map<std::string, Wire *> wire_dict);
    
    // todo: need to add destructor
    
    // getters
    double GetDt() const;
    double GetRho() const;
    bool GetCheckJacobian() const;
    std::unordered_map<std::string, Wire *> GetWireDict();
    double GetTime() const;
    
    // misc
    void UpdateTime();
    void UpdateYaf();
    
    // last here - need to make unit tests and test cases for EVERYTHING INSIDE ARGS
};

class LPNBlock {
    // private
    // these fields are inherently private, since they are declared above the "public" space
    std::string name;
    std::string type;
    std::vector<std::string> connecting_block_list;
    std::vector<int> flow_directions;
    int num_connections; // todo: need to create a function to set num_connections
    int neq; // todo: need to create a function to set neq
    int num_block_vars; // todo: need to create a function to set num_block_vars
    std::vector<std::string> connecting_wires_list; // todo: need to create a function to set connecting_wires_list
    std::vector<int> lpn_solution_ids; // solution IDs for the LPN block's internal solution variables // todo: need to create a function to set lpn_solution_ids

    std::unordered_map<char, std::vector<std::vector<double>>> mat; // todo: need to create a function to initialize map to have 'E', 'F', 'C', 'dE', 'dF', 'dC' keys with default values of empty vectors for the std::vector<std::vector<double>> values

    // row and column indices of block in global matrix
    std::vector<int> global_col_id; // todo: need to create a function to set global_col_id
    std::vector<int> global_row_id; // todo: need to create a function to set global_row_id

public:
    // constructors
    LPNBlock();
    LPNBlock(std::string name, std::string type, std::vector<std::string> connecting_block_list, std::vector<int> flow_directions);
    
    // destructor
    // virtual ~LPNBlock(); // https://stackoverflow.com/questions/3065154/undefined-reference-to-vtable ; https://www.geeksforgeeks.org/virtual-destructor/ // do I need to make this destructor virtual?

    // setters
    

    // getters
    std::string GetName() const;
    std::string GetType() const;
    std::vector<std::string> GetConnectingBlockList() const;
    std::vector<int> GetFlowDirections() const;
    int GetNumConnections() const;
    int GetNeq() const; 

    // misc
    void AddConnectingBlock(std::string block_name, int direction);
    void AddConnectingWire(std::string wire_name);
    int GetEquationId(std::unordered_map<std::string, Wire *>, int); // todo: write function for GetEquationId
    virtual void UpdateConstant(Args * args);
    virtual void UpdateTime(Args * args);
    virtual void UpdateSolution(Args * args);
    
    // todos: add unit tests for AddConnectingWire, AddConnectingBlock
};


//last here - continue adding the rest of blocks.py here
#endif