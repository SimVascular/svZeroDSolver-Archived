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

LPNVariable::LPNVariable() {}

LPNVariable::LPNVariable(std::string name, std::string type, double value) : name(name), type(type), value(value) {} // this is a member initialization/initializer list (between the colon and the squigly bracket)

std::string LPNVariable::GetName() const {
    return name;
}

std::string LPNVariable::GetType() const {
    return type;
}

double LPNVariable::GetValue() const {
    return value;
}

PressureVariable::PressureVariable() {}

FlowVariable::FlowVariable() {}

PressureVariable::PressureVariable(std::string name, std::string type, double value) : LPNVariable(name, type, value) {}

FlowVariable::FlowVariable(std::string name, std::string type, double value) : LPNVariable(name, type, value) {}

LPNBlock::LPNBlock() {}

LPNBlock::LPNBlock(std::string name, std::string type, std::vector<std::string> connecting_block_list, std::vector<int> flow_directions) : name(name), type(type), connecting_block_list(connecting_block_list), flow_directions(flow_directions) {}

std::string LPNBlock::GetName() const {
    return name;
}

std::string LPNBlock::GetType() const {
    return type;
}

std::vector<std::string> LPNBlock::GetConnectingBlockList() const {
    return connecting_block_list;
}

std::vector<int> LPNBlock::GetFlowDirections() const {
    return flow_directions;
}

int LPNBlock::GetNumConnections() const {
    return num_connections;
}

int LPNBlock::GetNeq() const {
    return neq;
}

void LPNBlock::AddConnectingBlock(std::string block_name, int direction) {
    // direction = +1 if flow sent to new block
    //           = -1 if flow received from new block
    connecting_block_list.push_back(block_name);
    num_connections = connecting_block_list.size();
    flow_directions.push_back(direction);
}

void LPNBlock::AddConnectingWire(std::string wire_name) {
    connecting_wires_list.push_back(wire_name);
}
// 
// int LPNBlock::GetEquationId(std::unordered_map<std::string, Wire *>, int) {
// 
// }

void LPNBlock::UpdateConstant(Args * args) {
    
}

void LPNBlock::UpdateTime(Args * args) {
    
}

void LPNBlock::UpdateSolution(Args * args) {
    
}

Wire::Wire(std::string name, std::array<LPNBlock *, 2> connecting_block_list) : name(name), connecting_block_list(connecting_block_list) {}

void Wire::SetLPNSolutionIds(std::array<int, 2> solution_ids) {
    lpn_solution_ids = solution_ids;
}

std::string Wire::GetName() const {
    return name;
}

std::array<int, 2> Wire::GetLPNSolutionIds() const {
    return lpn_solution_ids;
}

std::array<LPNBlock *, 2> Wire::GetConnectingBlockList() const {
    return connecting_block_list;
}