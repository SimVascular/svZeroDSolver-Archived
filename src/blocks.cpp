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

Wire::Wire(const std::string& name, const std::array<LPNBlock *, 2>& connecting_block_list) : name(name), connecting_block_list(connecting_block_list) {}

Wire::~Wire() {} // https://www.mycplus.com/tutorials/cplusplus-programming-tutorials/destructors/

void Wire::set_lpn_solution_ids(const std::array<int, 2>& solution_ids) {
  lpn_solution_ids = solution_ids;
}

std::string Wire::get_name() const {
  return name;
}

std::array<int, 2> Wire::get_lpn_solution_ids() const {
  return lpn_solution_ids;
}

std::array<LPNBlock *, 2> Wire::get_connecting_block_list() const {
  return connecting_block_list;
}

LPNBlock::LPNBlock() {}

LPNBlock::LPNBlock(const std::string& name, const std::string& type, const std::vector<std::string>& connecting_block_list, const std::vector<int>& flow_directions) : name(name), type(type), connecting_block_list(connecting_block_list), flow_directions(flow_directions) {}  // this is a member initialization/initializer list (between the colon and the squigly bracket)

LPNBlock::~LPNBlock() {} // https://stackoverflow.com/questions/3065154/undefined-reference-to-vtable

std::string LPNBlock::get_name() const {
  return name;
}

std::string LPNBlock::get_type() const {
  return type;
}

std::vector<std::string> LPNBlock::get_connecting_block_list() const {
  return connecting_block_list;
}

std::vector<int> LPNBlock::get_flow_directions() const {
  return flow_directions;
}

int LPNBlock::get_num_connections() const {
  return num_connections;
}

int LPNBlock::get_neq() const {
  return neq;
}

void LPNBlock::add_connecting_block(const std::string& block_name, int direction) {
  // direction = +1 if flow sent to new block
  //           = -1 if flow received from new block
  connecting_block_list.push_back(block_name);
  num_connections = connecting_block_list.size();
  flow_directions.push_back(direction);
}

void LPNBlock::add_connecting_wire(const std::string& wire_name) {
  connecting_wires_list.push_back(wire_name);
}
// 
// int LPNBlock::get_equation_id(std::unordered_map<std::string, Wire *>, int) {
// 
// }

void LPNBlock::update_constant(Args * args) {
  
}

void LPNBlock::update_time(Args * args) {
  
}

void LPNBlock::update_solution(Args * args) {
  
}