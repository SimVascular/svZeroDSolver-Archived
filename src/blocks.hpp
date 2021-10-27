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

class LPNVariable {
	// private
	// these fields are inherently private, since they are declared above the "public" space
	std::string name;
	std::string type;
	double value;

public:
    	// constructor
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
	PressureVariable(std::string name, std::string type, double value);
};

class FlowVariable : public LPNVariable {
public:
	FlowVariable(std::string name, std::string type, double value);
};

class Wire {
	// private
	// these fields are inherently private, since they are declared above the "public" space
	std::string name;
	PressureVariable P;
	FlowVariable Q;
	std::array<int, 2> lpn_solution_ids;
	std::array<LPNBlock, 2> connecting_block_list;

public:
	Wire(std::string name, std::array<LPNBlock, 2> connecting_block_list);
	void AddLPNSolutionIds(std::array<int, 2> lpn_solution_ids);
	void AddP(PressureVariable P);
	void AddQ(FlowVariable Q);
	std::string GetName() const;
	PressureVariable GetP() const;
	FlowVariable GetQ() const;
	std::array<int, 2> GetLPNSolutionIds() const;
};

#endif
