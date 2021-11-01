#include <string>
#include <iostream>
#include "blocks.hpp"

int main() {
	std::string name01 = "a";
	std::string type01 = "type_a";
	double value01 = 0.12;
	LPNVariable lpn_variable01 = LPNVariable(name01, type01, value01);
	
	if (name01.compare(lpn_variable01.GetName()) != 0) {
		std::cout << "Error. lpn_variable01.name, " << lpn_variable01.GetName() << ", does not match name01, " << name01 << ".\n";
	}
	if (type01.compare(lpn_variable01.GetType()) != 0) {
		std::cout << "Error. lpn_variable01.type, " << lpn_variable01.GetType() << ", does not match type01, " << type01 << ".\n";
	}
	if (value01 != lpn_variable01.GetValue()) {
		std::cout << "Error. lpn_variable01.value, " << lpn_variable01.GetValue() << ", does not match value01, " << value01 << ".\n";
	}
	
	std::string name02 = "b";
	std::string type02 = "type_b";
	double value02 = 1.0;
	PressureVariable pressure_variable01 = PressureVariable(name02, type02, value02);
	FlowVariable flow_variable01 = FlowVariable(name02, type02, value02);

	if (name02.compare(pressure_variable01.GetName()) != 0) {
		std::cout << "Error. pressure_variable01.name, " << pressure_variable01.GetName() << ", does not match name02, " << name02 << ".\n";
	}

	if (type02.compare(pressure_variable01.GetType()) != 0) {
		std::cout << "Error. pressure_variable01.type, " << pressure_variable01.GetType() << ", does not match type02, " << type02 << ".\n";
	}

	if (value02 != pressure_variable01.GetValue()) {
		std::cout << "Error. pressure_variable01.value, " << pressure_variable01.GetValue() << ", does not match value02, " << value02 << ".\n";
	}

	if (name02.compare(flow_variable01.GetName()) != 0) {
		std::cout << "Error. flow_variable01.name, " << flow_variable01.GetName() << ", does not match name02, " << name02 << ".\n";
	}

	if (type02.compare(flow_variable01.GetType()) != 0) {
		std::cout << "Error. flow_variable01.type, " << flow_variable01.GetType() << ", does not match type02, " << type02 << ".\n";
	}

	if (value02 != flow_variable01.GetValue()) {
		std::cout << "Error. flow_variable01.value, " << flow_variable01.GetValue() << ", does not match value02, " << value02 << ".\n";
	}

	std::string name03 = "c";
	std::string type03 = "type_c";
	std::vector<int> flow_directions03 {1};

	std::string name04 = "d";
	std::string type04 = "type_d";
	std::vector<int> flow_directions04 {-1, 1};

	std::string name05 = "e";
	std::string type05 = "type_e";
	std::vector<int> flow_directions05 {-1};

	std::vector<std::string> connecting_block_list03 {name04};
	std::vector<std::string> connecting_block_list04 {name03, name05};
	std::vector<std::string> connecting_block_list05 {name04};

	LPNBlock lpn_block03 = LPNBlock(name03, type03, connecting_block_list03, flow_directions03);
	LPNBlock lpn_block04 = LPNBlock(name04, type04, connecting_block_list04, flow_directions04);
	LPNBlock lpn_block05 = LPNBlock(name05, type05, connecting_block_list05, flow_directions05);

	if (lpn_block03.GetName().compare(name03) != 0) {
		std::cout << "Error. lpn_block03.name, " << lpn_block03.GetName() << ", does not match name03, " << name03 << ".\n";
	}

	if (lpn_block03.GetType().compare(type03) != 0) {
		std::cout << "Error. lpn_block03.type, " << lpn_block03.GetType() << ", does not match type03, " << type03 << ".\n";
	}

	if (lpn_block03.GetConnectingBlockList().size() != connecting_block_list03.size()) {
		std::cout << "Error. lpn_block03.connecting_block_list.size, " << lpn_block03.GetConnectingBlockList().size() << ", does not match connecting_block_list03.size, " << connecting_block_list03.size() << ".\n";
	} else {
		for (int i = 0; i < connecting_block_list03.size(); i++) {
			if (lpn_block03.GetConnectingBlockList()[i].compare(connecting_block_list03[i]) != 0) {
				std::cout << "Error. lpn_block03.connecting_block_list[i], " << lpn_block03.GetConnectingBlockList()[i] << ", does not match connecting_block_list03[i], " << connecting_block_list03[i] << ".\n";
			}
		}
	}

	if (lpn_block03.GetFlowDirections().size() != flow_directions03.size()) {
		std::cout << "Error. lpn_block03.flow_directions.size, " << lpn_block03.GetFlowDirections().size() << ", does not match flow_directions03.size, " << flow_directions03.size() << ".\n";
	} else {
		for (int i = 0; i < flow_directions03.size(); i++) {
			if (lpn_block03.GetFlowDirections()[i] != flow_directions03[i]) {
				std::cout << "Error. lpn_block03.flow_directions[i], " << lpn_block03.GetFlowDirections()[i] << ", does not match flow_directions03[i], " << flow_directions03[i] << ".\n";
			}
		}
	}

	if (lpn_block04.GetName().compare(name04) != 0) {
		std::cout << "Error. lpn_block04.name, " << lpn_block04.GetName() << ", does not match name04, " << name04 << ".\n";
	}

	if (lpn_block04.GetType().compare(type04) != 0) {
		std::cout << "Error. lpn_block04.type, " << lpn_block04.GetType() << ", does not match type04, " << type04 << ".\n";
	}

	if (lpn_block04.GetConnectingBlockList().size() != connecting_block_list04.size()) {
		std::cout << "Error. lpn_block04.connecting_block_list.size, " << lpn_block04.GetConnectingBlockList().size() << ", does not match connecting_block_list04.size, " << connecting_block_list04.size() << ".\n";
	} else {
		for (int i = 0; i < connecting_block_list04.size(); i++) {
			if (lpn_block04.GetConnectingBlockList()[i].compare(connecting_block_list04[i]) != 0) {
				std::cout << "Error. lpn_block04.connecting_block_list[i], " << lpn_block04.GetConnectingBlockList()[i] << ", does not match connecting_block_list04[i], " << connecting_block_list04[i] << ".\n";
			}
		}
	}

	if (lpn_block04.GetFlowDirections().size() != flow_directions04.size()) {
		std::cout << "Error. lpn_block04.flow_directions.size, " << lpn_block04.GetFlowDirections().size() << ", does not match flow_directions04.size, " << flow_directions04.size() << ".\n";
	} else {
		for (int i = 0; i < flow_directions04.size(); i++) {
			if (lpn_block04.GetFlowDirections()[i] != flow_directions04[i]) {
				std::cout << "Error. lpn_block04.flow_directions[i], " << lpn_block04.GetFlowDirections()[i] << ", does not match flow_directions04[i], " << flow_directions04[i] << ".\n";
			}
		}
	}

	if (lpn_block05.GetName().compare(name05) != 0) {
		std::cout << "Error. lpn_block05.name, " << lpn_block05.GetName() << ", does not match name05, " << name05 << ".\n";
	}

	if (lpn_block05.GetType().compare(type05) != 0) {
		std::cout << "Error. lpn_block05.type, " << lpn_block05.GetType() << ", does not match type05, " << type05 << ".\n";
	}

	if (lpn_block05.GetConnectingBlockList().size() != connecting_block_list05.size()) {
		std::cout << "Error. lpn_block05.connecting_block_list.size, " << lpn_block05.GetConnectingBlockList().size() << ", does not match connecting_block_list05.size, " << connecting_block_list05.size() << ".\n";
	} else {
		for (int i = 0; i < connecting_block_list05.size(); i++) {
			if (lpn_block05.GetConnectingBlockList()[i].compare(connecting_block_list05[i]) != 0) {
				std::cout << "Error. lpn_block05.connecting_block_list[i], " << lpn_block05.GetConnectingBlockList()[i] << ", does not match connecting_block_list05[i], " << connecting_block_list05[i] << ".\n";
			}
		}
	}

	if (lpn_block05.GetFlowDirections().size() != flow_directions05.size()) {
		std::cout << "Error. lpn_block05.flow_directions.size, " << lpn_block05.GetFlowDirections().size() << ", does not match flow_directions05.size, " << flow_directions05.size() << ".\n";
	} else {
		for (int i = 0; i < flow_directions05.size(); i++) {
			if (lpn_block05.GetFlowDirections()[i] != flow_directions05[i]) {
				std::cout << "Error. lpn_block05.flow_directions[i], " << lpn_block05.GetFlowDirections()[i] << ", does not match flow_directions05[i], " << flow_directions05[i] << ".\n";
			}
		}
	}

	std::array<LPNBlock *, 2> connecting_block_list06 {&lpn_block03, &lpn_block04};
	std::string name06 = "wire06";
	std::array<int, 2> lpn_solution_ids06 {1, 2};
	Wire wire06 = Wire(name06, connecting_block_list06);
	wire06.SetP(pressure_variable01);
	wire06.SetQ(flow_variable01);
	wire06.SetLPNSolutionIds(lpn_solution_ids06);

	if (wire06.GetName().compare(name06) != 0) {
		std::cout << "Error. wire06.name, " << wire06.GetName() << ", does not match name06, " << name06 << ".\n";
	}

	if (wire06.GetConnectingBlockList().size() != connecting_block_list06.size()) {
		std::cout << "Error. wire06.connecting_block_list.size, " << wire06.GetConnectingBlockList().size() << ", does not match connecting_block_list06.size, " << connecting_block_list06.size() << ".\n";
	} else {
		for (int i = 0; i < connecting_block_list06.size(); i++) {
			if ((*(wire06.GetConnectingBlockList()[i])).GetName().compare((*(connecting_block_list06[i])).GetName()) != 0) {
				std::cout << "Error. lpn_block06.connecting_block_list[i].name, " << (*(wire06.GetConnectingBlockList()[i])).GetName() << ", does not match connecting_block_list06[i].name, " << (*(connecting_block_list06[i])).GetName() << ".\n";
			}
		}
	}

	if ((wire06.GetP()).GetName() != pressure_variable01.GetName()) {
		std::cout << "Error. wire06.P.name, " << (wire06.GetP()).GetName() << ", does not match pressure_variable01.name, " << pressure_variable01.GetName() << ".\n";
	}

	if ((wire06.GetQ()).GetName() != flow_variable01.GetName()) {
		std::cout << "Error. wire06.Q.name, " << (wire06.GetQ()).GetName() << ", does not match flow_variable01.name, " << flow_variable01.GetName() << ".\n";
	}

	if (wire06.GetLPNSolutionIds().size() != lpn_solution_ids06.size()) {
		std::cout << "Error. wire06.lpn_solution_ids.size, " << wire06.GetLPNSolutionIds().size() << ", does not match lpn_solution_ids.size, " << lpn_solution_ids06.size() << ".\n";
	} else {
		for (int i = 0; i < lpn_solution_ids06.size(); i++) {
			if (wire06.GetLPNSolutionIds()[i] != lpn_solution_ids06[i]) {
				std::cout << "Error. wire06.lpn_solution_ids[i], " << wire06.GetLPNSolutionIds()[i] << ", does not match lpn_solution_ids06[i], " << lpn_solution_ids06[i] << ".\n";
			}
		}
	}

	return 0;
}
