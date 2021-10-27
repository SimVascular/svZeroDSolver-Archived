#include <string>
#include <iostream>
#include "blocks.hpp"

int main() {
	std::string name01 = "a";
	std::string type01 = "type_a";
	double value01 = 0.12;
	LPNVariable lpn_variable01 = LPNVariable(name01, type01, value01);
	std::cout << "name  = " << lpn_variable01.GetName()  << "\n";
	std::cout << "type  = " << lpn_variable01.GetType()  << "\n";
	std::cout << "value = " << lpn_variable01.GetValue() << "\n";
	
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
	return 0;
}
