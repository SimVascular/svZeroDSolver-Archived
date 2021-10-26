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
	
    return 0;
}
