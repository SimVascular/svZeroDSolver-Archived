#include <string>
#include <iostream>
#include "blocks.hpp"

int main() 
{
    std::string name03 = "c";
    std::vector<int> flow_directions03 {1};

    std::string name04 = "d";
    std::vector<int> flow_directions04 {-1, 1};

    std::string name05 = "e";
    std::vector<int> flow_directions05 {-1};

    std::vector<std::string> connecting_block_list03 {name04};
    std::vector<std::string> connecting_block_list04 {name03, name05};
    std::vector<std::string> connecting_block_list05 {name04};

    LPNBlock lpn_block03 = LPNBlock(name03, connecting_block_list03, flow_directions03);
    LPNBlock lpn_block04 = LPNBlock(name04, connecting_block_list04, flow_directions04);
    LPNBlock lpn_block05 = LPNBlock(name05, connecting_block_list05, flow_directions05);

    if (lpn_block03.get_name().compare(name03) != 0) {
        std::cerr << "Error. lpn_block03.name, " << lpn_block03.get_name() << ", does not match name03, " << name03 << ".\n";
        exit(1);
    }

    if (lpn_block03.get_connecting_block_list().size() != connecting_block_list03.size()) {
        std::cerr << "Error. lpn_block03.connecting_block_list.size, " << lpn_block03.get_connecting_block_list().size() << ", does not match connecting_block_list03.size, " << connecting_block_list03.size() << ".\n";
        exit(1);
    } else {
        for (int i = 0; i < connecting_block_list03.size(); i++) {
            if (lpn_block03.get_connecting_block_list()[i].compare(connecting_block_list03[i]) != 0) {
                std::cerr << "Error. lpn_block03.connecting_block_list[i], " << lpn_block03.get_connecting_block_list()[i] << ", does not match connecting_block_list03[i], " << connecting_block_list03[i] << ".\n";
                exit(1);
            }
        }
    }

    if (lpn_block03.get_flow_directions().size() != flow_directions03.size()) {
        std::cerr << "Error. lpn_block03.flow_directions.size, " << lpn_block03.get_flow_directions().size() << ", does not match flow_directions03.size, " << flow_directions03.size() << ".\n";
        exit(1);
    } else {
        for (int i = 0; i < flow_directions03.size(); i++) {
            if (lpn_block03.get_flow_directions()[i] != flow_directions03[i]) {
                std::cerr << "Error. lpn_block03.flow_directions[i], " << lpn_block03.get_flow_directions()[i] << ", does not match flow_directions03[i], " << flow_directions03[i] << ".\n";
                exit(1);
            }
        }
    }

    if (lpn_block04.get_name().compare(name04) != 0) {
        std::cerr << "Error. lpn_block04.name, " << lpn_block04.get_name() << ", does not match name04, " << name04 << ".\n";
        exit(1);
    }

    if (lpn_block04.get_connecting_block_list().size() != connecting_block_list04.size()) {
        std::cerr << "Error. lpn_block04.connecting_block_list.size, " << lpn_block04.get_connecting_block_list().size() << ", does not match connecting_block_list04.size, " << connecting_block_list04.size() << ".\n";
        exit(1);
    } else {
        for (int i = 0; i < connecting_block_list04.size(); i++) {
            if (lpn_block04.get_connecting_block_list()[i].compare(connecting_block_list04[i]) != 0) {
                std::cerr << "Error. lpn_block04.connecting_block_list[i], " << lpn_block04.get_connecting_block_list()[i] << ", does not match connecting_block_list04[i], " << connecting_block_list04[i] << ".\n";
                exit(1);
            }
        }
    }

    if (lpn_block04.get_flow_directions().size() != flow_directions04.size()) {
        std::cerr << "Error. lpn_block04.flow_directions.size, " << lpn_block04.get_flow_directions().size() << ", does not match flow_directions04.size, " << flow_directions04.size() << ".\n";
        exit(1);
    } else {
        for (int i = 0; i < flow_directions04.size(); i++) {
            if (lpn_block04.get_flow_directions()[i] != flow_directions04[i]) {
                std::cerr << "Error. lpn_block04.flow_directions[i], " << lpn_block04.get_flow_directions()[i] << ", does not match flow_directions04[i], " << flow_directions04[i] << ".\n";
                exit(1);
            }
        }
    }

    if (lpn_block05.get_name().compare(name05) != 0) {
        std::cerr << "Error. lpn_block05.name, " << lpn_block05.get_name() << ", does not match name05, " << name05 << ".\n";
        exit(1);
    }

    if (lpn_block05.get_connecting_block_list().size() != connecting_block_list05.size()) {
        std::cerr << "Error. lpn_block05.connecting_block_list.size, " << lpn_block05.get_connecting_block_list().size() << ", does not match connecting_block_list05.size, " << connecting_block_list05.size() << ".\n";
        exit(1);
    } else {
        for (int i = 0; i < connecting_block_list05.size(); i++) {
            if (lpn_block05.get_connecting_block_list()[i].compare(connecting_block_list05[i]) != 0) {
                std::cerr << "Error. lpn_block05.connecting_block_list[i], " << lpn_block05.get_connecting_block_list()[i] << ", does not match connecting_block_list05[i], " << connecting_block_list05[i] << ".\n";
                exit(1);
            }
        }
    }

    if (lpn_block05.get_flow_directions().size() != flow_directions05.size()) {
        std::cerr << "Error. lpn_block05.flow_directions.size, " << lpn_block05.get_flow_directions().size() << ", does not match flow_directions05.size, " << flow_directions05.size() << ".\n";
        exit(1);
    } else {
        for (int i = 0; i < flow_directions05.size(); i++) {
            if (lpn_block05.get_flow_directions()[i] != flow_directions05[i]) {
                std::cerr << "Error. lpn_block05.flow_directions[i], " << lpn_block05.get_flow_directions()[i] << ", does not match flow_directions05[i], " << flow_directions05[i] << ".\n";
                exit(1);
            }
        }
    }

    std::string name06 = "wire06";
    std::array<int, 2> lpn_solution_ids06 {1, 2};
    Wire wire06 = Wire(name06);
    wire06.set_lpn_solution_ids(lpn_solution_ids06);

    if (wire06.get_name().compare(name06) != 0) {
        std::cerr << "Error. wire06.name, " << wire06.get_name() << ", does not match name06, " << name06 << ".\n";
        exit(1);
    }

    if (wire06.get_lpn_solution_ids().size() != lpn_solution_ids06.size()) {
        std::cerr << "Error. wire06.lpn_solution_ids.size, " << wire06.get_lpn_solution_ids().size() << ", does not match lpn_solution_ids.size, " << lpn_solution_ids06.size() << ".\n";
        exit(1);
    } else {
        for (int i = 0; i < lpn_solution_ids06.size(); i++) {
            if (wire06.get_lpn_solution_ids()[i] != lpn_solution_ids06[i]) {
                std::cerr << "Error. wire06.lpn_solution_ids[i], " << wire06.get_lpn_solution_ids()[i] << ", does not match lpn_solution_ids06[i], " << lpn_solution_ids06[i] << ".\n";
                exit(1);
            }
        }
    }
  
    std::string name07 = "name07";
    std::vector<std::string> connecting_block_list07 {name03};
    std::vector<int> flow_directions07 {1};
    LPNBlock lpn_block07 = LPNBlock(name07, connecting_block_list07, flow_directions07);

    std::string name08 = "wire08";
    Wire wire08 = Wire(name08);

    lpn_block03.add_connecting_block(name07, -1);
    lpn_block03.add_connecting_wire(name08);
    lpn_block07.add_connecting_wire(name08);
  
    if (lpn_block03.get_connecting_block_list().size() != connecting_block_list03.size() + 1) {
        std::cerr << "Error. lpn_block03.connecting_block_list.size, " << lpn_block03.get_connecting_block_list().size() << ", does not match connecting_block_list03.size + 1 = " << connecting_block_list03.size() + 1 << ".\n";
        exit(1);
    } else {
        if (lpn_block03.get_connecting_block_list()[connecting_block_list03.size()].compare(name07) != 0) {
            std::cerr << "Error. lpn_block03.connecting_block_list[connecting_block_list03.size()], " << lpn_block03.get_connecting_block_list()[connecting_block_list03.size()] << ", does not match name07, " << name07 << ".\n";
            exit(1);
        }
    }
  
    if (lpn_block03.get_flow_directions().size() != flow_directions03.size() + 1) {
        std::cerr << "Error. lpn_block03.flow_directions.size, " << lpn_block03.get_flow_directions().size() << ", does not match flow_directions03.size + 1 = " << flow_directions03.size() + 1 << ".\n";
        exit(1);
    } else {
        if (lpn_block03.get_flow_directions()[flow_directions03.size()] != -1) {
            std::cerr << "Error. lpn_block03.flow_directions[flow_directions03.size()], " << lpn_block03.get_flow_directions()[flow_directions03.size()] << ", does not equate to 1.\n";
            exit(1);
        }
    }
    
    if (lpn_block03.get_connecting_wires_list().size() != 1) {
        std::cerr << "Error. lpn_block03.connecting_wire_list.size, " << lpn_block03.get_connecting_wires_list().size() << ", is not 1\n";
        exit(1);
    } else {
        if (lpn_block03.get_connecting_wires_list()[0].compare(name08) != 0) {
            std::cerr << "Error. lpn_block03.get_connecting_wires_list[0], " << lpn_block03.get_connecting_wires_list()[0] << ", does not match name08, " << name08 << ".\n";
            exit(1);
        }
    }
    
    double dt01 = 1.1;
    double rho01 = 2.2;
    bool check_jacobian01 = false;
    std::unordered_map<std::string, Wire *> wire_dict01;
    wire_dict01[name06] = &wire06;
    wire_dict01[name08] = &wire08;
    Args args01 = Args(dt01, rho01, check_jacobian01, wire_dict01);
    
    if (args01.get_dt() - dt01 != 0.0) {
        std::cerr << "args01.dt is not" << args01.get_dt() << ".\n";
        exit(1);
    }
    
    if (args01.get_rho() - rho01 != 0.0) {
        std::cerr << "args01.rho is not" << args01.get_rho() << ".\n";
        exit(1);
    }
    
    if (args01.get_check_jacobian() != false) {
        std::cerr <<  "args01.check_jacobian is not " << args01.get_check_jacobian() << ".\n";
        exit(1);
    }
    
    if (args01.get_wire_dict().size() != wire_dict01.size()) {
        std::cerr << "Error. args01.wire_dict.size, " << args01.get_wire_dict().size() << ", does not match wire_dict01.size, " << wire_dict01.size() << ".\n";
        exit(1);
    }
    
    // last here - make tests for Junction
    std::string junction_name01 = "j01";
    

    return 0;
}