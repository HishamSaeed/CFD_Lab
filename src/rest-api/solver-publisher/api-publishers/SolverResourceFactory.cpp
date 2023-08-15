#include "./SolverResourceFactory.hpp"
#include "../../../cfd-solver/cfd-solver.hpp"
#include <sstream>
#include <iomanip>
#include "../../../json.hpp"

SolverResourceFactory::SolverResourceFactory() {
    _resource = make_shared<restbed::Resource>();
    _resource->set_path(
        "/{operation: solve}"
        "/{num1: [-+]?[0-9]*\\.?[0-9]*}"
    );
    _resource->set_method_handler("GET", 
        [&](const shared_ptr<restbed::Session> session) {
            get_handler(session);
        });

}

std::string SolverResourceFactory::calculate(string problem) {
    char* args[] = {
            "./sim", // First element is typically the program name
            "-c",
        };

    int argn = sizeof(args) / sizeof(args[0]);
    if(problem == "0") {
        args[1] = "-c"; 
        solve_cfd(argn, args);
        return "solved";
    }
    else if(problem == "1") {
        args[1] = "-p"; 
        solve_cfd(argn, args);
        return "solved";
    }
    else if(problem == "2") {
        args[1] = "-k"; 
        solve_cfd(argn, args);
        return "solved";
    }
    else if(problem == "3") {
        args[1] = "-f"; 
        solve_cfd(argn, args);
        return "solved";
    }
    else if(problem == "4") {
        args[1] = "-n"; 
        solve_cfd(argn, args);
        return "solved";
    }
    else if(problem == "5") {
        args[1] = "-ft"; 
        solve_cfd(argn, args);
        return "solved";
    }
    else if(problem == "6") {
        args[1] = "-r"; 
        solve_cfd(argn, args);
        return "solved";
    }
    return "Not Solved";
}

string SolverResourceFactory::to_json(std::string result) {
    ostringstream str_stream;
    str_stream << result;
    nlohmann::json jsonResult = {
        {"result", str_stream.str()}
    };
    return jsonResult.dump();
}

std::string SolverResourceFactory::get_path_parameters 
(const shared_ptr<restbed::Session> session) {
    const auto& request = session->get_request();
    const auto operation = request->get_path_parameter("operation");
    auto num1 = request->get_path_parameter("num1").c_str();
    return num1;
}

shared_ptr<restbed::Resource> SolverResourceFactory::get_resource() const {
        return _resource;
}

void SolverResourceFactory::get_handler(const shared_ptr<restbed::Session> session) {
    const std::string num1 = get_path_parameters(session);
    auto result = calculate(num1);
    auto content = to_json(result);
    session->close(restbed::OK, content, 
        {{"Content-Length", to_string(content.size())}});
}