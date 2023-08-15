#include "./SolverResourceFactory.hpp"

#include <sstream>
#include <iomanip>
#include "../../../json.hpp"

SolverResourceFactory::SolverResourceFactory() {
    _resource = make_shared<restbed::Resource>();
    _resource->set_path(
        "/{operation: add|subtract|multiply|divide}"
        "/{num1: [-+]?[0-9]*\\.?[0-9]*}"
        "/{num2: [-+]?[0-9]*\\.?[0-9]*}"
    );
    _resource->set_method_handler("GET", 
        [&](const shared_ptr<restbed::Session> session) {
            get_handler(session);
        });

}

float SolverResourceFactory::calculate(float num1, float num2, string operation) {
    if(operation == "add") {
        return num1 + num2;
    }
    else if(operation == "subtract") {
        return num1 - num2;
    }
    else if(operation == "multiply") {
        return num1 * num2;
    }
    else if(operation == "divide") {
        return num1 / num2;
    }
    return 0;
}

string SolverResourceFactory::to_json(float result) {
    ostringstream str_stream;
    str_stream << result;
    nlohmann::json jsonResult = {
        {"result", str_stream.str()}
    };
    return jsonResult.dump();
}

tuple<float, float, string> SolverResourceFactory::get_path_parameters 
(const shared_ptr<restbed::Session> session) {

    const auto& request = session->get_request();
    const auto operation = request->get_path_parameter("operation");
    auto num1 = atof(request->get_path_parameter("num1").c_str());
    auto num2 = atof(request->get_path_parameter("num2").c_str());
    return make_tuple(num1, num2, operation);
}

shared_ptr<restbed::Resource> SolverResourceFactory::get_resource() const {
        return _resource;
}

void SolverResourceFactory::get_handler(const shared_ptr<restbed::Session> session) {
    const auto [num1, num2, operation] = get_path_parameters(session);
    auto result = calculate(num1, num2, operation);
    auto content = to_json(result);
    session->close(restbed::OK, content, 
        {{"Content-Length", to_string(content.size())}});
}