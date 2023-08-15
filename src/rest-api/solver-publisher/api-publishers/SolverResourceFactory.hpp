#pragma once

#include "../../IResourceFactory.hpp"
#include <string>
#include <tuple>

class SolverResourceFactory : public IResourceFactory {
    public: 
        SolverResourceFactory();
        std::shared_ptr<restbed::Resource> get_resource() const final;
    
    private:

        float calculate(float num1, float num2, string operation);
        tuple<float, float, string> get_path_parameters(const shared_ptr<restbed::Session> session);
        string to_json(float result);
        void get_handler(const shared_ptr<restbed::Session>);
        
        shared_ptr<restbed::Resource> _resource;

};