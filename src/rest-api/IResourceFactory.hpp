#pragma once

#include <memory>
#include <restbed>

using namespace std;

class IResourceFactory {
    
    public:
    virtual shared_ptr<restbed::Resource> get_resource() const = 0;
};