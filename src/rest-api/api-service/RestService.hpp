#pragma once

#include "../IService.hpp"

#include "../IResourceFactory.hpp"
#include "../IServiceSettingsFactory.hpp"

class RestService : public IService {

public:

    RestService(
        shared_ptr<IResourceFactory> resource_factory, 
        shared_ptr<IServiceSettingsFactory> settings_factory);
    void start() final;

private:

    restbed::Service _service;
    shared_ptr<IServiceSettingsFactory> _settings_factory;

};