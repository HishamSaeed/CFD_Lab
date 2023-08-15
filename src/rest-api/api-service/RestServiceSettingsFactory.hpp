#pragma once

#include "../IServiceSettingsFactory.hpp"

class RestServiceSettingsFactory : public IServiceSettingsFactory {

public:

    RestServiceSettingsFactory();
    shared_ptr<restbed::Settings> get_settings() const final;

private:

    std::shared_ptr<restbed::Settings> _settings;

};