#ifndef CFDLAB_CELL_HPP
#define CFDLAB_CELL_HPP
#include <array>
#include "enums.hpp"

class Cell {
public:
    // Constructors
    Cell();
    Cell(double& PI, double& UI, double& VI);

    // Get + Set pressure
    double& pressure();
    void set_pressure(double& value);

    // Get + Set velocity
    double& velocity(velocity_type type);
    void set_velocity(double& value, velocity_type type);

    // Get + Set vorder
    bool& border(border_position position);
    void set_border(border_position position);

    // Get + Set flag
    int& flag();
    void set_flag(int flag);

    // Get + Set temperature
    double& temperature();
    void set_temperature(double& value);


private:
    // one pressure value per call
    double _pressure = 0;

    // one temperature value per cell
    double _temperature = 0;


    // Flag for types of both cell and boundary
    int _flag;

    // Fixed size velocity
    std::array<double, 2> _velocity = {0};
    
    // Fixed number of borders
    std::array<bool, 4> _border = {false};
};

#endif //CFDLAB_CELL_HPP
