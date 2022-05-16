#include "grid.hpp"
#include <algorithm>
#include <iomanip>
#include <iostream>
#include <vector>

Grid::Grid(int imax_init, int jmax_init, int boundary_size_init, double &PI,
           double &UI, double &VI, double &TI, int **geometry)
    : _boundary_size(boundary_size_init),
      _imax_b(imax_init + 2 * boundary_size_init), // +2 for the total size of
                                                   // the vector / matrix
      _jmax_b(jmax_init + 2 * boundary_size_init), _imax(imax_init),
      _jmax(jmax_init) {

  // Resizing the grid cells
  // the boundary size is given as a single value for all borders
  _cells.resize(Grid::imaxb(),
                std::vector<Cell>(Grid::jmaxb(), Cell(PI, UI, VI)));

  // Resize Velocity Matrices
  _velocities[static_cast<int>(velocity_type::U)].resize(
      Grid::imaxb(), std::vector<double>(Grid::jmaxb(), UI));
  _velocities[static_cast<int>(velocity_type::V)].resize(
      Grid::imaxb(), std::vector<double>(Grid::jmaxb(), VI));

  // Do Not confuse between geometry and the bit poistion,fluid cell geometry is
  // 4 but bit 0 has to be set

  /**
   * Flag Setup 9 bit flag
   * bit 0 indicate fluid
   * bit 1 indicate No slip
   * bit 2 indicate Free slip
   * bit 3 indicate OutFlow
   * bit 4 indicate InFlow
   * bit 5 indicate boundary north
   * bit 6 indicate boundary south
   * bit 7 indicate boundary west
   * bit 8 indicate boundary east
   */

  /**Geometry
   * 0 No-Slip
   * 1 Free slip
   * 2 OutFlow
   * 3 InFlow
   * 4 Fluid
   */

  // Setting the flag
  for (int i = 0; i <= imaxb() - 1; i++) {
    for (int j = 0; j <= jmaxb() - 1; j++) {
      cell(i, j).flag() = 0;

      switch (geometry[i][j]) {
      // No slip cell, set bit 1
      case 0:
        cell(i, j).flag() = 1 << 1;
        break;

      // Free slip, set bit 2
      case 1:
        cell(i, j).flag() = 1 << 2;
        break;

      // OutFlow,set bit 2
      case 2:
        cell(i, j).flag() = 1 << 3;
        break;

      // InFlow,set bit 4
      case 3:
        cell(i, j).flag() = 1 << 4;
        break;

      // Fluid,set bit 0
      case 4:
        cell(i, j).flag() = 1 << 0;
        cell(i, j).velocity(velocity_type::U) = UI;
        cell(i, j).velocity(velocity_type::V) = VI;
        cell(i, j).pressure() = PI;
        cell(i, j).temperature() = TI;
        break;

      default:
        break;
      }

      if (!isfluid(geometry[i][j])) // set boundaries if not
      {
        // East boundary set bit 8
        if (i < imax() - 1 && geometry[i + 1][j] == 4) {
          cell(i, j).flag() |= 1 << 8;
        }
        // West boundary set bit 7
        if (i > 0 && geometry[i - 1][j] == 4) {
          cell(i, j).flag() |= 1 << 7;
        }
        // North boundary set bit 5
        if (j < jmax() - 1 && geometry[i][j + 1] == 4) {
          cell(i, j).flag() |= 1 << 5;
        }
        // South boundary set bit 6
        if (j > 0 && geometry[i][j - 1] == 4) {
          cell(i, j).flag() |= 1 << 6;
        }
      }
    }
  }

  // Setting the cells with a boundary cell
  // Top border
  std::vector<Cell> cells_temp;
  Grid::cells(cells_temp, matrix_selection::ROW, 0);
  std::for_each(cells_temp.begin(), cells_temp.end(),
                [](Cell &cell) { cell.set_border(border_position::BOTTOM); });
  Grid::set_cells(cells_temp, matrix_selection::ROW, 0);

  // Bottom
  Grid::cells(cells_temp, matrix_selection::ROW, jmaxb() - 1);
  std::for_each(cells_temp.begin(), cells_temp.end(),
                [](Cell &cell) { cell.set_border(border_position::TOP); });
  Grid::set_cells(cells_temp, matrix_selection::ROW, jmaxb() - 1);

  // Left
  Grid::cells(cells_temp, matrix_selection::COLUMN, 0);
  std::for_each(cells_temp.begin(), cells_temp.end(),
                [](Cell &cell) { cell.set_border(border_position::LEFT); });
  Grid::set_cells(cells_temp, matrix_selection::COLUMN, 0);

  // RIGHT
  Grid::cells(cells_temp, matrix_selection::COLUMN, Grid::imaxb() - 1);
  std::for_each(cells_temp.begin(), cells_temp.end(),
                [](Cell &cell) { cell.set_border(border_position::RIGHT); });
  Grid::set_cells(cells_temp, matrix_selection::COLUMN, Grid::imaxb() - 1);
};

int Grid::jmaxb() const { return _jmax_b; };

int Grid::imaxb() const { return _imax_b; };

int Grid::imax() const { return _imax; };

int Grid::jmax() const { return _jmax; };

int Grid::boundary_size() const { return _boundary_size; };

Cell &Grid::cell(int i, int j) { return _cells.at(i).at(j); };

void Grid::set_cell(Cell &cell, int i, int j) { _cells.at(i).at(j) = cell; };

void Grid::velocity(matrix<double> &vec, velocity_type type) {
  if (!(type == velocity_type::U || type == velocity_type::V)) {
    std::cerr << "Wrong velocity type" << std::endl;
  } else {
    // Resize vector and set all values to 0
    vec.resize(Grid::imaxb(), std::vector<double>(Grid::jmaxb(), 0));
    for (int x = 0; x < Grid::imaxb(); x++) {
      for (int y = 0; y < Grid::jmaxb(); y++) {
        // Accessing velocity
        vec.at(x).at(y) = _cells.at(x).at(y).velocity(type);
      }
    }
  }
}

void Grid::set_velocity(matrix<double> &vec, velocity_type type) {
  if (!(type == velocity_type::U || type == velocity_type::V)) {
    std::cerr << "Wrong velocity type" << std::endl;
  } else {
    for (int x = 0; x < Grid::imaxb(); x++) {
      for (int y = 0; y < Grid::jmaxb(); y++) {
        // Accessing velocity
        _cells.at(x).at(y).set_velocity(vec.at(x).at(y), type);
      }
    }
  }
}

void Grid::pressure(matrix<double> &vec) {
  // Resize vector and set all values to 0
  vec.resize(Grid::imaxb(), std::vector<double>(Grid::jmaxb(), 0));

  // Iterate over cells
  for (int x = 0; x < Grid::imaxb(); x++) {
    for (int y = 0; y < Grid::jmaxb(); y++) {
      // Accessing velocity
      vec.at(x).at(y) = _cells.at(x).at(y).pressure();
    }
  }
}

void Grid::set_pressure(matrix<double> &vec) {

  // Iterate over cells
  for (int x = 0; x < Grid::imaxb(); x++) {
    for (int y = 0; y < Grid::jmaxb(); y++) {
      // Accessing velocity
      _cells.at(x).at(y).set_pressure(vec.at(x).at(y));
    }
  }
}

void Grid::temperature(matrix<double> &vec) {
  // Resize vector and set all values to 0
  vec.resize(Grid::imaxb(), std::vector<double>(Grid::jmaxb(), 0));

  // Iterate over cells
  for (int x = 0; x < Grid::imaxb(); x++) {
    for (int y = 0; y < Grid::jmaxb(); y++) {
      // Accessing velocity
      vec.at(x).at(y) = _cells.at(x).at(y).temperature();
    }
  }
}

void Grid::set_temperature(matrix<double> &vec) {
  // Iterate over cells
  for (int x = 0; x < Grid::imaxb(); x++) {
    for (int y = 0; y < Grid::jmaxb(); y++) {
      // Accessing velocity
      _cells.at(x).at(y).set_temperature(vec.at(x).at(y));
    }
  }
}

void Grid::cells(std::vector<Cell> &cells, matrix_selection m, int index) {
  // Resizing cells vector
  switch (m) {
  case matrix_selection::ROW:
    cells.resize(Grid::imaxb());
    for (int i = 0; i < Grid::imaxb(); i++) {
      cells.at(i) = Grid::cell(i, index);
    }
    break;
  case matrix_selection::COLUMN:
    cells.resize(Grid::jmaxb());
    for (int i = 0; i < Grid::jmaxb(); i++) {
      cells.at(i) = Grid::cell(index, i);
    }
    break;
  default:
    break;
  }
}

void Grid::innercells(std::vector<std::vector<Cell>> &cells) {
  cells.resize(Grid::imax(), std::vector<Cell>(Grid::jmax(), Cell()));
  for (int y = 0; y < Grid::jmax(); y++) {
    for (int x = 0; x < Grid::imax(); x++) {
      // Accessing velocity
      cells.at(x).at(y) = Grid::_cells.at(x + Grid::_boundary_size)
                              .at(y + Grid::_boundary_size);
    }
  }
}

void Grid::set_innercells(std::vector<std::vector<Cell>> &cells) {
  for (int y = 0; y < Grid::jmax(); y++) {
    for (int x = 0; x < Grid::imax(); x++) {
      // Accessing velocity
      Grid::_cells.at(x + Grid::_boundary_size).at(y + Grid::_boundary_size) =
          cells.at(x).at(y);
    }
  }
}

void Grid::set_cells(std::vector<Cell> &cells, matrix_selection m, int index) {
  // Resizing cells vector
  switch (m) {
  case matrix_selection::ROW:
    cells.resize(Grid::imaxb());
    for (int i = 0; i < Grid::imaxb(); i++) {
      Grid::set_cell(cells.at(i), i, index);
    }
    break;
  case matrix_selection::COLUMN:
    cells.resize(Grid::jmaxb());
    for (int i = 0; i < Grid::jmaxb(); i++) {
      Grid::set_cell(cells.at(i), index, i);
    }
    break;
  default:
    break;
  }
}

void Grid::print_velocity(velocity_type type) {
  if (!(type == velocity_type::U || type == velocity_type::V)) {
    std::cerr << "Wrong velocity type" << std::endl;
  } else {
    matrix<double> vel;
    Grid::velocity(vel, type);
    for (int y = Grid::jmaxb() - 1; y >= 0; y--) {
      for (int x = 0; x < Grid::imaxb(); x++) {
        // Accessing velocity
        std::cout << std::right << std::setw(3);
        std::cout << std::fixed << std::setprecision(7) << vel.at(x).at(y)
                  << " ";
      }
      // Print new line
      std::cout << std::endl;
    }
  }
}

void Grid::print_pressure() {
  for (int y = Grid::jmaxb() - 1; y >= 0; y--) {
    for (int x = 0; x < Grid::imaxb(); x++) {
      // Accessing pressure
      std::cout << std::fixed << std::setprecision(7) << std::right
                << std::setw(2) << Grid::_cells.at(x).at(y).pressure() << " ";
    }
    // Print new line
    std::cout << std::endl;
  }
}

void Grid::print_flags() {
  for (int y = Grid::jmaxb() - 1; y >= 0; y--) {
    for (int x = 0; x < Grid::imaxb(); x++) {
      // Accessing flag
      std::cout << Grid::_cells.at(x).at(y).flag() << " ";
    }
    // Print new line
    std::cout << std::endl;
  }
}

void Grid::print_temperature() {
  for (int y = Grid::jmaxb() - 1; y >= 0; y--) {
    for (int x = 0; x < Grid::imaxb(); x++) {

      // Accessing pressure
      std::cout << std::fixed << std::setprecision(0) << std::right
                << std::setw(2) << Grid::_cells.at(x).at(y).temperature()
                << " ";
    }
    // Print new line
    std::cout << std::endl;
  }
}

int Grid::isfluid(int cell) {
  // check if it is inflow,outflow or fluid
  if ((cell == 2) || (cell == 3) || (cell == 4)) {
    return 1;
  } else {
    return 0;
  }
}
