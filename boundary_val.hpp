#ifndef __RANDWERTE_HPP__
#define __RANDWERTE_HPP__

#include "cstring"
#include "helper.hpp"
#include "datastructures.hpp"
#include "grid.hpp"

/**
 * The boundary values of the problem are set.
 */
void boundaryvalues_d(int imax, int jmax, Grid& grid); 

void boundaryvalues(int imax, int jmax, Grid& grid,double UT);

int B_N(int flag);

int B_S(int flag);

int B_W(int flag);

int B_E(int flag);

int B_NE(int flag);

int B_NW(int flag);

int B_SW(int flag);

int B_SE(int flag);

#endif
