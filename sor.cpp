#include "sor.hpp"
#include "boundary_val.hpp"
#include <cmath>
#include <iomanip>
#include <iostream>

void sor(double omg, double dx, double dy, int imax, int jmax, Grid &grid,
         matrix<double> &RS, double *res) {
  int i, j;
  double rloc;
  double coeff = omg / (2.0 * (1.0 / (dx * dx) + 1.0 / (dy * dy)));
  matrix<double> P;
  grid.pressure(P);

  /* set boundary values */
  for (i = 1; i <= imax; i++) {
    // P.at(i)[0] = P.at(i)[1];
    // P.at(i).at(jmax+1) = P.at(i).at(jmax);
    switch (grid.cell(i, 0).flag() &
            ((1 << 1) | (1 << 2) | (1 << 3) | (1 << 4))) {
    case 1 << 1: // No slip conditions
      if (B_N(grid.cell(i, 0).flag())) {
        P.at(i)[0] = P.at(i)[1];
      }

      break;
    }

    switch (grid.cell(i, jmax + 1).flag() &
            ((1 << 1) | (1 << 2) | (1 << 3) | (1 << 4))) {
    case 1 << 1: // No slip conditions
      if (B_S(grid.cell(i, jmax + 1).flag())) {
        P.at(i).at(jmax + 1) = P.at(i).at(jmax);
      }
      break;
    }
  }
  for (j = 1; j <= jmax; j++) {
    // P[0].at(j) = P[1].at(j);
    // P.at(imax+1).at(j) = P.at(imax).at(j);
    switch (grid.cell(0, j).flag() &
            ((1 << 1) | (1 << 2) | (1 << 3) | (1 << 4))) {
    case 1 << 1: // No slip conditions
      if (B_E(grid.cell(0, j).flag())) {
        //    P[0].at(j) = P[1].at(j);
        P[0].at(j) = P[1].at(j);
      }

      break;
    case (1 << 4):
      P[0].at(j) = P[1].at(j);
      break;
    }

    switch (grid.cell(imax + 1, j).flag() &
            ((1 << 1) | (1 << 2) | (1 << 3) | (1 << 4))) {
    case 1 << 1: // No slip conditions
      if (B_W(grid.cell(imax + 1, j).flag())) {
        P.at(imax + 1).at(j) = P.at(imax).at(j);
      }

      break;
    case 1 << 3:
      P.at(imax + 1).at(j) = 0.0;
      break;
    }
  }

  for (int i = 1; i <= imax; ++i) {
    for (int j = 1; j <= jmax; ++j) {

      if (B_E(grid.cell(i, j).flag()))
        P[i][j] = P[i + 1][j];

      if (B_W(grid.cell(i, j).flag()))
        P[i][j] = P[i - 1][j];

      if (B_N(grid.cell(i, j).flag()))
        P[i][j] = P[i][j + 1];

      if (B_S(grid.cell(i, j).flag()))
        P[i][j] = P[i][j - 1];
      // if ( B_NE(grid.cell(i,j).flag() ) ) P[i][j] = (P[i][j+1] +
      // P[i+1][j])/2;

      // if ( B_NW(grid.cell(i,j).flag() ) ) P[i][j] = (P[i][j+1] +
      // P[i-1][j])/2;

      // if ( B_SE(grid.cell(i,j).flag() ) ) P[i][j] = (P[i][j-1] +
      // P[i+1][j])/2;

      // if ( B_SW(grid.cell(i,j).flag() ) ) P[i][j] = (P[i][j-1] +
      // P[i-1][j])/2;
      if (B_NE(grid.cell(i, j).flag()))
        (P.at(i).at(j) = (P.at(i).at(j + 1) + P.at(i + 1).at(j)) / 2);
      // std::cout << std::fixed << std::setprecision(20) << P.at(i).at(j+1) <<
      // std::endl; std::cout << std::fixed << std::setprecision(20) <<
      // P.at(i+1).at(j) << std::endl; std::cout << (P.at(i).at(j+1) +
      // P.at(i+1).at(j))/2 << std::endl; std::cout <<  P.at(i).at(j) <<
      // std::endl;}

      if (B_NW(grid.cell(i, j).flag()))
        (P.at(i).at(j) = (P.at(i).at(j + 1) + P.at(i - 1).at(j)) / 2);
      // {std::cout << "i = " << i << " j = " << j << std::endl;
      // std::cout <<  P.at(i).at(j) << std::endl;

      // std::cout << std::fixed << std::setprecision(20) << P.at(i).at(j+1) <<
      // std::endl; std::cout << std::fixed << std::setprecision(20) <<
      // P.at(i-1).at(j) << std::endl; std::cout << (P.at(i).at(j+1) +
      // P.at(i-1).at(j))/2 << std::endl; std::cout <<  P.at(i).at(j) <<
      // std::endl;}

      if (B_SE(grid.cell(i, j).flag()))
        P.at(i).at(j) = (P.at(i).at(j - 1) + P.at(i + 1).at(j)) / 2;

      if (B_SW(grid.cell(i, j).flag()))
        P.at(i).at(j) = (P.at(i).at(j - 1) + P.at(i - 1).at(j)) / 2;
    }
  }
  // std::cout <<  P.at(38).at(30) << std::endl;
  /* SOR iteration */
  for (i = 1; i <= imax; i++) {
    for (j = 1; j <= jmax; j++) {
      if (grid.cell(i, j).flag() & (1 << 0)) {
        P.at(i).at(j) =
            (1.0 - omg) * P.at(i).at(j) +
            coeff * ((P.at(i + 1).at(j) + P.at(i - 1).at(j)) / (dx * dx) +
                     (P.at(i).at(j + 1) + P.at(i).at(j - 1)) / (dy * dy) -
                     RS.at(i).at(j));
        //     P[i][j] = (1.0-omg)*P[i][j]
        //   + coeff*( (P[i+1][j]+P[i-1][j])/(dx*dx) + (
        //   P[i][j+1]+P[i][j-1])/(dy*dy) - RS[i][j]);
      }
    }
  }

  /* compute the residual */
  rloc = 0;
  int fluid_cell_counter = 0;
  for (i = 1; i <= imax; i++) {
    for (j = 1; j <= jmax; j++) {
      if (grid.cell(i, j).flag() & (1 << 0)) {
        fluid_cell_counter++;
        rloc += ((P.at(i + 1).at(j) - 2.0 * P.at(i).at(j) + P.at(i - 1).at(j)) /
                     (dx * dx) +
                 (P.at(i).at(j + 1) - 2.0 * P.at(i).at(j) + P.at(i).at(j - 1)) /
                     (dy * dy) -
                 RS.at(i).at(j)) *
                ((P.at(i + 1).at(j) - 2.0 * P.at(i).at(j) + P.at(i - 1).at(j)) /
                     (dx * dx) +
                 (P.at(i).at(j + 1) - 2.0 * P.at(i).at(j) + P.at(i).at(j - 1)) /
                     (dy * dy) -
                 RS.at(i).at(j));

        //     rloc += ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + (
        //     P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j])*
        //   ( (P[i+1][j]-2.0*P[i][j]+P[i-1][j])/(dx*dx) + (
        //   P[i][j+1]-2.0*P[i][j]+P[i][j-1])/(dy*dy) - RS[i][j]);
      }
    }
  }
  // rloc = rloc/(imax*jmax);
  // std::cout << "counter " << fluid_cell_counter << std::endl;
  rloc = rloc / (fluid_cell_counter);
  rloc = sqrt(rloc);
  /* set residual */
  *res = rloc;

  /* set boundary values */
  // for(i = 1; i <= imax; i++) {
  //     // P.at(i)[0] = P.at(i)[1];
  //     // P.at(i).at(jmax+1) = P.at(i).at(jmax);
  //     switch(grid.cell(i,0).flag() & ((1<<1)|(1<<2)|(1<<3)|(1<<4))  ){
  //       case 1<<1://No slip conditions
  //         if ( B_N(grid.cell(i,0).flag()) )
  //         {
  //             P.at(i)[0] = P.at(i)[1];
  //         }

  //         break;
  //     }

  //     switch(grid.cell(i,jmax+1).flag() & ((1<<1)|(1<<2)|(1<<3)|(1<<4))  ){
  //       case 1<<1://No slip conditions
  //         if ( B_S(grid.cell(i,jmax+1).flag()) )
  //         {
  //           P.at(i).at(jmax+1) = P.at(i).at(jmax);
  //         }
  //         break;
  //     }
  // }
  // for(j = 1; j <= jmax; j++) {
  //     // P[0].at(j) = P[1].at(j);
  //     // P.at(imax+1).at(j) = P.at(imax).at(j);
  //     switch(grid.cell(0,j).flag() & ((1<<1)|(1<<2)|(1<<3)|(1<<4))  ){
  //       case 1<<1://No slip conditions
  //         if ( B_E(grid.cell(0,j).flag()) )
  //         {
  //            P[0].at(j) = P[1].at(j);
  //         }

  //           break;
  //       case (1<<4):
  //             P[0].at(j) = P[1].at(j);
  //         break;
  //     }

  //     switch(grid.cell(imax+1,j).flag() & ((1<<1)|(1<<2)|(1<<3)|(1<<4))  ){
  //       case 1<<1://No slip conditions
  //         if ( B_W(grid.cell(imax+1,j).flag()) )
  //         {
  //             P.at(imax+1).at(j) = P.at(imax).at(j);
  //         }

  //           break;
  //       case 1<<3:
  //             P.at(imax+1).at(j) = 0.0;
  //         break;
  //     }
  // }

  // for(int i = 1; i<=imax; ++i)
  //     {
  //     	for(int j = 1; j<=jmax; ++j)
  //     	{

  //     		if ( B_NE(grid.cell(i,j).flag() ) ) P[i][j] = (P[i][j+1]
  //     + P[i+1][j])/2;

  //     		if ( B_NW(grid.cell(i,j).flag() ) ) P[i][j] = (P[i][j+1]
  //     + P[i-1][j])/2;

  //     		if ( B_SE(grid.cell(i,j).flag() ) ) P[i][j] = (P[i][j-1]
  //     + P[i+1][j])/2;

  //     		if ( B_SW(grid.cell(i,j).flag() ) ) P[i][j] = (P[i][j-1]
  //     + P[i-1][j])/2;

  //     	}
  //     }
  grid.set_pressure(P);
}
