#include "boundary_val.hpp"
#include "datastructures.hpp"
#include "grid.hpp"

// Deprecated
void boundaryvalues_d(int imax, int jmax, Grid &grid) {

  // Setting the boundary values for the vertical boundaries

  // Setting the inner cells boundary conditions First
  // right Wall boundary conditions
  for (int j = 1; j <= jmax; j++) {
    // Horizontal Velocities
    grid.cell(imax, j).velocity(velocity_type::U) = 0.0;
  }

  for (int i = 1; i <= imax; i++) {
    // Vertical velocities
    grid.cell(i, jmax).velocity(velocity_type::V) = 0.0;
  }

  // Left and Right wall
  for (int j = 1; j <= jmax; j++) {
    // left Wall
    {
      // Vertical Velocities
      grid.cell(0, j).velocity(velocity_type::V) =
          (-1.0) * (grid.cell(1, j).velocity(velocity_type::V));

      // Horizontal Velocities
      grid.cell(0, j).velocity(velocity_type::U) = 0.0;
    }
    // right wall
    {
      // Vertical Velocities
      grid.cell(imax + 1, j).velocity(velocity_type::V) =
          (-1.0) * grid.cell(imax, j).velocity(velocity_type::V);
    }
  }
  // Top and Bottom wall
  for (int i = 1; i <= imax; i++) {
    // bottom wall
    {
      // Vertical velocities
      grid.cell(i, 0).velocity(velocity_type::V) = 0.0;

      // Horizontal Velocities
      grid.cell(i, 0).velocity(velocity_type::U) =
          (-1.0) * grid.cell(i, 1).velocity(velocity_type::U);
    }
    // moving wall
    {
      // Horizontal Velocitiy
      grid.cell(i, jmax + 1).velocity(velocity_type::U) =
          2.0 - grid.cell(i, jmax).velocity(velocity_type::U);
    }
  }
}

void boundaryvalues(int imax, int jmax, Grid &grid, double UT) {

  // Boundary Values that adjust inner cells Right Boundary
  for (int j = 1; j <= jmax; j++) {
    // Horizontal Velocities
    switch (grid.cell(imax + 1, j).flag() &
            ((1 << 1) | (1 << 2) | (1 << 3) | (1 << 4))) {
    case 1 << 1: // No slip conditions
      if (B_W(grid.cell(imax + 1, j).flag())) {
        grid.cell(imax, j).velocity(velocity_type::U) = 0;
      }

      break;
    }
  }
  // Boundary Values that adjust inner cells Top Boundary
  for (int i = 1; i <= imax; i++) {
    // Vertical velocities
    switch (grid.cell(i, jmax + 1).flag() &
            ((1 << 1) | (1 << 2) | (1 << 3) | (1 << 4))) {
    case 1 << 1: // No slip conditions
      if (B_S(grid.cell(i, jmax + 1).flag())) {
        grid.cell(i, jmax).velocity(velocity_type::V) = 0;
      }
      break;
    }
  }

  // Left and Right Wall
  for (int j = 1; j <= jmax; j++) {
    switch (grid.cell(0, j).flag() &
            ((1 << 1) | (1 << 2) | (1 << 3) | (1 << 4))) {
    case 1 << 1: // No slip conditions
      if (B_E(grid.cell(0, j).flag())) {
        grid.cell(0, j).velocity(velocity_type::U) = 0;
        grid.cell(0, j).velocity(velocity_type::V) =
            -(1.0) * grid.cell(1, j).velocity(velocity_type::V);
        // grid.cell(0,j).temperature() = 2* 1.0 - grid.cell(1,j).temperature();
      }

      break;
    case (1 << 4):
      grid.cell(0, j).velocity(velocity_type::U) = 1.0;
      grid.cell(0, j).velocity(velocity_type::V) = 0.0;
      break;
    }

    switch (grid.cell(imax + 1, j).flag() &
            ((1 << 1) | (1 << 2) | (1 << 3) | (1 << 4))) {
    case 1 << 1: // No slip conditions
      if (B_W(grid.cell(imax + 1, j).flag())) {
        grid.cell(imax, j).velocity(velocity_type::U) = 0;
        grid.cell(imax + 1, j).velocity(velocity_type::V) =
            -(1.0) * grid.cell(imax, j).velocity(velocity_type::V);
        // grid.cell(imax+1,j).temperature() = 2 * 0.0 -
        // grid.cell(imax,j).temperature();
      }

      break;
    case 1 << 3:
      grid.cell(imax + 1, j).velocity(velocity_type::U) =
          grid.cell(imax, j).velocity(velocity_type::U);
      ;
      grid.cell(imax + 1, j).velocity(velocity_type::V) =
          grid.cell(imax, j).velocity(velocity_type::V);
      break;
    }
  }
  // Top and Bottom Wall
  for (int i = 1; i <= imax; i++) {
    switch (grid.cell(i, 0).flag() &
            ((1 << 1) | (1 << 2) | (1 << 3) | (1 << 4))) {
    case 1 << 1: // No slip conditions
      if (B_N(grid.cell(i, 0).flag())) {
        grid.cell(i, 0).velocity(velocity_type::V) = 0;
        grid.cell(i, 0).velocity(velocity_type::U) =
            -(1.0) * grid.cell(i, 1).velocity(velocity_type::U);
        grid.cell(i, 0).temperature() = grid.cell(i, 1).temperature();
      }

      break;
    }

    switch (grid.cell(i, jmax + 1).flag() &
            ((1 << 1) | (1 << 2) | (1 << 3) | (1 << 4))) {
    case 1 << 1: // No slip conditions
      if (B_S(grid.cell(i, jmax + 1).flag())) {
        grid.cell(i, jmax).velocity(velocity_type::V) = 0;
        grid.cell(i, jmax + 1).velocity(velocity_type::U) =
            2 * UT - (1.0) * grid.cell(i, jmax).velocity(velocity_type::U);
        grid.cell(i, jmax + 1).temperature() = grid.cell(i, jmax).temperature();
      }
      break;
    }
  }

// Setting Corners
#pragma region
  // Left Bottom Corner
  grid.cell(0, 0).velocity(velocity_type::U) =
      -(1.0) * grid.cell(0, 1).velocity(velocity_type::U);
  grid.cell(0, 0).velocity(velocity_type::V) = 0.0;
  // Right Bottom Corner
  grid.cell(imax + 1, 0).velocity(velocity_type::U) =
      -(1.0) * grid.cell(imax + 1, 1).velocity(velocity_type::U);
  grid.cell(imax + 1, 0).velocity(velocity_type::V) = 0.0;
  // Left Top Corner
  grid.cell(0, jmax + 1).velocity(velocity_type::U) =
      -(1.0) * grid.cell(0, jmax).velocity(velocity_type::U);
  grid.cell(0, jmax + 1).velocity(velocity_type::V) = 0.0;
  // Right Top Corner
  grid.cell(imax + 1, jmax + 1).velocity(velocity_type::U) =
      -(1.0) * grid.cell(imax + 1, jmax).velocity(velocity_type::U);
  grid.cell(imax + 1, jmax + 1).velocity(velocity_type::V) = 0.0;
#pragma endregion

  for (int i = imax; i >= 0; i--) {
    for (int j = jmax; j >= 0; j--) {
      if (B_NE(grid.cell(i, j).flag())) {
        grid.cell(i, j).velocity(velocity_type::U) = 0;
        grid.cell(i - 1, j).velocity(velocity_type::U) =
            -(1.0) * grid.cell(i - 1, j + 1).velocity(velocity_type::U);
        grid.cell(i, j).velocity(velocity_type::V) = 0;
        grid.cell(i, j - 1).velocity(velocity_type::V) =
            -(1.0) * grid.cell(i + 1, j - 1).velocity(velocity_type::V);
      }

      if (B_NW(grid.cell(i, j).flag())) {
        grid.cell(i - 1, j).velocity(velocity_type::U) = 0;
        grid.cell(i, j).velocity(velocity_type::U) =
            -(1.0) * grid.cell(i, j + 1).velocity(velocity_type::U);
        grid.cell(i, j).velocity(velocity_type::V) = 0;
        grid.cell(i, j - 1).velocity(velocity_type::V) =
            -(1.0) * grid.cell(i - 1, j - 1).velocity(velocity_type::V);
      }
      if (B_SW(grid.cell(i, j).flag())) {
        grid.cell(i - 1, j).velocity(velocity_type::U) = 0;
        grid.cell(i, j).velocity(velocity_type::U) =
            -(1.0) * grid.cell(i, j - 1).velocity(velocity_type::U);
        grid.cell(i, j - 1).velocity(velocity_type::V) = 0;
        grid.cell(i, j).velocity(velocity_type::V) =
            -(1.0) * grid.cell(i - 1, j).velocity(velocity_type::V);
      }
      if (B_SE(grid.cell(i, j).flag())) {
        grid.cell(i, j).velocity(velocity_type::U) = 0;
        grid.cell(i - 1, j).velocity(velocity_type::U) =
            -(1.0) * grid.cell(i - 1, j - 1).velocity(velocity_type::U);
        grid.cell(i, j - 1).velocity(velocity_type::V) = 0;
        grid.cell(i, j).velocity(velocity_type::V) =
            -(1.0) * grid.cell(i + 1, j).velocity(velocity_type::V);
      }
    }
  }

  // for(int i = 0; i < imax; i++){
  //   for(int j = 0;j < jmax; j++){
  //     switch(grid.cell(i,j).flag() & ((1<<1)|(1<<2)|(1<<3))  ){
  //       case 1<<1://No slip conditions

  //         if ( B_N(grid.cell(i,j).flag()) )
  //           {
  //             grid.cell(i,j).velocity(velocity_type::V) = 0;
  //             grid.cell(i-1,j).velocity(velocity_type::U) = -(1.0) *
  //             grid.cell(i-1,j+1).velocity(velocity_type::U);
  //             grid.cell(i,j).velocity(velocity_type::U) = -(1.0) *
  //             grid.cell(i,j+1).velocity(velocity_type::U) ;
  //           }

  //         if ( B_S(grid.cell(i,j).flag()) )
  //           {
  //               grid.cell(i,j-1).velocity(velocity_type::V) = 0;
  //               grid.cell(i-1,j).velocity(velocity_type::U) = 2-(1.0) *
  //               grid.cell(i-1,j-1).velocity(velocity_type::U);
  //               grid.cell(i,j).velocity(velocity_type::U) = 2-(1.0) *
  //               grid.cell(i,j-1).velocity(velocity_type::U) ;
  //           }

  //         if ( B_E(grid.cell(i,j).flag()) )
  //         {
  //           grid.cell(i,j).velocity(velocity_type::U) = 0;
  //           grid.cell(i,j-1).velocity(velocity_type::V) = -(1.0) *
  //           grid.cell(i+1,j-1).velocity(velocity_type::V);
  //           grid.cell(i,j).velocity(velocity_type::V) = -(1.0) *
  //           grid.cell(i+1,j).velocity(velocity_type::V) ;

  //         }

  //         if ( B_W(grid.cell(i,j).flag()) )
  //             {
  //               grid.cell(i-1,j).velocity(velocity_type::U) = 0;
  //               grid.cell(i,j-1).velocity(velocity_type::V) = -(1.0) *
  //               grid.cell(i-1,j-1).velocity(velocity_type::V);
  //               grid.cell(i,j).velocity(velocity_type::V) = -(1.0) *
  //               grid.cell(i-1,j).velocity(velocity_type::V) ;
  //             // U[i-1][j] = 0;
  //             // V[i][j-1] = -V[i-1][j-1];
  //             // V[i][j] = -V[i-1][j];
  //             }

  //         if ( B_NE(grid.cell(i,j).flag())  )
  //             {
  //               grid.cell(i,j).velocity(velocity_type::U) = 0;
  //               grid.cell(i-1,j).velocity(velocity_type::U) = -(1.0) *
  //               grid.cell(i-1,j+1).velocity(velocity_type::U);
  //               grid.cell(i,j).velocity(velocity_type::V) = 0;
  //               grid.cell(i,j-1).velocity(velocity_type::V) = -(1.0) *
  //               grid.cell(i+1,j-1).velocity(velocity_type::V) ;
  //             // U[i][j] = 0;
  //             // U[i-1][j] = -U[i-1][j+1];
  //             // V[i][j] = 0;
  //             // V[i][j-1] = -V[i+1][j-1];
  //             }

  //         if ( B_NW(grid.cell(i,j).flag())  )
  //           {
  //               grid.cell(i-1,j).velocity(velocity_type::U) = 0;
  //               grid.cell(i,j).velocity(velocity_type::U) = -(1.0) *
  //               grid.cell(i,j+1).velocity(velocity_type::U);
  //               grid.cell(i,j).velocity(velocity_type::V) = 0;
  //               grid.cell(i,j-1).velocity(velocity_type::V) = -(1.0) *
  //               grid.cell(i-1,j-1).velocity(velocity_type::V) ;
  //           // U[i-1][j] = 0;
  //           // U[i][j] = - U[i][j+1];
  //           // V[i][j] = 0;
  //           // V[i][j-1] = -V[i-1][j-1];
  //           }

  //         if ( B_SE(grid.cell(i,j).flag()) ){
  //               grid.cell(i,j).velocity(velocity_type::U) = 0;
  //               grid.cell(i-1,j).velocity(velocity_type::U) = -(1.0) *
  //               grid.cell(i-1,j-1).velocity(velocity_type::U);
  //               grid.cell(i,j-1).velocity(velocity_type::V) = 0;
  //               grid.cell(i,j).velocity(velocity_type::V) = -(1.0) *
  //               grid.cell(i+1,j).velocity(velocity_type::V) ;
  //               // U[i][j]=0;
  //               // U[i-1][j] = -U[i-1][j-1];
  //               // V[i][j-1]=0;
  //               // V[i][j] = -V[i+1][j];
  //               }

  //         if ( B_SW(grid.cell(i,j).flag()) ){
  //               grid.cell(i-1,j).velocity(velocity_type::U) = 0;
  //               grid.cell(i,j).velocity(velocity_type::U) = -(1.0) *
  //               grid.cell(i,j-1).velocity(velocity_type::U);
  //               grid.cell(i,j-1).velocity(velocity_type::V) = 0;
  //               grid.cell(i,j).velocity(velocity_type::V) = -(1.0) *
  //               grid.cell(i-1,j).velocity(velocity_type::V) ;
  //             // U[i-1][j] = 0;
  //             // U[i][j] = -U[i][j-1];
  //             // V[i][j-1] = 0;
  //             // V[i][j] = -V[i-1][j];
  //             }
  //         break;

  //       case 1<<3:
  //         grid.cell(i,j).velocity(velocity_type::U) =
  //         grid.cell(i-1,j).velocity(velocity_type::U);;
  //         grid.cell(i,j).velocity(velocity_type::V) =
  //         grid.cell(i-1,j).velocity(velocity_type::V);
  //         grid.cell(i,j-1).velocity(velocity_type::V) =
  //         grid.cell(i-1,j-1).velocity(velocity_type::V) ; break;
  //       case (1<<4):
  //         grid.cell(i,j).velocity(velocity_type::U) = 1.0;
  //         grid.cell(i,j).velocity(velocity_type::V) =  0.0;
  //         grid.cell(i,j-1).velocity(velocity_type::V) = 0.0;
  //         break;
  //     }
  //   }
  // }
}

int B_N(int flag) {
  return ((flag & (1 << 5)) && ~(flag & ((1 << 7) | (1 << 8))));
}

int B_S(int flag) {
  return ((flag & (1 << 6)) && ~(flag & ((1 << 7) | (1 << 8))));
}

int B_W(int flag) {

  return ((flag & (1 << 7)) && ~(flag & ((1 << 5) | (1 << 6))));
}

int B_E(int flag) {
  return ((flag & (1 << 8)) && ~(flag & ((1 << 5) | (1 << 6))));
}

int B_NW(int flag) { return ((flag & (1 << 7)) && (flag & (1 << 5))); }

int B_NE(int flag) { return ((flag & (1 << 8)) && (flag & (1 << 5))); }

int B_SW(int flag) { return ((flag & (1 << 7)) && (flag & (1 << 6))); }

int B_SE(int flag) { return ((flag & (1 << 8)) && (flag & (1 << 6))); }
