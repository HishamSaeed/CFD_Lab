#include "uvp.hpp"
#include "boundary_val.hpp"
#include "datastructures.hpp"
#include "helper.hpp"
#include <cassert>
#include <cfloat>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
// ----------------------------------------------------------------------------------------------------
// Determines the value of F and G
void calculate_fg(double Re, double GX, double GY, double alpha, double dt,
                  double dx, double dy, int imax, int jmax, Grid &grid,
                  matrix<double> &F, matrix<double> &G) {
  // For the sake of readability: declaration/definition of variables needed in
  // calculation of F and G;
  double ui_j, uip1_j, uim1_j, ui_jp1, ui_jm1,
      uim1_jp1; // uip1_j = u[i+1][j], uim1_j = = u[i-1][j]
  double vi_j, vip1_j, vim1_j, vi_jp1, vi_jm1, vip1_jm1; // vi_jp1 = v[i][j+1]
  double Dx = 1 / dx;
  double Dy = 1 / dy;
  double kvis = 1 / Re; // kinematic viscosity

  // Extract Velocity matrices from Grid object
  // Thes matrices can be used to access the velocities instead of calling the
  // cell get function
  matrix<double> umatrix;
  matrix<double> vmatrix;

  // matrix<double> U;
  // matrix<double> V;
  grid.velocity(umatrix, velocity_type::U);
  grid.velocity(vmatrix, velocity_type::V);
  // grid.velocity(U, velocity_type::U);
  // grid.velocity(V, velocity_type::V);

  //        double a, b,c,d,du2x2,du2y2,du2dx,duvy,dv2y2,dv2x2,dv2dy,duvx;
  // Calculation of F
  // --------------------------------------------------------------------------------------
  for (int i = 1; i <= imax; i++) {
    for (int j = 1; j <= jmax; j++) {

      if (((grid.cell(i, j).flag() & (1 << 0)) & grid.cell(i + 1, j).flag()) ||
          ((grid.cell(i, j).flag() & (1 << 0)) &&
           (grid.cell(i + 1, j).flag() & (1 << 3)))) {
        // Reading the velocity in the x-direction of the current cell and the
        // sourrounding ones.
        ui_j = umatrix[i][j];
        uip1_j = umatrix[i + 1][j];
        uim1_j = umatrix[i - 1][j];
        ui_jp1 = umatrix[i][j + 1];
        ui_jm1 = umatrix[i][j - 1];

        uim1_jp1 = umatrix[i - 1][j + 1];

        // Reading the velocity in the x-direction of the current cell and the
        // sourrounding ones.
        vi_j = vmatrix[i][j];
        vip1_j = vmatrix[i + 1][j];
        vim1_j = vmatrix[i - 1][j];
        vi_jp1 = vmatrix[i][j + 1];
        vi_jm1 = vmatrix[i][j - 1];

        vip1_jm1 = vmatrix[i + 1][j - 1];

        // Calculation of F[i][j]:
        F[i][j] =
            ui_j +
            dt * (
                     // FD Aprroximation of U second order derivative in x and y
                     // directions
                     (kvis * ((Dx * Dx * (uip1_j - (2 * ui_j) + uim1_j)) +
                              (Dy * Dy * (ui_jp1 - (2 * ui_j) + ui_jm1)))) -
                     // d(u^2)/dx term in Doner scheme:
                     ((Dx * ((pow(0.5 * (ui_j + uip1_j), 2)) -
                             (pow(0.5 * (uim1_j + ui_j), 2)))) +
                      (alpha * Dx *
                       (((0.5 * abs(ui_j + uip1_j)) * (0.5 * (ui_j - uip1_j))) -
                        ((0.5 * abs(uim1_j + ui_j)) *
                         (0.5 * (uim1_j - ui_j)))))) -
                     // d(uv)/dy term in Doner scheme:
                     ((Dy *
                       (((0.5 * (vi_j + vip1_j)) * (0.5 * (ui_j + ui_jp1))) -
                        ((0.5 * (vi_jm1 + vip1_jm1)) *
                         (0.5 * (ui_jm1 + ui_j))))) +
                      (alpha * Dy *
                       (((0.5 * abs(vi_j + vip1_j)) * (0.5 * (ui_j - ui_jp1))) -
                        ((0.5 * abs(vi_jm1 + vip1_jm1)) *
                         (0.5 * (ui_jm1 - ui_j)))))) +
                     // Volume forces
                     GX);

        // du2x2= (U[i+1][j]-2*U[i][j]+U[i-1][j])/(dx*dx);
        // du2y2= (U[i][j+1]-2*U[i][j]+U[i][j-1])/(dy*dy);
        // a=(U[i][j]+U[i+1][j])/2;
        // b=(U[i-1][j]+U[i][j])/2;
        // du2dx=(a*a-b*b+
        // alpha*(fabs(a)*((U[i][j]-U[i+1][j])/2)-fabs(b)*((U[i-1][j]-U[i][j])/2)))/dx;
        // duvy=((V[i][j]+V[i+1][j])*(U[i][j]+U[i][j+1])-(V[i][j-1]+V[i+1][j-1])*(U[i][j-1]+U[i][j])+alpha*(fabs(V[i][j]+V[i+1][j])*(U[i][j]-U[i][j+1])-fabs(V[i][j-1]+V[i+1][j-1])*(U[i][j-1]-U[i][j])))/(4*dy);
        // F[i][j]=U[i][j]+dt*((du2x2+du2y2)*(1/Re)-du2dx-duvy+GX);
      }
    }
  }

  // Calculation of G
  // --------------------------------------------------------------------------------------
  for (int i = 1; i <= imax; i++) {
    for (int j = 1; j <= jmax; j++) {
      if ((grid.cell(i, j).flag() & (1 << 0)) & grid.cell(i, j + 1).flag()) {
        // Reading the velocity in the x-direction of the current cell and the
        // sourrounding ones.

        ui_j = umatrix[i][j];
        uip1_j = umatrix[i + 1][j];
        uim1_j = umatrix[i - 1][j];
        ui_jp1 = umatrix[i][j + 1];
        ui_jm1 = umatrix[i][j - 1];

        uim1_jp1 = umatrix[i - 1][j + 1];

        // Reading the velocity in the x-direction of the current cell and the
        // sourrounding ones.

        vi_j = vmatrix[i][j];
        vip1_j = vmatrix[i + 1][j];
        vim1_j = vmatrix[i - 1][j];
        vi_jp1 = vmatrix[i][j + 1];
        vi_jm1 = vmatrix[i][j - 1];

        vip1_jm1 = vmatrix[i + 1][j - 1];

        // Calculation of G[i][j]:
        G[i][j] =
            vi_j +
            dt * (
                     // FD Aprroximation of V second order derivative in x and y
                     // directions
                     (kvis * ((Dx * Dx * (vip1_j - (2 * vi_j) + vim1_j)) +
                              (Dy * Dy * (vi_jp1 - (2 * vi_j) + vi_jm1)))) -
                     // d(uv)/dx term in Doner scheme:
                     ((Dx *
                       (((0.5 * (ui_j + ui_jp1)) * (0.5 * (vi_j + vip1_j))) -
                        ((0.5 * (uim1_j + uim1_jp1)) *
                         (0.5 * (vim1_j + vi_j))))) +
                      (alpha * Dx *
                       (((0.5 * abs(ui_j + ui_jp1)) * (0.5 * (vi_j - vip1_j))) -
                        ((0.5 * abs(uim1_j + uim1_jp1)) *
                         (0.5 * (vim1_j - vi_j)))))) -
                     // d(v^2)/dy term in Doner scheme:
                     ((Dy * ((pow(0.5 * (vi_j + vi_jp1), 2)) -
                             (pow(0.5 * (vi_jm1 + vi_j), 2)))) +
                      (alpha * Dy *
                       (((0.5 * abs(vi_j + vi_jp1)) * (0.5 * (vi_j - vi_jp1))) -
                        ((0.5 * abs(vi_jm1 + vi_j)) *
                         (0.5 * (vi_jm1 - vi_j)))))) +
                     // volume forces
                     GY);

        // dv2y2= (V[i][j+1]-2*V[i][j]+V[i][j-1])/(dy*dy);
        // dv2x2= (V[i+1][j]-2*V[i][j]+V[i-1][j])/(dx*dx);
        // c=(V[i][j]+V[i][j+1])/2;
        // d=(V[i][j-1]+V[i][j])/2;
        // dv2dy=(c*c-d*d+
        // alpha*(fabs(c)*((V[i][j]-V[i][j+1])/2)-fabs(d)*((V[i][j-1]-V[i][j])/2)))/dy;
        // duvx=((U[i][j]+U[i][j+1])*(V[i][j]+V[i+1][j])-(U[i-1][j]+U[i-1][j+1])*(V[i-1][j]+V[i][j])+alpha*(fabs(U[i][j]+U[i][j+1])*(V[i][j]-V[i+1][j])-fabs(U[i-1][j]+U[i-1][j+1])*(V[i-1][j]-V[i][j])))/(4*dy);
        // G[i][j]=V[i][j]+dt*((dv2x2+dv2y2)*(1/Re)-dv2dy-duvx+GY);
      }
    }
  }

  // Updating the Boundary Conditions for F;
  for (int j = 1; j <= jmax; j++) {

    // F[0][j]    = umatrix[0][j];
    // F[imax][j] = umatrix[imax][j];
    switch (grid.cell(0, j).flag() &
            ((1 << 1) | (1 << 2) | (1 << 3) | (1 << 4))) {
    case 1 << 1: // No slip conditions
      if (B_E(grid.cell(0, j).flag())) {
        F[0][j] = umatrix[0][j];
        //        F[0][j]    = U[0][j];
      }

      break;
    case (1 << 4):
      F[0][j] = umatrix[0][j];
      break;
    }

    switch (grid.cell(imax + 1, j).flag() &
            ((1 << 1) | (1 << 2) | (1 << 3) | (1 << 4))) {
    case 1 << 1: // No slip conditions
      if (B_W(grid.cell(imax + 1, j).flag()))

      {
        F[imax][j] = umatrix[imax][j];
      }
      break;
    case 1 << 3:
      // F[imax+1][j] = umatrix[imax+1][j];
      // F[imax+1][j] = U[imax+1][j];
      // F[imax][j] = U[imax][j];
      break;
    }
  }

  for (int i = 1; i <= imax; i++) {
    // G[i][0]    = vmatrix[i][0];
    // G[i][jmax] = vmatrix[i][jmax];
    switch (grid.cell(i, 0).flag() &
            ((1 << 1) | (1 << 2) | (1 << 3) | (1 << 4))) {
    case 1 << 1: // No slip conditions
      if (B_N(grid.cell(i, 0).flag())) {
        G[i][0] = vmatrix[i][0];
        // G[i][0]    = V[i][0];
      }

      break;
    }

    switch (grid.cell(i, jmax + 1).flag() &
            ((1 << 1) | (1 << 2) | (1 << 3) | (1 << 4))) {
    case 1 << 1: // No slip conditions
      if (B_S(grid.cell(i, jmax + 1).flag())) {
        G[i][jmax] = vmatrix[i][jmax];
        //       G[i][jmax] = V[i][jmax];
      }
      break;
    }
  }

  for (int i = 1; i <= imax; ++i) {
    for (int j = 1; j <= jmax; ++j) {

      if (B_NE(grid.cell(i, j).flag())) {
        F[i][j] = umatrix[i][j];
        G[i][j] = vmatrix[i][j];
      }

      if (B_NW(grid.cell(i, j).flag())) {
        F[i - 1][j] = umatrix[i - 1][j];
        G[i][j] = vmatrix[i][j];
      }

      if (B_SE(grid.cell(i, j).flag())) {
        F[i][j] = umatrix[i][j];
        G[i][j - 1] = vmatrix[i][j - 1];
      }

      if (B_SW(grid.cell(i, j).flag())) {
        F[i - 1][j] = umatrix[i - 1][j];
        G[i][j - 1] = vmatrix[i][j - 1];
      }

      // if ( B_NE(grid.cell(i,j).flag()) ) { F[i][j] = U[i][j]; G[i][j] =
      // V[i][j]; }

      // if ( B_NW(grid.cell(i,j).flag()) ) { F[i-1][j] = U[i-1][j]; G[i][j] =
      // V[i][j]; }

      // if ( B_SE(grid.cell(i,j).flag()) ) { F[i][j] = U[i][j]; G[i][j-1] =
      // V[i][j-1]; }

      // if ( B_SW(grid.cell(i,j).flag()) ) { F[i-1][j] = U[i-1][j]; G[i][j-1] =
      // V[i][j-1]; }
    }
  }
}
// ----------------------------------------------------------------------------------------------------

// ----------------------------------------------------------------------------------------------------
// This operation computes the right hand side of the pressure poisson equation.
void calculate_rs(double dt, double dx, double dy, int imax, int jmax,
                  matrix<double> &F, matrix<double> &G, matrix<double> &RS,
                  Grid &grid) {

  for (int i = 1; i <= imax; i++) {
    for (int j = 1; j <= jmax; j++) {
      if (grid.cell(i, j).flag() & (1 << 0)) {
        RS[i][j] = (1 / dt) * (((F[i][j] - F[i - 1][j]) / dx) +
                               ((G[i][j] - G[i][j - 1]) / dy));
      }
    }
  }
}

// ----------------------------------------------------------------------------------------------------
// Determines the maximal time step size
void calculate_dt(double Re, double tau, double *dt, double dx, double dy,
                  int imax, int jmax, double PR, bool *temp_flag, Grid &grid) {
  // Extract Velocity matrices from Grid object
  // Thes matrices can be used to access the velocities instead of calling the
  // cell get function
  matrix<double> umatrix;
  matrix<double> vmatrix;

  grid.velocity(umatrix, velocity_type::U);
  grid.velocity(vmatrix, velocity_type::V);

  double umax_abs = 0, vmax_abs = 0;
  double Dx = 1 / dx;
  double Dy = 1 / dy;
  // std::cout << "umax init " << umax_abs << std::endl;
  // Search for the maximum absolute value
  for (int i = 1; i <= imax; i++) {
    for (int j = 1; j <= jmax; j++) {
      if (grid.cell(i, j).flag() & (1 << 0))
        // {

        //   std::cout << "U[" << i << "][" << j << "] = " << umatrix[i][j] <<
        //   std::endl;
        if (fabs(umatrix[i][j]) > umax_abs) {
          umax_abs = fabs(umatrix[i][j]);
          //   std::cout << "i = " << i << " j = " << j << std::endl;
          //   std::cout << "U = " << std::setprecision(15) << umatrix[i][j] <<
          //   std::endl; std::cout << "From grid " << std::setprecision(15) <<
          //   grid.cell(i,j).velocity(velocity_type::U) << std::endl; std::cout
          //   << "Umax = " << umax_abs << std::endl;
        }
      if (fabs(vmatrix[i][j]) > vmax_abs) {
        vmax_abs = fabs(vmatrix[i][j]);
      }
    }
  }
  // std::cout << "Umax = " << umax_abs << std::endl;
  // std::cout << "Vmax = " << vmax_abs << std::endl;
  // Checking for the minimum value based on stability conditions
  *dt = std::min((Re / 2) * pow((Dx * Dx) + (Dy * Dy), -1),
                 std::min(dx / umax_abs, dy / vmax_abs));
  // std::cout << "dt = " << *dt << std::endl;
  if (*temp_flag == true) {
    // *dt = std::min(*dt, ((1/alpha)*(1/(pow(Dx, 2) + pow(Dy, 2)))));
    // std::cout << "temp stability cond " << (0.5 * (PR) * (Re)*  pow( (
    // (Dx*Dx) + (Dy*Dy) ),-1)) << std::endl;
    *dt = std::min(*dt, (0.5 * (PR) * (Re)*pow(((Dx * Dx) + (Dy * Dy)), -1)));

    // std::cout << "dt after temp cond = " << *dt << std::endl;
  }
  *dt = tau * (*dt);

  // Checking for the minimum value based on stability conditions
  // *dt = tau * std::min( (Re/2)*pow((Dx*Dx)+(Dy*Dy),-1) , std::min(
  // dx/umax_abs , dy/vmax_abs ));
}
// ----------------------------------------------------------------------------------------------------
void calculate_uv(double dt, double dx, double dy, int imax, int jmax,
                  Grid &grid, matrix<double> &F, matrix<double> &G)

{

  matrix<double> pmatrix;
  grid.pressure(pmatrix);

  // Calculate u
  for (int i = 1; i <= imax; i++) {
    for (int j = 1; j <= jmax; j++) {
      if (((grid.cell(i, j).flag() & (1 << 0)) & grid.cell(i + 1, j).flag()) ||
          ((grid.cell(i, j).flag() & (1 << 0)) &&
           (grid.cell(i + 1, j).flag() & (1 << 3)))) {
        grid.cell(i, j).velocity(velocity_type::U) =
            F[i][j] - ((dt / dx) * (pmatrix[i + 1][j] - pmatrix[i][j]));
      }
    }
  }

  // Calculate v
  for (int i = 1; i <= imax; i++) {
    for (int j = 1; j <= jmax; j++) {
      if ((grid.cell(i, j).flag() & (1 << 0)) & grid.cell(i, j + 1).flag()) {
        grid.cell(i, j).velocity(velocity_type::V) =
            G[i][j] - ((dt / dy) * (pmatrix[i][j + 1] - pmatrix[i][j]));
      }
    }
  }
}
// ----------------------------------------------------------------------------------------------------
void calculate_temp(double dt, double dx, double dy, int imax, int jmax,
                    double Re, double PR, double alpha, double heat_source,
                    Grid &grid) {
  matrix<double> U;
  matrix<double> V;
  matrix<double> T;

  grid.velocity(U, velocity_type::U);
  grid.velocity(V, velocity_type::V);
  grid.temperature(T);

  double Dx = 1 / dx;
  double Dy = 1 / dy;
  double convection = 0.0;
  double diffusion = 0.0;

  for (int i = 0; i <= imax + 1; ++i) {
    for (int j = 0; j <= jmax + 1; ++j) {
      if (B_E(grid.cell(i, j).flag())) {
        T[i][j] = T[i + 1][j];
        grid.cell(i, j).temperature() = T[i + 1][j];
      }

      if (B_W(grid.cell(i, j).flag())) {
        T[i][j] = T[i - 1][j];
        grid.cell(i, j).temperature() = T[i - 1][j];
      }

      if (B_N(grid.cell(i, j).flag())) {
        T[i][j] = T[i][j + 1];
        grid.cell(i, j).temperature() = T[i][j + 1];
      }

      if (B_S(grid.cell(i, j).flag())) {
        T[i][j] = T[i][j - 1];
        grid.cell(i, j).temperature() = T[i][j - 1];
      }

      if (B_NE(grid.cell(i, j).flag())) {
        T[i][j] = (T[i][j + 1] + T[i + 1][j]) / 2;
        grid.cell(i, j).temperature() = (T[i][j + 1] + T[i + 1][j]) / 2;
      }

      if (B_NW(grid.cell(i, j).flag())) {
        T[i][j] = (T[i][j + 1] + T[i - 1][j]) / 2;
        grid.cell(i, j).temperature() = (T[i][j + 1] + T[i - 1][j]) / 2;
      }

      if (B_SE(grid.cell(i, j).flag())) {
        T[i][j] = (T[i][j - 1] + T[i + 1][j]) / 2;
        grid.cell(i, j).temperature() = (T[i][j - 1] + T[i + 1][j]) / 2;
      }

      if (B_SW(grid.cell(i, j).flag())) {
        T[i][j] = (T[i][j - 1] + T[i - 1][j]) / 2;
        grid.cell(i, j).temperature() = (T[i][j - 1] + T[i - 1][j]) / 2;
      }

      if (grid.cell(i, j).flag() & (1 << 3)) {
        T[i][j] = T[i - 1][j];
        grid.cell(i, j).temperature() = T[i - 1][j];
      }

      if (grid.cell(i, j).flag() & (1 << 4)) {
        T[i][j] = 0.0;
        grid.cell(i, j).temperature() = 0.0;
      }
    }
  }

  for (int j = 0; j <= jmax + 1; j++) {
    T[0][j] = 2 * 1.0 - T[1][j];
    grid.cell(0, j).temperature() = 2 * 1.0 - T[1][j];
    T[imax + 1][j] = (2 * (0.0)) - T[imax][j];
    grid.cell(imax + 1, j).temperature() = (2 * (0.0)) - T[imax][j];
  }

  // for(int i=0; i<=imax+1; i++)
  // {
  // 	T[i][0] = 2*294.78 - T[i][1];
  //         grid.cell(i,0).temperature() = 2*294.78 - T[i][0];
  // 	T[i][jmax+1] = (2*(291.20))- T[i][jmax];
  //         grid.cell(i,jmax+1).temperature() = (2*(291.20)) - T[i][jmax];

  // }

  double dut_dx;
  double dvt_dy;
  double dt2_dx2;
  double dt2_dy2;
  double Z;
  for (int i = 1; i <= imax; i++) {
    for (int j = 1; j <= jmax; j++) {
      if (grid.cell(i, j).flag() & (1 << 0)) {

        dut_dx = (1 / dx) * ((U[i][j] * (T[i][j] + T[i + 1][j]) * 0.5) -
                             (U[i - 1][j] * (T[i - 1][j] + T[i][j]) * 0.5)) +
                 (alpha / dx) *
                     ((fabs(U[i][j]) * (T[i][j] - T[i + 1][j]) * 0.5) -
                      (fabs(U[i - 1][j]) * (T[i - 1][j] - T[i][j]) * 0.5));

        dvt_dy = (1 / dy) * ((V[i][j] * (T[i][j] + T[i][j + 1]) * 0.5) -
                             (V[i][j - 1] * (T[i][j - 1] + T[i][j]) * 0.5)) +
                 (alpha / dy) *
                     ((fabs(V[i][j]) * (T[i][j] - T[i][j + 1]) * 0.5) -
                      (fabs(V[i][j - 1]) * (T[i][j - 1] - T[i][j]) * 0.5));

        dt2_dx2 = (T[i + 1][j] - 2 * T[i][j] + T[i - 1][j]) / (dx * dx);

        dt2_dy2 = (T[i][j + 1] - 2 * T[i][j] + T[i][j - 1]) / (dy * dy);

        Z = (1 / (Re * PR)) * (dt2_dx2 + dt2_dy2) - dut_dx - dvt_dy;

        // convection = (Dx * ((U[i][j] * ( (T[i][j] + T[i+1][j])/2) ) -
        // (U[i-1][j] * ((T[i-1][j] + T[i][j])/2) )) +
        //              gamma * Dx * ((abs(U[i][j]) * ( (T[i][j] - T[i+1][j])/2)
        //              ) - (abs(U[i-1][j]) * ( (T[i-1][j] - T[i][j])/2) ) ) ) +
        //         //      (Dy * ((V[i][j] * (T[i][j] + T[i][j+1])/2) -
        //         (V[i-1][j] * (T[i][j-1] + T[i][j])/2)) +
        //              (Dy * ((V[i][j] * ((T[i][j] + T[i][j+1])/2) ) -
        //              (V[i][j-1] * ((T[i][j-1] + T[i][j])/2) )) +
        //         //      gamma * Dy * ((abs(V[i][j]) * (T[i][j] -
        //         T[i][j+1])/2) - (abs(V[i-1][j]) * (T[i][j-1] - T[i][j])/2))
        //         );
        //              gamma * Dy * ((abs(V[i][j]) * ((T[i][j] - T[i][j+1])/2)
        //              ) - (abs(V[i][j-1]) * ((T[i][j-1] - T[i][j])/2) )) );
        // // diffusion = alpha * (pow(Dx, 2) * (T[i+1][j] -2*T[i][j] +
        // T[i-1][j]) +
        // //                      pow(Dy, 2) * (T[i][j+1] -2*T[i][j] +
        // T[i][j-1]));
        // // diffusion = (1/PR) * (1/Re) * ( (pow(Dx, 2) * (T[i+1][j]
        // -2*T[i][j] + T[i-1][j]) ) +
        // //                      (pow(Dy, 2) * (T[i][j+1] -2*T[i][j] +
        // T[i][j-1]))); diffusion = (1/(Re*PR)) * ( (pow(Dx, 2) * (T[i+1][j]
        // -2*T[i][j] + T[i-1][j]) ) +
        //                      (pow(Dy, 2) * (T[i][j+1] -2*T[i][j] +
        //                      T[i][j-1])));

        // std::cout << PR << std::endl;
        // T[i][j] = T[i][j] + dt * (-convection + diffusion) + heat_source;
        // grid.cell(i,j).temperature() = T[i][j] + dt * (-convection +
        // diffusion) + heat_source;
        grid.cell(i, j).temperature() = T[i][j] + (dt * Z);
      }
    }
  }

  // grid.set_temperature(T);
}

void calculate_fg_temp(double Re, double GX, double GY, double alpha, double dt,
                       double dx, double dy, int imax, int jmax, double beta,
                       Grid &grid, matrix<double> &F, matrix<double> &G) {
  // For the sake of readability: declaration/definition of variables needed in
  // calculation of F and G;
  double ui_j, uip1_j, uim1_j, ui_jp1, ui_jm1,
      uim1_jp1; // uip1_j = u[i+1][j], uim1_j = = u[i-1][j]
  double vi_j, vip1_j, vim1_j, vi_jp1, vi_jm1, vip1_jm1; // vi_jp1 = v[i][j+1]
  double Dx = 1 / dx;
  double Dy = 1 / dy;
  double kvis = 1 / Re; // kinematic viscosity

  // Extract Velocity matrices from Grid object
  // Thes matrices can be used to access the velocities instead of calling the
  // cell get function
  matrix<double> umatrix;
  matrix<double> vmatrix;
  matrix<double> T;

  // matrix<double> U;
  // matrix<double> V;
  grid.velocity(umatrix, velocity_type::U);
  grid.velocity(vmatrix, velocity_type::V);
  grid.temperature(T);
  // grid.velocity(U, velocity_type::U);
  // grid.velocity(V, velocity_type::V);

  //        double a, b,c,d,du2x2,du2y2,du2dx,duvy,dv2y2,dv2x2,dv2dy,duvx;
  // Calculation of F
  // --------------------------------------------------------------------------------------
  for (int i = 1; i <= imax; i++) {
    for (int j = 1; j <= jmax; j++) {

      if (((grid.cell(i, j).flag() & (1 << 0)) & grid.cell(i + 1, j).flag()) ||
          ((grid.cell(i, j).flag() & (1 << 0)) &&
           (grid.cell(i + 1, j).flag() & (1 << 3)))) {
        // Reading the velocity in the x-direction of the current cell and the
        // sourrounding ones.
        ui_j = umatrix[i][j];
        uip1_j = umatrix[i + 1][j];
        uim1_j = umatrix[i - 1][j];
        ui_jp1 = umatrix[i][j + 1];
        ui_jm1 = umatrix[i][j - 1];

        uim1_jp1 = umatrix[i - 1][j + 1];

        // Reading the velocity in the x-direction of the current cell and the
        // sourrounding ones.
        vi_j = vmatrix[i][j];
        vip1_j = vmatrix[i + 1][j];
        vim1_j = vmatrix[i - 1][j];
        vi_jp1 = vmatrix[i][j + 1];
        vi_jm1 = vmatrix[i][j - 1];

        vip1_jm1 = vmatrix[i + 1][j - 1];

        // Calculation of F[i][j]:
        F[i][j] =
            ui_j +
            dt * (
                     // FD Aprroximation of U second order derivative in x and y
                     // directions
                     (kvis * ((Dx * Dx * (uip1_j - (2 * ui_j) + uim1_j)) +
                              (Dy * Dy * (ui_jp1 - (2 * ui_j) + ui_jm1)))) -
                     // d(u^2)/dx term in Doner scheme:
                     ((Dx * ((pow(0.5 * (ui_j + uip1_j), 2)) -
                             (pow(0.5 * (uim1_j + ui_j), 2)))) +
                      (alpha * Dx *
                       (((0.5 * abs(ui_j + uip1_j)) * (0.5 * (ui_j - uip1_j))) -
                        ((0.5 * abs(uim1_j + ui_j)) *
                         (0.5 * (uim1_j - ui_j)))))) -
                     // d(uv)/dy term in Doner scheme:
                     ((Dy *
                       (((0.5 * (vi_j + vip1_j)) * (0.5 * (ui_j + ui_jp1))) -
                        ((0.5 * (vi_jm1 + vip1_jm1)) *
                         (0.5 * (ui_jm1 + ui_j))))) +
                      (alpha * Dy *
                       (((0.5 * abs(vi_j + vip1_j)) * (0.5 * (ui_j - ui_jp1))) -
                        ((0.5 * abs(vi_jm1 + vip1_jm1)) *
                         (0.5 * (ui_jm1 - ui_j)))))) +
                     // Volume forces
                     (GX * beta * (dt / 2) * (T[i][j] + T[i + 1][j]))

                 );
        // F[i][j] = F[i][j] - beta*(dt/2)*(T[i][j] + T[i+1][j])*GX;

        // du2x2= (U[i+1][j]-2*U[i][j]+U[i-1][j])/(dx*dx);
        // du2y2= (U[i][j+1]-2*U[i][j]+U[i][j-1])/(dy*dy);
        // a=(U[i][j]+U[i+1][j])/2;
        // b=(U[i-1][j]+U[i][j])/2;
        // du2dx=(a*a-b*b+
        // alpha*(fabs(a)*((U[i][j]-U[i+1][j])/2)-fabs(b)*((U[i-1][j]-U[i][j])/2)))/dx;
        // duvy=((V[i][j]+V[i+1][j])*(U[i][j]+U[i][j+1])-(V[i][j-1]+V[i+1][j-1])*(U[i][j-1]+U[i][j])+alpha*(fabs(V[i][j]+V[i+1][j])*(U[i][j]-U[i][j+1])-fabs(V[i][j-1]+V[i+1][j-1])*(U[i][j-1]-U[i][j])))/(4*dy);
        // F[i][j]=U[i][j]+dt*((du2x2+du2y2)*(1/Re)-du2dx-duvy+GX*((beta*dt)*(T[i][j]+T[i+1][j]))/2);
      }
    }
  }

  // Calculation of G
  // --------------------------------------------------------------------------------------
  for (int i = 1; i <= imax; i++) {
    for (int j = 1; j <= jmax; j++) {
      if ((grid.cell(i, j).flag() & (1 << 0)) & grid.cell(i, j + 1).flag()) {
        // Reading the velocity in the x-direction of the current cell and the
        // sourrounding ones.

        ui_j = umatrix[i][j];
        uip1_j = umatrix[i + 1][j];
        uim1_j = umatrix[i - 1][j];
        ui_jp1 = umatrix[i][j + 1];
        ui_jm1 = umatrix[i][j - 1];

        uim1_jp1 = umatrix[i - 1][j + 1];

        // Reading the velocity in the x-direction of the current cell and the
        // sourrounding ones.

        vi_j = vmatrix[i][j];
        vip1_j = vmatrix[i + 1][j];
        vim1_j = vmatrix[i - 1][j];
        vi_jp1 = vmatrix[i][j + 1];
        vi_jm1 = vmatrix[i][j - 1];

        vip1_jm1 = vmatrix[i + 1][j - 1];

        // Calculation of G[i][j]:
        G[i][j] =
            vi_j +
            dt * (
                     // FD Aprroximation of V second order derivative in x and y
                     // directions
                     (kvis * ((Dx * Dx * (vip1_j - (2 * vi_j) + vim1_j)) +
                              (Dy * Dy * (vi_jp1 - (2 * vi_j) + vi_jm1)))) -
                     // d(uv)/dx term in Doner scheme:
                     ((Dx *
                       (((0.5 * (ui_j + ui_jp1)) * (0.5 * (vi_j + vip1_j))) -
                        ((0.5 * (uim1_j + uim1_jp1)) *
                         (0.5 * (vim1_j + vi_j))))) +
                      (alpha * Dx *
                       (((0.5 * abs(ui_j + ui_jp1)) * (0.5 * (vi_j - vip1_j))) -
                        ((0.5 * abs(uim1_j + uim1_jp1)) *
                         (0.5 * (vim1_j - vi_j)))))) -
                     // d(v^2)/dy term in Doner scheme:
                     ((Dy * ((pow(0.5 * (vi_j + vi_jp1), 2)) -
                             (pow(0.5 * (vi_jm1 + vi_j), 2)))) +
                      (alpha * Dy *
                       (((0.5 * abs(vi_j + vi_jp1)) * (0.5 * (vi_j - vi_jp1))) -
                        ((0.5 * abs(vi_jm1 + vi_j)) *
                         (0.5 * (vi_jm1 - vi_j)))))) +
                     // volume forces
                     (GY * (beta * dt) * ((T[i][j] + T[i][j + 1]) / 2))

                 );
        // dv2y2= (V[i][j+1]-2*V[i][j]+V[i][j-1])/(dy*dy);
        // dv2x2= (V[i+1][j]-2*V[i][j]+V[i-1][j])/(dx*dx);
        // c=(V[i][j]+V[i][j+1])/2;
        // d=(V[i][j-1]+V[i][j])/2;
        // dv2dy=(c*c-d*d+
        // alpha*(fabs(c)*((V[i][j]-V[i][j+1])/2)-fabs(d)*((V[i][j-1]-V[i][j])/2)))/dy;
        // duvx=((U[i][j]+U[i][j+1])*(V[i][j]+V[i+1][j])-(U[i-1][j]+U[i-1][j+1])*(V[i-1][j]+V[i][j])+alpha*(fabs(U[i][j]+U[i][j+1])*(V[i][j]-V[i+1][j])-fabs(U[i-1][j]+U[i-1][j+1])*(V[i-1][j]-V[i][j])))/(4*dy);
        // G[i][j]=V[i][j]+dt*((dv2x2+dv2y2)*(1/Re)-dv2dy-duvx+GY*(((beta*dt)*(T[i][j]+T[i][j+1]))/2));
      }
    }
  }

  // Updating the Boundary Conditions for F;
  for (int j = 1; j <= jmax; j++) {

    // F[0][j]    = umatrix[0][j];
    // F[imax][j] = umatrix[imax][j];
    switch (grid.cell(0, j).flag() &
            ((1 << 1) | (1 << 2) | (1 << 3) | (1 << 4))) {
    case 1 << 1: // No slip conditions
      if (B_E(grid.cell(0, j).flag())) {
        F[0][j] = umatrix[0][j];
        //        F[0][j]    = U[0][j];
      }

      break;
    case (1 << 4):
      F[0][j] = umatrix[0][j];
      // F[0][j]    = U[0][j];
      break;
    }

    switch (grid.cell(imax + 1, j).flag() &
            ((1 << 1) | (1 << 2) | (1 << 3) | (1 << 4))) {
    case 1 << 1: // No slip conditions
      if (B_W(grid.cell(imax + 1, j).flag())) {
        F[imax][j] = umatrix[imax][j];
        // F[imax][j] = U[imax][j];
      }
      break;
    case 1 << 3:
      // F[imax+1][j] = umatrix[imax+1][j];
      // F[imax+1][j] = U[imax+1][j];
      // F[imax][j] = U[imax][j];
      break;
    }
  }

  for (int i = 1; i <= imax; i++) {
    // G[i][0]    = vmatrix[i][0];
    // G[i][jmax] = vmatrix[i][jmax];
    switch (grid.cell(i, 0).flag() &
            ((1 << 1) | (1 << 2) | (1 << 3) | (1 << 4))) {
    case 1 << 1: // No slip conditions
      if (B_N(grid.cell(i, 0).flag())) {
        G[i][0] = vmatrix[i][0];
        // G[i][0]    = V[i][0];
      }

      break;
    }

    switch (grid.cell(i, jmax + 1).flag() &
            ((1 << 1) | (1 << 2) | (1 << 3) | (1 << 4))) {
    case 1 << 1: // No slip conditions
      if (B_S(grid.cell(i, jmax + 1).flag())) {
        G[i][jmax] = vmatrix[i][jmax];
        //       G[i][jmax] = V[i][jmax];
      }
      break;
    }
  }

  for (int i = 1; i <= imax; ++i) {
    for (int j = 1; j <= jmax; ++j) {

      if (B_NE(grid.cell(i, j).flag())) {
        F[i][j] = umatrix[i][j];
        G[i][j] = vmatrix[i][j];
      }

      if (B_NW(grid.cell(i, j).flag())) {
        F[i - 1][j] = umatrix[i - 1][j];
        G[i][j] = vmatrix[i][j];
      }

      if (B_SE(grid.cell(i, j).flag())) {
        F[i][j] = umatrix[i][j];
        G[i][j - 1] = vmatrix[i][j - 1];
      }

      if (B_SW(grid.cell(i, j).flag())) {
        F[i - 1][j] = umatrix[i - 1][j];
        G[i][j - 1] = vmatrix[i][j - 1];
      }

      // if ( B_NE(grid.cell(i,j).flag()) ) { F[i][j] = U[i][j]; G[i][j] =
      // V[i][j]; }

      // if ( B_NW(grid.cell(i,j).flag()) ) { F[i-1][j] = U[i-1][j]; G[i][j] =
      // V[i][j]; }

      // if ( B_SE(grid.cell(i,j).flag()) ) { F[i][j] = U[i][j]; G[i][j-1] =
      // V[i][j-1]; }

      // if ( B_SW(grid.cell(i,j).flag()) ) { F[i-1][j] = U[i-1][j]; G[i][j-1] =
      // V[i][j-1]; }
    }
  }
}

// Test pre commit changed mode