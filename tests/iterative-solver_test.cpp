#define CATCH_CONFIG_MAIN
#include <math.h>
#include <iostream>
#include <catch/catch.hpp>
#include "../utilities.hpp"
#include "../src/cfd-solver/iterative-solver.hpp"
#include "../datastructures.hpp"
// Testing SOR method


/**
 * Uniform Grid, square domain from 0 to 1
 * N = Nx = Ny, uniform grid, is number of nodes per dim, total number of nodes in the 
 * grid is N * N 
 * n = nx = ny = (N-2) is the number of inner nodes per dim 
 * to be computed without boundaries, total number of inner nodes is n * n
 * s = n * n, matrix size, square matrix 
 * h = hx = hy = 1/(N-1)
 * ig = jg = 0 -> N, ig = jg = 1 -> n for inner nodes, g stands for grid. indicies of
 * the grid are different from indicies of matrix
*/

TEST_CASE("Test sor", "[test_sor]") {
    int N = 5;
    int n = N - 2;
    int num_rows = n * n;
    int num_cols = n * n;
    int num_elements = num_rows * num_cols;
    double h = 1.0/(N-1);
    double h_squared = std::pow(h,2);
    double eps = 0.0001;
    // Allocating matrix
    matrix<double> A(num_rows, std::vector<double>(num_elements, 0));

    double diag = -4 * (1/h_squared);
    double lu_diag = 1 * (1/h_squared);
    // Fill the diagonal
    for(int i = 0; i < num_rows; i++) {
      A[i][i] = diag; 
    }
    
    // Fill lower and upper diagonal below and below main diagonal n diagonals
    for(int i = n; i < num_rows; i++) {
      A[i][i-n] = lu_diag;
      A[i-n][i] = lu_diag;
    }
    
    // Fill lower and upper diagonal right exactly below and above main diagonal
    for(int i = 1; i < num_rows; i++) {
      A[i][i-1] = lu_diag;
      A[i-1][i] = lu_diag;
    }
    
    // Fill lower and upper diagonal with the zeros in the correct place
    // these diagonals are not completly filled with 1. 
    for(int i = n; i < num_rows; i=i+n) {
      A[i][i-1] = 0;
      A[i-1][i] = 0;
    }
    // grid points
    double x_i = 0.0;
    double y_i = 0.0;
    int k = 0;
    std::vector<double> x_analytical(num_rows, 0);
    std::vector<double> b(num_rows, 0);
    /**
     * fill analytical solution vector
     * from 2d to 1d (width * row + col)
     * from 1d to 2d row = index/width & col = index%width
     * make sure to check what is row and what is col, indices might switch
     * like the difference between the grid and matrix
    */
    for(int j = 1; j <= n; j++) {
      for(int i = 1; i <= n; i++) {
        x_i = i * h;
        y_i = j * h;
        k = n * (j-1) + (i-1);
        x_analytical[k] = std::sin(M_PI * x_i) * std::sin(M_PI * y_i);
        b[k] = -2 * (std::pow(M_PI,2)) * std::sin(M_PI * x_i) * std::sin(M_PI * y_i);
      }
    }

    std::vector<double> x(num_rows, 0);

    solve(A, b, x, eps, num_rows);

    double err = 0.0;
    for(int i = 0; i < num_rows; i++) {
      err += std::pow((x[i]-x_analytical[i]), 2);
    }
    err = err * (1.0/num_rows);
    err = std::sqrt(err);

    std::cout << "the error norm = " << err << std::endl;
}
