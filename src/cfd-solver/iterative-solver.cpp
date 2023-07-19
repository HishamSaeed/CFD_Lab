
#include "iterative-solver.hpp"



void solve(matrix<double> A, std::vector<double> b, std::vector<double> &x,const double eps,const int dim) {
    int itermax = 1;
    int i = 0;
    double res = 1.0;
    double omg = 1.5;

    while (res > eps) {
        sor(A, b, x, &res, dim, omg);
        // gs(A, b, x, &res, dim);
        i++;
    }
    std::cout << "number of iteration = " << i << std::endl;
}

void gs(matrix<double> A, std::vector<double> b, std::vector<double> &x,
        double *res, int dim) {

    double sum = 0.0;
    for(int i = 0; i < dim; i++) {
        sum = 0.0;
 
        // square brackets is faster to access data, but
        // .at() is safer and prohibits accessing elements out of array
        //  limit. But in release using square brackets is better
        for(int j = 0; j < i; j++) {
            sum += A[i][j] * x[j];
        }

        for(int j = i+1; j < dim; j++) {
            sum += A[i][j] * x[j];
        }

        x[i] = ((b[i] - sum)/A[i][i]);
    }
    calculateResidual(A, b, x, dim, res);   
}

void sor(matrix<double> A, std::vector<double> b, std::vector<double> &x, 
         double *res, int dim, double omg) {

    double sum = 0.0;

    for(int i = 0; i < dim; i++) {
        sum = 0.0;
 
        for(int j = 0; j < i; j++) {
            sum += A[i][j] * x[j];
        }

        for(int j = i+1; j < dim; j++) {
            sum += A[i][j] * x[j];
        }

        x[i] = (1 - omg) * x[i] + omg * ((b[i] - sum)/A[i][i]);
    }
    calculateResidual(A, b, x, dim, res);
}


void calculateResidual(matrix<double> A, std::vector<double> b, std::vector<double> &x, int dim, double* res) {

    double sum = 0.0;
    double rloc = 0.0;
    for(int k = 0; k < dim; k++) {
         sum = 0.0;
        for(int m = 0; m < dim; m++) {
            sum += A[k][m] * x[m];
        }
        rloc += std::pow((b[k] - sum),2);
    }
    rloc = rloc/dim;
    *res = std::sqrt(rloc);
}