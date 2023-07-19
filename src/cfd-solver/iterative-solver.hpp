#include <iostream>
#include <iostream>
#include <math.h>
#include "../../datastructures.hpp"


void solve(matrix<double> A, std::vector<double> b, std::vector<double> &x,const double eps,const int dim);

void gs(matrix<double> A, std::vector<double> b, std::vector<double> &x,
        double *res, int dim); 

void sor(matrix<double> A, std::vector<double> b, std::vector<double> &x, 
         double *res, int dim, double omg);

void calculateResidual(matrix<double>  A, std::vector<double> b, std::vector<double> &x, int dim, double* res);
