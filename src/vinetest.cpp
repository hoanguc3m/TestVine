#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::interfaces(r,cpp)]]

void (*Hfunc2) (int* family,int* n,double* v,double* u,double* theta,double* nu,double* out);

extern "C" void R_init_TestVine(DllInfo *dll) {
    Hfunc2 = (void (*) (int* ,int* ,double* ,double* ,double* ,double* ,double* )) R_GetCCallable("VineCopula", "Hfunc2");
}


// [[Rcpp::export]]
double Hfunc2_call() {

    double rho_value = 0.5;
    double u_val = 0.1;
    double v_val = 0.5;
    double u_double;

    int family_id = 1;
    int numb_obs = 1;
    double nu = 1;
    Hfunc2(&family_id,&numb_obs, &u_val, &v_val, &rho_value, &nu, &u_double);

    return u_double;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
timesTwo(42)
*/
