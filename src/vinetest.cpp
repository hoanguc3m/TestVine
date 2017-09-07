#include <Rcpp.h>
using namespace Rcpp;

#include "hfunc.h"

// [[Rcpp::interfaces(r,cpp)]]


// [[Rcpp::export]]
NumericVector timesTwo(NumericVector x) {

    double rho_value = 0.5;
    double u_val = 0.1;
    double v_val = 0.5;
    double u_double;

    int family_id = 1;
    int numb_obs = 1;
    double nu = 1;
    Hfunc2(&family_id,&numb_obs, &u_val, &v_val, &rho_value, &nu, &u_double);
    return x * 2;
    // return u_double;
}

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
timesTwo(42)
*/
