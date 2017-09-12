#include <Rcpp.h>
using namespace Rcpp;
#include <stan/math/rev/core.hpp>

// [[Rcpp::interfaces(r,cpp)]]

void (*Hfunc2) (int* family,int* n,double* v,double* u,double* theta,double* nu,double* out);

void (*diffhfunc_rho_tCopula) (double* u, double* v, int* n, double* param, int* copula, double* out);      // par
void (*diffhfunc_mod) (double* u, double* v, int* n, double* param, int* copula, double* out);              // par
void (*diffhfunc_nu_tCopula_new) (double* u, double* v, int* n, double* param, int* copula, double* out);   // par2
void (*diffhfunc_v_mod) (double* u, double* v, int* n, double* param, int* copula, double* out);            // u2



extern "C" void R_init_TestVine(DllInfo *dll) {
    Hfunc2 = (void (*) (int* ,int* ,double* ,double* ,double* ,double* ,double* )) R_GetCCallable("VineCopula", "Hfunc2");

    diffhfunc_rho_tCopula = (void (*) (double* , double* , int* , double* , int* , double* )) R_GetCCallable("VineCopula", "diffhfunc_rho_tCopula");
    diffhfunc_mod = (void (*) (double* , double* , int* , double* , int* , double* )) R_GetCCallable("VineCopula", "diffhfunc_mod");

    diffhfunc_nu_tCopula_new = (void (*) (double* , double* , int* , double* , int* , double* )) R_GetCCallable("VineCopula", "diffhfunc_nu_tCopula_new");

    diffhfunc_v_mod = (void (*) (double* , double* , int* , double* , int* , double* )) R_GetCCallable("VineCopula", "diffhfunc_v_mod");
}


using namespace stan::math;

// hfunc_trans with 3 vars arguement
stan::math::var hfunc_trans (int family,
                 double u,
                 const stan::math::var& v,
                 const stan::math::var& theta,
                 const stan::math::var& theta2){
    double v_val = v.val();
    double theta_val = theta.val();
    double theta2_val = theta2.val();
    double hfunc2_val;
    double hfunc2_dtheta;
    double hfunc2_dtheta2;
    double hfunc2_dv;


    int nn = 1;

    Hfunc2(&family, &nn, &u, &v_val, &theta_val, &theta2_val, &hfunc2_val);

    double *param = new double[2];
    param[0] = theta_val;
    param[1] = theta2_val;

    diffhfunc_v_mod(&u, &v_val, &nn, param, &family, &hfunc2_dv);

    if (family == 2){
        diffhfunc_rho_tCopula(&u, &v_val, &nn, param, &family, &hfunc2_dtheta);
        diffhfunc_nu_tCopula_new(&u, &v_val, &nn, param, &family, &hfunc2_dtheta2);
        return (new precomp_vvv_vari(hfunc2_val, v.vi_, theta.vi_,theta2.vi_, hfunc2_dv, hfunc2_dtheta, hfunc2_dtheta2) );
    }

    diffhfunc_mod(&u, &v_val, &nn, param, &family, &hfunc2_dtheta);
    return (new precomp_vv_vari(hfunc2_val, v.vi_, theta.vi_, hfunc2_dv, hfunc2_dtheta) );

}

// hfunc_trans with 2 vars arguement
stan::math::var hfunc_trans (int family,
                             double u,
                             const stan::math::var& v,
                             const stan::math::var& theta,
                             double theta2){
    double v_val = v.val();
    double theta_val = theta.val();
    double hfunc2_val;
    double hfunc2_dtheta;
    double hfunc2_dv;


    int nn = 1;

    Hfunc2(&family, &nn, &u, &v_val, &theta_val, &theta2, &hfunc2_val);

    double *param = new double[2];
    param[0] = theta_val;
    param[1] = theta2;

    if (family == 2){
        diffhfunc_rho_tCopula(&u, &v_val, &nn, param, &family, &hfunc2_dtheta);
        return (new precomp_vv_vari(hfunc2_val, v.vi_, theta.vi_,hfunc2_dv, hfunc2_dtheta) );
    }

    diffhfunc_mod(&u, &v_val, &nn, param, &family, &hfunc2_dtheta);
    return (new precomp_vv_vari(hfunc2_val, v.vi_, theta.vi_, hfunc2_dv, hfunc2_dtheta) );
}

// hfunc_trans with 1 vars arguement - v vari
stan::math::var hfunc_trans (int family,
                             double u,
                             const stan::math::var& v,
                             double theta,
                             double theta2){
    double v_val = v.val();
    double hfunc2_val;
    double hfunc2_dv;


    int nn = 1;

    Hfunc2(&family, &nn, &u, &v_val, &theta, &theta2, &hfunc2_val);

    double *param = new double[2];
    param[0] = theta;
    param[1] = theta2;

    diffhfunc_v_mod(&u, &v_val, &nn, param, &family, &hfunc2_dv);

    return (new precomp_v_vari(hfunc2_val, v.vi_, hfunc2_dv) );

}

// hfunc_trans with 1 vars arguement - theta vari
stan::math::var hfunc_trans (int family,
                             double u,
                             double v,
                             const stan::math::var& theta,
                             double theta2){
    double theta_val = theta.val();
    double hfunc2_val;
    double hfunc2_dtheta;


    int nn = 1;

    Hfunc2(&family, &nn, &u, &v, &theta_val, &theta2, &hfunc2_val);

    double *param = new double[2];
    param[0] = theta_val;
    param[1] = theta2;

    if (family == 2){
        diffhfunc_rho_tCopula(&u, &v, &nn, param, &family, &hfunc2_dtheta);
        return (new precomp_v_vari(hfunc2_val, theta.vi_,hfunc2_dtheta) );
    }

    diffhfunc_mod(&u, &v, &nn, param, &family, &hfunc2_dtheta);
    return (new precomp_v_vari(hfunc2_val, theta.vi_, hfunc2_dtheta) );
}

// hfunc_trans with 1 vars arguement - theta2 vari
stan::math::var hfunc_trans (int family,
                             double u,
                             double v,
                             double theta,
                             const stan::math::var& theta2){
    double theta2_val = theta2.val();
    double hfunc2_val;
    double hfunc2_dtheta2;
    double hfunc2_dv;


    int nn = 1;

    Hfunc2(&family, &nn, &u, &v, &theta, &theta2_val, &hfunc2_val);

    double *param = new double[2];
    param[0] = theta;
    param[1] = theta2_val;

    diffhfunc_nu_tCopula_new(&u, &v, &nn, param, &family, &hfunc2_dtheta2);
    return (new precomp_v_vari(hfunc2_val, theta2.vi_, hfunc2_dtheta2) );

}





// [[Rcpp::export]]
void Gaussian_hfunc() {

    stan::math::var rho_value = 0.6;
    stan::math::var theta2 = 0;
    double u_val = 0.2;
    stan::math::var v_val = 0.7;

    int family_id = 1;
    int numb_obs = 1;
    stan::math::var u_double = hfunc_trans(family_id, u_val, v_val, rho_value, theta2);
    u_double.grad();

    Rcpp::Rcout << " value of u trans " << u_double.val() << std::endl;
    Rcpp::Rcout << " deriv of u trans wrt theta " << rho_value.adj() << std::endl;
    Rcpp::Rcout << " deriv of u trans wrt v " << v_val.adj() << std::endl;
}

// [[Rcpp::export]]
void Student_hfunc() {

    stan::math::var rho_value = 0.6;
    stan::math::var theta2 = 5; // nu
    double u_val = 0.2;
    stan::math::var v_val = 0.7;

    int family_id = 2;
    int numb_obs = 1;

    stan::math::var u_double = hfunc_trans(family_id, u_val, v_val, rho_value, theta2);
    u_double.grad();

    Rcpp::Rcout << " value of u trans " << u_double.val() << std::endl;
    Rcpp::Rcout << " deriv of u trans wrt theta " << rho_value.adj() << std::endl;
    Rcpp::Rcout << " deriv of u trans wrt v " << v_val.adj() << std::endl;
    Rcpp::Rcout << " deriv of u trans wrt nu " << theta2.adj() << std::endl;
}

// [[Rcpp::export]]
void Clayton_hfunc() {

    stan::math::var rho_value = 0.6;
    stan::math::var theta2 = 0;
    double u_val = 0.2;
    stan::math::var v_val = 0.7;

    int family_id = 3;
    int numb_obs = 1;
    stan::math::var u_double = hfunc_trans(family_id, u_val, v_val, rho_value, theta2);
    u_double.grad();

    Rcpp::Rcout << " value of u trans " << u_double.val() << std::endl;
    Rcpp::Rcout << " deriv of u trans wrt theta " << rho_value.adj() << std::endl;
    Rcpp::Rcout << " deriv of u trans wrt v " << v_val.adj() << std::endl;
}

// [[Rcpp::export]]
void Gumbel_hfunc() {

    stan::math::var rho_value = 1.6;
    stan::math::var theta2 = 0;
    double u_val = 0.2;
    stan::math::var v_val = 0.7;

    int family_id = 4;
    int numb_obs = 1;
    stan::math::var u_double = hfunc_trans(family_id, u_val, v_val, rho_value, theta2);
    u_double.grad();

    Rcpp::Rcout << " value of u trans " << u_double.val() << std::endl;
    Rcpp::Rcout << " deriv of u trans wrt theta " << rho_value.adj() << std::endl;
    Rcpp::Rcout << " deriv of u trans wrt v " << v_val.adj() << std::endl;
}

// [[Rcpp::export]]
void Frank_hfunc() {

    stan::math::var rho_value = 0.6;
    stan::math::var theta2 = 0;
    double u_val = 0.2;
    stan::math::var v_val = 0.7;

    int family_id = 5;
    int numb_obs = 1;
    stan::math::var u_double = hfunc_trans(family_id, u_val, v_val, rho_value, theta2);
    u_double.grad();

    Rcpp::Rcout << " value of u trans " << u_double.val() << std::endl;
    Rcpp::Rcout << " deriv of u trans wrt theta " << rho_value.adj() << std::endl;
    Rcpp::Rcout << " deriv of u trans wrt v " << v_val.adj() << std::endl;
}


// [[Rcpp::export]]
void Joe_hfunc() {

    stan::math::var rho_value = 1.6;
    stan::math::var theta2 = 0;
    double u_val = 0.2;
    stan::math::var v_val = 0.7;

    int family_id = 6;
    int numb_obs = 1;
    stan::math::var u_double = hfunc_trans(family_id, u_val, v_val, rho_value, theta2);
    u_double.grad();

    Rcpp::Rcout << " value of u trans " << u_double.val() << std::endl;
    Rcpp::Rcout << " deriv of u trans wrt theta " << rho_value.adj() << std::endl;
    Rcpp::Rcout << " deriv of u trans wrt v " << v_val.adj() << std::endl;
}

