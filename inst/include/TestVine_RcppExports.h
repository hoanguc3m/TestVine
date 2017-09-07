// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#ifndef RCPP_TestVine_RCPPEXPORTS_H_GEN_
#define RCPP_TestVine_RCPPEXPORTS_H_GEN_

#include <Rcpp.h>

namespace TestVine {

    using namespace Rcpp;

    namespace {
        void validateSignature(const char* sig) {
            Rcpp::Function require = Rcpp::Environment::base_env()["require"];
            require("TestVine", Rcpp::Named("quietly") = true);
            typedef int(*Ptr_validate)(const char*);
            static Ptr_validate p_validate = (Ptr_validate)
                R_GetCCallable("TestVine", "_TestVine_RcppExport_validate");
            if (!p_validate(sig)) {
                throw Rcpp::function_not_exported(
                    "C++ function with signature '" + std::string(sig) + "' not found in TestVine");
            }
        }
    }

    inline NumericVector timesTwo(NumericVector x) {
        typedef SEXP(*Ptr_timesTwo)(SEXP);
        static Ptr_timesTwo p_timesTwo = NULL;
        if (p_timesTwo == NULL) {
            validateSignature("NumericVector(*timesTwo)(NumericVector)");
            p_timesTwo = (Ptr_timesTwo)R_GetCCallable("TestVine", "_TestVine_timesTwo");
        }
        RObject rcpp_result_gen;
        {
            RNGScope RCPP_rngScope_gen;
            rcpp_result_gen = p_timesTwo(Shield<SEXP>(Rcpp::wrap(x)));
        }
        if (rcpp_result_gen.inherits("interrupted-error"))
            throw Rcpp::internal::InterruptedException();
        if (rcpp_result_gen.inherits("try-error"))
            throw Rcpp::exception(as<std::string>(rcpp_result_gen).c_str());
        return Rcpp::as<NumericVector >(rcpp_result_gen);
    }

}

#endif // RCPP_TestVine_RCPPEXPORTS_H_GEN_
