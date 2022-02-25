#include <math.h>

#include <R.h>
#include <Rinternals.h>

SEXP FastPOI_C_sub(SEXP in_A, SEXP in_B, SEXP in_Delta, SEXP in_lambda, SEXP in_maxit) {
    int i, j, k, g;

    int p = nrows(in_Delta);
    int d = ncols(in_Delta);
    int maxit = asInteger(in_maxit);

    SEXP out_Z = PROTECT(allocMatrix(REALSXP, p, d));
    double* Z = REAL(out_Z);
    double* Zold = (double*)R_alloc(p * d, sizeof(double));
    double* Delta = REAL(in_Delta);
    double* a = (double*)R_alloc(d, sizeof(double));
    double* A = REAL(in_A);
    double* B = REAL(in_B);
    double a_norm;
    double lambda = asReal(in_lambda);
    double scale;
    double res;

    // Set initial value.
    for (j = 0; j < p * d; ++j) {
        Zold[j] = Z[j] = Delta[j];
    }

    for (i = 0; i < maxit; ++i) {
        // Store current value in Z 'old'.
        // Cyclic updating variables.
        for (g = 0; g < p; ++g) {
            for (j = 0; j < d; ++j) {
                a[j] = Delta[j * p + g];
                for (k = 0; k < p; ++k) {
                    if (k != g) {
                        a[j] -= B[k * p + g] * Z[j * p + k];
                    }
                }
            }

            a_norm = a[0] * a[0];
            for (j = 1; j < d; ++j) {
                a_norm += a[j] * a[j];
            }
            a_norm = sqrt(a_norm);

            if (a_norm > lambda) {
                scale = (1.0 - (lambda / a_norm)) / B[g * p + g];
            } else {
                scale = 0.0;
            }
            for (j = 0; j < d; ++j) {
                Z[j * p + g] = scale * a[j];
            }

        }

        // Copy Z to Zold and check break condition.
        res = 0;
        for (j = 0; j < p * d; ++j) {
            res += (Z[j] - Zold[j]) * (Z[j] - Zold[j]);
            Zold[j] = Z[j];
        }
        if (res < 1e-6) {
            break;
        }

    }

    UNPROTECT(1);
    return out_Z;
}
