#include <R.h>
#include <Rmath.h>

void rW(int *n, double *kappa, int *m, double *y)
{
    /* Follow Wood (1994). */
    double l = *kappa;
    double d = *m - 1;
    /* Algebraically equivalent to
         (-2. * l + sqrt(4. * l * l + d * d)) / d
       as in the reference, but numerically more stable:
    */
    double b = d / (sqrt(4. * l * l + d * d) + 2. * l);
    double x = (1. - b) / (1. + b);
    double c = l * x + d * log(1. - x * x);
    double u, w, z;

    Rboolean done;
    int i;

    GetRNGstate();
    for(i = 0; i < *n; i++) {
	done = FALSE;
	while(!done) {
	    z = rbeta(d / 2., d / 2.);
	    w = (1. - (1. + b) * z) / (1. - (1. - b) * z);
	    u = unif_rand();
	    if(l * w + d * log(1. - x * w) - c >= log(u)) {
		done = TRUE;
	    }
	}
	y[i] = w;
    }
    PutRNGstate();
}
