#include <R.h>
#include <Rmath.h>

void rW(int *n, double *kappa, int *m, double *y)
{
    /* Follow Wood (1994). */
    double l = *kappa;
    double d = *m - 1;
    double b = (- 2. * l + sqrt(4. * l * l + d * d)) / d;
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
	    u = runif(0., 1.);
	    if(l * w + d * log(1. - x * w) - c >= log(u)) {
		done = TRUE;
	    }
	}
	y[i] = w;
    }
    PutRNGstate();
}
