#include <R.h>

void mycfG(int *len, double *z, double *n, double *y)
{
    /* Compute
         I_{\nu}(z) / I_{\nu - 1}(z)
       via the Gauss continued fraction, implemented using Eqn 3.2'
       in Gautschi and Slavik (1978), modified to use a_0 = z/(2\nu)
       instead of a_0 = 1.
       See also "Forward Series Recurrence Algorithm" in
       http://dlmf.nist.gov/3.10#iii.
    */

    double tol = 1e-10;
    double zi, nu, xG, rho, p, s, t, u, v;
    int i, k;

    for(i = 0; i < *len; i++) {
	zi = z[i];
	nu = n[i];
	xG = zi * zi / 4.;
	p = s = zi / (2. * nu);
	rho = 0.;
	u = nu * (nu - 1.);
	v = 2. * nu;
	k = 1;
	while(fabs(p) > tol * s) {
	    u += v;
	    v += 2.;
	    t = xG * (1. + rho);
	    rho = - t / (u + t);
	    p *= rho;
	    s += p;
	    /* REprintf("k: %d p: %g s: %g\n", k, p, s); */
	    k++;
	}
	y[i] = s;
    }
}

void mycfP(int *len, double *z, double *n, double *y)
{
    /* Compute
         I_{\nu}(z) / I_{\nu - 1}(z)
       via the Perron continued fraction, implemented using Eqn 3.3'
       in Gautschi and Slavik (1978), modified to use a_0 = z/(z+2\nu)
       instead of a_0 = 1.
       See also "Forward Series Recurrence Algorithm" in
       http://dlmf.nist.gov/3.10#iii
    */

    double tol = 1e-10;    
    double zi, nu, xP, rho, p, s, t, u, v, w;
    int i, k;

    for(i = 0; i < *len; i++) {
	zi = z[i];
	nu = n[i];
	xP = zi / 2.;
	p = s = zi / (zi + 2. * nu);
	v = nu + zi + .5;
	u = (nu + zi) * v;
	w = xP * (nu + .5);
	rho = w / ((nu + xP) * v - w);
	p *= rho;
	s += p;
	k = 2;
	while(fabs(p) > tol * s) {
	    u += v;
	    v += .5;
	    w += xP;
	    t = w * (1. + rho);
	    rho = t / (u - t);
	    p *= rho;
	    s += p;
	    /* REprintf("k: %d p: %g s: %g\n", k, p, s); */
	    k++;
	}
	y[i] = s;
    }
}
