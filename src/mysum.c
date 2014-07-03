#include <Rmath.h>

void my0F1(int *n, double *x, double *nu, double *y0, double *y)
{
    double m, r, s, t, u, v, z;
    int i;

    for(i = 0; i < *n; i++) {

	z = x[i];
	t = y0[i];
	u = nu[i];
	v = 1.;
	r = 0.;
	s = t;

	/* Iterate until terms reach maximum. */
	m = (- u - 1 + sqrt((u - 1) * (u - 1)  + 4 * z)) / 2 + 1;
	while(v < m) {
	    r = s;
	    t *= z / (u * v);
	    s += t;	    
	    u++;
	    v++;
	}

	/* Iterate until numeric convergence if possible. */
	while(s > r) {
	    r = s;
	    t *= z / (u * v);
	    s += t;	    
	    u++;
	    v++;
	}

	y[i] = s;
	
    }
}
