void my0F1(int *n, double *x, double *nu, double *y0, double *y)
{
    double r, s, t, u, v, z;
    int i;

    for(i = 0; i < *n; i++) {

	z = x[i];
	t = y0[i];
	u = nu[i];
	v = 1.;
	r = 0.;
	s = t;

	/* Iterate until numeric convergence if possible. */

	while((s > r) && (v <= 500)) {
	    r = s;
	    t *= z / (u * v);
	    s += t;	    
	    u++;
	    v++;
	}

	y[i] = s;
	
    }
}
