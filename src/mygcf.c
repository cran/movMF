#include <Rmath.h>

void mygcf(int *len, double *x, double *s, double *y)
{
    /* Compute
         0F1(; s + 1; z) / 0F1( ; s; z)
       via the Gauss continued fraction
         0F1(; s + 1; z) / (s * 0F1( ; s; z))
	   = 1 / (s + z / (s + 1 + z / (s + 2 + z / ...)))
       using Steed's method (e.g., http://dlmf.nist.gov/3.10):
       Writing generalized continued fractions as
         b_0 + a_1 / (b_1 + a_2 / (b_2 + ...))
       we have
         a_1 = 1, a_2 = a_3 = ... = z;
	 b_0 = 0, b_1 = s, b_2 = s + 1, ...
       and hence the starting values
         C_0 = 0
	 D_1 = 1 / s
	 Delta_1 = 1 / s
	 C_1 = 1 / s
       and the recursion
         D_n = 1 / (D_{n-1} z + b_n)
	 Delta_n = (b_n D_n - 1) Delta_{n-1}
	 C_n = C_{n-1} + Delta_n
       for n >= 2.	 
    */

    double b, z, C, D, Delta;
    int i, n;

    for(i = 0; i < *len; i++) {
	z = x[i];
	b = s[i];
	D = Delta = C = 1 / b;
	n = 0;
	while((fabs(Delta / C) > 1e-12) && (n < 500)) {
	    n++;
	    b++;
	    D = 1 / (D * z + b);
	    Delta = (b * D - 1.) * Delta;
	    C = C + Delta;
	}
	y[i] = s[i] * C;
    }

}
