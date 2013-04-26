#include <math.h>
#include <matrix.h>
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// copy inputs: px_z, A, pi
	mxArray *px_z_in_m, *A_in_m, *pi_in_m; 
	px_z_in_m = mxDuplicateArray(prhs[0]);
	A_in_m = mxDuplicateArray(prhs[1]);
	pi_in_m = mxDuplicateArray(prhs[2]);

	// get pointers to inputs
	double *px_z, *A, *pi;
    px_z = mxGetPr(px_z_in_m);
    A = mxGetPr(A_in_m);
    pi = mxGetPr(pi_in_m);

	// dimensions: T, K
	const mwSize *dims;
	dims = mxGetDimensions(prhs[0]);
	int T = (int) dims[0];
	int K = (int) dims[1];

	// local variables: a, b, c
	double a[T*K];
	double b[T*K];
	double c[T]; 

	// initialize to zero
	int i;
	for (i=0; i<T*K; i++)
	{
		a[i] = 0;
		b[i] = 0;
	}
	for (i=0; i<T; i++)
	{
		c[i] = 0;
	}

	// Forward Sweep - Calculate
	//
	//   a(t, k)  =  sum_l px_z(t,k) A(l, k) alpha(t-1, l)  
	//   c(t)     =  sum_k a(t, k)
	//
	// and normalize 
	//
	//   a(t, k)  /=  c(t)

	// a(0, k)  =  px_z(0, k) pi(k)
	int k;
	for (k = 0; k < K; k++) 
	{
		a[k*T] = pi[k] * px_z[k*T];
		c[0] += a[k*T];
	}
	// normalize a(0,k) by c(k)
	for (k = 0; k < K; k++) 
	{
		a[k*T] /= c[0];
	}

	int t = 0;
	int l;
	for (t = 1; t < T; t++)
	{
		// a(t, k)  =  sum_l px_z(t,k) A(l, k) alpha(t-1, l)  
		for (k = 0; k < K; k++) 
		{
			for (l = 0; l < K; l++) 
			{
				// a(t,k) +=  px_z(t,k) A(l, k) alpha(t-1, l)  
				a[k*T + t] += px_z[k*T + t] * A[k*K + l] * a[l*T + t-1];
			}			
			// c(t) += a(t,k)
			c[t] += a[k*T + t];
		}
		// normalize a(t,k) by c(t)
		for (k = 0; k < K; k++) 
		{
			a[k*T + t] /= c[t];
		}
	}

		
	// Back sweep - calculate
	//
	// b(t,k)  =  1/c(t+1) sum_l px_z(t+1, l) A(k, l) beta(t+1, l) 

	// b(T-1,k) = 1
	for (k = 0; k < K; k++) 
	{
		b[k*T + (T-1)] = 1;
	}

	// t = T-2:0
	for (t = T-2; t >= 0; t--)
	{
		// b(t, k)  =  sum_l px_z(t+1,l) A(k, l) beta(t+1, l)  
		for (k = 0; k < K; k++) 
		{
			for (l = 0; l < K; l++) 
			{
				// b(t ,k) += px_z(t+1, l) A(k, l) betal(t+1, l)  
				b[k*T + t] += px_z[l*T + t+1] * A[l*K + k] * b[l*T + t+1];
			}			
			// normalize b(t,k) by c(t+1)
			b[k*T + t] /= c[t+1];
		}
	}

	// allocate outputs: g, xi, lnZ
	int xi_dims[3] = {T-1, K, K}; 
	mxArray *g_out_m, *xi_out_m, *lnZ_out_m; 
	plhs[0] = mxCreateDoubleMatrix(T, K, mxREAL);
	g_out_m = plhs[0];
	plhs[1] = mxCreateNumericArray(3, xi_dims, mxDOUBLE_CLASS, mxREAL);
	xi_out_m = plhs[1];
	plhs[2] = mxCreateDoubleMatrix(1, 1, mxREAL);
	lnZ_out_m = plhs[2];

	// get pointers to outputs
	double *g, *xi, *lnZ;
    g = mxGetPr(g_out_m);
    xi = mxGetPr(xi_out_m);
    lnZ = mxGetPr(lnZ_out_m);

	// g(t,k) = a(t,k) * b(t,k)
	for (i=0; i<T*K; i++)
	{
		g[i] = a[i] * b[i];
	}

	// xi(t, k, l) = alpha(t, k) A(k,l) px_z(t+1, l) beta(t+1, l) / c(t+1)
	for (t = 0; t < T-1; t++)
	{
		for (k = 0; k < K; k++) 
		{
			for (l = 0; l < K; l++) 
			{
				xi[l*K*(T-1) + k*(T-1) + t] = (a[k*T + t] \
				                               * A[l*K + k] \
				                               * px_z[l*T + t+1] \
				                               * b[l*T + t+1]) / c[t+1];
			}			
		}
	}

	// ln_Z = sum_t log(c[t])
	lnZ[0] = 0;
	for (t=0; t<T; t++)
	{
		lnZ[0] += log(c[t]);
	}

	return;
}
