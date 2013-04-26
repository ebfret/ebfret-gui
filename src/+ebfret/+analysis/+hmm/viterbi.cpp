#include <assert.h>
#include <math.h>
#include <matrix.h>
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	// copy inputs: px_z, A, pi
	mxArray *ln_px_z_in_m, *ln_A_in_m, *ln_pi_in_m; 
	ln_px_z_in_m = mxDuplicateArray(prhs[0]);
	ln_A_in_m = mxDuplicateArray(prhs[1]);
	ln_pi_in_m = mxDuplicateArray(prhs[2]);

	// get pointers to inputs
	double *ln_px_z, *ln_A, *ln_pi;
    ln_px_z = mxGetPr(ln_px_z_in_m);
    ln_A = mxGetPr(ln_A_in_m);
    ln_pi = mxGetPr(ln_pi_in_m);

	// dimensions: T, K
	const mwSize *dims;
	dims = mxGetDimensions(prhs[0]);
	long T = (long) dims[0];
	long K = (long) dims[1];

	// create output arrays
	plhs[0] = mxCreateDoubleMatrix(T, 1, mxREAL);
	// plhs[1] = mxCreateDoubleMatrix(T, K, mxREAL);
	// plhs[2] = mxCreateDoubleMatrix(T, K, mxREAL);
	double *z_hat = mxGetPr(plhs[0]);
	// double *z_max = mxGetPr(plhs[1]);
	// double *omega = mxGetPr(plhs[2]);

	// z_max(t, k) is most probable state at t-1, given state k at t
	double z_max[T*K];
	// log likehood of states:
	// omega(t, k) =  log(px_z(t,k)) + log(px_z(t-1, z_max(t,k))) 
	//                + omega(t-1, z_max(t,k))
	double omega[T*K];

	// Forward Sweep - Calculate
	//
	//   omega(t, k) = ln_px_z(t, k) + max_l { ln_A(l, k) + omega(t-1, l) }
	//   z_max(t, k) = argmax_l { ln_A(l, k) + omega(t-1, l) }

	// compute first timestep
	long k;
	for (k = 0; k < K; k++) 
	{
		// arbitrary value, since there is no predecessor to t=1
		z_max[k * T] = -1;
		// omega(1, k) = ln(pi(k)) + ln(px_z(1, k))
		omega[k * T] = ln_pi[k] + ln_px_z[k * T];
	}

	// remaining timesteps
	long t;
	long l;
	for (t = 1; t < T; t++)
	{
		for (k = 0; k < K; k++) 
		{
			z_max[k * T + t] = 0;
			omega[k * T + t] = ln_A[k * K] + omega[t - 1];
			double omega_l = 0;
			for (l = 1; l < K; l++)
			{
				omega_l = ln_A[k * K + l] + omega[l * T + t - 1];
				if (omega_l > omega[k * T + t]) 
				{
					z_max[k * T + t] = l;
					omega[k * T + t] = omega_l;
				}
			}
			omega[k * T + t] += ln_px_z[k * T + t];
		}
	}
	
	// Backward sweep 
    // get state for last time step
    double omega_max = omega[T-1]; 
    z_hat[T-1] = 0;
    for (k = 1; k < K; k++)
    {
    	if (omega[k * T + T-1] > omega_max)
    	{
    		omega_max = omega[k * T + T-1];
    		z_hat[T-1] = k;
    	}
    }
    z_hat[T-1];

    // remaining timesteps are determined by z_max
    for (t = T-2; t >= 0; t--)
    {
    	z_hat[t] = z_max[(long)(z_hat[t+1]) * T + t + 1];
    	if (z_hat[t] < 0)
    	{
    		mexPrintf("t=%d, z_hat[t]=%f", t, z_hat[t]);
    	}
    }

    // increment by one to turn into matlab index
    for (t = 0; t<T; t++)
    {
    	z_hat[t] += 1;
    }
	return;
}
