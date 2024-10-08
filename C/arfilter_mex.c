// Use Makefile in this directory to build .mex

/* WARNING: THIS FUNCTION PERFORMS NO ERROR CHECKING WHATSOEVER! Don't even
 * think of calling arfilter_mex itself! It is designed to be called via the
 * arfilter.m Matlab wrapper function (utils subfolder); all error checking and
 * documentation happens there. To compile, see Makefile in this folder */

#include "mex.h"
#include <string.h> // for memcpy

#define UNUSED __attribute__ ((unused))

void mexFunction(int UNUSED nlhs, mxArray *plhs[], int UNUSED nrhs, const mxArray *prhs[])
{
	const double* const a    = mxGetPr(prhs[0]);
	const bool          avec = (bool)(*mxGetLogicals(prhs[1]));
	const size_t        p    = (size_t)mxGetScalar(prhs[2]);
	const double* const x    = mxGetPr(prhs[3]);
	const bool          allt = (bool)(*mxGetLogicals(prhs[4]));

	const size_t n = mxGetM(prhs[3]);
	const size_t m = mxGetN(prhs[3]);

	double* const y = mxGetPr(plhs[0] = mxCreateDoubleMatrix(n,m,mxREAL)); /* allocate output array y  */

	const size_t nsq = n*n;

	memcpy(y,x,n*m*sizeof(double)); // copy input x to output y

	if (allt) {
		const size_t mmp = m < p ? m : p;
		if (avec) {
			for (size_t t=0; t<mmp; ++t) {
				double* const yt = y+n*t;
				for (size_t k=0; k<t; ++k) {
					const double* const ytlagk = yt-n*(k+1); // lag is k+1 !!!
					for (size_t i=0; i<n; ++i) yt[i] += a[k]*ytlagk[i];
				}
			}
			for (size_t t=mmp; t<m; ++t) {
				double* const yt = y+n*t;
				for (size_t k=0; k<p; ++k) {
					const double* const ytlagk = yt-n*(k+1); // lag is k+1 !!!
					for (size_t i=0; i<n; ++i) yt[i] += a[k]*ytlagk[i];
				}
			}
		}
		else { // ~avec
			for (size_t t=0; t<mmp; ++t) {
				double* const yt = y+n*t;
				for (size_t k=0; k<t; ++k) {
					const double* const ak = a+nsq*k;
					const double* const ytlagk = yt-n*(k+1); // lag is k+1 !!!
					for (size_t i=0; i<n; ++i) for (size_t j=0; j<n; ++j) yt[i] += ak[i+n*j]*ytlagk[j];
				}
			}
			for (size_t t=mmp; t<m; ++t) {
				double* const yt = y+n*t;
				for (size_t k=0; k<p; ++k) {
					const double* const ak = a+nsq*k;
					const double* const ytlagk = yt-n*(k+1); // lag is k+1 !!!
					for (size_t i=0; i<n; ++i) for (size_t j=0; j<n; ++j) yt[i] += ak[i+n*j]*ytlagk[j];
				}
			}
		}
	}
	else { // ~allt (can assume m > p)
		if (avec) {
			for (size_t t=p; t<m; ++t) {
				double* const yt = y+n*t;
				for (size_t k=0; k<p; ++k) {
					const double* const ytlagk = yt-n*(k+1); // lag is k+1 !!!
					for (size_t i=0; i<n; ++i) yt[i] += a[k]*ytlagk[i];
				}
			}
		}
		else { // ~avec
			for (size_t t=p; t<m; ++t) {
				double* const yt = y+n*t;
				for (size_t k=0; k<p; ++k) {
					const double* const ak = a+nsq*k;
					const double* const ytlagk = yt-n*(k+1); // lag is k+1 !!!
					for (size_t i=0; i<n; ++i) for (size_t j=0; j<n; ++j) yt[i] += ak[i+n*j]*ytlagk[j];
				}
			}
		}
	}
}
