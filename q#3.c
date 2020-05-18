#include <stdio.h>
#include <complex.h>
#include<math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])
#define numpts 256

// defining sinc function
double fx(double x)
{
	double f;
	if(x==0.)
		f=1.;
	else
		f=sin(x)/x;
	return f;
}

int main()
{
	double in[2*numpts];
	double xmin,xmax,dx;
	xmin=-50.;
	xmax=50.;
	dx=(xmax-xmin)/(numpts-1);
	double xarr[numpts],yarr[numpts];
	for (int i=0;i<numpts;++i)
	{
		xarr[i]=xmin+i*dx;
		yarr[i]=fx(xarr[i]);
		REAL(in,i)=yarr[i];
		IMAG(in,i)=0.;
	}
	size_t stride=1;
	int status;
	status=gsl_fft_complex_radix2_forward(in,stride,numpts);
	if(status!=GSL_SUCCESS)
		printf("error occured !");
	// normalising dft
	for (int i=0;i<numpts;++i)
	{
		REAL(in,i)*=1./(sqrt(numpts)); 
		IMAG(in,i)*=1./(sqrt(numpts));
	}
	double freq[numpts]; // array for storing values of k-space
	for (int i=0;i<(numpts/2);++i)
	{
		freq[i]=i/(dx);
	}
	for (int i=(numpts/2);i<numpts;++i)
	{
		freq[i]=-(numpts-i)/(dx);
	}
	double _Complex  factor;
	for (int i=0;i<numpts;++i)
	{
		factor=cexp(-2.*M_PI*(1./numpts)*freq[i]*xmin*I);
		REAL(in,i)=REAL(in,i)*creal(factor)-IMAG(in,i)*cimag(factor);
		REAL(in,i)*=dx*sqrt(numpts/(2.*M_PI));
		printf("%.20f\t%.20f\n",2.*M_PI*freq[i]/numpts,REAL(in,i));
	}
	
}
