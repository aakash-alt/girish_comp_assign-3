#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <complex.h>
#include <fftw3.h>
#include<math.h>
#define N 256
#define REAL 0
#define IMAG 1

double fx(double x) // defining sinc function
{
	double f;
	if (x==0.)
		f=1.;
	else
		f=sin(x)/x;
	return f;
}
double fk(double k) // defining Fourier transform of sinc function
{
	double f;
	if(k>-1&& k<1)
		f=sqrt(M_PI/2.);
	else
		f=0.;
	return f;
}
int main()
{
	double xmin,xmax,dx;
	xmin=-50.;
	xmax=50.;
	dx=(xmax-xmin)/(N-1);
	double x[N],f[N];
	for(int i=0;i<N;++i)
	{
		x[i]=xmin+i*dx;
		f[i]=fx(x[i]);
	}
	fftw_complex in[N], out[N];
	//double in[N];
	for (int i=0;i<N;++i)
	{
		in[i][REAL]=f[i];
		in[i][IMAG]=0.;
	}
	
	fftw_plan fft = fftw_plan_dft_1d(N,in,out,FFTW_FORWARD,FFTW_ESTIMATE);
	fftw_execute(fft);
	fftw_destroy_plan(fft);
	fftw_cleanup();
	for (int i=0;i<N;++i)
	{
		out[i][REAL]*=1./(sqrt(N)); // normalising dft
		out[i][IMAG]*=1./(sqrt(N));
	}
	double freq[N]; // array for storing values of k-space
	for (int i=0;i<int(N/2);++i)
	{
		freq[i]=i/(dx);
	}
	for (int i=int(N/2);i<N;++i)
	{
		freq[i]=-(N-i)/(dx);
	}
	
	double _Complex  factor;
	for(int i=0;i<N;i++)
	{ 
		factor=cexp(-2.*M_PI*(1./N)*freq[i]*xmin*I);
		out[i][REAL]=out[i][REAL]*creal(factor)-out[i][IMAG]*cimag(factor);
		out[i][REAL]*=dx*sqrt(N/(2.*M_PI));
		
		out[i][IMAG]=out[i][REAL]*cimag(factor)+out[i][IMAG]*creal(factor);
		out[i][IMAG]*=dx*sqrt(N/(2.*M_PI));
		printf("%.20f\t%.20f\n",2.*M_PI*freq[i]/N,out[i][REAL]);
	}
	
}
