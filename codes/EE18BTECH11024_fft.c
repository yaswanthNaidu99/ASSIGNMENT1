#include <stdio.h>
#include <complex.h>
#include <math.h>
typedef double complex cd;

void fft(complex *X, int n)
{
	if(n == 1) return;
	
	cd X_odd[n/2], X_even[n/2];
	
	for(int j = 0; j < n/2; j++) 
	{
		X_even[j] = X[2*j];
		X_odd[j] = X[2*j+1];
	}
	
	fft(X_even, n/2);	 	
	fft(X_odd, n/2);		
	
	cd exp;
    for (int j = 0;j < n/2;j++) 
    {
		exp = CMPLX(cos(2*M_PI*j/n),-sin(2*M_PI*j/n));
	    X[j] = X_even[j] + exp * X_odd[j];
        X[j+n/2] = X_even[j] - exp * X_odd[j];
      
    }
}

int main()
{	
	int n = 8;  
	double x[] = {1,2,3,4,1,2,3,4};
	
	cd X[n];
	
	for(int j = 0; j < n; j++) 
	{
		X[j] = x[j];
	}
	fft(X, n);
	
	for(int j = 0; j < n; j++) 
		printf("(%.5lf, %.5lf)\n", creal(X[j]), cimag(X[j]));
}
