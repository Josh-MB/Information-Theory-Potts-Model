#include "../include/psi.hpp"

#include <fmt/format.h>
#include <math.h>
#define el 0.5772156649015329

#if defined(___OSX___) || defined(WIN32)
	#define M_PI 3.14159265358979323846
#endif

double psi(double x)
{
    double s,ps,xa,x2;
    int n,k;
    static double a[] = {
        -0.8333333333333e-01,
         0.83333333333333333e-02,
        -0.39682539682539683e-02,
         0.41666666666666667e-02,
        -0.75757575757575758e-02,
         0.21092796092796093e-01,
        -0.83333333333333333e-01,
         0.4432598039215686};

    xa = fabs(x);
    s = 0.0;
    if ((x == (int)x) && (x <= 0.0)) {
        ps = 1e308;
        return ps;
    }
    if (xa == (int)xa) {
        n = (int)xa;
        for (k=1;k<n;k++) {
            s += 1.0/k;
        }
        ps =  s-el;
    }
    else if ((xa+0.5) == ((int)(xa+0.5))) {
        n = (int)(xa-0.5);
        for (k=1;k<=n;k++) {
            s += 1.0/(2.0*k-1.0);
        }
        ps = 2.0*s-el-1.386294361119891;
    }
    else {
        if (xa < 10.0) {
            n = 10-(int)xa;
            for (k=0;k<n;k++) {
                s += 1.0/(xa+k);
            }
            xa += n;
        }
        x2 = 1.0/(xa*xa);
        ps = log(xa)-0.5/xa+x2*(((((((a[7]*x2+a[6])*x2+a[5])*x2+
            a[4])*x2+a[3])*x2+a[2])*x2+a[1])*x2+a[0]);
        ps -= s;
    }
    if (x < 0.0)
        ps = ps - M_PI*cos(M_PI*x)/sin(M_PI*x)-1.0/x;
    return ps;
}


double digamma(double x) {
	double result = 0, xx, xx2, xx4;
	if(x <= 0) {
		fmt::print("digamma needs positive input\n");
		return 0;
	}
	for(; x < 7; ++x)
		result -= 1 / x;
	x -= 1.0 / 2.0;
	xx = 1.0 / x;
	xx2 = xx*xx;
	xx4 = xx2*xx2;
	result += log(x) + (1. / 24.)*xx2 - (7.0 / 960.0)*xx4 + (31.0 / 8064.0)*xx4*xx2 - (127.0 / 30720.0)*xx4*xx4;
	return result;
}