#pragma once

#include <complex>
#include <cmath>
#include <math.h>

double sinc(double x)
{
	if (x < 1E-10) return 1.;
	else return sin(x) / x;
}

// This implementation is following:
// https://granite.phys.s.u-tokyo.ac.jp/svn/LCGT/trunk/sensitivity/Matlab/official/@double/erfz.pdf
inline std::complex<double> erfz(std::complex<double> const &c)
{
	typedef std::complex<double> Complex;
	double const Pi = 3.141592653589793238462643383279502884197169399375;
	double const Re = c.real();
	double const R2 = Re * Re;
	double Im = c.imag();

	if (Im == 0) return std::complex<double>(erf(Re), 0);

	bool conj = Im < 0;
	if (conj) Im = -Im;

	double erfR = erf(Re);

    Complex E_ = exp(-R2) / (2 * Pi) * 2 * Im *
        Complex(sin(Re * Im) * sinc(Re * Im), sinc(2 * Re * Im));


    // 	exp(-0.25 * 12 ^ 2) / (0.25 * 12 ^ 2) = 6E-18
    int const N = 15;
    double F_R = 0;
    for (int n = N; n >= 1; --n) {
        F_R += exp(-0.25 * n * n) / (0.25 * n * n + R2);
    }
    F_R *= exp(-R2) * Re / Pi;


	int const M = (int)(2 * Im);
	int const M_N = M - N > 1 ? M - N : 1;

	Complex H(0., 0.);
	Complex G(0., 0.);

	for (int n = N; n >= 1; --n) {
        H += exp(-0.25 * n * n - n * Im - R2) / (0.25 * n * n + R2) * Complex(Re, 0.5 * n);
    }

	for (int n = M + N; n >= M_N; --n) {
        //G += exp(-0.25 * (n - 2 * Im) * (n - 2 * Im) + (Im * Im - R2)) / (0.25 * n * n + R2) * Complex(Re, -0.5 * n);
		G += exp(-0.25 * n * n + n * Im - R2) / (0.25 * n * n + R2) * Complex(Re, -0.5 * n);
	}

	H *= 1. / (2 * Pi);
	G *= 1. / (2 * Pi);

	Complex a(cos(-2 * Re * Im), sin(-2 * Re * Im));
	Complex b = a * (G + H);
	double real = erfR + E_.real() + F_R - b.real();
	double imag = 0 + E_.imag() + 0 - b.imag();

	if (conj) return Complex(real, -imag);
	return Complex(real, imag);
}
