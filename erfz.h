#pragma once

#include <complex>
#include <cmath>
#include <math.h>
#include <assert.h>

inline double sinc(double x)
{
    if (fabs(x) < 1E-10) return 1.;
    else return sin(x) / x;
}

inline double exp_minus(double a, double b)
{
    assert(b >= 0);
    if (b > 0.5) {
        return exp(a + b) - exp(a - b);
    } else {
        return exp(a - b) * expm1(2. * b);
    }
}
// This implementation is following:
// https://granite.phys.s.u-tokyo.ac.jp/svn/LCGT/trunk/sensitivity/Matlab/official/@double/erfz.pdf
// before you do some `smart` optimizations, think about if your optimization will cause an overflow!
inline std::complex<double> erfz(std::complex<double> const &c)
{
    // OK! I think one can stil improve the speed of the code by at least 50%, by some naive optimization.
    typedef std::complex<double> Complex;
    double const Pi = 3.141592653589793238462643383279502884197169399375;
    double const Re = c.real();
    double const R2 = Re * Re;

    double Im_ = c.imag();
    if (Im_ == 0) return std::complex<double>(erf(Re), 0);
    bool const conj = Im_ < 0;
    double const Im = conj ? (-Im_) : Im_;

    double const erfR = erf(Re);


    Complex const E_ = exp(-R2) / (2 * Pi) * 2 * Im *
        Complex(sin(Re * Im) * sinc(Re * Im), sinc(2 * Re * Im));


    // Note: exp(-0.25 * 12 ^ 2) / (0.25 * 12 ^ 2) = 6E-18
    int const N = 15;
    double F_R = 0;
    for (int n = N; n >= 1; --n) {
        F_R += exp(-0.25 * n * n) / (0.25 * n * n + R2);
    }
    F_R *= exp(-R2) * Re / Pi;


    int const M = (int)(2 * Im);
    int const M_N = std::max(M - N, 1);

    Complex HG(0, 0);
    Complex H(0, 0);
    Complex G(0, 0);

    // for H
    for (int n = std::min(M_N - 1, N); n >= 1; --n) {
        H += exp(-0.25 * n * n - n * Im - R2) / (0.25 * n * n + R2) * Complex(Re, 0.5 * n);
    }

    // overlap of H & G
    for (int n = N; n > M_N - 1; --n) {

        double HG_R = (exp(-0.25 * n * n - n * Im - R2) + exp(-0.25 * n * n + n * Im - R2)) / (0.25 * n * n + R2) * Re;
        double HG_I = exp_minus(-0.25 * n * n - R2, n * Im) / (0.25 * n * n + R2) * (-0.5 * n);

        HG.real(HG.real() + HG_R);
        HG.imag(HG.imag() + HG_I);
    }

    // for G
    for (int n = M + N; n > std::max(M_N - 1, N); --n) {
        G += exp(-0.25 * n * n + n * Im - R2) / (0.25 * n * n + R2) * Complex(Re, -0.5 * n);
    }

    H *= 1. / (2 * Pi);
    G *= 1. / (2 * Pi);
    HG *= 1. / (2 * Pi);

    Complex a(cos(-2 * Re * Im), sin(-2 * Re * Im));
    Complex b = a * (H + G + HG);
    double real = erfR + E_.real() + F_R - b.real();
    double imag = 0 + E_.imag() + 0 - b.imag();

    if (conj) return Complex(real, -imag);
    return Complex(real, imag);
}
