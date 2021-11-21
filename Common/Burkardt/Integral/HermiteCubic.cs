using System;
using Burkardt.Types;

namespace Burkardt.IntegralNS;

public static class HermiteCubic
{
    public static double hermite_cubic_integral(double x1, double f1, double d1, double x2,
            double f2, double d2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HERMITE_CUBIC_INTEGRAL returns the integral of a Hermite cubic polynomial.
        //
        //  Discussion:
        //
        //    The integral is taken over the definition interval [X1,X2].
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 February 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Fred Fritsch, Ralph Carlson,
        //    Monotone Piecewise Cubic Interpolation,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 17, Number 2, April 1980, pages 238-246.
        //
        //  Parameters:
        //
        //    Input, double X1, F1, D1, the left endpoint, function value
        //    and derivative.
        //
        //    Input, double X2, F2, D2, the right endpoint, function value
        //    and derivative.
        //
        //    Output, double HERMITE_CUBIC_INTEGRAL, the integral of the Hermite
        //    cubic polynomial over the interval X1 <= X <= X2.
        //
    {
        double h = x2 - x1;

        double q = 0.5 * h * (f1 + f2 + h * (d1 - d2) / 6.0);

        return q;
    }

    public static double hermite_cubic_integrate(double x1, double f1, double d1, double x2,
            double f2, double d2, double a, double b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HERMITE_CUBIC_INTEGRATE integrates a Hermite cubic polynomial from A to B.
        //
        //  Discussion:
        //
        //    A and B may be scalars, or one may be a vector, or both
        //    may be vectors of the same size.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 February 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Fred Fritsch, Ralph Carlson,
        //    Monotone Piecewise Cubic Interpolation,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 17, Number 2, April 1980, pages 238-246.
        //
        //  Parameters:
        //
        //    Input, double X1, F1, D1, the left endpoint, function value
        //    and derivative.
        //
        //    Input, double X2, F2, D2, the right endpoint, function value
        //    and derivative.
        //
        //    Input, double A, B, the left and right endpoints of the interval
        //    of integration.
        //
        //    Output, double HERMITE_CUBIC_INTEGRATE, the integral of the
        //    Hermite cubic polynomial over the interval A <= X <= B.
        //
    {
        double h = x2 - x1;

        double ta1 = (a - x1) / h;
        double ta2 = (x2 - a) / h;
        double tb1 = (b - x1) / h;
        double tb2 = (x2 - b) / h;

        double ua1 = ta1 * ta1 * ta1;
        double phia1 = ua1 * (2.0 - ta1);
        double psia1 = ua1 * (3.0 * ta1 - 4.0);

        double ua2 = ta2 * ta2 * ta2;
        double phia2 = ua2 * (2.0 - ta2);
        double psia2 = -ua2 * (3.0 * ta2 - 4.0);

        double ub1 = tb1 * tb1 * tb1;
        double phib1 = ub1 * (2.0 - tb1);
        double psib1 = ub1 * (3.0 * tb1 - 4.0);

        double ub2 = tb2 * tb2 * tb2;
        double phib2 = ub2 * (2.0 - tb2);
        double psib2 = -ub2 * (3.0 * tb2 - 4.0);

        double fterm = f1 * (phia2 - phib2) + f2 * (phib1 - phia1);
        double dterm = (d1 * (psia2 - psib2) + d2 * (psib1 - psia1)) * (h / 6.0);

        double q = 0.5 * h * (fterm + dterm);

        return q;
    }

    public static double[] hermite_cubic_lagrange_integral(double x1, double x2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HERMITE_CUBIC_LAGRANGE_INTEGRAL: Hermite cubic Lagrange integrals.
        //
        //  Discussion:
        //
        //    The Hermite cubic polynomial P(X) for interval (X1,X2) and data
        //    (F1,D1,F2,D2) satisfies:
        //
        //      P(X1) = F1,
        //      P'(X1) = D1,
        //      P(X2) = F2,
        //      P'(X2) = D2.
        //
        //    We can determine four Lagrange polynomials L1(X) through L4(X) so that
        //
        //      P(X) = F1 * L1(X) + D1 * L2(X) + F2 * L3(X) + D2 * L4(X).
        //
        //    This function returns the integrals of these four polynomials over
        //    the domain of definition [X1,X2].
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 February 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Fred Fritsch, Ralph Carlson,
        //    Monotone Piecewise Cubic Interpolation,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 17, Number 2, April 1980, pages 238-246.
        //
        //  Parameters:
        //
        //    Input, double X1, X2, the endpoints.
        //
        //    Output, double HERMITE_CUBIC_LAGRANGE_INTEGRAL[4], the integrals of the
        //    Hermite cubic Lagrange polynomials from X1 to X2.
        //
    {
        double[] q = new double[4];

        double h = x2 - x1;

        q[0] = h / 2.0;
        q[1] = h * h / 12.0;
        q[2] = h / 2.0;
        q[3] = -h * h / 12.0;

        return q;
    }

    public static double[] hermite_cubic_lagrange_integrate(double x1, double x2, double a,
            double b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HERMITE_CUBIC_LAGRANGE_INTEGRATE integrates Hermite cubic Lagrange polynomials.
        //
        //  Discussion:
        //
        //    A and B may be scalars, or one may be a vector, or both
        //    may be vectors of the same size.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 February 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Fred Fritsch, Ralph Carlson,
        //    Monotone Piecewise Cubic Interpolation,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 17, Number 2, April 1980, pages 238-246.
        //
        //  Parameters:
        //
        //    Input, double X1, X2, the endpoints of the interval of definition.
        //
        //    Input, double A, B, the left and right endpoints of the interval
        //    of integration.
        //
        //    Output, double HERMITE_CUBIC_LAGRANGE_INTEGRATE[4], the integrals of the
        //    Hermite cubic Lagrange polynomials over the interval A <= X <= B.
        //
    {
        double h = x2 - x1;
        double ta1 = (a - x1) / h;
        double ta2 = (x2 - a) / h;
        double tb1 = (b - x1) / h;
        double tb2 = (x2 - b) / h;

        double ua1 = ta1 * ta1 * ta1;
        double phia1 = ua1 * (2.0 - ta1);
        double psia1 = ua1 * (3.0 * ta1 - 4.0);

        double ua2 = ta2 * ta2 * ta2;
        double phia2 = ua2 * (2.0 - ta2);
        double psia2 = -ua2 * (3.0 * ta2 - 4.0);

        double ub1 = tb1 * tb1 * tb1;
        double phib1 = ub1 * (2.0 - tb1);
        double psib1 = ub1 * (3.0 * tb1 - 4.0);

        double ub2 = tb2 * tb2 * tb2;
        double phib2 = ub2 * (2.0 - tb2);
        double psib2 = -ub2 * (3.0 * tb2 - 4.0);

        double[] q = new double[4];

        q[0] = 0.5 * h * (phia2 - phib2);
        q[1] = 0.5 * h * (psia2 - psib2) * (h / 6.0);
        q[2] = 0.5 * h * (phib1 - phia1);
        q[3] = 0.5 * h * (psib1 - psia1) * (h / 6.0);

        return q;
    }

    public static void hermite_cubic_lagrange_value(double x1, double x2, int n, double[] x,
            ref double[] f, ref double[] d, ref double[] s, ref double[] t )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HERMITE_CUBIC_LAGRANGE_VALUE evaluates the Hermite cubic Lagrange polynomials.
        //
        //  Discussion:
        //
        //    The Hermite cubic polynomial P(X) for interval (X1,X2) and data
        //    (F1,D1,F2,D2) satisfies:
        //
        //      P(X1) = F1,
        //      P'(X1) = D1,
        //      P(X2) = F2,
        //      P'(X2) = D2.
        //
        //    We can determine four Lagrange polynomials L1(X) through L4(X) so that
        //
        //      P(X) = F1 * L1(X) + D1 * L2(X) + F2 * L3(X) + D2 * L4(X).
        //
        //    This function returns the values and derivatives of these four
        //    polynomials.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 February 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Fred Fritsch, Ralph Carlson,
        //    Monotone Piecewise Cubic Interpolation,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 17, Number 2, April 1980, pages 238-246.
        //
        //  Parameters:
        //
        //    Input, double X1, X2, the endpoints.
        //
        //    Input, int N, the number of sample points.
        //
        //    Input, double X[N], the sample points.
        //
        //    Output, double F[4*N], D[4*N], S[4*N], T[4*N], the value
        //    and first three derivatives of the Hermite cubic Lagrange polynomials
        //    at X.
        //
    {
        int j;

        double h = x2 - x1;

        for (j = 0; j < n; j++)
        {
            double dx = x[j] - x1;
            //
            //  F1.
            //
            f[0 + j * 4] = 1.0 + dx * dx / (h * h) * (-3.0 + dx / h * 2.0);
            d[0 + j * 4] = dx / (h * h) * (-6.0 + dx / h * 6.0);
            s[0 + j * 4] = 1.0 / (h * h) * (-6.0 + dx / h * 12.0);
            t[0 + j * 4] = 1.0 / (h * h * h) * 12.0;
            //
            //  D1
            //
            f[1 + j * 4] = dx + dx * dx / h * (-2.0 + dx / h);
            d[1 + j * 4] = 1.0 + dx / h * (-4.0 + dx / h * 3.0);
            s[1 + j * 4] = 1.0 / h * (-4.0 + dx / h * 6.0);
            t[1 + j * 4] = 1.0 / (h * h) * 6.0;
            //
            //  F2
            //
            f[2 + j * 4] = dx * dx / (h * h) * (3.0 - 2.0 * (dx / h));
            d[2 + j * 4] = dx / (h * h) * (6.0 - 6.0 * (dx / h));
            s[2 + j * 4] = 1.0 / (h * h) * (6.0 - 12.0 * (dx / h));
            t[2 + j * 4] = 1.0 / (h * h * h) * -12.0;
            //
            //  D2
            //
            f[3 + j * 4] = dx * dx / h * (-1.0 + dx / h);
            d[3 + j * 4] = dx / h * (-2.0 + dx / h * 3.0);
            s[3 + j * 4] = 1.0 / h * (-2.0 + dx / h * 6.0);
            t[3 + j * 4] = 1.0 / h * 6.0;
        }
    }

    public static double hermite_cubic_spline_integral(int nn, double[] xn, double[] fn,
            double[] dn )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HERMITE_CUBIC_SPLINE_INTEGRAL: Hermite cubic spline integral.
        //
        //  Discussion:
        //
        //    The integral is taken over the definition interval ( X[0], X[NN-1] ).
        //
        //    Note that if the intervals are equal in size, then the derivative
        //    information in DN has no effect on the integral value,
        //    except for the first and last entries.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 February 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Fred Fritsch, Ralph Carlson,
        //    Monotone Piecewise Cubic Interpolation,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 17, Number 2, April 1980, pages 238-246.
        //
        //  Parameters:
        //
        //    Input, int NN, the number of data points.
        //
        //    Input, double XN[NN], the coordinates of the data points.
        //    The entries in XN must be in strictly ascending order.
        //
        //    Input, double FN[NN], the function values.
        //
        //    Input, double DN[NN], the derivative values.
        //
        //    Output, double HERMITE_CUBIC_SPLINE_INTEGRAL, the integral of the
        //    Hermite cubic spline over the interval X[0] <= X <= X[NN-1].
        //
    {
        int i;

        double q = 0.0;

        for (i = 0; i < nn - 1; i++)
        {
            q += 0.5 * (xn[i + 1] - xn[i]) * (fn[i] + fn[i + 1]
                                                    + (xn[i + 1] - xn[i]) * (dn[i] - dn[i + 1]) / 6.0);
        }

        return q;
    }

    public static double[] hermite_cubic_spline_integrate(int nn, double[] xn, double[] fn,
            double[] dn, int n, double[] a, double[] b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HERMITE_CUBIC_SPLINE_INTEGRATE integrate Hermite cubic spline over [A,B].
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 February 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Fred Fritsch, Ralph Carlson,
        //    Monotone Piecewise Cubic Interpolation,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 17, Number 2, April 1980, pages 238-246.
        //
        //  Parameters:
        //
        //    Input, int NN, the number of data points.
        //
        //    Input, double XN[NN], the coordinates of the data points.
        //    The entries in XN must be in strictly ascending order.
        //
        //    Input, double FN[NN], the function values.
        //
        //    Input, double DN[NN], the derivative values.
        //
        //    Input, int N, the number of integration intervals.
        //
        //    Input, double A[N], B[N], the integration endpoints.
        //
        //    Output, double HERMITE_CUBIC_SPLINE_INTEGRATE[N], the integral
        //    over the interval [A,B].
        //
    {
        int ii;

        double[] q = new double[n];

        int i = n / 2;
        int j = n / 2;

        for (ii = 0; ii < n; ii++)
        {
            q[ii] = 0.0;

            double bb;
            double aa;
            double s;
            if (a[ii] <= b[ii])
            {
                aa = a[ii];
                bb = b[ii];
                s = +1.0;
            }
            else
            {
                aa = b[ii];
                bb = a[ii];
                s = -1.0;
            }

            typeMethods.r8vec_bracket3(nn, xn, aa, ref i);
            typeMethods.r8vec_bracket3(nn, xn, bb, ref j);
            //
            //  Evaluate the polynomial with the appropriate data.
            //
            if (i == j)
            {
                q[ii] = hermite_cubic_integrate(xn[i], fn[i], dn[i],
                    xn[i + 1], fn[i + 1], dn[i + 1], aa, bb);
            }
            else
            {
                q[ii] = hermite_cubic_integrate(xn[i], fn[i], dn[i],
                    xn[i + 1], fn[i + 1], dn[i + 1], aa, xn[i + 1]);

                int k;
                for (k = i + 1; k < j; k++)
                {
                    q[ii] += hermite_cubic_integral(xn[k], fn[k], dn[k],
                        xn[k + 1], fn[k + 1], dn[k + 1]);
                }

                q[ii] += hermite_cubic_integrate(xn[j], fn[j], dn[j],
                    xn[j + 1], fn[j + 1], dn[j + 1], xn[j], bb);
            }

            q[ii] = s * q[ii];
        }

        return q;
    }

    public static double[] hermite_cubic_spline_quad_rule(int nn, double[] xn)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HERMITE_CUBIC_SPLINE_QUAD_RULE: Hermite cubic spline quadrature rule.
        //
        //  Discussion:
        //
        //    The integral is taken over the definition interval ( X[0], X[NN-1] ).
        //
        //    Note that if the intervals are equal in size, then the derivative
        //    information in DN has no effect on the integral value,
        //    except for the first and last entries.
        //
        //    The quadrature rule is
        //
        //      Integral ( XN[0] <= X <= XN[NN-1] ) F(X) dX is approximately
        //
        //      Sum ( 0 <= I <= NN-1 ) W[0,I] * F(X[I]) + W[1,I] * F'(X[I])
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 March 2011
        //
        //  Author:
        //
        //    John Burkardt.
        //
        //  Reference:
        //
        //    Fred Fritsch, Ralph Carlson,
        //    Monotone Piecewise Cubic Interpolation,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 17, Number 2, April 1980, pages 238-246.
        //
        //  Parameters:
        //
        //    Input, int NN, the number of data points.
        //
        //    Input, double XN[NN], the coordinates of the data points.
        //    The entries in XN must be in strictly ascending order.
        //
        //    Output, double W[2*NN], the quadrature weights for F(1:NN)
        //    and DN(1:NN).
        //
    {
        int j;

        double[] w = new double[2 * nn];

        w[0 + 0 * 2] = 0.5 * (xn[1] - xn[0]);
        for (j = 1; j < nn - 1; j++)
        {
            w[0 + j * 2] = 0.5 * (xn[j + 1] - xn[j - 1]);
        }

        w[0 + (nn - 1) * 2] = 0.5 * (xn[nn - 1] - xn[nn - 2]);

        w[1 + 0 * 2] = Math.Pow(xn[1] - xn[0], 2) / 12.0;
        for (j = 1; j < nn - 1; j++)
        {
            w[1 + j * 2] = (xn[j + 1] - xn[j - 1])
                * (xn[j + 1] - 2.0 * xn[j] + xn[j - 1]) / 12.0;
        }

        w[1 + (nn - 1) * 2] = -Math.Pow(xn[nn - 2] - xn[nn - 1], 2) / 12.0;

        return w;
    }

    public static void hermite_cubic_spline_value(int nn, double[] xn, double[] fn,
            double[] dn, int n, double[] x, ref double[] f, ref double[] d, ref double[] s,
            ref double[] t )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HERMITE_CUBIC_SPLINE_VALUE evaluates a Hermite cubic spline.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 February 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Fred Fritsch, Ralph Carlson,
        //    Monotone Piecewise Cubic Interpolation,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 17, Number 2, April 1980, pages 238-246.
        //
        //  Parameters:
        //
        //    Input, int NN, the number of data points.
        //
        //    Input, double XN[NN], the coordinates of the data points.
        //    The entries in XN must be in strictly ascending order.
        //
        //    Input, double FN[NN], the function values.
        //
        //    Input, double DN[NN], the derivative values.
        //
        //    Input, int N, the number of sample points.
        //
        //    Input, double X[N], the coordinates of the sample points.
        //
        //    Output, double F[N], the function value at the sample points.
        //
        //    Output, double D[N], the derivative value at the sample points.
        //
        //    Output, double S[N], the second derivative value at the
        //    sample points.
        //
        //    Output, double T[N], the third derivative value at the
        //    sample points.
        //
    {
        int i;

        int left = n / 2;

        for (i = 0; i < n; i++)
        {
            typeMethods.r8vec_bracket3(nn, xn, x[i], ref left);

            hermite_cubic_value(xn[left], fn[left], dn[left], xn[left + 1],
                fn[left + 1], dn[left + 1], 1, x, ref f, ref d, ref s, ref t, xIndex: i, fIndex: i, dIndex: i, sIndex: i, tIndex: i);
        }
    }

    public static void hermite_cubic_to_power_cubic(double x1, double f1, double d1, double x2,
            double f2, double d2, ref double c0, ref double c1, ref double c2, ref double c3)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HERMITE_CUBIC_TO_POWER_CUBIC converts a Hermite cubic to power form.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 February 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Fred Fritsch, Ralph Carlson,
        //    Monotone Piecewise Cubic Interpolation,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 17, Number 2, April 1980, pages 238-246.
        //
        //  Parameters:
        //
        //    Input, double X1, F1, D1, the left endpoint, function value
        //    and derivative.
        //
        //    Input, double X2, F2, D2, the right endpoint, function value
        //    and derivative.
        //
        //    Output, double *C0, *C1, *C2, *C3, the power form of the polynomial.
        //
    {
        double h = x2 - x1;
        double df = (f2 - f1) / h;
        //
        //  Polynomial in terms of X - X1:
        //
        c0 = f1;
        c1 = d1;
        c2 = -(2.0 * d1 - 3.0 * df + d2) / h;
        c3 = (d1 - 2.0 * df + d2) / h / h;
        //
        //  Shift polynomial to X.
        //
        c2 -= x1 * c3;
        c1 -= x1 * c2;
        c0 -= x1 * c1;
        c2 -= x1 * c3;
        c1 -= x1 * c2;
        c2 -= x1 * c3;
    }

    public static void hermite_cubic_value(double x1, double f1, double d1, double x2,
            double f2, double d2, int n, double[] x, ref double[] f, ref double[] d,
            ref double[] s, ref double[] t, int xIndex = 0, int fIndex = 0, int dIndex = 0, int sIndex = 0, int tIndex = 0 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HERMITE_CUBIC_VALUE evaluates a Hermite cubic polynomial.
        //
        //  Discussion:
        //
        //    The input arguments can be vectors.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 February 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Fred Fritsch, Ralph Carlson,
        //    Monotone Piecewise Cubic Interpolation,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 17, Number 2, April 1980, pages 238-246.
        //
        //  Parameters:
        //
        //    Input, double X1, F1, D1, the left endpoint, function value
        //    and derivative.
        //
        //    Input, double X2, F2, D2, the right endpoint, function value
        //    and derivative.
        //
        //    Input, int N, the number of evaluation points.
        //
        //    Input, double X[N], the points at which the Hermite cubic
        //    is to be evaluated.
        //
        //    Output, double F[N], D[N], S[N], T[N], the value and first
        //    three derivatives of the Hermite cubic at X.
        //
    {
        int i;

        double h = x2 - x1;
        double df = (f2 - f1) / h;

        double c2 = -(2.0 * d1 - 3.0 * df + d2) / h;
        double c3 = (d1 - 2.0 * df + d2) / h / h;

        for (i = 0; i < n; i++)
        {
            f[fIndex + i] = f1 + (x[xIndex + i] - x1) * (d1
                                                         + (x[xIndex + i] - x1) * (c2
                                                             + (x[xIndex + i] - x1) * c3));
            d[dIndex + i] = d1 + (x[xIndex + i] - x1) * (2.0 * c2
                                                         + (x[xIndex + i] - x1) * 3.0 * c3);
            s[sIndex + i] = 2.0 * c2 + (x[xIndex + i] - x1) * 6.0 * c3;
            t[tIndex + i] = 6.0 * c3;
        }
    }

    public static void power_cubic_to_hermite_cubic(double c0, double c1, double c2, double c3,
            double x1, double x2, ref double f1, ref double d1, ref double f2, ref double d2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POWER_CUBIC_TO_HERMITE_CUBIC converts a power cubic to Hermite form.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 February 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Fred Fritsch, Ralph Carlson,
        //    Monotone Piecewise Cubic Interpolation,
        //    SIAM Journal on Numerical Analysis,
        //    Volume 17, Number 2, April 1980, pages 238-246.
        //
        //  Parameters:
        //
        //    Input, double C0, C1, C2, C3, the power form of the
        //    polynomial.
        //
        //    Input, double X1, X2, the left and right endpoints of
        //    the Hermite form.
        //
        //    Output, double *F1, *D1, the function and derivative values at X1.
        //
        //    Output, double *F2, *D2, the function and derivative values at X2.
        //
    {
        f1 = c0 + x1 * (c1 + x1 * (c2 + x1 * c3));
        d1 = c1 + x1 * (2.0 * c2 + x1 * 3.0 * c3);

        f2 = c0 + x2 * (c1 + x2 * (c2 + x2 * c3));
        d2 = c1 + x2 * (2.0 * c2 + x2 * 3.0 * c3);
    }
}