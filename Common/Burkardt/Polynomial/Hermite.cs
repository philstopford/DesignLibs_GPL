using System;
using Burkardt.MatrixNS;
using Burkardt.MonomialNS;
using Burkardt.Quadrature;

namespace Burkardt.PolynomialNS;

public static class Hermite
{
    public static void gen_hermite_poly ( int n, double x, double mu, ref double[] p )

        //******************************************************************************
        //
        //  Purpose:
        //
        //    GEN_HERMITE_POLY evaluates the generalized Hermite polynomials at X.
        //
        //  Discussion:
        //
        //    The generalized Hermite polynomials are orthogonal under the weight
        //    function:
        //
        //      w(x) = |x|^(2*MU) * exp ( - x^2 )
        //
        //    over the interval (-oo,+oo).
        //
        //    When MU = 0, the generalized Hermite polynomial reduces to the standard
        //    Hermite polynomial.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    26 February 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Theodore Chihara,
        //    An Introduction to Orthogonal Polynomials,
        //    Gordon and Breach, 1978,
        //    ISBN: 0677041500,
        //    LC: QA404.5 C44.
        //
        //  Parameters:
        //
        //    Input, int N, the highest order polynomial to compute.
        //
        //    Input, double X, the point at which the polynomials are 
        //    to be evaluated.
        //
        //    Input, double MU, the parameter.
        //    - 1 / 2 < MU.
        //
        //    Output, double P[N+1], the values of the first N+1
        //    polynomials at the point X.
        //
    {
        int i;

        switch (n)
        {
            case < 0:
                return;
        }

        p[0] = 1.0;

        switch (n)
        {
            case 0:
                return;
        }

        p[1] = 2.0 * x;
 
        for ( i = 1; i < n; i++ )
        {
            double theta = (i % 2) switch
            {
                0 => 0.0,
                _ => 2.0 * mu
            };

            p[i+1] = 2.0 * x * p[i] - 2.0 * ( i + theta ) * p[i-1];
        }
    }

    public static void hermite_poly_phys(int n, double x, ref double[] cx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HERMITE_POLY_PHYS evaluates the physicist's Hermite polynomials at X.
        //
        //  Differential equation:
        //
        //    Y'' - 2 X Y' + 2 N Y = 0
        //
        //  First terms:
        //
        //      1
        //      2 X
        //      4 X^2     -  2
        //      8 X^3     - 12 X
        //     16 X^4     - 48 X^2     + 12
        //     32 X^5    - 160 X^3    + 120 X
        //     64 X^6    - 480 X^4    + 720 X^2    - 120
        //    128 X^7   - 1344 X^5   + 3360 X^3   - 1680 X
        //    256 X^8   - 3584 X^6  + 13440 X^4  - 13440 X^2   + 1680
        //    512 X^9   - 9216 X^7  + 48384 X^5  - 80640 X^3  + 30240 X
        //   1024 X^10 - 23040 X^8 + 161280 X^6 - 403200 X^4 + 302400 X^2 - 30240
        //
        //  Recursion:
        //
        //    H(0,X) = 1,
        //    H(1,X) = 2*X,
        //    H(N,X) = 2*X * H(N-1,X) - 2*(N-1) * H(N-2,X)
        //
        //  Norm:
        //
        //    Integral ( -oo < X < +oo ) exp ( - X^2 ) * H(N,X)^2 dX
        //    = sqrt ( PI ) * 2^N * N!
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 May 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Milton Abramowitz, Irene Stegun,
        //    Handbook of Mathematical Functions,
        //    National Bureau of Standards, 1964,
        //    ISBN: 0-486-61272-4,
        //    LC: QA47.A34.
        //
        //  Parameters:
        //
        //    Input, int N, the highest order polynomial to compute.
        //    Note that polynomials 0 through N will be computed.
        //
        //    Input, double X, the point at which the polynomials are to be evaluated.
        //
        //    Output, double CX[N+1], the values of the first N+1 Hermite
        //    polynomials at the point X.
        //
    {
        int i;

        switch (n)
        {
            case < 0:
                return;
        }

        cx[0] = 1.0;

        switch (n)
        {
            case 0:
                return;
        }

        cx[1] = 2.0 * x;

        for (i = 2; i <= n; i++)
        {
            cx[i] = 2.0 * x * cx[i - 1] - 2.0 * (i - 1) * cx[i - 2];
        }
    }

    public static void hermite_poly_phys_coef(int n, ref double[] c)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HERMITE_POLY_PHYS_COEF: coefficients of the physicist's Hermite polynomial H(n,x).
        //
        //  First terms:
        //
        //    N/K     0     1      2      3       4     5      6    7      8    9   10
        //
        //     0      1
        //     1      0     2
        //     2     -2     0      4
        //     3      0   -12      0      8
        //     4     12     0    -48      0      16
        //     5      0   120      0   -160       0    32
        //     6   -120     0    720      0    -480     0     64
        //     7      0 -1680      0   3360       0 -1344      0   128
        //     8   1680     0 -13440      0   13440     0  -3584     0    256
        //     9      0 30240      0 -80640       0 48384      0 -9216      0 512
        //    10 -30240     0 302400      0 -403200     0 161280     0 -23040   0 1024
        //
        //  Recursion:
        //
        //    H(0,X) = 1,
        //    H(1,X) = 2*X,
        //    H(N,X) = 2*X * H(N-1,X) - 2*(N-1) * H(N-2,X)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 May 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Milton Abramowitz, Irene Stegun,
        //    Handbook of Mathematical Functions,
        //    National Bureau of Standards, 1964,
        //    ISBN: 0-486-61272-4,
        //    LC: QA47.A34.
        //
        //  Parameters:
        //
        //    Input, int N, the highest order polynomial to compute.
        //    Note that polynomials 0 through N will be computed.
        //
        //    Output, double C[(N+1)*(N+1)], the coefficients of the Hermite
        //    polynomials.
        //
    {
        int i;
        int j;

        switch (n)
        {
            case < 0:
                return;
        }

        for (i = 0; i <= n; i++)
        {
            for (j = 0; j <= n; j++)
            {
                c[i + j * (n + 1)] = 0.0;
            }
        }

        c[0 + 0 * (n + 1)] = 1.0;

        switch (n)
        {
            case 0:
                return;
        }

        c[1 + 1 * (n + 1)] = 2.0;

        for (i = 2; i <= n; i++)
        {
            c[i + 0 * (n + 1)] = -2.0 * (i - 1) * c[i - 2 + 0 * (n + 1)];
            for (j = 1; j <= i - 2; j++)
            {
                c[i + j * (n + 1)] = 2.0 * c[i - 1 + (j - 1) * (n + 1)]
                                     - 2.0 * (i - 1) * c[i - 2 + j * (n + 1)];
            }

            c[i + (i - 1) * (n + 1)] = 2.0 * c[i - 1 + (i - 2) * (n + 1)];
            c[i + i * (n + 1)] = 2.0 * c[i - 1 + (i - 1) * (n + 1)];
        }
    }

    public static double[] h_polynomial_coefficients(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    H_POLYNOMIAL_COEFFICIENTS: coefficients of H(i,x).
        //
        //  Discussion:
        //
        //    H(i,x) is the physicist's Hermite polynomial of degree I.
        //
        //  First terms:
        //
        //    N/K     0     1      2      3       4     5      6    7      8    9   10
        //
        //     0      1
        //     1      0     2
        //     2     -2     0      4
        //     3      0   -12      0      8
        //     4     12     0    -48      0      16
        //     5      0   120      0   -160       0    32
        //     6   -120     0    720      0    -480     0     64
        //     7      0 -1680      0   3360       0 -1344      0   128
        //     8   1680     0 -13440      0   13440     0  -3584     0    256
        //     9      0 30240      0 -80640       0 48384      0 -9216      0 512
        //    10 -30240     0 302400      0 -403200     0 161280     0 -23040   0 1024
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 March 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Milton Abramowitz, Irene Stegun,
        //    Handbook of Mathematical Functions,
        //    National Bureau of Standards, 1964,
        //    ISBN: 0-486-61272-4,
        //    LC: QA47.A34.
        //
        //  Parameters:
        //
        //    Input, int N, the highest order polynomial to compute.
        //    Note that polynomials 0 through N will be computed.
        //
        //    Output, double HERMITE_POLYNOMIAL_COEFFICIENTS[(N+1)*(N+1)], the 
        //    coefficients of the Hermite polynomials of orders 0 through N.
        //
    {
        int i;
        int j;

        switch (n)
        {
            case < 0:
                return null;
        }

        double[] c = new double[(n + 1) * (n + 1)];

        for (i = 0; i <= n; i++)
        {
            for (j = 0; j <= n; j++)
            {
                c[i + j * (n + 1)] = 0.0;
            }
        }

        c[0 + 0 * (n + 1)] = 1.0;

        switch (n)
        {
            case 0:
                return c;
        }

        c[1 + 1 * (n + 1)] = 2.0;

        for (i = 2; i <= n; i++)
        {
            c[i + 0 * (n + 1)] = -2.0 * (i - 1) * c[i - 2 + 0 * (n + 1)];
            for (j = 1; j <= i - 2; j++)
            {
                c[i + j * (n + 1)] = 2.0 * c[i - 1 + (j - 1) * (n + 1)]
                                     - 2.0 * (i - 1) * c[i - 2 + j * (n + 1)];
            }

            c[i + (i - 1) * (n + 1)] = 2.0 * c[i - 1 + (i - 2) * (n + 1)];
            c[i + i * (n + 1)] = 2.0 * c[i - 1 + (i - 1) * (n + 1)];
        }

        return c;
    }

    public static double[] h_polynomial_value(int m, int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    H_POLYNOMIAL_VALUE evaluates H(i,x).
        //
        //  Discussion:
        //
        //    H(i,x) is the physicist's Hermite polynomial of degree I.
        //
        //  Differential equation:
        //
        //    Y'' - 2 X Y' + 2 N Y = 0
        //
        //  First terms:
        //
        //      1
        //      2 X
        //      4 X^2     -  2
        //      8 X^3     - 12 X
        //     16 X^4     - 48 X^2     + 12
        //     32 X^5    - 160 X^3    + 120 X
        //     64 X^6    - 480 X^4    + 720 X^2    - 120
        //    128 X^7   - 1344 X^5   + 3360 X^3   - 1680 X
        //    256 X^8   - 3584 X^6  + 13440 X^4  - 13440 X^2   + 1680
        //    512 X^9   - 9216 X^7  + 48384 X^5  - 80640 X^3  + 30240 X
        //   1024 X^10 - 23040 X^8 + 161280 X^6 - 403200 X^4 + 302400 X^2 - 30240
        //
        //  Recursion:
        //
        //    H(0,X) = 1,
        //    H(1,X) = 2*X,
        //    H(N,X) = 2*X * H(N-1,X) - 2*(N-1) * H(N-2,X)
        //
        //  Norm:
        //
        //    Integral ( -oo < X < +oo ) exp ( - X^2 ) * H(N,X)^2 dX
        //    = sqrt ( PI ) * 2^N * N!
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 May 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Milton Abramowitz, Irene Stegun,
        //    Handbook of Mathematical Functions,
        //    National Bureau of Standards, 1964,
        //    ISBN: 0-486-61272-4,
        //    LC: QA47.A34.
        //
        //  Parameters:
        //
        //    Input, int M, the number of evaluation points.
        //
        //    Input, int N, the highest order polynomial to compute.
        //    Note that polynomials 0 through N will be computed.
        //
        //    Input, double X[M], the evaluation points.
        //
        //    Output, double H_POLYNOMIAL_VALUE[M*(N+1)], the values of the first 
        //    N+1 Hermite polynomials at the evaluation points.
        //
    {
        int i;
        int j;

        switch (n)
        {
            case < 0:
                return null;
        }

        double[] p = new double[m * (n + 1)];

        for (i = 0; i < m; i++)
        {
            p[i + 0 * m] = 1.0;
        }

        switch (n)
        {
            case 0:
                return p;
        }

        for (i = 0; i < m; i++)
        {
            p[i + 1 * m] = 2.0 * x[i];
        }

        for (j = 2; j <= n; j++)
        {
            for (i = 0; i < m; i++)
            {
                p[i + j * m] = 2.0 * x[i] * p[i + (j - 1) * m]
                               - 2.0 * (j - 1) * p[i + (j - 2) * m];
            }
        }

        return p;
    }


    public static double[] h_polynomial_zeros(int nt)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    H_POLYNOMIAL_ZEROS: zeros of H(i,x).
        //
        //  Discussion:
        //
        //    H(i,x) is the physicist's Hermite polynomial of degree I.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 February 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NT, the degree of the polynomial.
        //
        //    Output, double H_POLYNOMIAL_ZEROS[NT], the zeros of the polynomial.
        //
    {
        int i;

        double[] z = new double[nt];

        for (i = 0; i < nt; i++)
        {
            z[i] = 0.0;
        }

        double[] bj = new double[nt];

        for (i = 0; i < nt; i++)
        {
            bj[i] = Math.Sqrt((i + 1) / 2.0);
        }

        double[] wts = new double[nt];
        for (i = 0; i < nt; i++)
        {
            wts[i] = 0.0;
        }

        wts[0] = Math.Sqrt(Math.Sqrt(Math.PI));

        IMTQLX.imtqlx(nt, ref z, ref bj, ref wts);

        return z;
    }

    public static void hermite_recur(ref double p2, ref double dp2, ref double p1, double x, int order)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HERMITE_RECUR finds the value and derivative of a Hermite polynomial.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 April 2006
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Arthur Stroud, Don Secrest,
        //    Gaussian Quadrature Formulas,
        //    Prentice Hall, 1966,
        //    LC: QA299.4G3S7.
        //
        //  Parameters:
        //
        //    Output, double *P2, the value of H(ORDER)(X).
        //
        //    Output, double *DP2, the value of H'(ORDER)(X).
        //
        //    Output, double *P1, the value of H(ORDER-1)(X).
        //
        //    Input, double X, the point at which polynomials are evaluated.
        //
        //    Input, int ORDER, the order of the polynomial to be computed.
        //
    {
        int i;

        double q1 = 1.0;
        double dq1 = 0.0;

        double q2 = x;
        double dq2 = 1.0;

        for (i = 2; i <= order; i++)
        {
            double q0 = q1;
            double dq0 = dq1;

            q1 = q2;
            dq1 = dq2;

            q2 = x * q1 - 0.5 * (i - 1.0) * q0;
            dq2 = x * dq1 + q1 - 0.5 * (i - 1.0) * dq0;
        }

        p2 = q2;
        dp2 = dq2;
        p1 = q1;
    }

    public static void hermite_root(ref double x, int order, ref double dp2, ref double p1)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HERMITE_ROOT improves an approximate root of a Hermite polynomial.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 May 2006
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Arthur Stroud, Don Secrest.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Arthur Stroud, Don Secrest,
        //    Gaussian Quadrature Formulas,
        //    Prentice Hall, 1966,
        //    LC: QA299.4G3S7.
        //
        //  Parameters:
        //
        //    Input/output, double *X, the approximate root, which
        //    should be improved on output.
        //
        //    Input, int ORDER, the order of the Hermite polynomial.
        //
        //    Output, double *DP2, the value of H'(ORDER)(X).
        //
        //    Output, double *P1, the value of H(ORDER-1)(X).
        //
    {
        const double eps = 1.0E-12;
        double p2 = 0;
        int step;
        const int step_max = 10;

        for (step = 1; step <= step_max; step++)
        {
            hermite_recur(ref p2, ref dp2, ref p1, x, order);

            double d = p2 / dp2;
            x -= d;

            if (Math.Abs(d) <= eps * (Math.Abs(x) + 1.0))
            {
                return;
            }
        }
    }

    public static void hep_coefficients(int n, ref int o, ref double[] c, ref int[] f)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HEP_COEFFICIENTS: coefficients of Hermite polynomials He(n,x).
        //
        //  Discussion:
        //
        //    He(i,x) represents the probabilist's Hermite polynomial.
        //
        //  First terms:
        //
        //    N/K     0     1      2      3       4     5      6    7      8    9   10
        //
        //     0      1
        //     1      0     1
        //     2     -1     0      1
        //     3      0    -3      0      1
        //     4      3     0     -6      0       1
        //     5      0    15      0    -10       0     1
        //     6    -15     0     45      0     -15     0      1
        //     7      0  -105      0    105       0   -21      0     1
        //     8    105     0   -420      0     210     0    -28     0      1
        //     9      0   945      0  -1260       0   378      0   -36      0   1
        //    10   -945     0   4725      0   -3150     0    630     0    -45   0    1
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 October 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Milton Abramowitz, Irene Stegun,
        //    Handbook of Mathematical Functions,
        //    National Bureau of Standards, 1964,
        //    ISBN: 0-486-61272-4,
        //    LC: QA47.A34.
        //
        //    Daniel Zwillinger, editor,
        //    CRC Standard Mathematical Tables and Formulae,
        //    30th Edition,
        //    CRC Press, 1996.
        //
        //  Parameters:
        //
        //    Input, int N, the highest order polynomial to evaluate.
        //    Note that polynomials 0 through N will be evaluated.
        //
        //    Output, int &O, the number of coefficients.
        //
        //    Output, double C[(N+2)/2], the coefficients of the Legendre
        //    polynomial of degree N.
        //
        //    Output, int F[(N+2)/2], the exponents.
        //
    {
        int i;
        int j;

        double[] ct = new double[(n + 1) * (n + 1)];

        for (i = 0; i <= n; i++)
        {
            for (j = 0; j <= n; j++)
            {
                ct[i + j * (n + 1)] = 0.0;
            }
        }

        ct[0 + 0 * (n + 1)] = 1.0;

        switch (n)
        {
            case > 0:
            {
                ct[1 + 1 * (n + 1)] = 1.0;

                for (i = 2; i <= n; i++)
                {
                    ct[i + 0 * (n + 1)] = -(double) (i - 1) * ct[i - 2 + 0 * (n + 1)];
                    for (j = 1; j <= i - 2; j++)
                    {
                        ct[i + j * (n + 1)] =
                            ct[i - 1 + (j - 1) * (n + 1)] - (i - 1) * ct[i - 2 + j * (n + 1)];
                    }

                    ct[i + (i - 1) * (n + 1)] = ct[i - 1 + (i - 2) * (n + 1)];
                    ct[i + i * (n + 1)] = ct[i - 1 + (i - 1) * (n + 1)];
                }

                break;
            }
        }

        //
        //  Extract the nonzero data from the alternating columns of the last row.
        //
        o = (n + 2) / 2;

        int k = o;
        for (j = n; 0 <= j; j -= 2)
        {
            k -= 1;
            c[k] = ct[n + j * (n + 1)];
            f[k] = j;
        }
    }

    public static double[] hep_value(int n, int o, ref double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HEP_VALUE evaluates the Hermite polynomials He(i,x).
        //
        //  Discussion:
        //
        //    He(i,x) represents the probabilist's Hermite polynomial.
        //
        //    1
        //    X
        //    X^2  -  1
        //    X^3  -  3 X
        //    X^4  -  6 X^2 +   3
        //    X^5  - 10 X^3 +  15 X
        //    X^6  - 15 X^4 +  45 X^2 -   15
        //    X^7  - 21 X^5 + 105 X^3 -  105 X
        //    X^8  - 28 X^6 + 210 X^4 -  420 X^2 +  105
        //    X^9  - 36 X^7 + 378 X^5 - 1260 X^3 +  945 X
        //    X^10 - 45 X^8 + 630 X^6 - 3150 X^4 + 4725 X^2 - 945
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 October 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Milton Abramowitz, Irene Stegun,
        //    Handbook of Mathematical Functions,
        //    National Bureau of Standards, 1964,
        //    ISBN: 0-486-61272-4,
        //    LC: QA47.A34.
        //
        //    Daniel Zwillinger, editor,
        //    CRC Standard Mathematical Tables and Formulae,
        //    30th Edition,
        //    CRC Press, 1996.
        //
        //  Parameters:
        //
        //    Input, int N, the number of evaluation points.
        //
        //    Input, int O, the degree of the polynomial.
        //
        //    Input, double X[N], the evaluation points.
        //
        //    Output, double LP_VALUE[N], the value of the Hermite polynomial 
        //    of degree N at the points X.
        //
    {
        int i;

        double[] vtable = new double[n * (o + 1)];

        for (i = 0; i < n; i++)
        {
            vtable[i + 0 * n] = 1.0;
        }

        switch (o)
        {
            case >= 1:
            {
                for (i = 0; i < n; i++)
                {
                    vtable[i + 1 * n] = x[i];
                }

                int j;
                for (j = 2; j <= o; j++)
                {
                    for (i = 0; i < n; i++)
                    {
                        vtable[i + j * n] = x[i] * vtable[i + (j - 1) * n]
                                            - (j - 1) * vtable[i + (j - 2) * n];
                    }
                }

                break;
            }
        }

        double[] v = new double[n];

        for (i = 0; i < n; i++)
        {
            v[i] = vtable[i + o * n];
        }

        return v;
    }


    public static void hepp_to_polynomial(int m, int[] l, int o_max, int o, ref double[] c,
            ref int[] e)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HEPP_TO_POLYNOMIAL writes a Hermite Product Polynomial as a polynomial.
        //
        //  Discussion:
        //
        //    For example, if 
        //      M = 3,
        //      L = ( 1, 0, 2 ),
        //    then
        //      He(1,0,2)(X,Y,Z) 
        //      = He(1)(X) * He(0)(Y) * He(2)(Z)
        //      = X * 1 * ( Z^3-3Z)
        //      = - 3XZ + X Z^3
        //    so
        //      O = 2 (2 nonzero terms)
        //      C = -3.0
        //           1.0
        //      E =  8   <-- index in 3-space of exponent (1,0,1)
        //          23   <-- index in 3-space of exponent (1,0,3)
        //
        //    The output value of O is no greater than
        //      O_MAX = product ( 1 <= I <= M ) (L(I)+2)/2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 October 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the spatial dimension.
        //
        //    Input, int L[M], the index of each product 
        //    polynomial factor.  0 <= L(*).
        //
        //    Input, int O_MAX, an upper limit on the size of the 
        //    output arrays.
        //      O_MAX = product ( 1 <= I <= M ) (L(I)+2)/2.
        //
        //    Output, int &O, the "order" of the polynomial product.
        //
        //    Output, double C[O], the coefficients of the polynomial product.
        //
        //    Output, int E[O], the indices of the exponents of the 
        //    polynomial product.
        //
    {
        int i;
        int o2 = 0;
        int[] p = new int[1];

        double[] c1 = new double[o_max];
        double[] c2 = new double[o_max];
        int[] e1 = new int[o_max];
        int[] f2 = new int[o_max];
        int[] pp = new int[m];

        int o1 = 1;
        c1[0] = 1.0;
        e1[0] = 1;
        //
        //  Implicate one factor at a time.
        //
        for (i = 0; i < m; i++)
        {
            hep_coefficients(l[i], ref o2, ref c2, ref f2);

            o = 0;

            int j2;
            for (j2 = 0; j2 < o2; j2++)
            {
                int j1;
                for (j1 = 0; j1 < o1; j1++)
                {
                    c[o] = c1[j1] * c2[j2];
                    p = i switch
                    {
                        > 0 => Monomial.mono_unrank_grlex(i, e1[j1]),
                        _ => p
                    };

                    int i2;
                    for (i2 = 0; i2 < i; i2++)
                    {
                        if (p != null)
                        {
                            pp[i2] = p[i2];
                        }
                    }

                    pp[i] = f2[j2];
                    e[o] = Monomial.mono_rank_grlex(i + 1, pp);
                    o += 1;
                    p = i switch
                    {
                        > 0 => null,
                        _ => p
                    };
                }
            }

            Polynomial.polynomial_sort(o, ref c, ref e);
            Polynomial.polynomial_compress(o, c, e, ref o, ref c, ref e);

            o1 = o;
            int i1;
            for (i1 = 0; i1 < o; i1++)
            {
                c1[i1] = c[i1];
                e1[i1] = e[i1];
            }
        }
    }

    public static double[] hepp_value(int m, int n, int[] o, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HEPP_VALUE evaluates a Hermite Product Polynomial.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 October 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the spatial dimension.
        //
        //    Input, int N, the number of evaluation points.
        //
        //    Input, int O[M], the degree of the polynomial factors.
        //    0 <= O(*).
        //
        //    Input, double X[M*N], the evaluation points.
        //
        //    Output, double LPP_VALUE[N], the value of the Hermite Product 
        //    Polynomial of degree O at the points X.
        //
    {
        int i;
        int j;

        double[] v = new double[n];

        for (j = 0; j < n; j++)
        {
            v[j] = 1.0;
        }

        double[] xi = new double[n];

        for (i = 0; i < m; i++)
        {
            for (j = 0; j < n; j++)
            {
                xi[j] = x[i + j * m];
            }

            double[] vi = hep_value(n, o[i], ref xi);
            for (j = 0; j < n; j++)
            {
                v[j] *= vi[j];
            }
        }

        return v;
    }

    public static double[] he_polynomial_coefficients(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HE_POLYNOMIAL_COEFFICIENTS: coefficients of He(i,x).
        //
        //  Discussion:
        //
        //    He(i,x) represents the probabilist's Hermite polynomial.
        //
        //  First terms:
        //
        //    N/K     0     1      2      3       4     5      6    7      8    9   10
        //
        //     0      1
        //     1      0     1
        //     2     -1     0      1
        //     3      0    -3      0      1
        //     4      3     0     -6      0       1
        //     5      0    15      0    -10       0     1
        //     6    -15     0     45      0     -15     0      1
        //     7      0  -105      0    105       0   -21      0     1
        //     8    105     0   -420      0     210     0    -28     0      1
        //     9      0   945      0  -1260       0   378      0   -36      0   1
        //    10   -945     0   4725      0   -3150     0    630     0    -45   0    1
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 March 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Milton Abramowitz, Irene Stegun,
        //    Handbook of Mathematical Functions,
        //    National Bureau of Standards, 1964,
        //    ISBN: 0-486-61272-4,
        //    LC: QA47.A34.
        //
        //  Parameters:
        //
        //    Input, int N, the highest order polynomial to compute.
        //    Note that polynomials 0 through N will be computed.
        //
        //    Output, double HE_POLYNOMIAL_COEFFICIENTS[(N+1)*(N+1)], the coefficients 
        //    of the Hermite polynomials of orders 0 through N.
        //
    {
        int i;
        int j;

        switch (n)
        {
            case < 0:
                return null;
        }

        double[] c = new double[(n + 1) * (n + 1)];

        for (i = 0; i <= n; i++)
        {
            for (j = 0; j <= n; j++)
            {
                c[i + j * (n + 1)] = 0.0;
            }
        }

        c[0 + 0 * (n + 1)] = 1.0;

        switch (n)
        {
            case 0:
                return c;
        }

        c[1 + 1 * (n + 1)] = 1.0;

        for (i = 2; i <= n; i++)
        {
            c[i + 0 * (n + 1)] = -(double) (i - 1) * c[i - 2 + 0 * (n + 1)];
            for (j = 1; j <= i - 2; j++)
            {
                c[i + j * (n + 1)] = c[i - 1 + (j - 1) * (n + 1)] - (i - 1) * c[i - 2 + j * (n + 1)];
            }

            c[i + (i - 1) * (n + 1)] = c[i - 1 + (i - 2) * (n + 1)];
            c[i + i * (n + 1)] = c[i - 1 + (i - 1) * (n + 1)];
        }

        return c;
    }

    public static double[] he_polynomial_value(int m, int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HE_POLYNOMIAL_VALUE evaluates He(i,x).
        //
        //  Discussion:
        //
        //    He(i,x) represents the probabilist's Hermite polynomial.
        //
        //  Differential equation:
        //
        //    ( exp ( - 0.5 * x^2 ) * He(n,x)' )' + n * exp ( - 0.5 * x^2 ) * He(n,x) = 0
        //
        //  First terms:
        //
        //   1
        //   X
        //   X^2  -  1
        //   X^3  -  3 X
        //   X^4  -  6 X^2 +   3
        //   X^5  - 10 X^3 +  15 X
        //   X^6  - 15 X^4 +  45 X^2 -   15
        //   X^7  - 21 X^5 + 105 X^3 -  105 X
        //   X^8  - 28 X^6 + 210 X^4 -  420 X^2 +  105
        //   X^9  - 36 X^7 + 378 X^5 - 1260 X^3 +  945 X
        //   X^10 - 45 X^8 + 630 X^6 - 3150 X^4 + 4725 X^2 - 945
        //
        //  Recursion:
        //
        //    He(0,X) = 1,
        //    He(1,X) = X,
        //    He(N,X) = X * He(N-1,X) - (N-1) * He(N-2,X)
        //
        //  Orthogonality:
        //
        //    Integral ( -oo < X < +oo ) exp ( - 0.5 * X^2 ) * He(M,X) He(N,X) dX 
        //    = sqrt ( 2 * Math.PI ) * N// * delta ( N, M )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 February 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Milton Abramowitz, Irene Stegun,
        //    Handbook of Mathematical Functions,
        //    National Bureau of Standards, 1964,
        //    ISBN: 0-486-61272-4,
        //    LC: QA47.A34.
        //
        //    Frank Olver, Daniel Lozier, Ronald Boisvert, Charles Clark,
        //    NIST Handbook of Mathematical Functions,
        //    Cambridge University Press, 2010,
        //    ISBN: 978-0521192255,
        //    LC: QA331.N57.
        //
        //  Parameters:
        //
        //    Input, int M, the number of evaluation points.
        //
        //    Input, int N, the highest order polynomial to compute.
        //    Note that polynomials 0 through N will be computed.
        //
        //    Input, double X[M], the evaluation points.
        //
        //    Output, double HE_POLYNOMIAL_VALUE[M*(N+1)], the values of the
        //    probabilist's Hermite polynomials of index 0 through N.
        //
    {
        int i;
        int j;

        switch (n)
        {
            case < 0:
                return null;
        }

        double[] p = new double[m * (n + 1)];

        for (i = 0; i < m; i++)
        {
            p[i + 0 * m] = 1.0;
        }

        switch (n)
        {
            case 0:
                return p;
        }

        for (i = 0; i < m; i++)
        {
            p[i + 1 * m] = x[i];
        }

        for (j = 2; j <= n; j++)
        {
            for (i = 0; i < m; i++)
            {
                p[i + j * m] = x[i] * p[i + (j - 1) * m] - (j - 1) * p[i + (j - 2) * m];
            }
        }

        return p;
    }


    public static double[] he_polynomial_zeros(int nt)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HE_POLYNOMIAL_ZEROS: zeros of He(i,x).
        //
        //  Discussion:
        //
        //    He(i,x) represents the probabilist's Hermite polynomial.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 February 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NT, the degree of the polynomial.
        //
        //    Output, double HE_POLYNOMIAL_ZEROS[NT], the zeros of the polynomial.
        //
    {
        int i;

        double[] z = new double[nt];

        for (i = 0; i < nt; i++)
        {
            z[i] = 0.0;
        }

        double[] bj = new double[nt];

        for (i = 0; i < nt; i++)
        {
            bj[i] = Math.Sqrt((i + 1) / 2.0);
        }

        double[] wts = new double[nt];

        for (i = 0; i < nt; i++)
        {
            wts[i] = 0.0;
        }

        wts[0] = Math.Sqrt(Math.Sqrt(Math.PI));

        IMTQLX.imtqlx(nt, ref z, ref bj, ref wts);

        for (i = 0; i < nt; i++)
        {
            z[i] *= Math.Sqrt(2.0);
        }

        return z;
    }

    public static double[] hen_exponential_product(int p, double b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HEN_EXPONENTIAL_PRODUCT: exponential product exp(b*x)*Hen(i,x)*Hen(j,x).
        //
        //  Discussion:
        //
        //    Hen(i,x) is the normalized probabilist's Hermite polynomial of degree I.
        //
        //
        //    For polynomial chaos applications, it is of interest to know the
        //    value of the integrals of products of exp(B*X) with every possible pair
        //    of basis functions.  That is, we'd like to form
        //
        //      Tij = Integral ( -oo < X < +oo ) 
        //        exp(B*X) * Hen(I,X) * Hen(J,X) exp(-0.5*X*X) dx
        //
        //    We will estimate these integrals using Gauss-Hermite quadrature.
        //    Because of the exponential factor exp(B*X), the quadrature will not 
        //    be exact.
        //
        //    However, when B = 0, the quadrature is exact, and moreoever, the
        //    table will be the identity matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 February 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int P, the maximum degree of the 
        //    polyonomial factors.  0 <= P.
        //
        //    Input, double B, the coefficient of X in the exponential factor.
        //
        //    Output, double HEN_EXPONENTIAL_PRODUCT[(P+1)*(P+1)], the table of 
        //    integrals.  TABLE(I,J) represents the weighted integral of 
        //    exp(B*X) * Hen(I,X) * Hen(J,X).
        //
    {
        int i;
        int j;
        int k;

        double[] table = new double[(p + 1) * (p + 1)];
        for (j = 0; j <= p; j++)
        {
            for (i = 0; i <= p; i++)
            {
                table[i + j * (p + 1)] = 0.0;
            }
        }

        int order = (3 * p + 4) / 2;

        double[] x_table = new double[order];
        double[] w_table = new double[order];

        HermiteQuadrature.he_quadrature_rule(order, ref x_table, ref w_table);

        for (k = 0; k < order; k++)
        {
            double x = x_table[k];
            double[] h_table = hen_polynomial_value(1, p, x_table, +k);
            //
            //  The following formula is an outer product in H_TABLE.
            //
            for (j = 0; j <= p; j++)
            {
                for (i = 0; i <= p; i++)
                {
                    table[i + j * (p + 1)] += w_table[k] * Math.Exp(b * x) * h_table[i] * h_table[j];
                }
            }
        }

        return table;
    }

    public static double[] hen_polynomial_value(int m, int n, double[] x, int xIndex = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HEN_POLYNOMIAL_VALUE: evaluates Hen(i,x).
        //
        //  Discussion:
        //
        //    Hen(i,x) is the normalized probabilist's Hermite polynomial of degree I.
        //
        //    These polynomials satisfy the orthonormality condition:
        //
        //      Integral ( -oo < X < +oo ) exp ( - 0.5 * X^2 ) * Hen(M,X) Hen(N,X) dX 
        //      = delta ( N, M )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 February 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Milton Abramowitz, Irene Stegun,
        //    Handbook of Mathematical Functions,
        //    National Bureau of Standards, 1964,
        //    ISBN: 0-486-61272-4,
        //    LC: QA47.A34.
        //
        //    Frank Olver, Daniel Lozier, Ronald Boisvert, Charles Clark,
        //    NIST Handbook of Mathematical Functions,
        //    Cambridge University Press, 2010,
        //    ISBN: 978-0521192255,
        //    LC: QA331.N57.
        //
        //  Parameters:
        //
        //    Input, int M, the number of evaluation points.
        //
        //    Input, int N, the highest order polynomial to compute.
        //    Note that polynomials 0 through N will be computed.
        //
        //    Input, double X[M], the evaluation points.
        //
        //    Output, double HEN_POLYNOMIAL_VALUE[M*(N+1)], the values of the 
        //    polynomials of index 0 through N.
        //
    {
        int i;
        int j;


        switch (n)
        {
            case < 0:
                return null;
        }

        double[] p = new double[m * (n + 1)];

        for (i = 0; i < m; i++)
        {
            p[i + 0 * m] = 1.0;
        }

        switch (n)
        {
            case 0:
                return p;
        }

        for (i = 0; i < m; i++)
        {
            p[i + 1 * m] = x[xIndex + i];
        }

        for (j = 2; j <= n; j++)
        {
            for (i = 0; i < m; i++)
            {
                p[i + j * m] = x[xIndex + i] * p[i + (j - 1) * m] - (j - 1) * p[i + (j - 2) * m];
            }
        }

        //
        //  Normalize.
        //
        double fact = 1.0;
        for (j = 0; j <= n; j++)
        {
            for (i = 0; i < m; i++)
            {
                p[i + j * m] /= Math.Sqrt(fact * Math.Sqrt(2.0 * Math.PI));
            }

            fact *= j + 1;
        }

        return p;
    }

    public static double[] hen_power_product(int p, int e)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HEN_POWER_PRODUCT: power products, x^e*Hen(i,x)*Hen(j,x).
        //
        //  Discussion:
        //
        //    Hen(i,x) is the normalized probabilist's Hermite polynomial of degree I.
        //
        //    For polynomial chaos applications, it is of interest to know the
        //    value of the integrals of products of X with every possible pair
        //    of basis functions.  That is, we'd like to form
        //
        //      Tij = Integral ( -oo < X < +oo ) 
        //        X^E * Hen(I,X) * Hen(J,X) exp(-0.5*X*X) dx
        //
        //    We will estimate these integrals using Gauss-Hermite quadrature.
        //
        //    When E is 0, the computed table should be the identity matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 February 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int P, the maximum degree of the polyonomial 
        //    factors.  0 <= P.
        //
        //    Input, int E, the exponent of X in the integrand.
        //    0 <= E.
        //
        //    Output, double HEN_POWER_PRODUCT[(P+1)*(P+1)], the table of integrals.  
        //    TABLE(I,J) represents the weighted integral of 
        //    X^E * Hen(I,X) * Hen(J,X).
        //
    {
        int i;
        int j;
        int k;

        double[] table = new double[(p + 1) * (p + 1)];
        for (j = 0; j <= p; j++)
        {
            for (i = 0; i <= p; i++)
            {
                table[i + j * (p + 1)] = 0.0;
            }
        }

        int order = p + 1 + (e + 1) / 2;

        double[] x_table = new double[order];
        double[] w_table = new double[order];

        HermiteQuadrature.he_quadrature_rule(order, ref x_table, ref w_table);

        for (k = 0; k < order; k++)
        {
            double x = x_table[k];
            double[] h_table = hen_polynomial_value(1, p, x_table, +k);
            switch (e)
            {
                //
                //  The following formula is an outer product in H_TABLE.
                //
                case 0:
                {
                    for (j = 0; j <= p; j++)
                    {
                        for (i = 0; i <= p; i++)
                        {
                            table[i + j * (p + 1)] += w_table[k] * h_table[i] * h_table[j];
                        }
                    }

                    break;
                }
                default:
                {
                    for (j = 0; j <= p; j++)
                    {
                        for (i = 0; i <= p; i++)
                        {
                            table[i + j * (p + 1)] += w_table[k] * Math.Pow(x, e) * h_table[i] * h_table[j];
                        }
                    }

                    break;
                }
            }
        }

        return table;
    }

    public static double[] hf_exponential_product(int p, double b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HF_EXPONENTIAL_PRODUCT: exponential products, exp(b*x)*Hf(i,x)*Hf(j,x).
        //
        //  Discussion:
        //
        //    Hf(i,x) represents the Hermite function of "degree" I.   
        //
        //    For polynomial chaos applications, it is of interest to know the
        //    value of the integrals of products of exp(B*X) with every possible pair
        //    of basis functions.  That is, we'd like to form
        //
        //      Tij = Integral ( -oo < X < +oo ) exp(B*X) * Hf(I,X) * Hf(J,X) dx
        //
        //    We will estimate these integrals using Gauss-Hermite quadrature.
        //    Because of the exponential factor exp(B*X), the quadrature will not 
        //    be exact.
        //
        //    However, when B = 0, the quadrature is exact, and moreoever, the
        //    table will be the identity matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 February 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int P, the maximum degree of the polyonomial 
        //    factors.  0 <= P.
        //
        //    Input, double B, the coefficient of X in the exponential factor.
        //
        //    Output, double HF_EXPONENTIAL_PRODUCT[(P+1)*(P+1)], the table of
        //    integrals.  TABLE(I,J) represents the integral of 
        //    exp(B*X) * Hf(I,X) * Hf(J,X).
        //
    {
        int i;
        int j;
        int k;

        double[] table = new double[(p + 1) * (p + 1)];
        for (j = 0; j <= p; j++)
        {
            for (i = 0; i <= p; i++)
            {
                table[i + j * (p + 1)] = 0.0;
            }
        }

        int order = (3 * p + 4) / 2;

        double[] x_table = new double[order];
        double[] w_table = new double[order];

        HermiteQuadrature.hf_quadrature_rule(order, ref x_table, ref w_table);

        for (k = 0; k < order; k++)
        {
            double x = x_table[k];
            double[] h_table = hf_function_value(1, p, x_table, +k);
            //
            //  The following formula is an outer product in H_TABLE.
            //
            for (j = 0; j <= p; j++)
            {
                for (i = 0; i <= p; i++)
                {
                    table[i + j * (p + 1)] += w_table[k] * Math.Exp(b * x) * h_table[i] * h_table[j];
                }
            }
        }

        return table;
    }

    public static double[] hf_function_value(int m, int n, double[] x, int xIndex = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HF_FUNCTION_VALUE evaluates Hf(i,x).
        //
        //  Discussion:
        //
        //    Hf(i,x) represents the Hermite function of "degree" I.   
        //
        //    The Hermite function of degree n is related to the physicist's
        //    Hermite polynomial H(n,x):
        //
        //      Hf(n,x) = H(n,x) * exp ( - 0.5 * x^2 ) / sqrt ( 2^n n// sqrt ( Math.PI ) )
        //
        //    The Hermite functions are orthonormal:
        //
        //      Integral ( -oo < x < +oo ) Hf(m,x) Hf(n,x) dx = delta ( m, n )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 February 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Milton Abramowitz, Irene Stegun,
        //    Handbook of Mathematical Functions,
        //    National Bureau of Standards, 1964,
        //    ISBN: 0-486-61272-4,
        //    LC: QA47.A34.
        //
        //    Frank Olver, Daniel Lozier, Ronald Boisvert, Charles Clark,
        //    NIST Handbook of Mathematical Functions,
        //    Cambridge University Press, 2010,
        //    ISBN: 978-0521192255,
        //    LC: QA331.N57.
        //
        //  Parameters:
        //
        //    Input, int M, the number of evaluation points.
        //
        //    Input, int N, the highest order polynomial to compute.
        //    Note that polynomials 0 through N will be computed.
        //
        //    Input, double X[M], the point at which the polynomials are 
        //    to be evaluated.
        //
        //    Output, double HF_FUNCTION_VALUE[M*(N+1)], the values of the Hermite
        //    functions of index 0 through N at the evaluation points.
        //
    {
        int i;
        int j;
            

        double[] f = new double[m * (n + 1)];

        for (i = 0; i < m; i++)
        {
            f[i + 0 * m] = Math.Exp(-0.5 * x[i] * x[xIndex + i]) / Math.Sqrt(Math.Sqrt(Math.PI));
        }

        switch (n)
        {
            case 0:
                return f;
        }

        for (i = 0; i < m; i++)
        {
            f[i + 1 * m] = 2.0 * Math.Exp(-0.5 * x[i] * x[i]) * x[xIndex + i]
                           / Math.Sqrt(2.0 * Math.Sqrt(Math.PI));
        }

        for (j = 2; j <= n; j++)
        {
            for (i = 0; i < m; i++)
            {
                f[i + j * m] = (Math.Sqrt(2.0) * x[i] * f[i + (j - 1) * m]
                                - Math.Sqrt(j - 1) * f[i + (j - 2) * m])
                               / Math.Sqrt(j);
            }
        }

        return f;
    }


    public static double[] hf_power_product(int p, int e)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HF_POWER_PRODUCT: power products x^e*Hf(i,x)*Hf(j,x).
        //
        //  Discussion:
        //
        //    Hf(i,x) represents the Hermite function of "degree" I.   
        //
        //    For polynomial chaos applications, it is of interest to know the
        //    value of the integrals of products of X with every possible pair
        //    of basis functions.  That is, we'd like to form
        //
        //      Tij = Integral ( -oo < X < +oo ) X^E * Hf(I,X) * Hf(J,X) dx
        //
        //    We will estimate these integrals using Gauss-Hermite quadrature.
        //
        //    When E is 0, the computed table should be the identity matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 February 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int P, the maximum degree of the polyonomial 
        //    factors.  0 <= P.
        //
        //    Input, int E, the exponent of X in the integrand.
        //    0 <= E.
        //
        //    Output, double HF_POWER_PRODUCT[(P+1)*(P+1)], the table of integrals.  
        //    TABLE(I,J) represents the integral of X^E * Hf(I,X) * Hf(J,X).
        //
    {
        int i;
        int j;
        int k;

        double[] table = new double[(p + 1) * (p + 1)];
        for (j = 0; j <= p; j++)
        {
            for (i = 0; i <= p; i++)
            {
                table[i + j * (p + 1)] = 0.0;
            }
        }

        int order = p + 1 + (e + 1) / 2;

        double[] x_table = new double[order];
        double[] w_table = new double[order];

        HermiteQuadrature.hf_quadrature_rule(order, ref x_table, ref w_table);

        for (k = 0; k < order; k++)
        {
            double x = x_table[k];
            double[] h_table = hf_function_value(1, p, x_table, +k);
            switch (e)
            {
                //
                //  The following formula is an outer product in H_TABLE.
                //
                case 0:
                {
                    for (j = 0; j <= p; j++)
                    {
                        for (i = 0; i <= p; i++)
                        {
                            table[i + j * (p + 1)] += w_table[k] * h_table[i] * h_table[j];
                        }
                    }

                    break;
                }
                default:
                {
                    for (j = 0; j <= p; j++)
                    {
                        for (i = 0; i <= p; i++)
                        {
                            table[i + j * (p + 1)] += w_table[k] * Math.Pow(x, e) * h_table[i] * h_table[j];
                        }
                    }

                    break;
                }
            }
        }

        return table;
    }

    public static double[] hn_exponential_product(int p, double b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HN_EXPONENTIAL_PRODUCT: exponential products exp(b*x)*Hn(i,x)*Hn(j,x).
        //
        //  Discussion:
        //
        //    Hn(i,x) is the normalized physicist's Hermite polynomial of degree I.  
        //
        //    For polynomial chaos applications, it is of interest to know the
        //    value of the integrals of products of exp(B*X) with every possible pair
        //    of basis functions.  That is, we'd like to form
        //
        //      Tij = Integral ( -oo < X < +oo ) 
        //        exp(B*X) * Hn(I,X) * Hn(J,X) exp(-X*X) dx
        //
        //    We will estimate these integrals using Gauss-Hermite quadrature.
        //    Because of the exponential factor exp(B*X), the quadrature will not 
        //    be exact.
        //
        //    However, when B = 0, the quadrature is exact, and moreoever, the
        //    table will be the identity matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 February 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int P, the maximum degree of the polyonomial 
        //    factors.  0 <= P.
        //
        //    Input, double B, the coefficient of X in the exponential factor.
        //
        //    Output, double HN_EXPONENTIAL_PRODUCT[(P+1)*(P+1)], the table of 
        //    integrals.  TABLE(I,J) represents the weighted integral of 
        //    exp(B*X) * Hn(I,X) * Hn(J,X).
        //
    {
        int i;
        int j;
        int k;

        double[] table = new double[(p + 1) * (p + 1)];
        for (j = 0; j <= p; j++)
        {
            for (i = 0; i <= p; i++)
            {
                table[i + j * (p + 1)] = 0.0;
            }
        }

        int order = (3 * p + 4) / 2;

        double[] x_table = new double[order];
        double[] w_table = new double[order];

        HermiteQuadrature.h_quadrature_rule(order, ref x_table, ref w_table);

        for (k = 0; k < order; k++)
        {
            double x = x_table[k];
            double[] h_table = hn_polynomial_value(1, p, x_table, +k);
            //
            //  The following formula is an outer product in H_TABLE.
            //
            for (j = 0; j <= p; j++)
            {
                for (i = 0; i <= p; i++)
                {
                    table[i + j * (p + 1)] += w_table[k] * Math.Exp(b * x) * h_table[i] * h_table[j];
                }
            }
        }

        return table;
    }

    public static double[] hn_polynomial_value(int m, int n, double[] x, int xIndex = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HN_POLYNOMIAL_VALUE evaluates Hn(i,x).
        //
        //  Discussion:
        //
        //    Hn(i,x) is the normalized physicist's Hermite polynomial of degree I.  
        //
        //    These polynomials satisfy the orthonormality condition:
        //
        //      Integral ( -oo < X < +oo ) 
        //        exp ( - X^2 ) * Hn(M,X) Hn(N,X) dX = delta ( N, M )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    26 February 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Milton Abramowitz, Irene Stegun,
        //    Handbook of Mathematical Functions,
        //    National Bureau of Standards, 1964,
        //    ISBN: 0-486-61272-4,
        //    LC: QA47.A34.
        //
        //  Parameters:
        //
        //    Input, int M, the number of evaluation points.
        //
        //    Input, int N, the highest order polynomial to compute.
        //    Note that polynomials 0 through N will be computed.
        //
        //    Input, double X[M], the evaluation points.
        //
        //    Output, double HN_POLYNOMIAL_VALUE[M*(N+1)], the values of the first 
        //    N+1 Hermite polynomials at the evaluation points.
        //
    {
        int i;
        int j;

        switch (n)
        {
            case < 0:
                return null;
        }

        double[] p = new double[m * (n + 1)];

        for (i = 0; i < m; i++)
        {
            p[i + 0 * m] = 1.0;
        }

        switch (n)
        {
            case 0:
                return p;
        }

        for (i = 0; i < m; i++)
        {
            p[i + 1 * m] = 2.0 * x[xIndex + i];
        }

        for (j = 2; j <= n; j++)
        {
            for (i = 0; i < m; i++)
            {
                p[i + j * m] = 2.0 * x[xIndex + i] * p[i + (j - 1) * m]
                               - 2.0 * (j - 1) * p[i + (j - 2) * m];
            }
        }

        //
        //  Normalize.
        //
        double fact = 1.0;
        double two_power = 1.0;
        for (j = 0; j <= n; j++)
        {
            for (i = 0; i < m; i++)
            {
                p[i + j * m] /= Math.Sqrt(fact * two_power * Math.Sqrt(Math.PI));
            }

            fact *= j + 1;
            two_power *= 2.0;
        }

        return p;
    }

    public static double[] hn_power_product(int p, int e)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HN_POWER_PRODUCT: power products x^e*Hn(i,x)*Hn(j,x).
        //
        //  Discussion:
        //
        //    Hn(i,x) is the normalized physicist's Hermite polynomial of degree I.  
        //
        //    For polynomial chaos applications, it is of interest to know the
        //    value of the integrals of products of X with every possible pair
        //    of basis functions.  That is, we'd like to form
        //
        //      Tij = Integral ( -oo < X < +oo ) X^E * Hn(I,X) * Hn(J,X) exp(-X*X) dx
        //
        //    We will estimate these integrals using Gauss-Hermite quadrature.
        //
        //    When E is 0, the computed table should be the identity matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 February 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int P, the maximum degree of the polyonomial 
        //    factors.  0 <= P.
        //
        //    Input, int E, the exponent of X in the integrand.
        //    0 <= E.
        //
        //    Output, double HN_POWER_PRODUCT[(P+1)*(P+1)], the table of 
        //    integrals.  TABLE(I,J) represents the weighted integral of 
        //    X^E * Hn(I,X) * Hn(J,X).
        //
    {
        int i;
        int j;
        int k;

        double[] table = new double[(p + 1) * (p + 1)];
        for (j = 0; j <= p; j++)
        {
            for (i = 0; i <= p; i++)
            {
                table[i + j * (p + 1)] = 0.0;
            }
        }

        int order = p + 1 + (e + 1) / 2;

        double[] x_table = new double[order];
        double[] w_table = new double[order];

        HermiteQuadrature.h_quadrature_rule(order, ref x_table, ref w_table);

        for (k = 0; k < order; k++)
        {
            double x = x_table[k];
            double[] h_table = hn_polynomial_value(1, p, x_table, +k);
            switch (e)
            {
                //
                //  The following formula is an outer product in H_TABLE.
                //
                case 0:
                {
                    for (j = 0; j <= p; j++)
                    {
                        for (i = 0; i <= p; i++)
                        {
                            table[i + j * (p + 1)] += w_table[k] * h_table[i] * h_table[j];
                        }
                    }

                    break;
                }
                default:
                {
                    for (j = 0; j <= p; j++)
                    {
                        for (i = 0; i <= p; i++)
                        {
                            table[i + j * (p + 1)] += w_table[k] * Math.Pow(x, e) * h_table[i] * h_table[j];
                        }
                    }

                    break;
                }
            }
        }

        return table;
    }

    public static void hermite_genz_keister_lookup_points(int n, ref double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HERMITE_GENZ_KEISTER_LOOKUP_POINTS looks up Genz-Keister Hermite abscissas.
        //
        //  Discussion:
        //
        //    The integral:
        //
        //      integral ( -oo <= x <= +oo ) f(x) exp ( - x * x ) dx
        //
        //    The quadrature rule:
        //
        //      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
        //
        //    A nested family of rules for the Hermite integration problem
        //    was produced by Genz and Keister.  The structure of the nested
        //    family was denoted by 1+2+6+10+?, that is, it comprised rules
        //    of successive orders O = 1, 3, 9, 19, and a final rule of order
        //    35, 37, 41 or 43.
        //
        //    The precisions of these rules are P = 1, 5, 15, 29, 
        //    with the final rule of precision 51, 55, 63 or 67.
        //
        //    Three related families begin the same way, but end with a different final
        //    rule.  As a convenience, this function includes these final rules as well:
        //
        //    Designation  Orders       Precisions
        //
        //    1+2+6+10+16, 1,3,9,19,35  1,5,15,29,51
        //    1+2+6+10+18  1,3,9,19,37  1,5,15,29,55
        //    1+2+6+10+22  1,3,9,19,41  1,5,15,29,63
        //    1+2+6+10+24  1,3,9,19,43  1,5,15,29,67
        //
        //    Some of the data in this function was kindly supplied directly by
        //    Alan Genz on 24 April 2011.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Alan Genz, Bradley Keister,
        //    Fully symmetric interpolatory rules for multiple integrals
        //    over infinite regions with Gaussian weight,
        //    Journal of Computational and Applied Mathematics,
        //    Volume 71, 1996, pages 299-309
        //
        //    Florian Heiss, Viktor Winschel,
        //    Likelihood approximation by numerical integration on sparse grids,
        //    Journal of Econometrics,
        //    Volume 144, 2008, pages 62-80.
        //
        //    Thomas Patterson,
        //    The Optimal Addition of Points to Quadrature Formulae,
        //    Mathematics of Computation,
        //    Volume 22, Number 104, October 1968, pages 847-856.
        //
        //  Parameters:
        //
        //    Input, int N, the order.
        //    N must be 1, 3, 9, 19, 35, 37, 41, or 43.
        //
        //    Output, double X[N], the abscissas.
        //
    {
        switch (n)
        {
            case 1:
                x[0] = 0.0000000000000000E+00;
                break;
            case 3:
                x[0] = -1.2247448713915889E+00;
                x[1] = 0.0000000000000000E+00;
                x[2] = 1.2247448713915889E+00;
                break;
            case 9:
                x[0] = -2.9592107790638380E+00;
                x[1] = -2.0232301911005157E+00;
                x[2] = -1.2247448713915889E+00;
                x[3] = -5.2403354748695763E-01;
                x[4] = 0.0000000000000000E+00;
                x[5] = 5.2403354748695763E-01;
                x[6] = 1.2247448713915889E+00;
                x[7] = 2.0232301911005157E+00;
                x[8] = 2.9592107790638380E+00;
                break;
            case 19:
                x[0] = -4.4995993983103881E+00;
                x[1] = -3.6677742159463378E+00;
                x[2] = -2.9592107790638380E+00;
                x[3] = -2.2665132620567876E+00;
                x[4] = -2.0232301911005157E+00;
                x[5] = -1.8357079751751868E+00;
                x[6] = -1.2247448713915889E+00;
                x[7] = -8.7004089535290285E-01;
                x[8] = -5.2403354748695763E-01;
                x[9] = 0.0000000000000000E+00;
                x[10] = 5.2403354748695763E-01;
                x[11] = 8.7004089535290285E-01;
                x[12] = 1.2247448713915889E+00;
                x[13] = 1.8357079751751868E+00;
                x[14] = 2.0232301911005157E+00;
                x[15] = 2.2665132620567876E+00;
                x[16] = 2.9592107790638380E+00;
                x[17] = 3.6677742159463378E+00;
                x[18] = 4.4995993983103881E+00;
                break;
            case 35:
                x[0] = -6.3759392709822356E+00;
                x[1] = -5.6432578578857449E+00;
                x[2] = -5.0360899444730940E+00;
                x[3] = -4.4995993983103881E+00;
                x[4] = -4.0292201405043713E+00;
                x[5] = -3.6677742159463378E+00;
                x[6] = -3.3491639537131945E+00;
                x[7] = -2.9592107790638380E+00;
                x[8] = -2.5705583765842968E+00;
                x[9] = -2.2665132620567876E+00;
                x[10] = -2.0232301911005157E+00;
                x[11] = -1.8357079751751868E+00;
                x[12] = -1.5794121348467671E+00;
                x[13] = -1.2247448713915889E+00;
                x[14] = -8.7004089535290285E-01;
                x[15] = -5.2403354748695763E-01;
                x[16] = -1.7606414208200893E-01;
                x[17] = 0.0000000000000000E+00;
                x[18] = 1.7606414208200893E-01;
                x[19] = 5.2403354748695763E-01;
                x[20] = 8.7004089535290285E-01;
                x[21] = 1.2247448713915889E+00;
                x[22] = 1.5794121348467671E+00;
                x[23] = 1.8357079751751868E+00;
                x[24] = 2.0232301911005157E+00;
                x[25] = 2.2665132620567876E+00;
                x[26] = 2.5705583765842968E+00;
                x[27] = 2.9592107790638380E+00;
                x[28] = 3.3491639537131945E+00;
                x[29] = 3.6677742159463378E+00;
                x[30] = 4.0292201405043713E+00;
                x[31] = 4.4995993983103881E+00;
                x[32] = 5.0360899444730940E+00;
                x[33] = 5.6432578578857449E+00;
                x[34] = 6.3759392709822356E+00;
                break;
            case 37:
                x[0] = -6.853200069757519;
                x[1] = -6.124527854622158;
                x[2] = -5.521865209868350;
                x[3] = -4.986551454150765;
                x[4] = -4.499599398310388;
                x[5] = -4.057956316089741;
                x[6] = -3.667774215946338;
                x[7] = -3.315584617593290;
                x[8] = -2.959210779063838;
                x[9] = -2.597288631188366;
                x[10] = -2.266513262056788;
                x[11] = -2.023230191100516;
                x[12] = -1.835707975175187;
                x[13] = -1.561553427651873;
                x[14] = -1.224744871391589;
                x[15] = -0.870040895352903;
                x[16] = -0.524033547486958;
                x[17] = -0.214618180588171;
                x[18] = 0.000000000000000;
                x[19] = 0.214618180588171;
                x[20] = 0.524033547486958;
                x[21] = 0.870040895352903;
                x[22] = 1.224744871391589;
                x[23] = 1.561553427651873;
                x[24] = 1.835707975175187;
                x[25] = 2.023230191100516;
                x[26] = 2.266513262056788;
                x[27] = 2.597288631188366;
                x[28] = 2.959210779063838;
                x[29] = 3.315584617593290;
                x[30] = 3.667774215946338;
                x[31] = 4.057956316089741;
                x[32] = 4.499599398310388;
                x[33] = 4.986551454150765;
                x[34] = 5.521865209868350;
                x[35] = 6.124527854622158;
                x[36] = 6.853200069757519;
                break;
            case 41:
                x[0] = -7.251792998192644;
                x[1] = -6.547083258397540;
                x[2] = -5.961461043404500;
                x[3] = -5.437443360177798;
                x[4] = -4.953574342912980;
                x[5] = -4.4995993983103881;
                x[6] = -4.070919267883068;
                x[7] = -3.6677742159463378;
                x[8] = -3.296114596212218;
                x[9] = -2.9592107790638380;
                x[10] = -2.630415236459871;
                x[11] = -2.2665132620567876;
                x[12] = -2.043834754429505;
                x[13] = -2.0232301911005157;
                x[14] = -1.8357079751751868;
                x[15] = -1.585873011819188;
                x[16] = -1.2247448713915889;
                x[17] = -0.87004089535290285;
                x[18] = -0.52403354748695763;
                x[19] = -0.195324784415805;
                x[20] = 0.0000000000000000;
                x[21] = 0.195324784415805;
                x[22] = 0.52403354748695763;
                x[23] = 0.87004089535290285;
                x[24] = 1.2247448713915889;
                x[25] = 1.585873011819188;
                x[26] = 1.8357079751751868;
                x[27] = 2.0232301911005157;
                x[28] = 2.043834754429505;
                x[29] = 2.2665132620567876;
                x[30] = 2.630415236459871;
                x[31] = 2.9592107790638380;
                x[32] = 3.296114596212218;
                x[33] = 3.6677742159463378;
                x[34] = 4.070919267883068;
                x[35] = 4.4995993983103881;
                x[36] = 4.953574342912980;
                x[37] = 5.437443360177798;
                x[38] = 5.961461043404500;
                x[39] = 6.547083258397540;
                x[40] = 7.251792998192644;
                break;
            case 43:
                x[0] = -10.167574994881873;
                x[1] = -7.231746029072501;
                x[2] = -6.535398426382995;
                x[3] = -5.954781975039809;
                x[4] = -5.434053000365068;
                x[5] = -4.952329763008589;
                x[6] = -4.4995993983103881;
                x[7] = -4.071335874253583;
                x[8] = -3.6677742159463378;
                x[9] = -3.295265921534226;
                x[10] = -2.9592107790638380;
                x[11] = -2.633356763661946;
                x[12] = -2.2665132620567876;
                x[13] = -2.089340389294661;
                x[14] = -2.0232301911005157;
                x[15] = -1.8357079751751868;
                x[16] = -1.583643465293944;
                x[17] = -1.2247448713915889;
                x[18] = -0.87004089535290285;
                x[19] = -0.52403354748695763;
                x[20] = -0.196029453662011;
                x[21] = 0.0000000000000000;
                x[22] = 0.196029453662011;
                x[23] = 0.52403354748695763;
                x[24] = 0.87004089535290285;
                x[25] = 1.2247448713915889;
                x[26] = 1.583643465293944;
                x[27] = 1.8357079751751868;
                x[28] = 2.0232301911005157;
                x[29] = 2.089340389294661;
                x[30] = 2.2665132620567876;
                x[31] = 2.633356763661946;
                x[32] = 2.9592107790638380;
                x[33] = 3.295265921534226;
                x[34] = 3.6677742159463378;
                x[35] = 4.071335874253583;
                x[36] = 4.4995993983103881;
                x[37] = 4.952329763008589;
                x[38] = 5.434053000365068;
                x[39] = 5.954781975039809;
                x[40] = 6.535398426382995;
                x[41] = 7.231746029072501;
                x[42] = 10.167574994881873;
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("HERMITE_GENZ_KEISTER_LOOKUP_POINTS - Fatal error!");
                Console.WriteLine("  Illegal input value of N.");
                Console.WriteLine("  N must be 1, 3, 9, 19, 35, 37, 41 or 43.");
                break;
        }
    }

    public static void hermite_genz_keister_lookup_weights(int n, ref double[] w)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HERMITE_GENZ_KEISTER_LOOKUP_WEIGHTS looks up Genz-Keister Hermite weights.
        //
        //  Discussion:
        //
        //    The integral:
        //
        //      integral ( -oo <= x <= +oo ) f(x) exp ( - x * x ) dx
        //
        //    The quadrature rule:
        //
        //      sum ( 1 <= i <= n ) w(i) * f ( x(i) )
        //
        //    A nested family of rules for the Hermite integration problem
        //    was produced by Genz and Keister.  The structure of the nested
        //    family was denoted by 1+2+6+10+?, that is, it comprised rules
        //    of successive orders O = 1, 3, 9, 19, and a final rule of order
        //    35, 37, 41 or 43.
        //
        //    The precisions of these rules are P = 1, 5, 15, 29, 
        //    with the final rule of precision 51, 55, 63 or 67.
        //
        //    Three related families begin the same way, but end with a different final
        //    rule.  As a convenience, this function includes these final rules as well:
        //
        //    Designation  Orders       Precisions
        //
        //    1+2+6+10+16, 1,3,9,19,35  1,5,15,29,51
        //    1+2+6+10+18  1,3,9,19,37  1,5,15,29,55
        //    1+2+6+10+22  1,3,9,19,41  1,5,15,29,63
        //    1+2+6+10+24  1,3,9,19,43  1,5,15,29,67
        //
        //    Some of the data in this function was kindly supplied directly by
        //    Alan Genz on 24 April 2011.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Alan Genz, Bradley Keister,
        //    Fully symmetric interpolatory rules for multiple integrals
        //    over infinite regions with Gaussian weight,
        //    Journal of Computational and Applied Mathematics,
        //    Volume 71, 1996, pages 299-309
        //
        //    Florian Heiss, Viktor Winschel,
        //    Likelihood approximation by numerical integration on sparse grids,
        //    Journal of Econometrics,
        //    Volume 144, 2008, pages 62-80.
        //
        //    Thomas Patterson,
        //    The Optimal Addition of Points to Quadrature Formulae,
        //    Mathematics of Computation,
        //    Volume 22, Number 104, October 1968, pages 847-856.
        //
        //  Parameters:
        //
        //    Input, int N, the order.
        //    N must be 1, 3, 9, 19, 35, 37, 41, or 43.
        //
        //    Output, double W[N], the weights.
        //
    {
        switch (n)
        {
            case 1:
                w[0] = 1.7724538509055159E+00;
                break;
            case 3:
                w[0] = 2.9540897515091930E-01;
                w[1] = 1.1816359006036772E+00;
                w[2] = 2.9540897515091930E-01;
                break;
            case 9:
                w[0] = 1.6708826306882348E-04;
                w[1] = 1.4173117873979098E-02;
                w[2] = 1.6811892894767771E-01;
                w[3] = 4.7869428549114124E-01;
                w[4] = 4.5014700975378197E-01;
                w[5] = 4.7869428549114124E-01;
                w[6] = 1.6811892894767771E-01;
                w[7] = 1.4173117873979098E-02;
                w[8] = 1.6708826306882348E-04;
                break;
            case 19:
                w[0] = 1.5295717705322357E-09;
                w[1] = 1.0802767206624762E-06;
                w[2] = 1.0656589772852267E-04;
                w[3] = 5.1133174390883855E-03;
                w[4] = -1.1232438489069229E-02;
                w[5] = 3.2055243099445879E-02;
                w[6] = 1.1360729895748269E-01;
                w[7] = 1.0838861955003017E-01;
                w[8] = 3.6924643368920851E-01;
                w[9] = 5.3788160700510168E-01;
                w[10] = 3.6924643368920851E-01;
                w[11] = 1.0838861955003017E-01;
                w[12] = 1.1360729895748269E-01;
                w[13] = 3.2055243099445879E-02;
                w[14] = -1.1232438489069229E-02;
                w[15] = 5.1133174390883855E-03;
                w[16] = 1.0656589772852267E-04;
                w[17] = 1.0802767206624762E-06;
                w[18] = 1.5295717705322357E-09;
                break;
            case 35:
                w[0] = 1.8684014894510604E-18;
                w[1] = 9.6599466278563243E-15;
                w[2] = 5.4896836948499462E-12;
                w[3] = 8.1553721816916897E-10;
                w[4] = 3.7920222392319532E-08;
                w[5] = 4.3737818040926989E-07;
                w[6] = 4.8462799737020461E-06;
                w[7] = 6.3328620805617891E-05;
                w[8] = 4.8785399304443770E-04;
                w[9] = 1.4515580425155904E-03;
                w[10] = 4.0967527720344047E-03;
                w[11] = 5.5928828911469180E-03;
                w[12] = 2.7780508908535097E-02;
                w[13] = 8.0245518147390893E-02;
                w[14] = 1.6371221555735804E-01;
                w[15] = 2.6244871488784277E-01;
                w[16] = 3.3988595585585218E-01;
                w[17] = 9.1262675363737921E-04;
                w[18] = 3.3988595585585218E-01;
                w[19] = 2.6244871488784277E-01;
                w[20] = 1.6371221555735804E-01;
                w[21] = 8.0245518147390893E-02;
                w[22] = 2.7780508908535097E-02;
                w[23] = 5.5928828911469180E-03;
                w[24] = 4.0967527720344047E-03;
                w[25] = 1.4515580425155904E-03;
                w[26] = 4.8785399304443770E-04;
                w[27] = 6.3328620805617891E-05;
                w[28] = 4.8462799737020461E-06;
                w[29] = 4.3737818040926989E-07;
                w[30] = 3.7920222392319532E-08;
                w[31] = 8.1553721816916897E-10;
                w[32] = 5.4896836948499462E-12;
                w[33] = 9.6599466278563243E-15;
                w[34] = 1.8684014894510604E-18;
                break;
            case 37:
                w[0] = 0.337304188079177058E-20;
                w[1] = 0.332834739632930463E-16;
                w[2] = 0.323016866782871498E-13;
                w[3] = 0.809333688669950037E-11;
                w[4] = 0.748907559239519284E-09;
                w[5] = 0.294146671497083432E-07;
                w[6] = 0.524482423744884136E-06;
                w[7] = 0.586639457073896277E-05;
                w[8] = 0.571885531470621903E-04;
                w[9] = 0.41642095727577091E-03;
                w[10] = 0.174733389581099482E-02;
                w[11] = 0.313373786000304381E-02;
                w[12] = 0.768092665770660459E-02;
                w[13] = 0.274962713372148476E-01;
                w[14] = 0.783630990508037449E-01;
                w[15] = 0.16611584261479281E+00;
                w[16] = 0.253636910481387185E+00;
                w[17] = 0.261712932511430884E+00;
                w[18] = 0.171719680968980257E+00;
                w[19] = 0.261712932511430884E+00;
                w[20] = 0.253636910481387185E+00;
                w[21] = 0.16611584261479281E+00;
                w[22] = 0.783630990508037449E-01;
                w[23] = 0.274962713372148476E-01;
                w[24] = 0.768092665770660459E-02;
                w[25] = 0.313373786000304381E-02;
                w[26] = 0.174733389581099482E-02;
                w[27] = 0.41642095727577091E-03;
                w[28] = 0.571885531470621903E-04;
                w[29] = 0.586639457073896277E-05;
                w[30] = 0.524482423744884136E-06;
                w[31] = 0.294146671497083432E-07;
                w[32] = 0.748907559239519284E-09;
                w[33] = 0.809333688669950037E-11;
                w[34] = 0.323016866782871498E-13;
                w[35] = 0.332834739632930463E-16;
                w[36] = 0.337304188079177058E-20;
                break;
            case 41:
                w[0] = 0.117725656974405367E-22;
                w[1] = 0.152506745534300636E-18;
                w[2] = 0.202183949965101288E-15;
                w[3] = 0.724614869051195508E-13;
                w[4] = 0.103121966469463034E-10;
                w[5] = 0.710371395169350952E-09;
                w[6] = 0.264376044449260516E-07;
                w[7] = 0.558982787078644997E-06;
                w[8] = 0.675628907134744976E-05;
                w[9] = 0.512198007019776873E-04;
                w[10] = 0.335013114947200879E-03;
                w[11] = 0.249379691096933139E-02;
                w[12] = -0.25616995850607458E-01;
                w[13] = 0.317007878644325588E-01;
                w[14] = 0.125041498584003435E-02;
                w[15] = 0.293244560924894295E-01;
                w[16] = 0.799536390803302298E-01;
                w[17] = 0.164543666806555251E+00;
                w[18] = 0.258718519718241095E+00;
                w[19] = 0.293588795735908566E+00;
                w[20] = 0.997525375254611951E-01;
                w[21] = 0.293588795735908566E+00;
                w[22] = 0.258718519718241095E+00;
                w[23] = 0.164543666806555251E+00;
                w[24] = 0.799536390803302298E-01;
                w[25] = 0.293244560924894295E-01;
                w[26] = 0.125041498584003435E-02;
                w[27] = 0.317007878644325588E-01;
                w[28] = -0.25616995850607458E-01;
                w[29] = 0.249379691096933139E-02;
                w[30] = 0.335013114947200879E-03;
                w[31] = 0.512198007019776873E-04;
                w[32] = 0.675628907134744976E-05;
                w[33] = 0.558982787078644997E-06;
                w[34] = 0.264376044449260516E-07;
                w[35] = 0.710371395169350952E-09;
                w[36] = 0.103121966469463034E-10;
                w[37] = 0.724614869051195508E-13;
                w[38] = 0.202183949965101288E-15;
                w[39] = 0.152506745534300636E-18;
                w[40] = 0.117725656974405367E-22;
                break;
            case 43:
                w[0] = 0.968100020641528185E-37;
                w[1] = 0.15516931262860431E-22;
                w[2] = 0.175937309107750992E-18;
                w[3] = 0.217337608710893738E-15;
                w[4] = 0.747837010380540069E-13;
                w[5] = 0.104028132097205732E-10;
                w[6] = 0.70903573389336778E-09;
                w[7] = 0.263481722999966618E-07;
                w[8] = 0.560127964848432175E-06;
                w[9] = 0.680410934802210232E-05;
                w[10] = 0.508343873102544037E-04;
                w[11] = 0.32753080006610181E-03;
                w[12] = 0.267479828788552937E-02;
                w[13] = -0.687704270963253854E-02;
                w[14] = 0.119383201790913588E-01;
                w[15] = 0.248083722871002796E-02;
                w[16] = 0.29000335749726387E-01;
                w[17] = 0.798689557875757008E-01;
                w[18] = 0.164609842422580606E+00;
                w[19] = 0.258535954731607738E+00;
                w[20] = 0.292243810406117141E+00;
                w[21] = 0.102730713753441829E+00;
                w[22] = 0.292243810406117141E+00;
                w[23] = 0.258535954731607738E+00;
                w[24] = 0.164609842422580606E+00;
                w[25] = 0.798689557875757008E-01;
                w[26] = 0.29000335749726387E-01;
                w[27] = 0.248083722871002796E-02;
                w[28] = 0.119383201790913588E-01;
                w[29] = -0.687704270963253854E-02;
                w[30] = 0.267479828788552937E-02;
                w[31] = 0.32753080006610181E-03;
                w[32] = 0.508343873102544037E-04;
                w[33] = 0.680410934802210232E-05;
                w[34] = 0.560127964848432175E-06;
                w[35] = 0.263481722999966618E-07;
                w[36] = 0.70903573389336778E-09;
                w[37] = 0.104028132097205732E-10;
                w[38] = 0.747837010380540069E-13;
                w[39] = 0.217337608710893738E-15;
                w[40] = 0.175937309107750992E-18;
                w[41] = 0.15516931262860431E-22;
                w[42] = 0.968100020641528185E-37;
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("HERMITE_GENZ_KEISTER_LOOKUP_WEIGHTS - Fatal error!");
                Console.WriteLine("  Illegal input value of N.");
                Console.WriteLine("  N must be 1, 3, 9, 19, 35, 37, 41 or 43.");
                break;
        }
    }

}