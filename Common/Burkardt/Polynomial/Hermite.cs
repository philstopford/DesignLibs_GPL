using System;

namespace Burkardt.PolynomialNS
{
    public static class Hermite
    {
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
            double dq0;
            double dq1;
            double dq2;
            double q0;
            double q1;
            double q2;

            q1 = 1.0;
            dq1 = 0.0;

            q2 = x;
            dq2 = 1.0;

            for (i = 2; i <= order; i++)
            {
                q0 = q1;
                dq0 = dq1;

                q1 = q2;
                dq1 = dq2;

                q2 = x * q1 - 0.5 * ((double) (i) - 1.0) * q0;
                dq2 = x * dq1 + q1 - 0.5 * ((double) (i) - 1.0) * dq0;
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
            double d;
            double eps = 1.0E-12;
            double p2 = 0;
            int step;
            int step_max = 10;

            for (step = 1; step <= step_max; step++)
            {
                hermite_recur(ref p2, ref dp2, ref p1, x, order);

                d = p2 / (dp2);
                x = x - d;

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
            double[] ct;
            int i;
            int j;
            int k;

            ct = new double[(n + 1) * (n + 1)];

            for (i = 0; i <= n; i++)
            {
                for (j = 0; j <= n; j++)
                {
                    ct[i + j * (n + 1)] = 0.0;
                }
            }

            ct[0 + 0 * (n + 1)] = 1.0;

            if (0 < n)
            {
                ct[1 + 1 * (n + 1)] = 1.0;

                for (i = 2; i <= n; i++)
                {
                    ct[i + 0 * (n + 1)] = -(double) (i - 1) * ct[i - 2 + 0 * (n + 1)];
                    for (j = 1; j <= i - 2; j++)
                    {
                        ct[i + j * (n + 1)] =
                            ct[i - 1 + (j - 1) * (n + 1)] - (double) (i - 1) * ct[i - 2 + j * (n + 1)];
                    }

                    ct[i + (i - 1) * (n + 1)] = ct[i - 1 + (i - 2) * (n + 1)];
                    ct[i + i * (n + 1)] = ct[i - 1 + (i - 1) * (n + 1)];
                }
            }

            //
            //  Extract the nonzero data from the alternating columns of the last row.
            //
            o = (n + 2) / 2;

            k = o;
            for (j = n; 0 <= j; j = j - 2)
            {
                k = k - 1;
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
            int j;
            double[] v;
            double[] vtable;

            vtable = new double[n * (o + 1)];

            for (i = 0; i < n; i++)
            {
                vtable[i + 0 * n] = 1.0;
            }

            if (1 <= o)
            {
                for (i = 0; i < n; i++)
                {
                    vtable[i + 1 * n] = x[i];
                }

                for (j = 2; j <= o; j++)
                {
                    for (i = 0; i < n; i++)
                    {
                        vtable[i + j * n] = x[i] * vtable[i + (j - 1) * n]
                                            - (double) (j - 1) * vtable[i + (j - 2) * n];
                    }
                }
            }

            v = new double[n];

            for (i = 0; i < n; i++)
            {
                v[i] = vtable[i + o * n];
            }

            return v;
        }

        public static void hep_values(ref int n_data, ref int n, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HEP_VALUES returns values of the Hermite polynomials He(n,x).
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
            //    Stephen Wolfram,
            //    The Mathematica Book,
            //    Fourth Edition,
            //    Cambridge University Press, 1999,
            //    ISBN: 0-521-64314-7,
            //    LC: QA76.95.W65.
            //
            //  Parameters:
            //
            //    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, int &N, the order of the function.
            //
            //    Output, double &X, the point where the function is evaluated.
            //
            //    Output, double &FX, the value of the function.
            //
        {
            int N_MAX = 18;

            double[] fx_vec =
            {
                1.000000000000000E+00,
                5.000000000000000E+00,
                24.00000000000000E+00,
                110.0000000000000E+00,
                478.0000000000000E+00,
                1950.000000000000E+00,
                7360.000000000000E+00,
                25100.00000000000E+00,
                73980.00000000000E+00,
                169100.0000000000E+00,
                179680.0000000000E+00,
                -792600.0000000000E+00,
                -5939480.000000000E+00,
                0.000000000000000E+00,
                6.281250000000000E+00,
                6.000000000000000E+00,
                18.00000000000000E+00,
                90150.00000000000E+00
            };

            int[] n_vec =
            {
                0, 1, 2,
                3, 4, 5,
                6, 7, 8,
                9, 10, 11,
                12, 5, 5,
                5, 5, 5
            };

            double[] x_vec =
            {
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                0.0E+00,
                0.5E+00,
                1.0E+00,
                3.0E+00,
                1.0E+01
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                n = 0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                n = n_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
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
            double[] c1;
            double[] c2;
            int[] e1;
            int[] e2;
            int[] f2;
            int i;
            int i1;
            int i2;
            int j1;
            int j2;
            int o1;
            int o2 = 0;
            int[] p = new int[1];
            int[] pp;

            c1 = new double[o_max];
            c2 = new double[o_max];
            e1 = new int[o_max];
            e2 = new int[o_max];
            f2 = new int[o_max];
            pp = new int[m];

            o1 = 1;
            c1[0] = 1.0;
            e1[0] = 1;
            //
            //  Implicate one factor at a time.
            //
            for (i = 0; i < m; i++)
            {
                hep_coefficients(l[i], ref o2, ref c2, ref f2);

                o = 0;

                for (j2 = 0; j2 < o2; j2++)
                {
                    for (j1 = 0; j1 < o1; j1++)
                    {
                        c[o] = c1[j1] * c2[j2];
                        if (0 < i)
                        {
                            p = Monomial.mono_unrank_grlex(i, e1[j1]);
                        }

                        for (i2 = 0; i2 < i; i2++)
                        {
                            pp[i2] = p[i2];
                        }

                        pp[i] = f2[j2];
                        e[o] = Monomial.mono_rank_grlex(i + 1, pp);
                        o = o + 1;
                        if (0 < i)
                        {
                            p = null;
                        }
                    }
                }

                Polynomial.polynomial_sort(o, c, e);
                Polynomial.polynomial_compress(o, c, e, ref o, ref c, ref e);

                o1 = o;
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
            double[] v;
            double[] vi;
            double[] xi;

            v = new double[n];

            for (j = 0; j < n; j++)
            {
                v[j] = 1.0;
            }

            xi = new double[n];

            for (i = 0; i < m; i++)
            {
                for (j = 0; j < n; j++)
                {
                    xi[j] = x[i + j * m];
                }

                vi = hep_value(n, o[i], ref xi);
                for (j = 0; j < n; j++)
                {
                    v[j] = v[j] * vi[j];
                }
            }

            return v;
        }

    }
}