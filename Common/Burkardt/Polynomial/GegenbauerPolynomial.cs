using System;
using Burkardt.Types;

namespace Burkardt.PolynomialNS
{
    public static class GegenbauerPolynomial
    {
        public static bool gegenbauer_alpha_check(double alpha)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GEGENBAUER_ALPHA_CHECK checks the value of ALPHA.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    30 November 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double ALPHA, a parameter which is part of the definition of
            //    the Gegenbauer polynomials.  It must be greater than -0.5.
            //
            //    Output, bool GEGENBAUER_ALPHA_CHECK.
            //    TRUE, ALPHA is acceptable.
            //    FALSE, ALPHA is not acceptable. 
            //
        {
            bool check;
            bool squawk;

            squawk = false;

            if (-0.5 < alpha)
            {
                check = true;
            }
            else
            {
                check = false;
                if (squawk)
                {
                    Console.WriteLine("");
                    Console.WriteLine("GEGENBAUER_ALPHA_CHECK - Fatal error!");
                    Console.WriteLine("  Illegal value of ALPHA.");
                    Console.WriteLine("  ALPHA = " + alpha + "");
                    Console.WriteLine("  but ALPHA must be greater than -0.5.");
                }
            }

            return check;
        }

        public static void gegenbauer_ek_compute(int n, double alpha, ref double[] x, ref double[] w)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GEGENBAUER_EK_COMPUTE computes a Gauss-Gegenbauer quadrature rule.
            //
            //  Discussion:
            //
            //    The integral:
            //
            //      Integral ( -1 <= X <= 1 ) (1-X^2)^ALPHA * F(X) dX
            //
            //    The quadrature rule:
            //
            //      Sum ( 1 <= I <= N ) WEIGHT(I) * F ( XTAB(I) )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    20 November 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Sylvan Elhay, Jaroslav Kautsky,
            //    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of 
            //    Interpolatory Quadrature,
            //    ACM Transactions on Mathematical Software,
            //    Volume 13, Number 4, December 1987, pages 399-415.
            //
            //  Parameters:
            //
            //    Input, int N, the order of the quadrature rule.
            //
            //    Input, double ALPHA, the exponent of (1-X^2) in the weight.  
            //    -1.0 < ALPHA is required.
            //
            //    Input, double A, B, the left and right endpoints 
            //    of the interval.
            //
            //    Output, double X[N], the abscissas.
            //
            //    Output, double W[N], the weights.
            //
        {
            double abi;
            double[] bj;
            bool check;
            int i;
            double zemu;
            //
            //  Check N.
            //
            if (n < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("GEGENBAUER_EK_COMPUTE - Fatal error!");
                Console.WriteLine("  1 <= N is required.");
                return;
            }

            //
            //  Check ALPHA.
            //
            check = gegenbauer_alpha_check(alpha);
            if (!check)
            {
                Console.WriteLine("");
                Console.WriteLine("GEGENBAUER_EK_COMPUTE - Fatal error!");
                Console.WriteLine("  Illegal value of ALPHA.");
                return;
            }

            //
            //  Define the zero-th moment.
            //
            zemu = Math.Pow(2.0, 2.0 * alpha + 1.0)
                   * typeMethods.r8_gamma(alpha + 1.0)
                   * typeMethods.r8_gamma(alpha + 1.0)
                   / typeMethods.r8_gamma(2.0 * alpha + 2.0);
            //
            //  Define the Jacobi matrix.
            //
            for (i = 0; i < n; i++)
            {
                x[i] = 0.0;
            }

            bj = new double[n];

            bj[0] = 4.0 * Math.Pow(alpha + 1.0, 2)
                    / ((2.0 * alpha + 3.0) * Math.Pow(2.0 * alpha + 2.0, 2));

            for (i = 2; i <= n; i++)
            {
                abi = 2.0 * (alpha + (double) i);
                bj[i - 1] = 4.0 * (double) (i) * Math.Pow(alpha + i, 2) * (2.0 * alpha + i)
                            / ((abi - 1.0) * (abi + 1.0) * abi * abi);
            }

            for (i = 0; i < n; i++)
            {
                bj[i] = Math.Sqrt(bj[i]);
            }

            w[0] = Math.Sqrt(zemu);
            for (i = 1; i < n; i++)
            {
                w[i] = 0.0;
            }

            //
            //  Diagonalize the Jacobi matrix.
            //
            IMTQLX.imtqlx(n, ref x, ref bj, ref w);

            for (i = 0; i < n; i++)
            {
                w[i] = Math.Pow(w[i], 2);
            }
        }

        public static double[] gegenbauer_polynomial_value(int m, int n, double alpha, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GEGENBAUER_POLYNOMIAL_VALUE computes the Gegenbauer polynomials C(I,ALPHA)(X).
            //
            //  Differential equation:
            //
            //    (1-X*X) Y'' - (2 ALPHA + 1) X Y' + M (M + 2 ALPHA) Y = 0
            //
            //  Recursion:
            //
            //    C(0,ALPHA,X) = 1,
            //    C(1,ALPHA,X) = 2*ALPHA*X
            //    C(M,ALPHA,X) = (  ( 2*M-2+2*ALPHA) * X * C(M-1,ALPHA,X) 
            //                    + (  -M+2-2*ALPHA)   *   C(M-2,ALPHA,X) ) / M
            //
            //  Restrictions:
            //
            //    ALPHA must be greater than -0.5.
            //
            //  Special values:
            //
            //    If ALPHA = 1, the Gegenbauer polynomials reduce to the Chebyshev
            //    polynomials of the second kind.
            //
            //  Norm:
            //
            //    Integral ( -1 <= X <= 1 ) ( 1 - X^2 )^( ALPHA - 0.5 ) * C(M,ALPHA,X)^2 dX
            //
            //    = PI * 2^( 1 - 2 * ALPHA ) * Gamma ( M + 2 * ALPHA ) 
            //      / ( M! * ( M + ALPHA ) * ( Gamma ( ALPHA ) )^2 )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    01 December 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Stephen Wolfram,
            //    The Mathematica Book,
            //    Fourth Edition,
            //    Wolfram Media / Cambridge University Press, 1999.
            //
            //  Parameters:
            //
            //    Input, int M, the highest order polynomial to compute.
            //    Note that polynomials 0 through N will be computed.
            //
            //    Input, int N, the number of evaluation points.
            //
            //    Input, double ALPHA, a parameter which is part of the definition of
            //    the Gegenbauer polynomials.  It must be greater than -0.5.
            //
            //    Input, double X[N], the evaluation points.
            //
            //    Output, double GEGENBAUER_POLYNOMIAL_VALUE(1:M+1,N), the values of 
            //    Gegenbauer polynomials 0 through M
            //    at the N points X.  
            //
        {
            double[] c;
            bool check;
            int i;
            double i_r8;
            int j;

            check = gegenbauer_alpha_check(alpha);
            if (!check)
            {
                Console.WriteLine("");
                Console.WriteLine("GEGENBAUER_POLYNOMIAL_VALUE - Fatal error!");
                Console.WriteLine("  Illegal value of ALPHA.");
                return null;
            }

            c = new double[(m + 1) * n];

            if (m < 0)
            {
                return c;
            }

            if (n == 0)
            {
                return c;
            }

            for (j = 0; j < n; j++)
            {
                c[0 + j * (m + 1)] = 1.0;
            }

            if (m < 1)
            {
                return c;
            }

            for (j = 0; j < n; j++)
            {
                c[1 + j * (m + 1)] = 2.0 * alpha * x[j];
            }

            for (i = 2; i <= m; i++)
            {
                i_r8 = (double) i;
                for (j = 0; j < n; j++)
                {
                    c[i + j * (m + 1)] = ((2.0 * i_r8 - 2.0 + 2.0 * alpha) * x[j] * c[i - 1 + j * (m + 1)]
                                          + (-i_r8 + 2.0 - 2.0 * alpha) * c[i - 2 + j * (m + 1)])
                                         / i_r8;
                }
            }

            return c;
        }

        public static void gegenbauer_polynomial_values(ref int n_data, ref int n, ref double a, ref double x,
                ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GEGENBAUER_POLYNOMIAL_VALUES returns some values of the Gegenbauer polynomials.
            //
            //  Discussion:
            //
            //    The Gegenbauer polynomials are also known as the "spherical
            //    polynomials" or "ultraspherical polynomials".
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      GegenbauerC[n,m,x]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 August 2004
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
            //    Output, int &N, the order parameter of the function.
            //
            //    Output, double &A, the real parameter of the function.
            //
            //    Output, double &X, the argument of the function.
            //
            //    Output, double &FX, the value of the function.
            //
        {
            int N_MAX = 38;

            double[] a_vec =
            {
                0.5E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                0.0E+00,
                1.0E+00,
                2.0E+00,
                3.0E+00,
                4.0E+00,
                5.0E+00,
                6.0E+00,
                7.0E+00,
                8.0E+00,
                9.0E+00,
                10.0E+00,
                3.0E+00,
                3.0E+00,
                3.0E+00,
                3.0E+00,
                3.0E+00,
                3.0E+00,
                3.0E+00,
                3.0E+00,
                3.0E+00,
                3.0E+00,
                3.0E+00,
                3.0E+00,
                3.0E+00,
                3.0E+00,
                3.0E+00,
                3.0E+00
            };

            double[] fx_vec =
            {
                1.0000000000E+00,
                0.2000000000E+00,
                -0.4400000000E+00,
                -0.2800000000E+00,
                0.2320000000E+00,
                0.3075200000E+00,
                -0.0805760000E+00,
                -0.2935168000E+00,
                -0.0395648000E+00,
                0.2459712000E+00,
                0.1290720256E+00,
                0.0000000000E+00,
                -0.3600000000E+00,
                -0.0800000000E+00,
                0.8400000000E+00,
                2.4000000000E+00,
                4.6000000000E+00,
                7.4400000000E+00,
                10.9200000000E+00,
                15.0400000000E+00,
                19.8000000000E+00,
                25.2000000000E+00,
                -9.0000000000E+00,
                -0.1612800000E+00,
                -6.6729600000E+00,
                -8.3750400000E+00,
                -5.5267200000E+00,
                0.0000000000E+00,
                5.5267200000E+00,
                8.3750400000E+00,
                6.6729600000E+00,
                0.1612800000E+00,
                -9.0000000000E+00,
                -15.4252800000E+00,
                -9.6969600000E+00,
                22.4409600000E+00,
                100.8892800000E+00,
                252.0000000000E+00
            };

            int[] n_vec =
            {
                0, 1, 2,
                3, 4, 5,
                6, 7, 8,
                9, 10, 2,
                2, 2, 2,
                2, 2, 2,
                2, 2, 2,
                2, 5, 5,
                5, 5, 5,
                5, 5, 5,
                5, 5, 5,
                5, 5, 5,
                5, 5
            };

            double[] x_vec =
            {
                0.20E+00,
                0.20E+00,
                0.20E+00,
                0.20E+00,
                0.20E+00,
                0.20E+00,
                0.20E+00,
                0.20E+00,
                0.20E+00,
                0.20E+00,
                0.20E+00,
                0.40E+00,
                0.40E+00,
                0.40E+00,
                0.40E+00,
                0.40E+00,
                0.40E+00,
                0.40E+00,
                0.40E+00,
                0.40E+00,
                0.40E+00,
                0.40E+00,
                -0.50E+00,
                -0.40E+00,
                -0.30E+00,
                -0.20E+00,
                -0.10E+00,
                0.00E+00,
                0.10E+00,
                0.20E+00,
                0.30E+00,
                0.40E+00,
                0.50E+00,
                0.60E+00,
                0.70E+00,
                0.80E+00,
                0.90E+00,
                1.00E+00
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
                a = 0.0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                n = n_vec[n_data - 1];
                a = a_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }
    }
}