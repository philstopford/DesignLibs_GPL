using System;
using Burkardt.MatrixNS;
using Burkardt.Types;

namespace Burkardt.PolynomialNS
{
    public static class GegenbauerPolynomial
    {
        public static void gegenbauer_recur ( ref double p2, ref double dp2, ref double p1, double x,
                int order, double alpha, double[] c )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GEGENBAUER_RECUR evaluates a Gegenbauer polynomial.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 February 2008
            //
            //  Author:
            //
            //    John Burkardt
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
            //    Output, double *P2, the value of J(ORDER)(X).
            //
            //    Output, double *DP2, the value of J'(ORDER)(X).
            //
            //    Output, double *P1, the value of J(ORDER-1)(X).
            //
            //    Input, double X, the point at which polynomials are evaluated.
            //
            //    Input, int ORDER, the order of the polynomial.
            //
            //    Input, double ALPHA, the exponents of (1-X^2).
            //
            //    Input, double C[ORDER], the recursion coefficients.
            //
        {
            double dp0;
            double dp1;
            int i;
            double p0;

            p1 = 1.0;
            dp1 = 0.0;

            p2 = x;
            dp2 = 1.0;

            for ( i = 2; i <= order; i++ )
            {
                p0 = p1;
                dp0 = dp1;

                p1 = p2;
                dp1 = dp2;

                p2 = x *  ( p1 ) - c[i-1] * p0;
                dp2 = x * dp1 + ( p1 ) - c[i-1] * dp0;
            }
        }
        public static void gegenbauer_root ( ref double x, int order, double alpha, ref double dp2,
                ref double p1, double[] c )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GEGENBAUER_ROOT improves an approximate root of a Gegenbauer polynomial.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 February 2008
            //
            //  Author:
            //
            //    John Burkardt
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
            //    Input, int ORDER, the order of the polynomial.
            //
            //    Input, double ALPHA, the exponents of (1-X^2).
            //
            //    Output, double *DP2, the value of J'(ORDER)(X).
            //
            //    Output, double *P1, the value of J(ORDER-1)(X).
            //
            //    Input, double C[ORDER], the recursion coefficients.
            //
        {
            double d;
            double eps;
            double p2 = 0;
            int step;
            int step_max = 10;

            eps = typeMethods.r8_epsilon ( );

            for ( step = 1; step <= step_max; step++ )
            {
                gegenbauer_recur ( ref p2, ref dp2, ref p1, x, order, alpha, c );

                d = p2 / ( dp2 );
                x = x - d;

                if ( Math.Abs ( d ) <= eps * ( Math.Abs ( x ) + 1.0 ) )
                {
                    return;
                }
            }
        }
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

        public static void gegenbauer_poly(int n, double alpha, double x, double[] cx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GEGENBAUER_POLY computes the Gegenbauer polynomials C(0:N,ALPHA,X).
            //
            //  Discussion:
            //
            //    The Gegenbauer polynomial can be evaluated in Mathematica with
            //    the command
            //
            //      GegenbauerC[n,m,x]
            //
            //  Differential equation:
            //
            //    (1-X*X) Y'' - (2 ALPHA + 1) X Y' + N (N + 2 ALPHA) Y = 0
            //
            //  Recursion:
            //
            //    C(0,ALPHA,X) = 1,
            //    C(1,ALPHA,X) = 2*ALPHA*X
            //    C(N,ALPHA,X) = ( 2*(N-1+ALPHA)*X * C(N-1,ALPHA,X) - (N-2+2*ALPHA) * C(N-2,ALPHA,X) )/N
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
            //    Integral ( -1 <= X <= 1 )
            //      ( 1 - X^2 )^( ALPHA - 0.5 ) * C(N,ALPHA,X)^2 dX
            //
            //    = PI * 2^( 1 - 2 * ALPHA ) * Gamma ( N + 2 * ALPHA )
            //      / ( N! * ( N + ALPHA ) * ( Gamma ( ALPHA ) )^2 )
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
            //    Stephen Wolfram,
            //    The Mathematica Book,
            //    Fourth Edition,
            //    Cambridge University Press, 1999,
            //    ISBN: 0-521-64314-7,
            //    LC: QA76.95.W65.
            //
            //  Parameters:
            //
            //    Input, int N, the highest order polynomial to compute.
            //    Note that polynomials 0 through N will be computed.
            //
            //    Input, double ALPHA, a parameter which is part of the definition of
            //    the Gegenbauer polynomials.  It must be greater than -0.5.
            //
            //    Input, double X, the point at which the polynomials are to be evaluated.
            //
            //    Output, double CX[N+1], the values of the first N+1 Gegenbauer
            //    polynomials at the point X.
            //
        {
            int i;

            if (alpha <= -0.5)
            {
                Console.WriteLine("");
                Console.WriteLine("GEGENBAUER_POLY - Fatal error!");
                Console.WriteLine("  Illegal value of ALPHA = " + alpha + "");
                Console.WriteLine("  but ALPHA must be greater than -0.5.");
                return;
            }

            if (n < 0)
            {
                return;
            }

            cx[0] = 1.0;

            if (n == 0)
            {
                return;
            }

            cx[1] = 2.0 * alpha * x;

            for (i = 2; i <= n; i++)
            {
                cx[i] = (((double)(2 * i - 2) + 2.0 * alpha) * x * cx[i - 1]
                         + ((double)(-i + 2) - 2.0 * alpha) * cx[i - 2])
                        / (double)i;
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

    }
}