using System;

namespace Burkardt.PolynomialNS
{
    public static class Laguerre
    {
        public static void laguerre_recur(ref double p2, ref double dp2, ref double p1, double x,
            int order, double alpha, double[] b, double[] c )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGUERRE_RECUR finds the value and derivative of a Laguerre polynomial.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 May 2006
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
        //    Output, double *P2, the value of L(ORDER)(X).
        //
        //    Output, double *DP2, the value of L'(ORDER)(X).
        //
        //    Output, double *P1, the value of L(ORDER-1)(X).
        //
        //    Input, double X, the point at which polynomials are evaluated.
        //
        //    Input, int ORDER, the order of the polynomial to be computed.
        //
        //    Input, double ALPHA, the exponent of the X factor in the
        //    integrand.
        //
        //    Input, double B[ORDER], C[ORDER], the recursion coefficients.
        //
        {
            double dp0;
            double dp1;
            int i;
            double p0;

            p1 = 1.0;
            dp1 = 0.0;

            p2 = x - alpha - 1.0;
            dp2 = 1.0;

            for (i = 1; i < order; i++)
            {
                p0 = p1;
                dp0 = dp1;

                p1 = p2;
                dp1 = dp2;

                p2 = (x - b[i]) * (p1) - c[i] * p0;
                dp2 = (x - b[i]) * dp1 + (p1) - c[i] * dp0;
            }
        }

        public static void laguerre_root(ref double x, int order, double alpha, ref double dp2,
            ref double p1, double[] b, double[] c )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGUERRE_ROOT improves an approximate root of a Laguerre polynomial.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 May 2006
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
        //    Input, int ORDER, the order of the polynomial to be computed.
        //
        //    Input, double ALPHA, the exponent of the X factor.
        //
        //    Output, double *DP2, the value of L'(ORDER)(X).
        //
        //    Output, double *P1, the value of L(ORDER-1)(X).
        //
        //    Input, double B[ORDER], C[ORDER], the recursion coefficients.
        //
        {
            double d;
            double eps;
            double p2 = 0;
            int step;
            int step_max = 10;

            eps = double.Epsilon;

            for (step = 1; step <= step_max; step++)
            {
                laguerre_recur(ref p2, ref dp2, ref p1, x, order, alpha, b, c);

                d = p2 / (dp2);
                x = x - d;

                if (Math.Abs(d) <= eps * (Math.Abs(x) + 1.0))
                {
                    break;
                }
            }
        }
    }
}