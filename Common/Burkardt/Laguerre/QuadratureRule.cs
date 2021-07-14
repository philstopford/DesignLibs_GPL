using Burkardt.Types;

namespace Burkardt.Laguerre
{
    public static class QuadratureRule
    {
        public static void laguerre_compute(int order, ref double[] xtab, ref double[] weight,
                double alpha)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LAGUERRE_COMPUTE computes a Gauss-Laguerre quadrature rule.
            //
            //  Discussion:
            //
            //    In the simplest case, ALPHA is 0, and we are approximating the
            //    integral from 0 to +oo of EXP(-X) * F(X).  When this is so,
            //    it is easy to modify the rule to approximate the integral from
            //    A to +oo as well.
            //
            //    If ALPHA is nonzero, then there is no simple way to extend the
            //    rule to approximate the integral from A to +oo.  The simplest
            //    procedures would be to approximate the integral from 0 to A.
            //
            //    The integration interval is [ A, +oo ) or [ 0, +oo ).
            //
            //    The weight function is w(x) = exp ( -x ) or exp ( -x ) * x^alpha.
            //
            //
            //    If the integral to approximate is:
            //
            //        Integral ( A <= X < +oo ) EXP ( - X ) * F(X) dX
            //      or
            //        Integral ( 0 <= X < +oo ) EXP ( - X ) * X^ALPHA * F(X) dX
            //
            //    then the quadrature rule is:
            //
            //      EXP ( - A ) * Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( A+XTAB(I) )
            //    or
            //      sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
            //
            //    If the integral to approximate is:
            //
            //        Integral ( A <= X < +oo ) F(X) dX
            //      or
            //        Integral ( 0 <= X < +oo ) X^ALPHA * F(X) dX
            //
            //    then the quadrature rule is:
            //
            //      EXP ( - A ) * Sum ( 1 <= I <= ORDER ) 
            //        WEIGHT(I) * EXP(A+XTAB(I)) * F ( A+XTAB(I) )
            //    or
            //      sum ( 1 <= I <= ORDER ) WEIGHT(I) * EXP(XTAB(I)) * F ( XTAB(I) )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 May 2006
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
            //    Input, int ORDER, the order of the quadrature rule to be computed.
            //    ORDER must be at least 1.
            //
            //    Output, double XTAB[ORDER], the Gauss-Laguerre abscissas.
            //
            //    Output, double WEIGHT[ORDER], the Gauss-Laguerre weights.
            //
            //    Input, double ALPHA, the exponent of the X factor.
            //    Set ALPHA = 0.0 for the simplest rule.
            //    ALPHA must be nonnegative.
            //
        {
            double[] b;
            double[] c;
            double cc;
            double dp2 = 0;
            int i;
            double p1 = 0;
            double prod;
            double r1;
            double r2;
            double ratio;
            double x = 0;

            b = new double[order];
            c = new double[order];
            //
            //  Set the recursion coefficients.
            //
            for (i = 0; i < order; i++)
            {
                b[i] = (alpha + (double) (2 * i + 1));
            }

            for (i = 0; i < order; i++)
            {
                c[i] = (double) (i) * (alpha + (double) (i));
            }

            prod = 1.0;
            for (i = 1; i < order; i++)
            {
                prod = prod * c[i];
            }

            cc = typeMethods.r8_gamma(alpha + 1.0) * prod;

            for (i = 0; i < order; i++)
            {
                //
                //  Compute an estimate for the root.
                //
                if (i == 0)
                {
                    x = (1.0 + alpha) * (3.0 + 0.92 * alpha) /
                        (1.0 + 2.4 * (double) (order) + 1.8 * alpha);
                }
                else if (i == 1)
                {
                    x = x + (15.0 + 6.25 * alpha) /
                        (1.0 + 0.9 * alpha + 2.5 * (double) (order));
                }
                else
                {
                    r1 = (1.0 + 2.55 * (double) (i - 1))
                         / (1.9 * (double) (i - 1));

                    r2 = 1.26 * (double) (i - 1) * alpha /
                         (1.0 + 3.5 * (double) (i - 1));

                    ratio = (r1 + r2) / (1.0 + 0.3 * alpha);

                    x = x + ratio * (x - xtab[i - 2]);
                }

                //
                //  Use iteration to find the root.
                //
                PolynomialNS.Laguerre.laguerre_root(ref x, order, alpha, ref dp2, ref p1, b, c);
                //
                //  Set the abscissa and weight.
                //
                xtab[i] = x;
                weight[i] = (cc / dp2) / p1;
            }
        }
    }
}