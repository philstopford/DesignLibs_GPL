using System;

namespace Burkardt.Quadrature
{
    public static class GegenbauerQuadrature
    {
        public static void gegenbauer_ss_compute(int order, double alpha, ref double[] xtab,
        ref double[] weight )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GEGENBAUER_SS_COMPUTE computes a Gauss-Gegenbauer quadrature rule.
        //
        //  Discussion:
        //
        //    The integral:
        //
        //      Integral ( -1 <= X <= 1 ) (1-X^2)^ALPHA * F(X) dX
        //
        //    The quadrature rule:
        //
        //      Sum ( 1 <= I <= ORDER ) WEIGHT(I) * F ( XTAB(I) )
        //
        //    Thanks to Janiki Raman for pointing out a problem in an earlier
        //    version of the code that occurred when ALPHA was -0.5.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    24 June 2008
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
        //    Input, int ORDER, the order of the quadrature rule.
        //
        //    Input, double ALPHA, the exponent of (1-X^2) in the weight.  
        //    -1.0 < ALPHA is required.
        //
        //    Output, double XTAB[ORDER], the abscissas.
        //
        //    Output, double WEIGHT[ORDER], the weights.
        //
        {
            double an;
            double[] c;
            double cc;
            double delta;
            double dp2 = 0;
            int i;
            double p1 = 0;
            double prod;
            double r1;
            double r2;
            double r3;
            double temp;
            double x = 0;
            //
            //  Check ORDER.
            //
            if (order < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("GEGENBAUER_SS_COMPUTE - Fatal error!");
                Console.WriteLine("  1 <= ORDER is required.");
                return;
            }

            c = new double[order];
            //
            //  Check ALPHA.
            //
            if (alpha <= -1.0)
            {
                Console.WriteLine("");
                Console.WriteLine("GEGENBAUER_SS_COMPUTE - Fatal error!");
                Console.WriteLine("  -1.0 < ALPHA is required.");
                return;
            }

            //
            //  Set the recursion coefficients.
            //
            c[0] = 0.0;
            if (2 <= order)
            {
                c[1] = 1.0 / (2.0 * alpha + 3.0);
            }

            for (i = 3; i <= order; i++)
            {
                c[i - 1] = (double) (i - 1)
                           * (alpha + alpha + (double) (i - 1)) /
                           ((alpha + alpha + (double) (2 * i - 1))
                            * (alpha + alpha + (double) (2 * i - 3)));
            }

            delta = Helpers.Gamma(alpha + 1.0)
                    * Helpers.Gamma(alpha + 1.0)
                    / Helpers.Gamma(alpha + alpha + 2.0);

            prod = 1.0;
            for (i = 2; i <= order; i++)
            {
                prod = prod * c[i - 1];
            }

            cc = delta * Math.Pow(2.0, alpha + alpha + 1.0) * prod;

            for (i = 1; i <= order; i++)
            {
                if (i == 1)
                {
                    an = alpha / (double) (order);

                    r1 = (1.0 + alpha)
                         * (2.78 / (4.0 + (double) (order * order))
                            + 0.768 * an / (double) (order));

                    r2 = 1.0 + 2.44 * an + 1.282 * an * an;

                    x = (r2 - r1) / r2;
                }
                else if (i == 2)
                {
                    r1 = (4.1 + alpha) /
                         ((1.0 + alpha) * (1.0 + 0.156 * alpha));

                    r2 = 1.0 + 0.06 * ((double) (order) - 8.0) *
                        (1.0 + 0.12 * alpha) / (double) (order);

                    r3 = 1.0 + 0.012 * alpha *
                        (1.0 + 0.25 * Math.Abs(alpha)) / (double) (order);

                    x = x - r1 * r2 * r3 * (1.0 - x);
                }
                else if (i == 3)
                {
                    r1 = (1.67 + 0.28 * alpha) / (1.0 + 0.37 * alpha);

                    r2 = 1.0 + 0.22 * ((double) (order) - 8.0)
                        / (double) (order);

                    r3 = 1.0 + 8.0 * alpha /
                        ((6.28 + alpha) * (double) (order * order));

                    x = x - r1 * r2 * r3 * (xtab[0] - x);
                }
                else if (i < order - 1)
                {
                    x = 3.0 * xtab[i - 2] - 3.0 * xtab[i - 3] + xtab[i - 4];
                }
                else if (i == order - 1)
                {
                    r1 = (1.0 + 0.235 * alpha) / (0.766 + 0.119 * alpha);

                    r2 = 1.0 / (1.0 + 0.639
                        * ((double) (order) - 4.0)
                        / (1.0 + 0.71 * ((double) (order) - 4.0)));

                    r3 = 1.0 / (1.0 + 20.0 * alpha / ((7.5 + alpha) *
                                                      (double) (order * order)));

                    x = x + r1 * r2 * r3 * (x - xtab[i - 3]);
                }
                else if (i == order)
                {
                    r1 = (1.0 + 0.37 * alpha) / (1.67 + 0.28 * alpha);

                    r2 = 1.0 /
                         (1.0 + 0.22 * ((double) (order) - 8.0)
                             / (double) (order));

                    r3 = 1.0 / (1.0 + 8.0 * alpha /
                        ((6.28 + alpha) * (double) (order * order)));

                    x = x + r1 * r2 * r3 * (x - xtab[i - 3]);
                }

                gegenbauer_ss_root(ref x, order, alpha, ref dp2, ref p1, c);

                xtab[i - 1] = x;
                weight[i - 1] = cc / (dp2 * p1);
            }

            //
            //  Reverse the order of the values.
            //
            for (i = 1; i <= order / 2; i++)
            {
                temp = xtab[i - 1];
                xtab[i - 1] = xtab[order - i];
                xtab[order - i] = temp;
            }

            for (i = 1; i <= order / 2; i++)
            {
                temp = weight[i - 1];
                weight[i - 1] = weight[order - i];
                weight[order - i] = temp;
            }
        }

        public static void gegenbauer_ss_recur(ref double p2, ref double dp2, ref double p1, double x,
        int order, double alpha, double[] c )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GEGENBAUER_SS_RECUR: value and derivative of a Gegenbauer polynomial.
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
        //    Output, double &P2, the value of J(ORDER)(X).
        //
        //    Output, double &DP2, the value of J'(ORDER)(X).
        //
        //    Output, double &P1, the value of J(ORDER-1)(X).
        //
        //    Input, double X, the point at which polynomials are evaluated.
        //
        //    Input, int ORDER, the order of the polynomial to be computed.
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

            for (i = 2; i <= order; i++)
            {
                p0 = p1;
                dp0 = dp1;

                p1 = p2;
                dp1 = dp2;

                p2 = x * p1 - c[i - 1] * p0;
                dp2 = x * dp1 + p1 - c[i - 1] * dp0;
            }

            return;
        }

        public static void gegenbauer_ss_root(ref double x, int order, double alpha, ref double dp2,
        ref double p1, double[] c )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GEGENBAUER_SS_ROOT improves a root of a Gegenbauer polynomial.
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
        //    Input/output, double &X, the approximate root, which
        //    should be improved on output.
        //
        //    Input, int ORDER, the order of the polynomial to be computed.
        //
        //    Input, double ALPHA, the exponents of (1-X^2).
        //
        //    Output, double &DP2, the value of J'(ORDER)(X).
        //
        //    Output, double &P1, the value of J(ORDER-1)(X).
        //
        //    Input, double C[ORDER], the recursion coefficients.
        //
        {
            double d;
            double eps;
            double p2 = 0;
            int step;
            int step_max = 10;

            eps = typeMethods.r8_epsilon();

            for (step = 1; step <= step_max; step++)
            {
                gegenbauer_ss_recur(ref p2, ref dp2, ref p1, x, order, alpha, c);

                d = p2 / dp2;
                x = x - d;

                if (Math.Abs(d) <= eps * (Math.Abs(x) + 1.0))
                {
                    return;
                }
            }
        }
    }
}