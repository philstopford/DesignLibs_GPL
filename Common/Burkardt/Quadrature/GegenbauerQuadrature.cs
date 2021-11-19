using System;
using Burkardt.PolynomialNS;
using Burkardt.Types;

namespace Burkardt.Quadrature;

public static class GegenbauerQuadrature
{
    public static void gegenbauer_compute_np ( int order, int np, double[] p, ref double[] x,
            ref double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GEGENBAUER_COMPUTE_NP computes a Gegenbauer quadrature rule.
        //
        //  Discussion:
        //
        //    The integral:
        //
        //      Integral ( -1 <= X <= 1 ) (1-X^2)^ALPHA * F(X) dX
        //
        //    The quadrature rule:
        //
        //      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )
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
        //    22 June 2009
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
        //    Input, int ORDER, the order.
        //    1 <= ORDER.
        //
        //    Input, int NP, the number of parameters.
        //
        //    Input, double P[NP], contains parameters.
        //    P[0] = ALPHA = the exponent of (1-X^2).  -1.0 < ALPHA is required.
        //
        //    Output, double X[ORDER], the abscissas.
        //
        //    Output, double W[ORDER], the weights.
        //
    {
        double alpha = p[0];

        gegenbauer_compute ( order, alpha, ref x, ref w );
    }
        
    public static void gegenbauer_compute(int order, double alpha, ref double[] x, ref double[] w)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GEGENBAUER_COMPUTE computes a Gegenbauer quadrature rule.
        //
        //  Discussion:
        //
        //    The integral:
        //
        //      Integral ( -1 <= X <= 1 ) (1-X^2)^ALPHA * F(X) dX
        //
        //    The quadrature rule:
        //
        //      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )
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
        //    13 June 2009
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
        //    Input, int ORDER, the order.
        //    1 <= ORDER.
        //
        //    Input, double ALPHA, the exponent of (1-X^2).  -1.0 < ALPHA is required.
        //
        //    Output, double X[ORDER], the abscissas.
        //
        //    Output, double W[ORDER], the weights.
        //
    {
        double dp2 = 0;
        int i;
        double p1 = 0;
        double temp;
        double x0 = 0;
        switch (order)
        {
            //
            //  Check ORDER.
            //
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("GEGENBAUER_COMPUTE - Fatal error!");
                Console.WriteLine("  1 <= ORDER is required.");
                return;
        }

        double[] c = new double[order];
        switch (alpha)
        {
            //
            //  Check ALPHA.
            //
            case <= -1.0:
                Console.WriteLine("");
                Console.WriteLine("GEGENBAUER_COMPUTE - Fatal error!");
                Console.WriteLine("  -1.0 < ALPHA is required.");
                return;
        }

        //
        //  Set the recursion coefficients.
        //
        c[0] = 0.0;
        c[1] = order switch
        {
            >= 2 => 1.0 / (2.0 * alpha + 3.0),
            _ => c[1]
        };

        for (i = 3; i <= order; i++)
        {
            c[i - 1] = (i - 1)
                       * (alpha + alpha + (i - 1)) /
                       ((alpha + alpha + (2 * i - 1))
                        * (alpha + alpha + (2 * i - 3)));
        }

        double delta = typeMethods.r8_gamma(alpha + 1.0)
                       * typeMethods.r8_gamma(alpha + 1.0)
                       / typeMethods.r8_gamma(alpha + alpha + 2.0);

        double prod = 1.0;
        for (i = 2; i <= order; i++)
        {
            prod *= c[i - 1];
        }

        double cc = delta * Math.Pow(2.0, alpha + alpha + 1.0) * prod;

        for (i = 1; i <= order; i++)
        {
            double r2;
            double r1;
            double r3;
            switch (i)
            {
                case 1:
                    double an = alpha / order;

                    r1 = (1.0 + alpha)
                         * (2.78 / (4.0 + order * order)
                            + 0.768 * an / order);

                    r2 = 1.0 + 2.44 * an + 1.282 * an * an;

                    x0 = (r2 - r1) / r2;
                    break;
                case 2:
                    r1 = (4.1 + alpha) /
                         ((1.0 + alpha) * (1.0 + 0.156 * alpha));

                    r2 = 1.0 + 0.06 * (order - 8.0) *
                        (1.0 + 0.12 * alpha) / order;

                    r3 = 1.0 + 0.012 * alpha *
                        (1.0 + 0.25 * Math.Abs(alpha)) / order;

                    x0 -= r1 * r2 * r3 * (1.0 - x0);
                    break;
                case 3:
                    r1 = (1.67 + 0.28 * alpha) / (1.0 + 0.37 * alpha);

                    r2 = 1.0 + 0.22 * (order - 8.0)
                        / order;

                    r3 = 1.0 + 8.0 * alpha /
                        ((6.28 + alpha) * (order * order));

                    x0 -= r1 * r2 * r3 * (x[0] - x0);
                    break;
                default:
                {
                    if (i < order - 1)
                    {
                        x0 = 3.0 * x[i - 2] - 3.0 * x[i - 3] + x[i - 4];
                    }
                    else if (i == order - 1)
                    {
                        r1 = (1.0 + 0.235 * alpha) / (0.766 + 0.119 * alpha);

                        r2 = 1.0 / (1.0 + 0.639
                            * (order - 4.0)
                            / (1.0 + 0.71 * (order - 4.0)));

                        r3 = 1.0 / (1.0 + 20.0 * alpha / ((7.5 + alpha) *
                                                          (order * order)));

                        x0 += r1 * r2 * r3 * (x0 - x[i - 3]);
                    }
                    else if (i == order)
                    {
                        r1 = (1.0 + 0.37 * alpha) / (1.67 + 0.28 * alpha);

                        r2 = 1.0 /
                             (1.0 + 0.22 * (order - 8.0)
                                 / order);

                        r3 = 1.0 / (1.0 + 8.0 * alpha /
                            ((6.28 + alpha) * (order * order)));

                        x0 += r1 * r2 * r3 * (x0 - x[i - 3]);
                    }

                    break;
                }
            }

            GegenbauerPolynomial.gegenbauer_root(ref x0, order, alpha, ref dp2, ref p1, c);

            x[i - 1] = x0;
            w[i - 1] = cc / (dp2 * p1);
        }

        //
        //  Reverse the order of the values.
        //
        for (i = 1; i <= order / 2; i++)
        {
            temp = x[i - 1];
            x[i - 1] = x[order - i];
            x[order - i] = temp;
        }

        for (i = 1; i <= order / 2; i++)
        {
            temp = w[i - 1];
            w[i - 1] = w[order - i];
            w[order - i] = temp;
        }
    }

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
        double dp2 = 0;
        int i;
        double p1 = 0;
        double temp;
        double x = 0;
        switch (order)
        {
            //
            //  Check ORDER.
            //
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("GEGENBAUER_SS_COMPUTE - Fatal error!");
                Console.WriteLine("  1 <= ORDER is required.");
                return;
        }

        double[] c = new double[order];
        switch (alpha)
        {
            //
            //  Check ALPHA.
            //
            case <= -1.0:
                Console.WriteLine("");
                Console.WriteLine("GEGENBAUER_SS_COMPUTE - Fatal error!");
                Console.WriteLine("  -1.0 < ALPHA is required.");
                return;
        }

        //
        //  Set the recursion coefficients.
        //
        c[0] = 0.0;
        c[1] = order switch
        {
            >= 2 => 1.0 / (2.0 * alpha + 3.0),
            _ => c[1]
        };

        for (i = 3; i <= order; i++)
        {
            c[i - 1] = (i - 1)
                       * (alpha + alpha + (i - 1)) /
                       ((alpha + alpha + (2 * i - 1))
                        * (alpha + alpha + (2 * i - 3)));
        }

        double delta = Helpers.Gamma(alpha + 1.0)
                       * Helpers.Gamma(alpha + 1.0)
                       / Helpers.Gamma(alpha + alpha + 2.0);

        double prod = 1.0;
        for (i = 2; i <= order; i++)
        {
            prod *= c[i - 1];
        }

        double cc = delta * Math.Pow(2.0, alpha + alpha + 1.0) * prod;

        for (i = 1; i <= order; i++)
        {
            double r1;
            double r2;
            double r3;
            switch (i)
            {
                case 1:
                    double an = alpha / order;

                    r1 = (1.0 + alpha)
                         * (2.78 / (4.0 + order * order)
                            + 0.768 * an / order);

                    r2 = 1.0 + 2.44 * an + 1.282 * an * an;

                    x = (r2 - r1) / r2;
                    break;
                case 2:
                    r1 = (4.1 + alpha) /
                         ((1.0 + alpha) * (1.0 + 0.156 * alpha));

                    r2 = 1.0 + 0.06 * (order - 8.0) *
                        (1.0 + 0.12 * alpha) / order;

                    r3 = 1.0 + 0.012 * alpha *
                        (1.0 + 0.25 * Math.Abs(alpha)) / order;

                    x -= r1 * r2 * r3 * (1.0 - x);
                    break;
                case 3:
                    r1 = (1.67 + 0.28 * alpha) / (1.0 + 0.37 * alpha);

                    r2 = 1.0 + 0.22 * (order - 8.0)
                        / order;

                    r3 = 1.0 + 8.0 * alpha /
                        ((6.28 + alpha) * (order * order));

                    x -= r1 * r2 * r3 * (xtab[0] - x);
                    break;
                default:
                {
                    if (i < order - 1)
                    {
                        x = 3.0 * xtab[i - 2] - 3.0 * xtab[i - 3] + xtab[i - 4];
                    }
                    else if (i == order - 1)
                    {
                        r1 = (1.0 + 0.235 * alpha) / (0.766 + 0.119 * alpha);

                        r2 = 1.0 / (1.0 + 0.639
                            * (order - 4.0)
                            / (1.0 + 0.71 * (order - 4.0)));

                        r3 = 1.0 / (1.0 + 20.0 * alpha / ((7.5 + alpha) *
                                                          (order * order)));

                        x += r1 * r2 * r3 * (x - xtab[i - 3]);
                    }
                    else if (i == order)
                    {
                        r1 = (1.0 + 0.37 * alpha) / (1.67 + 0.28 * alpha);

                        r2 = 1.0 /
                             (1.0 + 0.22 * (order - 8.0)
                                 / order);

                        r3 = 1.0 / (1.0 + 8.0 * alpha /
                            ((6.28 + alpha) * (order * order)));

                        x += r1 * r2 * r3 * (x - xtab[i - 3]);
                    }

                    break;
                }
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
        int i;

        p1 = 1.0;
        double dp1 = 0.0;

        p2 = x;
        dp2 = 1.0;

        for (i = 2; i <= order; i++)
        {
            double p0 = p1;
            double dp0 = dp1;

            p1 = p2;
            dp1 = dp2;

            p2 = x * p1 - c[i - 1] * p0;
            dp2 = x * dp1 + p1 - c[i - 1] * dp0;
        }
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
        double p2 = 0;
        int step;
        const int step_max = 10;

        double eps = typeMethods.r8_epsilon();

        for (step = 1; step <= step_max; step++)
        {
            gegenbauer_ss_recur(ref p2, ref dp2, ref p1, x, order, alpha, c);

            double d = p2 / dp2;
            x -= d;

            if (Math.Abs(d) <= eps * (Math.Abs(x) + 1.0))
            {
                return;
            }
        }
    }
}