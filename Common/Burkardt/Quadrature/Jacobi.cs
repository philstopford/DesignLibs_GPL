using System;
using Burkardt.IntegralNS;
using Burkardt.MatrixNS;
using Burkardt.PolynomialNS;
using Burkardt.Types;

namespace Burkardt.Quadrature
{
    public static class JacobiQuadrature
    {
        public static void j_quadrature_rule(int n, double alpha, double beta, ref double[] x,
        ref double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    J_QUADRATURE_RULE: Gauss-Jacobi quadrature based on J(n,a,b,x).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 April 2012
        //
        //  Author:
        //
        //    John Burkardt.
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
        //    Input, int, N, the order.
        //
        //    Input, double, ALPHA, BETA, the parameters.
        //    -1 < ALPHA, BETA.
        //
        //    Output, double X[N], the abscissas.
        //
        //    Output, double W[N], the weights.
        //
        {
            double a2b2;
            double ab;
            double abi;
            double[] bj;
            int i;
            double i_r8;
            double zemu;

            ab = alpha + beta;
            abi = 2.0 + ab;
            //
            //  Define the zero-th moment.
            //
            zemu = Math.Pow(2.0, ab + 1.0) * Helpers.Gamma(alpha + 1.0)
                                      * Helpers.Gamma(beta + 1.0) / Helpers.Gamma(abi);
            //
            //  Define the Jacobi matrix.
            //
            x[0] = (beta - alpha) / abi;
            for (i = 1; i < n; i++)
            {
                x[i] = 0.0;
            }

            bj = new double[n];

            bj[0] = 4.0 * (1.0 + alpha) * (1.0 + beta)
                    / ((abi + 1.0) * abi * abi);
            for (i = 1; i < n; i++)
            {
                bj[i] = 0.0;
            }

            a2b2 = beta * beta - alpha * alpha;

            for (i = 1; i < n; i++)
            {
                i_r8 = (double) (i + 1);
                abi = 2.0 * i_r8 + ab;
                x[i] = a2b2 / ((abi - 2.0) * abi);
                abi = abi * abi;
                bj[i] = 4.0 * i_r8 * (i_r8 + alpha) * (i_r8 + beta)
                    * (i_r8 + ab) / ((abi - 1.0) * abi);
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
                w[i] = w[i] * w[i];
            }
        }

                public static double monomial_quadrature_jacobi(int expon, double alpha, double beta,
                int order, double[] w, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MONOMIAL_QUADRATURE_JACOBI applies a quadrature rule to a monomial.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    22 January 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int EXPON, the exponent.
            //
            //    Input, double ALPHA, the exponent of (1-X) in the weight factor.
            //
            //    Input, double BETA, the exponent of (1+X) in the weight factor.
            //
            //    Input, int ORDER, the number of points in the rule.
            //
            //    Input, double W[ORDER], the quadrature weights.
            //
            //    Input, double X[ORDER], the quadrature points.
            //
            //    Output, double MONOMIAL_QUADRATURE_JACOBI, the quadrature error.
            //
        {
            double exact;
            int i;
            double quad;
            double quad_error;
            //
            //  Get the exact value of the integral of the unscaled monomial.
            //
            exact = Integral.jacobi_integral(expon, alpha, beta);
            //
            //  Evaluate the unweighted monomial at the quadrature points.
            //
            quad = 0.0;
            for (i = 0; i < order; i++)
            {
                quad = quad + w[i] * Math.Pow(x[i], expon);
            }

            //
            //  Absolute error for cases where exact integral is zero,
            //  Relative error otherwise.
            //
            if (exact == 0.0)
            {
                quad_error = Math.Abs(quad);
            }
            else
            {
                quad_error = Math.Abs(quad - exact) / Math.Abs(exact);
            }

            return quad_error;
        }

        public static void jacobi_compute(int order, double alpha, double beta, ref double[] x,
                ref double[] w)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    JACOBI_COMPUTE computes a Jacobi quadrature rule.
            //
            //  Discussion:
            //
            //    The integration interval is [ -1, 1 ].
            //
            //    The weight function is w(x) = (1-X)^ALPHA * (1+X)^BETA.
            //
            //    The integral to approximate:
            //
            //      Integral ( -1 <= X <= 1 ) (1-X)^ALPHA * (1+X)^BETA * F(X) dX
            //
            //    The quadrature rule:
            //
            //      Sum ( 1 <= I <= ORDER ) W(I) * F ( X(I) )
            //
            //    Thanks to Xu Xiang of Fudan University for pointing out that
            //    an earlier implementation of this routine was incorrect!
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    18 February 2008
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
            //    Input, int ORDER, the order of the rule.
            //    1 <= ORDER.
            //
            //    Input, double ALPHA, BETA, the exponents of (1-X) and
            //    (1+X) in the quadrature rule.  For simple Legendre quadrature,
            //    set ALPHA = BETA = 0.0.  -1.0 < ALPHA and -1.0 < BETA are required.
            //
            //    Output, double X[ORDER], the abscissas.
            //
            //    Output, double W[ORDER], the weights.
            //
        {
            double an;
            double[] b;
            double bn;
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
            double x0 = 0;

            if (order < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("JACOBI_COMPUTE - Fatal error!");
                Console.WriteLine("  Illegal value of ORDER = " + order + "");
                return;
            }

            b = new double[order];
            c = new double[order];
            //
            //  Check ALPHA and BETA.
            //
            if (alpha <= -1.0)
            {
                Console.WriteLine("");
                Console.WriteLine("JACOBI_COMPUTE - Fatal error!");
                Console.WriteLine("  -1.0 < ALPHA is required.");
                return;
            }

            if (beta <= -1.0)
            {
                Console.WriteLine("");
                Console.WriteLine("JACOBI_COMPUTE - Fatal error!");
                Console.WriteLine("  -1.0 < BETA is required.");
                return;
            }

            //
            //  Set the recursion coefficients.
            //
            for (i = 1; i <= order; i++)
            {
                if (alpha + beta == 0.0 || beta - alpha == 0.0)
                {
                    b[i - 1] = 0.0;
                }
                else
                {
                    b[i - 1] = (alpha + beta) * (beta - alpha) /
                               ((alpha + beta + (double)(2 * i))
                                * (alpha + beta + (double)(2 * i - 2)));
                }

                if (i == 1)
                {
                    c[i - 1] = 0.0;
                }
                else
                {
                    c[i - 1] = 4.0 * (double)(i - 1)
                                   * (alpha + (double)(i - 1))
                                   * (beta + (double)(i - 1))
                                   * (alpha + beta + (double)(i - 1)) /
                               ((alpha + beta + (double)(2 * i - 1))
                                * Math.Pow(alpha + beta + (double)(2 * i - 2), 2)
                                * (alpha + beta + (double)(2 * i - 3)));
                }
            }

            delta = typeMethods.r8_gamma(alpha + 1.0)
                    * typeMethods.r8_gamma(beta + 1.0)
                    / typeMethods.r8_gamma(alpha + beta + 2.0);

            prod = 1.0;
            for (i = 2; i <= order; i++)
            {
                prod = prod * c[i - 1];
            }

            cc = delta * Math.Pow(2.0, alpha + beta + 1.0) * prod;

            for (i = 1; i <= order; i++)
            {
                if (i == 1)
                {
                    an = alpha / (double)(order);
                    bn = beta / (double)(order);

                    r1 = (1.0 + alpha)
                         * (2.78 / (4.0 + (double)(order * order))
                            + 0.768 * an / (double)(order));

                    r2 = 1.0 + 1.48 * an + 0.96 * bn
                         + 0.452 * an * an + 0.83 * an * bn;

                    x0 = (r2 - r1) / r2;
                }
                else if (i == 2)
                {
                    r1 = (4.1 + alpha) /
                         ((1.0 + alpha) * (1.0 + 0.156 * alpha));

                    r2 = 1.0 + 0.06 * ((double)(order) - 8.0) *
                        (1.0 + 0.12 * alpha) / (double)(order);

                    r3 = 1.0 + 0.012 * beta *
                        (1.0 + 0.25 * Math.Abs(alpha)) / (double)(order);

                    x0 = x0 - r1 * r2 * r3 * (1.0 - x0);
                }
                else if (i == 3)
                {
                    r1 = (1.67 + 0.28 * alpha) / (1.0 + 0.37 * alpha);

                    r2 = 1.0 + 0.22 * ((double)(order) - 8.0)
                        / (double)(order);

                    r3 = 1.0 + 8.0 * beta /
                        ((6.28 + beta) * (double)(order * order));

                    x0 = x0 - r1 * r2 * r3 * (x[0] - x0);
                }
                else if (i < order - 1)
                {
                    x0 = 3.0 * x[i - 2] - 3.0 * x[i - 3] + x[i - 4];
                }
                else if (i == order - 1)
                {
                    r1 = (1.0 + 0.235 * beta) / (0.766 + 0.119 * beta);

                    r2 = 1.0 / (1.0 + 0.639
                        * ((double)(order) - 4.0)
                        / (1.0 + 0.71 * ((double)(order) - 4.0)));

                    r3 = 1.0 / (1.0 + 20.0 * alpha / ((7.5 + alpha) *
                                                      (double)(order * order)));

                    x0 = x0 + r1 * r2 * r3 * (x0 - x[i - 3]);
                }
                else if (i == order)
                {
                    r1 = (1.0 + 0.37 * beta) / (1.67 + 0.28 * beta);

                    r2 = 1.0 /
                         (1.0 + 0.22 * ((double)(order) - 8.0)
                             / (double)(order));

                    r3 = 1.0 / (1.0 + 8.0 * alpha /
                        ((6.28 + alpha) * (double)(order * order)));

                    x0 = x0 + r1 * r2 * r3 * (x0 - x[i - 3]);
                }

                Jacobi.jacobi_root(ref x0, order, alpha, beta, ref dp2, ref p1, b, c);

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
    }
}