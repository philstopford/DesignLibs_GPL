using System;
using Burkardt.IntegralNS;
using Burkardt.MatrixNS;
using Burkardt.PolynomialNS;
using Burkardt.Types;

namespace Burkardt.Quadrature;

public static class JacobiQuadrature
{
    public class ParameterData
    {
        //
        //  Two global variables needed to support the "parameter" function.
        //
        public double[] P;
        public int[] NP;

        public double value = 0;
    }
    public static ParameterData jacobi_points(Func<ParameterData, int, int, ParameterData> parameter, ParameterData data, int n, int dim, ref double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    JACOBI_POINTS computes Jacobi quadrature points.
        //
        //  Discussion:
        //
        //    This function assumes the existence of a function:
        //      double parameter ( int dim, int offset )
        //    which can supply the values of the ALPHA and BETA parameters.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 April 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order.
        //
        //    Input, int DIM, the spatial dimension represented by this rule,
        //    in cases where a multidimensional product rule is being formed.
        //
        //    Output, double X[N], the abscissas.
        //
    {
        data = parameter(data, dim, 0);
        double alpha = data.value;
        data = parameter(data, dim, 1);
            
        double beta = data.value;

        jacobi_compute_points(n, alpha, beta, ref x);

        return data;
    }

    public static void jacobi_compute_points ( int order, double alpha, double beta,
            ref double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    JACOBI_COMPUTE_POINTS computes Jacobi quadrature points.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 October 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ORDER, the order.
        //
        //    Input, double ALPHA, BETA, the exponents of the (1-X) and (1+X) factors.
        //
        //    Output, double X[ORDER], the abscissas.
        //
    {
        double[] w = new double[order];

        jacobi_compute ( order, alpha, beta, ref x, ref w );
    }
        
    public static double[] jacobi_compute_points_np ( int order, int np, double[] p, double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    JACOBI_COMPUTE_POINTS_NP computes Jacobi quadrature points.
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
        //  Parameters:
        //
        //    Input, int ORDER, the order.
        //
        //    Input, int NP, the number of parameters.
        //
        //    Input, double P[NP], parameter values.
        //    P[0] = ALPHA, the exponent of (1-X)
        //    P[1] = BETA,  the exponent of (1+X).
        //    -1.0 < ALPHA and -1.0 < BETA are required.
        //
        //    Output, double X[ORDER], the abscissas.
        //
    {
        double alpha = p[0];
        double beta = p[1];

        jacobi_compute_points ( order, alpha, beta, ref x );

        return x;
    }
        
    public static ParameterData jacobi_weights(Func<ParameterData, int, int, ParameterData> parameter, ParameterData data, int n, int dim, ref double[] w)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    JACOBI_WEIGHTS computes Jacobi quadrature weights.
        //
        //  Discussion:
        //
        //    This function assumes the existence of a function:
        //      double parameter ( int dim, int offset )
        //    which can supply the values of the ALPHA and BETA parameters.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 April 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order.
        //
        //    Input, int DIM, the spatial dimension represented by this rule,
        //    in cases where a multidimensional product rule is being formed.
        //
        //    Output, double W[N], the weights.
        //
    {
        data = parameter(data, dim, 0);
        double alpha = data.value;
        data = parameter(data, dim, 1);
        double beta = data.value;

        jacobi_compute_weights(n, alpha, beta, ref w);

        return data;
    }
        
    public static void jacobi_compute_weights ( int order, double alpha, double beta,
            ref double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    JACOBI_COMPUTE_WEIGHTS computes Jacobi quadrature weights.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 October 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ORDER, the order.
        //
        //    Input, double ALPHA, BETA, the exponents of the (1-X) and (1+X) factors.
        //
        //    Output, double W[ORDER], the weights.
        //
    {
        double[] x = new double[order];

        jacobi_compute ( order, alpha, beta, ref x, ref w );
    }
        
    public static double[] jacobi_compute_weights_np ( int order, int np, double[] p, double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    JACOBI_COMPUTE_WEIGHTS_NP computes Jacobi quadrature weights.
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
        //  Parameters:
        //
        //    Input, int ORDER, the order.
        //
        //    Input, int NP, the number of parameters.
        //
        //    Input, double P[NP], parameter values.
        //    P[0] = ALPHA, the exponent of (1-X)
        //    P[1] = BETA,  the exponent of (1+X).
        //    -1.0 < ALPHA and -1.0 < BETA are required.
        //
        //    Output, double W[ORDER], the weights.
        //
    {
        double alpha = p[0];
        double beta = p[1];

        jacobi_compute_weights ( order, alpha, beta, ref w );

        return w;
    }

    public static void jacobi_handle(int order, double alpha, double beta, string output)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    JACOBI_HANDLE computes the requested Gauss-Jacobi rule and outputs it.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 January 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ORDER, the order of the rule.
        //
        //    Input, double ALPHA, BETA, the parameters.
        //
        //    Input, string OUTPUT_FILE, specifies the output.
        //    * "C++'", print as C++ code.
        //    * "F77", print as FORTRAN77 code.
        //    * "F90", print as FORTRAN90 code.
        //    * "MAT", print as MATLAB code.
        //    * file,  write files "file_w.txt", "file_x.txt", "file_r.txt" 
        //      defining weights, abscissas, and region.
        // 
    {
        int i;

        double[] r = new double[2];
        double[] w = new double[order];
        double[] x = new double[order];

        r[0] = -1.0;
        r[1] = +1.0;

        jacobi_compute(order, alpha, beta, ref x, ref w);

        switch (output)
        {
            case "C++":
            {
                Console.WriteLine("//");
                Console.WriteLine("//  Weights W, abscissas X and range R");
                Console.WriteLine("//  for a Gauss-Jacobi quadrature rule");
                Console.WriteLine("//  ORDER = " + order + "");
                Console.WriteLine("//  ALPHA = " + alpha + "");
                Console.WriteLine("//  BETA  = " + beta + "");
                Console.WriteLine("//");
                Console.WriteLine("//  Standard rule:");
                Console.WriteLine("//    Integral ( -1 <= x <= +1 ) (1-x)^ALPHA (1+x)^BETA f(x) dx");
                Console.WriteLine("//    is to be approximated by");
                Console.WriteLine("//    sum ( 1 <= I <= ORDER ) w(i) * f(x(i)).");
                Console.WriteLine("//");

                for (i = 0; i < order; i++)
                {
                    Console.WriteLine("  w[" + i + "] = "
                                      + w[i].ToString("0.################") + ";");
                }

                Console.WriteLine("");
                for (i = 0; i < order; i++)
                {
                    Console.WriteLine("  x[" + i + "] = "
                                      + x[i].ToString("0.################") + ";");
                }

                Console.WriteLine("");
                for (i = 0; i < 2; i++)
                {
                    Console.WriteLine("  r[" + i + "] = " + r[i] + ";");
                }

                break;
            }
            case "F77":
            {
                Console.WriteLine("c");
                Console.WriteLine("c  Weights W, abscissas X and range R");
                Console.WriteLine("c  for a Gauss-Jacobi quadrature rule");
                Console.WriteLine("c  ORDER = " + order + "");
                Console.WriteLine("c  ALPHA = " + alpha + "");
                Console.WriteLine("c  BETA  = " + beta + "");
                Console.WriteLine("c");
                Console.WriteLine("c  Standard rule:");
                Console.WriteLine("c    Integral ( -1 <= x <= +1 ) (1-x)^ALPHA (1+x)^BETA f(x) dx");
                Console.WriteLine("c    is to be approximated by");
                Console.WriteLine("c    sum ( 1 <= I <= ORDER ) w(i) * f(x(i)).");
                Console.WriteLine("c");

                for (i = 0; i < order; i++)
                {
                    Console.WriteLine("      w(" + i + 1 + ") = "
                                      + w[i].ToString("0.################") + "");
                }

                Console.WriteLine("");
                for (i = 0; i < order; i++)
                {
                    Console.WriteLine("      x(" + i + 1 + ") = "
                                      + x[i].ToString("0.################") + "");
                }

                Console.WriteLine("");
                for (i = 0; i < 2; i++)
                {
                    Console.WriteLine("      r(" + i + 1 + ") = " + r[i] + "");
                }

                break;
            }
            case "F90":
            {
                Console.WriteLine("!");
                Console.WriteLine("!  Weights W, abscissas X and range R");
                Console.WriteLine("!  for a Gauss-Jacobi quadrature rule");
                Console.WriteLine("!  ORDER = " + order + "");
                Console.WriteLine("!  ALPHA = " + alpha + "");
                Console.WriteLine("!  BETA  = " + beta + "");
                Console.WriteLine("!");
                Console.WriteLine("!  Standard rule:");
                Console.WriteLine("!    Integral ( -1 <= x <= +1 ) (1-x)^ALPHA (1+x)^BETA f(x) dx");
                Console.WriteLine("!    is to be approximated by");
                Console.WriteLine("!    sum ( 1 <= I <= ORDER ) w(i) * f(x(i)).");
                Console.WriteLine("!");

                for (i = 0; i < order; i++)
                {
                    Console.WriteLine("  w(" + i + 1 + ") = "
                                      + w[i].ToString("0.################") + "");
                }

                Console.WriteLine("");
                for (i = 0; i < order; i++)
                {
                    Console.WriteLine("  x(" + i + 1 + ") = "
                                      + x[i].ToString("0.################") + "");
                }

                Console.WriteLine("");
                for (i = 0; i < 2; i++)
                {
                    Console.WriteLine("  r(" + i + 1 + ") = " + r[i] + "");
                }

                break;
            }
            case "MAT":
            {
                Console.WriteLine("%");
                Console.WriteLine("%  Weights W, abscissas X and range R");
                Console.WriteLine("%  for a Gauss-Jacobi quadrature rule");
                Console.WriteLine("%  ORDER = " + order + "");
                Console.WriteLine("%  ALPHA = " + alpha + "");
                Console.WriteLine("%  BETA  = " + beta + "");
                Console.WriteLine("%");
                Console.WriteLine("%  Standard rule:");
                Console.WriteLine("%    Integral ( -1 <= x <= +1 ) (1-x)^ALPHA (1+x)^BETA f(x) dx");
                Console.WriteLine("%    is to be approximated by");
                Console.WriteLine("%    sum ( 1 <= I <= ORDER ) w(i) * f(x(i)).");
                Console.WriteLine("%");

                for (i = 0; i < order; i++)
                {
                    Console.WriteLine("  w(" + i + 1 + ") = "
                                      + w[i].ToString("0.################") + ";");
                }

                Console.WriteLine("");
                for (i = 0; i < order; i++)
                {
                    Console.WriteLine("  x(" + i + 1 + ") = "
                                      + x[i].ToString("0.################") + ";");
                }

                Console.WriteLine("");
                for (i = 0; i < 2; i++)
                {
                    Console.WriteLine("  r(" + i + 1 + ") = " + r[i] + ";");
                }

                break;
            }
            default:
                string output_w = output + "%s_w.txt";
                string output_x = output + "%s_x.txt";
                string output_r = output + "%s_r.txt";

                Console.WriteLine("");
                Console.WriteLine("  Creating quadrature files.");
                Console.WriteLine("");
                Console.WriteLine("  Root file name is     \"" + output + "\".");
                Console.WriteLine("");
                Console.WriteLine("  Weight file will be   \"" + output_w + "\".");
                Console.WriteLine("  Abscissa file will be \"" + output_x + "\".");
                Console.WriteLine("  Region file will be   \"" + output_r + "\".");

                typeMethods.r8mat_write(output_w, 1, order, w);
                typeMethods.r8mat_write(output_x, 1, order, x);
                typeMethods.r8mat_write(output_r, 1, 2, r);
                break;
        }
    }

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
        int i;

        double ab = alpha + beta;
        double abi = 2.0 + ab;
        //
        //  Define the zero-th moment.
        //
        double zemu = Math.Pow(2.0, ab + 1.0) * Helpers.Gamma(alpha + 1.0)
                                              * Helpers.Gamma(beta + 1.0) / Helpers.Gamma(abi);
        //
        //  Define the Jacobi matrix.
        //
        x[0] = (beta - alpha) / abi;
        for (i = 1; i < n; i++)
        {
            x[i] = 0.0;
        }

        double[] bj = new double[n];

        bj[0] = 4.0 * (1.0 + alpha) * (1.0 + beta)
                / ((abi + 1.0) * abi * abi);
        for (i = 1; i < n; i++)
        {
            bj[i] = 0.0;
        }

        double a2b2 = beta * beta - alpha * alpha;

        for (i = 1; i < n; i++)
        {
            double i_r8 = i + 1;
            abi = 2.0 * i_r8 + ab;
            x[i] = a2b2 / ((abi - 2.0) * abi);
            abi *= abi;
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
            w[i] *= w[i];
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
        int i;
        //
        //  Get the exact value of the integral of the unscaled monomial.
        //
        double exact = Integral.jacobi_integral(expon, alpha, beta);
        //
        //  Evaluate the unweighted monomial at the quadrature points.
        //
        double quad = 0.0;
        for (i = 0; i < order; i++)
        {
            quad += w[i] * Math.Pow(x[i], expon);
        }

        double quad_error = exact switch
        {
            //
            //  Absolute error for cases where exact integral is zero,
            //  Relative error otherwise.
            //
            0.0 => Math.Abs(quad),
            _ => Math.Abs(quad - exact) / Math.Abs(exact)
        };

        return quad_error;
    }
        
    public static void jacobi_compute_np ( int order, int np, double[] p, ref double[] x, ref double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    JACOBI_COMPUTE_NP computes a Jacobi quadrature rule.
        //
        //  Discussion:
        //
        //    The integral:
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
        //    22 June 2009
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
        //    Input, int ORDER, the order.
        //    1 <= ORDER.
        //
        //    Input, int NP, the number of parameters.
        //
        //    Input, double P[NP], parameter values.
        //    P[0] = ALPHA, the exponent of (1-X)
        //    P[1] = BETA,  the exponent of (1+X).
        //    -1.0 < ALPHA and -1.0 < BETA are required.
        //
        //    Output, double X[ORDER], the abscissas.
        //
        //    Output, double W[ORDER], the weights.
        //
    {
        double alpha = p[0];
        double beta = p[1];

        jacobi_compute ( order, alpha, beta, ref x, ref w );
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
        double dp2 = 0;
        int i;
        double p1 = 0;
        double temp;
        double x0 = 0;

        switch (order)
        {
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("JACOBI_COMPUTE - Fatal error!");
                Console.WriteLine("  Illegal value of ORDER = " + order + "");
                return;
        }

        double[] b = new double[order];
        double[] c = new double[order];
        switch (alpha)
        {
            //
            //  Check ALPHA and BETA.
            //
            case <= -1.0:
                Console.WriteLine("");
                Console.WriteLine("JACOBI_COMPUTE - Fatal error!");
                Console.WriteLine("  -1.0 < ALPHA is required.");
                return;
        }

        switch (beta)
        {
            case <= -1.0:
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
                           ((alpha + beta + 2 * i)
                            * (alpha + beta + (2 * i - 2)));
            }

            c[i - 1] = i switch
            {
                1 => 0.0,
                _ => 4.0 * (i - 1) * (alpha + (i - 1)) * (beta + (i - 1)) * (alpha + beta + (i - 1)) /
                     ((alpha + beta + (2 * i - 1)) * Math.Pow(alpha + beta + (2 * i - 2), 2) *
                      (alpha + beta + (2 * i - 3)))
            };
        }

        double delta = typeMethods.r8_gamma(alpha + 1.0)
                       * typeMethods.r8_gamma(beta + 1.0)
                       / typeMethods.r8_gamma(alpha + beta + 2.0);

        double prod = 1.0;
        for (i = 2; i <= order; i++)
        {
            prod *= c[i - 1];
        }

        double cc = delta * Math.Pow(2.0, alpha + beta + 1.0) * prod;

        for (i = 1; i <= order; i++)
        {
            double r2;
            double r3;
            double r1;
            switch (i)
            {
                case 1:
                    double an = alpha / order;
                    double bn = beta / order;

                    r1 = (1.0 + alpha)
                         * (2.78 / (4.0 + order * order)
                            + 0.768 * an / order);

                    r2 = 1.0 + 1.48 * an + 0.96 * bn
                         + 0.452 * an * an + 0.83 * an * bn;

                    x0 = (r2 - r1) / r2;
                    break;
                case 2:
                    r1 = (4.1 + alpha) /
                         ((1.0 + alpha) * (1.0 + 0.156 * alpha));

                    r2 = 1.0 + 0.06 * (order - 8.0) *
                        (1.0 + 0.12 * alpha) / order;

                    r3 = 1.0 + 0.012 * beta *
                        (1.0 + 0.25 * Math.Abs(alpha)) / order;

                    x0 -= r1 * r2 * r3 * (1.0 - x0);
                    break;
                case 3:
                    r1 = (1.67 + 0.28 * alpha) / (1.0 + 0.37 * alpha);

                    r2 = 1.0 + 0.22 * (order - 8.0)
                        / order;

                    r3 = 1.0 + 8.0 * beta /
                        ((6.28 + beta) * (order * order));

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
                        r1 = (1.0 + 0.235 * beta) / (0.766 + 0.119 * beta);

                        r2 = 1.0 / (1.0 + 0.639
                            * (order - 4.0)
                            / (1.0 + 0.71 * (order - 4.0)));

                        r3 = 1.0 / (1.0 + 20.0 * alpha / ((7.5 + alpha) *
                                                          (order * order)));

                        x0 += r1 * r2 * r3 * (x0 - x[i - 3]);
                    }
                    else if (i == order)
                    {
                        r1 = (1.0 + 0.37 * beta) / (1.67 + 0.28 * beta);

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