using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.MatrixNS;
using Burkardt.Types;

namespace Burkardt.ChebyshevPolynomialNS
{
    public static partial class ChebyshevPolynomial
    {
        public static double[] v_mass_matrix(int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    V_MASS_MATRIX computes the mass matrix for the Chebyshev V polynomial.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    20 July 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] a;
            int i;
            int j;
            double[] phi;
            double[] phiw;
            double[] w;
            double[] x;

            x = new double[n + 1];
            w = new double[n + 1];

            v_quadrature_rule(n + 1, ref x, ref w);

            phi = v_polynomial(n + 1, n, x);

            phiw = new double[(n + 1) * (n + 1)];

            for (i = 0; i <= n; i++)
            {
                for (j = 0; j <= n; j++)
                {
                    phiw[j + i * (n + 1)] = w[i] * phi[i + j * (n + 1)];
                }
            }

            a = typeMethods.r8mat_mm_new(n + 1, n + 1, n + 1, phiw, phi);

            return a;
        }
        //****************************************************************************80

        public static double v_moment(int e)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    V_MOMENT: integral ( -1 <= x <= +1 ) x^e sqrt(1+x) / sqrt(1-x) dx.
            //
            //  Discussion:
            //
            //     E    V_MOMENT
            //    --    --------------
            //     0      pi
            //     1      Math.PI / 2
            //     2      Math.PI / 2
            //     3    3 Math.PI / 8
            //     4    3 Math.PI / 8
            //     5    5 Math.PI / 16
            //     6    5 Math.PI / 16
            //     7   35 Math.PI / 128
            //     8   35 Math.PI / 128
            //     9   63 Math.PI / 256
            //    10   63 Math.PI / 256
            //    11  231 Math.PI / 1024
            //    12  231 Math.PI / 1024
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    18 July 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int E, the exponent of X.
            //    0 <= E.
            //
            //    Output, double V_MOMENT, the value of the integral.
            //
        {
            double f1;
            double f2;
            double f3;
            double f4;
            double f5;
            double f6;
            double f7;
            double f8;
            double r8_e;
            
            double value;

            r8_e = (double) (e);

            f1 = 1.0 / Helpers.Gamma(1.5 + r8_e);
            f2 = typeMethods.r8_mop(e);
            f3 = Math.PI * Helpers.Gamma(1.5 + r8_e);
            f4 = 2.0 * typeMethods.r8_hyper_2f1(0.5, -r8_e, 1.0, 2.0);
            f5 = (-1.0 + typeMethods.r8_mop(e)) * typeMethods.r8_hyper_2f1(0.5, -r8_e, 2.0, 2.0);
            f6 = Math.Sqrt(Math.PI) * typeMethods.r8_factorial(e);
            f7 = (-1.0 + typeMethods.r8_mop(e))
                 * typeMethods.r8_hyper_2f1(-0.5, 1.0 + r8_e, 1.5 + r8_e, -1.0);
            f8 = 2.0 * typeMethods.r8_hyper_2f1(0.5, 1.0 + r8_e, 1.5 + r8_e, -1.0);

            value = f1 * f2 * (f3 * (f4 + f5) - f6 * (f7 + f8));

            return value;
        }

        public static double[] v_polynomial(int m, int n, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    V_POLYNOMIAL evaluates Chebyshev polynomials V(n,x).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    24 April 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the number of evaluation points.
            //
            //    Input, int N, the highest polynomial to compute.
            //
            //    Input, double X[M], the evaluation points.
            //
            //    Output, double V_POLYNOMIAL[M*(N+1)], the values of the N+1 Chebyshev
            //    polynomials.
            //
        {
            int i;
            int j;
            double[] v;

            if (n < 0)
            {
                return null;
            }

            v = new double[m * (n + 1)];

            for (i = 0; i < m; i++)
            {
                v[i + 0 * m] = 1.0;
            }

            if (n < 1)
            {
                return v;
            }

            for (i = 0; i < m; i++)
            {
                v[i + 1 * m] = 2.0 * x[i] - 1.0;
            }

            for (i = 0; i < m; i++)
            {
                for (j = 2; j <= n; j++)
                {
                    v[i + j * m] = 2.0 * x[i] * v[i + (j - 1) * m] - v[i + (j - 2) * m];
                }
            }

            return v;
        }

        public static void v_polynomial_01_values(ref int n_data, ref int n, ref double x, ref double fx )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    V_POLYNOMIAL_01_VALUES: values of shifted Chebyshev polynomials V01(n,x).
        //
        //  Discussion:
        //
        //    V01(n,x) = V(n,2*x-1)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 July 2015
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
            int N_MAX = 25;

            double[] fx_vec =
            {
                0.0000000000000000,
                1.0000000000000000,
                0.4000000000000000,
                -0.4400000000000000,
                -1.0160000000000000,
                -0.9824000000000000,
                -0.3593600000000000,
                0.4792960000000000,
                1.0303744000000000,
                0.9632281600000000,
                0.3181450240000000,
                -0.5178251264000000,
                -1.0431002009600000,
                -0.9425151549440000,
                -15.000000000000000,
                3.1417984000000000,
                -1.3912448000000000,
                -1.2177792000000000,
                1.1837056000000000,
                1.0000000000000000,
                -0.8558976000000000,
                -0.8905088000000000,
                0.8752768000000000,
                0.1197696000000000,
                1.0000000000000000
            }
            ;

            int[] n_vec =
            {
                -1,
                0, 1, 2,
                3, 4, 5,
                6, 7, 8,
                9, 10, 11,
                12, 7, 7,
                7, 7, 7,
                7, 7, 7,
                7, 7, 7
            }
            ;

            double[] x_vec =
            {
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.85,
                0.00,
                0.10,
                0.20,
                0.30,
                0.40,
                0.50,
                0.60,
                0.70,
                0.80,
                0.90,
                1.00
            }
            ;

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

        public static double[] v_polynomial_ab(double a, double b, int m, int n, double[] xab)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    V_POLYNOMIAL_AB: Chebyshev polynomials VAB(n,x) in [A,B].
            //
            //  Discussion:
            //
            //    VAB(n,x) = V(n,(2*x-a-b)/(b-a))
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    20 July 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double A, B, the domain of definition.
            //
            //    Input, int M, the number of evaluation points.
            //
            //    Input, int N, the highest polynomial to compute.
            //
            //    Input, double XAB[M], the evaluation points.
            //    A <= XAB(*) <= B.
            //
            //    Output, double V_POLYNOMIAL_AB[M*(N+1)], the values.
            //
        {
            int i;
            double[] v;
            double[] x;

            x = new double[m];

            for (i = 0; i < m; i++)
            {
                x[i] = (2.0 * xab[i] - a - b) / (b - a);
            }

            v = v_polynomial(m, n, x);

            return v;
        }

        public static double v_polynomial_ab_value(double a, double b, int n, double xab)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    V_POLYNOMIAL_AB_VALUE: Chebyshev polynomial VAB(n,x) in [A,B].
            //
            //  Discussion:
            //
            //    VAB(n,x) = V(n,(2*x-a-b)/(b-a))
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    20 July 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double A, B, the domain of definition.
            //
            //    Input, int N, the order of the polynomial.
            //
            //    Input, double XAB, the evaluation point.
            //    A <= XAB(*) <= B.
            //
            //    Output, double V_POLYNOMIAL_AB_VALUE, the value.
            //
        {
            double v;
            double x;

            x = (2.0 * xab - a - b) / (b - a);

            v = v_polynomial_value(n, x);

            return v;
        }

        public static double[] v_polynomial_coefficients(int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    V_POLYNOMIAL_COEFFICIENTS: coefficients of the Chebyshev polynomial V(n,x).
            //
            //  First terms:
            //
            //    N/K     0     1      2      3       4     5      6    7      8    9   10
            //
            //     0      1
            //     1     -1     2
            //     2     -1    -2      4
            //     3      1    -4     -4      8
            //     4      1    +4    -12     -8      16
            //     5     -1     6    +12    -32     -16    32
            //     6     -1    -6     24    +32     -80   -32     64
            //     7     +1    -8    -24     80     +80  -192    -64   128
            //
            //  Recursion:
            //
            //    V(0,X) = 1,
            //    V(1,X) = 2 * X - 1,
            //    V(N,X) = 2 * X * V(N-1,X) - V(N-2,X)
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    16 July 2015
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
            //    Output, double V_POLYNOMIAL_COEFFICIENTS[(N+1)*(N+1)], the coefficients.
            //
        {
            double[] c;
            int i;
            int j;

            if (n < 0)
            {
                return null;
            }

            c = new double[(n + 1) * (n + 1)];

            for (i = 0; i <= n; i++)
            {
                for (j = 0; j <= n; j++)
                {
                    c[i + j * (n + 1)] = 0.0;
                }
            }

            c[0 + 0 * (n + 1)] = 1.0;

            if (n == 0)
            {
                return c;
            }

            c[1 + 0 * (n + 1)] = -1.0;
            c[1 + 1 * (n + 1)] = 2.0;

            for (i = 2; i <= n; i++)
            {
                c[i + 0 * (n + 1)] = -c[i - 2 + 0 * (n + 1)];
                for (j = 1; j <= i - 2; j++)
                {
                    c[i + j * (n + 1)] = 2.0 * c[i - 1 + (j - 1) * (n + 1)] - c[i - 2 + j * (n + 1)];
                }

                c[i + (i - 1) * (n + 1)] = 2.0 * c[i - 1 + (i - 2) * (n + 1)];
                c[i + i * (n + 1)] = 2.0 * c[i - 1 + (i - 1) * (n + 1)];
            }

            return c;
        }

        public static void v_polynomial_plot(int n_num, int[] n_val, string output_filename )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    V_POLYNOMIAL_PLOT plots Chebyshev polynomials V(n,x).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N_NUM, the number of polynomials to be plotted.
        //
        //    Input, int N_VAL[N_NUM], the degrees of 1 or more 
        //    Chebyshev polynomials to be plotted together.
        //
        //    Input, string OUTPUT_FILENAME, the name into which the 
        //    graphics information is to be stored.  Note that the PNG format will 
        //    be used.
        //
        {
            double a;
            double b;
            int column;
            string command_filename;
            List<string> command_unit = new List<string>();
            string data_filename;
            List<string> data_unit = new List<string>();
            int i;
            int j;
            int m = 501;
            int n;
            int n_max;
            double[] v;
            double[] x;

            a = -1.0;
            b = +1.0;

            x = typeMethods.r8vec_linspace_new(m, a, b);
            //
            //  Compute all the data.
            //
            n_max = typeMethods.i4vec_max(n_num, n_val);
            v = v_polynomial(m, n_max, x);
            //
            //  Create the data file.
            //
            data_filename = "v_polynomial_data.txt";
            for (i = 0; i < m; i++)
            {
                string line = x[i].ToString();
                for (j = 0; j < n_num; j++)
                {
                    n = n_val[j];
                    line += "  " + v[i + n * m];
                }

                data_unit.Add(line);
            }

            File.WriteAllLines(data_filename, data_unit);
            Console.WriteLine("");
            Console.WriteLine("  Created graphics data file '" + data_filename + "'.");
            //
            //  Plot the selected data.
            //
            command_filename = "v_polynomial_commands.txt";

            command_unit.Add("# " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("# Usage:");
            command_unit.Add("#  gnuplot < " + command_filename + "");
            command_unit.Add("#");
            command_unit.Add("set term png");
            command_unit.Add("set nokey");
            command_unit.Add("set output '" + output_filename + "'");
            command_unit.Add("set xlabel '<---X--->'");
            command_unit.Add("set ylabel '<---V(n,x)--->'");
            command_unit.Add("set title 'Chebyshev Polynomials V(n,x)'");
            command_unit.Add("set grid");
            command_unit.Add("set style data lines");
            for (j = 0; j < n_num; j++)
            {
                string line = "";
                column = n_val[j] + 1;
                if (j == 0)
                {
                    line += "plot ";
                }
                else
                {
                    line += "     ";
                }

                line += "'" + data_filename + "' using 1:" + column + " lw 3 linecolor rgb 'red'";
                if (j < n_num - 1)
                {
                    line += ", \\";
                }
                else
                {
                    line += "";
                }
                command_unit.Add(line);
            }

            File.WriteAllLines(command_filename, command_unit);
            Console.WriteLine("  Created graphics command file '" + command_filename + "'.");
        }

        public static double v_polynomial_value(int n, double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    V_POLYNOMIAL_VALUE: returns the single value V(n,x).
            //
            //  Discussion:
            //
            //    In cases where calling V_POLYNOMIAL is inconvenient, because it returns
            //    a vector of values for multiple arguments X, this simpler interface
            //    may be appropriate.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    11 July 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of the polynomial.
            //
            //    Input, double X, the argument of the polynomial.
            //
            //    Output, double V_POLYNOMIAL_VALUE, the value of V(n,x).
            //
        {
            int m;
            double[] v_vec;
            double value;
            double[] x_vec = new double[1];

            if (n < 0)
            {
                value = 0.0;
            }
            else
            {
                m = 1;
                x_vec[0] = x;

                v_vec = v_polynomial(m, n, x_vec);

                value = v_vec[n];

            }

            return value;
        }

        public static void v_polynomial_values(ref int n_data, ref int n, ref double x, ref double fx )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    V_POLYNOMIAL_VALUES returns values of Chebyshev polynomials V(n,x).
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      u = Sqrt[(x+1)/2],
        //      ChebyshevT[2*n+1,u] / u
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 July 2015
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
            int N_MAX = 14;

            double[] fx_vec =
            {
                0.0000000000000000E+00,
                1.0000000000000000E+00,
                0.6000000000000000E+00,
                -0.0400000000000000E+00,
                -0.6640000000000000E+00,
                -1.0224000000000000E+00,
                -0.9718400000000000E+00,
                -0.5325440000000000E+00,
                0.1197696000000000E+00,
                0.7241753600000000E+00,
                1.0389109760000000E+00,
                0.9380822016000000E+00,
                0.4620205465600000E+00,
                -0.1988493271040000E+00
            }
            ;

            int[] n_vec =
            {
                -1,
                0, 1, 2,
                3, 4, 5,
                6, 7, 8,
                9, 10, 11,
                12
            }
            ;

            double[] x_vec =
            {
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00,
                0.8E+00
            }
            ;

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

        public static double[] v_polynomial_zeros(int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    V_POLYNOMIAL_ZEROS returns zeroes of the Chebyshev polynomial V(n,x).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    26 April 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of the polynomial.
            //
            //    Output, double V_POLYNOMIAL_ZEROS[N], the zeroes.
            //
        {
            double angle;
            int i;
            
            double[] z;

            z = new double[n];

            for (i = 0; i < n; i++)
            {
                angle = (double) (2 * n - 2 * i - 1) * Math.PI / (double) (2 * n + 1);
                z[i] = Math.Cos(angle);
            }

            return z;
        }

        public static void v_quadrature_rule(int n, ref double[] t, ref double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    V_QUADRATURE_RULE: quadrature rule for V(n,x).
        //
        //  Discussion:
        //
        //    Thanks to Csaba Kertesz for pointing out some temporary memory that needed to be freed,
        //    06 November 2015.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 November 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the rule.
        //
        //    Output, double T[N], W[N], the points and weights of the rule.
        //
        {
            double[] bj;
            int i;
            

            for (i = 0; i < n; i++)
            {
                t[i] = 0.0;
            }

            t[0] = +0.5;

            bj = new double[n];
            for (i = 0; i < n; i++)
            {
                bj[i] = 0.5;
            }

            w[0] = Math.Sqrt(Math.PI);
            for (i = 1; i < n; i++)
            {
                w[i] = 0.0;
            }

            IMTQLX.imtqlx(n, ref t, ref bj, ref w);

            for (i = 0; i < n; i++)
            {
                w[i] = w[i] * w[i];
            }
        }

        public static double vv_product_integral(int i, int j)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    VV_PRODUCT_INTEGRAL: integral (-1<=x<=1) V(i,x)*V(j,x)*sqrt(1+x)/sqrt(1-x) dx
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    25 April 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int I, J, the polynomial indices.
            //    0 <= I, J.
            //
            //    Output, double VV_PRODUCT_INTEGRAL, the value of the integral.
            //
        {
            
            double value;

            if (i < 0)
            {
                Console.WriteLine("");
                Console.WriteLine("VV_PRODUCT_INTEGRAL - Fatal error!");
                Console.WriteLine("  0 <= I, is required.");
                return (1);
            }

            if (j < 0)
            {
                Console.WriteLine("");
                Console.WriteLine("VV_PRODUCT_INTEGRAL - Fatal error!");
                Console.WriteLine("  0 <= J is required.");
                return (1);
            }

            if (i != j)
            {
                value = 0.0;
            }
            else
            {
                value = Math.PI;
            }

            return value;
        }
    }
}