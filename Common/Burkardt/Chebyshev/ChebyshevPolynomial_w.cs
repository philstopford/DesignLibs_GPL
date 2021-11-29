using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using Burkardt.MatrixNS;
using Burkardt.Types;

namespace Burkardt.ChebyshevPolynomialNS;

public static partial class ChebyshevPolynomial
{
    public static double[] w_mass_matrix(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    W_MASS_MATRIX computes the mass matrix for the Chebyshev W polynomial.
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
        int i;

        double[] x = new double[n + 1];
        double[] w = new double[n + 1];

        w_quadrature_rule(n + 1, ref x, ref w);

        double[] phi = w_polynomial(n + 1, n, x);

        double[] phiw = new double[(n + 1) * (n + 1)];

        for (i = 0; i <= n; i++)
        {
            int j;
            for (j = 0; j <= n; j++)
            {
                phiw[j + i * (n + 1)] = w[i] * phi[i + j * (n + 1)];
            }
        }

        double[] a = typeMethods.r8mat_mm_new(n + 1, n + 1, n + 1, phiw, phi);

        return a;
    }

    public static double w_moment(int e)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    W_MOMENT: integral ( -1 <= x <= +1 ) x^e sqrt(1-x) / sqrt(1+x) dx.
        //
        //  Discussion:
        //
        //     E    W_MOMENT
        //    --    --------------
        //     0        pi
        //     1  -     Math.PI / 2
        //     2        Math.PI / 2
        //     3  -   3 Math.PI / 8
        //     4      3 Math.PI / 8
        //     5  -   5 Math.PI / 16
        //     6      5 Math.PI / 16
        //     7  -  35 Math.PI / 128
        //     8     35 Math.PI / 128
        //     9  -  63 Math.PI / 256
        //    10     63 Math.PI / 256
        //    11  - 231 Math.PI / 1024
        //    12    231 Math.PI / 1024
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
        double r8_e = e;

        double f1 = 1.0 / Helpers.Gamma(1.5 + r8_e);
        double f2 = typeMethods.r8_mop(e);
        double f3 = Math.PI * Helpers.Gamma(1.5 + r8_e);
        double f4 = 2.0 * typeMethods.r8_mop(e) * typeMethods.r8_hyper_2f1(0.5, -r8_e, 1.0, 2.0);
        double f5 = (-1.0 + typeMethods.r8_mop(e)) * typeMethods.r8_hyper_2f1(0.5, -r8_e, 2.0, 2.0);
        double f6 = Math.Sqrt(Math.PI) * typeMethods.r8_factorial(e);
        double f7 = (-1.0 + typeMethods.r8_mop(e))
                    * typeMethods.r8_hyper_2f1(-0.5, 1.0 + r8_e, 1.5 + r8_e, -1.0);
        double f8 = 2.0 * typeMethods.r8_mop(e) * typeMethods.r8_hyper_2f1(0.5, 1.0 + r8_e, 1.5 + r8_e, -1.0);

        double value = f1 * f2 * (f3 * (f4 - f5) + f6 * (f7 - f8));

        return value;
    }

    public static double[] w_polynomial(int m, int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    W_POLYNOMIAL evaluates Chebyshev polynomials W(n,x).
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
        //    Output, double W_POLYNOMIAL[M*(N+1)], the values of the N+1 
        //    Chebyshev polynomials.
        //
    {
        int i;

        switch (n)
        {
            case < 0:
                return null;
        }

        double[] v = new double[m * (n + 1)];

        for (i = 0; i < m; i++)
        {
            v[i + 0 * m] = 1.0;
        }

        switch (n)
        {
            case < 1:
                return v;
        }

        for (i = 0; i < m; i++)
        {
            v[i + 1 * m] = 2.0 * x[i] + 1.0;
        }

        for (i = 0; i < m; i++)
        {
            int j;
            for (j = 2; j <= n; j++)
            {
                v[i + j * m] = 2.0 * x[i] * v[i + (j - 1) * m] - v[i + (j - 2) * m];
            }
        }

        return v;
    }

    public static void w_polynomial_01_values(ref int n_data, ref int n, ref double x, ref double fx )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    W_POLYNOMIAL_01_VALUES: values of shifted Chebyshev polynomials W01(n,x).
        //
        //  Discussion:
        //
        //    W01(n,x) = W(n,2*x-1)
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
        const int N_MAX = 25;

        double[] fx_vec =
            {
                0.000000000000000,
                1.000000000000000,
                2.400000000000000,
                2.360000000000000,
                0.904000000000000,
                -1.094400000000000,
                -2.436160000000000,
                -2.316224000000000,
                -0.806553600000000,
                1.187048960000000,
                2.468422144000000,
                2.268742041600000,
                0.707816714240000,
                -1.277798641664000,
                -1.000000000000000,
                -0.119769600000000,
                -0.875276800000000,
                0.890508800000000,
                0.855897600000000,
                -1.000000000000000,
                -1.183705600000000,
                1.217779200000000,
                1.391244800000000,
                -3.141798400000000,
                15.00000000000000
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

        n_data = n_data switch
        {
            < 0 => 0,
            _ => n_data
        };

        n_data += 1;

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

    public static double[] w_polynomial_ab(double a, double b, int m, int n, double[] xab)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    W_POLYNOMIAL_AB: Chebyshev polynomials WAB(n,x) in [A,B].
        //
        //  Discussion:
        //
        //    WAB(n,x) = W(n,(2*x-a-b)/(b-a))
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
        //    Output, double W_POLYNOMIAL_AB[M*(N+1)], the values.
        //
    {
        int i;

        double[] x = new double[m];

        for (i = 0; i < m; i++)
        {
            x[i] = (2.0 * xab[i] - a - b) / (b - a);
        }

        double[] v = w_polynomial(m, n, x);

        return v;
    }

    public static double w_polynomial_ab_value(double a, double b, int n, double xab)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    W_POLYNOMIAL_AB_VALUE: Chebyshev polynomial WAB(n,x) in [A,B].
        //
        //  Discussion:
        //
        //    WAB(n,x) = W(n,(2*x-a-b)/(b-a))
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
        //    Output, double W_POLYNOMIAL_AB_VALUE, the value.
        //
    {
        double x = (2.0 * xab - a - b) / (b - a);

        double v = w_polynomial_value(n, x);

        return v;
    }

    public static double[] w_polynomial_coefficients(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    W_POLYNOMIAL_COEFFICIENTS: coefficients of the Chebyshev polynomial W(n,x).
        //
        //  First terms:
        //
        //    N/K     0     1      2      3       4     5      6    7      8    9   10
        //
        //     0      1
        //     1      1     2
        //     2     -1     2      4
        //     3     -1    -4      4      8
        //     4      1    -4    -12      8      16
        //     5      1     6    -12    -32     +16    32
        //     6     -1     6     24    -32     -80    32     64
        //     7     -1    -8    +24    +80     -80  -192     64   128
        //
        //  Recursion:
        //
        //    W(0,X) = 1,
        //    W(1,X) = 2 * X + 1,
        //    W(N,X) = 2 * X * T(N-1,X) - T(N-2,X)
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
        //    Output, double W_POLYNOMIAL_COEFFICIENTS[(N+1)*(N+1)], the coefficients.
        //
    {
        int i;
        int j;

        switch (n)
        {
            case < 0:
                return null;
        }

        double[] c = new double[(n + 1) * (n + 1)];

        for (i = 0; i <= n; i++)
        {
            for (j = 0; j <= n; j++)
            {
                c[i + j * (n + 1)] = 0.0;
            }
        }

        c[0 + 0 * (n + 1)] = 1.0;

        switch (n)
        {
            case 0:
                return c;
        }

        c[1] = 1.0;
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

    public static void w_polynomial_plot(int n_num, int[] n_val, string output_filename )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    W_POLYNOMIAL_PLOT plots Chebyshev polynomials W(n,x).
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
        List<string> command_unit = new();
        List<string> data_unit = new();
        int i;
        int j;
        const int m = 501;

        const double a = -1.0;
        const double b = +1.0;

        double[] x = typeMethods.r8vec_linspace_new(m, a, b);
        //
        //  Compute all the data.
        //
        int n_max = typeMethods.i4vec_max(n_num, n_val);
        double[] v = w_polynomial(m, n_max, x);
        //
        //  Create the data file.
        //
        const string data_filename = "w_polynomial_data.txt";
        for (i = 0; i < m; i++)
        {
            string line = x[i].ToString(CultureInfo.InvariantCulture);
            for (j = 0; j < n_num; j++)
            {
                int n = n_val[j];
                line += "  " + v[i + n * m].ToString(CultureInfo.InvariantCulture);
            }

            data_unit.Add(line);
        }

        File.WriteAllLines(data_filename, data_unit);
            
        Console.WriteLine("");
        Console.WriteLine("  Created graphics data file '" + data_filename + "'.");
        //
        //  Plot the selected data.
        //
        const string command_filename = "w_polynomial_commands.txt";

        command_unit.Add("# " + command_filename + "");
        command_unit.Add("#");
        command_unit.Add("# Usage:");
        command_unit.Add("#  gnuplot < " + command_filename + "");
        command_unit.Add("#");
        command_unit.Add("set term png");
        command_unit.Add("set nokey");
        command_unit.Add("set output '" + output_filename + "'");
        command_unit.Add("set xlabel '<---X--->'");
        command_unit.Add("set ylabel '<---W(n,x)--->'");
        command_unit.Add("set title 'Chebyshev Polynomials W(n,x)'");
        command_unit.Add("set grid");
        command_unit.Add("set style data lines");
        for (j = 0; j < n_num; j++)
        {
            string line = "";
            int column = n_val[j] + 1;
            line += j switch
            {
                0 => "plot ",
                _ => "     "
            };

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

    public static double w_polynomial_value(int n, double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    W_POLYNOMIAL_VALUE: returns the single value W(n,x).
        //
        //  Discussion:
        //
        //    In cases where calling W_POLYNOMIAL is inconvenient, because it returns
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
        //    Output, double W_POLYNOMIAL_VALUE, the value of W(n,x).
        //
    {
        double value;
        double[] x_vec = new double[1];

        switch (n)
        {
            case < 0:
                value = 0.0;
                break;
            default:
                const int m = 1;
                x_vec[0] = x;

                double[] v_vec = w_polynomial(m, n, x_vec);

                value = v_vec[n];
                break;
        }

        return value;
    }

    public static void w_polynomial_values(ref int n_data, ref int n, ref double x, ref double fx )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    W_POLYNOMIAL_VALUES returns values of Chebyshev polynomials W(n,x).
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      u = Sqrt[(x+1)/2],
        //      ChebyshevU[2*n,u]
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
        const int N_MAX = 14;

        double[] fx_vec =
            {
                0.000000000000000E+00,
                1.000000000000000E+00,
                2.600000000000000E+00,
                3.160000000000000E+00,
                2.456000000000000E+00,
                0.769600000000000E+00,
                -1.224640000000000E+00,
                -2.729024000000000E+00,
                -3.141798400000000E+00,
                -2.297853440000000E+00,
                -0.534767104000000E+00,
                1.442226073600000E+00,
                2.842328821760000E+00,
                3.105500041216000E+00
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

        n_data = n_data switch
        {
            < 0 => 0,
            _ => n_data
        };

        n_data += 1;

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

    public static double[] w_polynomial_zeros(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    W_POLYNOMIAL_ZEROS returns zeroes of the Chebyshev polynomial W(n,x).
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
        //    Output, double W_POLYNOMIAL_ZEROS[N], the zeroes.
        //
    {
        double[] z = new double[n];

        for (int i = 0; i < n; i++)
        {
            double angle = 2 * (n - i) * Math.PI / (2 * n + 1);
            z[i] = Math.Cos(angle);
        }

        return z;
    }

    public static void w_quadrature_rule(int n, ref double[] t, ref double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    W_QUADRATURE_RULE: quadrature rule for W(n,x).
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
        //    Input, int N, the order of the rule.
        //
        //    Output, double T[N], W[N], the points and weights of the rule.
        //
    {
        int i;
            

        for (i = 0; i < n; i++)
        {
            t[i] = 0.0;
        }

        t[0] = -0.5;

        double[] bj = new double[n];
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
            w[i] *= w[i];
        }
    }

    public static double ww_product_integral(int i, int j)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    WW_PRODUCT_INTEGRAL: integral (-1<=x<=1) W(i,x)*W(j,x)*sqrt(1-x)/sqrt(1+x) dx
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
        //    Output, double WW_PRODUCT_INTEGRAL, the value of the integral.
        //
    {
            
        double value = 0;

        switch (i)
        {
            case < 0:
                Console.WriteLine("");
                Console.WriteLine("WW_PRODUCT_INTEGRAL - Fatal error!");
                Console.WriteLine("  0 <= I, is required.");
                return 1;
        }

        switch (j)
        {
            case < 0:
                Console.WriteLine("");
                Console.WriteLine("WW_PRODUCT_INTEGRAL - Fatal error!");
                Console.WriteLine("  0 <= J is required.");
                return 1;
        }

        value = i != j ? 0.0 : Math.PI;

        return value;
    }
}