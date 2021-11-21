using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using Burkardt.MatrixNS;
using Burkardt.Types;

namespace Burkardt.ChebyshevPolynomialNS;

public static partial class ChebyshevPolynomial
{
    public static double[] cheby_u_zero ( int n )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHEBY_U_ZERO returns zeroes of the Chebyshev polynomial U(N)(X).
        //
        //  Discussion:
        //
        //    The I-th zero of U(N)(X) is cos((I-1)*PI/(N-1)), I = 1 to N
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 May 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the polynomial.
        //
        //    Output, double CHEBY_U_ZERO[N], the zeroes of U(N)(X).
        //
    {
        int i;

        double[] z = new double[n];

        for ( i = 0; i < n; i++ )
        {
            double angle = (i + 1) * Math.PI / (n + 1);
            z[i] = Math.Cos ( angle );
        }
        return z;
    }
    public static double[] u_mass_matrix(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    U_MASS_MATRIX computes the mass matrix for the Chebyshev U polynomial.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int k;

        double[] x = new double[n + 1];
        double[] w = new double[n + 1];

        u_quadrature_rule(n + 1, ref x, ref w);

        double[] phi = u_polynomial(n + 1, n, x);

        double[] phiw = new double[(n + 1) * (n + 1)];

        for (k = 0; k <= n; k++)
        {
            int i;
            for (i = 0; i <= n; i++)
            {
                phiw[i + k * (n + 1)] = w[k] * phi[k + i * (n + 1)];
            }
        }

        double[] a = typeMethods.r8mat_mm_new(n + 1, n + 1, n + 1, phiw, phi);

        return a;
    }

    public static double u_moment(int e)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    U_MOMENT: integral ( -1 <= x <= +1 ) x^e sqrt ( 1 - x^2 ) dx.
        //
        //  Discussion:
        //
        //     E    U_MOMENT
        //    --    -------------- 
        //     0         Math.PI /    2 
        //     2         Math.PI /    8
        //     4         Math.PI /   16
        //     6     5 * Math.PI /  128
        //     8     7 * Math.PI /  256
        //    10    21 * Math.PI / 1024
        //    12    33 * Math.PI / 2048
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 September 2013
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
        //    Output, double U_MOMENT, the value of the integral.
        //
    {
        double value;

        switch (e % 2)
        {
            case 1:
                value = 0.0;
                break;
            default:
                double arg1 = 0.5 * (1 + e);
                double arg2 = 2.0 + 0.5 * e;
                value = 0.5 * Math.Sqrt(Math.PI) * Helpers.Gamma(arg1) / Helpers.Gamma(arg2);
                break;
        }

        return value;
    }

    public static double[] u_polynomial(int m, int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    U_POLYNOMIAL evaluates Chebyshev polynomials U(n,x).
        //
        //  Differential equation:
        //
        //    (1-X*X) Y'' - 3 X Y' + N (N+2) Y = 0
        //
        //  First terms:
        //
        //    U(0,X) =   1
        //    U(1,X) =   2 X
        //    U(2,X) =   4 X^2 -   1
        //    U(3,X) =   8 X^3 -   4 X
        //    U(4,X) =  16 X^4 -  12 X^2 +  1
        //    U(5,X) =  32 X^5 -  32 X^3 +  6 X
        //    U(6,X) =  64 X^6 -  80 X^4 + 24 X^2 - 1
        //    U(7,X) = 128 X^7 - 192 X^5 + 80 X^3 - 8X
        //
        //  Recursion:
        //
        //    U(0,X) = 1,
        //    U(1,X) = 2 * X,
        //    U(N,X) = 2 * X * U(N-1,X) - U(N-2,X)
        //
        //  Norm:
        //
        //    Integral ( -1 <= X <= 1 ) ( 1 - X^2 ) * U(N,X)^2 dX = PI/2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 April 2012
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
        //    Output, double U_POLYNOMIAL[M*(N+1)], the values of the N+1 
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
            v[i + 1 * m] = 2.0 * x[i];
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

    public static void u_polynomial_01_values(ref int n_data, ref int n, ref double x, ref double fx )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    U_POLYNOMIAL_01_VALUES: values of shifted Chebyshev polynomials U01(n,x).
        //
        //  Discussion:
        //
        //    U01(n,x) = U(n,2*x-1)
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
                1.400000000000000,
                0.9600000000000000,
                -0.05600000000000000,
                -1.038400000000000,
                -1.397760000000000,
                -0.9184640000000000,
                0.1119104000000000,
                1.075138560000000,
                1.393283584000000,
                0.8754584576000000,
                -0.1676417433600000,
                -1.110156898304000,
                -8.000000000000000,
                1.511014400000000,
                -1.133260800000000,
                -0.1636352000000000,
                1.019801600000000,
                0.000000000000000,
                -1.019801600000000,
                0.1636352000000000,
                1.133260800000000,
                -1.511014400000000,
                8.000000000000000
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

    public static double[] u_polynomial_ab(double a, double b, int m, int n, double[] xab)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    U_POLYNOMIAL_AB: Chebyshev polynomials UAB(n,x) in [A,B].
        //
        //  Discussion:
        //
        //    UAB(n,x) = U(n,(2*x-a-b)/(b-a))
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
        //    Output, double U_POLYNOMIAL_AB[M*(N+1)], the values.
        //
    {
        int i;

        double[] x = new double[m];

        for (i = 0; i < m; i++)
        {
            x[i] = (2.0 * xab[i] - a - b) / (b - a);
        }

        double[] v = u_polynomial(m, n, x);
        return v;
    }

    public static double u_polynomial_ab_value(double a, double b, int n, double xab)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    U_POLYNOMIAL_AB_VALUE: Chebyshev polynomial UAB(n,x) in [A,B].
        //
        //  Discussion:
        //
        //    UAB(n,x) = U(n,(2*x-a-b)/(b-a))
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
        //    Output, double U_POLYNOMIAL_AB_VALUE, the value.
        //
    {
        double x = (2.0 * xab - a - b) / (b - a);

        double v = u_polynomial_value(n, x);

        return v;
    }

    public static double[] u_polynomial_coefficients(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    U_POLYNOMIAL_COEFFICIENTS evaluates coefficients of Chebyshev polynomials U(n,x).
        //
        //  First terms:
        //
        //    N/K     0     1      2      3       4     5      6    7      8    9   10
        //
        //     0      1
        //     1      0     2
        //     2     -1     0      4
        //     3      0    -4      0      8
        //     4      1     0    -12      0      16
        //     5      0     6      0    -32       0    32
        //     6     -1     0     24      0     -80     0     64
        //     7      0    -8      0     80       0  -192      0   128
        //
        //  Recursion:
        //
        //    U(0,X) = 1,
        //    U(1,X) = 2*X,
        //    U(N,X) = 2 * X * U(N-1,X) - U(N-2,X)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 February 2003
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
        //    Output, double U_POLYNOMIAL_COEFFICIENTS[(N+1)*((N+1)], the coefficients 
        //    of the Chebyshev U polynomials.
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

    public static void u_polynomial_plot(int n_num, int[] n_val, string output_filename )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    U_POLYNOMIAL_PLOT plots Chebyshev polynomials U(n,x).
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

        double a = -1.0;
        const double b = +1.0;

        double[] x = typeMethods.r8vec_linspace_new(m, a, b);
        //
        //  Compute all the data.
        //
        int n_max = typeMethods.i4vec_max(n_num, n_val);
        double[] v = u_polynomial(m, n_max, x);
        //
        //  Create the data file.
        //
        const string data_filename = "u_polynomial_data.txt";
        for (i = 0; i < m; i++)
        {
            string line = x[i].ToString(CultureInfo.InvariantCulture);
            for (j = 0; j < n_num; j++)
            {
                int n = n_val[j];
                line  += "  " + v[i + n * m];
            }

            data_unit.Add(line);
        }

        File.WriteAllLines(data_filename, data_unit);
        Console.WriteLine("");
        Console.WriteLine("  Created graphics data file '" + data_filename + "'.");
        //
        //  Plot the selected data.
        //
        string command_filename = "u_polynomial_commands.txt";

        command_unit.Add("# " + command_filename + "");
        command_unit.Add("#");
        command_unit.Add("# Usage:");
        command_unit.Add("#  gnuplot < " + command_filename + "");
        command_unit.Add("#");
        command_unit.Add("set term png");
        command_unit.Add("set nokey");
        command_unit.Add("set output '" + output_filename + "'");
        command_unit.Add("set xlabel '<---X--->'");
        command_unit.Add("set ylabel '<---U(n,x)--->'");
        command_unit.Add("set title 'Chebyshev Polynomials U(n,x)'");
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

        
    public static double u_polynomial_value(int n, double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    U_POLYNOMIAL_VALUE: returns the single value U(n,x).
        //
        //  Discussion:
        //
        //    In cases where calling U_POLYNOMIAL is inconvenient, because it returns
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
        //    Output, double U_POLYNOMIAL_VALUE, the value of U(n,x).
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
                int m = 1;
                x_vec[0] = x;

                double[] v_vec = u_polynomial(m, n, x_vec);

                value = v_vec[n];
                break;
        }

        return value;
    }

    public static void u_polynomial_values(ref int n_data, ref int n, ref double x, ref double fx )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    U_POLYNOMIAL_VALUES returns values of Chebyshev polynomials U(n,x).
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      ChebyshevU[n,x]
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
                0.0000000000000000E+00,
                0.1000000000000000E+01,
                0.1600000000000000E+01,
                0.1560000000000000E+01,
                0.8960000000000000E+00,
                -0.1264000000000000E+00,
                -0.1098240000000000E+01,
                -0.1630784000000000E+01,
                -0.1511014400000000E+01,
                -0.7868390400000000E+00,
                0.2520719360000000E+00,
                0.1190154137600000E+01,
                0.1652174684160000E+01,
                0.1453325357056000E+01
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

    public static double[] u_polynomial_zeros(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    U_POLYNOMIAL_ZEROS returns zeroes of Chebyshev polynomials U(n,x).
        //
        //  Discussion:
        //
        //    The I-th zero is cos((I-1)*PI/(N-1)), I = 1 to N
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 March 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the polynomial.
        //
        //    Output, double U_POLYNOMIAL_ZEROS[N], the zeroes.
        //
    {
        int i;

        double[] z = new double[n];

        for (i = 1; i <= n; i++)
        {
            double angle = i * Math.PI / (n + 1);
            z[i - 1] = Math.Cos(angle);
        }

        return z;
    }

    public static void u_quadrature_rule(int n, ref double[] t, ref double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    U_QUADRATURE_RULE: quadrature rule for U(n,x).
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
        int i;
            

        for (i = 0; i < n; i++)
        {
            t[i] = 0.0;
        }

        double[] bj = new double[n];
        for (i = 0; i < n; i++)
        {
            bj[i] = 0.5;
        }

        w[0] = Math.Sqrt(Math.PI / 2.0);
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

    public static double uu_product(int i, int j, double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UU_PRODUCT: evaluate U(i,x)*U(j,x)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I, J, the indices.
        //
        //    Input, double X, the argument.
        //
        //    Output, double UU_PRODUCT, the value.
        //
    {
        int k;
        double value = 0;

        value = 0.0;
        for (k = Math.Abs(i - j); k <= i + j; k += 2)
        {
            value += u_polynomial_value(k, x);
        }

        return value;
    }

    public static double uu_product_integral(int i, int j)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UU_PRODUCT_INTEGRAL: integral (-1<=x<=1) U(i,x)*U(j,x)*sqrt(1-x^2) dx
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
        //    Output, double UU_PRODUCT_INTEGRAL, the value of the integral.
        //
    {
            
        double value = 0;

        switch (i)
        {
            case < 0:
                Console.WriteLine("");
                Console.WriteLine("UU_PRODUCT_INTEGRAL - Fatal error!");
                Console.WriteLine("  0 <= I, is required.");
                return 1;
        }

        switch (j)
        {
            case < 0:
                Console.WriteLine("");
                Console.WriteLine("UU_PRODUCT_INTEGRAL - Fatal error!");
                Console.WriteLine("  0 <= J is required.");
                return 1;
        }

        if (i != j)
        {
            value = 0.0;
        }
        else
        {
            value = Math.PI / 2.0;
        }

        return value;
    }
}