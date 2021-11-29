using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using Burkardt.MatrixNS;
using Burkardt.SolveNS;
using Burkardt.Types;

namespace Burkardt.ChebyshevPolynomialNS;

public static partial class ChebyshevPolynomial
{
    public static double[] cheby_t_zero ( int n )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHEBY_T_ZERO returns zeroes of the Chebyshev polynomial T(N)(X).
        //
        //  Discussion:
        //
        //    The I-th zero of T(N)(X) is cos((2*I-1)*PI/(2*N)), I = 1 to N
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
        //    Output, double CHEBY_T_ZERO[N], the zeroes of T(N)(X).
        //
    {
        int i;

        double[] z = new double[n];

        for ( i = 0; i < n; i++ )
        {
            double angle = (2 * i + 1) * Math.PI / (2 * n);
            z[i] = Math.Cos ( angle );
        }
        return z;
    }
    public static double[] t_mass_matrix(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    T_MASS_MATRIX computes the mass matrix for the Chebyshev T polynomial.
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
        int i;

        double[] x = new double[n + 1];
        double[] w = new double[n + 1];

        t_quadrature_rule(n + 1, ref x, ref w);

        double[] phi = t_polynomial(n + 1, n, x);

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

    public static double t_moment(int e)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    T_MOMENT: integral ( -1 <= x <= +1 ) x^e dx / sqrt ( 1 - x^2 ).
        //
        //  Discussion:
        //
        //    Set 
        //      x = cos ( theta ), 
        //      dx = - sin ( theta ) d theta = - sqrt ( 1 - x^2 ) d theta
        //    to transform the integral to
        //      integral ( 0 <= theta <= Math.PI ) - ( cos ( theta ) )^e d theta
        //    which becomes
        //      0 if E is odd,
        //      (1/2^e) * choose ( e, e/2 ) * Math.PI if E is even.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 April 2012
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
        //    Output, double T_MOMENT, the value of the integral.
        //
    {
        double value = (e % 2) switch
        {
            1 => 0.0,
            _ => typeMethods.r8_choose(e, e / 2) * Math.PI / Math.Pow(2.0, e)
        };

        return value;
    }

    public static double[] t_polynomial(int m, int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    T_POLYNOMIAL evaluates Chebyshev polynomials T(n,x).
        //
        //  Discussion:
        //
        //    Chebyshev polynomials are useful as a basis for representing the
        //    approximation of functions since they are well conditioned, in the sense
        //    that in the interval [-1,1] they each have maximum absolute value 1.
        //    Hence an error in the value of a coefficient of the approximation, of
        //    size epsilon, is exactly reflected in an error of size epsilon between
        //    the computed approximation and the theoretical approximation.
        //
        //    Typical usage is as follows, where we assume for the moment
        //    that the interval of approximation is [-1,1].  The value
        //    of N is chosen, the highest polynomial to be used in the
        //    approximation.  Then the function to be approximated is
        //    evaluated at the N+1 points XJ which are the zeroes of the N+1-th
        //    Chebyshev polynomial.  Let these values be denoted by F(XJ).
        //
        //    The coefficients of the approximation are now defined by
        //
        //      C(I) = 2/(N+1) * sum ( 1 <= J <= N+1 ) F(XJ) T(I,XJ)
        //
        //    except that C(0) is given a value which is half that assigned
        //    to it by the above formula,
        //
        //    and the representation is
        //
        //    F(X) approximated by sum ( 0 <= J <= N ) C(J) T(J,X)
        //
        //    Now note that, again because of the fact that the Chebyshev polynomials
        //    have maximum absolute value 1, if the higher order terms of the
        //    coefficients C are small, then we have the option of truncating
        //    the approximation by dropping these terms, and we will have an
        //    exact value for maximum perturbation to the approximation that
        //    this will cause.
        //
        //    It should be noted that typically the error in approximation
        //    is dominated by the first neglected basis function (some multiple of
        //    T(N+1,X) in the example above).  If this term were the exact error,
        //    then we would have found the minimax polynomial, the approximating
        //    polynomial of smallest maximum deviation from the original function.
        //    The minimax polynomial is hard to compute, and another important
        //    feature of the Chebyshev approximation is that it tends to behave
        //    like the minimax polynomial while being easy to compute.
        //
        //    To evaluate a sum like
        //
        //      sum ( 0 <= J <= N ) C(J) T(J,X),
        //
        //    Clenshaw's recurrence formula is recommended instead of computing the
        //    polynomial values, forming the products and summing.
        //
        //    Assuming that the coefficients C(J) have been computed
        //    for J = 0 to N, then the coefficients of the representation of the
        //    indefinite integral of the function may be computed by
        //
        //      B(I) = ( C(I-1) - C(I+1))/2*(I-1) for I=1 to N+1,
        //
        //    with
        //
        //      C(N+1)=0
        //      B(0) arbitrary.
        //
        //    Also, the coefficients of the representation of the derivative of the
        //    function may be computed by:
        //
        //      D(I) = D(I+2)+2*I*C(I) for I=N-1, N-2, ..., 0,
        //
        //    with
        //
        //      D(N+1) = D(N)=0.
        //
        //    Some of the above may have to adjusted because of the irregularity of C(0).
        //
        //    The formula is:
        //
        //      T(N,X) = COS(N*ARCCOS(X))
        //
        //  Differential equation:
        //
        //    (1-X*X) Y'' - X Y' + N N Y = 0
        //
        //  First terms:
        //
        //    T(0,X) =  1
        //    T(1,X) =  1 X
        //    T(2,X) =  2 X^2 -   1
        //    T(3,X) =  4 X^3 -   3 X
        //    T(4,X) =  8 X^4 -   8 X^2 +  1
        //    T(5,X) = 16 X^5 -  20 X^3 +  5 X
        //    T(6,X) = 32 X^6 -  48 X^4 + 18 X^2 - 1
        //    T(7,X) = 64 X^7 - 112 X^5 + 56 X^3 - 7 X
        //
        //  Inequality:
        //
        //    abs ( T(N,X) ) <= 1 for -1 <= X <= 1
        //
        //  Orthogonality:
        //
        //    For integration over [-1,1] with weight
        //
        //      W(X) = 1 / sqrt(1-X*X),
        //
        //    if we write the inner product of T(I,X) and T(J,X) as
        //
        //      < T(I,X), T(J,X) > = integral ( -1 <= X <= 1 ) W(X) T(I,X) T(J,X) dX
        //
        //    then the result is:
        //
        //      0 if I /= J
        //      PI/2 if I == J /= 0
        //      PI if I == J == 0
        //
        //    A discrete orthogonality relation is also satisfied at each of
        //    the N zeroes of T(N,X):  sum ( 1 <= K <= N ) T(I,X) * T(J,X)
        //                              = 0 if I /= J
        //                              = N/2 if I == J /= 0
        //                              = N if I == J == 0
        //
        //  Recursion:
        //
        //    T(0,X) = 1,
        //    T(1,X) = X,
        //    T(N,X) = 2 * X * T(N-1,X) - T(N-2,X)
        //
        //    T'(N,X) = N * ( -X * T(N,X) + T(N-1,X) ) / ( 1 - X^2 )
        //
        //  Special values:
        //
        //    T(N,1) = 1
        //    T(N,-1) = (-1)^N
        //    T(2N,0) = (-1)^N
        //    T(2N+1,0) = 0
        //    T(N,X) = (-1)**N * T(N,-X)
        //
        //  Zeroes:
        //
        //    M-th zero of T(N,X) is cos((2*M-1)*PI/(2*N)), M = 1 to N
        //
        //  Extrema:
        //
        //    M-th extremum of T(N,X) is cos(PI*M/N), M = 0 to N
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 May 2003
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
        //    Output, double T_POLYNOMIAL[M*(N+1)], the values of the Chebyshev polynomials.
        //
    {
        int i;
        int j;

        switch (n)
        {
            case < 0:
                return null;
        }

        double[] v = new double[m * (n + 1)];

        for (i = 0; i < m; i++)
        {
            v[i] = 1.0;
        }

        switch (n)
        {
            case < 1:
                return v;
        }

        for (i = 0; i < m; i++)
        {
            v[i + 1 * m] = x[i];
        }

        for (j = 2; j <= n; j++)
        {
            for (i = 0; i < m; i++)
            {
                v[i + j * m] = 2.0 * x[i] * v[i + (j - 1) * m] - v[i + (j - 2) * m];
            }
        }

        return v;
    }

    public static void t_polynomial_01_values(ref int n_data, ref int n, ref double x, ref double fx )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    T_POLYNOMIAL_01_VALUES: values of shifted Chebyshev polynomials T01(n,x).
        //
        //  Discussion:
        //
        //    T01(n,x) = T(n,2*x-1)
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
                0.0000000000000000,
                1.0000000000000000,
                0.7000000000000000,
                -0.0200000000000000,
                -0.7280000000000000,
                -0.9992000000000000,
                -0.6708800000000000,
                0.0599680000000000,
                0.7548352000000000,
                0.9968012800000000,
                0.6406865920000000,
                -0.0998400512000000,
                -0.7804626636800000,
                -0.9928076779520000,
                -1.0000000000000000,
                0.2063872000000000,
                -0.9784704000000000,
                0.2580224000000000,
                0.9870208000000000,
                0.0000000000000000,
                -0.9870208000000000,
                -0.2580224000000000,
                0.9784704000000000,
                -0.2063872000000000,
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

    public static double[] t_polynomial_ab(double a, double b, int m, int n, ref double[] xab)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    T_POLYNOMIAL_AB: Chebyshev polynomials TAB(n,x) in [A,B].
        //
        //  Discussion:
        //
        //    TAB(n,x) = T(n,(2*x-a-b)/(b-a))
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 April 2012
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
        //    Output, double T_POLYNOMIAL_AB[M*(N+1)], the values.
        //
    {
        int i;

        double[] x = new double[m];

        for (i = 0; i < m; i++)
        {
            x[i] = (2.0 * xab[i] - a - b) / (b - a);
        }

        double[] v = t_polynomial(m, n, x);

        return v;
    }

    public static double t_polynomial_ab_value(double a, double b, int n, double xab)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    T_POLYNOMIAL_AB_VALUE: Chebyshev polynomial TAB(n,x) in [A,B].
        //
        //  Discussion:
        //
        //    TAB(n,x) = T(n,(2*x-a-b)/(b-a))
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
        //    Output, double T_POLYNOMIAL_AB_VALUE, the value.
        //
    {
        double x = (2.0 * xab - a - b) / (b - a);

        double v = t_polynomial_value(n, x);

        return v;
    }

    public static double[] t_polynomial_coefficients(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    T_POLYNOMIAL_COEFFICIENTS: coefficients of the Chebyshev polynomial T(n,x).
        //
        //  First terms:
        //
        //    N/K     0     1      2      3       4     5      6    7      8    9   10
        //
        //     0      1
        //     1      0     1
        //     2     -1     0      2
        //     3      0    -3      0      4
        //     4      1     0     -8      0       8
        //     5      0     5      0    -20       0    16
        //     6     -1     0     18      0     -48     0     32
        //     7      0    -7      0     56       0  -112      0    64
        //
        //  Recursion:
        //
        //    T(0,X) = 1,
        //    T(1,X) = X,
        //    T(N,X) = 2 * X * T(N-1,X) - T(N-2,X)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    22 April 2012
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
        //    Output, double T_POLYNOMIAL_COEFFICIENTS[(N+1)*(N+1)], the coefficients.
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

        c[1 + 1 * (n + 1)] = 1.0;

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

    public static void t_polynomial_plot(int n_num, int[] n_val, string output_filename )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    T_POLYNOMIAL_PLOT plots Chebyshev polynomials T(n,x).
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
        double[] v = t_polynomial(m, n_max, x);
        //
        //  Create the data file.
        //
        const string data_filename = "t_polynomial_data.txt";

        for (i = 0; i < m; i++)
        {
            string line = x[i].ToString(CultureInfo.InvariantCulture);
            for (j = 0; j < n_num; j++)
            {
                int n = n_val[j];
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
        const string command_filename = "t_polynomial_commands.txt";

        command_unit.Add("# " + command_filename + "");
        command_unit.Add("#");
        command_unit.Add("# Usage:");
        command_unit.Add("#  gnuplot < " + command_filename + "");
        command_unit.Add("#");
        command_unit.Add("set term png");
        command_unit.Add("set nokey");
        command_unit.Add("set output '" + output_filename + "'");
        command_unit.Add("set xlabel '<---X--->'");
        command_unit.Add("set ylabel '<---T(n,x)--->'");
        command_unit.Add("set title 'Chebyshev Polynomials T(n,x)'");
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

    public static double t_polynomial_value(int n, double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    T_POLYNOMIAL_VALUE: returns the single value T(n,x).
        //
        //  Discussion:
        //
        //    In cases where calling T_POLYNOMIAL is inconvenient, because it returns
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
        //    Output, double T_POLYNOMIAL_VALUE, the value of T(n,x).
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

                double[] v_vec = t_polynomial(m, n, x_vec);

                value = v_vec[n];
                break;
        }

        return value;
    }

    public static void t_polynomial_values(ref int n_data, ref int n, ref double x, ref double fx )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    T_POLYNOMIAL_VALUES returns values of the Chebyshev polynomial T(n,x).
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      ChebyshevT[n,x]
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
                0.8000000000000000E+00,
                0.2800000000000000E+00,
                -0.3520000000000000E+00,
                -0.8432000000000000E+00,
                -0.9971200000000000E+00,
                -0.7521920000000000E+00,
                -0.2063872000000000E+00,
                0.4219724800000000E+00,
                0.8815431680000000E+00,
                0.9884965888000000E+00,
                0.7000513740800000E+00,
                0.1315856097280000E+00
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

    public static double[] t_polynomial_zeros(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    T_POLYNOMIAL_ZEROS returns zeroes of the Chebyshev polynomial T(n,x).
        //
        //  Discussion:
        //
        //    The I-th zero is cos((2*I-1)*PI/(2*N)), I = 1 to N
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
        //    Output, double T_POLYNOMIAL_ZEROS[N], the zeroes.
        //
    {
        int i;

        double[] z = new double[n];

        for (i = 1; i <= n; i++)
        {
            double angle = (2 * i - 1) * Math.PI / (2 * n);
            z[i - 1] = Math.Cos(angle);
        }

        return z;
    }

    public static double[] t_project_coefficients(int n, Func<double, double> f)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    T_PROJECT_COEFFICIENTS: function projected onto Chebyshev polynomials T(n,x).
        //
        //  Discussion:
        //
        //    It is assumed that the interval of definition is -1 <= x <= +1.
        //
        //    Over this interval, f(x) will be well approximated by
        //
        //      f(x) approx sum ( 0 <= i <= n ) c(i) * T(i,x)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the highest order polynomial to compute.
        //
        //    Input, double F ( double X ), evaluates the function.
        //
        //    Output, double T_PROJECT_COEFFICIENTS[N+1], the projection coefficients 
        //    of f(x) onto T(0,x) through T(n,x).
        //
    {
        int j;
        int k;

        double[] d = new double[n + 1];

        for (k = 0; k <= n; k++)
        {
            double y = Math.Cos(Math.PI * (k + 0.5) / (n + 1));
            d[k] = f(y);
        }

        double fac = 2.0 / (n + 1);

        double[] c = new double[n + 1];

        for (j = 0; j <= n; j++)
        {
            double total = 0.0;
            for (k = 0; k <= n; k++)
            {
                total += d[k] * Math.Cos(Math.PI * j
                                                 * ((k + 0.5) / (n + 1)));
            }

            c[j] = fac * total;
        }

        c[0] /= 2.0;
        return c;
    }

    public static double[] t_project_coefficients_ab(int n, Func<double, double> f, double a,
            double b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    T_PROJECT_COEFFICIENTS_AB: function projected onto TAB(n,x) over [a,b]
        //
        //  Discussion:
        //
        //    TAB(n,x) = T(n,(2*x-a-b)/(b-a))
        //
        //    It is assumed that the interval of definition is a <= x <= b.
        //
        //    Over this interval, f(x) will be well approximated by
        //
        //      f(x) approx sum ( 0 <= i <= n ) c(i) * T(i,(2x-a-b)/(b-a))
        //
        //    where x* = ( x - b - a ) / ( b - a )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the highest order polynomial to compute.
        //
        //    Input, double F ( double X ), evaluates the function.
        //
        //    Input, double A, B, the interval of definition.
        //
        //    Output, double T_PROJECT_COEFFICIENTS_AB[N+1], the projection coefficients 
        //    of f(x) onto T(0,x) through T(n,x).
        //
    {
        int j;
        int k;

        double[] d = new double[n + 1];

        for (k = 0; k <= n; k++)
        {
            double t = Math.Cos(Math.PI * (k + 0.5) / (n + 1));

            double y = ((1.0 + t) * b
                        + (1.0 - t) * a)
                       / 2.0;

            d[k] = f(y);
        }

        double fac = 2.0 / (n + 1);

        double[] c = new double[n + 1];

        for (j = 0; j <= n; j++)
        {
            double total = 0.0;
            for (k = 0; k <= n; k++)
            {
                total += d[k] * Math.Cos(Math.PI * j
                                                 * ((k + 0.5) / (n + 1)));
            }

            c[j] = fac * total;
        }

        c[0] /= 2.0;

        return c;
    }

    public static double[] t_project_coefficients_data(double a, double b, int m, int n,
            double[] x, double[] d )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    T_PROJECT_COEFFICIENTS_DATA: project data onto Chebyshev polynomials T(n,x).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, the domain of definition.
        //
        //    Input, int M, the number of data values.
        //
        //    Input, int N, the desired order of the Chebyshev 
        //    expansion.
        //
        //    Input, double X[M], the data abscissas.  These need not 
        //    be sorted.  It must be the case that A <= X() <= B.
        //
        //    Input, double D[M], the data values.
        //
        //    Output, double T_PROJECT_COEFFICIENTS_DATA[N+1], the approximate 
        //    Chebshev coefficients.
        //
    {
        if (!typeMethods.r8vec_in_ab(m, x, a, b))
        {
            Console.WriteLine("");
            Console.WriteLine("T_PROJECT_COEFFICIENTS_DATA- Fatal error!");
            Console.WriteLine("  Some X not in [A,B].");
            return null;
        }

        //
        //  Compute the M by N+1 Chebyshev Vandermonde matrix V.
        //
        double[] v = t_polynomial_ab(a, b, m, n, ref x);
        //
        //  Compute the least-squares solution C.
        //
        double[] c = SVDSolve.svd_solve(m, n + 1, v, d);

        return c;
    }

    public static double[] t_project_value(int m, int n, double[] x, double[] c )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    T_PROJECT_VALUE evaluates an expansion in Chebyshev polynomials T(n,x).
        //
        //  Discussion:
        //
        //    The projection is assumed to be based on the interval [-1,+1].
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
        //    Input, int M, the number of evaluation points.
        //
        //    Input, int N, the highest order polynomial to compute.
        //
        //    Input, double X[M], the evaluation points.
        //
        //    Input, double C[N+1], the expansion coefficients.
        //
        //    Output, double T_PROJECT_VALUE[M], the value of the Chebyshev function.
        //
    {
        int i;
        int j;

        double[] b0 = new double[m];
        double[] b1 = new double[m];
        double[] b2 = new double[m];

        for (i = 0; i < m; i++)
        {
            b0[i] = 0.0;
        }

        for (i = 0; i < m; i++)
        {
            b1[i] = 0.0;
        }

        for (j = n; 0 <= j; j--)
        {
            for (i = 0; i < m; i++)
            {
                b2[i] = b1[i];
                b1[i] = b0[i];
                b0[i] = c[j] + 2.0 * x[i] * b1[i] - b2[i];
            }
        }

        double[] v = new double[m];

        for (i = 0; i < m; i++)
        {
            v[i] = 0.5 * (c[0] + b0[i] - b2[i]);
        }

        return v;
    }

    public static double[] t_project_value_ab(int m, int n, double[] x, double[] c, double a,
            double b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    T_PROJECT_VALUE_AB evaluates an expansion in Chebyshev polynomials TAB(n,x).
        //
        //  Discussion:
        //
        //    TAB(n,x) = T(n,(2*x-a-b)/(b-a))
        //
        //    The projection is assumed to be based on the interval [A,B].
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
        //    Input, int M, the number of evaluation points.
        //
        //    Input, int N, the highest order polynomial to compute.
        //
        //    Input, double X[M], the evaluation points.
        //
        //    Input, double C[N+1], the expansion coefficients.
        //
        //    Input, double A, B, the interval of definition.
        //
        //    Output, double T_PROJECT_VALUE_AB[M], the value of the Chebyshev function.
        //
    {
        int i;
        int j;

        double[] b0 = new double[m];
        double[] b1 = new double[m];
        double[] b2 = new double[m];

        for (i = 0; i < m; i++)
        {
            b0[i] = 0.0;
        }

        for (i = 0; i < m; i++)
        {
            b1[i] = 0.0;
        }

        for (j = n; 0 <= j; j--)
        {
            for (i = 0; i < m; i++)
            {
                b2[i] = b1[i];
                b1[i] = b0[i];
                b0[i] = c[j] + 2.0 / (b - a) * (2.0 * x[i] - a - b) * b1[i] - b2[i];
            }
        }

        double[] v = new double[m];

        for (i = 0; i < m; i++)
        {
            v[i] = 0.5 * (c[0] + b0[i] - b2[i]);
        }

        return v;
    }

    public static void t_quadrature_rule(int n, ref double[] t, ref double[] w )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    T_QUADRATURE_RULE: quadrature rule for T(n,x).
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
        bj[0] = Math.Sqrt(0.5);
        for (i = 1; i < n; i++)
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
        
    public static double tt_product(int i, int j, double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TT_PRODUCT: evaluate T(i,x)*T(j,x)
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
        //    Output, double TT_PRODUCT, the value.
        //
    {
        double value;

        if (i < 0 || j < 0)
        {
            value = 0.0;
        }
        else
        {
            int ipj = i + j;
            double tipj = t_polynomial_value(ipj, x);
            int imj = Math.Abs(i - j);
            double timj = t_polynomial_value(imj, x);
            value = 0.5 * (tipj + timj);
        }

        return value;
    }

    public static double tt_product_integral(int i, int j)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TT_PRODUCT_INTEGRAL: integral (-1<=x<=1) T(i,x)*T(j,x)/sqrt(1-x^2) dx
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 April 2012
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
        //    Output, double TT_PRODUCT_INTEGRAL, the value of the integral.
        //
    {
        double value = 0;

        switch (i)
        {
            case < 0:
                Console.WriteLine("");
                Console.WriteLine("TT_PRODUCT_INTEGRAL - Fatal error!");
                Console.WriteLine("  0 <= I, is required.");
                return 1;
        }

        switch (j)
        {
            case < 0:
                Console.WriteLine("");
                Console.WriteLine("TT_PRODUCT_INTEGRAL - Fatal error!");
                Console.WriteLine("  0 <= J is required.");
                return 1;
        }

        if (i != j)
        {
            value = 0.0;
        }
        else
        {
            value = i switch
            {
                0 => Math.PI,
                > 0 => Math.PI / 2.0,
                _ => value
            };
        }

        return value;
    }

    public static double ttt_product_integral(int i, int j, int k)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TTT_PRODUCT_INTEGRAL: integral (-1<=x<=1) T(i,x)*T(j,x)*T(k,x)/sqrt(1-x^2) dx
        //
        //  Discussion:
        //
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
        //  Reference:
        //
        //    John Mason, David Handscomb,
        //    Chebyshev Polynomials,
        //    CRC Press, 2002,
        //    ISBN: 0-8493-035509,
        //    LC: QA404.5.M37.
        //
        //  Parameters:
        //
        //    Input, int I, J, K, the polynomial indices.
        //    0 <= I, J.
        //
        //    Output, double TTT_PRODUCT_INTEGRAL, the integral.
        //
    {
        double value;

        switch (i)
        {
            case < 0:
                value = 0.0;
                break;
            default:
            {
                switch (j)
                {
                    case < 0:
                        value = 0.0;
                        break;
                    default:
                    {
                        value = k switch
                        {
                            < 0 => 0.0,
                            _ => 0.5 * (tt_product_integral(i + j, k) + +tt_product_integral(Math.Abs(i - j), k))
                        };

                        break;
                    }
                }

                break;
            }
        }

        return value;
    }

    public static double tu_product(int i, int j, double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TU_PRODUCT: evaluate T(i,x)*U(j,x)
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
        //    Output, double TU_PRODUCT, the value.
        //
    {
        double value;

        switch (i)
        {
            case < 0:
                value = 0.0;
                break;
            default:
            {
                switch (j)
                {
                    case < 0:
                        value = 0.0;
                        break;
                    default:
                    {
                        value = i switch
                        {
                            0 => u_polynomial_value(j, x),
                            _ => 0.5 * (uu_product(i, j, x) - uu_product(i - 2, j, x))
                        };

                        break;
                    }
                }

                break;
            }
        }

        return value;
    }
}