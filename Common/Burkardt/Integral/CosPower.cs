using System;

namespace Burkardt.IntegralNS;

public static class CosPower
{
    public static double cos_power_int(double a, double b, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COS_POWER_INT evaluates the cosine power integral.
        //
        //  Discussion:
        //
        //    The function is defined by
        //
        //      COS_POWER_INT(A,B,N) = Integral ( A <= T <= B ) ( cos ( t ))^n dt
        //
        //    The algorithm uses the following fact:
        //
        //      Integral cos^n ( t ) = -(1/n) * (
        //        cos^(n-1)(t) * sin(t) + ( n-1 ) * Integral cos^(n-2) ( t ) dt )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 March 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters
        //
        //    Input, double A, B, the limits of integration.
        //
        //    Input, integer N, the power of the sine function.
        //
        //    Output, double COS_POWER_INT, the value of the integral.
        //
    {
        int m;
        int mlo;
        double value;

        switch (n)
        {
            case < 0:
                Console.WriteLine("");
                Console.WriteLine("COS_POWER_INT - Fatal error!");
                Console.WriteLine("  Power N < 0.");
                return 1;
        }

        double sa = Math.Sin(a);
        double sb = Math.Sin(b);
        double ca = Math.Cos(a);
        double cb = Math.Cos(b);

        switch (n % 2)
        {
            case 0:
                value = b - a;
                mlo = 2;
                break;
            default:
                value = sb - sa;
                mlo = 3;
                break;
        }

        for (m = mlo; m <= n; m += 2)
        {
            value = ((m - 1) * value
                        - Math.Pow(ca, m - 1) * sa + Math.Pow(cb, m - 1) * sb)
                    / m;
        }

        return value;
    }

    public static void cos_power_int_values(ref int n_data, ref double a, ref double b, ref int n,
            ref double fx )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COS_POWER_INT_VALUES returns some values of the sine power integral.
        //
        //  Discussion:
        //
        //    The function has the form
        //
        //      COS_POWER_INT(A,B,N) = Integral ( A <= T <= B ) ( cos(T) )^N dt
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      Integrate [ ( Cos[x] )^n, { x, a, b } ]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 March 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
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
        //    Output, double &A, &B, the limits of integration.
        //
        //    Output, int &N, the power.
        //
        //    Output, double &FX, the value of the function.
        //
    {
        const int N_MAX = 11;

        double[] a_vec
                =
                {
                    0.00E+00,
                    0.00E+00,
                    0.00E+00,
                    0.00E+00,
                    0.00E+00,
                    0.00E+00,
                    0.00E+00,
                    0.00E+00,
                    0.00E+00,
                    0.00E+00,
                    0.00E+00
                }
            ;

        double[] b_vec
                =
                {
                    3.141592653589793,
                    3.141592653589793,
                    3.141592653589793,
                    3.141592653589793,
                    3.141592653589793,
                    3.141592653589793,
                    3.141592653589793,
                    3.141592653589793,
                    3.141592653589793,
                    3.141592653589793,
                    3.141592653589793
                }
            ;

        double[] fx_vec
                =
                {
                    3.141592653589793,
                    0.0,
                    1.570796326794897,
                    0.0,
                    1.178097245096172,
                    0.0,
                    0.9817477042468104,
                    0.0,
                    0.8590292412159591,
                    0.0,
                    0.7731263170943632
                }
            ;

        int[] n_vec
                =
                {
                    0,
                    1,
                    2,
                    3,
                    4,
                    5,
                    6,
                    7,
                    8,
                    9,
                    10
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
            a = 0.0;
            b = 0.0;
            n = 0;
            fx = 0.0;
        }
        else
        {
            a = a_vec[n_data - 1];
            b = b_vec[n_data - 1];
            n = n_vec[n_data - 1];
            fx = fx_vec[n_data - 1];
        }
    }
}