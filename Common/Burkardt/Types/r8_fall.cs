namespace Burkardt.Types;

public static partial class typeMethods
{
    public static double r8_fall(double x, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_FALL computes the falling factorial function [X]_N.
        //
        //  Discussion:
        //
        //    Note that the number of "injections" or 1-to-1 mappings from
        //    a set of N elements to a set of M elements is [M]_N.
        //
        //    The number of permutations of N objects out of M is [M]_N.
        //
        //    Moreover, the Stirling numbers of the first kind can be used
        //    to convert a falling factorial into a polynomial, as follows:
        //
        //      [X]_N = S^0_N + S^1_N * X + S^2_N * X^2 + ... + S^N_N X^N.
        //
        //    The formula is:
        //
        //      [X]_N = X * ( X - 1 ) * ( X - 2 ) * ... * ( X - N + 1 ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 May 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument of the falling factorial function.
        //
        //    Input, int N, the order of the falling factorial function.
        //    If N = 0, FALL = 1, if N = 1, FALL = X.  Note that if N is
        //    negative, a "rising" factorial will be computed.
        //
        //    Output, double R8_FALL, the value of the falling factorial function.
        //
    {
        int i;

        double value = 1.0;

        switch (n)
        {
            case > 0:
            {
                for (i = 1; i <= n; i++)
                {
                    value *= x;
                    x -= 1.0;
                }

                break;
            }
            case < 0:
            {
                for (i = -1; n <= i; i--)
                {
                    value *= x;
                    x += 1.0;
                }

                break;
            }
        }

        return value;
    }

    public static void r8_fall_values(ref int n_data, ref double x, ref int n, ref double f )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_FALL_VALUES returns some values of the falling factorial function.
        //
        //  Discussion:
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      FactorialPower[X,Y]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 December 2014
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
        //    Output, double &X, int &N, the arguments of the function.
        //
        //    Output, double &F, the value of the function.
        //
    {
        const int N_MAX = 15;

        double[] f_vec =
            {
                120.0000000000000,
                163.1601562500000,
                216.5625000000000,
                281.6601562500000,
                360.0000000000000,
                1.000000000000000,
                7.500000000000000,
                48.75000000000000,
                268.1250000000000,
                1206.562500000000,
                4222.968750000000,
                10557.42187500000,
                15836.13281250000,
                7918.066406250000,
                -3959.03320312500
            }
            ;

        int[] n_vec =
            {
                4,
                4,
                4,
                4,
                4,
                0,
                1,
                2,
                3,
                4,
                5,
                6,
                7,
                8,
                9
            }
            ;

        double[] x_vec =
            {
                5.00,
                5.25,
                5.50,
                5.75,
                6.00,
                7.50,
                7.50,
                7.50,
                7.50,
                7.50,
                7.50,
                7.50,
                7.50,
                7.50,
                7.50
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
            x = 0.0;
            n = 0;
            f = 0.0;
        }
        else
        {
            x = x_vec[n_data - 1];
            n = n_vec[n_data - 1];
            f = f_vec[n_data - 1];
        }

    }
}