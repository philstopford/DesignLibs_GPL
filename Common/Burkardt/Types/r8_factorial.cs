using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static double r8_factorial(int n)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_FACTORIAL computes the factorial of N.
            //
            //  Discussion:
            //
            //    factorial ( N ) = product ( 1 <= I <= N ) I
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 January 1999
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the argument of the factorial function.
            //    If N is less than 1, the function value is returned as 1.
            //
            //    Output, double R8_FACTORIAL, the factorial of N.
            //
        {
            double value = 1.0;

            for (int i = 1; i <= n; i++)
            {
                value = value * (double) (i);
            }

            return value;
        }

        public static double r8_factorial_stirling ( int n )
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_FACTORIAL_STIRLING computes Stirling's approximation to N!.
            //
            //  Discussion:
            //
            //    N! = Product ( 1 <= I <= N ) I
            //
            //    Stirling ( N ) = sqrt ( 2 * PI * N ) * ( N / E )^N * E^(1/(12*N) )
            //
            //    This routine returns the raw approximation for all nonnegative
            //    values of N.  If N is less than 0, the value is returned as 0,
            //    and if N is 0, the value of 1 is returned.  In all other cases,
            //    Stirling's formula is used.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 April 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the argument of the function.
            //
            //    Output, double R8_FACTORIAL_STIRLING, an approximation to N!.
            //
        {
            double value;

            if (n < 0)
            {
                value = 0.0;
            }
            else if (n == 0)
            {
                value = 1.0;
            }
            else
            {
                value = Math.Sqrt(2.0 * Math.PI * (double) (n))
                        * Math.Pow((double) (n) / Math.E, n)
                        * Math.Exp(1.0 / (double) (12 * n));
            }

            return value;
        }

        public static void r8_factorial_values(ref int n_data, ref int n, ref double fn )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_FACTORIAL_VALUES returns values of the real factorial function.
        //
        //  Discussion:
        //
        //    0! = 1
        //    I! = Product ( 1 <= J <= I ) J
        //
        //    Although the factorial is an int *valued function, it quickly
        //    becomes too large for an int *to hold.  This routine still accepts
        //    an int *as the input argument, but returns the function value
        //    as a real number.
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      n!
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 August 2004
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
        //    Output, int &N, the argument of the function.
        //
        //    Output, double &FN, the value of the function.
        //
        {
            int N_MAX = 25;

            double[] fn_vec =
            {
                0.1000000000000000E+01,
                0.1000000000000000E+01,
                0.2000000000000000E+01,
                0.6000000000000000E+01,
                0.2400000000000000E+02,
                0.1200000000000000E+03,
                0.7200000000000000E+03,
                0.5040000000000000E+04,
                0.4032000000000000E+05,
                0.3628800000000000E+06,
                0.3628800000000000E+07,
                0.3991680000000000E+08,
                0.4790016000000000E+09,
                0.6227020800000000E+10,
                0.8717829120000000E+11,
                0.1307674368000000E+13,
                0.2092278988800000E+14,
                0.3556874280960000E+15,
                0.6402373705728000E+16,
                0.1216451004088320E+18,
                0.2432902008176640E+19,
                0.1551121004333099E+26,
                0.3041409320171338E+65,
                0.9332621544394415E+158,
                0.5713383956445855E+263
            }
            ;

            int[] n_vec =
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
                10,
                11,
                12,
                13,
                14,
                15,
                16,
                17,
                18,
                19,
                20,
                25,
                50,
                100,
                150
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
                fn = 0.0;
            }
            else
            {
                n = n_vec[n_data - 1];
                fn = fn_vec[n_data - 1];
            }
        }

        public static double r8_factorial2(int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_FACTORIAL2 computes the double factorial function.
            //
            //  Discussion:
            //
            //    FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
            //                    = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
            //
            //  Example:
            //
            //     N Value
            //
            //     0     1
            //     1     1
            //     2     2
            //     3     3
            //     4     8
            //     5    15
            //     6    48
            //     7   105
            //     8   384
            //     9   945
            //    10  3840
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
            //    Input, int N, the argument of the double factorial
            //    function.  If N is less than 1, R8_FACTORIAL2 is returned as 1.0.
            //
            //    Output, double R8_FACTORIAL2, the value of Factorial2(N).
            //
        {
            int n_copy;
            double value;

            value = 1.0;

            if (n < 1)
            {
                return value;
            }

            n_copy = n;

            while (1 < n_copy)
            {
                value = value * (double) n_copy;
                n_copy = n_copy - 2;
            }

            return value;
        }

        public static void r8_factorial2_values(ref int n_data, ref int n, ref double f )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_FACTORIAL2_VALUES returns values of the double factorial function.
        //
        //  Formula:
        //
        //    FACTORIAL2( N ) = Product ( N * (N-2) * (N-4) * ... * 2 )  (N even)
        //                    = Product ( N * (N-2) * (N-4) * ... * 1 )  (N odd)
        //
        //    In Mathematica, the function can be evaluated by:
        //
        //      n!!
        //
        //  Example:
        //
        //     N    N!!
        //
        //     0     1
        //     1     1
        //     2     2
        //     3     3
        //     4     8
        //     5    15
        //     6    48
        //     7   105
        //     8   384
        //     9   945
        //    10  3840
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 February 2015
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
        //    Daniel Zwillinger,
        //    CRC Standard Mathematical Tables and Formulae,
        //    30th Edition,
        //    CRC Press, 1996, page 16.
        //
        //  Parameters:
        //
        //    Input/output, int &N_DATA.  The user sets N_DATA to 0 before the
        //    first call.  On each call, the routine increments N_DATA by 1, and
        //    returns the corresponding data; when there is no more data, the
        //    output value of N_DATA will be 0 again.
        //
        //    Output, int &N, the argument of the function.
        //
        //    Output, double &F, the value of the function.
        //
        {
            int N_MAX = 16;

            double[] f_vec =
            {
                1.0,
                1.0,
                2.0,
                3.0,
                8.0,
                15.0,
                48.0,
                105.0,
                384.0,
                945.0,
                3840.0,
                10395.0,
                46080.0,
                135135.0,
                645120.0,
                2027025.0
            }
            ;

            int[] n_vec =
            {
                0,
                1, 2, 3, 4, 5,
                6, 7, 8, 9, 10,
                11, 12, 13, 14, 15
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
                f = 0.0;
            }
            else
            {
                n = n_vec[n_data - 1];
                f = f_vec[n_data - 1];
            }
        }

    }
}