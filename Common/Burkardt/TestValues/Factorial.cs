namespace Burkardt.TestValues
{
    public static class Factorial
    {
        public static void i4_factorial_values(ref int n_data, ref int n, ref int fn)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4_FACTORIAL_VALUES returns values of the factorial function.
            //
            //  Discussion:
            //
            //    0! = 1
            //    I! = Product ( 1 <= J <= I ) I
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
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref int N, the argument of the function.
            //
            //    Output, ref int FN, the value of the function.
            //
        {
            int N_MAX = 13;

            int[] fn_vec =
            {
                1,
                1,
                2,
                6,
                24,
                120,
                720,
                5040,
                40320,
                362880,
                3628800,
                39916800,
                479001600
            };

            int[] n_vec =
            {
                0, 1, 2, 3,
                4, 5, 6, 7,
                8, 9, 10, 11,
                12
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                n = 0;
                fn = 0;
            }
            else
            {
                n = n_vec[n_data - 1];
                fn = fn_vec[n_data - 1];
            }
        }

        public static void i4_factorial2_values(ref int n_data, ref int n, ref int fn)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4_FACTORIAL2_VALUES returns values of the double factorial function.
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
            //    Daniel Zwillinger,
            //    CRC Standard Mathematical Tables and Formulae,
            //    30th Edition,
            //    CRC Press, 1996, page 16.
            //
            //  Parameters:
            //
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref int N, the argument of the function.
            //
            //    Output, ref int FN, the value of the function.
            //
        {
            int N_MAX = 16;

            int[] fn_vec =
            {
                1,
                1,
                2,
                3,
                8,
                15,
                48,
                105,
                384,
                945,
                3840,
                10395,
                46080,
                135135,
                645120,
                2027025
            };

            int[] n_vec =
            {
                0,
                1, 2, 3, 4, 5,
                6, 7, 8, 9, 10,
                11, 12, 13, 14, 15
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                n = 0;
                fn = 0;
            }
            else
            {
                n = n_vec[n_data - 1];
                fn = fn_vec[n_data - 1];
            }
        }

        public static void i4_fall_values(ref int n_data, ref int m, ref int n, ref int fmn)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4_FALL_VALUES returns values of the integer falling factorial function.
            //
            //  Discussion:
            //
            //    The definition of the falling factorial function is
            //
            //      (m)_n = (m)! / (m-n)!
            //            = ( m ) * ( m - 1 ) * ( m - 2 ) ... * ( m - n + 1 )
            //            = Gamma ( m + 1 ) / Gamma ( m - n + 1 )
            //
            //    We assume 0 <= N <= M.
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      FactorialPower[m,n]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 December 2014
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
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref int M, &N, the arguments of the function.
            //
            //    Output, ref int FMN, the value of the function.
            //
        {
            int N_MAX = 15;

            int[] fmn_vec =
            {
                1, 5, 20, 60, 120,
                120, 0, 1, 10, 4000,
                90, 4896, 24, 912576, 0
            };

            int[] m_vec =
            {
                5, 5, 5, 5, 5,
                5, 5, 50, 10, 4000,
                10, 18, 4, 98, 1
            };

            int[] n_vec =
            {
                0, 1, 2, 3, 4,
                5, 6, 0, 1, 1,
                2, 3, 4, 3, 7
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                m = 0;
                n = 0;
                fmn = 0;
            }
            else
            {
                m = m_vec[n_data - 1];
                n = n_vec[n_data - 1];
                fmn = fmn_vec[n_data - 1];
            }
        }

        public static void i4_rise_values(ref int n_data, ref int m, ref int n, ref int fmn)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4_RISE_VALUES returns values of the integer rising factorial function.
            //
            //  Discussion:
            //
            //    The integer rising factorial function is sometimes symbolized by (m)_n.
            //
            //    The definition is
            //
            //      (m)_n = (m-1+n)! / (m-1)!
            //            = ( m ) * ( m + 1 ) * ( m + 2 ) ... * ( m - 1 + n )
            //            = Gamma ( m + n ) / Gamma ( m )
            //
            //    We assume 0 <= N <= M.
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      Pochhammer[m,n]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 December 2014
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
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref int M, &N, the arguments of the function.
            //
            //    Output, ref int FMN, the value of the function.
            //
        {
            int N_MAX = 15;

            int[] fmn_vec =
            {
                1, 5, 30, 210, 1680,
                15120, 151200, 1, 10, 4000,
                110, 6840, 840, 970200, 5040
            };

            int[] m_vec =
            {
                5, 5, 5, 5, 5,
                5, 5, 50, 10, 4000,
                10, 18, 4, 98, 1
            };

            int[] n_vec =
            {
                0, 1, 2, 3, 4,
                5, 6, 0, 1, 1,
                2, 3, 4, 3, 7
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                m = 0;
                n = 0;
                fmn = 0;
            }
            else
            {
                m = m_vec[n_data - 1];
                n = n_vec[n_data - 1];
                fmn = fmn_vec[n_data - 1];
            }
        }

        public static void r8_factorial_values(ref int n_data, ref int n, ref double fn)

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
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref int N, the argument of the function.
            //
            //    Output, ref double FN, the value of the function.
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
            };

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
            };

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

        public static void r8_factorial_log_values(ref int n_data, ref int n, ref double fn)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_FACTORIAL_LOG_VALUES returns values of log(n!).
            //
            //  Discussion:
            //
            //    The function log(n!) can be written as
            //
            //     log(n!) = sum ( 1 <= i <= n ) log ( i )
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      Log[n!]
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
            //    Stephen Wolfram,
            //    The Mathematica Book,
            //    Fourth Edition,
            //    Cambridge University Press, 1999,
            //    ISBN: 0-521-64314-7,
            //    LC: QA76.95.W65.
            //
            //    Daniel Zwillinger, editor,
            //    CRC Standard Mathematical Tables and Formulae,
            //    30th Edition,
            //    CRC Press, 1996.
            //
            //  Parameters:
            //
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref int N, the argument of the function.
            //
            //    Output, ref double FN, the value of the function.
            //
        {
            int N_MAX = 27;

            double[] fn_vec =
            {
                0.0000000000000000E+00,
                0.0000000000000000E+00,
                0.6931471805599453E+00,
                0.1791759469228055E+01,
                0.3178053830347946E+01,
                0.4787491742782046E+01,
                0.6579251212010101E+01,
                0.8525161361065414E+01,
                0.1060460290274525E+02,
                0.1280182748008147E+02,
                0.1510441257307552E+02,
                0.1750230784587389E+02,
                0.1998721449566189E+02,
                0.2255216385312342E+02,
                0.2519122118273868E+02,
                0.2789927138384089E+02,
                0.3067186010608067E+02,
                0.3350507345013689E+02,
                0.3639544520803305E+02,
                0.3933988418719949E+02,
                0.4233561646075349E+02,
                0.5800360522298052E+02,
                0.1484777669517730E+03,
                0.3637393755555635E+03,
                0.6050201058494237E+03,
                0.2611330458460156E+04,
                0.5912128178488163E+04
            };

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
                150,
                500,
                1000
            };

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

        public static void r8_factorial2_values(ref int n_data, ref int n, ref double f)

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
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref int N, the argument of the function.
            //
            //    Output, ref double F, the value of the function.
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
            };

            int[] n_vec =
            {
                0,
                1, 2, 3, 4, 5,
                6, 7, 8, 9, 10,
                11, 12, 13, 14, 15
            };

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

        public static void r8_fall_values(ref int n_data, ref double x, ref int n, ref double f)

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
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref double X, ref int N, the arguments of the function.
            //
            //    Output, ref double F, the value of the function.
            //
        {
            int N_MAX = 15;

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
            };

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
            };

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
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

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

        public static void r8_rise_values(ref int n_data, ref double x, ref int n, ref double f)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_RISE_VALUES returns some values of the rising factorial function.
            //
            //  Discussion:
            //
            //    Pochhammer(X,Y) = Gamma(X+Y) / Gamma(X)
            //
            //    For integer arguments, Pochhammer(M,N) = ( M + N - 1 )! / ( N - 1 )!
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      Pochhammer[X,Y]
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
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref double X, ref int N, the arguments of the function.
            //
            //    Output, ref double F, the value of the function.
            //
        {
            int N_MAX = 15;

            double[] f_vec =
            {
                1680.000000000000,
                1962.597656250000,
                2279.062500000000,
                2631.972656250000,
                3024.000000000000,
                1.000000000000000,
                7.500000000000000,
                63.75000000000000,
                605.6250000000000,
                6359.062500000000,
                73129.21875000000,
                914115.2343750000,
                1.234055566406250E+07,
                1.789380571289063E+08,
                2.773539885498047E+09
            };

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
            };

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
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

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


        public static void subfactorial_values(ref int n_data, ref int n, ref int fn)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SUBFACTORIAL_VALUES returns values of the subfactorial function.
            //
            //  Discussion:
            //
            //    The subfactorial function Subfactorial(N) counts the number of
            //    permutations of N objects which leave no object unchanged.
            //
            //    Such a permutation is known as a derangement.
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      << DiscreteMath`CombinatorialFunctions`
            //      Subfactorial[n]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    21 March 2007
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
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref int N, the argument of the function.
            //
            //    Output, ref int FN, the value of the function.
            //
        {
            int N_MAX = 13;

            int[] fn_vec =
            {
                1,
                0,
                1,
                2,
                9,
                44,
                265,
                1854,
                14833,
                133496,
                1334961,
                14684570,
                176214841
            };

            int[] n_vec =
            {
                0, 1, 2, 3,
                4, 5, 6, 7,
                8, 9, 10, 11,
                12
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                n = 0;
                fn = 0;
            }
            else
            {
                n = n_vec[n_data - 1];
                fn = fn_vec[n_data - 1];
            }
        }


    }
}