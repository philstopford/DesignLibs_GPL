namespace TestValues
{
    public static class Sine
    {
        public static void si_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SI_VALUES returns some values of the sine integral function.
            //
            //  Discussion:
            //
            //    SI(X) = integral ( 0 <= T <= X ) sin ( T ) / T dt
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      SinIntegral[x]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    12 August 2004
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
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 16;

            double[] fx_vec =
            {
                0.4931074180430667E+00,
                0.5881288096080801E+00,
                0.6812222391166113E+00,
                0.7720957854819966E+00,
                0.8604707107452929E+00,
                0.9460830703671830E+00,
                0.1108047199013719E+01,
                0.1256226732779218E+01,
                0.1389180485870438E+01,
                0.1505816780255579E+01,
                0.1605412976802695E+01,
                0.1778520173443827E+01,
                0.1848652527999468E+01,
                0.1833125398665997E+01,
                0.1758203138949053E+01,
                0.1654140414379244E+01
            };

            double[] x_vec =
            {
                0.5E+00,
                0.6E+00,
                0.7E+00,
                0.8E+00,
                0.9E+00,
                1.0E+00,
                1.2E+00,
                1.4E+00,
                1.6E+00,
                1.8E+00,
                2.0E+00,
                2.5E+00,
                3.0E+00,
                3.5E+00,
                4.0E+00,
                4.5E+00
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
                fx = 0.0;
            }
            else
            {
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void sin_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SIN_VALUES returns some values of the sine function.
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      Sin[x]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    12 June 2007
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
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 13;

            double[] fx_vec =
            {
                0.00000000000000000000,
                0.25881904510252076235,
                0.47942553860420300027,
                0.50000000000000000000,
                0.70710678118654752440,
                0.84147098480789650665,
                0.86602540378443864676,
                1.00000000000000000000,
                0.90929742682568169540,
                0.14112000805986722210,
                0.00000000000000000000,
                -0.75680249530792825137,
                -0.95892427466313846889
            };

            double[] x_vec =
            {
                0.0000000000000000000,
                0.26179938779914943654,
                0.50000000000000000000,
                0.52359877559829887308,
                0.78539816339744830962,
                1.0000000000000000000,
                1.0471975511965977462,
                1.5707963267948966192,
                2.0000000000000000000,
                3.0000000000000000000,
                3.1415926535897932385,
                4.0000000000000000000,
                5.0000000000000000000
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
                fx = 0.0;
            }
            else
            {
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void sin_degree_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SIN_DEGREE_VALUES: the sine function with argument in degrees.
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      Sin[x Degree]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 March 2010
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
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 22;

            double[] fx_vec =
            {
                -0.087155742747658173558,
                0.000000000000000000000,
                0.017452406437283512819,
                0.034899496702500971646,
                0.052335956242943832722,
                0.069756473744125300776,
                0.087155742747658173558,
                0.17364817766693034885,
                0.25881904510252076235,
                0.50000000000000000000,
                0.70710678118654752440,
                0.86602540378443864676,
                0.96592582628906828675,
                0.99619469809174553230,
                0.99756405025982424761,
                0.99862953475457387378,
                0.99939082701909573001,
                0.99984769515639123916,
                1.0000000000000000000,
                0.99984769515639123916,
                0.96592582628906828675,
                0.00000000000000000000
            };
            double[] x_vec =
            {
                -5.0,
                0.0,
                1.0,
                2.0,
                3.0,
                4.0,
                5.0,
                10.0,
                15.0,
                30.0,
                45.0,
                60.0,
                75.0,
                85.0,
                86.0,
                87.0,
                88.0,
                89.0,
                90.0,
                91.0,
                105.0,
                180.0
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
                fx = 0.0;
            }
            else
            {
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void sin_power_int_values(ref int n_data, ref double a, ref double b, ref int n,
                ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SIN_POWER_INT_VALUES returns some values of the sine power integral.
            //
            //  Discussion:
            //
            //    The function has the form
            //
            //      SIN_POWER_INT(A,B,N) = Integral ( A <= T <= B ) ( sin(T) )^N dt
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      Integrate [ ( Sin[x] )^n, { x, a, b } ]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    02 September 2004
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
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref double A, &B, the limits of integration.
            //
            //    Output, ref int N, the power.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 10;

            double[] a_vec =
            {
                0.10E+02,
                0.00E+00,
                0.00E+00,
                0.00E+00,
                0.00E+00,
                0.00E+00,
                0.00E+00,
                0.10E+01,
                0.00E+00,
                0.00E+00
            };

            double[] b_vec =
            {
                0.20E+02,
                0.10E+01,
                0.10E+01,
                0.10E+01,
                0.10E+01,
                0.10E+01,
                0.20E+01,
                0.20E+01,
                0.10E+01,
                0.10E+01
            };

            double[] fx_vec =
            {
                0.10000000000000000000E+02,
                0.45969769413186028260E+00,
                0.27267564329357957615E+00,
                0.17894056254885809051E+00,
                0.12402556531520681830E+00,
                0.88974396451575946519E-01,
                0.90393123848149944133E+00,
                0.81495684202992349481E+00,
                0.21887522421729849008E-01,
                0.17023439374069324596E-01
            };

            int[] n_vec =
            {
                0,
                1,
                2,
                3,
                4,
                5,
                5,
                5,
                10,
                11
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

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

        public static void sinh_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SINH_VALUES returns some values of the hyperbolic sine function.
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      Sinh[x]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 June 2007
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
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 18;

            double[] fx_vec =
            {
                -74.203210577788758977,
                -1.1752011936438014569,
                0.00000000000000000000,
                0.10016675001984402582,
                0.20133600254109398763,
                0.30452029344714261896,
                0.41075232580281550854,
                0.52109530549374736162,
                0.63665358214824127112,
                0.75858370183953350346,
                0.88810598218762300657,
                1.0265167257081752760,
                1.1752011936438014569,
                3.6268604078470187677,
                10.017874927409901899,
                27.289917197127752449,
                74.203210577788758977,
                11013.232874703393377
            };

            double[] x_vec =
            {
                -5.0,
                -1.0,
                0.0,
                0.1,
                0.2,
                0.3,
                0.4,
                0.5,
                0.6,
                0.7,
                0.8,
                0.9,
                1.0,
                2.0,
                3.0,
                4.0,
                5.0,
                10.0
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
                fx = 0.0;
            }
            else
            {
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

    }
}