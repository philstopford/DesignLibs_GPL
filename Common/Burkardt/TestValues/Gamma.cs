namespace Burkardt.TestValues
{
    public static class Gamma
    {

        public static void gamma_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GAMMA_VALUES returns some values of the Gamma function.
            //
            //  Discussion:
            //
            //    The Gamma function is defined as:
            //
            //      Gamma(Z) = Integral ( 0 <= T < +oo ) T^(Z-1) exp(-T) dT
            //
            //    It satisfies the recursion:
            //
            //      Gamma(X+1) = X * Gamma(X)
            //
            //    Gamma is undefined for nonpositive integral X.
            //    Gamma(0.5) = sqrt(PI)
            //    For N a positive integer, Gamma(N+1) = N!, the standard factorial.
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      Gamma[x]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    20 May 2007
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
            int N_MAX = 25;

            double[] fx_vec =
            {
                -0.3544907701811032E+01,
                -0.1005871979644108E+03,
                0.9943258511915060E+02,
                0.9513507698668732E+01,
                0.4590843711998803E+01,
                0.2218159543757688E+01,
                0.1772453850905516E+01,
                0.1489192248812817E+01,
                0.1164229713725303E+01,
                0.1000000000000000E+01,
                0.9513507698668732E+00,
                0.9181687423997606E+00,
                0.8974706963062772E+00,
                0.8872638175030753E+00,
                0.8862269254527580E+00,
                0.8935153492876903E+00,
                0.9086387328532904E+00,
                0.9313837709802427E+00,
                0.9617658319073874E+00,
                0.1000000000000000E+01,
                0.2000000000000000E+01,
                0.6000000000000000E+01,
                0.3628800000000000E+06,
                0.1216451004088320E+18,
                0.8841761993739702E+31
            };

            double[] x_vec =
            {
                -0.50E+00,
                -0.01E+00,
                0.01E+00,
                0.10E+00,
                0.20E+00,
                0.40E+00,
                0.50E+00,
                0.60E+00,
                0.80E+00,
                1.00E+00,
                1.10E+00,
                1.20E+00,
                1.30E+00,
                1.40E+00,
                1.50E+00,
                1.60E+00,
                1.70E+00,
                1.80E+00,
                1.90E+00,
                2.00E+00,
                3.00E+00,
                4.00E+00,
                10.00E+00,
                20.00E+00,
                30.00E+00
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

        public static void gamma_01_pdf_values(ref int n_data, ref double alpha, ref double x,
                ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GAMMA_01_PDF_VALUES returns some values of the standard Gamma PDF.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    28 July 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref double ALPHA, the shape parameter.
            //
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 10;

            double[] alpha_vec =
            {
                1.092091484911879,
                4.147546169663503,
                2.076535407379806,
                1.287888961910225,
                0.2191449888955355,
                0.3086361453280091,
                2.006531407488083,
                3.986434770531281,
                4.487520304498656,
                0.472723751058401
            };

            double[] fx_vec =
            {
                0.00009260811963612823,
                0.1260335478747823,
                0.1363536772414351,
                0.5114450139194701,
                0.0001230139468263628,
                0.001870342832511005,
                0.004476000451227789,
                0.0,
                0.2056668486524041,
                0.0
            };

            double[] x_vec =
            {
                9.541334553343761,
                5.39780214905239,
                0.1942467166183289,
                0.6545463320909413,
                6.156639979175331,
                4.220159083225351,
                7.424071607424807,
                -0.4806971028367454,
                3.18289954879574,
                -0.3570226383736496
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                alpha = 0.0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                alpha = alpha_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void gamma_cdf_values(ref int n_data, ref double mu, ref double sigma, ref double x,
                ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GAMMA_CDF_VALUES returns some values of the Gamma CDF.
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      Needs["Statistics`ContinuousDistributions`"]
            //      dist = GammaDistribution [ mu, sigma ]
            //      CDF [ dist, x ]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 September 2004
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
            //    Output, ref double MU, the mean of the distribution.
            //
            //    Output, ref double SIGMA, the standard deviation of the distribution.
            //
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 12;

            double[] fx_vec =
            {
                0.8646647167633873E+00,
                0.9816843611112658E+00,
                0.9975212478233336E+00,
                0.9996645373720975E+00,
                0.6321205588285577E+00,
                0.4865828809674080E+00,
                0.3934693402873666E+00,
                0.3296799539643607E+00,
                0.4421745996289254E+00,
                0.1911531694619419E+00,
                0.6564245437845009E-01,
                0.1857593622214067E-01
            };

            double[] mu_vec =
            {
                0.1000000000000000E+01,
                0.1000000000000000E+01,
                0.1000000000000000E+01,
                0.1000000000000000E+01,
                0.1000000000000000E+01,
                0.1000000000000000E+01,
                0.1000000000000000E+01,
                0.1000000000000000E+01,
                0.2000000000000000E+01,
                0.3000000000000000E+01,
                0.4000000000000000E+01,
                0.5000000000000000E+01
            };

            double[] sigma_vec =
            {
                0.5000000000000000E+00,
                0.5000000000000000E+00,
                0.5000000000000000E+00,
                0.5000000000000000E+00,
                0.2000000000000000E+01,
                0.3000000000000000E+01,
                0.4000000000000000E+01,
                0.5000000000000000E+01,
                0.2000000000000000E+01,
                0.2000000000000000E+01,
                0.2000000000000000E+01,
                0.2000000000000000E+01
            };

            double[] x_vec =
            {
                0.1000000000000000E+01,
                0.2000000000000000E+01,
                0.3000000000000000E+01,
                0.4000000000000000E+01,
                0.2000000000000000E+01,
                0.2000000000000000E+01,
                0.2000000000000000E+01,
                0.2000000000000000E+01,
                0.3000000000000000E+01,
                0.3000000000000000E+01,
                0.3000000000000000E+01,
                0.3000000000000000E+01
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                mu = 0.0;
                sigma = 0.0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                mu = mu_vec[n_data - 1];
                sigma = sigma_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void gamma_inc_p_values(ref int n_data, ref double a, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GAMMA_INC_P_VALUES: values of normalized incomplete Gamma function P(A,X).
            //
            //  Discussion:
            //
            //    The (normalized) incomplete Gamma function is defined as:
            //
            //      P(A,X) = 1/Gamma(A) * Integral ( 0 <= T <= X ) T^(A-1) * exp(-T) dT.
            //
            //    With this definition, for all A and X,
            //
            //      0 <= P(A,X) <= 1
            //
            //    and
            //
            //      P(A,oo) = 1.0
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      1 - GammaRegularized[A,X]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 April 2010
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
            //    Output, ref double A, the parameter of the function.
            //
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 20;

            double[] a_vec =
            {
                0.10E+00,
                0.10E+00,
                0.10E+00,
                0.50E+00,
                0.50E+00,
                0.50E+00,
                0.10E+01,
                0.10E+01,
                0.10E+01,
                0.11E+01,
                0.11E+01,
                0.11E+01,
                0.20E+01,
                0.20E+01,
                0.20E+01,
                0.60E+01,
                0.60E+01,
                0.11E+02,
                0.26E+02,
                0.41E+02
            };

            double[] fx_vec =
            {
                0.7382350532339351E+00,
                0.9083579897300343E+00,
                0.9886559833621947E+00,
                0.3014646416966613E+00,
                0.7793286380801532E+00,
                0.9918490284064973E+00,
                0.9516258196404043E-01,
                0.6321205588285577E+00,
                0.9932620530009145E+00,
                0.7205974576054322E-01,
                0.5891809618706485E+00,
                0.9915368159845525E+00,
                0.1018582711118352E-01,
                0.4421745996289254E+00,
                0.9927049442755639E+00,
                0.4202103819530612E-01,
                0.9796589705830716E+00,
                0.9226039842296429E+00,
                0.4470785799755852E+00,
                0.7444549220718699E+00
            };

            double[] x_vec =
            {
                0.30E-01,
                0.30E+00,
                0.15E+01,
                0.75E-01,
                0.75E+00,
                0.35E+01,
                0.10E+00,
                0.10E+01,
                0.50E+01,
                0.10E+00,
                0.10E+01,
                0.50E+01,
                0.15E+00,
                0.15E+01,
                0.70E+01,
                0.25E+01,
                0.12E+02,
                0.16E+02,
                0.25E+02,
                0.45E+02
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
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                a = a_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void gamma_inc_q_values(ref int n_data, ref double a, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GAMMA_INC_Q_VALUES: values of normalized incomplete Gamma function Q(A,X).
            //
            //  Discussion:
            //
            //    The (normalized) incomplete Gamma function is defined as:
            //
            //      Q(A,X) = 1/Gamma(A) * Integral ( X <= T < oo ) T^(A-1) * exp(-T) dT.
            //
            //    With this definition, for all A and X,
            //
            //      0 <= Q(A,X) <= 1
            //
            //    and
            //
            //      Q(A,oo) = 0.0
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      GammaRegularized[A,X]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 April 2010
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
            //    Output, ref double A, the parameter of the function.
            //
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 20;

            double[] a_vec =
            {
                0.10E+00,
                0.10E+00,
                0.10E+00,
                0.50E+00,
                0.50E+00,
                0.50E+00,
                0.10E+01,
                0.10E+01,
                0.10E+01,
                0.11E+01,
                0.11E+01,
                0.11E+01,
                0.20E+01,
                0.20E+01,
                0.20E+01,
                0.60E+01,
                0.60E+01,
                0.11E+02,
                0.26E+02,
                0.41E+02
            };

            double[] fx_vec =
            {
                0.2617649467660649E+00,
                0.09164201026996572E+00,
                0.01134401663780527E+00,
                0.6985353583033387E+00,
                0.2206713619198468E+00,
                0.008150971593502700E+00,
                0.9048374180359596E+00,
                0.3678794411714423E+00,
                0.006737946999085467E+00,
                0.9279402542394568E+00,
                0.4108190381293515E+00,
                0.008463184015447498E+00,
                0.9898141728888165E+00,
                0.5578254003710746E+00,
                0.007295055724436130E+00,
                0.9579789618046939E+00,
                0.02034102941692837E+00,
                0.07739601577035708E+00,
                0.5529214200244148E+00,
                0.2555450779281301E+00
            };

            double[] x_vec =
            {
                0.30E-01,
                0.30E+00,
                0.15E+01,
                0.75E-01,
                0.75E+00,
                0.35E+01,
                0.10E+00,
                0.10E+01,
                0.50E+01,
                0.10E+00,
                0.10E+01,
                0.50E+01,
                0.15E+00,
                0.15E+01,
                0.70E+01,
                0.25E+01,
                0.12E+02,
                0.16E+02,
                0.25E+02,
                0.45E+02
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
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                a = a_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void gamma_inc_tricomi_values(ref int n_data, ref double a, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GAMMA_INC_TRICOMI_VALUES: values of Tricomi's incomplete Gamma function.
            //
            //  Discussion:
            //
            //    Tricomi's incomplete Gamma function is defined as:
            //
            //      1/Gamma(A) * 1/X^A * Integral ( 0 <= T <= X ) T^(A-1) * exp(-T) dT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 April 2010
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
            //    Output, ref double A, the parameter of the function.
            //
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 20;

            double[] a_vec =
            {
                0.10E+00,
                0.10E+00,
                0.10E+00,
                0.50E+00,
                0.50E+00,
                0.50E+00,
                0.10E+01,
                0.10E+01,
                0.10E+01,
                0.11E+01,
                0.11E+01,
                0.11E+01,
                0.20E+01,
                0.20E+01,
                0.20E+01,
                0.60E+01,
                0.60E+01,
                0.11E+02,
                0.26E+02,
                0.41E+02
            };

            double[] fx_vec =
            {
                1.048292641463504E+00,
                1.024577737369574E+00,
                0.9493712443185374E+00,
                1.100793230316492E+00,
                0.8998911979655218E+00,
                0.5301656062431039E+00,
                0.9516258196404043E+00,
                0.6321205588285577E+00,
                0.1986524106001829E+00,
                0.9071784510537487E+00,
                0.5891809618706485E+00,
                0.1688269752193589E+00,
                0.4527034271637121E+00,
                0.1965220442795224E+00,
                0.02025928457705232E+00,
                0.0001721181724479739E+00,
                3.280858070850586E-07,
                5.244396471821590E-14,
                2.013462926183376E-37,
                1.230623887499875E-68
            };

            double[] x_vec =
            {
                0.30E-01,
                0.30E+00,
                0.15E+01,
                0.75E-01,
                0.75E+00,
                0.35E+01,
                0.10E+00,
                0.10E+01,
                0.50E+01,
                0.10E+00,
                0.10E+01,
                0.50E+01,
                0.15E+00,
                0.15E+01,
                0.70E+01,
                0.25E+01,
                0.12E+02,
                0.16E+02,
                0.25E+02,
                0.45E+02
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
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                a = a_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void gamma_inc_values(ref int n_data, ref double a, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GAMMA_INC_VALUES returns some values of the incomplete Gamma function.
            //
            //  Discussion:
            //
            //    The incomplete Gamma function is defined as:
            //
            //      Integral ( X <= T < oo ) T^(A-1) * exp(-T) dT.
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      Gamma[A,X]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 April 2010
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
            //    Output, ref double A, the parameter of the function.
            //
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 20;

            double[] a_vec =
            {
                0.10E+00,
                0.10E+00,
                0.10E+00,
                0.50E+00,
                0.50E+00,
                0.50E+00,
                0.10E+01,
                0.10E+01,
                0.10E+01,
                0.11E+01,
                0.11E+01,
                0.11E+01,
                0.20E+01,
                0.20E+01,
                0.20E+01,
                0.60E+01,
                0.60E+01,
                0.11E+02,
                0.26E+02,
                0.41E+02
            };

            double[] fx_vec =
            {
                2.490302836300570E+00,
                0.8718369702247978E+00,
                0.1079213896175866E+00,
                1.238121685818417E+00,
                0.3911298052193973E+00,
                0.01444722098952533E+00,
                0.9048374180359596E+00,
                0.3678794411714423E+00,
                0.006737946999085467E+00,
                0.8827966752611692E+00,
                0.3908330082003269E+00,
                0.008051456628620993E+00,
                0.9898141728888165E+00,
                0.5578254003710746E+00,
                0.007295055724436130E+00,
                114.9574754165633E+00,
                2.440923530031405E+00,
                280854.6620274718E+00,
                8.576480283455533E+24,
                2.085031346403364E+47
            };

            double[] x_vec =
            {
                0.30E-01,
                0.30E+00,
                0.15E+01,
                0.75E-01,
                0.75E+00,
                0.35E+01,
                0.10E+00,
                0.10E+01,
                0.50E+01,
                0.10E+00,
                0.10E+01,
                0.50E+01,
                0.15E+00,
                0.15E+01,
                0.70E+01,
                0.25E+01,
                0.12E+02,
                0.16E+02,
                0.25E+02,
                0.45E+02
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
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                a = a_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }
        
        public static void gamma_pdf_values(ref int n_data, ref double beta, ref double alpha, ref double x,
                ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GAMMA_PDF_VALUES returns some values of a Gamma PDF.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    29 July 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref double BETA, the rate parameter.
            //
            //    Output, ref double ALPHA, the shape parameter.
            //
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 10;

            double[] alpha_vec =
            {
                4.781587882544648,
                2.076535407379806,
                0.549783967662353,
                0.3086361453280091,
                3.773367432107051,
                4.487520304498656,
                0.06808445791730976,
                0.6155195788227712,
                4.562418534907164,
                4.114436583429598
            };

            double[] beta_vec =
            {
                1.092091484911879,
                2.808477213834471,
                1.287888961910225,
                3.169828561512062,
                2.006531407488083,
                0.009191855792026001,
                0.472723751058401,
                4.204237253278341,
                1.301514988836825,
                1.758143299519481
            };

            double[] fx_vec =
            {
                0.1672017697220646,
                0.8522122814089312,
                2.122272611165834,
                0.00006993771842317114,
                0.01679379733182281,
                0.0000000006687464259463117,
                0.001295436045931343,
                0.0,
                0.01189893036865762,
                0.3658836103539945
            };

            double[] x_vec =
            {
                4.942957250382744,
                0.2099361564793942,
                0.07173978623046406,
                2.587141553904492,
                4.743179115458788,
                1.974664495479389,
                5.126400502735112,
                -0.1534233427414219,
                0.5047170879434957,
                1.456220075613112
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                beta = 0.0;
                alpha = 0.0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                beta = beta_vec[n_data - 1];
                alpha = alpha_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void inverse_gamma_pdf_values(ref int n_data, ref double alpha, ref double beta,
                ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    INVERSE_GAMMA_PDF_VALUES returns values of the inverse Gamma PDF.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 August 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref double ALPHA, &BETA, the parameters.
            //
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 12;

            double[] alpha_vec =
            {
                1.00,
                1.00,
                1.00,
                1.00,
                1.00,
                1.00,
                1.00,
                1.00,
                2.00,
                3.00,
                4.00,
                5.00
            };

            double[] beta_vec =
            {
                0.50,
                0.50,
                0.50,
                0.50,
                2.00,
                3.00,
                4.00,
                5.00,
                2.00,
                2.00,
                2.00,
                2.00
            };

            double[] fx_vec =
            {
                0.3032653298563167,
                0.09735009788392561,
                0.04702676249392300,
                0.02757802820576861,
                0.1839397205857212,
                0.1673476201113224,
                0.1353352832366127,
                0.1026062482798735,
                0.07606179541223586,
                0.02535393180407862,
                0.005634207067573026,
                0.0009390345112621711
            };

            double[] x_vec =
            {
                1.00,
                2.00,
                3.00,
                4.00,
                2.00,
                2.00,
                2.00,
                2.00,
                3.00,
                3.00,
                3.00,
                3.00
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                alpha = 0.0;
                beta = 0.0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                alpha = alpha_vec[n_data - 1];
                beta = beta_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }


    }
}