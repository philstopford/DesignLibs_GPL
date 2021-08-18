namespace Burkardt.Values
{
    public static class Beta
    {
        public static void beta_cdf_values(ref int n_data, ref double a, ref double b, ref double x,
                ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BETA_CDF_VALUES returns some values of the Beta CDF.
            //
            //  Discussion:
            //
            //    The incomplete Beta function may be written
            //
            //      BETA_INC(A,B,X) = Integral (0 <= t <= X) T^(A-1) * (1-T)^(B-1) dT
            //                      / Integral (0 <= t <= 1) T^(A-1) * (1-T)^(B-1) dT
            //
            //    Thus,
            //
            //      BETA_INC(A,B,0.0) = 0.0;
            //      BETA_INC(A,B,1.0) = 1.0
            //
            //    The incomplete Beta function is also sometimes called the
            //    "modified" Beta function, or the "normalized" Beta function
            //    or the Beta CDF (cumulative density function.
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      BETA[X,A,B] / BETA[A,B]
            //
            //    The function can also be evaluated by using the Statistics package:
            //
            //      Needs["Statistics`ContinuousDistributions`"]
            //      dist = BetaDistribution [ a, b ]
            //      CDF [ dist, x ]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    28 April 2013
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
            //    Karl Pearson,
            //    Tables of the Incomplete Beta Function,
            //    Cambridge University Press, 1968.
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
            //    Output, ref double A, &B, the parameters of the function.
            //
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 45;

            double[] a_vec =
            {
                0.5E+00,
                0.5E+00,
                0.5E+00,
                1.0E+00,
                1.0E+00,
                1.0E+00,
                1.0E+00,
                1.0E+00,
                2.0E+00,
                2.0E+00,
                2.0E+00,
                2.0E+00,
                2.0E+00,
                2.0E+00,
                2.0E+00,
                2.0E+00,
                2.0E+00,
                5.5E+00,
                10.0E+00,
                10.0E+00,
                10.0E+00,
                10.0E+00,
                20.0E+00,
                20.0E+00,
                20.0E+00,
                20.0E+00,
                20.0E+00,
                30.0E+00,
                30.0E+00,
                40.0E+00,
                0.1E+01,
                0.1E+01,
                0.1E+01,
                0.1E+01,
                0.1E+01,
                0.1E+01,
                0.1E+01,
                0.1E+01,
                0.2E+01,
                0.3E+01,
                0.4E+01,
                0.5E+01,
                1.30625,
                1.30625,
                1.30625
            };

            double[] b_vec =
            {
                0.5E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                1.0E+00,
                2.0E+00,
                2.0E+00,
                2.0E+00,
                2.0E+00,
                2.0E+00,
                2.0E+00,
                2.0E+00,
                2.0E+00,
                2.0E+00,
                5.0E+00,
                0.5E+00,
                5.0E+00,
                5.0E+00,
                10.0E+00,
                5.0E+00,
                10.0E+00,
                10.0E+00,
                20.0E+00,
                20.0E+00,
                10.0E+00,
                10.0E+00,
                20.0E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                0.2E+01,
                0.3E+01,
                0.4E+01,
                0.5E+01,
                0.2E+01,
                0.2E+01,
                0.2E+01,
                0.2E+01,
                11.7562,
                11.7562,
                11.7562
            };

            double[] fx_vec =
            {
                0.6376856085851985E-01,
                0.2048327646991335E+00,
                0.1000000000000000E+01,
                0.0000000000000000E+00,
                0.5012562893380045E-02,
                0.5131670194948620E-01,
                0.2928932188134525E+00,
                0.5000000000000000E+00,
                0.2800000000000000E-01,
                0.1040000000000000E+00,
                0.2160000000000000E+00,
                0.3520000000000000E+00,
                0.5000000000000000E+00,
                0.6480000000000000E+00,
                0.7840000000000000E+00,
                0.8960000000000000E+00,
                0.9720000000000000E+00,
                0.4361908850559777E+00,
                0.1516409096347099E+00,
                0.8978271484375000E-01,
                0.1000000000000000E+01,
                0.5000000000000000E+00,
                0.4598773297575791E+00,
                0.2146816102371739E+00,
                0.9507364826957875E+00,
                0.5000000000000000E+00,
                0.8979413687105918E+00,
                0.2241297491808366E+00,
                0.7586405487192086E+00,
                0.7001783247477069E+00,
                0.5131670194948620E-01,
                0.1055728090000841E+00,
                0.1633399734659245E+00,
                0.2254033307585166E+00,
                0.3600000000000000E+00,
                0.4880000000000000E+00,
                0.5904000000000000E+00,
                0.6723200000000000E+00,
                0.2160000000000000E+00,
                0.8370000000000000E-01,
                0.3078000000000000E-01,
                0.1093500000000000E-01,
                0.918884684620518,
                0.21052977489419,
                0.1824130512500673
            };

            double[] x_vec =
            {
                0.01E+00,
                0.10E+00,
                1.00E+00,
                0.00E+00,
                0.01E+00,
                0.10E+00,
                0.50E+00,
                0.50E+00,
                0.10E+00,
                0.20E+00,
                0.30E+00,
                0.40E+00,
                0.50E+00,
                0.60E+00,
                0.70E+00,
                0.80E+00,
                0.90E+00,
                0.50E+00,
                0.90E+00,
                0.50E+00,
                1.00E+00,
                0.50E+00,
                0.80E+00,
                0.60E+00,
                0.80E+00,
                0.50E+00,
                0.60E+00,
                0.70E+00,
                0.80E+00,
                0.70E+00,
                0.10E+00,
                0.20E+00,
                0.30E+00,
                0.40E+00,
                0.20E+00,
                0.20E+00,
                0.20E+00,
                0.20E+00,
                0.30E+00,
                0.30E+00,
                0.30E+00,
                0.30E+00,
                0.225609,
                0.0335568,
                0.0295222
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
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                a = a_vec[n_data - 1];
                b = b_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void beta_inc_values(ref int n_data, ref double a, ref double b, ref double x,
                ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BETA_INC_VALUES returns some values of the incomplete Beta function.
            //
            //  Discussion:
            //
            //    The incomplete Beta function may be written
            //
            //      BETA_INC(A,B,X) = Integral (0 <= t <= X) T^(A-1) * (1-T)^(B-1) dT
            //                      / Integral (0 <= t <= 1) T^(A-1) * (1-T)^(B-1) dT
            //
            //    Thus,
            //
            //      BETA_INC(A,B,0.0) = 0.0;
            //      BETA_INC(A,B,1.0) = 1.0
            //
            //    The incomplete Beta function is also sometimes called the
            //    "modified" Beta function, or the "normalized" Beta function
            //    or the Beta CDF (cumulative density function.
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      BETA[X,A,B] / BETA[A,B]
            //
            //    The function can also be evaluated by using the Statistics package:
            //
            //      Needs["Statistics`ContinuousDistributions`"]
            //      dist = BetaDistribution [ a, b ]
            //      CDF [ dist, x ]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    28 April 2013
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
            //    Karl Pearson,
            //    Tables of the Incomplete Beta Function,
            //    Cambridge University Press, 1968.
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
            //    Output, ref double A, &B, the parameters of the function.
            //
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 45;

            double[] a_vec =
            {
                0.5E+00,
                0.5E+00,
                0.5E+00,
                1.0E+00,
                1.0E+00,
                1.0E+00,
                1.0E+00,
                1.0E+00,
                2.0E+00,
                2.0E+00,
                2.0E+00,
                2.0E+00,
                2.0E+00,
                2.0E+00,
                2.0E+00,
                2.0E+00,
                2.0E+00,
                5.5E+00,
                10.0E+00,
                10.0E+00,
                10.0E+00,
                10.0E+00,
                20.0E+00,
                20.0E+00,
                20.0E+00,
                20.0E+00,
                20.0E+00,
                30.0E+00,
                30.0E+00,
                40.0E+00,
                0.1E+01,
                0.1E+01,
                0.1E+01,
                0.1E+01,
                0.1E+01,
                0.1E+01,
                0.1E+01,
                0.1E+01,
                0.2E+01,
                0.3E+01,
                0.4E+01,
                0.5E+01,
                1.30625,
                1.30625,
                1.30625
            };

            double[] b_vec =
            {
                0.5E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                1.0E+00,
                2.0E+00,
                2.0E+00,
                2.0E+00,
                2.0E+00,
                2.0E+00,
                2.0E+00,
                2.0E+00,
                2.0E+00,
                2.0E+00,
                5.0E+00,
                0.5E+00,
                5.0E+00,
                5.0E+00,
                10.0E+00,
                5.0E+00,
                10.0E+00,
                10.0E+00,
                20.0E+00,
                20.0E+00,
                10.0E+00,
                10.0E+00,
                20.0E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                0.5E+00,
                0.2E+01,
                0.3E+01,
                0.4E+01,
                0.5E+01,
                0.2E+01,
                0.2E+01,
                0.2E+01,
                0.2E+01,
                11.7562,
                11.7562,
                11.7562
            };

            double[] fx_vec =
            {
                0.6376856085851985E-01,
                0.2048327646991335E+00,
                0.1000000000000000E+01,
                0.0000000000000000E+00,
                0.5012562893380045E-02,
                0.5131670194948620E-01,
                0.2928932188134525E+00,
                0.5000000000000000E+00,
                0.2800000000000000E-01,
                0.1040000000000000E+00,
                0.2160000000000000E+00,
                0.3520000000000000E+00,
                0.5000000000000000E+00,
                0.6480000000000000E+00,
                0.7840000000000000E+00,
                0.8960000000000000E+00,
                0.9720000000000000E+00,
                0.4361908850559777E+00,
                0.1516409096347099E+00,
                0.8978271484375000E-01,
                0.1000000000000000E+01,
                0.5000000000000000E+00,
                0.4598773297575791E+00,
                0.2146816102371739E+00,
                0.9507364826957875E+00,
                0.5000000000000000E+00,
                0.8979413687105918E+00,
                0.2241297491808366E+00,
                0.7586405487192086E+00,
                0.7001783247477069E+00,
                0.5131670194948620E-01,
                0.1055728090000841E+00,
                0.1633399734659245E+00,
                0.2254033307585166E+00,
                0.3600000000000000E+00,
                0.4880000000000000E+00,
                0.5904000000000000E+00,
                0.6723200000000000E+00,
                0.2160000000000000E+00,
                0.8370000000000000E-01,
                0.3078000000000000E-01,
                0.1093500000000000E-01,
                0.918884684620518,
                0.21052977489419,
                0.1824130512500673
            };

            double[] x_vec =
            {
                0.01E+00,
                0.10E+00,
                1.00E+00,
                0.00E+00,
                0.01E+00,
                0.10E+00,
                0.50E+00,
                0.50E+00,
                0.10E+00,
                0.20E+00,
                0.30E+00,
                0.40E+00,
                0.50E+00,
                0.60E+00,
                0.70E+00,
                0.80E+00,
                0.90E+00,
                0.50E+00,
                0.90E+00,
                0.50E+00,
                1.00E+00,
                0.50E+00,
                0.80E+00,
                0.60E+00,
                0.80E+00,
                0.50E+00,
                0.60E+00,
                0.70E+00,
                0.80E+00,
                0.70E+00,
                0.10E+00,
                0.20E+00,
                0.30E+00,
                0.40E+00,
                0.20E+00,
                0.20E+00,
                0.20E+00,
                0.20E+00,
                0.30E+00,
                0.30E+00,
                0.30E+00,
                0.30E+00,
                0.225609,
                0.0335568,
                0.0295222
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
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                a = a_vec[n_data - 1];
                b = b_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void beta_log_values(ref int n_data, ref double x, ref double y, ref double fxy)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BETA_LOG_VALUES returns some values of the logarithm of the Beta function.
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      Log[Beta[x]]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 August 2004
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
            //    Output, ref double X, &Y, the arguments of the function.
            //
            //    Output, ref double FXY, the value of the function.
            //
        {
            int N_MAX = 17;

            double[] fxy_vec =
            {
                0.1609437912434100E+01,
                0.9162907318741551E+00,
                0.5108256237659907E+00,
                0.2231435513142098E+00,
                0.1609437912434100E+01,
                0.9162907318741551E+00,
                0.0000000000000000E+00,
                -0.1791759469228055E+01,
                -0.3401197381662155E+01,
                -0.4941642422609304E+01,
                -0.6445719819385578E+01,
                -0.3737669618283368E+01,
                -0.5123963979403259E+01,
                -0.6222576268071369E+01,
                -0.7138866999945524E+01,
                -0.7927324360309794E+01,
                -0.9393661429103221E+01
            };

            double[] x_vec =
            {
                0.2E+00,
                0.4E+00,
                0.6E+00,
                0.8E+00,
                1.0E+00,
                1.0E+00,
                1.0E+00,
                2.0E+00,
                3.0E+00,
                4.0E+00,
                5.0E+00,
                6.0E+00,
                6.0E+00,
                6.0E+00,
                6.0E+00,
                6.0E+00,
                7.0E+00
            };

            double[] y_vec =
            {
                1.0E+00,
                1.0E+00,
                1.0E+00,
                1.0E+00,
                0.2E+00,
                0.4E+00,
                1.0E+00,
                2.0E+00,
                3.0E+00,
                4.0E+00,
                5.0E+00,
                2.0E+00,
                3.0E+00,
                4.0E+00,
                5.0E+00,
                6.0E+00,
                7.0E+00
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
                y = 0.0;
                fxy = 0.0;
            }
            else
            {
                x = x_vec[n_data - 1];
                y = y_vec[n_data - 1];
                fxy = fxy_vec[n_data - 1];
            }
        }

        public static void beta_noncentral_cdf_values(ref int n_data, ref double a, ref double b,
                ref double lambda, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BETA_NONCENTRAL_CDF_VALUES returns some values of the noncentral Beta CDF.
            //
            //  Discussion:
            //
            //    The values presented here are taken from the reference, where they
            //    were given to a limited number of decimal places.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    24 January 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    R Chattamvelli, R Shanmugam,
            //    Algorithm AS 310:
            //    Computing the Non-central Beta Distribution Function,
            //    Applied Statistics,
            //    Volume 46, Number 1, 1997, pages 146-156.
            //
            //  Parameters:
            //
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0
            //    before the first call.  On each call, the routine increments N_DATA by 1,
            //    and returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref double A, &B, the shape parameters.
            //
            //    Output, ref double LAMBDA, the noncentrality parameter.
            //
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 25;

            double[] a_vec =
            {
                5.0,
                5.0,
                5.0,
                10.0,
                10.0,
                10.0,
                20.0,
                20.0,
                20.0,
                10.0,
                10.0,
                15.0,
                20.0,
                20.0,
                20.0,
                30.0,
                30.0,
                10.0,
                10.0,
                10.0,
                15.0,
                10.0,
                12.0,
                30.0,
                35.0
            };
            double[] b_vec =
            {
                5.0,
                5.0,
                5.0,
                10.0,
                10.0,
                10.0,
                20.0,
                20.0,
                20.0,
                20.0,
                10.0,
                5.0,
                10.0,
                30.0,
                50.0,
                20.0,
                40.0,
                5.0,
                10.0,
                30.0,
                20.0,
                5.0,
                17.0,
                30.0,
                30.0
            };
            double[] fx_vec =
            {
                0.4563021,
                0.1041337,
                0.6022353,
                0.9187770,
                0.6008106,
                0.0902850,
                0.9998655,
                0.9925997,
                0.9641112,
                0.9376626573,
                0.7306817858,
                0.1604256918,
                0.1867485313,
                0.6559386874,
                0.9796881486,
                0.1162386423,
                0.9930430054,
                0.0506899273,
                0.1030959706,
                0.9978417832,
                0.2555552369,
                0.0668307064,
                0.0113601067,
                0.7813366615,
                0.8867126477
            };
            double[] lambda_vec =
            {
                54.0,
                140.0,
                170.0,
                54.0,
                140.0,
                250.0,
                54.0,
                140.0,
                250.0,
                150.0,
                120.0,
                80.0,
                110.0,
                65.0,
                130.0,
                80.0,
                130.0,
                20.0,
                54.0,
                80.0,
                120.0,
                55.0,
                64.0,
                140.0,
                20.0
            };
            double[] x_vec =
            {
                0.8640,
                0.9000,
                0.9560,
                0.8686,
                0.9000,
                0.9000,
                0.8787,
                0.9000,
                0.9220,
                0.868,
                0.900,
                0.880,
                0.850,
                0.660,
                0.720,
                0.720,
                0.800,
                0.644,
                0.700,
                0.780,
                0.760,
                0.795,
                0.560,
                0.800,
                0.670
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
                lambda = 0.0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                a = a_vec[n_data - 1];
                b = b_vec[n_data - 1];
                lambda = lambda_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void beta_pdf_values(ref int n_data, ref double alpha, ref double beta, ref double x,
                ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BETA_PDF_VALUES returns some values of the Beta PDF.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    30 July 2015
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
            //    Output, ref double ALPHA, &BETA, the parameters of the function.
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

            double[] beta_vec =
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

            double[] fx_vec =
            {
                0.002826137156803199,
                0.04208950342768649,
                0.2184064957817208,
                0.1335142301445414,
                0.1070571849830009,
                0.005796394377470491,
                0.5518796772414584,
                0.0,
                2.87907465409348,
                2.126992854611924
            };

            double[] x_vec =
            {
                0.8667224264776531,
                0.04607764003473368,
                0.02211617261254013,
                0.4582543823302144,
                0.8320834756642252,
                0.3520587633290876,
                0.898529119425846,
                -0.01692420862048847,
                0.09718884992568674,
                0.2621671905296927
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

        public static void beta_values(ref int n_data, ref double x, ref double y, ref double fxy)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BETA_VALUES returns some values of the Beta function.
            //
            //  Discussion:
            //
            //    Beta(X,Y) = ( Gamma(X) * Gamma(Y) ) / Gamma(X+Y)
            //
            //    Both X and Y must be greater than 0.
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      Beta[X,Y]
            //
            //  Properties:
            //
            //    Beta(X,Y) = Beta(Y,X).
            //    Beta(X,Y) = Integral ( 0 <= T <= 1 ) T^(X-1) (1-T)^(Y-1) dT.
            //    Beta(X,Y) = Gamma(X) * Gamma(Y) / Gamma(X+Y)
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 August 2004
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
            //    Output, ref double X, &Y, the arguments of the function.
            //
            //    Output, ref double FXY, the value of the function.
            //
        {
            int N_MAX = 17;

            double[] b_vec =
            {
                0.5000000000000000E+01,
                0.2500000000000000E+01,
                0.1666666666666667E+01,
                0.1250000000000000E+01,
                0.5000000000000000E+01,
                0.2500000000000000E+01,
                0.1000000000000000E+01,
                0.1666666666666667E+00,
                0.3333333333333333E-01,
                0.7142857142857143E-02,
                0.1587301587301587E-02,
                0.2380952380952381E-01,
                0.5952380952380952E-02,
                0.1984126984126984E-02,
                0.7936507936507937E-03,
                0.3607503607503608E-03,
                0.8325008325008325E-04
            };

            double[] x_vec =
            {
                0.2E+00,
                0.4E+00,
                0.6E+00,
                0.8E+00,
                1.0E+00,
                1.0E+00,
                1.0E+00,
                2.0E+00,
                3.0E+00,
                4.0E+00,
                5.0E+00,
                6.0E+00,
                6.0E+00,
                6.0E+00,
                6.0E+00,
                6.0E+00,
                7.0E+00
            };

            double[] y_vec =
            {
                1.0E+00,
                1.0E+00,
                1.0E+00,
                1.0E+00,
                0.2E+00,
                0.4E+00,
                1.0E+00,
                2.0E+00,
                3.0E+00,
                4.0E+00,
                5.0E+00,
                2.0E+00,
                3.0E+00,
                4.0E+00,
                5.0E+00,
                6.0E+00,
                7.0E+00
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
                y = 0.0;
                fxy = 0.0;
            }
            else
            {
                x = x_vec[n_data - 1];
                y = y_vec[n_data - 1];
                fxy = b_vec[n_data - 1];
            }
        }

    }
}