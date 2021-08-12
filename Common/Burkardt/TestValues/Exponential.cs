namespace Burkardt.TestValues
{
    public static class Exponential
    {
        public static void exp_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    EXP_VALUES returns some values of the exponential function.
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      Exp[x]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 March 2008
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
            int N_MAX = 24;

            double[] fx_vec =
            {
                0.000045399929762484851536E+00,
                0.0067379469990854670966E+00,
                0.36787944117144232160E+00,
                1.0000000000000000000E+00,
                1.0000000100000000500E+00,
                1.0001000050001666708E+00,
                1.0010005001667083417E+00,
                1.0100501670841680575E+00,
                1.1051709180756476248E+00,
                1.2214027581601698339E+00,
                1.3498588075760031040E+00,
                1.4918246976412703178E+00,
                1.6487212707001281468E+00,
                1.8221188003905089749E+00,
                2.0137527074704765216E+00,
                2.2255409284924676046E+00,
                2.4596031111569496638E+00,
                2.7182818284590452354E+00,
                7.3890560989306502272E+00,
                23.140692632779269006E+00,
                148.41315910257660342E+00,
                22026.465794806716517E+00,
                4.8516519540979027797E+08,
                2.3538526683701998541E+17
            };

            double[] x_vec =
            {
                -10.0E+00,
                -5.0E+00,
                -1.0E+00,
                0.0E+00,
                0.00000001E+00,
                0.0001E+00,
                0.001E+00,
                0.01E+00,
                0.1E+00,
                0.2E+00,
                0.3E+00,
                0.4E+00,
                0.5E+00,
                0.6E+00,
                0.7E+00,
                0.8E+00,
                0.9E+00,
                1.0E+00,
                2.0E+00,
                3.1415926535897932385E+00,
                5.0E+00,
                10.0E+00,
                20.0E+00,
                40.0E+00
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

        public static void exp3_int_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    EXP3_INT_VALUES returns some values of the EXP3 integral function.
            //
            //  Discussion:
            //
            //    The function is defined by:
            //
            //      EXP3_INT(x) = Integral ( 0 <= t <= x ) exp ( -t^3 ) dt
            //
            //    The data was reported by McLeod.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    28 August 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Allan McLeod,
            //    Algorithm 757:
            //    MISCFUN: A software package to compute uncommon special functions,
            //    ACM Transactions on Mathematical Software,
            //    Volume 22, Number 3, September 1996, pages 288-301.
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
            int N_MAX = 20;

            double[] fx_vec =
            {
                0.19531249963620212007E-02,
                0.78124990686775522671E-02,
                0.31249761583499728667E-01,
                0.12493899888803079984E+00,
                0.48491714311363971332E+00,
                0.80751118213967145286E+00,
                0.86889265412623270696E+00,
                0.88861722235357162648E+00,
                0.89286018500218176869E+00,
                0.89295351429387631138E+00,
                0.89297479112737843939E+00,
                0.89297880579798112220E+00,
                0.89297950317496621294E+00,
                0.89297951152951902903E+00,
                0.89297951156918122102E+00,
                0.89297951156924734716E+00,
                0.89297951156924917298E+00,
                0.89297951156924921121E+00,
                0.89297951156924921122E+00,
                0.89297951156924921122E+00
            };

            double[] x_vec =
            {
                0.0019531250E+00,
                0.0078125000E+00,
                0.0312500000E+00,
                0.1250000000E+00,
                0.5000000000E+00,
                1.0000000000E+00,
                1.2500000000E+00,
                1.5000000000E+00,
                1.8750000000E+00,
                2.0000000000E+00,
                2.1250000000E+00,
                2.2500000000E+00,
                2.5000000000E+00,
                2.7500000000E+00,
                3.0000000000E+00,
                3.1250000000E+00,
                3.2500000000E+00,
                3.5000000000E+00,
                3.7500000000E+00,
                4.0000000000E+00
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

        public static void exponential_01_pdf_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    EXPONENTIAL_01_PDF_VALUES returns some values of the standard exponential PDF.
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
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 10;

            double[] fx_vec =
            {
                0.4959398481993681,
                0.00856777959135697,
                0.01720937842266235,
                0.07507070056996956,
                0.1679332083261492,
                0.0,
                0.399845179478639,
                0.9005384971416223,
                0.0,
                0.05044803826563792
            };

            double[] x_vec =
            {
                0.7013006334030669,
                4.759746670799113,
                4.062300786629853,
                2.589324935217918,
                1.784188948117787,
                -0.1363469579618277,
                0.9166778581012469,
                0.1047623644285883,
                -0.2589405122149109,
                2.986811417663269
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

        public static void exponential_cdf_values(ref int n_data, ref double lambda, ref double x,
                ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    EXPONENTIAL_CDF_VALUES returns some values of the Exponential CDF.
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      Needs["Statistics`ContinuousDistributions`"]
            //      dist = ExponentialDistribution [ lambda ]
            //      CDF [ dist, x ]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    29 August 2004
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
            //    Output, ref double LAMBDA, the parameter of the distribution.
            //
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 9;

            double[] fx_vec =
            {
                0.3934693402873666E+00,
                0.6321205588285577E+00,
                0.7768698398515702E+00,
                0.8646647167633873E+00,
                0.8646647167633873E+00,
                0.9816843611112658E+00,
                0.9975212478233336E+00,
                0.9996645373720975E+00,
                0.9999546000702375E+00
            };

            double[] lambda_vec =
            {
                0.5000000000000000E+00,
                0.5000000000000000E+00,
                0.5000000000000000E+00,
                0.5000000000000000E+00,
                0.1000000000000000E+01,
                0.2000000000000000E+01,
                0.3000000000000000E+01,
                0.4000000000000000E+01,
                0.5000000000000000E+01
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
                0.2000000000000000E+01
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                lambda = 0.0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                lambda = lambda_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void exponential_pdf_values(ref int n_data, ref double lambda, ref double x,
                ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    EXPONENTIAL_PDF_VALUES returns some values of the Exponential PDF.
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
            //    Output, ref double LAMBDA, the parameter of the distribution.
            //
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 10;

            double[] fx_vec =
            {
                0.0001446999730194618,
                0.06289850821824726,
                0.3663607831924032,
                0.3542787877169571,
                1.472582451176006E-12,
                1.829637907028298E-06,
                0.01173398427218792,
                0.0,
                0.1034724689882351,
                1.95394780436833
            };

            double[] lambda_vec =
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

            double[] x_vec =
            {
                9.558807522740191,
                5.573123971945631,
                0.5677992226519164,
                1.010563614677953,
                6.303053694254367,
                4.440343499102481,
                7.522202212856243,
                -0.08143245130010748,
                3.442598613603521,
                0.03753060499296568
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                lambda = 0.0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                lambda = lambda_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

    }
}