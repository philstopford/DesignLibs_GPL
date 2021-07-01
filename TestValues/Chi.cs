namespace TestValues
{
    public static class Chi
    {
        public static void chi_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CHI_VALUES returns some values of the hyperbolic cosine integral function.
            //
            //  Discussion:
            //
            //    The hyperbolic cosine integral is defined by
            //
            //      CHI(X) = gamma + log ( x )
            //        + integral ( 0 <= T < X ) ( cosh ( T ) - 1 ) / T  dT
            //
            //    where gamma is Euler's constant.
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      CoshIntegral[x]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 June 2007
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
                -0.05277684495649362,
                0.1577508933739787,
                0.3455691756953907,
                0.5183999848333915,
                0.6813138871854339,
                0.8378669409802082,
                1.141841924170595,
                1.445494075789644,
                1.759505807660965,
                2.092577214062032,
                2.452666922646915,
                3.524425488354165,
                4.960392094765610,
                6.959191927647393,
                9.813547558823186,
                13.96581164859243
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

        public static void chi_square_cdf_values(ref int n_data, ref int a, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CHI_SQUARE_CDF_VALUES returns some values of the Chi-Square CDF.
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      Needs["Statistics`ContinuousDistributions`"]
            //      dist = ChiSquareDistribution [ df ]
            //      CDF [ dist, x ]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    09 August 2004
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
            //    Output, ref int A, the parameter of the function.
            //
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            //
            int N_MAX = 21;

            int[] a_vec =
            {
                1, 2, 1, 2,
                1, 2, 3, 4,
                1, 2, 3, 4,
                5, 3, 3, 3,
                3, 3, 10, 10,
                10
            };

            double[] fx_vec =
            {
                0.7965567455405796E-01,
                0.4987520807317687E-02,
                0.1124629160182849E+00,
                0.9950166250831946E-02,
                0.4729107431344619E+00,
                0.1812692469220181E+00,
                0.5975750516063926E-01,
                0.1752309630642177E-01,
                0.6826894921370859E+00,
                0.3934693402873666E+00,
                0.1987480430987992E+00,
                0.9020401043104986E-01,
                0.3743422675270363E-01,
                0.4275932955291202E+00,
                0.6083748237289110E+00,
                0.7385358700508894E+00,
                0.8282028557032669E+00,
                0.8883897749052874E+00,
                0.1721156299558408E-03,
                0.3659846827343712E-02,
                0.1857593622214067E-01
            };

            double[] x_vec =
            {
                0.01E+00,
                0.01E+00,
                0.02E+00,
                0.02E+00,
                0.40E+00,
                0.40E+00,
                0.40E+00,
                0.40E+00,
                1.00E+00,
                1.00E+00,
                1.00E+00,
                1.00E+00,
                1.00E+00,
                2.00E+00,
                3.00E+00,
                4.00E+00,
                5.00E+00,
                6.00E+00,
                1.00E+00,
                2.00E+00,
                3.00E+00
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                a = 0;
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

        public static void chi_square_pdf_values(ref int n_data, ref double df, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CHI_SQUARE_PDF_VALUES returns some values of the Chi-Square PDF.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    01 August 2015
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
            //    Output, ref double DF, the degrees of freedom.
            //
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            //
            int N_MAX = 21;

            double[] df_vec =
            {
                1.0, 2.0, 1.0, 2.0,
                1.0, 2.0, 3.0, 4.0,
                1.0, 2.0, 3.0, 4.0,
                5.0, 3.0, 3.0, 3.0,
                3.0, 3.0, 10.0, 10.0,
                10.0
            };

            double[] fx_vec =
            {
                3.969525474770117,
                0.4975062395963412,
                2.792879016972342,
                0.4950249168745841,
                0.5164415474672784,
                0.4093653765389909,
                0.2065766189869113,
                0.08187307530779819,
                0.2419707245191434,
                0.3032653298563167,
                0.2419707245191434,
                0.1516326649281584,
                0.08065690817304777,
                0.2075537487102974,
                0.1541803298037693,
                0.1079819330263761,
                0.07322491280963248,
                0.04865217332964145,
                0.0007897534631674914,
                0.00766415502440505,
                0.02353325907815472
            };

            double[] x_vec =
            {
                0.01E+00,
                0.01E+00,
                0.02E+00,
                0.02E+00,
                0.40E+00,
                0.40E+00,
                0.40E+00,
                0.40E+00,
                1.00E+00,
                1.00E+00,
                1.00E+00,
                1.00E+00,
                1.00E+00,
                2.00E+00,
                3.00E+00,
                4.00E+00,
                5.00E+00,
                6.00E+00,
                1.00E+00,
                2.00E+00,
                3.00E+00
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                df = 0.0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                df = df_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void chi_square_noncentral_cdf_values(ref int n_data, ref int df, ref double lambda,
                ref double x, ref double cdf)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CHI_SQUARE_NONCENTRAL_CDF_VALUES returns values of the noncentral chi CDF.
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      Needs["Statistics`ContinuousDistributions`"]
            //      dist = NoncentralChiSquareDistribution [ df, lambda ]
            //      CDF [ dist, x ]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    30 August 2004
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
            //    Output, ref int DF, the number of degrees of freedom.
            //
            //    Output, ref double LAMBDA, the noncentrality parameter.
            //
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double CDF, the noncentral chi CDF.
            //
        {
            int N_MAX = 28;

            double[] cdf_vec =
            {
                0.8399444269398261E+00,
                0.6959060300435139E+00,
                0.5350879697078847E+00,
                0.7647841496310313E+00,
                0.6206436532195436E+00,
                0.4691667375373180E+00,
                0.3070884345937569E+00,
                0.2203818092990903E+00,
                0.1500251895581519E+00,
                0.3071163194335791E-02,
                0.1763982670131894E-02,
                0.9816792594625022E-03,
                0.1651753140866208E-01,
                0.2023419573950451E-03,
                0.4984476352854074E-06,
                0.1513252400654827E-01,
                0.2090414910614367E-02,
                0.2465021206048452E-03,
                0.2636835050342939E-01,
                0.1857983220079215E-01,
                0.1305736595486640E-01,
                0.5838039534819351E-01,
                0.4249784402463712E-01,
                0.3082137716021596E-01,
                0.1057878223400849E+00,
                0.7940842984598509E-01,
                0.5932010895599639E-01,
                0.2110395656918684E+00
            };

            int[] df_vec =
            {
                1, 2, 3,
                1, 2, 3,
                1, 2, 3,
                1, 2, 3,
                60, 80, 100,
                1, 2, 3,
                10, 10, 10,
                10, 10, 10,
                10, 10, 10,
                8
            };

            double[] lambda_vec =
            {
                0.5E+00,
                0.5E+00,
                0.5E+00,
                1.0E+00,
                1.0E+00,
                1.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                20.0E+00,
                20.0E+00,
                20.0E+00,
                30.0E+00,
                30.0E+00,
                30.0E+00,
                5.0E+00,
                5.0E+00,
                5.0E+00,
                2.0E+00,
                3.0E+00,
                4.0E+00,
                2.0E+00,
                3.0E+00,
                4.0E+00,
                2.0E+00,
                3.0E+00,
                4.0E+00,
                0.5E+00
            };

            double[] x_vec =
            {
                3.000E+00,
                3.000E+00,
                3.000E+00,
                3.000E+00,
                3.000E+00,
                3.000E+00,
                3.000E+00,
                3.000E+00,
                3.000E+00,
                3.000E+00,
                3.000E+00,
                3.000E+00,
                60.000E+00,
                60.000E+00,
                60.000E+00,
                0.050E+00,
                0.050E+00,
                0.050E+00,
                4.000E+00,
                4.000E+00,
                4.000E+00,
                5.000E+00,
                5.000E+00,
                5.000E+00,
                6.000E+00,
                6.000E+00,
                6.000E+00,
                5.000E+00
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
                lambda = 0.0;
                df = 0;
                cdf = 0.0;
            }
            else
            {
                x = x_vec[n_data - 1];
                lambda = lambda_vec[n_data - 1];
                df = df_vec[n_data - 1];
                cdf = cdf_vec[n_data - 1];
            }
        }

        public static void inverse_chi_square_pdf_values(ref int n_data, ref double df, ref double x,
                ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    INVERSE_CHI_SQUARE_PDF_VALUES returns values of the inverse Chi-Square PDF.
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
            //    Output, ref double DF, the degrees of freedom.
            //
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 21;

            double[] df_vec =
            {
                1.0,
                2.0,
                1.0,
                2.0,
                1.0,
                2.0,
                3.0,
                4.0,
                1.0,
                2.0,
                3.0,
                4.0,
                5.0,
                3.0,
                3.0,
                3.0,
                3.0,
                3.0,
                10.0,
                10.0,
                10.0
            };

            double[] fx_vec =
            {
                0.08500366602520342,
                0.3368973499542734,
                0.3661245640481622,
                1.026062482798735,
                0.4518059816704532,
                0.8953274901880941,
                1.129514954176133,
                1.119159362735118,
                0.2419707245191433,
                0.3032653298563167,
                0.2419707245191433,
                0.1516326649281584,
                0.08065690817304778,
                0.05492391118346530,
                0.02166329508030457,
                0.01100204146138436,
                0.006457369034861447,
                0.004162370481945731,
                0.0007897534631674914,
                0.00001584474249412852,
                1.511920090468204E-06
            };

            double[] x_vec =
            {
                0.10,
                0.10,
                0.20,
                0.20,
                0.40,
                0.40,
                0.40,
                0.40,
                1.00,
                1.00,
                1.00,
                1.00,
                1.00,
                2.00,
                3.00,
                4.00,
                5.00,
                6.00,
                1.00,
                2.00,
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
                df = 0.0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                df = df_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void scaled_inverse_chi_square_pdf_values(ref int n_data, ref double df, ref double xi,
                ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SCALED_INVERSE_CHI_SQUARE_PDF_VALUES: scaled inverse Chi-Square PDF values.
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
            //    Output, ref double DF, the degrees of freedom.
            //
            //    Output, ref double XI, the scale parameter.
            //
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
            //
        {
            int N_MAX = 18;

            double[] df_vec =
            {
                1.0,
                2.0,
                1.0,
                2.0,
                1.0,
                2.0,
                1.0,
                2.0,
                1.0,
                2.0,
                1.0,
                2.0,
                1.0,
                2.0,
                1.0,
                2.0,
                1.0,
                2.0
            };

            double[] fx_vec =
            {
                0.7322491280963244,
                0.3368973499542734,
                0.9036119633409063,
                1.026062482798735,
                0.5968580144169457,
                0.8953274901880941,
                0.08500366602520342,
                0.004539992976248485,
                0.3661245640481622,
                0.1684486749771367,
                0.4518059816704532,
                0.5130312413993675,
                0.0008099910956089117,
                4.122307244877116E-07,
                0.04250183301260171,
                0.002269996488124243,
                0.1830622820240811,
                0.08422433748856834
            };

            double[] x_vec =
            {
                0.10,
                0.10,
                0.20,
                0.20,
                0.40,
                0.40,
                0.10,
                0.10,
                0.20,
                0.20,
                0.40,
                0.40,
                0.10,
                0.10,
                0.20,
                0.20,
                0.40
            };

            double[] xi_vec =
            {
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                0.50,
                1.00,
                1.00,
                1.00,
                1.00,
                1.00,
                1.00,
                2.00,
                2.00,
                2.00,
                2.00,
                2.00,
                2.00
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                df = 0.0;
                xi = 0.0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                df = df_vec[n_data - 1];
                xi = xi_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }


    }
}