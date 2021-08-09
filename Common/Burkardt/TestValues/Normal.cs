namespace TestValues
{
    public static class Normal
    {
        public static void normal_cdf_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NORMAL_CDF_VALUES returns some values of the Normal CDF.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 February 2021
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Milton Abramowitz and Irene Stegun,
            //    Handbook of Mathematical Functions,
            //    US Department of Commerce, 1964.
            //
            //  Parameters:
            //
            //    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, double *X, the argument of the function.
            //
            //    Output double *FX, the value of the function.
            //
        {
            int N_MAX = 13;

            double[] fx_vec =
            {
                0.500000000000000E+00, 0.539827837277029E+00, 0.579259709439103E+00,
                0.617911422188953E+00, 0.655421741610324E+00, 0.691462461274013E+00,
                0.725746882249927E+00, 0.758036347776927E+00, 0.788144601416604E+00,
                0.815939874653241E+00, 0.841344746068543E+00, 0.933192798731142E+00,
                0.977249868051821E+00
            };
            double[] x_vec =
            {
                0.00E+00, 0.10E+00, 0.20E+00,
                0.30E+00, 0.40E+00, 0.50E+00,
                0.60E+00, 0.70E+00, 0.80E+00,
                0.90E+00, 1.00E+00, 1.50E+00,
                2.00E+00
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                x = 0.0E+00;
                fx = 0.0E+00;
            }
            else
            {
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void normal_01_cdf_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NORMAL_01_CDF_VALUES returns some values of the Normal 01 CDF.
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      Needs["Statistics`ContinuousDistributions`"]
            //      dist = NormalDistribution [ 0, 1 ]
            //      CDF [ dist, x ]
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
            int N_MAX = 17;

            double[] fx_vec =
            {
                0.5000000000000000E+00,
                0.5398278372770290E+00,
                0.5792597094391030E+00,
                0.6179114221889526E+00,
                0.6554217416103242E+00,
                0.6914624612740131E+00,
                0.7257468822499270E+00,
                0.7580363477769270E+00,
                0.7881446014166033E+00,
                0.8159398746532405E+00,
                0.8413447460685429E+00,
                0.9331927987311419E+00,
                0.9772498680518208E+00,
                0.9937903346742239E+00,
                0.9986501019683699E+00,
                0.9997673709209645E+00,
                0.9999683287581669E+00
            };

            double[] x_vec =
            {
                0.0000000000000000E+00,
                0.1000000000000000E+00,
                0.2000000000000000E+00,
                0.3000000000000000E+00,
                0.4000000000000000E+00,
                0.5000000000000000E+00,
                0.6000000000000000E+00,
                0.7000000000000000E+00,
                0.8000000000000000E+00,
                0.9000000000000000E+00,
                0.1000000000000000E+01,
                0.1500000000000000E+01,
                0.2000000000000000E+01,
                0.2500000000000000E+01,
                0.3000000000000000E+01,
                0.3500000000000000E+01,
                0.4000000000000000E+01
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

        public static void normal_01_pdf_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NORMAL_01_PDF_VALUES returns some values of the Normal 01 PDF.
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      Needs["Statistics`ContinuousDistributions`"]
            //      dist = NormalDistribution [ 0, 1 ]
            //      PDF [ dist, x ]
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
            int N_MAX = 10;

            double[] fx_vec =
            {
                0.03155059887555709,
                0.0005094586261557538,
                0.01235886992552887,
                0.353192862601275,
                0.3171212685764107,
                0.0009653372813755943,
                0.06083856556197816,
                0.003066504313116445,
                0.0005116437388114821,
                0.2246444116615346
            };

            double[] x_vec =
            {
                -2.252653624140994,
                3.650540612071437,
                2.636073871461605,
                0.4935635421351536,
                -0.6775433481923101,
                -3.471050120671749,
                -1.939377660943641,
                -3.120345651740235,
                -3.649368017767143,
                1.0717256984193
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

        public static void normal_cdf_values(ref int n_data, ref double mu, ref double sigma, ref double x,
                ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NORMAL_CDF_VALUES returns some values of the Normal CDF.
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      Needs["Statistics`ContinuousDistributions`"]
            //      dist = NormalDistribution [ mu, sigma ]
            //      CDF [ dist, x ]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 August 2004
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
                0.5000000000000000E+00,
                0.9772498680518208E+00,
                0.9999683287581669E+00,
                0.9999999990134124E+00,
                0.6914624612740131E+00,
                0.6305586598182364E+00,
                0.5987063256829237E+00,
                0.5792597094391030E+00,
                0.6914624612740131E+00,
                0.5000000000000000E+00,
                0.3085375387259869E+00,
                0.1586552539314571E+00
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

        public static void normal_pdf_values(ref int n_data, ref double mu, ref double sigma, ref double x,
                ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NORMAL_PDF_VALUES returns some values of the Normal PDF.
            //
            //  Discussion:
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      Needs["Statistics`ContinuousDistributions`"]
            //      dist = NormalDistribution [ mu, sigma ]
            //      PDF [ dist, x ]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    27 July 2015
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
            int N_MAX = 10;

            double[] fx_vec =
            {
                0.01180775937213258,
                0.006307849174478944,
                0.0147514774470322,
                0.9468437743011001,
                0.02140312299941794,
                0.05939959967353488,
                0.2348929157422787,
                0.007207515678571277,
                0.005944396897656727,
                0.03637663165771322
            };

            double[] mu_vec =
            {
                -56.31634060352484,
                12.33908855337884,
                -48.48444152359102,
                26.7931424604825,
                -19.73874370047668,
                -99.63232576831896,
                -81.09104995766396,
                68.16949013113364,
                -47.93940044652702,
                -29.67426801922078
            };

            double[] sigma_vec =
            {
                4.785956124893755,
                2.13500469923221,
                0.6387882883091059,
                0.4024634224214489,
                3.79790008346491,
                4.497769898408682,
                0.1667227687589636,
                0.7032091872463158,
                4.57117016420902,
                4.132147851761006
            };

            double[] x_vec =
            {
                -46.85424018542929,
                6.781057314200307,
                -50.23282168570062,
                26.67129012408019,
                -12.9643468135976,
                -103.6600156181528,
                -80.73183222587458,
                66.09155915000321,
                -58.53544475210675,
                -35.44773135435396
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

        public static void truncated_normal_b_pdf_values(ref int n_data, ref double mu, ref double sigma,
                ref double b, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRUNCATED_NORMAL_B_PDF_VALUES: values of the Upper Truncated Normal PDF.
            //
            //  Discussion:
            //
            //    The Normal distribution, with mean Mu and standard deviation Sigma,
            //    is truncated to the interval (-oo,B].
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 September 2013
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
            //    Output, double &MU, the mean of the distribution.
            //
            //    Output, double &SIGMA, the standard deviation of the distribution.
            //
            //    Output, double &B, the upper truncation limit.
            //
            //    Output, double &X, the argument of the function.
            //
            //    Output, double &FX, the value of the function.
            //
        {
            int N_MAX = 11;

            double[] b_vec
                    =
                    {
                        150.0,
                        150.0,
                        150.0,
                        150.0,
                        150.0,
                        150.0,
                        150.0,
                        150.0,
                        150.0,
                        150.0,
                        150.0
                    }
                ;

            double[] fx_vec
                    =
                    {
                        0.01507373507401876,
                        0.01551417047139894,
                        0.01586560931024694,
                        0.01612150073158793,
                        0.01627701240029317,
                        0.01632918226724295,
                        0.01627701240029317,
                        0.01612150073158793,
                        0.01586560931024694,
                        0.01551417047139894,
                        0.01507373507401876
                    }
                ;

            double[] mu_vec
                    =
                    {
                        100.0,
                        100.0,
                        100.0,
                        100.0,
                        100.0,
                        100.0,
                        100.0,
                        100.0,
                        100.0,
                        100.0,
                        100.0
                    }
                ;

            double[] sigma_vec
                    =
                    {
                        25.0,
                        25.0,
                        25.0,
                        25.0,
                        25.0,
                        25.0,
                        25.0,
                        25.0,
                        25.0,
                        25.0,
                        25.0
                    }
                ;

            double[] x_vec
                    =
                    {
                        90.0,
                        92.0,
                        94.0,
                        96.0,
                        98.0,
                        100.0,
                        102.0,
                        104.0,
                        106.0,
                        108.0,
                        110.0
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
                b = 0.0;
                mu = 0.0;
                sigma = 0.0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                b = b_vec[n_data - 1];
                mu = mu_vec[n_data - 1];
                sigma = sigma_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void truncated_normal_ab_cdf_values(ref int n_data, ref double mu, ref double sigma,
                ref double a, ref double b, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRUNCATED_NORMAL_AB_CDF_VALUES: values of the Truncated Normal CDF.
            //
            //  Discussion:
            //
            //    The Normal distribution, with mean Mu and standard deviation Sigma,
            //    is truncated to the interval [A,B].
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 September 2013
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
            //    Output, double &MU, the mean of the distribution.
            //
            //    Output, double &SIGMA, the standard deviation of the distribution.
            //
            //    Output, double &A, &B, the lower and upper truncation limits.
            //
            //    Output, double &X, the argument of the function.
            //
            //    Output, double &FX, the value of the function.
            //
        {
            int N_MAX = 11;

            double[] a_vec =
            {
                50.0,
                50.0,
                50.0,
                50.0,
                50.0,
                50.0,
                50.0,
                50.0,
                50.0,
                50.0,
                50.0
            };

            double[] b_vec =
            {
                150.0,
                150.0,
                150.0,
                150.0,
                150.0,
                150.0,
                150.0,
                150.0,
                150.0,
                150.0,
                150.0
            };

            double[] fx_vec =
            {
                0.3371694242213513,
                0.3685009225506048,
                0.4006444233448185,
                0.4334107066903040,
                0.4665988676496338,
                0.5000000000000000,
                0.5334011323503662,
                0.5665892933096960,
                0.5993555766551815,
                0.6314990774493952,
                0.6628305757786487
            };

            double[] mu_vec =
            {
                100.0,
                100.0,
                100.0,
                100.0,
                100.0,
                100.0,
                100.0,
                100.0,
                100.0,
                100.0,
                100.0
            };

            double[] sigma_vec =
            {
                25.0,
                25.0,
                25.0,
                25.0,
                25.0,
                25.0,
                25.0,
                25.0,
                25.0,
                25.0,
                25.0
            };

            double[] x_vec =
            {
                90.0,
                92.0,
                94.0,
                96.0,
                98.0,
                100.0,
                102.0,
                104.0,
                106.0,
                108.0,
                110.0
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
                mu = 0.0;
                sigma = 0.0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                a = a_vec[n_data - 1];
                b = b_vec[n_data - 1];
                mu = mu_vec[n_data - 1];
                sigma = sigma_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }

        }

        public static void truncated_normal_a_cdf_values(ref int n_data, ref double mu, ref double sigma,
                ref double a, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRUNCATED_NORMAL_A_CDF_VALUES: values of the Lower Truncated Normal CDF.
            //
            //  Discussion:
            //
            //    The Normal distribution, with mean Mu and standard deviation Sigma,
            //    is truncated to the interval [A,+oo).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 September 2013
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
            //    Output, double &MU, the mean of the distribution.
            //
            //    Output, double &SIGMA, the standard deviation of the distribution.
            //
            //    Output, double &A, the lower truncation limit.
            //
            //    Output, double &X, the argument of the function.
            //
            //    Output, double &FX, the value of the function.
            //
        {
            int N_MAX = 11;

            double[] a_vec =
            {
                50.0,
                50.0,
                50.0,
                50.0,
                50.0,
                50.0,
                50.0,
                50.0,
                50.0,
                50.0,
                50.0
            };

            double[] fx_vec =
            {
                0.3293202045481688,
                0.3599223134505957,
                0.3913175216041539,
                0.4233210140873113,
                0.4557365629792204,
                0.4883601253415709,
                0.5209836877039214,
                0.5533992365958304,
                0.5854027290789878,
                0.6167979372325460,
                0.6474000461349729
            };

            double[] mu_vec =
            {
                100.0,
                100.0,
                100.0,
                100.0,
                100.0,
                100.0,
                100.0,
                100.0,
                100.0,
                100.0,
                100.0
            };

            double[] sigma_vec =
            {
                25.0,
                25.0,
                25.0,
                25.0,
                25.0,
                25.0,
                25.0,
                25.0,
                25.0,
                25.0,
                25.0
            };

            double[] x_vec =
            {
                90.0,
                92.0,
                94.0,
                96.0,
                98.0,
                100.0,
                102.0,
                104.0,
                106.0,
                108.0,
                110.0
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
                mu = 0.0;
                sigma = 0.0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                a = a_vec[n_data - 1];
                mu = mu_vec[n_data - 1];
                sigma = sigma_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void truncated_normal_a_pdf_values(ref int n_data, ref double mu, ref double sigma,
                ref double a, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRUNCATED_NORMAL_A_PDF_VALUES: values of the Lower Truncated Normal PDF.
            //
            //  Discussion:
            //
            //    The Normal distribution, with mean Mu and standard deviation Sigma,
            //    is truncated to the interval [A,+oo).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 September 2013
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
            //    Output, double &MU, the mean of the distribution.
            //
            //    Output, double &SIGMA, the standard deviation of the distribution.
            //
            //    Output, double &A, the lower truncation limit.
            //
            //    Output, double &X, the argument of the function.
            //
            //    Output, double &FX, the value of the function.
            //
        {
            int N_MAX = 11;

            double[] a_vec =
            {
                50.0,
                50.0,
                50.0,
                50.0,
                50.0,
                50.0,
                50.0,
                50.0,
                50.0,
                50.0,
                50.0
            };

            double[] fx_vec =
            {
                0.01507373507401876,
                0.01551417047139894,
                0.01586560931024694,
                0.01612150073158793,
                0.01627701240029317,
                0.01632918226724295,
                0.01627701240029317,
                0.01612150073158793,
                0.01586560931024694,
                0.01551417047139894,
                0.01507373507401876
            };

            double[] mu_vec =
            {
                100.0,
                100.0,
                100.0,
                100.0,
                100.0,
                100.0,
                100.0,
                100.0,
                100.0,
                100.0,
                100.0
            };

            double[] sigma_vec =
            {
                25.0,
                25.0,
                25.0,
                25.0,
                25.0,
                25.0,
                25.0,
                25.0,
                25.0,
                25.0,
                25.0
            };

            double[] x_vec =
            {
                90.0,
                92.0,
                94.0,
                96.0,
                98.0,
                100.0,
                102.0,
                104.0,
                106.0,
                108.0,
                110.0
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
                mu = 0.0;
                sigma = 0.0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                a = a_vec[n_data - 1];
                mu = mu_vec[n_data - 1];
                sigma = sigma_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

        public static void truncated_normal_b_cdf_values(ref int n_data, ref double mu, ref double sigma,
                ref double b, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRUNCATED_NORMAL_B_CDF_VALUES: values of the upper Truncated Normal CDF.
            //
            //  Discussion:
            //
            //    The Normal distribution, with mean Mu and standard deviation Sigma,
            //    is truncated to the interval (-oo,B].
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 September 2013
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
            //    Output, double &MU, the mean of the distribution.
            //
            //    Output, double &SIGMA, the standard deviation of the distribution.
            //
            //    Output, double &B, the upper truncation limit.
            //
            //    Output, double &X, the argument of the function.
            //
            //    Output, double &FX, the value of the function.
            //
        {
            int N_MAX = 11;

            double[] b_vec =
            {
                150.0,
                150.0,
                150.0,
                150.0,
                150.0,
                150.0,
                150.0,
                150.0,
                150.0,
                150.0,
                150.0
            };

            double[] fx_vec =
            {
                0.3525999538650271,
                0.3832020627674540,
                0.4145972709210122,
                0.4466007634041696,
                0.4790163122960786,
                0.5116398746584291,
                0.5442634370207796,
                0.5766789859126887,
                0.6086824783958461,
                0.6400776865494043,
                0.6706797954518312
            };

            double[] mu_vec =
            {
                100.0,
                100.0,
                100.0,
                100.0,
                100.0,
                100.0,
                100.0,
                100.0,
                100.0,
                100.0,
                100.0
            };

            double[] sigma_vec =
            {
                25.0,
                25.0,
                25.0,
                25.0,
                25.0,
                25.0,
                25.0,
                25.0,
                25.0,
                25.0,
                25.0
            };

            double[] x_vec =
            {
                90.0,
                92.0,
                94.0,
                96.0,
                98.0,
                100.0,
                102.0,
                104.0,
                106.0,
                108.0,
                110.0
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                b = 0.0;
                mu = 0.0;
                sigma = 0.0;
                x = 0.0;
                fx = 0.0;
            }
            else
            {
                b = b_vec[n_data - 1];
                mu = mu_vec[n_data - 1];
                sigma = sigma_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }

    }
}