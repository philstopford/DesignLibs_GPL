namespace Burkardt.TestValues
{
    public static class Truncated
    {
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
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref double MU, the mean of the distribution.
            //
            //    Output, ref double SIGMA, the standard deviation of the distribution.
            //
            //    Output, ref double A, &B, the lower and upper truncation limits.
            //
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
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

        public static void truncated_normal_ab_pdf_values(ref int n_data, ref double mu, ref double sigma,
                ref double a, ref double b, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRUNCATED_NORMAL_AB_PDF_VALUES: values of the Truncated Normal PDF.
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
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref double MU, the mean of the distribution.
            //
            //    Output, ref double SIGMA, the standard deviation of the distribution.
            //
            //    Output, ref double A, &B, the lower and upper truncation limits.
            //
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
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
                0.01543301171801836,
                0.01588394472270638,
                0.01624375997031919,
                0.01650575046469259,
                0.01666496869385951,
                0.01671838200940538,
                0.01666496869385951,
                0.01650575046469259,
                0.01624375997031919,
                0.01588394472270638,
                0.01543301171801836
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
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref double MU, the mean of the distribution.
            //
            //    Output, ref double SIGMA, the standard deviation of the distribution.
            //
            //    Output, ref double A, the lower truncation limit.
            //
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
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
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref double MU, the mean of the distribution.
            //
            //    Output, ref double SIGMA, the standard deviation of the distribution.
            //
            //    Output, ref double A, the lower truncation limit.
            //
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
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
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref double MU, the mean of the distribution.
            //
            //    Output, ref double SIGMA, the standard deviation of the distribution.
            //
            //    Output, ref double B, the upper truncation limit.
            //
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
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
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, ref double MU, the mean of the distribution.
            //
            //    Output, ref double SIGMA, the standard deviation of the distribution.
            //
            //    Output, ref double B, the upper truncation limit.
            //
            //    Output, ref double X, the argument of the function.
            //
            //    Output, ref double FX, the value of the function.
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