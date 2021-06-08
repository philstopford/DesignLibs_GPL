namespace Burkardt.CDFLib
{
    public static partial class CDF
    {
        public static void chi_noncentral_cdf_values(ref int n_data, ref double x, ref double lambda,
                ref int df, ref double cdf)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CHI_NONCENTRAL_CDF_VALUES returns values of the noncentral chi CDF.
            //
            //  Discussion:
            //
            //    The CDF of the noncentral chi square distribution can be evaluated
            //    within Mathematica by commands such as:
            //
            //      Needs["Statistics`ContinuousDistributions`"]
            //      CDF [ NoncentralChiSquareDistribution [ DF, LAMBDA ], X ]
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
            //    Stephen Wolfram,
            //    The Mathematica Book,
            //    Fourth Edition,
            //    Wolfram Media / Cambridge University Press, 1999.
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
            //    Output, double *LAMBDA, the noncentrality parameter.
            //
            //    Output, int *DF, the number of degrees of freedom.
            //
            //    Output, double *CDF, the noncentral chi CDF.
            //
        {
            int N_MAX = 27;

            double[] cdf_vec =  {
                0.839944E+00, 0.695906E+00, 0.535088E+00,
                0.764784E+00, 0.620644E+00, 0.469167E+00,
                0.307088E+00, 0.220382E+00, 0.150025E+00,
                0.307116E-02, 0.176398E-02, 0.981679E-03,
                0.165175E-01, 0.202342E-03, 0.498448E-06,
                0.151325E-01, 0.209041E-02, 0.246502E-03,
                0.263684E-01, 0.185798E-01, 0.130574E-01,
                0.583804E-01, 0.424978E-01, 0.308214E-01,
                0.105788E+00, 0.794084E-01, 0.593201E-01
            }
            ;
            int[] df_vec =  {
                1, 2, 3,
                1, 2, 3,
                1, 2, 3,
                1, 2, 3,
                60, 80, 100,
                1, 2, 3,
                10, 10, 10,
                10, 10, 10,
                10, 10, 10
            }
            ;
            double[] lambda_vec =  {
                0.5E+00, 0.5E+00, 0.5E+00,
                1.0E+00, 1.0E+00, 1.0E+00,
                5.0E+00, 5.0E+00, 5.0E+00,
                20.0E+00, 20.0E+00, 20.0E+00,
                30.0E+00, 30.0E+00, 30.0E+00,
                5.0E+00, 5.0E+00, 5.0E+00,
                2.0E+00, 3.0E+00, 4.0E+00,
                2.0E+00, 3.0E+00, 4.0E+00,
                2.0E+00, 3.0E+00, 4.0E+00
            }
            ;
            double[] x_vec =  {
                3.000E+00, 3.000E+00, 3.000E+00,
                3.000E+00, 3.000E+00, 3.000E+00,
                3.000E+00, 3.000E+00, 3.000E+00,
                3.000E+00, 3.000E+00, 3.000E+00,
                60.000E+00, 60.000E+00, 60.000E+00,
                0.050E+00, 0.050E+00, 0.050E+00,
                4.000E+00, 4.000E+00, 4.000E+00,
                5.000E+00, 5.000E+00, 5.000E+00,
                6.000E+00, 6.000E+00, 6.000E+00
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
                x = 0.0E+00;
                lambda = 0.0E+00;
                df = 0;
                cdf = 0.0E+00;
            }
            else
            {
                x = x_vec[n_data - 1];
                lambda = lambda_vec[n_data - 1];
                df = df_vec[n_data - 1];
                cdf = cdf_vec[n_data - 1];
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
            //    The value of CHI_CDF ( DF, X ) can be evaluated in Mathematica by
            //    commands like:
            //
            //      Needs["Statistics`ContinuousDistributions`"]
            //      CDF[ChiSquareDistribution[DF], X ]
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
            //    Stephen Wolfram,
            //    The Mathematica Book,
            //    Fourth Edition,
            //    Wolfram Media / Cambridge University Press, 1999.
            //
            //  Parameters:
            //
            //    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, int *A, the parameter of the function.
            //
            //    Output, double *X, the argument of the function.
            //
            //    Output, double *FX, the value of the function.
            //
        {
            int N_MAX = 21;

            int[] a_vec =  {
                1, 2, 1, 2,
                1, 2, 3, 4,
                1, 2, 3, 4,
                5, 3, 3, 3,
                3, 3, 10, 10,
                10
            }
            ;
            double[] fx_vec =  {
                0.0796557E+00, 0.00498752E+00, 0.112463E+00, 0.00995017E+00,
                0.472911E+00, 0.181269E+00, 0.0597575E+00, 0.0175231E+00,
                0.682689E+00, 0.393469E+00, 0.198748E+00, 0.090204E+00,
                0.0374342E+00, 0.427593E+00, 0.608375E+00, 0.738536E+00,
                0.828203E+00, 0.88839E+00, 0.000172116E+00, 0.00365985E+00,
                0.0185759E+00
            }
            ;
            double[] x_vec =  {
                0.01E+00, 0.01E+00, 0.02E+00, 0.02E+00,
                0.40E+00, 0.40E+00, 0.40E+00, 0.40E+00,
                1.00E+00, 1.00E+00, 1.00E+00, 1.00E+00,
                1.00E+00, 2.00E+00, 3.00E+00, 4.00E+00,
                5.00E+00, 6.00E+00, 1.00E+00, 2.00E+00,
                3.00E+00
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
                a = 0;
                x = 0.0E+00;
                fx = 0.0E+00;
            }
            else
            {
                a = a_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }
        }
    }
}