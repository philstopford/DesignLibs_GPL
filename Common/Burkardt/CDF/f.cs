namespace Burkardt.CDFLib
{
    public static partial class CDF
    {
        public static void f_cdf_values(ref int n_data, ref int a, ref int b, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    F_CDF_VALUES returns some values of the F CDF test function.
            //
            //  Discussion:
            //
            //    The value of F_CDF ( DFN, DFD, X ) can be evaluated in Mathematica by
            //    commands like:
            //
            //      Needs["Statistics`ContinuousDistributions`"]
            //      CDF[FRatioDistribution[ DFN, DFD ], X ]
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
            //    Output, int *A, int *B, the parameters of the function.
            //
            //    Output, double *X, the argument of the function.
            //
            //    Output, double *FX, the value of the function.
            //
        {
            int N_MAX = 20;

            int[] a_vec =  {
                1, 1, 5, 1,
                2, 4, 1, 6,
                8, 1, 3, 6,
                1, 1, 1, 1,
                2, 3, 4, 5
            }
            ;
            int[] b_vec =  {
                1, 5, 1, 5,
                10, 20, 5, 6,
                16, 5, 10, 12,
                5, 5, 5, 5,
                5, 5, 5, 5
            }
            ;
            double[] fx_vec =  {
                0.500000E+00, 0.499971E+00, 0.499603E+00, 0.749699E+00,
                0.750466E+00, 0.751416E+00, 0.899987E+00, 0.899713E+00,
                0.900285E+00, 0.950025E+00, 0.950057E+00, 0.950193E+00,
                0.975013E+00, 0.990002E+00, 0.994998E+00, 0.999000E+00,
                0.568799E+00, 0.535145E+00, 0.514343E+00, 0.500000E+00
            }
            ;
            double[] x_vec =  {
                1.00E+00, 0.528E+00, 1.89E+00, 1.69E+00,
                1.60E+00, 1.47E+00, 4.06E+00, 3.05E+00,
                2.09E+00, 6.61E+00, 3.71E+00, 3.00E+00,
                10.01E+00, 16.26E+00, 22.78E+00, 47.18E+00,
                1.00E+00, 1.00E+00, 1.00E+00, 1.00E+00
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
                b = 0;
                x = 0.0E+00;
                fx = 0.0E+00;
            }
            else
            {
                a = a_vec[n_data - 1];
                b = b_vec[n_data - 1];
                x = x_vec[n_data - 1];
                fx = fx_vec[n_data - 1];
            }

        }

        public static void f_noncentral_cdf_values(ref int n_data, ref int a, ref int b, ref double lambda,
                ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    F_NONCENTRAL_CDF_VALUES returns some values of the F CDF test function.
            //
            //  Discussion:
            //
            //    The value of NONCENTRAL_F_CDF ( DFN, DFD, LAMDA, X ) can be evaluated
            //    in Mathematica by commands like:
            //
            //      Needs["Statistics`ContinuousDistributions`"]
            //      CDF[NoncentralFRatioDistribution[ DFN, DFD, LAMBDA ], X ]
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
            //    Output, int *A, int *B, double *LAMBDA, the
            //    parameters of the function.
            //
            //    Output, double *X, the argument of the function.
            //
            //    Output, double *FX, the value of the function.
            //
        {
            int N_MAX = 22;

            int[] a_vec =  {
                1, 1, 1, 1,
                1, 1, 1, 1,
                1, 1, 2, 2,
                3, 3, 4, 4,
                5, 5, 6, 6,
                8, 16
            }
            ;
            int[] b_vec =  {
                1, 5, 5, 5,
                5, 5, 5, 5,
                5, 5, 5, 10,
                5, 5, 5, 5,
                1, 5, 6, 12,
                16, 8
            }
            ;
            double[] fx_vec =  {
                0.500000E+00, 0.636783E+00, 0.584092E+00, 0.323443E+00,
                0.450119E+00, 0.607888E+00, 0.705928E+00, 0.772178E+00,
                0.819105E+00, 0.317035E+00, 0.432722E+00, 0.450270E+00,
                0.426188E+00, 0.337744E+00, 0.422911E+00, 0.692767E+00,
                0.363217E+00, 0.421005E+00, 0.426667E+00, 0.446402E+00,
                0.844589E+00, 0.816368E+00
            }
            ;
            double[] lambda_vec =  {
                0.00E+00, 0.000E+00, 0.25E+00, 1.00E+00,
                1.00E+00, 1.00E+00, 1.00E+00, 1.00E+00,
                1.00E+00, 2.00E+00, 1.00E+00, 1.00E+00,
                1.00E+00, 2.00E+00, 1.00E+00, 1.00E+00,
                0.00E+00, 1.00E+00, 1.00E+00, 1.00E+00,
                1.00E+00, 1.00E+00
            }
            ;
            double[] x_vec =  {
                1.00E+00, 1.00E+00, 1.00E+00, 0.50E+00,
                1.00E+00, 2.00E+00, 3.00E+00, 4.00E+00,
                5.00E+00, 1.00E+00, 1.00E+00, 1.00E+00,
                1.00E+00, 1.00E+00, 1.00E+00, 2.00E+00,
                1.00E+00, 1.00E+00, 1.00E+00, 1.00E+00,
                2.00E+00, 2.00E+00
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
                b = 0;
                lambda = 0.0E+00;
                x = 0.0E+00;
                fx = 0.0E+00;
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
    }
}