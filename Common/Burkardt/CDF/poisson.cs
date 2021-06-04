namespace Burkardt.CDFLib
{
    public static partial class CDF
    {
        public static void poisson_cdf_values(ref int n_data, ref double a, ref int x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POISSON_CDF_VALUES returns some values of the Poisson CDF.
            //
            //  Discussion:
            //
            //    CDF(X)(A) is the probability of at most X successes in unit time,
            //    given that the expected mean number of successes is A.
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
            //    Daniel Zwillinger,
            //    CRC Standard Mathematical Tables and Formulae,
            //    30th Edition, CRC Press, 1996, pages 653-658.
            //
            //  Parameters:
            //
            //    Input/output, int *N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, double *A, the parameter of the function.
            //
            //    Output, int *X, the argument of the function.
            //
            //    Output, double *FX, the value of the function.
            //
        {
            int N_MAX = 21;

            double[] a_vec =  {
                0.02E+00, 0.10E+00, 0.10E+00, 0.50E+00,
                0.50E+00, 0.50E+00, 1.00E+00, 1.00E+00,
                1.00E+00, 1.00E+00, 2.00E+00, 2.00E+00,
                2.00E+00, 2.00E+00, 5.00E+00, 5.00E+00,
                5.00E+00, 5.00E+00, 5.00E+00, 5.00E+00,
                5.00E+00
            }
            ;
            double[] fx_vec =  {
                0.980E+00, 0.905E+00, 0.995E+00, 0.607E+00,
                0.910E+00, 0.986E+00, 0.368E+00, 0.736E+00,
                0.920E+00, 0.981E+00, 0.135E+00, 0.406E+00,
                0.677E+00, 0.857E+00, 0.007E+00, 0.040E+00,
                0.125E+00, 0.265E+00, 0.441E+00, 0.616E+00,
                0.762E+00
            }
            ;
            int[] x_vec =  {
                0, 0, 1, 0,
                1, 2, 0, 1,
                2, 3, 0, 1,
                2, 3, 0, 1,
                2, 3, 4, 5,
                6
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
                a = 0.0E+00;
                x = 0;
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