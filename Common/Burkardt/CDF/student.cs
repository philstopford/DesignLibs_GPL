namespace Burkardt.CDFLib
{
    public static partial class CDF
    {
        public static void student_cdf_values(ref int n_data, ref int a, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    STUDENT_CDF_VALUES returns some values of the Student CDF.
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
            //    Output, int *A, the parameter of the function.
            //
            //    Output, double *X, the argument of the function.
            //
            //    Output, double *FX, the value of the function.
            //
        {
            int N_MAX = 13;

            int[] a_vec =  {
                1, 2, 3, 4,
                5, 2, 5, 2,
                5, 2, 3, 4,
                5
            }
            ;
            double[] fx_vec =  {
                0.60E+00, 0.60E+00, 0.60E+00, 0.60E+00,
                0.60E+00, 0.75E+00, 0.75E+00, 0.95E+00,
                0.95E+00, 0.99E+00, 0.99E+00, 0.99E+00,
                0.99E+00
            }
            ;
            double[] x_vec =  {
                0.325E+00, 0.289E+00, 0.277E+00, 0.271E+00,
                0.267E+00, 0.816E+00, 0.727E+00, 2.920E+00,
                2.015E+00, 6.965E+00, 4.541E+00, 3.747E+00,
                3.365E+00
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