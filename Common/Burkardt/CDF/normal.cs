namespace Burkardt.CDFLib
{
    public static partial class CDF
    {
        public static void normal_cdf_values ( ref int n_data, ref double x, ref double fx )

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

            double[] fx_vec = {
                0.500000000000000E+00, 0.539827837277029E+00, 0.579259709439103E+00,
                0.617911422188953E+00, 0.655421741610324E+00, 0.691462461274013E+00,
                0.725746882249927E+00, 0.758036347776927E+00, 0.788144601416604E+00,
                0.815939874653241E+00, 0.841344746068543E+00, 0.933192798731142E+00,
                0.977249868051821E+00 };
            double[] x_vec = {
                0.00E+00, 0.10E+00, 0.20E+00,
                0.30E+00, 0.40E+00, 0.50E+00,
                0.60E+00, 0.70E+00, 0.80E+00,
                0.90E+00, 1.00E+00, 1.50E+00,
                2.00E+00 };

            if ( n_data < 0 )
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if ( N_MAX < n_data )
            {
                n_data = 0;
                x = 0.0E+00;
                fx = 0.0E+00;
            }
            else
            {
                x = x_vec[n_data-1];
                fx = fx_vec[n_data-1];
            }
        }
    }
}