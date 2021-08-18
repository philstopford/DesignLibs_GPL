namespace Burkardt.Values
{
    public static class Intgr
    {
        public static void int_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    INT_VALUES returns some values of the "integer part" function.
            //
            //  Discussion:
            //
            //    INT(X) = returns the integer part of a real number.
            //
            //    The result is returned as a real number.
            //
            //    The result is computed by rounding the absolute value of the
            //   input towards 0, and then restoring the sign.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    02 October 2011
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
            //    Output double FX, the value of the function.
            //
        {
            int N_MAX = 25;

            double[] fx_vec =
            {
                -2.00E+00,
                -1.00E+00,
                -1.00E+00,
                -1.00E+00,
                -1.00E+00,
                -1.00E+00,
                0.00E+00,
                0.00E+00,
                0.00E+00,
                0.00E+00,
                0.00E+00,
                0.00E+00,
                0.00E+00,
                0.00E+00,
                0.00E+00,
                0.00E+00,
                0.00E+00,
                0.00E+00,
                0.00E+00,
                1.00E+00,
                1.00E+00,
                1.00E+00,
                1.00E+00,
                1.00E+00,
                2.00E+00
            };

            double[] x_vec =
            {
                -2.01E+00,
                -1.99E+00,
                -1.50E+00,
                -1.10E+00,
                -1.01E+00,
                -1.00E+00,
                -0.99E+00,
                -0.90E+00,
                -0.51E+00,
                -0.50E+00,
                -0.49E+00,
                -0.01E+00,
                0.00E+00,
                0.01E+00,
                0.49E+00,
                0.50E+00,
                0.51E+00,
                0.90E+00,
                0.99E+00,
                1.00E+00,
                1.01E+00,
                1.10E+00,
                1.50E+00,
                1.99E+00,
                2.01E+00
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
    }
}