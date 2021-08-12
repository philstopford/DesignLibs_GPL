namespace Burkardt.TestValues
{
    public static class Sqrt
    {

        public static void sqrt_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SQRT_VALUES returns some values of the square root function.
            //
            //  Discussion:
            //
            //    SQRT(X) = positive real number Y such that Y * Y = X.
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      Sqrt[x]
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
            //    Output, ref double X, the argument of the function.
            //
            //    Output double FX, the value of the function.
            //
        {
            int N_MAX = 14;

            double[] fx_vec =
            {
                0.0000000000000000E+00,
                0.9000000040950000E-04,
                0.3000000000000000E+00,
                0.3162277660168379E+00,
                0.6324555320336759E+00,
                0.1000000000000000E+01,
                0.1414213562373095E+01,
                0.1732050807568877E+01,
                0.1772453850905516E+01,
                0.4358898943540674E+01,
                0.5385164807134504E+01,
                0.8426149773176359E+01,
                0.9848857801796105E+01,
                0.1111111106055556E+05
            };

            double[] x_vec =
            {
                0.0000000000000000E+00,
                0.8100000073710001E-08,
                0.9000000000000000E-01,
                0.1000000000000000E+00,
                0.4000000000000000E+00,
                0.1000000000000000E+01,
                0.2000000000000000E+01,
                0.3000000000000000E+01,
                0.3141592653589793E+01,
                0.1900000000000000E+02,
                0.2900000000000000E+02,
                0.7100000000000000E+02,
                0.9700000000000000E+02,
                0.1234567890000000E+09
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