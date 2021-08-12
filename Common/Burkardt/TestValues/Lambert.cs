namespace Burkardt.TestValues
{
    public static class Lambert
    {

        public static void lambert_w_values(ref int n_data, ref double x, ref double fx)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LAMBERT_W_VALUES returns some values of the Lambert W function.
            //
            //  Discussion:
            //
            //    The function W(X) is defined implicitly by:
            //
            //      W(X) * e^W(X) = X
            //
            //    The function is also known as the "Omega" function.
            //
            //    In Mathematica, the function can be evaluated by:
            //
            //      W = ProductLog [ X ]
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    24 February 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Brian Hayes,
            //    "Why W?",
            //    The American Scientist,
            //    Volume 93, March-April 2005, pages 104-108.
            //
            //    Eric Weisstein,
            //    "Lambert's W-Function",
            //    CRC Concise Encyclopedia of Mathematics,
            //    CRC Press, 1998.
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
            int N_MAX = 22;

            double[] fx_vec =
            {
                0.0000000000000000E+00,
                0.3517337112491958E+00,
                0.5671432904097839E+00,
                0.7258613577662263E+00,
                0.8526055020137255E+00,
                0.9585863567287029E+00,
                0.1000000000000000E+01,
                0.1049908894964040E+01,
                0.1130289326974136E+01,
                0.1202167873197043E+01,
                0.1267237814307435E+01,
                0.1326724665242200E+01,
                0.1381545379445041E+01,
                0.1432404775898300E+01,
                0.1479856830173851E+01,
                0.1524345204984144E+01,
                0.1566230953782388E+01,
                0.1605811996320178E+01,
                0.1745528002740699E+01,
                0.3385630140290050E+01,
                0.5249602852401596E+01,
                0.1138335808614005E+02
            };
            double[] x_vec =
            {
                0.0000000000000000E+00,
                0.5000000000000000E+00,
                0.1000000000000000E+01,
                0.1500000000000000E+01,
                0.2000000000000000E+01,
                0.2500000000000000E+01,
                0.2718281828459045E+01,
                0.3000000000000000E+01,
                0.3500000000000000E+01,
                0.4000000000000000E+01,
                0.4500000000000000E+01,
                0.5000000000000000E+01,
                0.5500000000000000E+01,
                0.6000000000000000E+01,
                0.6500000000000000E+01,
                0.7000000000000000E+01,
                0.7500000000000000E+01,
                0.8000000000000000E+01,
                0.1000000000000000E+02,
                0.1000000000000000E+03,
                0.1000000000000000E+04,
                0.1000000000000000E+07
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