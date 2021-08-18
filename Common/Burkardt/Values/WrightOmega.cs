using System.Numerics;

namespace Burkardt.Values
{
    public static class WrightOmega
    {
        public static void wright_omega_values(ref int n_data, ref Complex z,
                ref Complex fz)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    WRIGHT_OMEGA_VALUES returns some values of the Wright Omega function.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 May 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Robert Corless, David Jeffrey,
            //    The Wright Omega Function,
            //    in Artificial Intelligence, Automated Reasoning, and Symbolic Computation,
            //    ed J Calmet, B Benhamou, O Caprotti, L Henocque, V Sorge,
            //    Lecture Notes in Artificial Intelligence, volume 2385,
            //    Springer, 2002, pages 76-89.
            //
            //  Parameters:
            //
            //    Input/output, ref int N_DATA.  The user sets N_DATA to 0 before the
            //    first call.  On each call, the routine increments N_DATA by 1, and
            //    returns the corresponding data; when there is no more data, the
            //    output value of N_DATA will be 0 again.
            //
            //    Output, new Complex &Z, the argument of the function.
            //
            //    Output, new Complex &FZ, the value of the function.
            //
        {
            int N_MAX = 10;

            Complex[] fz_vec =
            {
                new Complex(+0.5671432904097838, +0.0000000000000000),
                new Complex(+1.000000000000000, +0.0000000000000000),
                new Complex(+2.718281828459045, +0.0000000000000000),
                new Complex(-1.000000000000000, +0.0000000000000000),
                new Complex(-1.000000000000000, +0.0000000000000000),
                new Complex(-2.000000000000000, +0.0000000000000000),
                new Complex(-0.40637573995996, +0.0000000000000000),
                new Complex(+0.000000000000000, +1.0000000000000000),
                new Complex(-0.3181315052047641, +1.337235701430689),
                new Complex(+0.9372082083733697, +0.5054213160131512)
            };

            Complex[] z_vec =
            {
                new Complex(+0.000000000000000, +0.000000000000000),
                new Complex(+1.000000000000000, +0.000000000000000),
                new Complex(+3.718281828459045, +0.000000000000000),
                new Complex(-1.000000000000000, +3.141592653589793),
                new Complex(-1.000000000000000, -3.141592653589793),
                new Complex(-1.306852819440055, +3.141592653589793),
                new Complex(-1.306852819440055, -3.141592653589793),
                new Complex(+0.000000000000000, +2.570796326794897),
                new Complex(+0.000000000000000, +3.141592653589793),
                new Complex(+1.000000000000000, +1.000000000000000)
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            n_data = n_data + 1;

            if (N_MAX < n_data)
            {
                n_data = 0;
                z = new Complex(0.0, 0.0);
                fz = new Complex(0.0, 0.0);
            }
            else
            {
                z = z_vec[n_data - 1];
                fz = fz_vec[n_data - 1];
            }
        }

    }
}