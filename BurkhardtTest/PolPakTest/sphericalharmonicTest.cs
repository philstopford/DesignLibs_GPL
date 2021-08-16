using System;
using Burkardt.Function;

namespace PolPakTest
{
    public static class sphericalharmonicTest
    {
        public static void spherical_harmonic_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERICAL_HARMONIC_TEST tests SPHERICAL_HARMONIC.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 June 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N_MAX = 20;

            double[] c = new double[N_MAX + 1];
            int l = 0;
            int m = 0;
            int n_data = 0;
            double phi = 0;
            double[] s = new double[N_MAX + 1];
            double theta = 0;
            double yi = 0;
            double yi2 = 0;
            double yr = 0;
            double yr2 = 0;

            Console.WriteLine("");
            Console.WriteLine("SPHERICAL_HARMONIC_TEST:");
            Console.WriteLine("  SPHERICAL_HARMONIC evaluates spherical harmonic functions.");
            Console.WriteLine("");
            Console.WriteLine(
                "         N         M    THETA      PHI            YR            YI");
            Console.WriteLine("");

            n_data = 0;

            for (;;)
            {
                Burkardt.TestValues.SphericalHarmonic.spherical_harmonic_values(ref n_data, ref l, ref m, ref theta,
                    ref phi, ref yr, ref yi);

                if (n_data == 0)
                {
                    break;
                }

                SphericalHarmonic.spherical_harmonic(l, m, theta, phi, ref c, ref s);

                yr2 = c[l];
                yi2 = s[l];

                Console.WriteLine("  "
                                  + l.ToString().PadLeft(8) + "  "
                                  + m.ToString().PadLeft(8) + "  "
                                  + theta.ToString().PadLeft(8) + "  "
                                  + phi.ToString().PadLeft(8) + "  "
                                  + yr.ToString().PadLeft(14) + "  "
                                  + yi.ToString().PadLeft(14) + "");

                Console.WriteLine("  "
                                  + "        " + "  "
                                  + "        " + "  "
                                  + "        " + "  "
                                  + "        " + "  "
                                  + yr2.ToString().PadLeft(14) + "  "
                                  + yi2.ToString().PadLeft(14) + "");
            }

        }

    }
}