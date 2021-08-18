using System;
using Burkardt.Values;

namespace TestValuesTest
{
    public class SphericalHarmonicTest
    {
        public static void spherical_harmonic_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SPHERICAL_HARMONIC_VALUES_TEST tests SPHERICAL_HARMONIC_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    08 February 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int l = 0;
            int m = 0;
            int n_data;
            double phi = 0;
            double theta = 0;
            double yi = 0;
            double yr = 0;
            Console.WriteLine("");
            Console.WriteLine("SPHERICAL_HARMONIC_VALUES_TEST:");
            Console.WriteLine("  SPHERICAL_HARMONIC_VALUES stores values of");
            Console.WriteLine("  the spherical harmonic function.");
            Console.WriteLine("");
            Console.WriteLine("   L   M    THETA       PHI           Yr                    Yi");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                SphericalHarmonic.spherical_harmonic_values(ref n_data, ref l, ref m, ref theta, ref phi, ref yr,
                    ref yi);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + l.ToString().PadLeft(2) + "  "
                                  + m.ToString().PadLeft(2) + "  "
                                  + theta.ToString("0.####").PadLeft(8) + "  "
                                  + phi.ToString("0.####").PadLeft(8) + "  "
                                  + yr.ToString("0.################").PadLeft(24) + "  "
                                  + yi.ToString("0.################").PadLeft(24) + "");
            }
        }
    }
}