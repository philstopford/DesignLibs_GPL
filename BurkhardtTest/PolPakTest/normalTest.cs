using System;
using Burkardt.CDFLib;

namespace PolPakTest
{
    public static class normalTest
    {
        public static void normal_01_cdf_inverse_test ( )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    NORMAL_01_CDF_INVERSE_TEST tests NORMAL_01_CDF_INVERSE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 February 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double fx = 0;
            int n_data = 0;
            double x = 0;
            double x2 = 0;

            Console.WriteLine("");
            Console.WriteLine("NORMAL_01_CDF_INVERSE_TEST:");
            Console.WriteLine("  NORMAL_01_CDF_INVERSE inverts the normal 01 CDF.");
            Console.WriteLine("");
            Console.WriteLine("    FX      X    NORMAL_01_CDF_INVERSE(FX)");
            Console.WriteLine("");

            n_data = 0;

            for ( ; ; )
            {
                Burkardt.Values.Normal.normal_01_cdf_values ( ref n_data, ref x, ref fx );

                if ( n_data == 0 )
                {
                    break;
                }

                x2 = CDF.normal_01_cdf_inv ( fx );

                Console.WriteLine("  "
                                        + fx.ToString().PadLeft(8)   + "  "
                                        + x.ToString().PadLeft(14)  + "  "
                                        + x2.ToString().PadLeft(14) + "");
            }

        }
    }
}