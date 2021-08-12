using System;
using Burkardt.TestValues;

namespace TestValuesTest
{
    public class StirlingTest
    {
        public static void stirling1_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    STIRLING1_VALUES_TEST tests STIRLING1_VALUES.
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
            int m = 0;
            int n = 0;
            int n_data;
            int s1 = 0;
            Console.WriteLine("");
            Console.WriteLine("STIRLING1_VALUES_TEST:");
            Console.WriteLine("  STIRLING1_VALUES returns values of");
            Console.WriteLine("  the Stirling numbers of the first kind.");
            Console.WriteLine("");
            Console.WriteLine("     N     N        S1");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Stirling.stirling1_values(ref n_data, ref n, ref m, ref s1);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + n.ToString().PadLeft(6) + "  "
                                  + m.ToString().PadLeft(6) + "  "
                                  + s1.ToString().PadLeft(12) + "");
            }
        }

        public static void stirling2_values_test()
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    STIRLING2_VALUES_TEST tests STIRLING2_VALUES.
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
            int m = 0;
            int n = 0;
            int n_data;
            int s2 = 0;
            Console.WriteLine("");
            Console.WriteLine("STIRLING2_VALUES_TEST:");
            Console.WriteLine("  STIRLING2_VALUES returns values of");
            Console.WriteLine("  the Stirling numbers of the second kind.");
            Console.WriteLine("");
            Console.WriteLine("     N     N        S2");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Stirling.stirling1_values(ref n_data, ref n, ref m, ref s2);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + n.ToString().PadLeft(6) + "  "
                                  + m.ToString().PadLeft(6) + "  "
                                  + s2.ToString().PadLeft(12) + "");
            }
        }

    }
}