using System;
using Burkardt.TestValues;

namespace TestValuesTest
{
    public class PsiTest
    {
        public static void psi_values_test()
            //****************************************************************************80
            //
            //  Purpose: 
            //
            //    PSI_VALUES_TEST tests PSI_VALUES.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    09 February 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double fx = 0;
            int n_data;
            double x = 0;
            Console.WriteLine("");
            Console.WriteLine("PSI_VALUES_TEST");
            Console.WriteLine("  PSI_VALUES stores values of");
            Console.WriteLine("  the PSI function.");
            Console.WriteLine("");
            Console.WriteLine("      X            PSI(X)");
            Console.WriteLine("");
            n_data = 0;
            for (;;)
            {
                Psi.psi_values(ref n_data, ref x, ref fx);
                if (n_data == 0)
                {
                    break;
                }

                Console.WriteLine("  "
                                  + x.ToString().PadLeft(12) + "  "
                                  + fx.ToString("0.################").PadLeft(24) + "");
            }
        }
    }
}