using System;
using Burkardt.Probability;

namespace Burkardt.ProbabilityTest
{
    partial class Program
    {
        static void bessel_i0_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST0105 tests BESSEL_I0.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 August 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double fx = 0;
            int n_data = 0;
            double x = 0;

            Console.WriteLine("");
            Console.WriteLine("BESSEL_I0_TEST:");
            Console.WriteLine("  BESSEL_I0 evaluates the Bessel function of the");
            Console.WriteLine("  first kind and order 0;");
            Console.WriteLine("");
            Console.WriteLine("      X       Exact F       BESSEL_I0(X)");
            Console.WriteLine("");

            n_data = 0;

            for (;;)
            {
                Bessel.bessel_i0_values(ref n_data, ref x, ref fx);

                if (n_data == 0)
                {
                    break;
                }

                double fx2 = Bessel.bessel_i0(x);

                Console.WriteLine("  "
                                  + x.ToString().PadLeft(8) + "  "
                                  + fx.ToString().PadLeft(16) + "  "
                                  + fx2.ToString().PadLeft(16) + "");
            }

            return;
        }

        static void bessel_i1_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BESSEL_I1_TEST tests BESSEL_I1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 August 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double fx = 0;
            double x = 0;

            Console.WriteLine("");
            Console.WriteLine("BESSEL_I1_TEST:");
            Console.WriteLine("  BESSEL_I1 evaluates the Bessel function of the");
            Console.WriteLine("  first kind and order 1;");
            Console.WriteLine("");
            Console.WriteLine("      X       Exact F       BESSEL_I1(X)");
            Console.WriteLine("");

            int n_data = 0;

            for (;;)
            {
                Bessel.bessel_i1_values(ref n_data, ref x, ref fx);

                if (n_data == 0)
                {
                    break;
                }

                double fx2 = Bessel.bessel_i1(x);

                Console.WriteLine("  "
                                  + x.ToString().PadLeft(8) + "  "
                                  + fx.ToString().PadLeft(16) + "  "
                                  + fx2.ToString().PadLeft(16) + "");
            }

            return;
        }
    }
}