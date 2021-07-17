using System;
using Burkardt.PolynomialNS;

namespace LegendreShiftedPolynomialTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for LEGENDRE_SHIFTED_POLYNOMIAL_TEST.
            //
            //  Discussion:
            //
            //    LEGENDRE_SHIFTED_POLYNOMIAL_TEST tests the LEGENDRE_SHIFTED_POLYNOMIAL library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 March 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("LEGENDRE_SHIFTED_POLYNOMIAL_TEST:");
            Console.WriteLine("  Test the LEGENDRE_SHIFTED_POLYNOMIAL library.");

            p01_polynomial_value_test();
            //
            //  Terminate.
            //
            Console.WriteLine("");
            Console.WriteLine("LEGENDRE_SHIFTED_POLYNOMIAL_TEST:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void p01_polynomial_value_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    P01_POLYNOMIAL_VALUE_TEST tests P01_POLYNOMIAL_VALUE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 March 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int n_data;
            double e;
            double fx1 = 0;
            double fx2;
            double[] fx2_vec;
            int n = 0;
            double x = 0;
            double[] x_vec = new double[1];

            Console.WriteLine("");
            Console.WriteLine("P01_POLYNOMIAL_VALUE_TEST:");
            Console.WriteLine("  P01_POLYNOMIAL_VALUE evaluates the shifted Legendre polynomial P01(n,x).");
            Console.WriteLine("");
            Console.WriteLine("                        Tabulated                 Computed");
            Console.WriteLine("     N        X          P01(N,X)                 P01(N,X)                     Error");
            Console.WriteLine("");

            n_data = 0;

            for (;;)
            {
                LegendreShifted.p01_polynomial_values(ref n_data, ref n, ref x, ref fx1);

                if (n_data == 0)
                {
                    break;
                }

                x_vec[0] = x;
                fx2_vec = LegendreShifted.p01_polynomial_value(1, n, x_vec);
                fx2 = fx2_vec[n];

                e = fx1 - fx2;

                Console.WriteLine("  " + n.ToString().PadLeft(4)
                                       + "  " + x.ToString().PadLeft(12)
                                       + "  " + fx1.ToString().PadLeft(24)
                                       + "  " + fx2.ToString().PadLeft(24)
                                       + "  " + e.ToString().PadLeft(8) + "");
            }
        }
    }
}