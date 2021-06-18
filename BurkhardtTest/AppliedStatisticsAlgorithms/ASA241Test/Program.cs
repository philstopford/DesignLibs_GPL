using System;
using Burkardt.AppliedStatistics;

namespace ASA241Test
{
    class Program
    {
        static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for ASA241_TEST.
        //
        //  Discussion:
        //
        //    ASA241_TEST tests the ASA241 library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 January 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            Console.WriteLine("");
            Console.WriteLine("ASA241_TEST:");
            Console.WriteLine("  Test the ASA241 library.");

            test01();
            test02();

            Console.WriteLine("");
            Console.WriteLine("ASA241_TEST:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");

        }

        static void test01()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests R4_NORMAL_01_CDF_INVERSE, NORMAL_01_CDF_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 February 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            double fx = 0;
            double x = 0;

            Console.WriteLine("");
            Console.WriteLine("TEST01:");
            Console.WriteLine("  Let FX = NormalCDF ( X ).");
            Console.WriteLine("");
            Console.WriteLine("  NORMAL_01_CDF_VALUES returns some values of ( X, FX ).");
            Console.WriteLine("");
            Console.WriteLine("  R4_NORMAL_01_CDF_INVERSE takes the value of FX, and computes an");
            Console.WriteLine("    estimate X2, of the corresponding input argument,");
            Console.WriteLine("    accurate to about 7 decimal places.");
            Console.WriteLine("");
            Console.WriteLine("          FX                        X                        X2          DIFF");
            Console.WriteLine("");

            int n_data = 0;

            for (;;)
            {
                Algorithms.normal_01_cdf_values(ref n_data, ref x, ref fx);

                if (n_data == 0)
                {
                    break;
                }

                float fx2 = (float) (fx);
                float x2 = Algorithms.r4_normal_01_cdf_inverse(fx2);

                Console.WriteLine("  " + fx.ToString("0.################").PadLeft(24)
                    + "  " + x.ToString("0.################").PadLeft(24)
                    + "  " + x2.ToString("0.################").PadLeft(24)
                    + "  " + Math.Abs(x - x2).ToString("0.################") + "");
            }
        }

        static void test02()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests R8_NORMAL_01_CDF_INVERSE, NORMAL_01_CDF_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 February 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            double fx = 0;
            double x = 0;

            Console.WriteLine("");
            Console.WriteLine("TEST02:");
            Console.WriteLine("  Let FX = NormalCDF ( X ).");
            Console.WriteLine("");
            Console.WriteLine("  NORMAL_01_CDF_VALUES returns some values of ( X, FX ).");
            Console.WriteLine("");
            Console.WriteLine("  R8_NORMAL_01_CDF_INVERSE takes the value of FX, and computes an");
            Console.WriteLine("    estimate X2, of the corresponding input argument,");
            Console.WriteLine("    accurate to about 16 decimal places.");
            Console.WriteLine("");
            Console.WriteLine("          FX                        X                        X2          DIFF");
            Console.WriteLine("");

            int n_data = 0;

            for (;;)
            {
                Algorithms.normal_01_cdf_values(ref n_data, ref x, ref fx);

                if (n_data == 0)
                {
                    break;
                }

                double x2 = Algorithms.r8_normal_01_cdf_inverse(fx);

                Console.WriteLine("  " + fx.ToString("0.################").PadLeft(24)
                    + "  " + x.ToString("0.################").PadLeft(24)
                    + "  " + x2.ToString("0.################").PadLeft(24)
                    + "  " + Math.Abs(x - x2).ToString("0.################") + "");
            }

            return;
        }

    }
}