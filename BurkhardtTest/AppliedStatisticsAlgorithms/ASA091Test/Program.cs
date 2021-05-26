using System;
using Burkardt.AppliedStatistics;

namespace Burkardt.ASA091Test
{
    class Program
    {
        static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for ASA091_TEST.
        //
        //  Discussion:
        //
        //    ASA091_TEST tests the ASA091 library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 February 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            Console.WriteLine("");
            Console.WriteLine("ASA091_TEST:");
            Console.WriteLine("  C++ version");
            Console.WriteLine("  Test the ASA091 library.");

            test01();
            test02();

            Console.WriteLine("");
            Console.WriteLine("ASA091_TEST:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void test01()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 makes a single simple calculation with PPCHI2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 February 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            double g;
            int ifault = 0;
            double p;
            double v;
            double value;
            double value_correct = 0.4;

            p = 0.017523;
            v = 4.0;

            Console.WriteLine("");
            Console.WriteLine("TEST01:");
            Console.WriteLine("  Perform a simple sample calculation using");
            Console.WriteLine("  PPCHI2 to invert the Chi-Squared CDF.");

            g = Helpers.LogGamma(v / 2.0);

            Console.WriteLine("");
            Console.WriteLine("  P =                  "
                 + p.ToString("0.################").PadLeft(24) + "");
            Console.WriteLine("  V =                  "
                 + v.ToString("0.################").PadLeft(24) + "");
            Console.WriteLine("  G Log(Gamma(V/2)) =  "
                 + g.ToString("0.################").PadLeft(24) + "");

            value = Algorithms.ppchi2(p, v, g, ref ifault);

            Console.WriteLine("  VALUE =              "
                 + value.ToString("0.################").PadLeft(24) + "");
            Console.WriteLine("  VALUE (correct) =    "
                 + value_correct.ToString("0.################").PadLeft(24) + "");

            Console.WriteLine("");
            Console.WriteLine("  Error flag IFAULT = " + ifault + "");

            return;
        }

        static void test02()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 compares PPCHI2 against tabulated values.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 February 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            int a = 0;
            double fx = 0;
            double g;
            int ifault = 0;
            int n_data;
            double v;
            double x = 0;
            double x2;

            Console.WriteLine("");
            Console.WriteLine("TEST02:");
            Console.WriteLine("  PPCHI2 computes percentage points of the Chi-Square CDF.");
            Console.WriteLine("  Compare to tabulated values.");
            Console.WriteLine("");
            Console.WriteLine("         N        CDF       X                        "
                + " X2                    DIFF");
            Console.WriteLine("                           (tabulated)               "
                + "(PPCHI2)");
            Console.WriteLine("");

            n_data = 0;

            for (;;)
            {
                Algorithms.chi_square_cdf_values(ref n_data, ref a, ref x, ref fx);

                if (n_data == 0)
                {
                    break;
                }

                v = a;

                g = Helpers.LogGamma(v / 2.0);

                x2 = Algorithms.ppchi2(fx, v, g, ref ifault);

                Console.WriteLine("  " + a.ToString("0.####").PadLeft(10)
                    + "  " + fx.ToString("0.####").PadLeft(10)
                    + "  " + x.ToString("0.################").PadLeft(24)
                    + "  " + x2.ToString("0.################").PadLeft(24)
                    + "  " + Math.Abs((x - x2)).ToString("0.####").PadLeft(10) + "");
            }
        }

    }
}