using System;
using Burkardt;
using Burkardt.AppliedStatistics;

namespace BetaNCTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for BETA_NC_TEST.
            //
            //  Discussion:
            //
            //    BETA_NC_TEST tests the BETA_NC library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    29 January 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("BETA_NC_TEST:");
            Console.WriteLine("  Test the BETA_NC library.");

            test01();
            //
            //  Terminate.
            //
            Console.WriteLine("");
            Console.WriteLine("BETA_NC_TEST:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");

        }

        static void test01()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST01 demonstrates the use of BETA_NONCENTRAL_CDF.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    29 January 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double a = 0;
            double b = 0;
            double fx = 0;
            double lambda = 0;
            double x = 0;

            double error_max = 1.0E-10;

            Console.WriteLine("");
            Console.WriteLine("TEST01:");
            Console.WriteLine("  BETA_NONCENTRAL_CDF computes the noncentral incomplete ");
            Console.WriteLine("  Beta function.");
            Console.WriteLine("  Compare to tabulated values.");
            Console.WriteLine("");
            Console.WriteLine("      A      B     LAMBDA        X      "
                + "    FX                        FX2");
            Console.WriteLine("                                        "
                + "    (Tabulated)               (computed)          DIFF");
            Console.WriteLine("");

            int n_data = 0;

            for (;;)
            {
                Algorithms.beta_noncentral_cdf_values(ref n_data, ref a, ref b, ref lambda, ref x, ref fx);

                if (n_data == 0)
                {
                    break;
                }

                double fx2 = BetaNC.beta_noncentral_cdf(a, b, lambda, x, error_max);

                Console.WriteLine("  " + a.ToString("0.##").PadLeft(5)
                    + "  " + b.ToString("0.##").PadLeft(5)
                    + "  " + lambda.ToString("0.###").PadLeft(7)
                    + "  " + x.ToString("0.####").PadLeft(10)
                    + "  " + fx.ToString("0.################").PadLeft(24)
                    + "  " + fx2.ToString("0.################").PadLeft(24)
                    + "  " + Math.Abs(fx - fx2).ToString("0.####").PadLeft(10) + "");
            }
        }
    }
}