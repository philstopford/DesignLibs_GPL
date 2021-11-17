using System;
using Burkardt.AppliedStatistics;

namespace TOMS291Test;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for TOMS291_TEST.
        //
        //  Discussion:
        //
        //    TOMS291_TEST tests the TOMS291 library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 January 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("TOMS291_TEST:");
        Console.WriteLine("  Test the TOMS291 library.");

        test01();

        Console.WriteLine("");
        Console.WriteLine("TOMS291_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 demonstrates the use of ALOGAM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 January 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        double fx2;
        int ifault = 0;
        int n_data;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST01:");
        Console.WriteLine("  ALOGAM computes the logarithm of the ");
        Console.WriteLine("  Gamma function.  We compare the result");
        Console.WriteLine("  to tabulated values.");
        Console.WriteLine("");
        Console.WriteLine("          X                     "
                          + "FX                        FX2");
        Console.WriteLine("                                "
                          + "(Tabulated)               (ALOGAM)                DIFF");
        Console.WriteLine("");

        n_data = 0;

        for (;;)
        {
            Burkardt.Values.Gamma.gamma_log_values(ref n_data, ref x, ref fx);

            if (n_data == 0)
            {
                break;
            }

            fx2 = Algorithms.alogam(x, ref ifault);

            Console.WriteLine("  " + x.ToString("0.################").PadLeft(24)
                                   + "  " + fx.ToString("0.################").PadLeft(24)
                                   + "  " + fx2.ToString("0.################").PadLeft(24)
                                   + "  " + Math.Abs(fx - fx2).ToString("0.####").PadLeft(10) + "");
        }

    }
}