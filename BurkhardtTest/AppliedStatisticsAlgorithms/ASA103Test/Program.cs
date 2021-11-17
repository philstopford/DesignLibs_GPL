using System;
using Burkardt.AppliedStatistics;

namespace ASA103Test;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for ASA103_TEST.
        //
        //  Discussion:
        //
        //    ASA103_TEST tests the ASA103 library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 January 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("ASA103_TEST:");
        Console.WriteLine("  Test the ASA103 library.");

        digamma_test();

        Console.WriteLine("");
        Console.WriteLine("ASA103_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }


    private static void digamma_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIGAMMA_TEST tests DIGAMA.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 January 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double fx = 0;
        int ifault = 0;
        double x = 0;

        Console.WriteLine("");
        Console.WriteLine("DIGAMMA_TEST:");
        Console.WriteLine("  DIGAMMA computes the Digamma or Psi function. ");
        Console.WriteLine("  Compare the result to tabulated values.");
        Console.WriteLine("");
        Console.WriteLine("          X       "
                          + "FX                        FX2");
        Console.WriteLine("                  "
                          + "(Tabulated)               (DIGAMMA)               DIFF");
        Console.WriteLine("");

        int n_data = 0;

        for (;;)
        {
            Algorithms.psi_values(ref n_data, ref x, ref fx);

            if (n_data == 0)
            {
                break;
            }

            double fx2 = Algorithms.digamma(x, ref ifault);

            Console.WriteLine("  " + x.ToString("0.####").PadLeft(10)
                                   + "  " + fx.ToString("0.################").PadLeft(24)
                                   + "  " + fx2.ToString("0.################").PadLeft(24)
                                   + "  " + Math.Abs(fx - fx2).ToString("0.####").PadLeft(10) + "");
        }
    }

}