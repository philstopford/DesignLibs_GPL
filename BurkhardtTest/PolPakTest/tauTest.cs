using System;
using Burkardt.Function;

namespace PolPakTest;

public static class tauTest
{
    public static void tau_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TAU_TEST tests TAU.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 May 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int c = 0;
        int n = 0;
        int n_data;

        Console.WriteLine("");
        Console.WriteLine("TAU_TEST");
        Console.WriteLine("  TAU computes the Tau function.");
        Console.WriteLine("");
        Console.WriteLine("  N  exact C(I)  computed C(I)");
        Console.WriteLine("");

        n_data = 0;

        for (;;)
        {
            Burkardt.Values.Tau.tau_values(ref n_data, ref n, ref c);

            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + n.ToString().PadLeft(4) + "  "
                              + c.ToString().PadLeft(10) + "  "
                              + Tau.tau(n).ToString().PadLeft(10) + "");
        }

    }

}