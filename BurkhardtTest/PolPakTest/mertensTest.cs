using System;
using Burkardt.Function;

namespace PolPakTest;

public static class mertensTest
{
    public static void mertens_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MERTENS_TEST tests MERTENS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 October 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int c = 0;
        int n = 0;
        int n_data = 0;

        Console.WriteLine("");
        Console.WriteLine("MERTENS_TEST");
        Console.WriteLine("  MERTENS computes the Mertens function.");
        Console.WriteLine("");
        Console.WriteLine("      N   Exact   MERTENS(N)");
        Console.WriteLine("");

        n_data = 0;

        for (;;)
        {
            Burkardt.Values.Mertens.mertens_values(ref n_data, ref n, ref c);

            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + n.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                              + c.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  "
                              + Mertens.mertens(n).ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

    }
}