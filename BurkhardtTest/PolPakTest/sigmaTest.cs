using System;
using Burkardt.Function;

namespace PolPakTest;

public static class sigmaTest
{
    public static void sigma_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SIGMA_TEST tests SIGMA.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 June 2007
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
        Console.WriteLine("SIGMA_TEST");
        Console.WriteLine("  SIGMA computes the SIGMA function.");
        Console.WriteLine("");
        Console.WriteLine("  N   Exact   SIGMA(N)");
        Console.WriteLine("");

        n_data = 0;

        for (;;)
        {
            Burkardt.Values.Sigma.sigma_values(ref n_data, ref n, ref c);

            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + n.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  "
                              + c.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "  "
                              + Sigma.sigma(n).ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

    }

}