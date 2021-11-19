using System;
using Burkardt.SolveNS;
using Burkardt.Types;

namespace SubsetTestNS;

public static class PellTest
{
    public static void pell_basic_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PELL_BASIC_TEST tests PELL_BASIC.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 March 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int d = 0;
        int q = 0;
        int r = 0;
        int x0 = 0;
        int y0 = 0;

        Console.WriteLine("");
        Console.WriteLine("PELL_BASIC_TEST");
        Console.WriteLine("  PELL_BASIC solves the basic Pell equation.");

        Console.WriteLine("");
        Console.WriteLine("       D       X        Y         R");
        Console.WriteLine("");

        for (d = 2; d <= 20; d++)
        {
            typeMethods.i4_sqrt(d, ref q, ref r);

            if (r != 0)
            {
                Pell.pell_basic(d, ref x0, ref y0);

                r = x0 * x0 - d * y0 * y0;

                Console.WriteLine(d.ToString().PadLeft(9) + "  "
                                                          + x0.ToString().PadLeft(9) + "  "
                                                          + y0.ToString().PadLeft(9) + "  "
                                                          + r.ToString().PadLeft(9) + "");
            }
        }
    }

    public static void pell_next_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PELL_NEXT_TEST tests PELL_NEXT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 March 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int d = 0;
        int q = 0;
        int r = 0;
        int x0 = 0;
        int x1 = 0;
        int y0 = 0;
        int y1 = 0;

        Console.WriteLine("");
        Console.WriteLine("PELL_NEXT_TEST");
        Console.WriteLine("  PELL_NEXT computes the next solution.");
        Console.WriteLine("");
        Console.WriteLine("       D       X        Y         R");
        Console.WriteLine("");

        for (d = 2; d <= 20; d++)
        {
            typeMethods.i4_sqrt(d, ref q, ref r);

            if (r != 0)
            {
                Pell.pell_basic(d, ref x0, ref y0);

                r = x0 * x0 - d * y0 * y0;

                Console.WriteLine(d.ToString().PadLeft(9) + "  "
                                                          + x0.ToString().PadLeft(9) + "  "
                                                          + y0.ToString().PadLeft(9) + "  "
                                                          + r.ToString().PadLeft(9) + "");

                Pell.pell_next(d, x0, y0, x0, y0, ref x1, ref y1);

                r = x1 * x1 - d * y1 * y1;

                Console.WriteLine("         " + "  "
                                              + x1.ToString().PadLeft(9) + "  "
                                              + y1.ToString().PadLeft(9) + "  "
                                              + r.ToString().PadLeft(9) + "");

            }

        }
    }

}