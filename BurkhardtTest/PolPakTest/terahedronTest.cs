using System;
using Burkardt.TetrahedronNS;

namespace PolPakTest;

public static class terahedronTest
{
    public static void tetrahedron_num_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TETRAHEDRON_NUM_TEST tests TETRAHEDRON_NUM.
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
        int n;

        Console.WriteLine("");
        Console.WriteLine("TETRAHEDRON_NUM_TEST");
        Console.WriteLine("  TETRAHEDRON_NUM computes the tetrahedron numbers.");
        Console.WriteLine("");

        for (n = 1; n <= 10; n++)
        {
            Console.WriteLine("  "
                              + n.ToString().PadLeft(4) + "  "
                              + Tetrahedron.tetrahedron_num(n).ToString().PadLeft(6) + "");
        }

    }
}