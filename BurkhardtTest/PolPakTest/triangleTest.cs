using System;
using Burkardt.TriangleNS;
using Burkardt.Types;

namespace PolPakTest;

public static class triangleTest
{
    public static void triangle_num_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_NUM_TEST tests TRIANGLE_NUM.
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
        Console.WriteLine("TRIANGLE_NUM_TEST");
        Console.WriteLine("  TRIANGLE_NUM computes the triangular numbers.");
        Console.WriteLine("");

        for (n = 1; n <= 10; n++)
        {
            Console.WriteLine("  "
                              + n.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "  "
                              + Triangle.triangle_num(n).ToString(CultureInfo.InvariantCulture).PadLeft(6) + "");
            ;
        }

    }

    public static void triangle_lower_to_i4_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_LOWER_TO_I4_TEST tests TRIANGLE_LOWER_TO_I4.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 April 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int j;
        int k;

        Console.WriteLine("");
        Console.WriteLine("TRIANGLE_LOWER_TO_I4_TEST");
        Console.WriteLine("  TRIANGLE_LOWER_TO_I4 converts a lower triangular index to a");
        Console.WriteLine("  linear one.");
        Console.WriteLine("");
        Console.WriteLine("     I     J ==>   K");
        Console.WriteLine("");

        for (i = 0; i <= 4; i++)
        {
            for (j = 0; j <= i; j++)
            {
                k = typeMethods.triangle_lower_to_i4(i, j);

                Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                       + "  " + j.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                       + "    " + k.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "");
            }
        }

    }

}