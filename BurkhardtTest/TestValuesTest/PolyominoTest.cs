using System;
using Burkardt.Values;

namespace TestValuesTest;

public class PolyominoTest
{
    public static void polyomino_chiral_count_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYOMINO_CHIRAL_COUNT_VALUES_TEST tests POLYOMINO_CHIRAL_COUNT_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 May 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n_data;
        long number = 0;
        int order = 0;
        Console.WriteLine("");
        Console.WriteLine("POLYOMINO_CHIRAL_COUNT_VALUES_TEST:");
        Console.WriteLine("  POLYOMINO_CHIRAL_COUNT_VALUES returns the number of chiral");
        Console.WriteLine("  polyominoes of a given order;");
        n_data = 0;
        Console.WriteLine("");
        Console.WriteLine("   Order      Number");
        Console.WriteLine("");
        for (;;)
        {
            Polyomino.polyomino_chiral_count_values(ref n_data, ref order, ref number);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  " + order.ToString(CultureInfo.InvariantCulture).PadLeft(4) + order
                              + "  " + number.ToString(CultureInfo.InvariantCulture).PadLeft(24) + number + "");
        }
    }

    public static void polyomino_fixed_count_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYOMINO_FIXED_COUNT_VALUES_TEST tests POLYOMINO_FIXED_COUNT_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 April 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n_data;
        long number = 0;
        int order = 0;
        Console.WriteLine("");
        Console.WriteLine("POLYOMINO_FIXED_COUNT_VALUES_TEST:");
        Console.WriteLine("  POLYOMINO_FIXED_COUNT_VALUES returns the number of fixed");
        Console.WriteLine("  polyominoes of a given order;");
        n_data = 0;
        Console.WriteLine("");
        Console.WriteLine("   Order      Number");
        Console.WriteLine("");
        for (;;)
        {
            Polyomino.polyomino_fixed_count_values(ref n_data, ref order, ref number);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  " + order.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + number.ToString(CultureInfo.InvariantCulture).PadLeft(24) + "");
        }
    }

    public static void polyomino_free_count_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYOMINO_FREE_COUNT_VALUES_TEST tests POLYOMINO_FREE_COUNT_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 April 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n_data;
        long number = 0;
        int order = 0;
        Console.WriteLine("");
        Console.WriteLine("POLYOMINO_FREE_COUNT_VALUES_TEST:");
        Console.WriteLine("  POLYOMINO_FREE_COUNT_VALUES returns the number of free");
        Console.WriteLine("  polyominoes of a given order;");
        n_data = 0;
        Console.WriteLine("");
        Console.WriteLine("   Order      Number");
        Console.WriteLine("");
        for (;;)
        {
            Polyomino.polyomino_free_count_values(ref n_data, ref order, ref number);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  " + order.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + number.ToString(CultureInfo.InvariantCulture).PadLeft(24) + "");
        }
    }
}