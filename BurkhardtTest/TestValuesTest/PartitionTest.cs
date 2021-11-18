using System;
using Burkardt.Values;

namespace TestValuesTest;

public class PartitionTest
{
    public static void partition_count_values_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PARTITION_COUNT_VALUES_TEST tests PARTITION_COUNT_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int fn = 0;
        int n = 0;
        int n_data;
        Console.WriteLine("");
        Console.WriteLine("PARTITION_COUNT_VALUES_TEST:");
        Console.WriteLine("  PARTITION_COUNT_VALUES returns values of ");
        Console.WriteLine("  the integer partition count function.");
        Console.WriteLine("");
        Console.WriteLine("     N         P(N)");
        Console.WriteLine("");
        n_data = 0;
        for (;;)
        {
            Partition.partition_count_values(ref n_data, ref n, ref fn);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + n.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                              + fn.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    public static void partition_distinct_count_values_test()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PARTITION_DISTINCT_COUNT_VALUES_TEST tests PARTITION_DISTINCT_COUNT_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int fn = 0;
        int n = 0;
        int n_data;
        Console.WriteLine("");
        Console.WriteLine("PARTITION_DISTINCT_COUNT_VALUES_TEST:");
        Console.WriteLine("  PARTITION_DISTINCT_COUNT_VALUES returns values of ");
        Console.WriteLine("  the integer distinct partition count function.");
        Console.WriteLine("");
        Console.WriteLine("     N         Q(N)");
        Console.WriteLine("");
        n_data = 0;
        for (;;)
        {
            Partition.partition_distinct_count_values(ref n_data, ref n, ref fn);
            if (n_data == 0)
            {
                break;
            }

            Console.WriteLine("  "
                              + n.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                              + fn.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }
}