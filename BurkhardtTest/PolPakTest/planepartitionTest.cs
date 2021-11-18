using System;
using Burkardt.Function;

namespace PolPakTest;

public static class planepartitionTest
{
    public static void plane_partition_num_test ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PLANE_PARTITION_NUM_TEST tests PLANE_PARTITION_NUM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    24 February 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n;

        Console.WriteLine("");
        Console.WriteLine("PLANE_PARTITION_NUM_TEST");
        Console.WriteLine("  PLANE_PARTITION_NUM computes the number of plane");
        Console.WriteLine("  partitions of an integer.");
        Console.WriteLine("");

        for ( n = 1; n <= 10; n++ )
        {
            Console.WriteLine("  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + Partition.plane_partition_num ( n ).ToString(CultureInfo.InvariantCulture).PadLeft(6) + "");
        }

    }
}