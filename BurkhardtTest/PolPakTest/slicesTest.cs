using System;
using Burkardt.Function;
using Burkardt.Types;

namespace PolPakTest;

public static class slicesTest
{
    public static void slices_test ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SLICES_TEST tests SLICES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 August 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int DIM_MAX = 5;
        const int SLICE_MAX = 8;

        int dim_num;
        int[] p = new int[DIM_MAX*SLICE_MAX];

        Console.WriteLine("");
        Console.WriteLine("SLICES_TEST:");
        Console.WriteLine("  SLICES determines the maximum number of pieces created");
        Console.WriteLine("  by SLICE_NUM slices in a DIM_NUM space.");

        for ( dim_num = 1; dim_num <= DIM_MAX; dim_num++ )
        {
            int slice_num;
            for ( slice_num = 1; slice_num <= SLICE_MAX; slice_num++ )
            {
                int piece_num = Slices.slices ( dim_num, slice_num );
                p[dim_num-1+(slice_num-1)*DIM_MAX] = piece_num;
            }
        }

        typeMethods.i4mat_print ( DIM_MAX, SLICE_MAX, p, "  Slice Array:" );
    }

}