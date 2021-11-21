using System;
using Burkardt.Types;

namespace Burkardt.Function;

public static class Slices
{
    public static int slices ( int dim_num, int slice_num )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SLICES: maximum number of pieces created by a given number of slices.
        //
        //  Discussion:
        //
        //    If we imagine slicing a pizza, each slice produce more pieces.  
        //    The position of the slice affects the number of pieces created, but there
        //    is a maximum.  
        //
        //    This function determines the maximum number of pieces created by a given
        //    number of slices, applied to a space of a given dimension.
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
        //  Reference:
        //
        //    Robert Banks,
        //    Slicing Pizzas, Racing Turtles, and Further Adventures in 
        //    Applied Mathematics,
        //    Princeton, 1999,
        //    ISBN13: 9780691059471,
        //    LC: QA93.B358.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int SLICE_NUM, the number of slices.
        //
        //    Input, int SLICES, the maximum number of pieces that can
        //    be created by the given number of slices applied in the given dimension.
        //
    {
        int j;

        int piece_num = 0;
        for ( j = 0; j <= Math.Min ( dim_num, slice_num ); j++ )
        {
            piece_num += typeMethods.i4_choose ( slice_num, j );
        }

        return piece_num;
    }
}