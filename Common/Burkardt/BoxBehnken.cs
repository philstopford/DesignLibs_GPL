using System;

namespace Burkardt.BoxBehnkenNS;

public static class BoxBehnken
{
    public static double[] box_behnken ( int dim_num, int x_num, double[] range )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BOX_BEHNKEN returns a Box-Behnken design for the given number of factors.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 October 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    George Box, Donald Behnken,
        //    Some new three level designs for the study of quantitative variables,
        //    Technometrics,
        //    Volume 2, pages 455-475, 1960.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int X_NUM, the number of elements of the design.
        //    X_NUM should be equal to DIM_NUM * 2^(DIM_NUM-1) + 1.
        //
        //    Input, double RANGE[DIM_NUM*2], the minimum and maximum
        //    value for each component.
        //
        //    Output, double BOX_BEHNKEN[DIM_NUM*X_NUM], the elements of the design.
        //
    {
        //
        //  Ensure that the range is legal.
        //
        for (int i = 0; i < dim_num; i++ )
        {
            if (!(range[i + 1 * dim_num] <= range[i + 0 * dim_num]))
            {
                continue;
            }

            Console.WriteLine();
            Console.WriteLine("BOX_BEHNKEN - Fatal error!");
            Console.WriteLine("  RANGE[" + i + ",1] <= RANGE[" + i + ",0].");
            return null;
        }

        double[] x = new double[dim_num*x_num];
        //
        //  The first point is the center.
        //
        int j = 0;

        for (int i = 0; i < dim_num; i++ )
        {
            x[i+j*dim_num] = ( range[i+0*dim_num] + range[i+1*dim_num] ) / 2.0;
        }
        //
        //  For subsequent elements, one entry is fixed at the middle of the range.
        //  The others are set to either extreme.
        //
        for (int i = 0; i < dim_num; i++ )
        {
            j += 1;
            for (int i2 = 0; i2 < dim_num; i2++ )
            {
                x[i2+j*dim_num] = range[i2+0*dim_num];
            }
            x[i+j*dim_num] = ( range[i+0*dim_num] + range[i+1*dim_num] ) / 2.0;
            //
            //  The next element is made by finding the last low value, making it
            //  high, and all subsequent high values low.
            //
            for ( ; ; )
            {
                int last_low = -1;

                for (int i2 = 0; i2 < dim_num; i2++ )
                {
                    if ( Math.Abs(x[i2+j*dim_num] - range[i2+0*dim_num]) <= double.Epsilon )
                    {
                        last_low = i2;
                    }
                }

                if ( last_low == -1 )
                {
                    break;
                }

                j += 1;
                for (int i2 = 0; i2 < dim_num; i2++ )
                {
                    x[i2+j*dim_num] = x[i2+(j-1)*dim_num];
                }
                x[last_low+j*dim_num] = range[last_low+1*dim_num];

                for (int i2 = last_low + 1; i2 < dim_num; i2++ )
                {
                    if ( Math.Abs(x[i2+j*dim_num] - range[i2+1*dim_num]) <= double.Epsilon )
                    {
                        x[i2+j*dim_num] = range[i2+0*dim_num];
                    }
                }
            }
        }
        return x;
    }

        
    public static int box_behnken_size ( int dim_num )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BOX_BEHNKEN_SIZE returns the size of a Box-Behnken design.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 October 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    George Box, Donald Behnken,
        //    Some new three level designs for the study of quantitative variables,
        //    Technometrics,
        //    Volume 2, pages 455-475, 1960.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Output, int X_NUM, the number of elements of the design.
        //    X_NUM will be equal to DIM_NUM * 2^(DIM_NUM-1) + 1.
        //
    {
        int x_num = dim_num switch
        {
            >= 1 => 1 + dim_num * (int) Math.Pow(2, dim_num - 1),
            _ => -1
        };

        return x_num;
    }
        
}