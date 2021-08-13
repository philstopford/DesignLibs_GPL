using System;

namespace Burkardt.Function
{
    public static class Hilbert
    {
        public static void d2xy ( int m, int d, ref int x, ref int y )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    D2XY converts a 1D Hilbert coordinate to a 2D Cartesian coordinate.
        //
        //  Modified:
        //
        //    24 December 2015
        //
        //  Parameters:
        //
        //    Input, int M, the index of the Hilbert curve.
        //    The number of cells is N=2^M.
        //    0 < M.
        //
        //    Input, int D, the Hilbert coordinate of the cell.
        //    0 <= D < N * N.
        //
        //    Output, int &X, &Y, the Cartesian coordinates of the cell.
        //    0 <= X, Y < N.
        //
        {
            int n;
            int rx;
            int ry;
            int s;
            int t = d;

            n = (int)Math.Pow ( 2, m );

            x = 0;
            y = 0;
            for ( s = 1; s < n; s = s * 2 )
            {
                rx = ( 1 & ( t / 2 ) );
                ry = ( 1 & ( t ^ rx ) );
                rot ( s, ref x, ref y, rx, ry );
                x = x + s * rx;
                y = y + s * ry;
                t = t / 4;
            }
        }
        
        public static void rot ( int n, ref int x, ref int y, int rx, int ry ) 

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ROT rotates and flips a quadrant appropriately
        //
        //  Modified:
        //
        //    24 December 2015
        //
        //  Parameters:
        //
        //    Input, int N, the length of a side of the square.  N must be a power of 2.
        //
        //    Input/output, int &X, &Y, the input and output coordinates of a point.
        //
        //    Input, int RX, RY, ???
        //
        {
            int t;

            if ( ry == 0 )
            {
                //
                //  Reflect.
                //
                if ( rx == 1 )
                {
                    x = n - 1 - x;
                    y = n - 1 - y;
                }
                //
                //  Flip.
                //
                t = x;
                x = y;
                y = t;
            }
        }
        
        public static int xy2d ( int m, int x, int y )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    XY2D converts a 2D Cartesian coordinate to a 1D Hilbert coordinate.
            //
            //  Discussion:
            //
            //    It is assumed that a square has been divided into an NxN array of cells,
            //    where N is a power of 2.
            //
            //    Cell (0,0) is in the lower left corner, and (N-1,N-1) in the upper 
            //    right corner.
            //
            //  Modified:
            //
            //    24 December 2015
            //
            //  Parameters:
            //
            //    Input, int M, the index of the Hilbert curve.
            //    The number of cells is N=2^M.
            //    0 < M.
            //
            //    Input, int X, Y, the Cartesian coordinates of a cell.
            //    0 <= X, Y < N.
            //
            //    Output, int XY2D, the Hilbert coordinate of the cell.
            //    0 <= D < N * N.
            //
        {
            int d = 0;
            int n;
            int rx;
            int ry;
            int s;

            n = (int)Math.Pow ( 2, m );

            for ( s = n / 2; s > 0; s = s / 2 )
            {
                rx = ( x & s ) > 0 ? 1 : 0;
                ry = ( y & s ) > 0 ? 1 : 0;
                d = d + s * s * ( ( 3 * rx ) ^ ry );
                rot ( s, ref x, ref y, rx, ry );
            }
            return d;
        }
    }
}