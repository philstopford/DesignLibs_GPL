namespace Burkardt.FEM
{
    public static class XY
    {
        public static void xy_set ( int nx, int ny, int node_num, double xl, double xr, double yb,
                double yt, ref double[] node_xy )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    XY_SET sets the XY coordinates of the nodes.
            //
            //  Discussion:
            //
            //    The nodes are laid out in an evenly spaced grid, in the unit square.
            //
            //    The first node is at the origin.  More nodes are created to the
            //    right until the value of X = 1 is reached, at which point
            //    the next layer is generated starting back at X = 0, and an
            //    increased value of Y.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    19 March 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NX, NY, the number of elements in the X and
            //    Y direction.
            //
            //    Input, int NODE_NUM, the number of nodes.
            //
            //    Input, double XL, XR, YB, YT, the X coordinates of
            //    the left and right sides of the rectangle, and the Y coordinates
            //    of the bottom and top of the rectangle.
            //
            //    Output, double NODE_XY[2*NODE_NUM], the nodes.
            //
        {
            int i;
            int j;

            for ( j = 1; j <= 2*ny-1; j++ )
            {
                for ( i = 1; i <= 2*nx - 1; i++ )
                {
                    node_xy[0+(i-1+(j-1)*(2*nx-1))*2] =
                        ( ( double ) ( 2 * nx - i - 1 ) * xl
                          + ( double ) (          i - 1 ) * xr )
                        / ( double ) ( 2 * nx     - 2 );

                    node_xy[1+(i-1+(j-1)*(2*nx-1))*2] =
                        ( ( double ) ( 2 * ny - j - 1 ) * yb
                          + ( double ) (          j - 1 ) * yt )
                        / ( double ) ( 2 * ny     - 2 );

                }
            }
        }
    }
}