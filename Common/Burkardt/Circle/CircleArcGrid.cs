using System;

namespace Burkardt.CircleNS;

public static class Circle
{
    public static double[] circle_arc_grid ( double r, double[] c, double[] a, int n )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CIRCLE_ARC_GRID computes grid points along a circular arc.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 November 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the circle.
        //
        //    Input, double C[2], the coordinates of the center.
        //
        //    Input, double A[2], the angle of the first and last
        //    points, in DEGREES.
        //
        //    Input, int N, the number of points.
        //
        //    Output, double CIRCLE_ARC_GRID[2*N], the grid points.
        //
    {
        double aj;
        int j;
        double[] xy;

        xy = new double[2*n];

        for ( j = 0; j < n; j++ )
        {
            aj = ( (n - j - 1) * a[0]
                   + j * a[1] ) 
                 / (n     - 1);

            xy[0+j*2] = c[0] + r * Math.Cos ( aj * Math.PI / 180.0 );
            xy[1+j*2] = c[1] + r * Math.Sin ( aj * Math.PI / 180.0 );
        }

        return xy;
    }
}