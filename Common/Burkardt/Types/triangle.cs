using System;

namespace Burkardt.Types
{
    public static class triangle
    {
        public static double triangle_area_2d ( double[] t )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_AREA_2D computes the area of a triangle in 2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double T[2*3], the vertices of the triangle.
        //
        //    Output, double TRIANGLE_AREA_2D, the area of the triangle.  AREA will
        //    be nonnegative.
        //
        {
            double area = Math.Abs ( 0.5 * (
                t[0+0*2] * ( t[1+2*2] - t[1+1*2] ) +
                t[0+1*2] * ( t[1+0*2] - t[1+2*2] ) +
                t[0+2*2] * ( t[1+1*2] - t[1+0*2] ) ) );

            return area;
        }
    }
}
