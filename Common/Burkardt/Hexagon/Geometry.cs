using System;
using Burkardt.Types;

namespace Burkardt.Hexagon
{
    public static class Geometry
    {
        public static double hexagon_area_2d(double r)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HEXAGON_AREA_2D returns the area of a regular hexagon in 2D.
            //
            //  Discussion:
            //
            //    The radius of a regular hexagon is the distance from the center
            //    of the hexagon to any vertex.  This happens also to equal the
            //    length of any side.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    01 September 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double R, the radius of the hexagon.
            //
            //    Output, double HEXAGON_AREA_2D, the area of the hexagon.
            //
        {
            double value;

            value = r * r * hexagon_unit_area_2d();

            return value;
        }

        public static bool hexagon_contains_point_2d(double[] v, double[] p)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HEXAGON_CONTAINS_POINT_2D finds if a point is inside a hexagon in 2D.
            //
            //  Discussion:
            //
            //    This test is only valid if the hexagon is convex.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    20 June 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double V[2*6], the vertics, in counter clockwise order.
            //
            //    Input, double P[2], the point to be tested.
            //
            //    Output, bool HEXAGON_CONTAINS_POINT_2D, is TRUE if X is in the hexagon.
            //
        {
            int N = 6;

            int i;
            int j;
            //
            //  A point is inside a convex hexagon if and only if it is "inside"
            //  each of the 6 halfplanes defined by lines through consecutive
            //  vertices.
            //
            for (i = 0; i < N; i++)
            {
                j = (i + 1) % N;

                if (v[0 + i * 2] * (v[1 + j * 2] - p[1])
                    + v[0 + j * 2] * (p[1] - v[1 + i * 2])
                    + p[0] * (v[1 + i * 2] - v[1 + j * 2]) < 0.0)
                {
                    return false;
                }
            }

            return true;
        }

        public static void hexagon_shape_2d(double angle, ref double[] p)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HEXAGON_SHAPE_2D returns points on the unit regular hexagon in 2D.
            //
            //  Discussion:
            //
            //    The unit regular hexagon has radius 1. The radius is the distance from
            //    the center to any vertex, and it is also the length of any side.
            //    An example of a unit hexagon is the convex hull of the points
            //
            //      (   1,              0 ),
            //      (   0.5,   Math.Sqrt (3)/2 ),
            //      ( - 0.5,   Math.Sqrt (3)/2 ),
            //      ( - 1,              0 ),
            //      ( - 0.5, - Math.Sqrt (3)/2 ),
            //      (   0.5, - Math.Sqrt (3)/2 ).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    08 September 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double ANGLE, the angle, in degrees, of the point.
            //
            //    Output, double P[2], the coordinates of the point.
            //
        {
            //
            //  Ensure that 0.0 <= ANGLE2 < 360.
            //
            angle = typeMethods.r8_modp(angle, 360.0);
            //
            //  y = - sqrt(3) * x + sqrt(3)
            //
            if (0.0 <= angle && angle <= 60.0)
            {
                p[0] = Math.Sqrt(3.0) / (typeMethods.r8_tand(angle) + Math.Sqrt(3.0));
                p[1] = typeMethods.r8_tand(angle) * p[0];
            }
            //
            //  y = sqrt(3) / 2
            //
            else if (angle <= 120.0)
            {
                p[1] = Math.Sqrt(3.0) / 2.0;
                p[0] = typeMethods.r8_cotd(angle) * p[1];
            }
            //
            //  y = sqrt(3) * x + sqrt(3)
            //
            else if (angle <= 180.0)
            {
                p[0] = Math.Sqrt(3.0) / (typeMethods.r8_tand(angle) - Math.Sqrt(3.0));
                p[1] = typeMethods.r8_tand(angle) * p[0];
            }
            //
            //  y = - sqrt(3) * x - sqrt(3)
            //
            else if (angle <= 240.0)
            {
                p[0] = -Math.Sqrt(3.0) / (typeMethods.r8_tand(angle) + Math.Sqrt(3.0));
                p[1] = typeMethods.r8_tand(angle) * p[0];
            }
            //
            //  y = - sqrt(3) / 2
            //
            else if (angle <= 300.0)
            {
                p[1] = -Math.Sqrt(3.0) / 2.0;
                p[0] = typeMethods.r8_cotd(angle) * p[1];
            }
            //
            //  y = sqrt(3) * x - sqrt(3)
            //
            else if (angle <= 360.0)
            {
                p[0] = -Math.Sqrt(3.0) / (typeMethods.r8_tand(angle) - Math.Sqrt(3.0));
                p[1] = typeMethods.r8_tand(angle) * p[0];
            }

        }

        public static double hexagon_unit_area_2d()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HEXAGON_UNIT_AREA_2D returns the area of a unit regular hexagon in 2D.
            //
            //  Discussion:
            //
            //    A "unit" regular hexagon has both a "radius" of 1 (distance
            //    from the center to any vertex), and a side length of 1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    01 September 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, double HEXAGON_UNIT_AREA_2D, the area of the hexagon.
            //
        {
            double value;

            value = 3.0 * Math.Sqrt(3.0) / 2.0;

            return value;
        }

        public static void hexagon_vertices_2d(ref double[] h)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    HEXAGON_VERTICES_2D returns the vertices of the unit hexagon in 2D.
            //
            //  Discussion:
            //
            //    The unit hexagon has maximum radius 1, and is the hull of the points
            //
            //      (   1,              0 ),
            //      (   0.5,   Math.Sqrt (3)/2 ),
            //      ( - 0.5,   Math.Sqrt (3)/2 ),
            //      ( - 1,              0 ),
            //      ( - 0.5, - Math.Sqrt (3)/2 ),
            //      (   0.5, - Math.Sqrt (3)/2 ).
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
            //    Output, double H[2*6], the coordinates of the vertices.
            //
        {
            double A = 0.8660254037844386;
            int DIM_NUM = 2;

            h[0 + 0 * 2] = 1.0;
            h[0 + 1 * 2] = 0.5;
            h[0 + 2 * 2] = -0.5;
            h[0 + 3 * 2] = -1.0;
            h[0 + 4 * 2] = -0.5;
            h[0 + 5 * 2] = 0.5;

            h[1 + 0 * 2] = 0.0;
            h[1 + 1 * 2] = A;
            h[1 + 2 * 2] = A;
            h[1 + 3 * 2] = 0.0;
            h[1 + 4 * 2] = -A;
            h[1 + 5 * 2] = -A;

        }
    }
}