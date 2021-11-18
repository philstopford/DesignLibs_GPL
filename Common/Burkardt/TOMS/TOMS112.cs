namespace Burkardt.TOMSNS;

public static partial class TOMS
{
    public static bool point_in_polygon ( int n, double[] x, double[] y, double x0, double y0 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POINT_IN_POLYGON determines if a point is inside a polygon
        //
        //  Discussion:
        //
        //    If the points ( x(i), y(i) ) ( i = 1, 2, ..., n ) are,
        //    in this cyclic order, the vertices of a simple closed polygon and
        //    (x0,y0) is a point not on any side of the polygon, then the
        //    procedure determines, by setting "point_in_polygon" to TRUE or FALSE,
        //    whether (x0,y0) lies in the interior of the polygon.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 November 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Moshe Shimrat,
        //    ACM Algorithm 112,
        //    Position of Point Relative to Polygon,
        //    Communications of the ACM,
        //    Volume 5, Number 8, page 434, August 1962.
        //
        //    Richard Hacker,
        //    Certification of Algorithm 112,
        //    Communications of the ACM,
        //    Volume 5, Number 12, page  606, December 1962.
        //
        //  Parameters:
        //
        //    Input, int N, the number of nodes or vertices in 
        //    the polygon.  N must be at least 3.
        //
        //    Input, double X[N], Y[N], the vertices of the polygon.
        //
        //    Input, double X0, Y0, the coordinates of the point to be tested.
        //
        //    Output, bool POINT_IN_POLYGON, is TRUE if the point is
        //    inside the polygon.
        //
    {
        int i;

        bool value = false;

        for ( i = 0; i < n; i++ )
        {
            int ip1 = ( i + 1 ) % n;

            if (y[ip1] < y0 != y0 <= y[i])
            {
                continue;
            }

            double t = x0 - x[i] - ( y0 - y[i] ) * ( x[ip1] - x[i] ) / ( y[ip1] - y[i] );
            value = t switch
            {
                < 0.0 => !value,
                _ => value
            };
        }

        return value;
    }
}