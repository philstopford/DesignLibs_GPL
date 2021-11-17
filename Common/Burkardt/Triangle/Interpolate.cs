namespace Burkardt.TriangleNS;

public static class Interpolate
{
    public static double[] triangle_interpolate_linear ( int m, int n, double[] p1, double[] p2, 
            double[] p3, double[] p, double[] v1, double[] v2, double[] v3 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_INTERPOLATE_LINEAR interpolates data given on a triangle's vertices.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 January 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the dimension of the quantity.
        //
        //    Input, int N, the number of points.
        //
        //    Input, double P1[2], P2[2], P3[2], the vertices of the triangle,
        //    in counterclockwise order.
        //
        //    Input, double P[2*N], the point at which the interpolant is desired.
        //
        //    Input, double V1[M], V2[M], V3[M], the value of some quantity at the vertices.
        //
        //    Output, double TRIANGLE_INTERPOLATE_LINEAR[M,N], the interpolated value 
        //    of the quantity at P.
        //
    {
        double abc;
        double apc;
        double abp;
        int i;
        int j;
        double pbc;
        double[] v;

        v = new double[m*n];

        abc = triangle_area ( p1[0], p1[1], p2[0], p2[1], p3[0], p3[1] );

        for ( j = 0; j < n; j++ )
        {
            pbc = triangle_area ( p[0+j*2], p[1+j*2], p2[0],    p2[1],    p3[0],    p3[1] );
            apc = triangle_area ( p1[0],    p1[1],    p[0+j*2], p[1+j*2], p3[0],    p3[1] );
            abp = triangle_area ( p1[0],    p1[1],    p2[0],    p2[1],    p[0+j*2], p[1+j*2] );
            for ( i = 0; i < m; i++ )
            {
                v[i+j*m] =
                    ( pbc * v1[i]
                      + apc * v2[i]
                      + abp * v3[i] )
                    / abc;
            }
        }

        return v;
    }
        
    public static double triangle_area ( double p1x, double p1y, double p2x, double p2y, 
            double p3x, double p3y )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_AREA computes the area of a triangle in 2D.
        //
        //  Discussion:
        //
        //    If the triangle's vertices are given in counter clockwise order,
        //    the area will be positive.  If the triangle's vertices are given
        //    in clockwise order, the area will be negative!
        //
        //    An earlier version of this routine always returned the absolute
        //    value of the computed area.  I am convinced now that that is
        //    a less useful result!  For instance, by returning the signed 
        //    area of a triangle, it is possible to easily compute the area 
        //    of a nonconvex polygon as the sum of the (possibly negative) 
        //    areas of triangles formed by node 1 and successive pairs of vertices.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1X, P1Y, P2X, P2Y, P3X, P3Y, the coordinates
        //    of the vertices P1, P2, and P3.
        //
        //    Output, double TRIANGLE_AREA, the area of the triangle.
        //
    {
        double area;

        area = 0.5 * ( 
            p1x * ( p2y - p3y ) + 
            p2x * ( p3y - p1y ) + 
            p3x * ( p1y - p2y ) );
 
        return area;
    }
}