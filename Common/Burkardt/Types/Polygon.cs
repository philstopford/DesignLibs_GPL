namespace Burkardt.Types
{
    public static class Polygon
    {
        public static double[] polygon_centroid_2d(int n, double[] v)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYGON_CENTROID_2D computes the centroid of a polygon in 2D.
            //
            //  Formula:
            //
            //    Denoting the centroid coordinates by CENTROID, then
            //
            //      CENTROID(1) = Integral ( Polygon interior ) x dx dy / Area ( Polygon )
            //      CENTROID(2) = Integral ( Polygon interior ) y dx dy / Area ( Polygon ).
            //
            //    Green's theorem states that
            //
            //      Integral ( Polygon boundary ) ( M dx + N dy ) =
            //      Integral ( Polygon interior ) ( dN/dx - dM/dy ) dx dy.
            //
            //    Using M = 0 and N = x * x / 2, we get:
            //
            //      CENTROID(1) = 0.5 * Integral ( Polygon boundary ) x * x dy,
            //
            //    which becomes
            //
            //      CENTROID(1) = 1/6 Sum ( 1 <= I <= N )
            //        ( X(I+1) + X(I) ) * ( X(I) * Y(I+1) - X(I+1) * Y(I))
            //
            //    where, when I = N, the index "I+1" is replaced by 1.
            //
            //    A similar calculation gives us a formula for CENTROID(2).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    17 August 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Gerard Bashein and Paul Detmer,
            //    Centroid of a Polygon,
            //    Graphics Gems IV, edited by Paul Heckbert,
            //    AP Professional, 1994.
            //
            //  Parameters:
            //
            //    Input, int N, the number of sides of the polygonal shape.
            //
            //    Input, double V[2*N], the coordinates of the vertices
            //    of the shape.
            //
            //    Output, double POLYGON_CENTROID_2D[2], the coordinates of the
            //    centroid of the shape.
            //
        {
            double area;
            double[] centroid;
            int i;
            int ip1;
            double temp;
            //
            area = 0.0;
            centroid = new double[2];
            centroid[0] = 0.0;
            centroid[1] = 0.0;

            for (i = 0; i < n; i++)
            {
                if (i < n - 1)
                {
                    ip1 = i + 1;
                }
                else
                {
                    ip1 = 0;
                }

                temp = (v[0 + i * 2] * v[1 + ip1 * 2] - v[0 + ip1 * 2] * v[1 + i * 2]);

                area = area + temp;

                centroid[0] = centroid[0] + (v[0 + ip1 * 2] + v[0 + i * 2]) * temp;
                centroid[1] = centroid[1] + (v[1 + ip1 * 2] + v[1 + i * 2]) * temp;

            }

            area = area / 2.0;

            centroid[0] = centroid[0] / (6.0 * area);
            centroid[1] = centroid[1] / (6.0 * area);

            return centroid;
        }
    }
}