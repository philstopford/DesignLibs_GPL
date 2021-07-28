namespace Burkardt.SimplexNS
{
    public static class Simplex
    {
        public static double[] simplex_to_triangle ( double[] tvert1, double[] tvert2, 
        double[] tvert3, double[] s )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SIMPLEX_TO_TRIANGLE maps points from the simplex to a triangle.
        //
        //  Discussion:
        //
        //    The simplex has vertices:
        //
        //      (  0, 0 )
        //      (  1, 0 )
        //      (  0, 1 )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double TVERT1[2], TVERT2[2], TVERT3[2], the coordinates
        //    of the vertices of the triangle.  These vertices will be taken
        //    to be the images of (0,0), (1,0) and (0,1) respectively.
        //
        //    Input, double S[2], the coordinates of the point in the simplex.
        //
        //    Output, double SIMPLEX_TO_TRIANGLE[2], the coordinates of the point in
        //    the triangle.
        //
        {
            int i;
            double[] t;

            t = new double[2];

            for ( i = 0; i < 2; i++ )
            {
                t[i] = tvert1[i] * ( 1.0 - s[0] - s[1] ) 
                       + tvert2[i] * s[0]
                       + tvert3[i] * s[1];
            }

            return t;
        }
    }
}