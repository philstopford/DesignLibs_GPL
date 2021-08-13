namespace Burkardt.SimplexNS
{
    public static class Simplex
    {
        public static int simplex_num ( int m, int n )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SIMPLEX_NUM evaluates the N-th Simplex number in M dimensions.
            //
            //  Discussion:
            //
            //     N\M: 1    2    3    4    5
            //    --   --   --   --   --   --
            //     0    0    0    0    0    0
            //     1    1    1    1    1    1
            //     2    2    3    4    5    6
            //     3    3    6   10   15   21
            //     4    4   10   20   35   56
            //     5    5   15   35   70  126
            //     6    6   21   56  126  252
            //     7    7   28   84  210  462
            //     8    8   36  120  330  792
            //     9    9   45  165  495 1287
            //    10   10   55  220  715 2002
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    27 February 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters
            //
            //    Input, int M, the spatial dimension.
            //
            //    Input, int N, the index of the number.
            //
            //    Output, int SIMPLEX_NUM, the desired value.
            //
        {
            int i;
            int value;
  
            value = 1;
            for ( i = 1; i <= m; i++ )
            {
                value = ( value * ( n + i - 1 ) ) / i;
            }

            return value;
        }

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