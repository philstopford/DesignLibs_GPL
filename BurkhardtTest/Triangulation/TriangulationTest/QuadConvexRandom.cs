using Burkardt.PointsNS;
using Burkardt.Uniform;

namespace TriangulationTest;

public static partial class TriangulationSampleData
{
    public static void quad_convex_random ( ref int seed, ref double[] xy )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    QUAD_CONVEX_RANDOM returns a random convex quadrilateral.
        //
        //  Description:
        //
        //    The quadrilateral is constrained in that the vertices must all lie
        //    with the unit square.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input/output, int *SEED, a seed for the random number
        //    generator.
        //
        //    Output, double XY[2*NODE_NUM], the coordinates of the
        //    nodes of the quadrilateral, given in counterclockwise order.
        //
    {
        int[] hull = new int[4];
        int hull_num = 0;
        int i;
        int j;
        double[] xy_random = new double[2*4];

        for ( ; ; )
        {
            //
            //  Generate 4 random points.
            //
            xy_random = UniformRNG.r8mat_uniform_01 ( 2, 4, ref seed );
            //
            //  Determine the convex hull.
            //
            Points.points_hull_2d ( 4, xy_random, ref hull_num, ref hull );
            //
            //  If HULL_NUM < 4, then our convex hull is a triangle.
            //  Try again.
            //
            if ( hull_num == 4 )
            {
                break;
            }
        }
        //
        //  Make an ordered copy of the random points.
        //
        for ( j = 0; j < 4; j++ )
        {
            for ( i = 0; i < 2; i++ )
            {
                xy[i+j*2] = xy_random[i+(hull[j]-1)*2];
            }
        }
    }
}