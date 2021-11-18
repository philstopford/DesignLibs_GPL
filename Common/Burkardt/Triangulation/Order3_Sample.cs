using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.TriangulationNS;

public static class Sample
{
    public static void triangulation_order3_sample(int node_num, double[] node_xy,
            int triangle_num, int[] triangle_node, int num_ran, ref int seed,
            ref double[] xd, ref int[] td)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGULATION_ORDER3_SAMPLE returns random points in a triangulation.
        //
        //  Discussion:
        //
        //    It is assumed that the triangulation consists of a set of non-overlapping
        //    triangles.
        //
        //    The point is chosen uniformly in the area covered by the triangulation.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 June 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
        //
        //    Input, int TRIANGLE_NUM, the number of triangles.
        //
        //    Input, int TRIANGLE_NODE[3*TRIANGLE_NUM], the nodes that make up the
        //    triangles.
        //
        //    Input, int NUM_RAN, the number of points to sample.
        //
        //    Input/output, int *SEED, a seed for the random number generator.
        //
        //    Output, double XD[2*NUM_RAN], the sample points.
        //
        //    Output, int TD[NUM_RAN], the triangle to which each sample point
        //    belongs.
        //
    {
        int i;
        int i1;
        int i2;
        int i3;
        int left = 0;
        int right = 0;
        double[] t = new double[2 * 3];
        //
        //  Compute the areas of the triangles.
        //  Build a cumulative area vector.
        //  Convert it to a relative cumulative area vector.
        //
        double[] area_cum = new double[triangle_num + 1];
        area_cum[0] = 0.0;

        for (i = 0; i < triangle_num; i++)
        {
            i1 = triangle_node[0 + i * 3];
            t[0 + 0 * 2] = node_xy[0 + i1 * 2];
            t[1 + 0 * 2] = node_xy[1 + i1 * 2];

            i2 = triangle_node[1 + i * 3];
            t[0 + 1 * 2] = node_xy[0 + i2 * 2];
            t[1 + 1 * 2] = node_xy[1 + i2 * 2];

            i3 = triangle_node[2 + i * 3];
            t[0 + 2 * 2] = node_xy[0 + i3 * 2];
            t[1 + 2 * 2] = node_xy[1 + i3 * 2];

            area_cum[i + 1] = area_cum[i] + typeMethods.triangle_area_2d(t);
        }

        double area_total = area_cum[triangle_num];

        for (i = 0; i <= triangle_num; i++)
        {
            area_cum[i] /= area_total;
        }

        //
        //  Pick random values.  A random value R indicates the corresponding triangle
        //  whose cumulative relative area contains R.
        //
        //  Bracket the random value in the cumulative relative areas,
        //  indicating a triangle.
        //
        //  Pick a random point in the triangle.
        //
        for (i = 0; i < num_ran; i++)
        {
            double r = UniformRNG.r8_uniform_01(ref seed);

            typeMethods.r8vec_bracket(triangle_num + 1, area_cum, r, ref left, ref right);

            td[i] = right - 1;

            i1 = triangle_node[0 + (td[i] - 1) * 3];
            t[0 + 0 * 2] = node_xy[0 + i1 * 2];
            t[1 + 0 * 2] = node_xy[1 + i1 * 2];

            i2 = triangle_node[1 + (td[i] - 1) * 3];
            t[0 + 1 * 2] = node_xy[0 + i2 * 2];
            t[1 + 1 * 2] = node_xy[1 + i2 * 2];

            i3 = triangle_node[2 + (td[i] - 1) * 3];
            t[0 + 2 * 2] = node_xy[0 + i3 * 2];
            t[1 + 2 * 2] = node_xy[1 + i3 * 2];

            typeMethods.triangle_sample(t, 1, ref seed, ref xd, pIndex: + i * 2);
        }
    }
}