using System;
using Burkardt.TriangleNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.QuadMesh
{
    public static class Sample
    {
        public static void sample_q4_mesh(int node_num, double[] node_xy, int element_num,
                int[] element_node, int sample_num, ref int seed, ref double[] sample_xy,
                ref int[] sample_element)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SAMPLE_Q4_MESH returns random points in a Q4 mesh.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    15 March 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NODE_NUM, the number of nodes.
            //
            //    Input, double NODE_XY(2,NODE_NUM), the coordinates of the nodes.
            //
            //    Input, int ELEMENT_NUM, the number of elements.
            //
            //    Input, int ELEMENT_NODE(4,ELEMENT_NUM), the nodes
            //    that form the elements.
            //
            //    Input, int SAMPLE_NUM, the number of points to sample.
            //
            //    Input/output, int *SEED, a seed for the random
            //     number generator.
            //
            //    Output, double SAMPLE_XY(2,SAMPLE_NUM), the sample points.
            //
            //    Output, int SAMPLE_ELEMENT(SAMPLE_NUM), the elements from
            //    which each point was drawn.
            //
        {
            double area;
            double[] area_cum;
            double area_total;
            int element;
            int i1;
            int i2;
            int i3;
            int i4;
            int left = 0;
            double[] quad_xy = new double[2 * 4];
            double r;
            int right = 0;
            int sample;
            //
            //  Compute the areas of the quadrilaterals.
            //
            area_cum = new double[element_num + 1];

            area_cum[0] = 0.0;

            for (element = 1; element <= element_num; element++)
            {
                i1 = element_node[0 + (element - 1) * 4];
                i2 = element_node[1 + (element - 1) * 4];
                i3 = element_node[2 + (element - 1) * 4];
                i4 = element_node[3 + (element - 1) * 4];

                quad_xy[0 + 0 * 2] = node_xy[0 + i1 * 2];
                quad_xy[1 + 0 * 2] = node_xy[1 + i1 * 2];
                quad_xy[0 + 1 * 2] = node_xy[0 + i2 * 2];
                quad_xy[1 + 1 * 2] = node_xy[1 + i2 * 2];
                quad_xy[0 + 2 * 2] = node_xy[0 + i3 * 2];
                quad_xy[1 + 2 * 2] = node_xy[1 + i3 * 2];
                quad_xy[0 + 3 * 2] = node_xy[0 + i4 * 2];
                quad_xy[1 + 3 * 2] = node_xy[1 + i4 * 2];

                area = Area.area_quad(quad_xy);

                area_cum[element] = area_cum[element - 1] + area;
            }

            area_total = area_cum[element_num];

            for (element = 0; element <= element_num; element++)
            {
                area_cum[element] = area_cum[element] / area_total;
            }

            //
            //  A random value R indicates the corresponding quadrilateral whose
            //  cumulative relative area first includes the number R.
            //
            for (sample = 0; sample < sample_num; sample++)
            {
                r = UniformRNG.r8_uniform_01(ref seed);

                typeMethods.r8vec_bracket(element_num + 1, area_cum, r, ref left, ref right);

                element = right - 1;

                i1 = element_node[0 + (element - 1) * 4];
                i2 = element_node[1 + (element - 1) * 4];
                i3 = element_node[2 + (element - 1) * 4];
                i4 = element_node[3 + (element - 1) * 4];

                quad_xy[0 + 0 * 2] = node_xy[0 + i1 * 2];
                quad_xy[1 + 0 * 2] = node_xy[1 + i1 * 2];
                quad_xy[0 + 1 * 2] = node_xy[0 + i2 * 2];
                quad_xy[1 + 1 * 2] = node_xy[1 + i2 * 2];
                quad_xy[0 + 2 * 2] = node_xy[0 + i3 * 2];
                quad_xy[1 + 2 * 2] = node_xy[1 + i3 * 2];
                quad_xy[0 + 3 * 2] = node_xy[0 + i4 * 2];
                quad_xy[1 + 3 * 2] = node_xy[1 + i4 * 2];

                sample_quad(quad_xy, 1, ref seed, ref sample_xy, xyIndex: +sample * 2);

                sample_element[sample] = element;
            }
        }

        public static void sample_quad(double[] quad_xy, int n, ref int seed, ref double[] xy, int xyIndex = 0)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SAMPLE_QUAD returns random points in a quadrilateral.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    22 February 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double QUAD_XY[2*4], the coordinates of the nodes.
            //
            //    Input, int N, the number of points to sample.
            //
            //    Input/output, int *SEED, a seed for the random 
            //     number generator.
            //
            //    Output, double XY[2*N], the sample points.
            //
        {
            double area1;
            double area2;
            double area_total;
            int i;
            double r;
            double[] t1 = new double[2 * 3];
            double[] t2 = new double[2 * 3];

            t1[0 + 0 * 2] = quad_xy[0 + 0 * 2];
            t1[1 + 0 * 2] = quad_xy[1 + 0 * 2];
            t1[0 + 1 * 2] = quad_xy[0 + 1 * 2];
            t1[1 + 1 * 2] = quad_xy[1 + 1 * 2];
            t1[0 + 2 * 2] = quad_xy[0 + 2 * 2];
            t1[1 + 2 * 2] = quad_xy[1 + 2 * 2];

            area1 = Integrals.triangle_area(t1);

            t2[0 + 0 * 2] = quad_xy[0 + 2 * 2];
            t2[1 + 0 * 2] = quad_xy[1 + 2 * 2];
            t2[0 + 1 * 2] = quad_xy[0 + 3 * 2];
            t2[1 + 1 * 2] = quad_xy[1 + 3 * 2];
            t2[0 + 2 * 2] = quad_xy[0 + 0 * 2];
            t2[1 + 2 * 2] = quad_xy[1 + 0 * 2];

            area2 = Integrals.triangle_area(t2);

            if (area1 < 0.0 || area2 < 0.0)
            {
                t1[0 + 0 * 2] = quad_xy[0 + 1 * 2];
                t1[1 + 0 * 2] = quad_xy[1 + 1 * 2];
                t1[0 + 1 * 2] = quad_xy[0 + 2 * 2];
                t1[1 + 1 * 2] = quad_xy[1 + 2 * 2];
                t1[0 + 2 * 2] = quad_xy[0 + 3 * 2];
                t1[1 + 2 * 2] = quad_xy[1 + 3 * 2];

                area1 = Integrals.triangle_area(t1);

                t2[0 + 0 * 2] = quad_xy[0 + 3 * 2];
                t2[1 + 0 * 2] = quad_xy[1 + 3 * 2];
                t2[0 + 1 * 2] = quad_xy[0 + 0 * 2];
                t2[1 + 1 * 2] = quad_xy[1 + 0 * 2];
                t2[0 + 2 * 2] = quad_xy[0 + 1 * 2];
                t2[1 + 2 * 2] = quad_xy[1 + 1 * 2];

                area2 = Integrals.triangle_area(t2);

                if (area1 < 0.0 || area2 < 0.0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("SAMPLE_QUAD - Fatal error!");
                    Console.WriteLine("  The quadrilateral nodes seem to be listed in");
                    Console.WriteLine("  the wrong order, or the quadrilateral is");
                    Console.WriteLine("  degenerate.");
                    return;
                }
            }

            area_total = area1 + area2;
            //
            //  Choose a triangle at random, weighted by the areas.
            //  Then choose a point in that triangle.
            //
            for (i = 0; i < n; i++)
            {
                r = UniformRNG.r8_uniform_01(ref seed);

                if (r * area_total < area1)
                {
                    typeMethods.triangle_sample(t1, 1, ref seed, ref xy, pIndex: +i * 2);
                }
                else
                {
                    typeMethods.triangle_sample(t2, 1, ref seed, ref xy, pIndex: +i * 2);
                }
            }
        }

        public static double[] sample_quad_new(double[] quad_xy, int n, ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SAMPLE_QUAD_NEW returns random points in a quadrilateral.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    22 February 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double QUAD_XY[2*4], the coordinates of the nodes.
            //
            //    Input, int N, the number of points to sample.
            //
            //    Input/output, int *SEED, a seed for the random 
            //     number generator.
            //
            //    Output, double SAMPLE_QUAD[2*N], the sample points.
            //
        {
            double area1;
            double area2;
            double area_total;
            int i;
            double r;
            double[] t1 = new double[2 * 3];
            double[] t2 = new double[2 * 3];
            double[] xy;

            t1[0 + 0 * 2] = quad_xy[0 + 0 * 2];
            t1[1 + 0 * 2] = quad_xy[1 + 0 * 2];
            t1[0 + 1 * 2] = quad_xy[0 + 1 * 2];
            t1[1 + 1 * 2] = quad_xy[1 + 1 * 2];
            t1[0 + 2 * 2] = quad_xy[0 + 2 * 2];
            t1[1 + 2 * 2] = quad_xy[1 + 2 * 2];

            area1 = Integrals.triangle_area(t1);

            t2[0 + 0 * 2] = quad_xy[0 + 2 * 2];
            t2[1 + 0 * 2] = quad_xy[1 + 2 * 2];
            t2[0 + 1 * 2] = quad_xy[0 + 3 * 2];
            t2[1 + 1 * 2] = quad_xy[1 + 3 * 2];
            t2[0 + 2 * 2] = quad_xy[0 + 0 * 2];
            t2[1 + 2 * 2] = quad_xy[1 + 0 * 2];

            area2 = Integrals.triangle_area(t2);

            if (area1 < 0.0 || area2 < 0.0)
            {
                t1[0 + 0 * 2] = quad_xy[0 + 1 * 2];
                t1[1 + 0 * 2] = quad_xy[1 + 1 * 2];
                t1[0 + 1 * 2] = quad_xy[0 + 2 * 2];
                t1[1 + 1 * 2] = quad_xy[1 + 2 * 2];
                t1[0 + 2 * 2] = quad_xy[0 + 3 * 2];
                t1[1 + 2 * 2] = quad_xy[1 + 3 * 2];

                area1 = Integrals.triangle_area(t1);

                t2[0 + 0 * 2] = quad_xy[0 + 3 * 2];
                t2[1 + 0 * 2] = quad_xy[1 + 3 * 2];
                t2[0 + 1 * 2] = quad_xy[0 + 0 * 2];
                t2[1 + 1 * 2] = quad_xy[1 + 0 * 2];
                t2[0 + 2 * 2] = quad_xy[0 + 1 * 2];
                t2[1 + 2 * 2] = quad_xy[1 + 1 * 2];

                area2 = Integrals.triangle_area(t2);

                if (area1 < 0.0 || area2 < 0.0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("SAMPLE_QUAD - Fatal error!");
                    Console.WriteLine("  The quadrilateral nodes seem to be listed in");
                    Console.WriteLine("  the wrong order, or the quadrilateral is");
                    Console.WriteLine("  degenerate.");
                    return null;
                }
            }

            area_total = area1 + area2;
            //
            //  Choose a triangle at random, weighted by the areas.
            //  Then choose a point in that triangle.
            //
            xy = new double[2 * n];

            for (i = 0; i < n; i++)
            {
                r = UniformRNG.r8_uniform_01(ref seed);

                if (r * area_total < area1)
                {
                    typeMethods.triangle_sample(t1, 1, ref seed, ref xy, +i * 2);
                }
                else
                {
                    typeMethods.triangle_sample(t2, 1, ref seed, ref xy, +i * 2);
                }
            }

            return xy;
        }
    }
}