using System;

namespace Burkardt.Polygon
{
    public static class Triangulate
    {
                public static int[] polygon_triangulate(int n, double[] x, double[] y)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYGON_TRIANGULATE determines a triangulation of a polygon.
            //
            //  Discussion:
            //
            //    There are N-3 triangles in the triangulation.
            //
            //    For the first N-2 triangles, the first edge listed is always an
            //    internal diagonal.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 May 2014
            //
            //  Author:
            //
            //    Original C version by Joseph ORourke.
            //    This C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Joseph ORourke,
            //    Computational Geometry in C,
            //    Cambridge, 1998,
            //    ISBN: 0521649765,
            //    LC: QA448.D38.
            //
            //  Parameters:
            //
            //    Input, int N, the number of vertices.
            //
            //    Input, double X[N], Y[N], the coordinates of each vertex.
            //
            //    Output, int TRIANGLES[3*(N-2)], the triangles of the 
            //    triangulation.
            //
        {
            double area;
            bool[] ear;
            int i;
            int i0;
            int i1;
            int i2;
            int i3;
            int i4;
            int[] next;
            int node;
            int node_m1;
            int[] prev;
            int triangle_num;
            int[] triangles;
            //
            //  We must have at least 3 vertices.
            //
            if (n < 3)
            {
                Console.WriteLine("");
                Console.WriteLine("POLYGON_TRIANGULATE - Fatal error!");
                Console.WriteLine("  N < 3.");
                return null;
            }

            //
            //  Consecutive vertices cannot be equal.
            //
            node_m1 = n - 1;
            for (node = 0; node < n; node++)
            {
                if (x[node_m1] == x[node] && y[node_m1] == y[node])
                {
                    Console.WriteLine("");
                    Console.WriteLine("POLYGON_TRIANGULATE - Fatal error!");
                    Console.WriteLine("  Two consecutive nodes are identical.");
                    return null;
                }

                node_m1 = node;
            }

            //
            //  Area must be positive.
            //
            area = 0.0;
            for (node = 0; node < n - 2; node++)
            {
                area = area + 0.5 *
                (
                    (x[node + 1] - x[node]) * (y[node + 2] - y[node])
                    - (x[node + 2] - x[node]) * (y[node + 1] - y[node])
                );
            }

            if (area <= 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("POLYGON_TRIANGULATE - Fatal error!");
                Console.WriteLine("  Polygon has zero or negative area.");
                return null;
            }

            triangles = new int[3 * (n - 2)];
            //
            //  PREV and NEXT point to the previous and next nodes.
            //
            prev = new int[n];
            next = new int[n];

            i = 0;
            prev[i] = n - 1;
            next[i] = i + 1;

            for (i = 1; i < n - 1; i++)
            {
                prev[i] = i - 1;
                next[i] = i + 1;
            }

            i = n - 1;
            prev[i] = i - 1;
            next[i] = 0;
            //
            //  EAR indicates whether the node and its immediate neighbors form an ear
            //  that can be sliced off immediately.
            //
            ear = new bool[n];
            for (i = 0; i < n; i++)
            {
                ear[i] = Properties.diagonal(prev[i], next[i], n, prev, next, x, y);
            }

            triangle_num = 0;

            i2 = 0;

            while (triangle_num < n - 3)
            {
                //
                //  If I2 is an ear, gather information necessary to carry out
                //  the slicing operation and subsequent "healing".
                //
                if (ear[i2])
                {
                    i3 = next[i2];
                    i4 = next[i3];
                    i1 = prev[i2];
                    i0 = prev[i1];
                    //
                    //  Make vertex I2 disappear.
                    //
                    next[i1] = i3;
                    prev[i3] = i1;
                    //
                    //  Update the earity of I1 and I3, because I2 disappeared.
                    //
                    ear[i1] = Properties.diagonal(i0, i3, n, prev, next, x, y);
                    ear[i3] = Properties.diagonal(i1, i4, n, prev, next, x, y);
                    //
                    //  Add the diagonal [I3, I1, I2] to the list.
                    //
                    triangles[0 + triangle_num * 3] = i3;
                    triangles[1 + triangle_num * 3] = i1;
                    triangles[2 + triangle_num * 3] = i2;
                    triangle_num = triangle_num + 1;
                }

                //
                //  Try the next vertex.
                //
                i2 = next[i2];
            }

            //
            //  The last triangle is formed from the three remaining vertices.
            //
            i3 = next[i2];
            i1 = prev[i2];

            triangles[0 + triangle_num * 3] = i3;
            triangles[1 + triangle_num * 3] = i1;
            triangles[2 + triangle_num * 3] = i2;
            triangle_num = triangle_num + 1;

            return triangles;
        }

    }
}