using Burkardt.Types;

namespace Burkardt.TriangulationNS
{
    public static class VBEDG
    {
        public static void vbedg(double x, double y, int node_num, double[] node_xy,
        int triangle_num, int[] triangle_node, int[] triangle_neighbor, ref int ltri,
        ref int ledg, ref int rtri, ref int redg )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VBEDG determines which boundary edges are visible to a point.
        //
        //  Discussion:
        //
        //    The point (X,Y) is assumed to be outside the convex hull of the
        //    region covered by the 2D triangulation.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 September 2008
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Barry Joe.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Barry Joe,
        //    GEOMPACK - a software package for the generation of meshes
        //    using geometric algorithms,
        //    Advances in Engineering Software,
        //    Volume 13, pages 325-331, 1991.
        //
        //  Parameters:
        //
        //    Input, double X, Y, the coordinates of a point outside the convex hull
        //    of the current triangulation.
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
        //
        //    Input, int TRIANGLE_NUM, the number of triangles.
        //
        //    Input, int TRIANGLE_NODE[3*TRIANGLE_NUM], the triangle incidence list.
        //
        //    Input, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], the triangle neighbor list;
        //    negative values are used for links of a counter clockwise linked list
        //    of boundary edges;
        //      LINK = -(3*I + J-1) where I, J = triangle, edge index.
        //
        //    Input/output, int *LTRI, *LEDG.  If LTRI != 0 then these values are
        //    assumed to be already computed and are not changed, else they are updated.
        //    On output, LTRI is the index of boundary triangle to the left of the
        //    leftmost boundary triangle visible from (X,Y), and LEDG is the boundary
        //    edge of triangle LTRI to the left of the leftmost boundary edge visible
        //    from (X,Y).  1 <= LEDG <= 3.
        //
        //    Input/output, int *RTRI.  On input, the index of the boundary triangle
        //    to begin the search at.  On output, the index of the rightmost boundary
        //    triangle visible from (X,Y).
        //
        //    Input/output, int *REDG, the edge of triangle RTRI that is visible
        //    from (X,Y).  1 <= REDG <= 3.
        //
        {
            int a;
            double ax;
            double ay;
            int b;
            double bx;
            double by;
            bool done;
            int e;
            int l;
            int lr;
            int t;
            //
            //  Find the rightmost visible boundary edge using links, then possibly
            //  leftmost visible boundary edge using triangle neighbor information.
            //
            if (ltri == 0)
            {
                done = false;
                ltri = rtri;
                ledg = redg;
            }
            else
            {
                done = true;
            }

            for (;;)
            {
                l = -triangle_neighbor[3 * ((rtri) - 1) + (redg) - 1];
                t = l / 3;
                e = 1 + l % 3;
                a = triangle_node[3 * (t - 1) + e - 1];

                if (e <= 2)
                {
                    b = triangle_node[3 * (t - 1) + e];
                }
                else
                {
                    b = triangle_node[3 * (t - 1) + 0];
                }

                ax = node_xy[2 * (a - 1) + 0];
                ay = node_xy[2 * (a - 1) + 1];

                bx = node_xy[2 * (b - 1) + 0];
                by = node_xy[2 * (b - 1) + 1];

                lr = Helpers.lrline(x, y, ax, ay, bx, by, 0.0);

                if (lr <= 0)
                {
                    break;
                }

                rtri = t;
                redg = e;

            }

            if (done)
            {
                return;
            }

            t = ltri;
            e = ledg;

            for (;;)
            {
                b = triangle_node[3 * (t - 1) + e - 1];
                e = typeMethods.i4_wrap(e - 1, 1, 3);

                while (0 < triangle_neighbor[3 * (t - 1) + e - 1])
                {
                    t = triangle_neighbor[3 * (t - 1) + e - 1];

                    if (triangle_node[3 * (t - 1) + 0] == b)
                    {
                        e = 3;
                    }
                    else if (triangle_node[3 * (t - 1) + 1] == b)
                    {
                        e = 1;
                    }
                    else
                    {
                        e = 2;
                    }

                }

                a = triangle_node[3 * (t - 1) + e - 1];
                ax = node_xy[2 * (a - 1) + 0];
                ay = node_xy[2 * (a - 1) + 1];

                bx = node_xy[2 * (b - 1) + 0];
                by = node_xy[2 * (b - 1) + 1];

                lr = Helpers.lrline(x, y, ax, ay, bx, by, 0.0);

                if (lr <= 0)
                {
                    break;
                }

            }

            ltri = t;
            ledg = e;
        }
    }
}