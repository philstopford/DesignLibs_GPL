using System;

namespace Burkardt.TriangulationNS
{
    public class DelaunaySearchData
    {
        public int triangle_index_save = -1;
    }
    
    public static class Search
    {
        public static void triangulation_search_delaunay(ref DelaunaySearchData data, int node_num, double[] node_xy,
        int triangle_order, int triangle_num, int[] triangle_node,
        int[] triangle_neighbor, double[] p, ref int triangle_index,
        ref double alpha, ref double beta, ref double gamma, ref int edge,
        ref int step_num, int pIndex = 0 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGULATION_SEARCH_DELAUNAY searches a triangulation for a point.
        //
        //  Discussion:
        //
        //    The algorithm "walks" from one triangle to its neighboring triangle,
        //    and so on, until a triangle is found containing point P, or P is found
        //    to be outside the convex hull.
        //
        //    The algorithm computes the barycentric coordinates of the point with
        //    respect to the current triangle.  If all three quantities are positive,
        //    the point is contained in the triangle.  If the I-th coordinate is
        //    negative, then (X,Y) lies on the far side of edge I, which is opposite
        //    from vertex I.  This gives a hint as to where to search next.
        //
        //    For a Delaunay triangulation, the search is guaranteed to terminate.
        //    For other triangulations, a cycle may occur.
        //
        //    Note the surprising fact that, even for a Delaunay triangulation of
        //    a set of nodes, the nearest point to (X,Y) need not be one of the
        //    vertices of the triangle containing (X,Y).
        //
        //    The code can be called for triangulations of any order, but only
        //    the first three nodes in each triangle are considered.  Thus, if
        //    higher order triangles are used, and the extra nodes are intended
        //    to give the triangle a polygonal shape, these will have no effect,
        //    and the results obtained here might be misleading.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 October 2012
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
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
        //
        //    Input, int TRIANGLE_ORDER, the order of the triangles.
        //
        //    Input, int TRIANGLE_NUM, the number of triangles in the triangulation.
        //
        //    Input, int TRIANGLE_NODE[TRIANGLE_ORDER*TRIANGLE_NUM],
        //    the nodes of each triangle.
        //
        //    Input, int TRIANGLE_NEIGHBOR[3*TRIANGLE_NUM], the triangle neighbor list.
        //
        //    Input, double P[2], the coordinates of a point.
        //
        //    Output, int *TRIANGLE_INDEX, the index of the triangle where the search ended.
        //    If a cycle occurred, then TRIANGLE_INDEX = -1.
        //
        //    Output, double *ALPHA, *BETA, *GAMMA, the barycentric coordinates
        //    of the point with respect to triangle *TRIANGLE_INDEX.
        //
        //    Output, int *EDGE, indicates the position of the point (X,Y) in
        //    triangle TRIANGLE:
        //    0, the interior or boundary of the triangle;
        //    -1, outside the convex hull of the triangulation, past edge 1;
        //    -2, outside the convex hull of the triangulation, past edge 2;
        //    -3, outside the convex hull of the triangulation, past edge 3.
        //
        //    Output, int *STEP_NUM, the number of steps.
        {
            int a;
            int b;
            int c;
            double det;
            double dxp;
            double dxa;
            double dxb;
            double dyp;
            double dya;
            double dyb;

            step_num = -1;
            edge = 0;

            if (data.triangle_index_save < 0 || triangle_num <= data.triangle_index_save)
            {
                triangle_index = (triangle_num + 1) / 2;
            }
            else
            {
                triangle_index = data.triangle_index_save;
            }

            for (;;)
            {
                step_num = step_num + 1;

                if (triangle_num < step_num)
                {
                    Console.WriteLine("");
                    Console.WriteLine("TRIANGULATION_SEARCH_DELAUNAY - Fatal error!");
                    Console.WriteLine("  The algorithm seems to be cycling.");
                    Console.WriteLine("  Current triangle is " + triangle_index + "");
                    triangle_index = -1;
                    alpha = -1.0;
                    beta = -1.0;
                    gamma = -1.0;
                    edge = -1;
                    return;
                }

                //
                //  Get the vertices of triangle TRIANGLE.
                //
                a = triangle_node[0 + (triangle_index - 1) * triangle_order];
                b = triangle_node[1 + (triangle_index - 1) * triangle_order];
                c = triangle_node[2 + (triangle_index - 1) * triangle_order];
                //
                //  Using vertex C as a base, compute the distances to vertices A and B,
                //  and the point (X,Y).
                //
                dxa = node_xy[0 + a * 2] - node_xy[0 + c * 2];
                dya = node_xy[1 + a * 2] - node_xy[1 + c * 2];

                dxb = node_xy[0 + b * 2] - node_xy[0 + c * 2];
                dyb = node_xy[1 + b * 2] - node_xy[1 + c * 2];

                dxp = p[pIndex + 0] - node_xy[0 + c * 2];
                dyp = p[pIndex + 1] - node_xy[1 + c * 2];

                det = dxa * dyb - dya * dxb;
                //
                //  Compute the barycentric coordinates of the point (X,Y) with respect
                //  to this triangle.
                //
                alpha = (dxp * dyb - dyp * dxb) / det;
                beta = (dxa * dyp - dya * dxp) / det;
                gamma = 1.0 - alpha - beta;
                //
                //  If the barycentric coordinates are all positive, then the point
                //  is inside the triangle and we're done.
                //
                if (0.0 <= alpha &&
                    0.0 <= beta &&
                    0.0 <= gamma)
                {
                    break;
                }

                //
                //  At least one barycentric coordinate is negative.
                //
                //  If there is a negative barycentric coordinate for which there exists
                //  an opposing triangle neighbor closer to the point, move to that triangle.
                //
                //  (Two coordinates could be negative, in which case we could go for the
                //  most negative one, or the most negative one normalized by the actual
                //  distance it represents).
                //
                if (alpha < 0.0 && 0 <= triangle_neighbor[1 + (triangle_index - 1) * 3])
                {
                    triangle_index = triangle_neighbor[1 + (triangle_index - 1) * 3];
                    continue;
                }
                else if (beta < 0.0 && 0 <= triangle_neighbor[2 + (triangle_index - 1) * 3])
                {
                    triangle_index = triangle_neighbor[2 + (triangle_index - 1) * 3];
                    continue;
                }
                else if (gamma < 0.0 && 0 <= triangle_neighbor[0 + (triangle_index - 1) * 3])
                {
                    triangle_index = triangle_neighbor[0 + (triangle_index - 1) * 3];
                    continue;
                }

                //
                //  All negative barycentric coordinates correspond to vertices opposite
                //  sides on the convex hull.
                //
                //  Note the edge and exit.
                //
                if (alpha < 0.0)
                {
                    edge = -2;
                    break;
                }
                else if (beta < 0.0)
                {
                    edge = -3;
                    break;
                }
                else if (gamma < 0.0)
                {
                    edge = -1;
                    break;
                }
                else
                {
                    Console.WriteLine("");
                    Console.WriteLine("TRIANGULATION_ORDER3_SEARCH - Fatal error!");
                    Console.WriteLine("  The algorithm seems to have reached a dead end");
                    Console.WriteLine("  after " + step_num + " steps.");
                    triangle_index = -1;
                    edge = -1;
                    return;
                }
            }

            data.triangle_index_save = triangle_index;
        }

        public static int triangulation_search_naive(int node_num, double[] node_xy,
        int triangle_order, int triangle_num, int[] triangle_node, double[] p, int pIndex = 0 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGULATION_SEARCH_NAIVE naively searches a triangulation for a point.
        //
        //  Discussion:
        //
        //    The algorithm simply checks each triangle to see if point P is
        //    contained in it.  Surprisingly, this is not the fastest way to
        //    do the check, at least if the triangulation is Delaunay.
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
        //    Input, int TRIANGLE_ORDER, the order of the triangles.
        //
        //    Input, int TRIANGLE_NUM, the number of triangles in the triangulation.
        //
        //    Input, int TRIANGLE_NODE[TRIANGLE_ORDER*TRIANGLE_NUM],
        //    the nodes of each triangle.
        //
        //    Input, double P[2], the coordinates of a point.
        //
        //    Output, int TRIANGULATION_SEARCH_NAIVE, the index of the triangle
        //    containing the point, or -1 if no triangle was found containing
        //    the point.
        //
        {
            int a;
            double alpha;
            int b;
            double beta;
            int c;
            double det;
            double dxp;
            double dxa;
            double dxb;
            double dyp;
            double dya;
            double dyb;
            double gamma;
            int triangle;
            int triangle_index;

            triangle_index = -1;

            for (triangle = 0; triangle < triangle_num; triangle++)
            {
                //
                //  Get the vertices of triangle TRIANGLE.
                //
                a = triangle_node[0 + triangle * triangle_order];
                b = triangle_node[1 + triangle * triangle_order];
                c = triangle_node[2 + triangle * triangle_order];
                //
                //  Using vertex C as a base, compute the distances to vertices A and B,
                //  and the point (X,Y).
                //
                dxa = node_xy[0 + a * 2] - node_xy[0 + c * 2];
                dya = node_xy[1 + a * 2] - node_xy[1 + c * 2];

                dxb = node_xy[0 + b * 2] - node_xy[0 + c * 2];
                dyb = node_xy[1 + b * 2] - node_xy[1 + c * 2];

                dxp = p[pIndex + 0] - node_xy[0 + c * 2];
                dyp = p[pIndex + 1] - node_xy[1 + c * 2];

                det = dxa * dyb - dya * dxb;
                //
                //  Compute the barycentric coordinates of the point (X,Y) with respect
                //  to this triangle.
                //
                alpha = (dxp * dyb - dyp * dxb) / det;
                beta = (dxa * dyp - dya * dxp) / det;
                gamma = 1.0 - alpha - beta;
                //
                //  If the barycentric coordinates are all positive, then the point
                //  is inside the triangle and we're done.
                //
                if (0.0 <= alpha &&
                    0.0 <= beta &&
                    0.0 <= gamma)
                {
                    triangle_index = triangle + 1;
                    break;
                }
            }

            return triangle_index;
        }

    }
}