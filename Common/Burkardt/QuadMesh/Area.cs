using System;
using Burkardt.TriangleNS;
using Burkardt.Types;

namespace Burkardt.QuadMesh;

public static class Area
{
    public static void area_q4_mesh(int node_num, int element_num, double[] node_xy,
            int[] element_node, double[] element_area, ref double mesh_area)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    AREA_Q4_MESH computes areas of elements in a Q4 mesh.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    24 February 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, int ELEMENT_NUM, the number of elements.
        //
        //    Input, double NODE_XY[2*NODE_NUM], the node coordinates.
        //
        //    Input, int ELEMENT_NODE[4*ELEMENT_NUM], lists the 
        //    nodes that make up each element, in counterclockwise order.
        //
        //    Output, double ELEMENT_AREA[ELEMENT_NUM], the element areas.
        //
        //    Output, double *MESH_AREA, the mesh area.
        //
    {
        int element;
        double[] q4 = new double[2 * 4];

        for (element = 0; element < element_num; element++)
        {
            int node;
            for (node = 0; node < 4; node++)
            {
                int dim;
                for (dim = 0; dim < 2; dim++)
                {
                    q4[dim + 2 * node] = node_xy[dim + 2 * element_node[node + 4 * element]];
                }
            }

            element_area[element] = area_quad(q4);
        }

        mesh_area = typeMethods.r8vec_sum(element_num, element_area);
    }

    public static double area_quad(double[] quad_xy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    AREA_QUAD returns the area of a quadrilateral.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 February 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double QUAD_XY[2*4], the coordinates of the nodes.
        //
        //    Output, double AREA_QUAD, the area.
        //
    {
        double[] t1 = new double[2 * 3];
        double[] t2 = new double[2 * 3];

        t1[0 + 0 * 2] = quad_xy[0 + 0 * 2];
        t1[1 + 0 * 2] = quad_xy[1 + 0 * 2];
        t1[0 + 1 * 2] = quad_xy[0 + 1 * 2];
        t1[1 + 1 * 2] = quad_xy[1 + 1 * 2];
        t1[0 + 2 * 2] = quad_xy[0 + 2 * 2];
        t1[1 + 2 * 2] = quad_xy[1 + 2 * 2];

        double area1 = Integrals.triangle_area(t1);

        t2[0 + 0 * 2] = quad_xy[0 + 2 * 2];
        t2[1 + 0 * 2] = quad_xy[1 + 2 * 2];
        t2[0 + 1 * 2] = quad_xy[0 + 3 * 2];
        t2[1 + 1 * 2] = quad_xy[1 + 3 * 2];
        t2[0 + 2 * 2] = quad_xy[0 + 0 * 2];
        t2[1 + 2 * 2] = quad_xy[1 + 0 * 2];

        double area2 = Integrals.triangle_area(t2);

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
                Console.WriteLine("AREA_QUAD - Fatal error!");
                Console.WriteLine("  The quadrilateral nodes seem to be listed in");
                Console.WriteLine("  the wrong order, or the quadrilateral is");
                Console.WriteLine("  degenerate.");
                return 0;
            }
        }

        double area = area1 + area2;

        return area;
    }
}