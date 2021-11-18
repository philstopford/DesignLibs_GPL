namespace Burkardt.TriangulationNS;

public static partial class Triangulation
{
    public static double triangulation_area(int node_num, double[] node_xy, int element_order,
            int element_num, int[] element_node )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGULATION_AREA computes the area of a triangulation.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 December 2011
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
        //    Input, int ELEMENT_ORDER, the order of the triangles.
        //
        //    Input, int ELEMENT_NUM, the number of triangles.
        //
        //    Input, int ELEMENT_NODE[ELEMENT_ORDER*ELEMENT_NUM],
        //    the nodes making up each triangle.
        //
        //    Output, double TRIANGULATION_AREA, the area.
        //
    {
        int element;
        double[] element_xy = new double[2 * 3];

        double value = 0.0;

        for (element = 0; element < element_num; element++)
        {
            int j;
            for (j = 0; j < 3; j++)
            {
                int nj = element_node[j + element * element_order];
                element_xy[0 + j * 2] = node_xy[0 + nj * 2];
                element_xy[1 + j * 2] = node_xy[1 + nj * 2];
            }

            double element_area = 0.5 * (
                element_xy[0 + 0 * 2] * (element_xy[1 + 1 * 2] - element_xy[1 + 2 * 2]) +
                element_xy[0 + 1 * 2] * (element_xy[1 + 2 * 2] - element_xy[1 + 0 * 2]) +
                element_xy[0 + 2 * 2] * (element_xy[1 + 0 * 2] - element_xy[1 + 1 * 2]));

            value += element_area;
        }

        return value;
    }

    public static double triangulation_areas(int node_num, double[] node_xy, int triangle_order,
            int triangle_num, int[] triangle_node, double[] triangle_area )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGULATION_AREAS computes triangle and triangulation areas.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 September 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NODE_NUM, the number of nodes in the
        //    triangulation.
        //
        //    Input, double NODE_XY[2*NODE_NUM], the coordinates of the nodes.
        //
        //    Input, int TRIANGLE_ORDER, the order of triangles in
        //    the triangulation.
        //
        //    Input, int TRIANGLE_NUM, the number of triangles in
        //    the triangulation.
        //
        //    Input, int TRIANGLE_NODE[TRIANGLE_ORDER*TRIANGLE_NUM],
        //    the nodes making up each triangle.
        //
        //    Output, double TRIANGLE_AREA[TRIANGLE_NUM], the area of
        //    the triangles.
        //
        //    Output, double TRIANGULATION_AREAS, the area of the triangulation.
        //
    {
        double[] t_xy = new double[2 * 3];
        int triangle;

        double triangulation_area = 0.0;

        for (triangle = 0; triangle < triangle_num; triangle++)
        {
            int j;
            for (j = 0; j < 3; j++)
            {
                int nj = triangle_node[j + triangle * triangle_order];
                t_xy[0 + j * 2] = node_xy[0 + nj * 2];
                t_xy[1 + j * 2] = node_xy[1 + nj * 2];
            }

            triangle_area[triangle] = 0.5 * (
                t_xy[0 + 0 * 2] * (t_xy[1 + 1 * 2] - t_xy[1 + 2 * 2]) +
                t_xy[0 + 1 * 2] * (t_xy[1 + 2 * 2] - t_xy[1 + 0 * 2]) +
                t_xy[0 + 2 * 2] * (t_xy[1 + 0 * 2] - t_xy[1 + 1 * 2]));

            triangulation_area += triangle_area[triangle];
        }

        return triangulation_area;
    }
}