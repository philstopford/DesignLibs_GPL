using System;
using Burkardt.Types;

namespace Burkardt.TriangulationNS
{
    public static partial class Quad
    {
        public static void triangulation_order3_quad(int node_num, double[] node_xy,
        int triangle_order, int triangle_num, int[] triangle_node,
        Func<int, double[], double[], double[] > quad_fun, int quad_num,
        double[] quad_xy, double[] quad_w, ref double quad_value, ref double region_area )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGULATION_ORDER3_QUAD approximates an integral over a triangulation.
        //
        //  Discussion:
        //
        //    The routine will accept triangulations of order higher than 3.
        //    However, only the first three nodes (the vertices) of each
        //    triangle will be used.  This will still produce correct results
        //    for higher order triangulations, as long as the sides of the
        //    triangle are straight.
        //
        //    We assume that the vertices of each triangle are listed first
        //    in the description of higher order triangles, and we assume that
        //    the vertices are listed in counterclockwise order.
        //
        //    The approximation of the integral is made using a quadrature rule
        //    defined on the unit triangle, and supplied by the user.
        //
        //    The user also supplies the name of a subroutine, here called "QUAD_FUN",
        //    which evaluates the integrand at a set of points.  The form is:
        //
        //      void quad_fun ( int n, double xy_vec[], double f_vec[] )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NODE_NUM, the number of nodes in the triangulation.
        //
        //    Input, double NODE_XY(2,NODE_NUM), the coordinates of the nodes.
        //
        //    Input, int TRIANGLE_ORDER, the order of triangles in the triangulation.
        //
        //    Input, int TRIANGLE_NUM, the number of triangles in the triangulation.
        //
        //    Input, int TRIANGLE_NODE[TRIANGLE_ORDER*TRIANGLE_NUM],
        //    the nodes making up each triangle.
        //
        //    Input, void QUAD_FUN ( int N, double XY_VEC[], double F_VEC[] ),
        //    the name of the function that evaluates the integrand.
        //
        //    Input, int QUAD_NUM, the order of the quadrature rule.
        //
        //    Input, double QUAD_XY(2,QUAD_NUM), the abscissas of the
        //    quadrature rule, in the unit triangle.
        //
        //    Input, double QUAD_W(QUAD_NUM), the weights of the
        //    quadrature rule.
        //
        //    Output, double *QUAD_VALUE, the estimate of the integral
        //    of F(X,Y) over the region covered by the triangulation.
        //
        //    Output, double *REGION_AREA, the area of the region.
        //
        {
            int i;
            int j;
            int quad;
            double[] quad_f;
            double[] quad2_xy;
            double temp;
            int triangle;
            double triangle_area;
            double[] triangle_xy = new double[2 * 3];

            quad_f = new double[quad_num];
            quad2_xy = new double[2 * quad_num];

            quad_value = 0.0;
            region_area = 0.0;

            for (triangle = 0; triangle < triangle_num; triangle++)
            {
                for (j = 0; j < 3; j++)
                {
                    for (i = 0; i < 2; i++)
                    {
                        triangle_xy[i + j * 2] = node_xy[i + (triangle_node[j + triangle * 3] - 1) * 2];
                    }
                }

                triangle_area = typeMethods.triangle_area_2d(triangle_xy);

                Triangulation.triangle_order3_reference_to_physical(triangle_xy,
                    quad_num, quad_xy, ref quad2_xy);

                quad_fun(quad_num, quad2_xy, quad_f);

                temp = 0.0;
                for (quad = 0; quad < quad_num; quad++)
                {
                    temp = temp + quad_w[quad] * quad_f[quad];
                }

                quad_value = quad_value + triangle_area * temp;

                region_area = region_area + triangle_area;
            }
        }

    }
}