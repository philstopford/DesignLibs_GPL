using Burkardt.SolveNS;
using Burkardt.Types;

namespace Burkardt.FEM;

public static class FEM_2D_Transfer
{
    public static double[] fem2d_transfer(int sample_node_num, int sample_element_order,
            int sample_element_num, int sample_value_dim, int sample_value_num,
            double[] sample_node_xy, int[] sample_element_node,
            int[] sample_element_neighbor, double[] sample_value, int fem_node_num,
            int fem_element_order, int fem_element_num, int fem_value_dim,
            int fem_value_num, double[] fem_node_xy, int[] fem_element_node)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FEM2D_TRANSFER "transfers" from one finite element mesh to another.
        //
        //  BAD THINGS:
        //
        //    1) the linear system A*X=B is defined with A being a full storage matrix.
        //    2) the quadrature rule used is low order.
        //    3) the triangular elements are assumed to be linear.
        //
        //  Discussion:
        //
        //    We are also given a set of "sample" finite element function defined
        //    by SAMPLE_NODE_XY, SAMPLE_ELEMENT, and SAMPLE_VALUE.
        //
        //    We are given a second finite element mesh, FEM_NODE_XY and 
        //    FEM_ELEMENT_NODE.
        //
        //    Our aim is to "project" the sample data values into the finite element 
        //    space, that is, to come up with a finite element function FEM_VALUE which
        //    well approximates the sample data.
        //
        //    Now let W(x,y) represent a function interpolating the sample data, and
        //    let Vk(x,y) represent the finite element basis function associated with
        //    node K.
        //
        //    Then we seek the coefficient vector U corresponding to a finite element 
        //    function U(x,y) of the form:
        //
        //      U(x,y) = sum ( 1 <= K <= N ) Uk * Vk(x,y)
        //
        //    To determine the coefficent vector entries U, we form a set of 
        //    projection equations.  For node K at grid point (I,J), the associated 
        //    basis function Vk(x,y) is used to pose the equation:
        //
        //      Integral U(x,y) Vk(x,y) dx dy = Integral W(x,y) Vk(x,y) dx dy
        //
        //    The left hand side is the usual stiffness matrix times the desired
        //    coefficient vector U.  To complete the system, we simply need to
        //    determine the right hand side, that is, the integral of the data function
        //    W against the basis function Vk.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 August 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int SAMPLE_NODE_NUM, the number of nodes.
        //
        //    Input, int SAMPLE_ELEMENT_ORDER, the element order.
        //
        //    Input, int SAMPLE_ELEMENT_NUM, the number of elements.
        //
        //    Input, int SAMPLE_VALUE_DIM, the value dimension.
        //
        //    Input, int SAMPLE_VALUE_NUM, the number of values.
        //
        //    Input, double SAMPLE_NODE_XY[2*SAMPLE_NODE_NUM], the nodes.
        //
        //    Input, int SAMPLE_ELEMENT_NODE[SAMPLE_ELEMENT_ORDER*SAMPLE_ELEMENT_NUM],
        //    the nodes that make up each element.
        //
        //    Input, int SAMPLE_ELEMENT_NEIGHBOR[3*SAMPLE_ELEMENT_NUM],
        //    the neighbor triangles.
        //
        //    Input, double SAMPLE_VALUE[SAMPLE_VALUE_DIM*SAMPLE_NODE_NUM],
        //    the values.
        //
        //    Input, int FEM_NODE_NUM, the number of nodes.
        //
        //    Input, int FEM_ELEMENT_ORDER, the element order.
        //
        //    Input, int FEM_ELEMENT_NUM, the number of elements.
        //
        //    Input, int FEM_VALUE_DIM, the value dimension.
        //
        //    Input, int FEM_VALUE_NUM, the number of values.
        //
        //    Input, double FEM_NODE_XY[2*FEM_NODE_NUM], the nodes.  
        //
        //    Input, int FEM_ELEMENT_NODE[FEM_ELEMENT_ORDER*FEM_ELEMENT_NUM],
        //    the nodes that make up each element.
        //
        //    Output, double FEM2D_TRANSFER[FEM_VALUE_DIM*FEM_VALUE_NUM], 
        //    the values.
        //
    {
        int element;
        int j;
        const int project_node_num = 1;
        double[] project_node_xy = new double[2 * 1];
        const int quad_num = 3;
        //
        //  Assemble the coefficient matrix A and the right-hand side B.
        //
        double[] b = typeMethods.r8mat_zero_new(fem_node_num, fem_value_dim);
        double[] a = typeMethods.r8mat_zero_new(fem_node_num, fem_node_num);

        for (element = 0; element < fem_element_num; element++)
        {
            int i1 = fem_element_node[0 + element * 3];
            int i2 = fem_element_node[1 + element * 3];
            int i3 = fem_element_node[2 + element * 3];

            double area = 0.5 *
                          (fem_node_xy[0 + i1 * 2] * (fem_node_xy[1 + i2 * 2] - fem_node_xy[1 + i3 * 2])
                           + fem_node_xy[0 + i2 * 2] * (fem_node_xy[1 + i3 * 2] - fem_node_xy[1 + i1 * 2])
                           + fem_node_xy[0 + i3 * 2] * (fem_node_xy[1 + i1 * 2] - fem_node_xy[1 + i2 * 2]));
            //
            //  Consider each quadrature point.
            //  Here, we use the midside nodes as quadrature points.
            //
            int quad;
            for (quad = 0; quad < quad_num; quad++)
            {
                int q2 = (quad + 1) % quad_num;

                int nq1 = fem_element_node[quad + element * fem_element_order];
                int nq2 = fem_element_node[q2 + element * fem_element_order];

                double xq = 0.5 * (fem_node_xy[0 + nq1 * 2] + fem_node_xy[0 + nq2 * 2]);
                double yq = 0.5 * (fem_node_xy[1 + nq1 * 2] + fem_node_xy[1 + nq2 * 2]);
                const double wq = 1.0 / 3.0;
                //
                //  Consider each test function in the element.
                //
                int ti1;
                for (ti1 = 0; ti1 < 3; ti1++)
                {
                    int ti2 = (ti1 + 1) % 3;
                    int ti3 = (ti1 + 2) % 3;

                    int nti1 = fem_element_node[ti1 + element * fem_element_order];
                    int nti2 = fem_element_node[ti2 + element * fem_element_order];
                    int nti3 = fem_element_node[ti3 + element * fem_element_order];

                    double qi = 0.5 * (
                        (fem_node_xy[0 + nti3 * 2] - fem_node_xy[0 + nti2 * 2])
                        * (yq - fem_node_xy[1 + nti2 * 2])
                        - (fem_node_xy[1 + nti3 * 2] - fem_node_xy[1 + nti2 * 2])
                        * (xq - fem_node_xy[0 + nti2 * 2])) / area;
                    //
                    //  The projection takes place here.  The finite element code needs the value
                    //  of the sample function at the point (XQ,YQ).  The call to PROJECTION
                    //  locates (XQ,YQ) in the triangulated mesh of sample data, and returns a
                    //  value produced by piecewise linear interpolation.
                    //
                    project_node_xy[0 + 0 * 2] = xq;
                    project_node_xy[1 + 0 * 2] = yq;

                    double[] project_value = FEM_2D_Projection.projection(sample_node_num, sample_node_xy,
                        sample_element_order, sample_element_num, sample_element_node,
                        sample_element_neighbor, sample_value_dim, sample_value,
                        project_node_num, project_node_xy);

                    for (j = 0; j < fem_value_dim; j++)
                    {
                        b[nti1 + j * fem_node_num] += area * wq * (project_value[j + 0 * fem_value_dim] * qi);
                    }

                    //
                    //  Consider each basis function in the element.
                    //
                    int tj1;
                    for (tj1 = 0; tj1 < 3; tj1++)
                    {
                        int tj2 = (tj1 + 1) % 3;
                        int tj3 = (tj1 + 2) % 3;

                        int ntj1 = fem_element_node[tj1 + element * fem_element_order];
                        int ntj2 = fem_element_node[tj2 + element * fem_element_order];
                        int ntj3 = fem_element_node[tj3 + element * fem_element_order];

                        double qj = 0.5 * (
                            (fem_node_xy[0 + ntj3 * 2] - fem_node_xy[0 + ntj2 * 2])
                            * (yq - fem_node_xy[1 + ntj2 * 2])
                            - (fem_node_xy[1 + ntj3 * 2] - fem_node_xy[1 + ntj2 * 2])
                            * (xq - fem_node_xy[0 + ntj2 * 2])) / area;

                        a[nti1 + ntj1 * fem_node_num] += area * wq * (qi * qj);
                    }
                }
            }
        }

        //
        //  SOLVE the linear system A * X = B.
        //
        double[] x = Solve.r8ge_fss_new(fem_node_num, a, fem_value_dim, b);
        //
        //  Copy solution.
        //
        double[] fem_value = new double[fem_value_dim * fem_value_num];

        for (j = 0; j < fem_value_num; j++)
        {
            int i;
            for (i = 0; i < fem_value_dim; i++)
            {
                fem_value[i + j * fem_value_dim] = x[j + i * fem_value_num];
            }
        }

        return fem_value;
    }
}