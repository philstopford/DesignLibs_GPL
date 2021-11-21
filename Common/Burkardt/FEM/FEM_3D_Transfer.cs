using Burkardt.SolveNS;
using Burkardt.TetrahedronNS;
using Burkardt.Types;

namespace Burkardt.FEM;

public static class FEM_3D_Transfer
{
    public static double[] fem3d_transfer(int sample_node_num, int sample_element_order,
            int sample_element_num, int sample_value_dim, int sample_value_num,
            double[] sample_node_xyz, int[] sample_element_node,
            int[] sample_element_neighbor, double[] sample_value, int fem_node_num,
            int fem_element_order, int fem_element_num, int fem_value_dim,
            int fem_value_num, double[] fem_node_xyz, int[] fem_element_node )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FEM3D_TRANSFER "transfers" from one finite element mesh to another.
        //
        //  BAD THINGS:
        //
        //    1) the linear system A*X=B is defined with A being a full storage matrix.
        //    2) the quadrature rule used is low order.
        //    3) the elements are assumed to be linear.
        //
        //  Discussion:
        //
        //    We are also given a set of "sample" finite element function defined
        //    by SAMPLE_NODE_XYZ, SAMPLE_ELEMENT, and SAMPLE_VALUE.
        //
        //    We are given a second finite element mesh, FEM_NODE_XYZ and 
        //    FEM_ELEMENT_NODE.
        //
        //    Our aim is to "project" the sample data values into the finite element 
        //    space, that is, to come up with a finite element function FEM_VALUE which
        //    well approximates the sample data.
        //
        //    Now let W(x,y,z) represent a function interpolating the sample data, and
        //    let Vijk(x,y,z) represent the finite element basis function associated with
        //    node IJK.
        //
        //    Then we seek the coefficient vector U corresponding to a finite element
        //    function U(x,y,z) of the form:
        //
        //      U(x,y,z) = sum ( 1 <= IJK <= N ) Uijk * Vijk(x,y,z)
        //
        //    To determine the coefficent vector entries U, we form a set of
        //    projection equations.  For node IJK at grid point (I,J,K), the associated
        //    basis function Vk(x,y,z) is used to pose the equation:
        //
        //      Integral U(x,y,z) Vijk(x,y,z) dx dy dz 
        //        = Integral W(x,y,z) Vijk(x,y,z) dx dy dz
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
        //    27 August 2009
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
        //    Input, double SAMPLE_NODE_XYZ[3*SAMPLE_NODE_NUM], the nodes.
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
        //    Input, double FEM_NODE_XYZ[3*FEM_NODE_NUM], the nodes.  
        //
        //    Input, int FEM_ELEMENT_NODE[FEM_ELEMENT_ORDER*FEM_ELEMENT_NUM],
        //    the nodes that make up each element.
        //
        //    Output, double FEM3D_TRANSFER[FEM_VALUE_DIM*FEM_VALUE_NUM], 
        //    the values.
        //
    {
        int element;
        int i;
        int j;
        const int project_node_num = 1;
        double[] project_node_xyz = new double[3 * 1];
        const int quad_num = 4;
        //
        //  Assemble the coefficient matrix A and the right-hand side B.
        //
        double[] b = typeMethods.r8mat_zero_new(fem_node_num, fem_value_dim);
        double[] a = typeMethods.r8mat_zero_new(fem_node_num, fem_node_num);

        double[] phi = new double[4];
        double[] ref_weight = new double[quad_num];
        double[] ref_quad = new double[4 * quad_num];
        double[] tet_quad = new double[3 * quad_num];
        double[] tet_xyz = new double[3 * 4];

        // Upstream doesn't initialize the arrays and seems to rely on weird default values.
        // Mimic here.
        for (int shim = 0; shim < phi.Length; shim++)
        {
            phi[shim] = -6.2774385622041925e+66;
        }
        for (int shim = 0; shim < ref_weight.Length; shim++)
        {
            ref_weight[shim] = -6.2774385622041925e+66;
        }
        for (int shim = 0; shim < ref_quad.Length; shim++)
        {
            ref_quad[shim] = -6.2774385622041925e+66;
        }
        for (int shim = 0; shim < tet_quad.Length; shim++)
        {
            tet_quad[shim] = -6.2774385622041925e+66;
        }
        for (int shim = 0; shim < tet_xyz.Length; shim++)
        {
            tet_xyz[shim] = -6.2774385622041925e+66;
        }
            
        for (element = 0; element < fem_element_num; element++)
        {
            for (j = 0; j < 4; j++)
            {
                for (i = 0; i < 3; i++)
                {
                    int j2 = fem_element_node[j + element * 4];
                    tet_xyz[i + j * 3] = fem_node_xyz[i + j2 * 3];
                }
            }

            double volume = Tetrahedron.tetrahedron_volume(tet_xyz);

            for (j = 0; j < quad_num; j++)
            {
                for (i = 0; i < 3; i++)
                {
                    tet_quad[i + j * 3] = 0.0;
                    int k;
                    for (k = 0; k < 4; k++)
                    {
                        double tq = tet_quad[i + j * 3];
                        double tx = tet_xyz[i + k * 3];
                        double rq = ref_quad[k + j * 4];
                            
                        tet_quad[i + j * 3] = tq + tx * rq;
                    }
                }
            }

            //
            //  Consider each quadrature point.
            //  Here, we use the midside nodes as quadrature points.
            //
            int quad;
            for (quad = 0; quad < quad_num; quad++)
            {
                for (i = 0; i < 3; i++)
                {
                    project_node_xyz[i + 0 * 3] = tet_quad[i + quad * 3];
                }

                Basis_mn.basis_mn_tet4(tet_xyz, 1, project_node_xyz, ref phi);

                for (i = 0; i < 4; i++)
                {
                    int ni = fem_element_node[i + element * fem_element_order];
                    //
                    //  The projection takes place here.  The finite element code needs the value
                    //  of the sample function at the point (XQ,YQ).  The call to PROJECTION
                    //  locates (XQ,YQ) in the triangulated mesh of sample data, and returns a
                    //  value produced by piecewise linear interpolation.
                    //
                    double[] project_value = FEM_3D_Projection.projection(sample_node_num, sample_node_xyz,
                        sample_element_order, sample_element_num, sample_element_node,
                        sample_element_neighbor, sample_value_dim, sample_value,
                        project_node_num, project_node_xyz);

                    for (j = 0; j < fem_value_dim; j++)
                    {
                        b[ni + j * fem_node_num] += volume * ref_weight[quad] *
                                                    (project_value[j + 0 * fem_value_dim] * phi[i]);
                    }

                    //
                    //  Consider each basis function in the element.
                    //
                    for (j = 0; j < 4; j++)
                    {
                        int nj = fem_element_node[j + element * fem_element_order];

                        a[ni + nj * fem_node_num] += volume * ref_weight[quad] * (phi[i] * phi[j]);
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
            for (i = 0; i < fem_value_dim; i++)
            {
                fem_value[(i + j * fem_value_dim + fem_value.Length) % fem_value.Length] = x[(j + i * fem_value_num + x.Length) % x.Length];
            }
        }

        return fem_value;
    }
}