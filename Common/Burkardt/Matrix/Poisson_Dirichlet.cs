using System;
using Burkardt.FEM;
using Burkardt.Quadrature;
using Burkardt.Types;

namespace Burkardt.MatrixNS;

public static partial class Matrix
{
    public static void assemble_poisson(int node_num, double[] node_xy,
            int element_num, int[] element_node, int quad_num, int ib, ref double[] a,
            ref double[] f, Func<int, double[], double[]> rhs, Func<int, double[], double[]> h_coef,
            Func<int, double[], double[]> k_coef)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ASSEMBLE_POISSON assembles the system for the Poisson equation.
        //
        //  Discussion:
        //
        //    The matrix is known to be banded.  A special matrix storage format
        //    is used to reduce the space required.  Details of this format are
        //    discussed in the routine DGB_FA.
        //
        //    Note that a 3 point quadrature rule, which is sometimes used to
        //    assemble the matrix and right hand side, is just barely accurate
        //    enough for simple problems.  If you want better results, you
        //    should use a quadrature rule that is more accurate.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 December 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, double NODE_XY[2*NODE_NUM], the coordinates of nodes.
        //
        //    Input, int ELEMENT_NUM, the number of elements.
        //
        //    Input, int ELEMENT_NODE[3*ELEMENT_NUM];
        //    ELEMENT_NODE(I,J) is the global index of local node I in element J.
        //
        //    Input, int QUAD_NUM, the number of quadrature points used in assembly.
        //
        //    Input, int IB, the half-bandwidth of the matrix.
        //
        //    Output, double A(3*IB+1,NODE_NUM), the NODE_NUM by NODE_NUM
        //    coefficient matrix, stored in a compressed format.
        //
        //    Output, double F(NODE_NUM), the right hand side.
        //
        //  Local parameters:
        //
        //    Local, double BI, DBIDX, DBIDY, the value of some basis function
        //    and its first derivatives at a quadrature point.
        //
        //    Local, double BJ, DBJDX, DBJDY, the value of another basis
        //    function and its first derivatives at a quadrature point.
        //
    {
        double bi = 0;
        double bj = 0;
        double dbidx = 0;
        double dbidy = 0;
        double dbjdx = 0;
        double dbjdy = 0;
        int element;
        int i;
        int node;
        double[] p = new double[2];
        double[] t3 = new double[2 * 3];

        double[] phys_h = new double[quad_num];
        double[] phys_k = new double[quad_num];
        double[] phys_rhs = new double[quad_num];
        double[] phys_xy = new double[2 * quad_num];
        double[] quad_w = new double[quad_num];
        double[] quad_xy = new double[2 * quad_num];
        double[] w = new double[quad_num];
        //
        //  Initialize the arrays to zero.
        //
        for (node = 0; node < node_num; node++)
        {
            f[node] = 0.0;
        }

        for (node = 0; node < node_num; node++)
        {
            for (i = 0; i < 3 * ib + 1; i++)
            {
                a[i + node * (3 * ib + 1)] = 0.0;
            }
        }

        //
        //  Get the quadrature weights and nodes.
        //
        QuadratureRule.quad_rule(quad_num, ref quad_w, ref quad_xy);
        //
        //  Add up all quantities associated with the ELEMENT-th element.
        //
        for (element = 0; element < element_num; element++)
        {
            //
            //  Make a copy of the element.
            //
            int j;
            for (j = 0; j < 3; j++)
            {
                for (i = 0; i < 2; i++)
                {
                    t3[i + j * 2] = node_xy[i + (element_node[j + element * 3] - 1) * 2];
                }
            }

            //
            //  Map the quadrature points QUAD_XY to points PHYS_XY in the physical element.
            //
            Reference.reference_to_physical_t3(t3, quad_num, quad_xy, ref phys_xy);

            double area = Math.Abs(typeMethods.triangle_area_2d(t3));

            int quad;
            for (quad = 0; quad < quad_num; quad++)
            {
                w[quad] = quad_w[quad] * area;
            }

            phys_rhs = rhs(quad_num, phys_xy);
            phys_h = h_coef(quad_num, phys_xy);
            phys_k = k_coef(quad_num, phys_xy);
            //
            //  Consider the QUAD-th quadrature point.
            //
            for (quad = 0; quad < quad_num; quad++)
            {
                p[0] = phys_xy[0 + quad * 2];
                p[1] = phys_xy[1 + quad * 2];
                //
                //  Consider the TEST-th test function.
                //
                //  We generate an integral for every node associated with an unknown.
                //  But if a node is associated with a boundary condition, we do nothing.
                //
                int test;
                for (test = 1; test <= 3; test++)
                {
                    i = element_node[test - 1 + element * 3];

                    Basis11.basis_one_t3(t3, test, p, ref bi, ref dbidx, ref dbidy);

                    f[i - 1] += w[quad] * phys_rhs[quad] * bi;
                    //
                    //  Consider the BASIS-th basis function, which is used to form the
                    //  value of the solution function.
                    //
                    int basis;
                    for (basis = 1; basis <= 3; basis++)
                    {
                        j = element_node[basis - 1 + element * 3];

                        Basis11.basis_one_t3(t3, basis, p, ref bj, ref dbjdx, ref dbjdy);

                        a[i - j + 2 * ib + (j - 1) * (3 * ib + 1)] += w[quad] * (
                            phys_h[quad] *
                            (dbidx * dbjdx + dbidy * dbjdy)
                            + phys_k[quad] * bj * bi);
                    }
                }
            }
        }
    }

    public static void dirichlet_apply(int node_num, double[] node_xy, int[] node_condition,
            int ib, ref double[] a, ref double[] f, Func<int, double[], double[]> dirichlet_condition)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIRICHLET_APPLY accounts for Dirichlet boundary conditions.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 July 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, double NODE_XY[2*NODE_NUM], the coordinates of nodes.
        //
        //    Input, int NODE_CONDITION[NODE_NUM], reports the condition
        //    used to set the unknown associated with the node.
        //    0, unknown.
        //    1, finite element equation.
        //    2, Dirichlet condition;
        //    3, Neumann condition.
        //
        //    Input, int IB, the half-bandwidth of the matrix.
        //
        //    Input/output, double A[(3*IB+1)*NODE_NUM], the NODE_NUM by
        //    NODE_NUM coefficient matrix, stored in a compressed format; on output,
        //    the matrix has been adjusted for Dirichlet boundary conditions.
        //
        //    Input/output, double F[NODE_NUM], the right hand side.
        //    On output, the right hand side has been adjusted for Dirichlet
        //    boundary conditions.
        //
    {
        const int DIRICHLET = 2;
        int node;

        double[] node_bc = dirichlet_condition(node_num, node_xy);

        for (node = 0; node < node_num; node++)
        {
            if (node_condition[node] == DIRICHLET)
            {
                int column_low = Math.Max(node + 1 - ib, 1);
                int column_high = Math.Min(node + 1 + ib, node_num);

                int column;
                for (column = column_low; column <= column_high; column++)
                {
                    a[node + 1 - column + 2 * ib + (column - 1) * (3 * ib + 1)] = 0.0;
                }

                a[2 * ib + node * (3 * ib + 1)] = 1.0;

                f[node] = node_bc[node];
            }
        }
    }

    public static double[] residual_poisson(int node_num, double[] node_xy, int[] node_condition,
            int element_num, int[] element_node, int quad_num, int ib, double[] a,
            double[] f, double[] node_u, Func<int, double[], double[]> rhs, Func<int, double[], double[]> h_coef,
            Func<int, double[], double[]> k_coef, Func<int, double[], double[]> dirichlet_condition)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RESIDUAL_POISSON evaluates the residual for the Poisson equation.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 November 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, double NODE_XY[2*NODE_NUM], the
        //    coordinates of nodes.
        //
        //    Input, int NODE_CONDITION[NODE_NUM], reports the condition
        //    used to set the unknown associated with the node.
        //    0, unknown.
        //    1, finite element equation.
        //    2, Dirichlet condition;
        //    3, Neumann condition.
        //
        //    Input, int ELEMENT_NUM, the number of elements.
        //
        //    Input, int ELEMENT_NODE[3*ELEMENT_NUM];
        //    ELEMENT_NODE(I,J) is the global index of local node I in element J.
        //
        //    Input, int QUAD_NUM, the number of quadrature points used in assembly.
        //
        //    Input, int IB, the half-bandwidth of the matrix.
        //
        //    Workspace, double A[(3*IB+1)*NODE_NUM], the NODE_NUM by NODE_NUM
        //    coefficient matrix, stored in a compressed format.
        //
        //    Workspace, double F[NODE_NUM], the right hand side.
        //
        //    Input, double NODE_U[NODE_NUM], the value of the solution
        //    at each node.
        //
        //    Output, double NODE_R[NODE_NUM], the finite element
        //    residual at each node.
        //
        //  Local parameters:
        //
        //    Local, double BI, DBIDX, DBIDY, the value of some basis function
        //    and its first derivatives at a quadrature point.
        //
        //    Local, double BJ, DBJDX, DBJDY, the value of another basis
        //    function and its first derivatives at a quadrature point.
        //
    {
        double bi = 0;
        double bj = 0;
        double dbidx = 0;
        double dbidy = 0;
        double dbjdx = 0;
        double dbjdy = 0;
        int element;
        int i;
        int node;
        double[] p = new double[2];
        double[] t3 = new double[2 * 3];

        double[] phys_h = new double[quad_num];
        double[] phys_k = new double[quad_num];
        double[] phys_rhs = new double[quad_num];
        double[] phys_xy = new double[2 * quad_num];

        double[] quad_w = new double[quad_num];
        double[] quad_xy = new double[2 * quad_num];

        double[] w = new double[quad_num];
        //
        //  Initialize the arrays to zero.
        //
        for (node = 0; node < node_num; node++)
        {
            f[node] = 0.0;
        }

        for (node = 0; node < node_num; node++)
        {
            for (i = 0; i < 3 * ib + 1; i++)
            {
                a[i + node * (3 * ib + 1)] = 0.0;
            }
        }

        //
        //  Get the quadrature weights and nodes.
        //
        QuadratureRule.quad_rule(quad_num, ref quad_w, ref quad_xy);
        //
        //  The actual values of A and F are determined by summing up
        //  contributions from all the elements.
        //
        for (element = 0; element < element_num; element++)
        {
            //
            //  Make a copy of the element.
            //
            int j;
            for (j = 0; j < 3; j++)
            {
                for (i = 0; i < 2; i++)
                {
                    t3[i + j * 2] = node_xy[i + (element_node[j + element * 3] - 1) * 2];
                }
            }

            //
            //  Map the quadrature points QUAD_XY to points XY in the physical element.
            //
            Reference.reference_to_physical_t3(t3, quad_num, quad_xy, ref phys_xy);

            double area = Math.Abs(typeMethods.triangle_area_2d(t3));

            int quad;
            for (quad = 0; quad < quad_num; quad++)
            {
                w[quad] = quad_w[quad] * area;
            }

            phys_rhs = rhs(quad_num, phys_xy);
            phys_h = h_coef(quad_num, phys_xy);
            phys_k = k_coef(quad_num, phys_xy);
            //
            //  Consider a quadrature point QUAD, with coordinates (X,Y).
            //
            for (quad = 0; quad < quad_num; quad++)
            {
                p[0] = phys_xy[0 + quad * 2];
                p[1] = phys_xy[1 + quad * 2];
                //
                //  Consider one of the basis functions, which will play the
                //  role of test function in the integral.
                //
                //  We generate an integral for every node associated with an unknown.
                //  But if a node is associated with a boundary condition, we do nothing.
                //
                int test;
                for (test = 1; test <= 3; test++)
                {
                    i = element_node[test - 1 + element * 3];

                    Basis11.basis_one_t3(t3, test, p, ref bi, ref dbidx, ref dbidy);

                    f[i - 1] += w[quad] * phys_rhs[quad] * bi;
                    //
                    //  Consider another basis function, which is used to form the
                    //  value of the solution function.
                    //
                    int basis;
                    for (basis = 1; basis <= 3; basis++)
                    {
                        j = element_node[basis - 1 + element * 3];

                        Basis11.basis_one_t3(t3, basis, p, ref bj, ref dbjdx, ref dbjdy);

                        a[i - j + 2 * ib + (j - 1) * (3 * ib + 1)] += w[quad]
                                                                      * (phys_h[quad] *
                                                                         (dbidx * dbjdx + dbidy * dbjdy)
                                                                         + phys_k[quad] * bj * bi);
                    }
                }
            }
        }

        //
        //  Apply boundary conditions.
        //
        dirichlet_apply(node_num, node_xy, node_condition, ib, ref a, ref f, dirichlet_condition);
        //
        //  Compute A*U.
        //
        double[] node_r = dgb_mxv(node_num, node_num, ib, ib, a, node_u);
        //
        //  Set RES = A * U - F.
        //
        for (node = 0; node < node_num; node++)
        {
            node_r[node] -= f[node];
        }

        return node_r;
    }
}