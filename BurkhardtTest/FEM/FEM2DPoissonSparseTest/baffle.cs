using System;

namespace FEM2DPoissonSparseTest;

public static class baffle
{
    public static double[] dirichlet_condition(int node_num, double[] node_xy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIRICHLET_CONDITION sets the value of a Dirichlet boundary condition.
        //
        //  Discussion:
        //
        //    This routine specifies that the solution is zero on the boundary.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    25 January 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, double NODE_XY[2*NODE_NUM], the coordinates of the points.
        //
        //    Output, double NODE_BC[NODE_NUM], the value of the 
        //    Dirichlet boundary conditions at the points.
        //
    {
        int node;
        double[] node_bc = new double[node_num];
        for (node = 0; node < node_num; node++)
        {
            node_bc[node] = 0.0;
        }

        return node_bc;
    }

    public static double[] h_coef(int node_num, double[] node_xy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    H_COEF evaluates the coefficient H(X,Y) of DEL U in the Poisson equation.
        //
        //  Discussion:
        //
        //    The equation is 
        //
        //      -DEL H(X,Y) DEL U(X,Y) + K(X,Y) * U(X,Y) = F(X,Y)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    25 January 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, double NODE_XY[2*NODE_NUM], the coordinates of the points.
        //
        //    Output, double NODE_H[NODE_NUM], the value of the 
        //    H function at the points.
        //
    {
        int node;
        double[] node_h = new double[node_num];
        for (node = 0; node < node_num; node++)
        {
            node_h[node] = 1.0;
        }

        return node_h;
    }

    public static double[] k_coef(int node_num, double[] node_xy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    K_COEF evaluates the coefficient K(X,Y) of U in the Poisson equation.
        //
        //  Discussion:
        //
        //    The equation is 
        //
        //      -DEL H(X,Y) DEL U(X,Y) + K(X,Y) * U(X,Y) = F(X,Y)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    25 January 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, double NODE_XY[2*NODE_NUM], the coordinates of the points.
        //
        //    Output, double NODE_K[NODE_NUM], the value of the 
        //    K function at the points.
        //
    {
        int node;
        double[] node_k = new double[node_num];

        for (node = 0; node < node_num; node++)
        {
            node_k[node] = 0.0;
        }

        return node_k;
    }

    public static double[] rhs(int node_num, double[] node_xy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RHS gives the right-hand side of the differential equation.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    25 January 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, double NODE_XY[2*NODE_NUM], the coordinates of the points.
        //
        //    Output, double NODE_RHS[NODE_NUM], the value of the 
        //    right hand side function at the points.
        //
    {
        int node;
        double z = 0;
        double[] node_rhs = new double[node_num];

        for (node = 0; node < node_num; node++)
        {
            z = z switch
            {
                < 0.0 => 0.0,
                _ => 4.0 - Math.Sqrt(
                    Math.Pow(node_xy[0 + node * 2] - 2.0, 2) + Math.Pow(node_xy[1 + node * 2] - 2.0, 2))
            };

            node_rhs[node] = z;
        }

        return node_rhs;
    }
}