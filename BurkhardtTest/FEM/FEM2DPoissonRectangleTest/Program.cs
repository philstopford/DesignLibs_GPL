using System;
using Burkardt;
using Burkardt.FEM;
using Burkardt.MatrixNS;
using Burkardt.Quadrature;
using Burkardt.TriangulationNS;
using Burkardt.Types;

namespace FEM2DPoissonRectangleTest;

internal class Program
{
    private static void Main(string[] args)
        /****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for FEM2D_POISSON_RECTANGLE.
        //
        //  Discussion:
        //
        //    FEM2D_POISSON_RECTANGLE solves
        //
        //      -Laplacian U(X,Y) = F(X,Y)
        //
        //    in a rectangular region in the plane.  Along the boundary,
        //    Dirichlet boundary conditions are imposed.
        //
        //      U(X,Y) = G(X,Y)
        //
        //    The code uses continuous piecewise quadratic basis functions on
        //    triangles determined by a uniform grid of NX by NY points.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 September 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Local parameters:
        //
        //    Local, double A[(3*IB+1)*NUNK], the coefficient matrix.
        //
        //    Local, double ELEMENT_AREA[ELEMENT_NUM], the area of each element.
        //
        //    Local, double C[NUNK], the finite element coefficients, solution of A * C = F.
        //
        //    Local, double EH1, the H1 seminorm error.
        //
        //    Local, double EL2, the L2 error.
        //
        //    Local, int ELEMENT_NODE[ELEMENT_NUM*NNODES]; ELEMENT_NODE(I,J) is the
        //    global node index of the local node J in element I.
        //
        //    Local, int ELEMENT_NUM, the number of elements.
        //
        //    Local, double F[NUNK], the right hand side.
        //
        //    Local, int IB, the half-bandwidth of the matrix.
        //
        //    Local, int INDX[NODE_NUM], gives the index of the unknown quantity
        //    associated with the given node.
        //
        //    Local, int NNODES, the number of nodes used to form one element.
        //
        //    Local, double NODE_XY[2*NODE_NUM], the X and Y coordinates of nodes.
        //
        //    Local, int NQ, the number of quadrature points used for assembly.
        //
        //    Local, int NUNK, the number of unknowns.
        //
        //    Local, int NX, the number of points in the X direction.
        //
        //    Local, int NY, the number of points in the Y direction.
        //
        //    Local, double WQ[NQ], quadrature weights.
        //
        //    Local, double XL, XR, YB, YT, the X coordinates of
        //    the left and right sides of the rectangle, and the Y coordinates
        //    of the bottom and top of the rectangle.
        //
        //    Local, double XQ[NQ*ELEMENT_NUM], YQ[NQ*ELEMENT_NUM], the X and Y
        //    coordinates of the quadrature points in each element.
        */
    {
        int NNODES = 6;
        int NQ = 3;
        int NX = 7;
        int NY = 7;
        int ELEMENT_NUM = (NX - 1) * (NY - 1) * 2;
        int NODE_NUM = (2 * NX - 1) * (2 * NY - 1);

        double[] a;
        double[] c;
        double eh1 = 0;
        double el2 = 0;
        double[] element_area = new double[ELEMENT_NUM];
        int[] element_node = new int[NNODES * ELEMENT_NUM];
        double[] f;
        int ib;
        int ierr;
        int[] indx = new int[NODE_NUM];
        int job;
        string node_eps_file_name = "fem2d_poisson_rectangle_nodes.eps";
        string node_txt_file_name = "fem2d_poisson_rectangle_nodes.txt";
        bool node_label;
        int node_show;
        double[] node_xy = new double[2 * NODE_NUM];
        int nunk = 0;
        int[] pivot;
        string solution_txt_file_name = "fem2d_poisson_rectangle_solution.txt";
        int triangle_show;
        string triangulation_eps_file_name = "fem2d_poisson_rectangle_elements.eps";
        string triangulation_txt_file_name = "fem2d_poisson_rectangle_elements.txt";
        double[] wq = new double[NQ];
        double xl = 0.0E+00;
        double[] xq = new double[NQ * ELEMENT_NUM];
        double xr = 1.0E+00;
        double yb = 0.0E+00;
        double[] yq = new double[NQ * ELEMENT_NUM];
        double yt = 1.0E+00;

        Console.WriteLine("");
        Console.WriteLine("FEM2D_POISSON_RECTANGLE:");
        Console.WriteLine("");
        Console.WriteLine("  Solution of the Poisson equation on a unit box");
        Console.WriteLine("  in 2 dimensions.");
        Console.WriteLine("");
        Console.WriteLine("  - Uxx - Uyy = F(x,y) in the box");
        Console.WriteLine("       U(x,y) = G(x,y) on the boundary.");
        Console.WriteLine("");
        Console.WriteLine("  The finite element method is used, with piecewise");
        Console.WriteLine("  quadratic basis functions on 6 node triangular");
        Console.WriteLine("  elements.");
        Console.WriteLine("");
        Console.WriteLine("  The corner nodes of the triangles are generated by an");
        Console.WriteLine("  underlying grid whose dimensions are");
        Console.WriteLine("");
        Console.WriteLine("  NX =                 " + NX + "");
        Console.WriteLine("  NY =                 " + NY + "");
        Console.WriteLine("");
        Console.WriteLine("  Number of nodes    = " + NODE_NUM + "");
        Console.WriteLine("  Number of elements = " + ELEMENT_NUM + "");
        //
        //  Set the coordinates of the nodes.
        //
        XY.xy_set(NX, NY, NODE_NUM, xl, xr, yb, yt, ref node_xy);
        //
        //  Organize the nodes into a grid of 6-node triangles.
        //
        Grid.grid_t6(NX, NY, NNODES, ELEMENT_NUM, ref element_node);
        //
        //  Set the quadrature rule for assembly.
        //
        QuadratureRule.quad_a(node_xy, element_node, ELEMENT_NUM, NODE_NUM,
            NNODES, ref wq, ref xq, ref yq);
        //
        //  Determine the areas of the elements.
        //
        Element.area_set(NODE_NUM, node_xy, NNODES, ELEMENT_NUM,
            element_node, element_area);
        //
        //  Determine which nodes are boundary nodes and which have a
        //  finite element unknown.  Then set the boundary values.
        //
        Burkardt.FEM.Boundary.indx_set(NX, NY, NODE_NUM, ref indx, ref nunk);

        Console.WriteLine("  Number of unknowns =       " + nunk + "");
        //
        //  Determine the bandwidth of the coefficient matrix.
        //
        ib = Matrix.bandwidth(NNODES, ELEMENT_NUM, element_node, NODE_NUM, indx);

        Console.WriteLine("");
        Console.WriteLine("  Total bandwidth is " + 3 * ib + 1 + "");
        switch (NX)
        {
            //
            //  Make an EPS picture of the nodes.
            //
            case <= 10 when NY <= 10:
                node_label = true;
                typeMethods.nodes_plot(node_eps_file_name, NODE_NUM, node_xy, node_label);

                Console.WriteLine("");
                Console.WriteLine("FEM2D_POISSON_RECTANGLE:");
                Console.WriteLine("  Wrote an EPS file");
                Console.WriteLine("    \"" + node_eps_file_name + "\".");
                Console.WriteLine("  containing a picture of the nodes.");
                break;
        }

        //
        //  Write the nodes to an ASCII file that can be read into MATLAB.
        //
        typeMethods.nodes_write(NODE_NUM, node_xy, node_txt_file_name);

        Console.WriteLine("");
        Console.WriteLine("FEM2D_POISSON_RECTANGLE:");
        Console.WriteLine("  Wrote an ASCII node file");
        Console.WriteLine("    " + node_txt_file_name + "");
        Console.WriteLine("  of the form");
        Console.WriteLine("    X(I), Y(I)");
        Console.WriteLine("  which can be used for plotting.");
        switch (NX)
        {
            //
            //  Make a picture of the elements.
            //
            case <= 10 when NY <= 10:
                node_show = 1;
                triangle_show = 2;

                Plot.triangulation_order6_plot(triangulation_eps_file_name, NODE_NUM,
                    node_xy, ELEMENT_NUM, element_node, node_show, triangle_show);

                Console.WriteLine("");
                Console.WriteLine("FEM2D_POISSON_RECTANGLE:");
                Console.WriteLine("  Wrote an EPS file");
                Console.WriteLine("    \"" + triangulation_eps_file_name + "\".");
                Console.WriteLine("  containing a picture of the elements.");
                break;
        }

        //
        //  Write the elements to a file that can be read into MATLAB.
        //
        Element.element_write(NNODES, ELEMENT_NUM, element_node,
            triangulation_txt_file_name);

        Console.WriteLine("");
        Console.WriteLine("FEM2D_POISSON_RECTANGLE:");
        Console.WriteLine("  Wrote an ASCII element file");
        Console.WriteLine("    \"" + triangulation_txt_file_name + "\".");
        Console.WriteLine("  of the form");
        Console.WriteLine("    Node(1) Node(2) Node(3) Node(4) Node(5) Node(6)");
        Console.WriteLine("  which can be used for plotting.");
        //
        //  Allocate space for the coefficient matrix A and right hand side F.
        //
        a = new double[(3 * ib + 1) * nunk];
        f = new double[nunk];
        pivot = new int[nunk];
        //
        //  Assemble the coefficient matrix A and the right-hand side F of the
        //  finite element equations.
        //
        Matrix.assemble(NODE_NUM, node_xy, NNODES,
            ELEMENT_NUM, element_node, NQ,
            wq, xq, yq, element_area, indx, ib, nunk, ref a, ref f);
        //
        //  Print a tiny portion of the matrix.
        //
        Matrix.dgb_print_some(nunk, nunk, ib, ib, a, 1, 1, 5, 5,
            "  Initial 5 x 5 block of coefficient matrix A:");

        typeMethods.r8vec_print_some(nunk, f, 10, "  Part of the right hand side F:");
        //
        //  Modify the coefficient matrix and right hand side to account for
        //  boundary conditions.
        //
        Burkardt.FEM.Boundary.boundary(NX, NY, NODE_NUM, node_xy, indx, ib, nunk, ref a, ref f, exact);
        //
        //  Print a tiny portion of the matrix.
        //
        Matrix.dgb_print_some(nunk, nunk, ib, ib, a, 1, 1, 5, 5,
            "  A after boundary adjustment:");

        typeMethods.r8vec_print_some(nunk, f, 10, "  F after boundary adjustment:");
        //
        //  Solve the linear system using a banded solver.
        //
        ierr = Matrix.dgb_fa(nunk, ib, ib, ref a, ref pivot);

        if (ierr != 0)
        {
            Console.WriteLine("");
            Console.WriteLine("FEM2D_POISSON_RECTANGLE - Error!");
            Console.WriteLine("  DGB_FA returned an error condition.");
            Console.WriteLine("");
            Console.WriteLine("  The linear system was not factored, and the");
            Console.WriteLine("  algorithm cannot proceed.");
            return;
        }

        job = 0;
        c = Matrix.dgb_sl(nunk, ib, ib, a, pivot, f, job);

        typeMethods.r8vec_print_some(nunk, c, 10, "  Part of the solution vector:");
        //
        //  Calculate error using 13 point quadrature rule.
        //
        QuadratureRule.errors(element_area, element_node, indx, node_xy, c,
            ELEMENT_NUM, NNODES, nunk, NODE_NUM, ref el2, ref eh1, exact);
        //
        //  Compare the exact and computed solutions just at the nodes.
        //
        typeMethods.compare(NODE_NUM, node_xy, indx, nunk, c, exact);
        //
        //  Write an ASCII file that can be read into MATLAB.
        //
        typeMethods.solution_write(c, indx, NODE_NUM, nunk, solution_txt_file_name,
            node_xy, exact);

        Console.WriteLine("");
        Console.WriteLine("FEM2D_POISSON_RECTANGLE:");
        Console.WriteLine("  Wrote an ASCII solution file");
        Console.WriteLine("    " + solution_txt_file_name + "");
        Console.WriteLine("  of the form");
        Console.WriteLine("    U( X(I), Y(I) )");
        Console.WriteLine("  which can be used for plotting.");

        Console.WriteLine("");
        Console.WriteLine("FEM2D_POISSON_RECTANGLE:");
        Console.WriteLine("  Normal end of execution.");

        Console.WriteLine("");
    }

    private static ExactResult exact(double x, double y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT calculates the exact solution and its first derivatives.
        //
        //  Discussion:
        //
        //    The function specified here depends on the problem being
        //    solved.  This is one of the routines that a user will
        //    normally want to change.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 September 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, Y, the coordinates of a point
        //    in the region, at which the right hand side of the
        //    differential equation is to be evaluated.
        //
        //    Output, double *U, *DUDX, *DUDY, the value of
        //    the exact solution U and its derivatives dUdX
        //    and dUdY at the point (X,Y).
        //
    {

        ExactResult res = new()
        {
            u = Math.Sin(Math.PI * x) * Math.Sin(Math.PI * y) + x,
            dudx = Math.PI * Math.Cos(Math.PI * x) * Math.Sin(Math.PI * y) + 1.0E+00,
            dudy = Math.PI * Math.Sin(Math.PI * x) * Math.Cos(Math.PI * y)
        };

        return res;
    }
}