using System;
using Burkardt.MatrixNS;
using Burkardt.Table;
using Burkardt.TriangulationNS;
using Burkardt.Types;

namespace FEM2DPoissonTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for FEM2D_POISSON.
            //
            //  Discussion:
            //
            //    FEM2D_POISSON solves the Poisson equation
            //
            //      -DEL H(X,Y) DEL U(X,Y) + K(X,Y) * U(X,Y) = F(X,Y)
            //
            //    in a triangulated region in the plane.
            //
            //    Along the boundary of the region, Dirichlet conditions
            //    are imposed:
            //
            //      U(X,Y) = G(X,Y)
            //
            //    The code uses continuous piecewise linear basis functions on
            //    triangles.
            //
            //  Problem specification:
            //
            //    The user defines the geometry by supplying two data files
            //    which list the node coordinates, and list the nodes that make up
            //    each element.
            //
            //    The user specifies the right hand side of the Dirichlet boundary
            //    conditions by supplying a function
            //
            //      void dirichlet_condition ( int node_num, double node_xy[2*node_num],
            //        double node_bc[node_num] )
            //
            //    The user specifies the coefficient function K(X,Y) of the Poisson
            //    equation by supplying a routine of the form
            //
            //      void k_coef ( int node_num, double node_xy[2*node_num],
            //        double node_k[node_num] )
            //
            //    The user specifies the right hand side of the Poisson equation
            //    by supplying a routine of the form
            //
            //      void rhs ( int node_num, double node_xy[2*node_num],
            //        double node_f[node_num] )
            //
            //  Usage:
            //
            //    fem2d_poisson 'prefix'
            //
            //    where 'prefix' is the common filename prefix so that:
            //
            //    * prefix_nodes.txt contains the coordinates of the nodes;
            //    * prefix_elements.txt contains the indices of nodes forming each element.
            //
            //    Files created include:
            //
            //    * prefix_nodes.eps, an image of the nodes;
            //    * prefix_elements.eps, an image of the elements;
            //    * prefix_solution.txt, the value of the solution at every node.
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
            //  Local parameters:
            //
            //    Local, double A[(3*IB+1)*NODE_NUM], the coefficient matrix.
            //
            //    Local, double EH1, the H1 seminorm error.
            //
            //    Local, double EL2, the L2 error.
            //
            //    Local, int ELEMENT_NODE[3*ELEMENT_NUM];
            //    ELEMENT_NODE(I,J) is the global index of local node I in element J.
            //
            //    Local, int ELEMENT_NUM, the number of elements.
            //
            //    Local, integer ELEMENT_ORDER, the element order.
            //
            //    Local, double F[NODE_NUM], the right hand side.
            //
            //    Local, int IB, the half-bandwidth of the matrix.
            //
            //    Local, bool NODE_BOUNDARY[NODE_NUM], is TRUE if the node is
            //    found to lie on the boundary of the region.
            //
            //    Local, int NODE_CONDITION[NODE_NUM],
            //    indicates the condition used to determine the variable at a node.
            //    0, there is no condition (and no variable) at this node.
            //    1, a finite element equation is used;
            //    2, a Dirichlet condition is used.
            //    3, a Neumann condition is used.
            //
            //    Local, int NODE_NUM, the number of nodes.
            //
            //    Local, double NODE_R[NODE_NUM], the residual error.
            //
            //    Local, double NODE_U[NODE_NUM], the finite element coefficients.
            //
            //    Local, double NODE_XY[2*NODE_NUM], the coordinates of nodes.
            //
            //    Local, integer QUAD_NUM, the number of quadrature points used for
            //    assembly.  This is currently set to 3, the lowest reasonable value.
            //    Legal values are 1, 3, 4, 6, 7, 9, 13, and for some problems, a value
            //    of QUAD_NUM greater than 3 may be appropriate.
            //
        {
            double[] a;
            bool debug = false;
            int dim_num;
            string element_eps_filename;
            string element_filename;
            int[] element_node;
            int element_num;
            int element_order;
            int element_show;
            double[] f;
            int ib;
            int ierr;
            int job;
            int node;
            bool[] node_boundary;
            int[] node_condition;
            string node_eps_filename;
            string node_filename;
            bool node_label;
            int node_num;
            int node_show;
            double[] node_r = null;
            double[] node_u;
            double[] node_xy;
            int[] pivot;
            string prefix;
            int quad_num = 7;
            string solution_filename;
            double temp;

            Console.WriteLine("");
            Console.WriteLine("FEM2D_POISSON:");
            Console.WriteLine("");
            Console.WriteLine("  Solution of the Poisson equation in an arbitrary region");
            Console.WriteLine("  in 2 dimensions.");
            Console.WriteLine("");
            Console.WriteLine("  - DEL H(x,y) DEL U(x,y) + K(x,y) * U(x,y) = F(x,y) in the region");
            Console.WriteLine("");
            Console.WriteLine("                                     U(x,y) = G(x,y) on the boundary.");
            Console.WriteLine("");
            Console.WriteLine("  The finite element method is used,");
            Console.WriteLine("  with triangular elements,");
            Console.WriteLine("  which must be a 3 node linear triangle.");
            //
            //  Get the filename prefix.
            //
            try
            {
                prefix = args[0];
            }
            catch
            {
                Console.WriteLine("");
                Console.WriteLine("  Please enter the filename prefix:");

                prefix = Console.ReadLine();
            }

            if ((prefix != "ell") && (prefix != "lake"))
            {
                Console.WriteLine("Supported prefix value in this test is one of : baffle, ell, lake");
                return;
            }

            //
            //  Create the file names.
            //
            node_filename = prefix + "_nodes.txt";
            element_filename = prefix + "_elements.txt";
            node_eps_filename = prefix + "_nodes.eps";
            element_eps_filename = prefix + "_elements.eps";
            solution_filename = prefix + "_solution.txt";

            Console.WriteLine("");
            Console.WriteLine("  Node file is \"" + node_filename + "\".");
            Console.WriteLine("  Element file is \"" + element_filename + "\".");
            //
            //  Read the node coordinate file.
            //
            TableHeader h = typeMethods.r8mat_header_read(node_filename);

            dim_num = h.m;
            node_num = h.n;

            Console.WriteLine("  Number of nodes =          " + node_num + "");

            node_condition = new int[node_num];

            node_xy = typeMethods.r8mat_data_read(node_filename, dim_num, node_num);

            typeMethods.r8mat_transpose_print_some(dim_num, node_num, node_xy, 1, 1, 2, 10,
                "  First 10 nodes");
            //
            //  Read the element description file.
            //
            h = typeMethods.i4mat_header_read(element_filename);
            element_order = h.m;
            element_num = h.n;

            Console.WriteLine("");
            Console.WriteLine("  Element order =            " + element_order + "");
            Console.WriteLine("  Number of elements =       " + element_num + "");

            if (element_order != 3)
            {
                Console.WriteLine("");
                Console.WriteLine("FEM2D_POISSON - Fatal error!");
                Console.WriteLine("  The input triangulation has order " + element_order + "");
                Console.WriteLine("  However, a triangulation of order 3 is required.");
                return;
            }

            element_node = typeMethods.i4mat_data_read(element_filename, element_order,
                element_num);

            typeMethods.i4mat_transpose_print_some(3, element_num,
                element_node, 1, 1, 3, 10, "  First 10 elements");

            Console.WriteLine("");
            Console.WriteLine("  Quadrature order =          " + quad_num + "");
            //
            //  Determine which nodes are boundary nodes and which have a
            //  finite element unknown.  Then set the boundary values.
            //
            node_boundary = Boundary.triangulation_order3_boundary_node(node_num, element_num,
                element_node);

            if (debug)
            {
                typeMethods.lvec_print(node_num, node_boundary, "    Node  Boundary?");
            }

            //
            //  Determine the node conditions.
            //  For now, we'll just assume all boundary nodes are Dirichlet.
            //
            for (node = 0; node < node_num; node++)
            {
                if (node_boundary[node])
                {
                    node_condition[node] = 2;
                }
                else
                {
                    node_condition[node] = 1;
                }
            }

            //
            //  Determine the bandwidth of the coefficient matrix.
            //
            ib = Matrix.bandwidth(element_num, element_node);

            Console.WriteLine("");
            Console.WriteLine("  The matrix half bandwidth is " + ib + "");
            Console.WriteLine("  The matrix bandwidth is      " + 2 * ib + 1 + "");
            Console.WriteLine("  The storage bandwidth is     " + 3 * ib + 1 + "");
            //
            //  Make a picture of the nodes.
            //
            if (node_num <= 1000)
            {
                if (node_num <= 100)
                {
                    node_label = true;
                }
                else
                {
                    node_label = false;
                }

                typeMethods.points_plot(node_eps_filename, node_num, node_xy, node_label);
            }

            //
            //  Make a picture of the elements.
            //
            if (node_num <= 1000)
            {
                if (node_num <= 100)
                {
                    node_show = 2;
                    element_show = 2;
                }
                else
                {
                    node_show = 0;
                    element_show = 1;
                }

                Plot.triangulation_order3_plot(element_eps_filename, node_num,
                    node_xy, element_num, element_node, node_show, element_show);
            }

            //
            //  Allocate space for the coefficient matrix A and right hand side F.
            //
            a = new double[(3 * ib + 1) * node_num];
            f = new double[node_num];
            node_u = new double[node_num];
            pivot = new int[node_num];
            //
            //  Assemble the finite element coefficient matrix A and the right-hand side F.
            //
            switch (prefix)
            {
                case "ell":
                    Matrix.assemble_poisson(node_num, node_xy, element_num,
                    element_node, quad_num, ib, ref a, ref f, ell.rhs, ell.h_coef, ell.k_coef);
                    break;
                case "lake":
                    Matrix.assemble_poisson(node_num, node_xy, element_num,
                        element_node, quad_num, ib, ref a, ref f, lake.rhs, lake.h_coef, lake.k_coef);
                    break;
            }

            //
            //  Print a tiny portion of the matrix.
            //
            if (debug)
            {
                Matrix.dgb_print_some(node_num, node_num, ib, ib, a, 1, 1, 10, 10,
                    "  Initial block of Finite Element matrix A:");

                typeMethods.r8vec_print_some(node_num, f, 1, 10,
                    "  Part of right hand side vector:");
            }

            //
            //  Adjust the linear system to account for Dirichlet boundary conditions.
            //
            switch (prefix)
            {
                case "ell":
                    Matrix.dirichlet_apply(node_num, node_xy, node_condition, ib, ref a, ref f, ell.dirichlet_condition);
                    break;
                case "lake":
                    Matrix.dirichlet_apply(node_num, node_xy, node_condition, ib, ref a, ref f, lake.dirichlet_condition);
                    break;
            }

            if (debug)
            {
                Matrix.dgb_print_some(node_num, node_num, ib, ib, a, 1, 1, 10, 10,
                    "  Finite Element matrix A after boundary adjustments:");

                typeMethods.r8vec_print_some(node_num, f, 1, 10,
                    "  Part of right hand side vector:");
            }

            //
            //  Solve the linear system using a banded solver.
            //
            ierr = Matrix.dgb_fa(node_num, ib, ib, ref a, ref pivot);

            if (ierr != 0)
            {
                Console.WriteLine("");
                Console.WriteLine("FEM2D_POISSON - Fatal error!");
                Console.WriteLine("  DGB_FA returned the error condition IERR = " + ierr + ".");
                Console.WriteLine("");
                Console.WriteLine("  The linear system was not factored, and the");
                Console.WriteLine("  algorithm cannot proceed.");
                return;
            }

            job = 0;

            node_u = Matrix.dgb_sl(node_num, ib, ib, a, pivot, f, job);

            typeMethods.r8vec_print_some(node_num, node_u, 1, 10,
                "  Part of the solution vector U:");
            //
            //  Write an ASCII file that can be read into MATLAB.
            //
            typeMethods.r8mat_write(solution_filename, 1, node_num, node_u);

            Console.WriteLine("");
            Console.WriteLine("FEM2D_POISSON:");
            Console.WriteLine("  Wrote an ASCII file");
            Console.WriteLine("    \"" + solution_filename + "\".");
            Console.WriteLine("  of the form");
            Console.WriteLine("    U ( X(I), Y(I) )");
            Console.WriteLine("  which can be used for plotting.");

            if (debug)
            {
                typeMethods.r8vec_print_some(node_num, node_u, 1, 10,
                    "  Part of the solution vector:");
            }

            switch (prefix)
            {
                case "ell":
                    node_r = Matrix.residual_poisson(node_num, node_xy, node_condition,
                        element_num, element_node, quad_num, ib, a, f, node_u , ell.rhs, ell.h_coef, ell.k_coef, ell.dirichlet_condition);
                    break;
                case "lake":
                    node_r = Matrix.residual_poisson(node_num, node_xy, node_condition,
                        element_num, element_node, quad_num, ib, a, f, node_u , lake.rhs, lake.h_coef, lake.k_coef, lake.dirichlet_condition);
                    break;
            }

            temp = typeMethods.r8vec_amax(node_num, node_r);

            Console.WriteLine("");
            Console.WriteLine("  Maximum absolute residual = " + temp + "");

            Console.WriteLine("");
            Console.WriteLine("FEM2D_POISSON:");
            Console.WriteLine("  Normal end of execution.");

            Console.WriteLine("");
        }
    }
}