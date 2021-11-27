using System;
using Burkardt.MatrixNS;
using Burkardt.SolveNS;
using Burkardt.Table;
using Burkardt.TriangulationNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace FEM2DPoissonSparseTest;

internal static class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for FEM2D_POISSON_SPARSE.
        //
        //  Discussion:
        //
        //    This program uses a sparse matrix storage format and an iterative solver,
        //    which allow it to solve larger problems faster.
        //
        //    This program solves the Poisson equation
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
        //    The user specifies the coefficient function H(X,Y) of the Poisson
        //    equation by supplying a routine of the form
        //
        //      void h_coef ( int node_num, double node_xy[2*node_num],
        //        double node_h[node_num] )
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
        //    fem2d_poisson_sparse prefix
        //
        //    where 'prefix' is the common filename prefix so that:
        //
        //    * prefix_nodes.txt contains the coordinates of the nodes;
        //    * prefix_elements.txt contains the indices of nodes forming each element.
        //
        //    Files created include:
        //
        //    * prefix_values.txt, the value of the solution at every node.
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
        //  Local parameters:
        //
        //    Local, double A[NZ_NUM], the coefficient matrix.
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
        //    Local, int IA[NZ_NUM], the row indices of the nonzero entries
        //    of the coefficient matrix.
        //
        //    Local, int JA[NZ_NUM], the column indices of the nonzero entries
        //    of the coefficient matrix.
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
        //    Local, double NODE_U[NODE_NUM], the finite element coefficients.
        //
        //    Local, double NODE_XY[2*NODE_NUM], the coordinates of nodes.
        //
        //    Local, int NZ_NUM, the number of nonzero entries
        //    in the coefficient matrix.
        //
        //    Local, integer QUAD_NUM, the number of quadrature points used for
        //    assembly.  This is currently set to 3, the lowest reasonable value.
        //    Legal values are 1, 3, 4, 6, 7, 9, 13, and for some problems, a value
        //    of QUAD_NUM greater than 3 may be appropriate.
        //
    {
        const bool debug = false;
        int node;
        double[] node_u = null;
        string prefix;
        const int quad_num = 3;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("FEM2D_POISSON_SPARSE:");
        Console.WriteLine("");
        Console.WriteLine("  A finite element method solver for the Poisson problem");
        Console.WriteLine("  in an arbitrary triangulated region in 2 dimensions,");
        Console.WriteLine("  using sparse storage and an iterative solver.");
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

        if (prefix != "baffle" && prefix != "ell" && prefix != "lake")
        {
            Console.WriteLine("Supported prefix value in this test is one of : baffle, ell, lake");
            return;
        }

        //
        //  Create the file names.
        //
        string node_filename = prefix + "_nodes.txt";
        string element_filename = prefix + "_elements.txt";
        string solution_filename = prefix + "_values.txt";

        Console.WriteLine("");
        Console.WriteLine("  Node file is \"" + node_filename + "\".");
        Console.WriteLine("  Element file is \"" + element_filename + "\".");
        //
        //  Read the node coordinate file.
        //
        TableHeader h = typeMethods.r8mat_header_read(node_filename);
        int dim_num = h.m;
        int node_num = h.n;

        Console.WriteLine("  Number of nodes =          " + node_num + "");

        int[] node_condition = new int[node_num];

        double[] node_xy = typeMethods.r8mat_data_read(node_filename, dim_num, node_num);

        typeMethods.r8mat_transpose_print_some(dim_num, node_num, node_xy, 1, 1, 2, 10,
            "  First 10 nodes");
        //
        //  Read the triangle description file.
        //
        h = typeMethods.i4mat_header_read(element_filename);
        int element_order = h.m;
        int element_num = h.n;

        Console.WriteLine("");
        Console.WriteLine("  Element order =            " + element_order + "");
        Console.WriteLine("  Number of elements =       " + element_num + "");

        if (element_order != 3)
        {
            Console.WriteLine("");
            Console.WriteLine("FEM2D_POISSON_SPARSE - Fatal error!");
            Console.WriteLine("  The input triangulation has order " + element_order + "");
            Console.WriteLine("  However, a triangulation of order 3 is required.");
            return;
        }

        int[] element_node = typeMethods.i4mat_data_read(element_filename, element_order,
            element_num);

        typeMethods.i4mat_transpose_print_some(3, element_num,
            element_node, 1, 1, 3, 10, "  First 10 elements");

        Console.WriteLine("");
        Console.WriteLine("  Quadrature order =          " + quad_num + "");
        //
        //  Determine which nodes are boundary nodes and which have a
        //  finite element unknown.  Then set the boundary values.
        //
        bool[] node_boundary = Boundary.triangulation_order3_boundary_node(node_num, element_num,
            element_node);
        //
        //  Determine the node conditions.
        //  For now, we'll just assume all boundary nodes are Dirichlet.
        //
        for (node = 0; node < node_num; node++)
        {
            node_condition[node] = node_boundary[node] switch
            {
                true => 2,
                _ => 1
            };
        }

        //
        //  Determine the element neighbor array, just so we can estimate
        //  the nonzeros.
        //
        int[] element_neighbor = new int[3 * element_num];

        element_neighbor = Neighbor.triangulation_order3_neighbor_triangles(element_num, element_node);
        //
        //  Count the number of nonzeros.
        //
        int[] adj_col = new int[node_num + 1];

        int nz_num = Adjacency.triangulation_order3_adj_count(node_num, element_num,
            element_node, element_neighbor, adj_col);

        Console.WriteLine("");
        Console.WriteLine("  Number of nonzero coefficients NZ_NUM = " + nz_num + "");
        //
        //  Set up the sparse row and column index vectors.
        //
        int[] ia = new int[nz_num];
        int[] ja = new int[nz_num];

        Adjacency.triangulation_order3_adj_set2(node_num, element_num, element_node,
            element_neighbor, nz_num, adj_col, ia, ja);

        //
        //  Allocate space for the coefficient matrix A and right hand side F.
        //
        double[] a = new double[nz_num];
        double[] f = new double[node_num];
        //
        //  Assemble the finite element coefficient matrix A and the right-hand side F.
        //
        switch (prefix)
        {
            case "baffle":
                DSP.assemble_poisson_dsp(node_num, node_xy, element_num,
                    element_node, quad_num, nz_num, ia, ja, ref a, ref f, baffle.rhs, baffle.h_coef, baffle.k_coef);
                break;
            case "ell":
                DSP.assemble_poisson_dsp(node_num, node_xy, element_num,
                    element_node, quad_num, nz_num, ia, ja, ref a, ref f, ell.rhs, ell.h_coef, ell.k_coef);
                break;
            case "lake":
                DSP.assemble_poisson_dsp(node_num, node_xy, element_num,
                    element_node, quad_num, nz_num, ia, ja, ref a, ref f, lake.rhs, lake.h_coef, lake.k_coef);
                break;
        }

        switch (debug)
        {
            //
            //  Print a portion of the matrix.
            //
            case true:
                DSP.dsp_print_some(node_num, node_num, nz_num, ia, ja, a, 1, 1, 10, 10,
                    "  Part of Finite Element matrix A:");

                typeMethods.r8vec_print_some(node_num, f, 1, 10,
                    "  Part of right hand side vector F:");
                break;
        }

        //
        //  Adjust the linear system to account for Dirichlet boundary conditions.
        //
        switch (prefix)
        {
            case "baffle":
                DSP.dirichlet_apply_dsp(node_num, node_xy, node_condition, nz_num, ia, ja,
                    ref a, ref f, baffle.dirichlet_condition);
                break;
            case "ell":
                DSP.dirichlet_apply_dsp(node_num, node_xy, node_condition, nz_num, ia, ja,
                    ref a, ref f, ell.dirichlet_condition);
                break;
            case "lake":
                DSP.dirichlet_apply_dsp(node_num, node_xy, node_condition, nz_num, ia, ja,
                    ref a, ref f, lake.dirichlet_condition);
                break;
        }

        switch (debug)
        {
            case true:
                DSP.dsp_print_some(node_num, node_num, nz_num, ia, ja, a, 1, 1, 10, 10,
                    "  Part of finite Element matrix A after boundary adjustments:");

                typeMethods.r8vec_print_some(node_num, f, 1, 10,
                    "  Part of right hand side vector F:");
                break;
        }

        //
        //  Solve the linear system using an iterative solver.
        //
        UniformRNG.r8vec_uniform_01(node_num, ref seed, ref node_u);

        int itr_max = 20;
        int mr = 20;
        double tol_abs = 0.000001;
        double tol_rel = 0.000001;

        RestartedGeneralizedMinimumResidual.mgmres(a, ia, ja, ref node_u, f, node_num, nz_num, itr_max, mr, tol_abs,
            tol_rel);

        typeMethods.r8vec_print_some(node_num, node_u, 1, 10,
            "  Part of the solution vector vector U:");
        //
        //  Write an ASCII file that can be read into MATLAB.
        //
        typeMethods.r8mat_write(solution_filename, 1, node_num, node_u);

        Console.WriteLine("");
        Console.WriteLine("FEM2D_POISSON_SPARSE:");
        Console.WriteLine("  Wrote an ASCII file");
        Console.WriteLine("    \"" + solution_filename + "\".");
        Console.WriteLine("  of the form");
        Console.WriteLine("    U ( X(I), Y(I) )");
        Console.WriteLine("  which can be used for plotting.");

        switch (debug)
        {
            case true:
                typeMethods.r8vec_print_some(node_num, node_u, 1, 10,
                    "  Part of the solution vector:");
                break;
        }

        Console.WriteLine("");
        Console.WriteLine("FEM2D_POISSON_SPARSE:");
        Console.WriteLine("  Normal end of execution.");

        Console.WriteLine("");
    }
}