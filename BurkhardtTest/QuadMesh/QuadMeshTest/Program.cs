using System;
using System.Globalization;
using Burkardt.QuadMesh;
using Burkardt.Types;

namespace QuadMeshTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for QUAD_MESH_TEST.
        //
        //  Discussion:
        //
        //    QUAD_MESH_TEST tests the QUAD_MESH library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 March 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {

        Console.WriteLine("");
        Console.WriteLine("QUAD_MESH_TEST");
        Console.WriteLine("  Test the QUAD_MESH library.");

        adj_size_q4_mesh_test();
        area_q4_mesh_test();
        area_quad_test();
        boundary_edge_count_q4_mesh_test();
        boundary_edge_count_euler_q4_mesh_test();
        example1_q4_mesh_test();
        example2_q4_mesh_test();
        test08();
        test09();
        test10();
        test105();
        sample_quad_test();
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("QUAD_MESH_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void adj_size_q4_mesh_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ADJ_SIZE_Q4_MESH_TEST tests ADJ_SIZE_Q4_MESH.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 February 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int element_num = 0;
        int hole_num = 0;
        int node;
        int node_num = 0;

        Console.WriteLine("");
        Console.WriteLine("ADJ_SIZE_Q4_MESH_TEST");
        Console.WriteLine("  ADJ_SIZE_Q4_MEXH counts the node adjacencies.");
        //
        //  Get the sizes of the example.
        //
        Burkardt.Values.QuadMesh.example1_q4_mesh_size(ref node_num, ref element_num, ref hole_num);

        int[] element_neighbor = new int[4 * element_num];
        int[] element_node = new int[4 * element_num];
        double[] node_xy = new double[2 * node_num];

        Burkardt.Values.QuadMesh.example1_q4_mesh(node_num, element_num, ref node_xy, ref element_node,
            ref element_neighbor);
        //
        //  Get the count of the node adjacencies.
        //
        int[] adj_col = new int[node_num + 1];

        int adj_num = Adjacency.adj_size_q4_mesh(node_num, element_num, element_node,
            element_neighbor, ref adj_col);

        Console.WriteLine("");
        Console.WriteLine("  Number of adjacency entries is " + adj_num + "");

        Console.WriteLine("");
        Console.WriteLine("  Adjacency pointers:");
        Console.WriteLine("");
        for (node = 0; node < node_num; node++)
        {
            Console.WriteLine("  " + node.ToString().PadLeft(8)
                                   + "  " + adj_col[node].ToString().PadLeft(8)
                                   + "  " + (adj_col[node + 1] - 1).ToString().PadLeft(8) + "");
        }
    }

    private static void area_q4_mesh_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    AREA_Q4_MESH_TEST tests AREA_Q4_MESH.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 February 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int element_num = 0;
        int hole_num = 0;
        double mesh_area = 0;
        int node_num = 0;

        Console.WriteLine("");
        Console.WriteLine("AREA_Q4_MESH_TEST");
        Console.WriteLine("  AREA_Q4_MESH computes the area of each element");
        Console.WriteLine("  in a Q4 mesh.");

        Burkardt.Values.QuadMesh.example1_q4_mesh_size(ref node_num, ref element_num, ref hole_num);

        int[] element_neighbor = new int[4 * element_num];
        int[] element_node = new int[4 * element_num];
        double[] node_xy = new double[2 * node_num];

        Burkardt.Values.QuadMesh.example1_q4_mesh(node_num, element_num, ref node_xy, ref element_node,
            ref element_neighbor);

        double[] element_area = new double[element_num];

        Area.area_q4_mesh(node_num, element_num, node_xy, element_node,
            element_area, ref mesh_area);

        typeMethods.r8vec_print(element_num, element_area, "  Element areas:");

        Console.WriteLine("");
        Console.WriteLine("   Mesh =   " + mesh_area + "");
    }

    private static void area_quad_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    AREA_QUAD_TEST demonstrates AREA_QUAD.
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
    {
        double[] quad_xy =
        {
            1.0, 2.0,
            5.0, 2.0,
            5.0, 3.0,
            4.0, 4.0
        };

        Console.WriteLine("");
        Console.WriteLine("AREA_QUAD_TEST");
        Console.WriteLine("  AREA_QUAD computes the area of a quadrilateral.");

        double area = Area.area_quad(quad_xy);

        Console.WriteLine("");
        Console.WriteLine("  Area = " + area + "");

    }

    private static void boundary_edge_count_q4_mesh_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BOUNDARY_EDGE_COUNT_Q4_MESH_TEST tests BOUNDARY_EDGE_COUNT_Q4_MESH.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 February 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int element_num = 0;
        int hole_num = 0;
        int node_num = 0;

        Console.WriteLine("");
        Console.WriteLine("BOUNDARY_EDGE_COUNT_Q4_MESH_TEST");
        Console.WriteLine("  BOUNDARY_EDGE_COUNT_Q4_MESH counts the");
        Console.WriteLine("  boundary edges.");

        Burkardt.Values.QuadMesh.example1_q4_mesh_size(ref node_num, ref element_num, ref hole_num);

        int[] element_neighbor = new int[4 * element_num];
        int[] element_node = new int[4 * element_num];
        double[] node_xy = new double[2 * node_num];

        Burkardt.Values.QuadMesh.example1_q4_mesh(node_num, element_num, ref node_xy, ref element_node,
            ref element_neighbor);

        int boundary_edge_num = Boundary.boundary_edge_count_q4_mesh(element_num, element_node);

        Console.WriteLine("");
        Console.WriteLine("  Number of boundary edges = " + boundary_edge_num + "");
        Console.WriteLine("  Correct number =           " + 22 + "");

    }

    private static void boundary_edge_count_euler_q4_mesh_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BOUNDARY_EDGE_COUNT_EULER_Q4_MESH_TEST tests BOUNDARY_EDGE_COUNT_EULER_Q4_MESH.
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
    {
        int boundary_edge_num = 0;
        int element_num = 0;
        int hole_num = 0;
        int node_num = 0;

        Console.WriteLine("");
        Console.WriteLine("BOUNDARY_EDGE_COUNT_EULER_Q4_MESH_TEST");
        Console.WriteLine("  BOUNDARY_EDGE_COUNT_EULER_Q4_MESH counts the");
        Console.WriteLine("  boundary edges using Euler's formula.");

        Burkardt.Values.QuadMesh.example1_q4_mesh_size(ref node_num, ref element_num, ref hole_num);

        boundary_edge_num = Boundary.boundary_edge_count_euler_q4_mesh(node_num,
            element_num, hole_num);

        Console.WriteLine("");
        Console.WriteLine("  Number of boundary edges = " + boundary_edge_num + "");
        Console.WriteLine("  Correct number =           " + 22 + "");
    }

    private static void example1_q4_mesh_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXAMPLE1_Q4_MESH_TEST tests EXAMPLE1_Q4_MESH.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 December 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int element_num = 0;
        int hole_num = 0;
        int node_num = 0;

        Console.WriteLine("");
        Console.WriteLine("EXAMPLE1_Q4_MESH_TEST");
        Console.WriteLine("  EXAMPLE1_Q4_MESH sets up example #1 Q4 mesh.");

        Burkardt.Values.QuadMesh.example1_q4_mesh_size(ref node_num, ref element_num, ref hole_num);

        Console.WriteLine("");
        Console.WriteLine("  Number of nodes =    " + node_num + "");
        Console.WriteLine("  Number of elements = " + element_num + "");
        Console.WriteLine("  Number of holes =    " + hole_num + "");

        int[] element_neighbor = new int[4 * element_num];
        int[] element_node = new int[4 * element_num];
        double[] node_xy = new double[2 * node_num];

        Burkardt.Values.QuadMesh.example1_q4_mesh(node_num, element_num, ref node_xy, ref element_node,
            ref element_neighbor);
        //
        //  Print the mesh.
        //
        typeMethods.r8mat_transpose_print(2, node_num, node_xy, "  Node coordinates:");
        typeMethods.i4mat_transpose_print(4, element_num, element_node, "  Elements:");
        typeMethods.i4mat_transpose_print(4, element_num, element_neighbor,
            "  Element neighbors");
        //
        //  Plot the mesh.
        //
        const int node_show = 2;
        const int element_show = 2;
        const string output_filename = "q4_mesh_ex1.eps";

        Plot.plot_q4_mesh(node_num, element_num, node_xy, element_node,
            node_show, element_show, output_filename);

    }

    private static void example2_q4_mesh_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXAMPLE2_Q4_MESH_TEST tests EXAMPLE2_Q4_MESH.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 December 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int element_num = 0;
        int hole_num = 0;
        int node_num = 0;

        Console.WriteLine("");
        Console.WriteLine("EXAMPLE2_Q4_MESH_TEST");
        Console.WriteLine("  EXAMPLE2_Q4_MESH sets up example #2 Q4 mesh.");

        Burkardt.Values.QuadMesh.example2_q4_mesh_size(ref node_num, ref element_num, ref hole_num);

        Console.WriteLine("");
        Console.WriteLine("  Number of nodes =    " + node_num + "");
        Console.WriteLine("  Number of elements = " + element_num + "");
        Console.WriteLine("  Number of holes =    " + hole_num + "");

        int[] element_neighbor = new int[4 * element_num];
        int[] element_node = new int[4 * element_num];
        double[] node_xy = new double[2 * node_num];

        Burkardt.Values.QuadMesh.example2_q4_mesh(node_num, element_num, ref node_xy, ref element_node,
            ref element_neighbor);
        //
        //  Print the mesh.
        //
        typeMethods.r8mat_transpose_print(2, node_num, node_xy, "  Node coordinates:");
        typeMethods.i4mat_transpose_print(4, element_num, element_node, "  Elements:");
        typeMethods.i4mat_transpose_print(4, element_num, element_neighbor,
            "  Element neighbors");
        //
        //  Plot the mesh.
        //
        const int node_show = 2;
        const int element_show = 2;
        const string output_filename = "q4_mesh_ex2.eps";

        Plot.plot_q4_mesh(node_num, element_num, node_xy, element_node,
            node_show, element_show, output_filename);
    }

    private static void test08()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST08 tests NEIGHBOR_ELEMENTS_Q4_MESH.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 February 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int element_num = 0;
        int hole_num = 0;
        int node_num = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST08");
        Console.WriteLine("  NEIGHBOR_ELEMENTS_Q4_MESH determines the");
        Console.WriteLine("  adjacency relationships between elements.");

        Burkardt.Values.QuadMesh.example1_q4_mesh_size(ref node_num, ref element_num, ref hole_num);

        int[] element_neighbor = new int[4 * element_num];
        int[] element_node = new int[4 * element_num];
        double[] node_xy = new double[2 * node_num];

        Burkardt.Values.QuadMesh.example1_q4_mesh(node_num, element_num, ref node_xy, ref element_node,
            ref element_neighbor);

        typeMethods.i4mat_transpose_print(4, element_num, element_neighbor,
            "  Element neighbors as reported by EXAMPLE1_Q4_MESH:");

        int[] element_neighbor2 = Neighbors.neighbor_elements_q4_mesh(element_num, element_node);

        typeMethods.i4mat_transpose_print(4, element_num, element_neighbor2,
            "  Element neighbors computed by NEIGHBOR_ELEMENTS_Q4_MESH:");
    }

    private static void test09()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST09 writes data to files.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 February 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int element_num = 0;
        int hole_num = 0;
        int node_num = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST09");
        Console.WriteLine("  Write Q4 Mesh Example #1 to files.");

        Burkardt.Values.QuadMesh.example2_q4_mesh_size(ref node_num, ref element_num, ref hole_num);

        int[] element_neighbor = new int[4 * element_num];
        int[] element_node = new int[4 * element_num];
        double[] node_xy = new double[2 * node_num];

        Burkardt.Values.QuadMesh.example2_q4_mesh(node_num, element_num, ref node_xy, ref element_node,
            ref element_neighbor);

        string output_filename = "q4_mesh_ex2_element_neighbors.txt";
        typeMethods.i4mat_write(output_filename, 4, element_num, element_neighbor);
        Console.WriteLine("");
        Console.WriteLine("  Element neighbors written to \"" + output_filename + "\".");

        output_filename = "q4_mesh_ex2_elements.txt";
        typeMethods.i4mat_write(output_filename, 4, element_num, element_node);
        Console.WriteLine("  Elements written to \"" + output_filename + "\".");

        output_filename = "q4_mesh_ex2_xy.txt";
        typeMethods.r8mat_write(output_filename, 2, node_num, node_xy);
        Console.WriteLine("  Node coordinates written to \"" + output_filename + "\".");
    }

    private static void test10()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST10 tests NODE_ORDER_Q4_MESH.
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
    {
        int element_num = 0;
        int hole_num = 0;
        int node_num = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST10");
        Console.WriteLine("  NODE_ORDER_4_MESH computes the order");
        Console.WriteLine("  of the nodes in a Q4 mesh.");

        Burkardt.Values.QuadMesh.example1_q4_mesh_size(ref node_num, ref element_num, ref hole_num);

        int[] element_neighbor = new int[4 * element_num];
        int[] element_node = new int[4 * element_num];
        double[] node_xy = new double[2 * node_num];

        Burkardt.Values.QuadMesh.example1_q4_mesh(node_num, element_num, ref node_xy, ref element_node,
            ref element_neighbor);

        int[] node_order = NodeOrder.node_order_q4_mesh(element_num, element_node, node_num);

        typeMethods.i4vec_print(node_num, node_order, "      NODE         ORDER");
    }

    private static void test105()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST105 tests SAMPLE_Q4_MESH.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 March 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int element_num = 0;
        int hole_num = 0;
        int node_num = 0;
        int sample;

        Console.WriteLine("");
        Console.WriteLine("TEST105");
        Console.WriteLine("  SAMPLE_Q4_MESH returns uniform sample points from");
        Console.WriteLine("  a Q4 mesh.");

        Burkardt.Values.QuadMesh.example1_q4_mesh_size(ref node_num, ref element_num, ref hole_num);

        int[] element_neighbor = new int[4 * element_num];
        int[] element_node = new int[4 * element_num];
        double[] node_xy = new double[2 * node_num];

        Burkardt.Values.QuadMesh.example1_q4_mesh(node_num, element_num, ref node_xy, ref element_node,
            ref element_neighbor);

        const int sample_num = 20;

        double[] sample_xy = new double[2 * sample_num];
        int[] sample_element = new int[sample_num];

        int seed = 123456789;

        Sample.sample_q4_mesh(node_num, node_xy, element_num, element_node,
            sample_num, ref seed, ref sample_xy, ref sample_element);

        Console.WriteLine("");
        Console.WriteLine("             X        Y     Element");
        Console.WriteLine("");

        for (sample = 0; sample < sample_num; sample++)
        {
            Console.WriteLine("  " + sample.ToString().PadLeft(8)
                                   + "  " + sample_xy[0 + sample * 2].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + sample_xy[1 + sample * 2].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + sample_element[sample].ToString().PadLeft(8) + "");
        }
    }

    private static void sample_quad_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SAMPLE_QUAD_TEST demonstrates SAMPLE_QUAD.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 February 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const string output_filename = "sample_quad.txt";
        double[] quad_xy =
        {
            1.0, 2.0,
            5.0, 2.0,
            5.0, 3.0,
            4.0, 4.0
        };

        Console.WriteLine("");
        Console.WriteLine("SAMPLE_QUAD_TEST");
        Console.WriteLine("  SAMPLE_QUAD computes N random points in a quadrilateral.");
        Console.WriteLine("  Write them to a file.");

        const int n = 5000;

        int seed = 123456789;

        double[] xy = Sample.sample_quad_new(quad_xy, n, ref seed);

        typeMethods.r8mat_write(output_filename, 2, n, xy);

        Console.WriteLine("");
        Console.WriteLine("  Point coordinates written to \"" + output_filename + "\".");
    }
}