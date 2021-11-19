using System;
using Burkardt.GMesh;
using Burkardt.Types;

namespace GMeshIOTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for GMSH_IO_TEST.
        //
        //  Discussion:
        //
        //    GMSH_IO_TEST tests the GMSH_IO library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 October 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("GMSH_IO_TEST");
        Console.WriteLine("  Test the GMSH_IO library.");

        test01();
        test02();

        Console.WriteLine("");
        Console.WriteLine("GMSH_IO_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 gets the example 2D data and writes it to a file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 October 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] element_node;
        int element_num = 0;
        int element_order = 0;
        string gmsh_filename = "example_2d.msh";
        int m = 0;
        int node_num = 0;
        double[] node_x;

        Console.WriteLine("");
        Console.WriteLine("TEST01:");
        Console.WriteLine("  Get example 2D data, write to a file.");
        //
        //  Get sizes.
        //
        Mesh2D.gmsh_mesh2d_node_size_example(ref node_num, ref m);

        Mesh2D.gmsh_mesh2d_element_size_example(ref element_num, ref element_order);
        //
        //  Print the sizes.
        //
        Console.WriteLine("");
        Console.WriteLine("  Number of nodes = " + node_num + "");
        Console.WriteLine("  Spatial dimension = " + m + "");
        Console.WriteLine("  Number of elements = " + element_num + "");
        Console.WriteLine("  Order of elements = " + element_order + "");
        //
        //  Get the data.
        //
        node_x = Mesh2D.gmsh_mesh2d_node_data_example(node_num, m);

        element_node = Mesh2D.gmsh_mesh2d_element_data_example(element_num, element_order);
        //
        //  Print some of the data.
        //
        typeMethods.r8mat_transpose_print_some(m, node_num, node_x,
            1, 1, m, 10, "  Coordinates for first 10 nodes:");

        typeMethods.i4mat_transpose_print_some(element_order, element_num, element_node,
            1, 1, element_order, 10, "  Node connectivity of first 10 elements:");
        //
        //  Write the GMSH file.
        //
        Mesh2D.gmsh_mesh2d_write(gmsh_filename, m, node_num, node_x,
            element_order, element_num, element_node);

        Console.WriteLine("");
        Console.WriteLine("  Wrote example data to file \"" + gmsh_filename + "\"");
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 reads the example data from a file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 October 2014
        //
        //  Author:
        //
        //   John Burkardt
        //
    {
        int[] element_node;
        int element_num = 0;
        int element_order = 0;
        string gmsh_filename = "example_2d.msh";
        int m = 0;
        int node_num = 0;
        double[] node_x;

        Console.WriteLine("");
        Console.WriteLine("TEST02:");
        Console.WriteLine("  Read data from a file.");
        //
        //  Get the data size.
        //
        IO.gmsh_size_read(gmsh_filename, ref node_num, ref m, ref element_num,
            ref element_order);
        //
        //  Print the sizes.
        //
        Console.WriteLine("");
        Console.WriteLine("  Node data read from file \"" + gmsh_filename + "\"");
        Console.WriteLine("");
        Console.WriteLine("  Number of nodes = " + node_num + "");
        Console.WriteLine("  Spatial dimension = " + m + "");
        Console.WriteLine("  Number of elements = " + element_num + "");
        Console.WriteLine("  Element order = " + element_order + "");
        //
        //  Allocate memory.
        //
        node_x = new double[m * node_num];
        element_node = new int[element_order * element_num];
        //
        //  Get the data.
        //
        IO.gmsh_data_read(gmsh_filename, m, node_num, node_x,
            element_order, element_num, element_node);
        //
        //  Print some of the data.
        //
        typeMethods.r8mat_transpose_print_some(m, node_num, node_x,
            1, 1, m, 10, "  Coordinates for first 10 nodes:");

        typeMethods.i4mat_transpose_print_some(element_order, element_num, element_node,
            1, 1, element_order, 10, "  Connectivity for first 10 elements:");
    }
}