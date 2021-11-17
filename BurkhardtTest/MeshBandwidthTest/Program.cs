using System;
using Burkardt;
using Burkardt.MeshNS;
using Burkardt.Table;
using Burkardt.Types;

namespace MeshBandwidthTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for MESH_BANDWIDTH.
        //
        //  Discussion:
        //
        //    MESH_BANDWIDTH determines the geometric bandwidth of a mesh of elements.
        //
        //    The user supplies an element file, listing the indices of the nodes that
        //    make up each element.
        //
        //    The program computes the geometric bandwidth associated with the mesh.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 September 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Usage:
        //
        //    mesh_bandwidth element_filename
        //
    {
        string element_filename;
        int[] element_node;
        int element_num = 0;
        int element_order = 0;
        int m = 0;
        int ml = 0;
        int mu = 0;

        Console.WriteLine("");
        Console.WriteLine("MESH_BANDWIDTH");
        Console.WriteLine("  Read a mesh file which defines");
        Console.WriteLine("  a \"triangulation\" of a region in the plane,");
        Console.WriteLine("  or a \"tetrahedronization\" of a region in space,");
        Console.WriteLine("  or any division of a regino in ND space into elements,");
        Console.WriteLine("  using a mesh of elements of uniform order.");
        Console.WriteLine("");
        Console.WriteLine("  Determine the geometric mesh bandwidth.");
        Console.WriteLine("");
        Console.WriteLine("    M = ML + 1 + MU.");
        Console.WriteLine("");
        Console.WriteLine("  which is the bandwidth of the vertex connectivity");
        Console.WriteLine("  matrix.");
        Console.WriteLine("");
        Console.WriteLine("  Note that a matrix associated with variables defined");
        Console.WriteLine("  at the  nodes could have a greater bandwidth than M,");
        Console.WriteLine("  since you might have multiple variables at a vertex,");
        Console.WriteLine("  or the variable might be a vector quantity,");
        Console.WriteLine("  or physical effects might link two variables that are");
        Console.WriteLine("  not associated with vertices that are connected.");
        //
        //  If at least one command line argument, it's the element file.
        //
        try
        {
            element_filename = args[0];

        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("MESH_BANDWIDTH:");
            Console.WriteLine("  Please enter the name of the element file.");

            element_filename = Console.ReadLine();
        }

        //
        //  Read the triangulation data.
        //
        TableMisc.readHeader(element_filename, ref element_order, ref element_num);

        Console.WriteLine("");
        Console.WriteLine("  Read the header of \"" + element_filename + "\".");
        Console.WriteLine("");
        Console.WriteLine("  Element order ELEMENT_ORDER =    " + element_order + "");
        Console.WriteLine("  Number of element ELEMENT_NUM  = " + element_num + "");

        element_node = typeMethods.i4mat_data_read(element_filename, element_order,
            element_num);

        Console.WriteLine("");
        Console.WriteLine("  Read the data in \"" + element_filename + "\".");

        typeMethods.i4mat_transpose_print_some(element_order, element_num, element_node, 1, 1,
            element_order, 10, "  Portion of data read from file:");
        //
        //  Compute the bandwidth.
        //
        Mesh.bandwidth_mesh(element_order, element_num, element_node, ref ml, ref mu, ref m);

        Console.WriteLine("");
        Console.WriteLine("  Lower bandwidth ML = " + ml + "");
        Console.WriteLine("  Upper bandwidth MU = " + mu + "");
        Console.WriteLine("  Total bandwidth M  = " + m + "");

        Console.WriteLine("");
        Console.WriteLine("MESH_BANDWIDTH");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }
}