using System;
using Burkardt.Table;
using Burkardt.TetrahedronNS;
using Burkardt.Types;

namespace TetMeshQ2LTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for TET_MESH_Q2L.
            //
            //  Discussion:
            //
            //    TET_MESH_Q2L makes a linear tet mesh from a quadratic one.
            //
            //    A quadratic tet mesh is assumed to consist of 10-node
            //    tetrahedrons.  This routine rearranges information so as to 
            //    define a 4-node tet mesh.
            //
            //  Usage:
            //
            //    tet_mesh_q2l prefix
            //
            //    where prefix is the common file prefix:
            //
            //    * prefix_nodes.txt,    the node coordinates (not needed by this program);
            //    * prefix_elements.txt,    the linear element definitions.
            //    * prefix_q2l_elements.txt,    the quadratic element definitions,
            //                                created by the program.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    01 October 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            string input_element_filename;
            int node_num1 = 0;
            int node_num2 = 0;
            string output_element_filename;
            string prefix;
            int[] tetra_node1;
            int[] tetra_node2;
            int tetra_num1 = 0;
            int tetra_num2 = 0;
            int tetra_order1;
            int tetra_order2 = 4;

            Console.WriteLine("");
            Console.WriteLine("TET_MESH_Q2L");
            Console.WriteLine("  Read a \"quadratic\" tet mesh and");
            Console.WriteLine("  write out a \"linear\" one.");
            Console.WriteLine("");
            Console.WriteLine("  Read a tet mesh of TETRA_NUM1 tetrahedrons");
            Console.WriteLine("  using 10 nodes.");
            Console.WriteLine("");
            Console.WriteLine("  Create a 4 node tet mesh by breaking");
            Console.WriteLine("  every 10 node tetrahedron into 8 smaller ones.");
            Console.WriteLine("  Write the new linear tet mesh to a file.");
            Console.WriteLine("");
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
                Console.WriteLine("TET_MESH_Q2L:");
                Console.WriteLine("  Please enter the filename prefix.");

                prefix = Console.ReadLine();
            }

            //
            //  Create the filenames.
            //
            input_element_filename = prefix + "_elements.txt";
            output_element_filename = prefix + "_q2l_elements.txt";
            //
            //  Read the tet mesh data.
            //
            TableHeader h = typeMethods.i4mat_header_read(input_element_filename);
            tetra_order1 = h.m;
            tetra_num1 = h.n;

            if (tetra_order1 != 10)
            {
                Console.WriteLine("");
                Console.WriteLine("TET_MESH_Q2L - Fatal error!");
                Console.WriteLine("  The tet mesh must have order 10.");
                return;
            }

            Console.WriteLine("");
            Console.WriteLine("  Read the header of \"" + input_element_filename + "\".");
            Console.WriteLine("");
            Console.WriteLine("  Tetrahedron order = " + tetra_order1 + "");
            Console.WriteLine("  Number of tetras  = " + tetra_num1 + "");

            tetra_node1 = typeMethods.i4mat_data_read(input_element_filename, tetra_order1,
                tetra_num1);

            Console.WriteLine("");
            Console.WriteLine("  Read the data in \"" + input_element_filename + "\".");

            typeMethods.i4mat_transpose_print_some(tetra_order1, tetra_num1,
                tetra_node1, 1, 1, tetra_order1, 5, "  First 5 tetrahedrons:");
            //
            //  Set the number of linear tetrahedrons:
            //  We didn't read in the node data, so set that count to 0.
            //
            node_num1 = 0;

            TetMesh_Q2L.tet_mesh_order10_to_order4_size(node_num1, tetra_num1,
                ref node_num2, ref tetra_num2);

            Console.WriteLine("  Number of linear tetrahedrons =    " + tetra_num2 + "");
            //
            //  Allocate space.
            //
            tetra_node2 = new int[tetra_order2 * tetra_num2];
            //
            //  Convert the data.
            //
            TetMesh_Q2L.tet_mesh_order10_to_order4_compute(tetra_num1, tetra_node1,
                tetra_num2, ref tetra_node2);

            typeMethods.i4mat_transpose_print_some(tetra_order2, tetra_num2, tetra_node2,
                1, 1, tetra_order2, 5, "  First 5 linear tetras");
            //
            //  Write out the tetrahedron data for the quadratic mesh
            //
            typeMethods.i4mat_write(output_element_filename, tetra_order2, tetra_num2,
                tetra_node2);

            Console.WriteLine("");
            Console.WriteLine("  Created the file \"" + output_element_filename + "\".");

            Console.WriteLine("");
            Console.WriteLine("TET_MESH_Q2L:");
            Console.WriteLine("  Normal end of execution.");

            Console.WriteLine("");
        }
    }
}