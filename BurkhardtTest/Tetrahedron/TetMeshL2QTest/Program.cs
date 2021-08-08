using System;
using Burkardt.Table;
using Burkardt.TetrahedronNS;
using Burkardt.Types;

namespace TetMeshL2QTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for TET_MESH_L2Q.
            //
            //  Discussion:
            //
            //    TET_MESH_L2Q makes a quadratic tet mesh from a linear one.
            //
            //  Usage:
            //
            //    tet_mesh_l2q prefix
            //
            //    where prefix is the common file prefix:
            //
            //    * prefix_nodes.txt,    the node coordinates;
            //    * prefix_elements.txt,    the linear element definitions.
            //    * prefix_l2q_nodes.txt, the quadratic node coordinates,
            //    * prefix_l2q_elements.txt,    the quadratic element definitions.
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
            int dim_num;
            int[] edge_data;
            int[] element_node1;
            int[] element_node2;
            int element_num;
            int element_order1;
            int element_order2 = 10;
            string input_node_filename;
            string input_element_filename;
            int node_num1;
            int node_num2 = 0;
            double[] node_xyz1;
            double[] node_xyz2;
            string output_node_filename;
            string output_element_filename;
            string prefix;

            Console.WriteLine("");
            Console.WriteLine("TET_MESH_L2Q");
            Console.WriteLine("  Read a \"linear\" tet mesh and");
            Console.WriteLine("  write out a \"quadratic\" one.");
            Console.WriteLine("");
            Console.WriteLine("  Read a node file of NODE_NUM1 nodes in 3 dimensions.");
            Console.WriteLine("  Read an associated tet mesh of TETRA_NUM");
            Console.WriteLine("  tetrahedrons, using 4 nodes per tetrahedron.");
            Console.WriteLine("");
            Console.WriteLine("  Create new nodes which are midpoints of sides,");
            Console.WriteLine("  generate new node and tet mesh data for");
            Console.WriteLine("  quadratic 10-node tetrahedrons, and write them out.");
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
                Console.WriteLine("TET_MESH_L2Q:");
                Console.WriteLine("  Please enter the filename prefix.");

                prefix = Console.ReadLine();
            }

            //
            //  Create the filenames.
            //
            input_node_filename = prefix + "_nodes.txt";
            input_element_filename = prefix + "_elements.txt";
            output_node_filename = prefix + "_l2q_nodes.txt";
            output_element_filename = prefix + "_l2q_elements.txt";
            //
            //  Read the node data.
            //
            TableHeader h = typeMethods.r8mat_header_read(input_node_filename );
            dim_num = h.m;
            node_num1 = h.n;

            if (dim_num != 3)
            {
                Console.WriteLine("");
                Console.WriteLine("TET_MESH_L2Q - Fatal error!");
                Console.WriteLine("  The spatial dimension must be 3.");
                Console.WriteLine("  This data has dimension = " + dim_num + "");
                return;
            }

            Console.WriteLine("");
            Console.WriteLine("  Read the header of \"" + input_node_filename + "\".");
            Console.WriteLine("");
            Console.WriteLine("  Spatial dimension = " + dim_num + "");
            Console.WriteLine("  Number of nodes   = " + node_num1 + "");

            node_xyz1 = typeMethods.r8mat_data_read(input_node_filename, dim_num,
                node_num1);

            Console.WriteLine("");
            Console.WriteLine("  Read the data in \"" + input_node_filename + "\".");

            typeMethods.r8mat_transpose_print_some(dim_num, node_num1,
                node_xyz1, 1, 1, dim_num, 5, "  First 5 nodes:");
            //
            //  Read the tet mesh data.
            //
            h = typeMethods.i4mat_header_read(input_element_filename );
            element_order1 = h.m;
            element_num = h.n;

            if (element_order1 != 4)
            {
                Console.WriteLine("");
                Console.WriteLine("TET_MESH_L2Q - Fatal error!");
                Console.WriteLine("  The tet mesh must have order 4.");
                Console.WriteLine("  This mesh has order " + element_order1 + "");
                return;
            }

            Console.WriteLine("");
            Console.WriteLine("  Read the header of \"" + input_element_filename + "\".");
            Console.WriteLine("");
            Console.WriteLine("  Tetrahedron order = " + element_order1 + "");
            Console.WriteLine("  Number of tetras  = " + element_num + "");

            element_node1 = typeMethods.i4mat_data_read(input_element_filename, element_order1,
                element_num);

            Console.WriteLine("");
            Console.WriteLine("  Read the data in \"" + input_element_filename + "\".");

            typeMethods.i4mat_transpose_print_some(element_order1, element_num,
                element_node1, 1, 1, element_order1, 5, "  First 5 tetrahedrons:");
            //
            //  If the element information is 1-based, make it 0-based.
            //
            TetMesh.tet_mesh_base_zero(node_num1, element_order1, element_num, ref element_node1);
            //
            //  Compute the quadratic mesh.
            //
            edge_data = new int[5 * 6 * element_num];

            TetMesh_L2Q.tet_mesh_order4_to_order10_size(element_num, element_node1, node_num1,
                ref edge_data, ref node_num2);

            Console.WriteLine("  Number of quadratic nodes = " + node_num2 + "");

            node_xyz2 = new double[dim_num * node_num2];
            element_node2 = new int[element_order2 * element_num];

            TetMesh_L2Q.tet_mesh_order4_to_order10_compute(element_num, element_node1, node_num1,
                node_xyz1, edge_data, element_node2, node_num2, node_xyz2);
            //
            //  Print a small amount of the quadratic data.
            //
            typeMethods.r8mat_transpose_print_some(dim_num, node_num2, node_xyz2,
                1, 1, dim_num, 5, "  First 5 quadratic nodes:");

            typeMethods.i4mat_transpose_print_some(element_order2, element_num, element_node2,
                1, 1, element_order2, 5, "  First 5 quadratic tetras");
            //
            //  Write out the node and tetra data for the quadratic mesh
            //
            typeMethods.r8mat_write(output_node_filename, dim_num, node_num2, node_xyz2);

            Console.WriteLine("");
            Console.WriteLine("  Wrote the file \"" + output_node_filename + "\".");

            typeMethods.i4mat_write(output_element_filename, element_order2, element_num, element_node2);

            Console.WriteLine("  Wrote the file \"" + output_element_filename + "\".");
 
            Console.WriteLine("");
            Console.WriteLine("TET_MESH_L2Q:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }
    }
}