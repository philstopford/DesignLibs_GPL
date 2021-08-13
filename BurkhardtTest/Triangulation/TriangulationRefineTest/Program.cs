using System;
using Burkardt;
using Burkardt.MeshNS;
using Burkardt.Table;
using Burkardt.TriangulationNS;
using Burkardt.Types;

namespace TriangulationRefineTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for TRIANGULATION_REFINE.
            //
            //  Discussion:
            //
            //    TRIANGULATION_REFINE refines a triangulation by doubling.
            //
            //  Usage:
            //
            //    triangulation_refine prefix
            //
            //    where 'prefix' is the common filename prefix:
            //
            //    * prefix_nodes.txt contains the node coordinates,
            //    * prefix_elements.txt contains the element definitions.
            //    * prefix_ref_nodes.txt will contain the refined nodes;
            //    * prefix_ref_elements.txt will contain the refined elements.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    04 October 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            bool debug = true;
            int dim_num;
            int[] edge_data;
            string element_filename;
            string element_ref_filename;
            string node_filename;
            int node_num1 = 0;
            int node_num2 = 0;
            string node_ref_filename;
            double[] node_xy1;
            double[] node_xy2;
            string prefix;
            int triangle_num1 = 0;
            int triangle_num2 = 0;
            int triangle_order;
            int[] triangle_node1;
            int[] triangle_node2;

            Console.WriteLine("");
            Console.WriteLine("TRIANGULATION_REFINE");
            Console.WriteLine("  Read a \"linear\" or \"quadratic\" triangulation");
            Console.WriteLine("  and write out a refined triangulation.");
            Console.WriteLine("");
            Console.WriteLine("  In particular:");
            Console.WriteLine("");
            Console.WriteLine("  Read a dataset of NODE_NUM1 points in 2 dimensions.");
            Console.WriteLine("");
            Console.WriteLine("  Read an associated triangulation dataset of TRIANGLE_NUM1 ");
            Console.WriteLine("  triangles which use 3 or 6 nodes per triangle.");
            Console.WriteLine("");
            Console.WriteLine("  Subdivide each triangle into 4 triangles,");
            Console.WriteLine("  generate new nodes as midpoints of current nodes.");
            Console.WriteLine("");
            Console.WriteLine("  Write out the new node and triangulation data.");
            Console.WriteLine("");
            Console.WriteLine("  If the input triangulation was Delaunay, then");
            Console.WriteLine("  the output triangulation will be Delaunay.");
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
                Console.WriteLine("TRIANGULATION_REFINE:");
                Console.WriteLine("  Please enter the filename prefix.");

                prefix = Console.ReadLine();
            }

            //
            //  Create the filenames.
            //
            node_filename = prefix + "_nodes.txt";
            element_filename = prefix + "_elements.txt";
            node_ref_filename = prefix + "_ref_nodes.txt";
            element_ref_filename = prefix + "_ref_elements.txt";
            //
            //  Read the node data.
            //
            TableHeader h = typeMethods.r8mat_header_read(node_filename);
            dim_num = h.m;
            node_num1 = h.n;

            Console.WriteLine("");
            Console.WriteLine("  Read the header of \"" + node_filename + "\".");
            Console.WriteLine("");
            Console.WriteLine("  Spatial dimension DIM_NUM = " + dim_num + "");
            Console.WriteLine("  Number of nodes NODE_NUM1 = " + node_num1 + "");

            node_xy1 = typeMethods.r8mat_data_read(node_filename, dim_num, node_num1);

            Console.WriteLine("");
            Console.WriteLine("  Read the data in \"" + node_filename + "\".");

            typeMethods.r8mat_transpose_print_some(dim_num, node_num1, node_xy1, 1, 1, dim_num, 5,
                "  First 5 nodes:");
            //
            //  Read the element data.
            //
            h = typeMethods.i4mat_header_read(element_filename);
            triangle_order = h.m;
            triangle_num1 = h.n;

            if (triangle_order != 3 && triangle_order != 6)
            {
                Console.WriteLine("");
                Console.WriteLine("TRIANGULATION_REFINE - Fatal error!");
                Console.WriteLine("  Data is not for a 3-node or 6-node triangulation.");
                return;
            }

            Console.WriteLine("");
            Console.WriteLine("  Read the header of \"" + element_filename + "\".");
            Console.WriteLine("");
            Console.WriteLine("  Triangle order TRIANGLE_ORDER = " + triangle_order + "");
            Console.WriteLine("  Number of triangles TRIANGLE_NUM1  = " + triangle_num1 + "");

            triangle_node1 = typeMethods.i4mat_data_read(element_filename,
                triangle_order, triangle_num1);

            Console.WriteLine("");
            Console.WriteLine("  Read the data in \"" + element_filename + "\".");

            typeMethods.i4mat_transpose_print_some(triangle_order, triangle_num1, triangle_node1,
                1, 1, triangle_order, 5, "  First 5 triangles:");
            //
            //  Detect and correct 1-based node indexing.
            //
            Mesh.mesh_base_zero(node_num1, triangle_order, triangle_num1, ref triangle_node1);
            //
            //  Determine the size of the refined mesh.
            //
            edge_data = new int[5 * (3 * triangle_num1)];

            if (triangle_order == 3)
            {
                Refine.triangulation_order3_refine_size(node_num1, triangle_num1,
                    triangle_node1, ref node_num2, ref triangle_num2, ref edge_data);
            }
            else if (triangle_order == 6)
            {
                Refine.triangulation_order6_refine_size(node_num1, triangle_num1,
                    triangle_node1, ref node_num2, ref triangle_num2, ref edge_data);
            }

            Console.WriteLine("");
            Console.WriteLine("  Number of nodes in refined mesh =      " + node_num2 + "");
            Console.WriteLine("  Number of triangles in refined mesh =  "
                              + triangle_num2 + "");

            node_xy2 = new double[dim_num * node_num2];
            triangle_node2 = new int[triangle_order * triangle_num2];
            //
            //  Compute the refined mesh.
            //
            if (triangle_order == 3)
            {
                Refine.triangulation_order3_refine_compute(node_num1, triangle_num1,
                    node_xy1, triangle_node1, node_num2, triangle_num2, edge_data, ref node_xy2,
                    ref triangle_node2);
            }
            else if (triangle_order == 6)
            {
                Refine.triangulation_order6_refine_compute(node_num1, triangle_num1,
                    node_xy1, triangle_node1, node_num2, triangle_num2, edge_data, ref node_xy2,
                    ref triangle_node2);
            }

            if (debug)
            {
                typeMethods.r8mat_transpose_print_some(dim_num, node_num2, node_xy2,
                    1, 1, dim_num, 5, "  First 5 output nodes:");

                typeMethods.i4mat_transpose_print_some(triangle_order, triangle_num2, triangle_node2,
                    1, 1, triangle_order, 5, "  First 5 output triangles");
            }

            //
            //  Write out the node and triangle data.
            //
            typeMethods.r8mat_write(node_ref_filename, dim_num, node_num2, node_xy2);

            Console.WriteLine("");
            Console.WriteLine("  Wrote the refined node data to \"" + node_ref_filename + "\".");

            typeMethods.i4mat_write(element_ref_filename, triangle_order, triangle_num2,
                triangle_node2);

            Console.WriteLine("  Wrote the refined element data to \"" + element_ref_filename + "\".");

            Console.WriteLine("");
            Console.WriteLine("TRIANGULATION_REFINE:");
            Console.WriteLine("  Normal end of execution.");

        }
    }
}