using System;
using Burkardt;
using Burkardt.MeshNS;
using Burkardt.Table;
using Burkardt.TetrahedronNS;
using Burkardt.Types;

namespace TetMeshBoundaryTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for TET_MESH_BOUNDARY.
            //
            //  Discussion:
            //
            //    TET_MESH_BOUNDARY reads files defining a tetrahedral mesh, and
            //    creates new files that define the 3D triangular mesh formed
            //    by the boundary of the tetrahedral mesh.
            //
            //  Usage:
            //
            //    tet_mesh_boundary prefix
            //
            //    where
            //
            //    * prefix_nodes.txt contains the node coordinates;
            //    * prefix_elements.txt contains the element definitions;
            //    * 'prefix'_boundary_node_mask.txt is 0 for interior nodes,
            //       1 for boundary nodes;
            //    * prefix_boundary_nodes.txt contains the boundary nodes.
            //    * prefix_boundary_elements.txt contains the boundary elements (faces).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    25 December 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            string boundary_element_filename;
            int[] boundary_element_node;
            int boundary_element_num = 0;
            int boundary_element_order = 0;
            string boundary_node_filename;
            int[] boundary_node_index;
            int[] boundary_node_mask;
            string boundary_node_mask_filename;
            int boundary_node_num = 0;
            double[] boundary_node_xyz;
            int dim_num;
            string element_filename;
            int[] element_node;
            int element_num;
            int element_order;
            int i;
            int j;
            int j2;
            int node_num;
            string node_filename;
            double[] node_xyz;
            string prefix;

            Console.WriteLine("");
            Console.WriteLine("TET_MESH_BOUNDARY");
            Console.WriteLine("  Read files defining a tet mesh dataset.");
            Console.WriteLine("");
            Console.WriteLine("  Determine the faces that form the boundary;");
            Console.WriteLine("  write a boundary node mask file;");
            Console.WriteLine("  write a boundary node file;");
            Console.WriteLine("  write a boundary element file,");
            Console.WriteLine("  defining the boundary as a TRI_SURFACE.");
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
                Console.WriteLine("TET_MESH_BOUNDARY:");
                Console.WriteLine("  Please enter the filename prefix.");

                prefix = Console.ReadLine();
            }

            //
            //  Create the filenames.
            //
            node_filename = prefix + "_nodes.txt";
            element_filename = prefix + "_elements.txt";
            boundary_node_mask_filename = prefix + "_boundary_node_mask.txt";
            boundary_node_filename = prefix + "_boundary_nodes.txt";
            boundary_element_filename = prefix + "_boundary_elements.txt";
            //
            //  Read the nodes.
            //
            TableHeader h = typeMethods.r8mat_header_read(node_filename);
            dim_num = h.m;
            node_num = h.n;

            Console.WriteLine("");
            Console.WriteLine("  Read the header of \"" + node_filename + "\".");

            if (dim_num != 3)
            {
                Console.WriteLine("");
                Console.WriteLine("TET_MESH_BOUNDARY - Fatal error!");
                Console.WriteLine("  The spatial dimension must be 3.");
                return;
            }

            Console.WriteLine("");
            Console.WriteLine("  Spatial dimension = " + dim_num + "");
            Console.WriteLine("  Number of nodes   = " + node_num + "");

            node_xyz = typeMethods.r8mat_data_read(node_filename, dim_num, node_num);

            Console.WriteLine("");
            Console.WriteLine("  Read the data in \"" + node_filename + "\".");

            typeMethods.r8mat_transpose_print_some(dim_num, node_num,
                node_xyz, 1, 1, dim_num, 5, "  First 5 nodes:");
            //
            //  Read the elements.
            //
            h = typeMethods.i4mat_header_read(element_filename);
            element_order = h.m;
            element_num = h.n;

            Console.WriteLine("");
            Console.WriteLine("  Read the header of \"" + element_filename + "\".");

            if (element_order != 4 && element_order != 10)
            {
                Console.WriteLine("");
                Console.WriteLine("TET_MESH_BOUNDARY - Fatal error!");
                Console.WriteLine("  The elements must have order 4 or 10.");
                return;
            }

            Console.WriteLine("");
            Console.WriteLine("  Element order = " + element_order + "");
            Console.WriteLine("  Number of tetras  = " + element_num + "");

            element_node = typeMethods.i4mat_data_read(element_filename, element_order,
                element_num);

            Console.WriteLine("");
            Console.WriteLine("  Read the data in \"" + element_filename + "\".");

            typeMethods.i4mat_transpose_print_some(element_order, element_num,
                element_node, 1, 1, element_order, 5, "  First 5 elements:");
            //
            //  If the element information is 1-based, make it 0-based.
            //
            Mesh.mesh_base_zero(node_num, element_order, element_num, ref element_node);
            //
            //  Count the boundary faces and nodes.
            //
            boundary_node_mask = new int[node_num];

            TetMesh_Boundary.tet_mesh_boundary_count(element_order, element_num, element_node,
                node_num, ref boundary_node_num, ref boundary_element_num, ref boundary_node_mask);

            Console.WriteLine("");
            Console.WriteLine("  Number of faces is " + 4 * element_num + "");
            Console.WriteLine("  Number of boundary faces is " + boundary_element_num + "");
            Console.WriteLine("  Number of boundary nodes is " + boundary_node_num + "");
            //
            //  Set the boundary nodes and write them out.
            //
            boundary_node_xyz = new double[dim_num * boundary_node_num];

            j2 = 0;
            for (j = 0; j < node_num; j++)
            {
                if (boundary_node_mask[j] == 1)
                {
                    for (i = 0; i < 3; i++)
                    {
                        boundary_node_xyz[i + j2 * 3] = node_xyz[i + j * 3];
                    }

                    j2 = j2 + 1;
                }
            }

            typeMethods.r8mat_write(boundary_node_filename, dim_num, boundary_node_num,
                boundary_node_xyz);

            Console.WriteLine("");
            Console.WriteLine("  Wrote the boundary node coordinates to \""
                              + boundary_node_filename + "\".");
            //
            //  Write the boundary node mask data.
            //
            typeMethods.i4mat_write(boundary_node_mask_filename, 1, node_num,
                boundary_node_mask);

            Console.WriteLine("  Wrote the boundary node mask to \""
                              + boundary_node_mask_filename + "\".");
            //
            //  Compute the reduced indices for the boundary nodes.
            //
            boundary_node_index = typeMethods.i4vec_cum(node_num, boundary_node_mask);
            //
            //  Subtract 1 so boundary node indices are zero-based.
            //
            for (i = 0; i < node_num; i++)
            {
                boundary_node_index[i] = boundary_node_index[i] - 1;
            }

            //
            //  Set the boundary faces, applying the reduced index labels, and write them out.
            //
            if (element_order == 4)
            {
                boundary_element_order = 3;
            }
            else if (element_order == 10)
            {
                boundary_element_order = 6;
            }

            boundary_element_node = TetMesh_Boundary.tet_mesh_boundary_set(element_order, element_num,
                element_node, boundary_element_order, boundary_element_num);

            for (j = 0; j < boundary_element_num; j++)
            {
                for (i = 0; i < boundary_element_order; i++)
                {
                    boundary_element_node[i + j * boundary_element_order] =
                        boundary_node_index[boundary_element_node[i + j * boundary_element_order]];
                }
            }

            typeMethods.i4mat_write(boundary_element_filename, boundary_element_order,
                boundary_element_num, boundary_element_node);

            Console.WriteLine("  Wrote the boundary face coordinates to \""
                              + boundary_element_filename + "\".");

            Console.WriteLine("");
            Console.WriteLine("TET_MESH_BOUNDARY:");
            Console.WriteLine("  Normal end of execution.");
        }
    }
}