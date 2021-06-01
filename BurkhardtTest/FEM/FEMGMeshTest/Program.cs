using System;
using Burkardt.Table;
using Burkardt.FEM;

namespace Burkardt.FEMGMeshTest
{
    class Program
    {
        static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for FEM_TO_GMSH.
        //
        //  Discussion:
        //
        //    FEM_TO_GMSH converts a 1D, 2D or 3D mesh from FEM to GMSH format.
        //
        //  Usage:
        //
        //    fem_to_gmsh prefix
        //
        //    where 'prefix' is the common filename prefix:
        //
        //    * 'prefix'_nodes.txt contains the node coordinates,
        //    * 'prefix'_elements.txt contains the element node connectivity.
        //    * 'prefix'.msh will contain the Gmsh version of the data.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 October 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            string element_filename;
            int[] element_node = new int[1];
            int element_num = 0;
            int element_order = 0;
            string gmsh_filename;
            int m = 0;
            string node_filename;
            int node_num = 0;
            double[] node_x = new double[1];
            string prefix;

            Console.WriteLine("");

            Console.WriteLine("");
            Console.WriteLine("FEM_TO_GMSH");
            Console.WriteLine("  Convert a 1D, 2D or 3D mesh from FEM to GMSH format.");
            Console.WriteLine("");
            Console.WriteLine("  Read \"prefix\"_nodes.txt, node coordinates.");
            Console.WriteLine("  Read \"prefix\"_elements.txt, element node connectivity.");
            Console.WriteLine("");
            Console.WriteLine("  Create \"prefix\".msh, a corresponding Gmsh mesh file.");

            if (args.Length < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("FEM_TO_GMSH:");
                Console.WriteLine("  Please enter the filename prefix.");

                prefix = Console.ReadLine();
            }
            else
            {
                prefix = args[0];
            }

            //
            //  Create the filenames.
            //
            node_filename = prefix + "_nodes.txt";
            element_filename = prefix + "_elements.txt";
            gmsh_filename = prefix + ".msh";
            //
            //  Read the node data.
            //
            TableHeader h = TableReader.r8mat_header_read(node_filename);
            node_num = h.n;
            m = h.m;

            Console.WriteLine("");
            Console.WriteLine("  Read the header of \"" + node_filename + "\".");
            Console.WriteLine("");
            Console.WriteLine("  Spatial dimension = " + m + "");
            Console.WriteLine("  Number of nodes  = " + node_num + "");

            node_x = TableReader.r8mat_data_read(node_filename, m, node_num);

            Console.WriteLine("");
            Console.WriteLine("  Read the data in \"" + node_filename + "\".");

            TableMisc.r8mat_transpose_print_some(m, node_num, node_x, 1, 1, m, 5,
                "  Portion of node coordinate data:");
            //
            //  Read the element data.
            //
            h = TableReader.i4mat_header_read(element_filename);
            element_order = h.m;
            element_num = h.n;

            if (m == 1)
            {
                if (element_order == 2)
                {
                }
                else
                {
                }
            }
            else if (m == 2)
            {
                if (element_order == 3)
                {
                }
                else if (element_order == 6)
                {
                }
                else
                {
                    Console.Write("");
                    Console.Write("FEM_TO_GMSH - Fatal error!");
                    Console.Write("2D mesh data must use 3 or 6 nodes.");
                    return;
                }
            }
            else if (m == 3)
            {
                if (element_order == 4)
                {
                }
                else if (element_order == 10)
                {
                }
                else if (element_order == 20)
                {
                }
                else
                {
                    Console.Write("");
                    Console.Write("FEM_TO_GMSH - Fatal error!");
                    Console.Write("  3D mesh data must use 4, 10, or 20 nodes.");
                    return;
                }
            }

            Console.WriteLine("");
            Console.WriteLine("  Read the header of \"" + element_filename + "\".");
            Console.WriteLine("");
            Console.WriteLine("  Element order = " + element_order + "");
            Console.WriteLine("  Number of elements  = " + element_num + "");

            element_node = TableReader.i4mat_data_read(element_filename, element_order,
                element_num);

            Console.WriteLine("");
            Console.WriteLine("  Read the data in \"" + element_filename + "\".");

            TableMisc.i4mat_transpose_print_some(element_order, element_num, element_node,
                1, 1, element_order, 10, "  Initial portion of element data:");
            //
            //  Write out the Gmsh version of the data.
            //
            if (m == 1)
            {
                GMesh.gmsh_mesh1d_write(gmsh_filename, m, node_num, node_x, element_order,
                    element_num, element_node);
            }
            else if (m == 2)
            {
                GMesh.gmsh_mesh2d_write(gmsh_filename, m, node_num, node_x, element_order,
                    element_num, element_node);
            }
            else if (m == 3)
            {
                GMesh.gmsh_mesh3d_write(gmsh_filename, m, node_num, node_x, element_order,
                    element_num, element_node);
            }

            Console.WriteLine("");
            Console.WriteLine("  Created the GMSH file \"" + gmsh_filename + "\".");
            //
            //  Terminate.
            //
            Console.WriteLine("");
            Console.WriteLine("FEM_TO_GMSH:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }
    }
}