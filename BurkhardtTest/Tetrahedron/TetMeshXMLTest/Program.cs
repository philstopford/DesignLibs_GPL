using System;
using Burkardt;
using Burkardt.FEM;
using Burkardt.Table;
using Burkardt.TetrahedronNS;
using Burkardt.Types;

namespace TetMeshXMLTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for TET_MESH_TO_XML.
            //
            //  Discussion:
            //
            //    TET_MESH_TO_XML converts a tet mesh, using 4 or 10 node 
            //    tetrahedral elements, to a DOLFIN XML mesh file.
            //
            //  Usage:
            //
            //    tet_mesh_to_xml prefix
            //
            //    where prefix is the common file prefix:
            //
            //    * prefix_nodes.txt,       the node coordinates;
            //    * prefix_elements.txt,    the linear element definitions.
            //    * prefix.xml,             the DOLFIN XML mesh file.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    04 June 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Anders Logg, Kent-Andre Mardal, Garth Wells,
            //    Automated Solution of Differential Equations by the Finite Element
            //    Method: The FEniCS Book,
            //    Lecture Notes in Computational Science and Engineering,
            //    Springer, 2011,
            //    ISBN13: 978-3642230981
            //
        {
            int dim_num;
            string element_filename;
            int[] element_node;
            int element_num;
            int element_order;
            string node_filename;
            int node_num;
            double[] node_xyz;
            string prefix;
            string xml_filename;

            Console.WriteLine("");
            Console.WriteLine("TET_MESH_TO_XML;");
            Console.WriteLine("  Convert a tet mesh to DOLFIN XML format");
            Console.WriteLine("");
            Console.WriteLine("  Read \"prefix\"_nodes.txt, node coordinates.");
            Console.WriteLine("  Read \"prefix\"_elements.txt, 4, 10 or 20 node element definitions.");
            Console.WriteLine("");
            Console.WriteLine("  Create \"prefix\".xml, a corresponding DOLFIN XML mesh file.");
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
                Console.WriteLine("TET_MESH_TO_XML:");
                Console.WriteLine("  Please enter the filename prefix.");

                prefix = Console.ReadLine();
            }

            //
            //  Create the filenames.
            //
            node_filename = prefix + "_nodes.txt";
            element_filename = prefix + "_elements.txt";
            xml_filename = prefix + ".xml";
            //
            //  Read the node data.
            //
            TableHeader h = typeMethods.r8mat_header_read(node_filename);
            dim_num = h.m;
            node_num = h.n;

            if (dim_num != 3)
            {
                Console.WriteLine("");
                Console.WriteLine("TET_MESH_TO_XML; - Fatal error!");
                Console.WriteLine("  The spatial dimension must be 3.");
                return;
            }

            Console.WriteLine("");
            Console.WriteLine("  Read the header of \"" + node_filename + "\".");
            Console.WriteLine("");
            Console.WriteLine("  Spatial dimension = " + dim_num + "");
            Console.WriteLine("  Number of nodes   = " + node_num + "");

            node_xyz = typeMethods.r8mat_data_read(node_filename, dim_num, node_num);

            Console.WriteLine("");
            Console.WriteLine("  Read the data in \"" + node_filename + "\".");

            typeMethods.r8mat_transpose_print_some(dim_num, node_num,
                node_xyz, 1, 1, dim_num, 5, "  First 5 nodes:");
            //
            //  Read the element data.
            //
            h = typeMethods.i4mat_header_read(element_filename);
            element_order = h.m;
            element_num = h.n;

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
                Console.WriteLine("");
                Console.WriteLine("TET_MESH_TO_XML - Fatal error!");
                Console.WriteLine("  The tet mesh must have order 4, 10, or 20.");
                return;
            }

            Console.WriteLine("");
            Console.WriteLine("  Read the header of \"" + element_filename + "\".");
            Console.WriteLine("");
            Console.WriteLine("  Element order = " + element_order + "");
            Console.WriteLine("  Number of elements  = " + element_num + "");

            element_node = typeMethods.i4mat_data_read(element_filename, element_order,
                element_num);

            Console.WriteLine("");
            Console.WriteLine("  Read the data in \"" + element_filename + "\".");

            typeMethods.i4mat_transpose_print_some(element_order, element_num,
                element_node, 1, 1, element_order, 5, "  First 5 elements:");
            //
            //  Use 0-based indexing.
            //
            TetMesh.tet_mesh_base_zero(node_num, element_order, element_num, ref element_node);
            //
            //  Write the DOLFIN XML file.
            //
            XML.xml_write(xml_filename, dim_num, node_num, node_xyz, element_order,
                element_num, element_node);

            Console.WriteLine("");
            Console.WriteLine("  Created XML file \"" + xml_filename + "\".");

            Console.WriteLine("");
            Console.WriteLine("TET_MESH_TO_XML;");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }
    }
}