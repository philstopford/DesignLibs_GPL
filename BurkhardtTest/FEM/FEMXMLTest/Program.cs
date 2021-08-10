using System;
using Burkardt;
using Burkardt.FEM;
using Burkardt.Table;
using Burkardt.Types;

namespace FEMXMLTest
{
    class Program
    {
        static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for FEM_TO_XML.
        //
        //  Discussion:
        //
        //    FEM_TO_XML converts mesh data from FEM to DOLFIN XML format.
        //
        //  Usage:
        //
        //    fem_to_xml prefix
        //
        //    where 'prefix' is the common filename prefix:
        //
        //    * 'prefix'_nodes.txt contains the node coordinates,
        //    * 'prefix'_elements.txt contains the element definitions.
        //    * 'prefix'.xml will contain the DOLFIN XML version of the data.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 October 2014
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
        //    ISBN13: 978-364223098
        //
        {
            string element_filename;
            int[] element_node;
            int element_num;
            int element_order;
            int m;
            string node_filename;
            int node_num;
            double[] node_x;
            string prefix;
            string xml_filename;

            Console.WriteLine("");
            Console.WriteLine("FEM_TO_XML");
            Console.WriteLine("  Convert 1D, 2D or 3D mesh data from FEM to DOLFIN XML format.");
            Console.WriteLine("");
            Console.WriteLine("  Read \"prefix\"_nodes.txt, node coordinates.");
            Console.WriteLine("  Read \"prefix\"_elements.txt, element node connectivity.");
            Console.WriteLine("");
            Console.WriteLine("  Create \"prefix\".xml, a corresponding DOLFIN XML file.");
            //
            //  Get the filename prefix.
            //
            if (args.Length < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("FEM_TO_XML:");
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
            xml_filename = prefix + ".xml";
//
//  Read the node data.
//
            TableHeader h = typeMethods.r8mat_header_read(node_filename);
            m = h.m;
            node_num = h.n;

            Console.WriteLine("");
            Console.WriteLine("  Read the header of \"" + node_filename + "\".");
            Console.WriteLine("");
            Console.WriteLine("  Spatial dimension = " + m + "");
            Console.WriteLine("  Number of nodes  = " + node_num + "");

            node_x = typeMethods.r8mat_data_read(node_filename, m, node_num);

            Console.WriteLine("");
            Console.WriteLine("  Read the data in \"" + node_filename + "\".");

            typeMethods.r8mat_transpose_print_some(m, node_num, node_x, 1, 1, m, 5,
                "  Portion of node coordinate data:");
//
//  Read the element data.
//
            h = typeMethods.i4mat_header_read(element_filename);
            element_order = h.m;
            element_num = h.n;

            if (m == 1)
            {
                if (element_order == 2)
                {
                }
                else
                {
                    Console.WriteLine("");
                    Console.WriteLine("FEM_TO_XML - Fatal error!");
                    Console.WriteLine("  1D element data must use 2 vertices.");
                    return;
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
                    Console.WriteLine("");
                    Console.WriteLine("FEM_TO_XML - Fatal error!");
                    Console.WriteLine("  2D element data must use 3 vertices.");
                    return;
                }
            }
            else if (m == 3)
            {
                if (element_order == 4)
                {
                }
                else
                {
                    Console.WriteLine("");
                    Console.WriteLine("FEM_TO_XML - Fatal error!");
                    Console.WriteLine("  3D element data must use 4 vertices.");
                    return;
                }
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

            typeMethods.i4mat_transpose_print_some(element_order, element_num, element_node,
                1, 1, element_order, 10, "  Initial portion of element data:");
            //
            //  Write the XML version of the data.
            //
            if (m == 1)
            {
                XML.xml_mesh1d_write(xml_filename, m, node_num, node_x, element_order,
                    element_num, element_node);
            }
            else if (m == 2)
            {
                XML.xml_mesh2d_write(xml_filename, m, node_num, node_x, element_order,
                    element_num, element_node);
            }
            else if (m == 3)
            {
                XML.xml_mesh3d_write(xml_filename, m, node_num, node_x, element_order,
                    element_num, element_node);
            }

            Console.WriteLine("");
            Console.WriteLine("  Created XML file \"" + xml_filename + "\".");

            Console.WriteLine("");
            Console.WriteLine("FEM_TO_XML:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }
    }
}
