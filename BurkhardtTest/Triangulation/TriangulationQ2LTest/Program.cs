using System;
using Burkardt.Table;
using Burkardt.TriangulationNS;
using Burkardt.Types;

namespace TriangulationQ2LTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for TRIANGULATION_Q2L.
            //
            //  Discussion:
            //
            //    TRIANGULATION_Q2L makes a linear triangulation from a quadratic one.
            //
            //    A quadratic triangulation is assumed to consist of 6-node triangles,
            //    as in the following:
            //
            //    11-12-13-14-15
            //     |\    |\    |
            //     | \   | \   |
            //     6  7  8  9 10
            //     |   \ |   \ |
            //     |    \|    \|
            //     1--2--3--4--5
            //
            //    This routine rearranges information so as to define the 3-node
            //    triangulation:
            //
            //    11-12-13-14-15
            //     |\ |\ |\ |\ |
            //     | \| \| \| \|
            //     6--7--8--9-10
            //     |\ |\ |\ |\ |
            //     | \| \| \| \|
            //     1--2--3--4--5
            //
            //  Usage:
            //
            //    triangulation_q2l prefix
            //
            //    where 'prefix' is the common filename prefix:
            //
            //    * prefix_nodes.txt contains the node coordinates (not needed by this program),
            //    * prefix_elements.txt contains the element definitions.
            //    * prefix_q2l_elements.txt will contain the linearized element definitions.
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
            string prefix;
            string element_filename;
            string element_q2l_filename;
            int[] triangle_node1;
            int[] triangle_node2;
            int triangle_num1;
            int triangle_num2;
            int triangle_order1;
            int triangle_order2;

            Console.WriteLine("");

            Console.WriteLine("");
            Console.WriteLine("TRIANGULATION_Q2L");
            Console.WriteLine("  Read a \"quadratic\" triangulation and");
            Console.WriteLine("  write out a \"linear\" triangulation.");
            Console.WriteLine("");
            Console.WriteLine("  Read a triangulation dataset of TRI_NUM1");
            Console.WriteLine("  triangles using 6 nodes..");
            Console.WriteLine("");
            Console.WriteLine("  Create a 3 node triangulation by breaking");
            Console.WriteLine("  every 6 node triangle into 4 smaller ones.");
            Console.WriteLine("  Write the new linear triangulation to a file.");
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
                Console.WriteLine("TRIANGULATION_Q2L:");
                Console.WriteLine("  Please enter the filename prefix.");

                prefix = Console.ReadLine();
            }

            //
            //  Create the filenames.
            //
            element_filename = prefix + "_elements.txt";
            element_q2l_filename = prefix + "_q2l_elements.txt";
            //
            //  Read the data.
            //
            TableHeader h = typeMethods.i4mat_header_read(element_filename);
            triangle_order1 = h.m;
            triangle_num1 = h.n;

            if (triangle_order1 != 6)
            {
                Console.WriteLine("");
                Console.WriteLine("TRIANGULATION_Q2L - Fatal error!");
                Console.WriteLine("  Data is not for a 6-node triangulation.");
                return;
            }

            Console.WriteLine("");
            Console.WriteLine("  Read the header of \"" + element_filename + "\".");
            Console.WriteLine("");
            Console.WriteLine("  Triangle order TRIANGLE_ORDER1 = " + triangle_order1 + "");
            Console.WriteLine("  Number of triangles TRIANGLE_NUM1 = " + triangle_num1 + "");

            triangle_node1 = typeMethods.i4mat_data_read(element_filename,
                triangle_order1, triangle_num1);

            Console.WriteLine("");
            Console.WriteLine("  Read the data in \"" + element_filename + "\".");

            typeMethods.i4mat_transpose_print_some(triangle_order1, triangle_num1, triangle_node1, 1, 1,
                triangle_order1, 10, "  Portion of data read from file:");
            //
            //  Set the number of linear triangles:
            //
            triangle_num2 = 4 * triangle_num1;
            triangle_order2 = 3;
            //
            //  Convert the data.
            //
            triangle_node2 = Conversion.triangulation_order6_to_order3(triangle_num1,
                triangle_node1);

            typeMethods.i4mat_transpose_print(triangle_order2, triangle_num2,
                triangle_node2, "  TRIANGLE_NODE2");
            //
            //  Write out the node and triangle data for the quadratic mesh
            //
            typeMethods.i4mat_write(element_q2l_filename, triangle_order2,
                triangle_num2, triangle_node2);

            Console.WriteLine("");
            Console.WriteLine("  Wrote the linearized element data in \"" + element_q2l_filename + "\".");

            Console.WriteLine("");
            Console.WriteLine("TRIANGULATION_Q2L:");
            Console.WriteLine("  Normal end of execution.");

            Console.WriteLine("");
        }
    }
}