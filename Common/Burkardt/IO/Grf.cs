using System;
using System.Collections.Generic;
using System.IO;
using System.Linq;
using Burkardt.Types;

namespace Burkardt.IO;

public static class Grf
{
    public static void grf_data_print(int node_num, int edge_num, int[] edge_pointer,
            int[] edge_data, double[] xy )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRF_DATA_PRINT prints the data of a GRF file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 January 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Stephen Skiena,
        //    Implementing Discrete Mathematics,
        //    Combinatorics and Graph Theory with Mathematica,
        //    Addison-Wesley, 1990.
        //
        //  Parameters:
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, int EDGE_NUM, the number of edges.
        //
        //    Input, int EDGE_POINTER[NODE_NUM+1], pointers to 
        //    the beginning of edge data for each node.
        //
        //    Input, int EDGE_DATA[EDGE_NUM], the edge data.
        //
        //    Input, double XY[2*NODE_NUM], the node coordinates.
        //
    {
        int edge;
        int node;

        Console.WriteLine("");
        Console.WriteLine("  Edge pointers:");
        Console.WriteLine("");
        Console.WriteLine("  Node     First      Last");
        Console.WriteLine("");
        for (node = 0; node < node_num; node++)
        {
            Console.WriteLine("  " + node.ToString().PadLeft(4)
                                   + "  " + edge_pointer[node].ToString().PadLeft(8)
                                   + "  " + (edge_pointer[node + 1] - 1).ToString().PadLeft(8) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  Edge data:");
        Console.WriteLine("");
        Console.WriteLine("  Node     Adjacent nodes");
        Console.WriteLine("");

        for (node = 0; node < node_num; node++)
        {
            string tmp = "  " + node.ToString().PadLeft(4);
            for (edge = edge_pointer[node]; edge <= edge_pointer[node + 1] - 1; edge++)
            {
                tmp += "  " + edge_data[edge].ToString().PadLeft(8);
            }

            Console.WriteLine(tmp);
        }

        Console.WriteLine("");
        Console.WriteLine("  Node        X          Y");
        Console.WriteLine("");

        for (node = 0; node < node_num; node++)
        {
            Console.WriteLine("  " + node.ToString().PadLeft(4)
                                   + "  " + xy[0 + node * 2].ToString().PadLeft(10)
                                   + "  " + xy[1 + node * 2].ToString().PadLeft(10) + "");
        }
    }

    public static void grf_data_read(string input_filename, int node_num, int edge_num,
            ref int[] edge_pointer, ref int[] edge_data, ref double[] xy )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRF_DATA_READ reads the data of a GRF file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 January 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Stephen Skiena,
        //    Implementing Discrete Mathematics,
        //    Combinatorics and Graph Theory with Mathematica,
        //    Addison-Wesley, 1990.
        //
        //  Parameters:
        //
        //    Input, string INPUT_FILENAME, the name of the file.
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, int EDGE_NUM, the number of edges.
        //
        //    Output, int EDGE_POINTER[NODE_NUM+1], pointers to 
        //    the beginning of edge data for each node.
        //
        //    Output, int EDGE_DATA[EDGE_NUM], the edge data.
        //
        //    Output, double XY[2*NODE_NUM], the node coordinates.
        //
    {
        int edge;
        int i;
        string[] input_unit;
        int n;
        int node;
        int node_j;
        double xval;
        double yval;

        for (edge = 0; edge < edge_num; edge++)
        {
            edge_data[edge] = -1;
        }

        for (node = 0; node < node_num + 1; node++)
        {
            edge_pointer[node] = -1;
        }

        try
        {
            input_unit = File.ReadAllLines(input_filename);
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("GRF_DATA_READ - Fatal error!");
            Console.WriteLine("  Cannot open the input file \"" + input_filename + "\".");
            return;
        }

        //
        //  Read a line.  If it's a blank or comment, skip it.
        //  Otherwise, count the number of "words", and then reread it.
        //
        edge = 0;
        node = 0;
        edge_pointer[0] = 0;
            
        string text = "";

        for (node = 0; node < node_num; node++)
        {
            try
            {
                text = input_unit[node];
            }
            catch
            {
                Console.WriteLine("");
                Console.WriteLine("GRF_DATA_READ - Fatal error!");
                Console.WriteLine("  Unexpected end of information;");
                return;
            }

            if (typeMethods.s_len_trim(text) <= 0)
            {
                continue;
            }

            switch (text[0])
            {
                case '#':
                    continue;
            }

            string[] tokens = Helpers.splitStringByWhitespace(text);
                
            n = tokens.Length;

            switch (n)
            {
                case < 3:
                    Console.WriteLine("");
                    Console.WriteLine("GRF_DATA_READ - Fatal error!");
                    Console.WriteLine("  Record has less than 3 items.");
                    return;
            }

            xval = Convert.ToDouble(tokens[1]);
            yval = Convert.ToDouble(tokens[2]);

            edge_pointer[node + 1] = edge_pointer[node] + n - 3;

            xy[0 + node * 2] = xval;
            xy[1 + node * 2] = yval;

            for (i = n - 3; i < n; i++)
            {
                node_j = Convert.ToInt32(tokens[i]);
                edge_data[edge] = node_j;
                edge += 1;
            }
        }
    }
        
    public static void grf_data_write(ref List<string> output_unit, int node_num, int edge_num,
            int[] edge_pointer, int[] edge_data, double[] xy )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRF_DATA_WRITE writes the data of a GRF file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 January 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Stephen Skiena,
        //    Implementing Discrete Mathematics,
        //    Combinatorics and Graph Theory with Mathematica,
        //    Addison-Wesley, 1990.
        //
        //  Parameters:
        //
        //    Input, ofstream &OUTPUT_UNIT, a pointer to the GRF file.
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, int EDGE_NUM, the number of edges.
        //
        //    Input, int EDGE_POINTER[NODE_NUM+1], pointers to 
        //    the beginning of edge data for each node.
        //
        //    Input, int EDGE_DATA[EDGE_NUM], the edge data.
        //
        //    Input, double XY[2*NODE_NUM], the node coordinates.
        //
    {
        int edge;
        int node;

        for (node = 0; node < node_num; node++)
        {
            string tmp = "  " + (node + 1).ToString().PadLeft(4)
                              + "  " + xy[0 + node * 2].ToString().PadLeft(10)
                              + "  " + xy[1 + node * 2].ToString().PadLeft(10);
            for (edge = edge_pointer[node]; edge <= edge_pointer[node + 1] - 1; edge++)
            {
                tmp += "  " + edge_data[edge].ToString().PadLeft(8);
            }

            output_unit.Add(tmp);
        }
    }

    public static void grf_example(int node_num, int edge_num, ref int[] edge_pointer,
            ref int[] edge_data, ref double[] xy )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRF_EXAMPLE sets up a GRF example.
        //
        //  Discussion:
        //
        //    The example is known as the Coxeter graph.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 January 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Stephen Skiena,
        //    Implementing Discrete Mathematics,
        //    Combinatorics and Graph Theory with Mathematica,
        //    Addison-Wesley, 1990.
        //
        //  Parameters:
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, int EDGE_NUM, the number of edges.
        //
        //    Output, int EDGE_POINTER[NODE_NUM+1], pointers to 
        //    the beginning of edge data for each node.
        //
        //    Output, int EDGE_DATA[EDGE_NUM], the edge data.
        //
        //    Output, double XY[2*NODE_NUM], the node coordinates.
        //
    {
        int[] EDGE_DATA =  {
                8, 2, 3,
                14, 1, 5,
                9, 4, 1,
                10, 7, 3,
                13, 2, 6,
                12, 5, 7,
                11, 6, 4,
                25, 20, 1,
                24, 21, 3,
                23, 15, 4,
                22, 16, 7,
                28, 17, 6,
                27, 18, 5,
                26, 19, 2,
                10, 18, 19,
                11, 19, 20,
                12, 21, 20,
                13, 15, 21,
                14, 16, 15,
                8, 17, 16,
                9, 18, 17,
                11, 27, 24,
                10, 28, 25,
                9, 26, 22,
                8, 23, 27,
                14, 24, 28,
                13, 25, 22,
                12, 26, 23
            }
            ;

        int[] EDGE_POINTER =  {
                0, 3, 6, 9, 12, 15, 18, 21, 24, 27,
                30, 33, 36, 39, 42, 45, 48, 51, 54, 57,
                60, 63, 66, 69, 72, 75, 78, 81, 84
            }
            ;

        double[] XY =  {
                0.412, 0.984,
                0.494, 0.984,
                0.366, 0.926,
                0.388, 0.862,
                0.546, 0.926,
                0.518, 0.860,
                0.458, 0.818,
                0.152, 0.684,
                0.264, 0.682,
                0.354, 0.680,
                0.458, 0.670,
                0.554, 0.672,
                0.658, 0.668,
                0.774, 0.692,
                0.164, 0.450,
                0.228, 0.448,
                0.274, 0.390,
                0.242, 0.330,
                0.194, 0.278,
                0.146, 0.328,
                0.102, 0.390,
                0.668, 0.472,
                0.638, 0.416,
                0.656, 0.334,
                0.714, 0.270,
                0.798, 0.326,
                0.830, 0.408,
                0.754, 0.466
            }
            ;

        typeMethods.i4vec_copy(edge_num, EDGE_DATA, ref edge_data);
        typeMethods.i4vec_copy(node_num + 1, EDGE_POINTER, ref edge_pointer);
        typeMethods.r8vec_copy(2 * node_num, XY, ref xy);
    }

    public static void grf_example_size(ref int node_num, ref int edge_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRF_EXAMPLE_SIZE sizes a GRF example.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 January 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Stephen Skiena,
        //    Implementing Discrete Mathematics,
        //    Combinatorics and Graph Theory with Mathematica,
        //    Addison-Wesley, 1990.
        //
        //  Parameters:
        //
        //    Output, int *NODE_NUM, the number of nodes.
        //
        //    Output, int *EDGE_NUM, the number of edges.
        //
    {
        node_num = 28;
        edge_num = 84;
    }

    public static void grf_header_print(int node_num, int edge_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRF_HEADER_PRINT prints the header of a GRF file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 January 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Stephen Skiena,
        //    Implementing Discrete Mathematics,
        //    Combinatorics and Graph Theory with Mathematica,
        //    Addison-Wesley, 1990.
        //
        //  Parameters:
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, int EDGE_NUM, the number of edges.
        //
    {
        Console.WriteLine("");
        Console.WriteLine("  The number of nodes NODE_NUM = " + node_num + "");
        Console.WriteLine("  The number of edges EDGE_NUM = " + edge_num + "");
    }

    public static void grf_header_read(string input_filename, ref int node_num, ref int edge_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRF_HEADER_READ reads the header of a GRF file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 January 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Stephen Skiena,
        //    Implementing Discrete Mathematics,
        //    Combinatorics and Graph Theory with Mathematica,
        //    Addison-Wesley, 1990.
        //
        //  Parameters:
        //
        //    Input, string INPUT_FILENAME, the name of the file.
        //
        //    Output, int *NODE_NUM, the number of nodes.
        //
        //    Output, int *EDGE_NUM, the number of edges.
        //
    {
        int n;
        string[] input;

        edge_num = 0;
        node_num = 0;

        try
        {
            input = File.ReadAllLines(input_filename);
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("GRF_HEADER_READ - Fatal error!");
            Console.WriteLine("  Cannot open the input file \"" + input_filename + "\".");
            return;
        }

        foreach (string text in input)
        {
            if (text[0] == '#' || typeMethods.s_len_trim(text) == 0)
            {
                continue;
            }

            string[] tokens = text.Split(' ');
                
            n = tokens.Length;

            if (n < 3)
            {
                Console.WriteLine("");
                Console.WriteLine("GRF_HEADER_READ - Fatal error!");
                Console.WriteLine("  Illegal record has less than 3 data items");
                break;
            }

            edge_num = edge_num + n - 3;
            node_num += 1;
        }
    }

    public static void grf_header_write(string output_filename, ref List<string> output_unit,
            int node_num, int edge_num )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRF_HEADER_WRITE writes the header of a GRF file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 January 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string OUTPUT_FILENAME, the name of the output file.
        //
        //    Input, ofstream &OUTPUT_UNIT, the output file unit number.
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, int EDGE_NUM, the number of edges.
        //
    {
        //
        //  Write the header.
        //
        output_unit.Add("#  " + output_filename + "");
        output_unit.Add("#  created by grf_io::grf_header_write.C");
        output_unit.Add("#");
        output_unit.Add("#  Number of nodes  = " + node_num + "");
        output_unit.Add("#  Number of edges =  " + edge_num + "");
        output_unit.Add("#");
    }

    public static void grf_write(string output_filename, int node_num, int edge_num,
            int[] edge_pointer, int[] edge_data, double[] xy )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GRF_WRITE writes a GRF file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 January 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string OUTPUT_FILENAME, the name of the file
        //    to which the data should be written.
        //
        //    Input, int NODE_NUM, the number of nodes.
        //
        //    Input, int EDGE_NUM, the number of edges.
        //
        //    Input, int EDGE_POINTER[NODE_NUM+1], pointers to the
        //    first edge item for each node.
        //
        //    Input, int EDGE_DATA[EDGE_NUM], indices of adjacent nodes.
        //
        //    Input, double XY[2*NODE_NUM], the node coordinates.
        //
    {
        List<string> output_unit = new();

        switch (false)
        {
            //
            //  Write the header.
            //
            case true:
                grf_header_write(output_filename, ref output_unit, node_num, edge_num);
                break;
        }

        //
        //  Write the data.
        //
        grf_data_write(ref output_unit, node_num, edge_num, edge_pointer,
            edge_data, xy);
            
        //
        //  Open the output file.
        //
        try
        {
            File.WriteAllLines(output_filename, output_unit);
        }
        catch
        {
            Console.WriteLine("");
            Console.WriteLine("GRF_WRITE - Fatal error!");
            Console.WriteLine("  Cannot open the output file \"" + output_filename + "\".");
        }
    }
}