using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.Quadrature;
using Burkardt.Types;

namespace Burkardt.TriangleNS;

public static class SymmetricQuadrature
{
    public static void triasymq(int n, double[] vert1, double[] vert2, double[] vert3,
            ref double[] rnodes, ref double[] weights, int numnodes)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIASYMQ returns a symmetric quadrature formula for a user triangle.
        //
        //  Discussion:
        //
        //    This routine constructs (or rather retrieves) a symmetric
        //    quadrature formula on the user-defined triangle in the plane
        //    specified by its vertices (vert1,vert2,vert3).
        //
        //    The total number of nodes and weights to be created is numnodes.
        //
        //      n       1     2     3     4     5     6     7     8     9    10
        //     -----------------------------------------------------------------
        //    nodes     1     3     6     6     7    12    15    16    19    25
        //
        //
        //      n      11    12    13    14    15    16    17    18    19    20
        //     -----------------------------------------------------------------
        //    nodes    28    33    37    42    49    55    60    67    73    79
        //
        //
        //      n      21    22    23    24    25    26    27    28    29    30
        //     -----------------------------------------------------------------
        //    nodes    87    96   103   112   120   130   141   150   159   171
        //
        //
        //      n      31    32    33    34    35    36    37    38    39    40
        //     -----------------------------------------------------------------
        //    nodes   181   193   204   214   228   243   252   267   282   295
        //
        //
        //      n      41    42    43    44    45    46    47    48    49    50
        //     -----------------------------------------------------------------
        //    nodes   309   324   339   354   370   385   399   423   435   453
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU GPL license.
        //
        //  Modified:
        //
        //    30 June 2014
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Hong Xiao, Zydrunas Gimbutas,
        //    A numerical algorithm for the construction of efficient quadrature
        //    rules in two and higher dimensions,
        //    Computers and Mathematics with Applications,
        //    Volume 59, 2010, pages 663-676.
        //
        //  Parameters:
        //
        //    Input, int N, the degree of the quadrature (must not exceed 50).
        //    Note that the total number of nodes to be created is numnodes.
        //
        //    Input, double VERT1[2], VERT2[2], VERT3[2], the vertices of
        //    the triangle on which the quadrature rule is to be constructed.
        //
        //    Output, double RNODES[2*NUMNODES], the nodes in the plane
        //    (all inside the user-supplied triangle)
        //
        //    Output, double WEIGHTS[NUMNODES], the quadrature weights.
        //
        //    Input, int NUMNODES, the number of nodes.
        //
    {
        int itype = 0;

        Quae.quaequad(itype, n, ref rnodes, ref weights, numnodes);

        typeMethods.trianmap(numnodes, vert1, vert2, vert3, ref rnodes, ref weights);
    }

    public static void triasymq_gnuplot(double[] vert1, double[] vert2, double[] vert3,
            int numnodes, double[] rnodes, string header)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIASYMQ_GNUPLOT: set up a GNUPLOT plot of the triangle quadrature rule.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU GPL license.
        //
        //  Modified:
        //
        //    11 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Hong Xiao, Zydrunas Gimbutas,
        //    A numerical algorithm for the construction of efficient quadrature
        //    rules in two and higher dimensions,
        //    Computers and Mathematics with Applications,
        //    Volume 59, 2010, pages 663-676.
        //
        //  Parameters:
        //
        //    Input, double VERT1[2], VERT2[2], VERT3[2], the vertices
        //    of the triangle.
        //
        //    Input, double NUMNODES, the number of nodes.
        //
        //    Input, double RNODES[2*NUMNODES], the coordinates of
        //    the nodes.
        //
        //    Input, string HEADER, a string to be used to identify
        //    the files created.
        //
    {
        List<string> command_unit = new();
        int j;
        List<string> node_unit = new();
        List<string> vertex_unit = new();
        //
        //  Create the vertex file.
        //
        string vertex_filename = header + "_vertices.txt";
        vertex_unit.Add(vert1[0] + "  "
                                 + vert1[1] + "");
        vertex_unit.Add(vert2[0] + "  "
                                 + vert2[1] + "");
        vertex_unit.Add(vert3[0] + "  "
                                 + vert3[1] + "");
        vertex_unit.Add(vert1[0] + "  "
                                 + vert1[1] + "");
        File.WriteAllLines(vertex_filename, vertex_unit);

        Console.WriteLine("");
        Console.WriteLine("  Created vertex file '" + vertex_filename + "'");
        //
        //  Create node file.
        //
        string node_filename = header + "_nodes.txt";
        for (j = 0; j < numnodes; j++)
        {
            node_unit.Add(rnodes[0 + j * 2] + "  "
                                            + rnodes[1 + j * 2] + "");
        }

        File.WriteAllLines(node_filename, node_unit);
        Console.WriteLine("  Created node file '" + node_filename + "'");
        //
        //  Create graphics command file.
        //
        string command_filename = header + "_commands.txt";
        command_unit.Add("# " + command_filename + "");
        command_unit.Add("#");
        command_unit.Add("# Usage:");
        command_unit.Add("#  gnuplot < " + command_filename + "");
        command_unit.Add("#");
        command_unit.Add("set term png");
        string plot_filename = header + ".png";
        command_unit.Add("set output '" + plot_filename + "'");
        command_unit.Add("set xlabel '<--- X --->'");
        command_unit.Add("set ylabel '<--- Y --->'");
        command_unit.Add("set title '" + header + "'");
        command_unit.Add("set grid");
        command_unit.Add("set key off");
        command_unit.Add("set size ratio -1");
        command_unit.Add("set style data lines");
        command_unit.Add("set timestamp");
        command_unit.Add("plot '" + vertex_filename + "' with lines lw 3, \\");
        command_unit.Add("     '" + node_filename + "' with points pt 7 lt 0");
        File.WriteAllLines(command_filename, command_unit);

        Console.WriteLine("  Created command file '" + command_filename + "'");

    }
}