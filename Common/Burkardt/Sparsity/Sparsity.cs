using System;
using System.Collections.Generic;
using System.IO;
using System.Net;
using Burkardt.Types;

namespace Burkardt.SparsityNS;

public static class Sparsity
{
    public static void spy_file(string header, string data_filename)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPY_FILE plots a sparsity pattern stored in a file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 September 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string HEADER, the name to be used for the
        //    title of the plot, and as part of the names of the command
        //    and plot files.
        //
        //    Input, string DATA_FILENAME, the name of the file
        //    containing the indices of nonzero matrix entries.
        //
    {
        int n0 = +typeMethods.i4_huge();
        int n1 = -typeMethods.i4_huge();
        int m0 = +typeMethods.i4_huge();
        int m1 = -typeMethods.i4_huge();
        int nz_num = 0;

        string[] data_unit = File.ReadAllLines(data_filename);
            
        foreach (string line in data_unit)
        {
            i4vec vals = typeMethods.s_to_i4vec(line, 2);
            int i = vals.ivec[0];
            int j = vals.ivec[1];
                
            nz_num += 1;
            m0 = Math.Min(m0, i);
            m1 = Math.Max(m1, i);
            n0 = Math.Min(n0, j);
            n1 = Math.Max(n1, j);
        }

        //
        //  Create command file.
        //
        string command_filename = header + "_commands.txt";

        List<string> command_unit = new()
        {
            "# " + command_filename,
            "#",
            "# Usage:",
            "#  gnuplot < " + command_filename + "",
            "#",
            "unset key",
            "set term png"
        };

        string png_filename = header + ".png";
        command_unit.Add("set output '" + png_filename + "'");
        command_unit.Add("set size ratio -1");
        command_unit.Add("set xlabel '<--- J --->'");
        command_unit.Add("set ylabel '<--- I --->'");

        command_unit.Add("set title '"
                         + nz_num + " nonzeros for \""
                         + header + "\"'");
        command_unit.Add("set timestamp");
        command_unit.Add("plot [y="
                         + m0 + ":"
                         + m1 + "] [x="
                         + n0 + ":"
                         + n1 + "] '"
                         + data_filename + "' with points pt 5");

        File.WriteAllLines(command_filename, command_unit);
        Console.WriteLine("  Created graphics command file '" + command_filename + "'");
    }

    public static void spy_ge(int m, int n, double[] a, string header )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPY_GE plots a sparsity pattern for a general storage (GE) matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 September 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns
        //    in the matrix.
        //
        //    Input, double A[M*N], the matrix.
        //
        //    Input, string HEADER, the name to be used for the
        //    title of the plot, and as part of the names of the data, command
        //    and plot files.
        //
    {
        //
        //  Create data file.
        //
        string data_filename = header + "_data.txt";

        List<string> data_unit = new();
            
        int nz_num = 0;
        for (int j = 0; j < n; j++)
        {
            for (int i = 0; i < m; i++)
            {
                if (a[i + j * m] != 0.0)
                {
                    data_unit.Add(j.ToString() + "  " + i.ToString());
                    nz_num += 1;
                }
            }
        }
            
        File.WriteAllLines(data_filename, data_unit);

        Console.WriteLine("");
        Console.WriteLine("  Created sparsity data file '" + data_filename + "'");
        //
        //  Create command file.
        //
        string command_filename = header + "_commands.txt";

        List<string> command_unit = new()
        {
            "# " + command_filename + "",
            "#",
            "# Usage:",
            "#  gnuplot < " + command_filename + "",
            "#",
            "unset key",
            "set term png"
        };

        string png_filename = header + ".png";
        command_unit.Add("set output '" + png_filename + "'");
        command_unit.Add("set size ratio -1");
        command_unit.Add("set xlabel '<--- J --->'");
        command_unit.Add("set ylabel '<--- I --->'");
        command_unit.Add("set title '" + nz_num + " nonzeros for \""
                         + header + "\"'");
        command_unit.Add("set timestamp");
        command_unit.Add("plot [y=0:" + (n - 1).ToString() + "] [x="
                         + (m - 1).ToString() + ":0] '"
                         + data_filename + "' with points pt 5");

        File.WriteAllLines(command_filename, command_unit);
            
        Console.WriteLine("  Created graphics command file '" + command_filename + "'");
    }
}