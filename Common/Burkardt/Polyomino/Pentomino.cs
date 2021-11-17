using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.Types;

namespace Burkardt.PolyominoNS;

public static class Pentomino
{
    public static void pentomino_matrix(string name, ref int p_m, ref int p_n, ref int[] p)

        //*****************************************************************************/
        //
        //  Purpose:
        //
        //    PENTOMINO_MATRIX returns a 0/1 matrix defining a particular pentomino.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 April 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, char NAME, is "f", "i", "l", "n", "p", "t", "u", "v", "w", "x", 
        //    "y" or "z".
        //
        //    Output, int &P_M, &P_N, the number of rows and columns of the
        //    representation.
        //
        //    Output, int P[P_M*P_N], a matrix of 0"s and 1"s that indicates
        //    the shape of the pentomino.
        //
    {
        int[] f_mat =
        {
            0, 1, 1,
            1, 1, 0,
            0, 1, 0
        };
        int[] i_mat =
        {
            1,
            1,
            1,
            1,
            1
        };
        int[] l_mat =
        {
            1, 0,
            1, 0,
            1, 0,
            1, 1
        };
        int[] n_mat =
        {
            1, 1, 0, 0,
            0, 1, 1, 1
        };
        int[] p_mat =
        {
            1, 1,
            1, 1,
            1, 0
        };
        int[] t_mat =
        {
            1, 1, 1,
            0, 1, 0,
            0, 1, 0
        };
        int[] u_mat =
        {
            1, 0, 1,
            1, 1, 1
        };
        int[] v_mat =
        {
            1, 0, 0,
            1, 0, 0,
            1, 1, 1
        };
        int[] w_mat =
        {
            1, 0, 0,
            1, 1, 0,
            0, 1, 1
        };
        int[] x_mat =
        {
            0, 1, 0,
            1, 1, 1,
            0, 1, 0
        };
        int[] y_mat =
        {
            0, 0, 1, 0,
            1, 1, 1, 1
        };
        int[] z_mat =
        {
            1, 1, 0,
            0, 1, 0,
            0, 1, 1
        };

        switch (name)
        {
            case "f":
            case "F":
                p_m = 3;
                p_n = 3;
                p = typeMethods.i4rows_copy_new(p_m, p_n, f_mat);
                break;
            case "i":
            case "I":
                p_m = 5;
                p_n = 1;
                p = typeMethods.i4rows_copy_new(p_m, p_n, i_mat);
                break;
            case "l":
            case "L":
                p_m = 4;
                p_n = 2;
                p = typeMethods.i4rows_copy_new(p_m, p_n, l_mat);
                break;
            case "n":
            case "N":
                p_m = 2;
                p_n = 4;
                p = typeMethods.i4rows_copy_new(p_m, p_n, n_mat);
                break;
            case "p":
            case "P":
                p_m = 3;
                p_n = 2;
                p = typeMethods.i4rows_copy_new(p_m, p_n, p_mat);
                break;
            case "t":
            case "T":
                p_m = 3;
                p_n = 3;
                p = typeMethods.i4rows_copy_new(p_m, p_n, t_mat);
                break;
            case "u":
            case "U":
                p_m = 2;
                p_n = 3;
                p = typeMethods.i4rows_copy_new(p_m, p_n, u_mat);
                break;
            case "v":
            case "V":
                p_m = 3;
                p_n = 3;
                p = typeMethods.i4rows_copy_new(p_m, p_n, v_mat);
                break;
            case "w":
            case "W":
                p_m = 3;
                p_n = 3;
                p = typeMethods.i4rows_copy_new(p_m, p_n, w_mat);
                break;
            case "x":
            case "X":
                p_m = 3;
                p_n = 3;
                p = typeMethods.i4rows_copy_new(p_m, p_n, x_mat);
                break;
            case "y":
            case "Y":
                p_m = 2;
                p_n = 4;
                p = typeMethods.i4rows_copy_new(p_m, p_n, y_mat);
                break;
            case "z":
            case "Z":
                p_m = 3;
                p_n = 3;
                p = typeMethods.i4rows_copy_new(p_m, p_n, z_mat);
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("PENTOMINO_MATRIX - Fatal error!");
                Console.WriteLine("  Illegal name = '" + name + "'");
                Console.WriteLine("  Legal names: f, i, l, n, p, t, u, v, w, x, y, z.");
                break;
        }
    }

    public static void pentomino_plot(int p_m, int p_n, int[] p, string label)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PENTOMINO_PLOT plots a particular pentomino in a 5x5 grid.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 April 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int P_M, P_N, the number of rows and columns of the
        //    representation.
        //
        //    Input, int P[P_M*P_N], a matrix of 0"s and 1"s.
        //    1 <= P_M, P_N <= 5.  There should be exactly 5 values of one.
        //
        //    Input, string LABEL, a title for the plot.
        //
    {
        string color = "";
        int[] color_index;
        string command_filename;
        List<string> command_unit = new();
        int i;
        int i_reverse;
        int j;
        int k;
        int m = 5;
        int n = 5;
        string plot_filename;

        command_filename = label + "_commands.txt";
        plot_filename = label + ".png";
        //
        //  Initially, the grid is entirely white (color 0)
        //
        color_index = typeMethods.i4rows_zeros_new(m, n);
        //
        //  Place the pentomino on the grid, so that it is "snug" in the upper left corner.
        //
        for (i = 0; i < p_m; i++)
        {
            for (j = 0; j < p_n; j++)
            {
                color_index[i * n + j] = p[i * p_n + j];
            }
        }

        //
        //  Create the command file.
        //
        command_unit.Add("# " + command_filename + "");
        command_unit.Add("#");
        command_unit.Add("# Usage:");
        command_unit.Add("#  gnuplot < " + command_filename + "");
        command_unit.Add("#");
        command_unit.Add("set term png");
        command_unit.Add("set output '" + plot_filename + "'");
        command_unit.Add("set title '" + label + "'");
        //
        //  Get a plot of TRUE SQUARES.
        //
        command_unit.Add("set xrange [ 0 : 5 ]");
        command_unit.Add("set yrange [ 0 : 5 ]");
        command_unit.Add("set size square");
        command_unit.Add("unset border");
        command_unit.Add("unset tics");
        command_unit.Add("set nokey");

        k = 0;
        for (i = 0; i < m; i++)
        {
            for (j = 0; j < n; j++)
            {
                k += 1;

                color = color_index[i * n + j] switch
                {
                    0 => "white",
                    1 => "black",
                    _ => color
                };

                i_reverse = m - 1 - i;
                command_unit.Add("set object " + k
                                               + " rect from " + j
                                               + ", " + i_reverse
                                               + " to " + (j + 1)
                                               + ", " + (i_reverse + 1) + " back");
                command_unit.Add("set object " + k
                                               + " rect fc rgb '" + color
                                               + "' fillstyle solid 1.0");
            }
        }

        //
        //  If you don"t have some bogus PLOT command here, all the previous work
        //  results in no plot all.  Way to go, gnuplot!
        //  Here, we plot the function y = -1, which is out of range and won"t show up.
        //
        command_unit.Add("plot -1 with lines");
        File.WriteAllLines(command_filename, command_unit);
        Console.WriteLine("  PENTOMINO_PLOT created command file '" + command_filename + "'");
    }
}