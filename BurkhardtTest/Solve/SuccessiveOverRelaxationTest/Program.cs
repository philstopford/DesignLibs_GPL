﻿using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.MatrixNS;
using Burkardt.SolveNS;
using Burkardt.Types;

namespace SuccessiveOverRelaxationTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for SOR_TEST.
        //
        //  Discussion:
        //
        //    SOR_TEST tests the SOR library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 May 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("SOR_TEST");
            
        Console.WriteLine("  Test the SOR library.");

        double w = 0.5;
        sor_test01(w);

        w = 1.0;
        sor_test01(w);

        w = 1.5;
        sor_test01(w);

        Console.WriteLine("");
        Console.WriteLine("SOR_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void sor_test01(double w)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SOR_TEST01 tests SOR1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 June 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double W, the relaxation factor.
        //    0 < W < 2 is required.
        //
    {
        double[] b = new double[1];
        List<string> command = new();
        List<string> data = new();
        int i;
        int it;
        int j;

        Console.WriteLine("");
        Console.WriteLine("SOR1_TEST01:");
        Console.WriteLine("  Relaxation parameter W = " + w + "");

        const int it_num = 2000;
        const int n = 33;
        //
        //  Set the matrix A.
        //
        double[] a = Matrix.dif2(n, n);
        //
        //  Determine the right hand side vector B.
        //
        double[] x_exact = new double[n];
        for (i = 0; i < n; i++)
        {
            double t = i / (double)(n - 1);
            x_exact[i] = Math.Exp(t) * (t - 1) * t;
            //   x_exact[i] = ( double ) ( i + 1 );
        }

        typeMethods.r8mat_mv(n, n, a, x_exact, ref b);
        //
        //  Set the initial estimate for the solution.
        //
        double[] x = new double[n];
        for (i = 0; i < n; i++)
        {
            x[i] = 0.0;
        }

        //
        //  Allocate plot arrays.
        //
        double[] m_plot = new double[it_num + 1];
        double[] r_plot = new double[it_num + 1];
        double[] s_plot = new double[it_num + 1];
        double[] x_plot = new double[n * (it_num + 1)];
        //
        //  Initialize plot arrays.
        //
        r_plot[0] = typeMethods.r8mat_residual_norm(n, n, a, x, b);
        m_plot[0] = 1.0;
        for (i = 0; i < n; i++)
        {
            x_plot[i + 0 * n] = x[i];
        }

        for (j = 0; j <= it_num; j++)
        {
            s_plot[j] = j;
        }

        //
        //  Carry out the iteration.
        //
        for (it = 1; it <= it_num; it++)
        {
            double[] x_new = SuccessiveOverRelaxation.sor1(n, a, b, x, w);

            r_plot[it] = typeMethods.r8mat_residual_norm(n, n, a, x_new, b);
            //
            //  Compute the average point motion.
            //
            m_plot[it] = typeMethods.r8vec_diff_norm_squared(n, x, x_new) / n;
            //
            //  Update the solution
            //
            typeMethods.r8vec_copy(n, x_new, ref x);

            for (i = 0; i < n; i++)
            {
                x_plot[i + 0 * n] = x[i];
            }

        }

        typeMethods.r8vec_print(n, x, "Solution");
        //
        //  Plot the residual.
        //
        double[] rl_plot = new double[it_num + 1];
        for (j = 0; j <= it_num; j++)
        {
            rl_plot[j] = Math.Log(r_plot[j]);
        }

        //
        //  Create the data file.
        //
        string data_filename = "residual_data.txt";

        for (j = 0; j <= it_num; j++)
        {
            data.Add(j + "  "
                       + rl_plot[j] + "");
        }

        File.WriteAllLines(data_filename, data);
        data.Clear();

        Console.WriteLine(" ");
        Console.WriteLine("  Data stored in \"" + data_filename + "\".");
        //
        //  Create the command file.
        //
        string command_filename = "residual_commands.txt";

        command.Add("# residual_commands.txt");
        command.Add("#");
        command.Add("# Usage:");
        command.Add("#  gnuplot < residual_commands.txt");
        command.Add("#");
        command.Add("set term png");
        command.Add("set output 'residual.png'");
        command.Add("set style data lines");
        command.Add("set xlabel 'Iteration'");
        command.Add("set ylabel 'Residual'");
        command.Add("set title 'Log(Residual) over Iterations'");
        command.Add("set grid");
        command.Add("plot 'residual_data.txt' using 1:2 lw 2");
        command.Add("quit");

        File.WriteAllLines(command_filename, command);
        command.Clear();

        Console.WriteLine("  Plot commands stored in \"" + command_filename + "\".");
        //
        //  Plot the average point motion.
        //
        double[] ml_plot = new double[it_num + 1];
        for (j = 0; j <= it_num; j++)
        {
            ml_plot[j] = Math.Log(m_plot[j]);
        }

        //
        //  Create the data file.
        //
        data_filename = "motion_data.txt";

        for (j = 0; j <= it_num; j++)
        {
            data.Add(j + "  "
                       + ml_plot[j] + "");
        }

        File.WriteAllLines(data_filename, data);

        Console.WriteLine(" ");
        Console.WriteLine("  Data stored in \"" + data_filename + "\".");
        //
        //  Create the command file.
        //
        command_filename = "motion_commands.txt";

        command.Add("# motion_commands.txt");
        command.Add("#");
        command.Add("# Usage:");
        command.Add("#  gnuplot < motion_commands.txt");
        command.Add("#");
        command.Add("set term png");
        command.Add("set output 'motion.png'");
        command.Add("set style data lines");
        command.Add("set xlabel 'Iteration'");
        command.Add("set ylabel 'Motion'");
        command.Add("set title 'Log(Motion) over Iterations'");
        command.Add("set grid");
        command.Add("plot 'motion_data.txt' using 1:2 lw 2");
        command.Add("quit");

        File.WriteAllLines(command_filename, command);

        Console.WriteLine("  Plot commands stored in \"" + command_filename + "\".");

    }
}