using System;
using System.Collections.Generic;
using System.IO;
using Burkardt;
using Burkardt.MatrixNS;
using Burkardt.Sequence;
using Burkardt.Types;

namespace JacobiTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for JACOBI_TEST.
        //
        //  Discussion:
        //
        //    JACOBI_TEST tests the JACOBI library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 March 2020
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("JACOBI_TEST");
        Console.WriteLine("  Test JACOBI.");

        jacobi_test01();
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("JACOBI_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void jacobi_test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    JACOBI_TEST01 tests JACOBI1.
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
    {
        double[] a;
        double[] b;
        string command_filename;
        List<string> command = new();
        string data_filename;
        List<string> data = new();
        int i;
        int it;
        int it_num;
        int j;
        double[] m_plot;
        double[] ml_plot;
        int n;
        double[] r_plot;
        double[] rl_plot;
        double[] s_plot;
        double t;
        double w = 0.5;
        double[] x;
        double[] x_exact;
        double[] x_new;
        double[] x_plot;

        Console.WriteLine("");
        Console.WriteLine("JACOBI_TEST01:");

        it_num = 2000;
        n = 33;
        //
        //  Set the matrix A.
        //
        a = Matrix.dif2(n, n);
        //
        //  Determine the right hand side vector B.
        //
        x_exact = new double[n];
        for (i = 0; i < n; i++)
        {
            t = i / (double) (n - 1);
            x_exact[i] = Math.Exp(t) * (t - 1) * t;
            //   x_exact[i] = ( double ) ( i + 1 );
        }

        b = typeMethods.r8mat_mv_new(n, n, a, x_exact);
        //
        //  Set the initial estimate for the solution.
        //
        x = new double[n];
        for (i = 0; i < n; i++)
        {
            x[i] = 0.0;
        }

        //
        //  Allocate plot arrays.
        //
        m_plot = new double[it_num + 1];
        r_plot = new double[it_num + 1];
        s_plot = new double[it_num + 1];
        x_plot = new double[n * (it_num + 1)];
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
            x_new = Jacobi.jacobi1(n, a, b, x);

            r_plot[it] = typeMethods.r8mat_residual_norm(n, n, a, x_new, b);
            //
            //  Compute the average point motion.
            //
            m_plot[it] = typeMethods.r8vec_diff_norm_squared(n, x, x_new) / n;
            //
            //  Update the solution
            //
            for (i = 0; i < n; i++)
            {
                x[i] = (1.0 - w) * x[i] + w * x_new[i];
            }

            //  r8vec_copy ( n, x_new, x );

            for (i = 0; i < n; i++)
            {
                x_plot[i + 0 * n] = x[i];
            }
        }

        typeMethods.r8vec_print(n, x, "Solution");
        //
        //  Plot the residual.
        //
        rl_plot = new double[it_num + 1];
        for (j = 0; j <= it_num; j++)
        {
            rl_plot[j] = Math.Log(r_plot[j]);
        }

        //
        //  Create the data file.
        //
        data_filename = "residual_data.txt";
            
        for (j = 0; j <= it_num; j++)
        {
            data.Add(j + "  "
                       + rl_plot[j] + "");
        }

        File.WriteAllLines(data_filename, data);
            
        Console.WriteLine(" ");
        Console.WriteLine("  Data stored in \"" + data_filename + "\".");
        //
        //  Create the command file.
        //
        command_filename = "residual_commands.txt";

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
            
        Console.WriteLine("  Plot commands stored in \"" + command_filename + "\".");
        //
        //  Plot the average point motion.
        //
        ml_plot = new double[it_num + 1];
        for (j = 0; j <= it_num; j++)
        {
            ml_plot[j] = Math.Log(m_plot[j]);
        }

        //
        //  Create the data file.
        //
        data_filename = "motion_data.txt";
        data.Clear();

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
        command.Clear();

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
            
        //
        //  Plot the evolution of the locations of the generators.
        //
        //figure ( 3 )

        //y = ( 0 : it_num );
        //for k = 1 : n
        //  plot ( x_plot(k,1:it_num+1), y )
        //  hold on;
        //end
        //grid on
        //hold off;

        //title ( "Generator evolution." );
        //xlabel ( "Generator positions" );
        //ylabel ( "Iterations" ); 
    }
}