using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.Types;

namespace Burkardt.StochasticDifferentialEquations;

public static class ChainRule
{
    public static void chain(ref typeMethods.r8vecNormalData data, ref int seed, int n, ref double[] xem, ref double[] vem, ref double diff)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHAIN tests the stochastic Chain Rule.
        //
        //  Discussion:
        //
        //    This function solves a stochastic differential equation for 
        //
        //      V = sqrt(X) 
        //
        //    where X satisfies the stochastic differential equation:
        // 
        //      dX = ( alpha - X ) * dt + beta * sqrt(X) dW,
        //      X(0) = Xzero,
        //
        //    with 
        //
        //      alpha = 2,
        //      beta = 1,
        //      Xzero = 1.
        //
        //    From the stochastic Chain Rule, the SDE for V is therefore:
        //
        //      dV = ( ( 4 * alpha - beta^2 ) / ( 8 * V ) - 1/2 V ) dt + 1/2 beta dW
        //      V(0) = sqrt ( Xzero ).
        //
        //    Xem is the Euler-Maruyama solution for X. 
        //
        //    Vem is the Euler-Maruyama solution of the SDE for V from
        //    the stochastic Chain Rule.
        //
        //    Hence, we compare sqrt(Xem) and Vem.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 September 2012
        //
        //  Author:
        //
        //    Original Matlab version by Desmond Higham.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Desmond Higham,
        //    An Algorithmic Introduction to Numerical Simulation of
        //    Stochastic Differential Equations,
        //    SIAM Review,
        //    Volume 43, Number 3, September 2001, pages 525-546.
        //
        //  Parameters:
        //
        //    Input, int &SEED, a seed for the random number generator.
        //
        //    Input, int N, the number of time steps.
        //
        //    Output, double XEM[N+1], the computed value of X.
        //
        //    Output, double VEM[N+1], the computed value of V.
        //
        //    Output, double &DIFF, the maximum value of |sqrt(XEM)-V|.
        //
    {
        double alpha;
        double beta;
        double dt;
        double dt2;
        double[] dw;
        int i;
        int j;
        double tmax;
        //
        //  Set problem parameters.
        //
        alpha = 2.0;
        beta = 1.0;
        //
        //  Stepping parameters.
        //  dt2 is the size of the Euler-Maruyama steps.
        //
        tmax = 1.0;
        dt = tmax / n;
        dt2 = dt;
        //
        //  Define the increments dW.
        //
        dw = typeMethods.r8vec_normal_01_new(n, ref data, ref seed);

        for (i = 0; i < n; i++)
        {
            dw[i] = Math.Sqrt(dt) * dw[i];
        }

        //
        //  Solve for X(t).
        //
        xem[0] = 1.0;
        for (j = 1; j <= n; j++)
        {
            xem[j] = xem[j - 1] + (alpha - xem[j - 1]) * dt2
                                + beta * Math.Sqrt(xem[j - 1]) * dw[j - 1];
        }

        //
        //  Solve for V(t).
        //
        vem[0] = Math.Sqrt(xem[0]);
        for (j = 1; j <= n; j++)
        {
            vem[j] = vem[j - 1]
                     + ((4.0 * alpha - beta * beta) / (8.0 * vem[j - 1])
                        - 0.5 * vem[j - 1]) * dt2
                     + 0.5 * beta * dw[j - 1];
        }

        //
        //  Compare sqrt(X) and V.
        //
        diff = 0.0;
        for (i = 0; i <= n; i++)
        {
            diff = Math.Max(diff, Math.Abs(Math.Sqrt(xem[i]) - vem[i]));
        }
    }

    public static void chain_gnuplot(int n, double[] x, double[] v)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHAIN_GNUPLOT writes a GNUPLOT input file to plot CHAIN data.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 September 2012
        //
        //  Author:
        //
        //    John Burkardt.
        //
        //  Reference:
        //
        //    Desmond Higham,
        //    An Algorithmic Introduction to Numerical Simulation of
        //    Stochastic Differential Equations,
        //    SIAM Review,
        //    Volume 43, Number 3, September 2001, pages 525-546.
        //
        //  Parameters:
        //
        //    Input, int N, the number of steps.
        //
        //    Input, double X[N+1], the value of X.
        //
        //    Input, double V[N+1], the value of V.
        //
    {
        string command_filename = "chain_commands.txt";
        List<string> command = new();
        string data_filename = "chain_data.txt";
        List<string> data = new();
        int i;
        double t;

        for (i = 0; i <= n; i++)
        {
            t = i / (double) n;
            data.Add("  " + t
                          + "  " + Math.Sqrt(x[i])
                          + "  " + v[i] + "");
        }

        File.WriteAllLines(data_filename, data);

        Console.WriteLine("");
        Console.WriteLine("  CHAIN data stored in \"" + data_filename + "\".");

        command.Add("# chain_commands.txt");
        command.Add("# created by sde::chain_gnuplot.");
        command.Add("#");
        command.Add("# Usage:");
        command.Add("#  gnuplot < chain_commands.txt");
        command.Add("#");
        command.Add("set term png");
        command.Add("set output 'chain.png'");
        command.Add("set xlabel 't'");
        command.Add("set ylabel 'Sqrt(X(t)) vs V(X(t))'");
        command.Add("set title 'V(X(t)) from X(t) and from Chain Rule'");
        command.Add("set grid");
        command.Add("set style data lines");
        command.Add("plot 'chain_data.txt' using 1:2 title 'Sqrt(X(t))', \\");
        command.Add("     'chain_data.txt' using 1:3 title 'V(X(t))'");
        command.Add("quit");

        File.WriteAllLines(command_filename, command);

        Console.WriteLine("  CHAIN plot commands stored in \"" + command_filename + "\".");

    }
}