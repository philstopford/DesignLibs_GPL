using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.Types;

namespace Burkardt.StochasticDifferentialEquations;

public static class EulerMaruyama
{
    public static void em(ref typeMethods.r8vecNormalData data, ref int seed, int n, ref double[] t, ref double[] xtrue, ref double[] t2,
            ref double[] xem, ref double error)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EM applies the Euler-Maruyama method to a linear SDE.
        //
        //  Discussion:
        //
        //    The SDE is 
        //
        //      dX = lambda * X dt + mu * X dW,   
        //      X(0) = Xzero,
        //
        //    where 
        //
        //      lambda = 2,
        //      mu = 1,
        //      Xzero = 1.
        //
        //    The discretized Brownian path over [0,1] uses
        //    a stepsize dt = 2^(-8).
        //
        //    The Euler-Maruyama method uses a larger timestep Dt = R*dt,
        //    where R is an integer.  For an SDE of the form
        //
        //      dX = f(X(t)) dt + g(X(t)) dW(t)
        //
        //    it has the form
        //
        //      X(j) = X(j-1) + f(X(j-1)) * Dt + g(X(j-1)) * ( W(j*Dt) - W((j-1)*Dt) )
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
        //    Volume 43, Number 3, September 2001, pages 525-546
        //
        //  Parameters:
        //
        //    Input/output, int &SEED, a seed for the random 
        //    number generator.
        //
        //    Input, int N, the number of time steps.  A typical
        //    value is 2^8.  N should be a multiple of 4.
        //
        //    Output, double T[N+1], the time values for the exact solution.
        //
        //    Output, double XTRUE[N+1], the exact solution.
        //
        //    Output, double T2[N/4+1], the time values for the 
        //    Euler-Maruyama solution.
        //
        //    Output, double XEM[N/4+1], the Euler-Maruyama solution.
        //
        //    Output, double &ERROR, the value of | XEM(T) - XTRUE(T) |.
        //
    {
        int j;
        //
        //  Set problem parameters.
        //
        const double lambda = 2.0;
        const double mu = 1.0;
        const double xzero = 1.0;
        //
        //  Set stepping parameters.
        //
        const double tmax = 1.0;
        double dt = tmax / n;
        //
        //  Define the increments dW.
        //
        double[] dw = typeMethods.r8vec_normal_01_new(n, ref data, ref seed);

        for (j = 0; j < n; j++)
        {
            dw[j] = Math.Sqrt(dt) * dw[j];
        }

        //
        //  Sum the Brownian increments.
        //
        double[] w = new double[n + 1];
        w[0] = 0.0;
        for (j = 1; j <= n; j++)
        {
            w[j] = w[j - 1] + dw[j - 1];
        }

        for (j = 0; j <= n; j++)
        {
            t[j] = j * tmax / n;
        }

        //
        //  Compute the discretized Brownian path.
        //
        for (j = 0; j <= n; j++)
        {
            xtrue[j] = xzero * Math.Exp((lambda - 0.5 * mu * mu) * (t[j] + mu * w[j]));
        }

        //
        //  Set:
        //  R, the multiplier for the EM step, 
        //  Dt, the EM stepsize,
        //  L, the number of EM steps (we need N to be a multiple of R!)
        //
        const int r = 4;
        double dt2 = r * dt;
        int l = n / r;

        for (j = 0; j <= l; j++)
        {
            t2[j] = j * tmax / l;
        }

        //
        //  Compute XEM.
        //
        xem[0] = xzero;
        for (j = 1; j <= l; j++)
        {
            double dw2 = 0.0;
            int i;
            for (i = r * (j - 1); i < r * j; i++)
            {
                dw2 += dw[i];
            }

            xem[j] = xem[j - 1] + dt2 * lambda * xem[j - 1] + mu * xem[j - 1] * dw2;
        }

        error = Math.Abs(xem[l] - xtrue[n]);
    }

    public static void em_gnuplot(int n, double[] t, double[] xtrue, double[] t2, double[] xem)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EM_GNUPLOT writes a GNUPLOT input file to plot EM data.
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
        //    Input, double T[N+1], the time values for the exact solution.
        //
        //    Input, double XTRUE[N+1], the exact solution.
        //
        //    Input, double T2[N/4+1], the time values for the 
        //    Euler-Maruyama solution.
        //
        //    Input, double XEM[N/4+1], the Euler-Maruyama solution.
        //
    {
        const string command_filename = "em_commands.txt";
        List<string> command = new();
        const string data1_filename = "em1_data.txt";
        const string data2_filename = "em2_data.txt";
        List<string> data = new();
        int i;
        //
        //  Create data file #1.
        //

        for (i = 0; i <= n; i++)
        {
            data.Add("  " + t[i]
                          + "  " + xtrue[i] + "");
        }

        File.WriteAllLines(data1_filename, data);

        Console.WriteLine("");
        Console.WriteLine("  EM data #1 stored in \"" + data1_filename + "\".");

        data.Clear();
        for (i = 0; i <= n / 4; i++)
        {
            data.Add("  " + t2[i]
                          + "  " + xem[i] + "");
        }

        File.WriteAllLines(data2_filename, data);
        Console.WriteLine("");
        Console.WriteLine("  EM data #2 stored in \"" + data2_filename + "\".");

        command.Add("# em_commands.txt");
        command.Add("# created by sde::em_gnuplot.");
        command.Add("#");
        command.Add("# Usage:");
        command.Add("#  gnuplot < em_commands.txt");
        command.Add("#");
        command.Add("set term png");
        command.Add("set output 'em.png'");
        command.Add("set xlabel 't'");
        command.Add("set ylabel 'X(t)'");
        command.Add("set title 'Exact X(t) and Euler-Maruyama Estimate'");
        command.Add("set grid");
        command.Add("set style data lines");
        command.Add("plot 'em1_data.txt' using 1:2 title 'Exact X(t))', \\");
        command.Add("     'em2_data.txt' using 1:2 title 'EM X(t)'");
        command.Add("quit");

        File.WriteAllLines(command_filename, command);

        Console.WriteLine("  EM plot commands stored in \"" + command_filename + "\".");

    }
}