using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.SolveNS;
using Burkardt.Types;

namespace Burkardt.StochasticDifferentialEquations;

public static class MilsteinStrong
{
    public static void milstrong(ref typeMethods.r8vecNormalData data, ref int seed, int p_max, ref double[] dtvals, ref double[] xerr)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MILSTRONG tests the strong convergence of the Milstein method.
        //
        //  Discussion:
        //
        //    This function solves the stochastic differential equation
        //
        //      dX = sigma * X * ( k - X ) dt + beta * X dW,  
        //      X(0) = Xzero,
        //
        //    where 
        //
        //       sigma = 2, 
        //       k = 1, 
        //       beta = 1,
        //       Xzero = 0.5.
        //
        //    The discretized Brownian path over [0,1] has dt = 2^(-11).
        //
        //    The Milstein method uses timesteps 128*dt, 64*dt, 32*dt, 16*dt 
        //    (also dt for reference).
        //
        //    We examine strong convergence at T=1:  
        //
        //      E | X_L - X(T) |.
        //
        //    The code is vectorized: all paths computed simultaneously.
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
        //    Original MATLAB version by Desmond Higham.
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
        //    Input/output, int &SEED, a seed for the random 
        //    number generator.
        //
        //    Input, int P_MAX, the number of time step sizes to use.
        //    A typical value is 4.
        //
        //    Output, double DTVALS[P_MAX], the time steps used.
        //
        //    Output, double XERR[P_MAX], the averaged absolute error in the
        //    solution estimate at the final time.
        //
    {
        int i;
        int j;
        int p;
        //
        //  Set problem parameters.
        //
        const double sigma = 2.0;
        const double k = 1.0;
        const double beta = 0.25;
        const double xzero = 0.5;
        //
        //  Set stepping parameters.
        //
        double tmax = 1.0;
        int n = (int) Math.Pow(2, 11);
        double dt = tmax / n;
        //
        //  Number of paths sampled.
        //
        int m = 500;
        //
        //  Define the increments dW.
        //
        double[] dw = typeMethods.r8mat_normal_01_new(m, n, ref data, ref seed);
        for (j = 0; j < n; j++)
        {
            for (i = 0; i < m; i++)
            {
                dw[i + j * m] = Math.Sqrt(dt) * dw[i + j * m];
            }
        }

        //
        //  Estimate the reference solution at time T M times.
        //
        double[] xref = new double[m];

        for (i = 0; i < m; i++)
        {
            xref[i] = xzero;
        }

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < m; i++)
            {
                xref[i] = xref[i]
                          + dt * sigma * xref[i] * (k - xref[i])
                          + beta * xref[i] * dw[i + j * m]
                          + 0.5 * beta * beta * xref[i] * (dw[i + j * m] * dw[i + j * m] - dt);
            }
        }

        //
        //  Now compute M Milstein approximations at each of 4 timesteps,
        //  and record the average errors.
        //
        for (p = 0; p < p_max; p++)
        {
            dtvals[p] = dt * 8.0 * Math.Pow(2.0, p + 1);
        }

        for (p = 0; p < p_max; p++)
        {
            xerr[p] = 0.0;
        }

        double[] xtemp = new double[m];
        for (p = 0; p < p_max; p++)
        {
            int r = 8 * (int) Math.Pow(2, p + 1);
            double dtp = dtvals[p];
            int l = n / r;
            for (i = 0; i < m; i++)
            {
                xtemp[i] = xzero;
            }

            for (j = 0; j < l; j++)
            {
                for (i = 0; i < m; i++)
                {
                    double winc = 0.0;
                    int i2;
                    for (i2 = r * j; i2 < r * (j + 1); i2++)
                    {
                        winc += dw[i + i2 * m];
                    }

                    xtemp[i] = xtemp[i]
                               + dtp * sigma * xtemp[i] * (k - xtemp[i])
                               + beta * xtemp[i] * winc
                               + 0.5 * beta * beta * xtemp[i] * (winc * winc - dtp);
                }
            }

            xerr[p] = 0.0;
            for (i = 0; i < m; i++)
            {
                xerr[p] += Math.Abs(xtemp[i] - xref[i]);
            }

            xerr[p] /= m;
        }

        //
        //  Least squares fit of error = C * dt^q
        //
        double[] a = new double[p_max * 2];
        double[] rhs = new double[p_max];
        for (p = 0; p < p_max; p++)
        {
            a[p + 0 * p_max] = 1.0;
            a[p + 1 * p_max] = Math.Log(dtvals[p]);
            rhs[p] = Math.Log(xerr[p]);
        }

        double[] sol = QRSolve.qr_solve(p_max, 2, a, rhs);

        Console.WriteLine("");
        Console.WriteLine("MILSTEIN:");
        Console.WriteLine("  Least squares solution to Error = c * dt ^ q");
        Console.WriteLine("  Expecting Q to be about 1.");
        Console.WriteLine("  Computed Q = " + sol[1] + "");

        double resid = 0.0;
        for (i = 0; i < p_max; i++)
        {
            double e = a[i + 0 * p_max] * sol[0] + a[i + 1 * p_max] * sol[1] - rhs[i];
            resid += e * e;
        }

        resid = Math.Sqrt(resid);
        Console.WriteLine("  Residual is " + resid + "");

    }

    public static void milstrong_gnuplot(int p_max, double[] dtvals, double[] xerr)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MILSTRONG_GNUPLOT writes a GNUPLOT input file to plot MILSTRONG data.
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
        //    Input, int P_MAX, the number of time step sizes to use.
        //
        //    Input, double DTVALS(P_MAX), the time steps used.
        //
        //    Input, double XERR(P_MAX), the averaged absolute error in the
        //    solution estimate at the final time.
        //
    {
        const string command_filename = "milstrong_commands.txt";
        List<string> command = new();
        const string data_filename = "milstrong_data.txt";
        List<string> data = new();
        int i;
        for (i = 0; i < p_max; i++)
        {
            data.Add("  " + dtvals[i]
                          + "  " + xerr[i] + "");
        }

        File.WriteAllLines(data_filename, data);

        Console.WriteLine("");
        Console.WriteLine("  MILSTRONG data stored in \"" + data_filename + "\".");

        command.Add("# milstrong_commands.txt");
        command.Add("# created by sde::milstrong_gnuplot.");
        command.Add("#");
        command.Add("# Usage:");
        command.Add("#  gnuplot < milstrong_commands.txt");
        command.Add("#");
        command.Add("set term png");
        command.Add("set output 'milstrong.png'");
        command.Add("set xlabel 'Log(dt)'");
        command.Add("set ylabel 'Log(Averaged Error at final T)'");
        command.Add("set logscale xy 10");
        command.Add("set title 'Milstein Error as function of DT'");
        command.Add("set grid");
        command.Add("set style data linespoints");
        command.Add("plot 'milstrong_data.txt' using 1:2 title 'Error', \\");
        command.Add("     'milstrong_data.txt' using 1:1 title 'Slope = 1'");
        command.Add("quit");

        File.WriteAllLines(command_filename, command);
        Console.WriteLine("  MILSTRONG plot commands stored in \"" + command_filename + "\".");

    }
}