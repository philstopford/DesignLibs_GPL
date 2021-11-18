﻿using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.SolveNS;
using Burkardt.Types;

namespace Burkardt.StochasticDifferentialEquations;

public static class EulerMaruyamaWeak
{
    public static void emweak(ref typeMethods.r8vecNormalData data, ref int seed, int method, int m, int p_max, ref double[] dtvals,
            ref double[] xerr)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EMWEAK tests the weak convergence of the Euler-Maruyama method.
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
        //    The discretized Brownian path over [0,1] has dt = 2^(-9).
        //
        //    The Euler-Maruyama method will use 5 different timesteps:
        //
        //      2^(p-10),  p = 1,2,3,4,5.
        //
        //    We examine weak convergence at T=1:
        //
        //      | E (X_L) - E (X(T)) |.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 September 2012
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
        //    Volume 43, Number 3, September 2001, pages 525-546
        //
        //  Parameters:
        //
        //    Input, int &SEED, a seed for the random number generator.
        //
        //    Input, int METHOD.
        //    0, use the standard Euler-Maruyama method;
        //    1, use the weak Euler-Maruyama method.
        //
        //    Input, int M, the number of simulations to perform.
        //    A typical value is M = 1000.
        //
        //    Input, int P_MAX, the number of time step sizes to use.
        //    A typical value is 5.
        //
        //    Output, double DTVALS[P_MAX], the time steps used.
        //
        //    Output, double XERR[P_MAX], the averaged absolute error in the
        //    solution estimate at the final time.
        //
    {
        int i;
        int p;
        //
        //  Problem parameters;
        //
        const double lambda = 2.0;
        const double mu = 0.1;
        const double xzero = 1.0;
        //
        //  Stepping parameters.
        //
        for (p = 0; p < p_max; p++)
        {
            dtvals[p] = Math.Pow(2.0, p - 9);
        }

        //
        //  Take various Euler timesteps.
        //  For stepsize dt, we will need to take L Euler steps to reach time TMAX.
        //
        double[] xtemp = new double[m];
        double[] xem = new double[p_max];

        for (p = 0; p < p_max; p++)
        {
            int l = (int) Math.Pow(2, 9 - p);
            double dt = dtvals[p];

            for (i = 0; i < m; i++)
            {
                xtemp[i] = xzero;
            }

            int j;
            for (j = 0; j < l; j++)
            {
                double[] winc = typeMethods.r8vec_normal_01_new(m, ref data, ref seed);
                switch (method)
                {
                    case 0:
                    {
                        for (i = 0; i < m; i++)
                        {
                            winc[i] = Math.Sqrt(dt) * winc[i];
                        }

                        break;
                    }
                    default:
                    {
                        for (i = 0; i < m; i++)
                        {
                            winc[i] = Math.Sqrt(dt) * typeMethods.r8_sign(winc[i]);
                        }

                        break;
                    }
                }

                for (i = 0; i < m; i++)
                {
                    xtemp[i] = xtemp[i] + dt * lambda * xtemp[i]
                                        + mu * xtemp[i] * winc[i];
                }
            }

            //
            //  Average the M results for this stepsize.
            //
            xem[p] = typeMethods.r8vec_mean(m, xtemp);
        }

        //
        //  Compute the error in the estimates for each stepsize.
        //
        for (p = 0; p < p_max; p++)
        {
            xerr[p] = Math.Abs(xem[p] - Math.Exp(lambda));
        }

        //
        //  Least squares fit of error = c * dt^q.
        //
        double[] a = new double[p_max * 2];
        double[] rhs = new double[p_max];

        for (i = 0; i < p_max; i++)
        {
            a[i + 0 * p_max] = 1.0;
            a[i + 1 * p_max] = Math.Log(dtvals[i]);
            rhs[i] = Math.Log(xerr[i]);
        }

        double[] sol = QRSolve.qr_solve(p_max, 2, a, rhs);

        Console.WriteLine("");
        Console.WriteLine("EMWEAK:");
        switch (method)
        {
            case 0:
                Console.WriteLine("  Using standard Euler-Maruyama method.");
                break;
            default:
                Console.WriteLine("  Using weak Euler-Maruyama method.");
                break;
        }

        Console.WriteLine("  Least squares solution to Error = c * dt ^ q");
        Console.WriteLine("  (Expecting Q to be about 1.)");
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

    public static void emweak_gnuplot(int p_max, double[] dtvals, double[] xerr, int method)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EMWEAK_GNUPLOT writes a GNUPLOT input file to plot EMWEAK data.
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
        //    Input, double DTVALS[P_MAX], the time steps used.
        //
        //    Input, double XERR[P_MAX], the averaged absolute error in the
        //    solution estimate at the final time.
        //
        //    Input, int METHOD.
        //    0, use the standard Euler-Maruyama method;
        //    1, use the weak Euler-Maruyama method.
        //
    {
        List<string> command = new();
        List<string> data = new();
        int i;
        string data_filename = method switch
        {
            //
            //  Create data file.
            //
            0 => "emweak0_data.txt",
            _ => "emweak1_data.txt"
        };

        for (i = 0; i < p_max; i++)
        {
            data.Add("  " + dtvals[i]
                          + "  " + xerr[i] + "");
        }

        File.WriteAllLines(data_filename, data);
        Console.WriteLine("");
        Console.WriteLine("  EMWEAK data stored in \"" + data_filename + "\".");
        string command_filename = method switch
        {
            //
            //  Create the command file.
            //
            0 => "emweak0_commands.txt",
            _ => "emweak1_commands.txt"
        };

        command.Add("# " + command_filename + "");
        command.Add("# created by sde::emweak_gnuplot.");
        command.Add("#");
        command.Add("# Usage:");
        command.Add("#  gnuplot < " + command_filename + "");
        command.Add("#");
        command.Add("set term png");
        switch (method)
        {
            case 0:
                command.Add("set output 'emweak0.png'");
                break;
            default:
                command.Add("set output 'emweak1.png'");
                break;
        }

        command.Add("set xlabel 'Log(dt)'");
        command.Add("set ylabel 'Log(Averaged Error at final T)'");
        command.Add("set logscale xy 10");
        switch (method)
        {
            case 0:
                command.Add("set title 'Standard Euler-Maruyama Error as function of DT'");
                break;
            default:
                command.Add("set title 'Weak Euler-Maruyama Error as function of DT'");
                break;
        }

        command.Add("set grid");
        command.Add("set style data linespoints");
        switch (method)
        {
            case 0:
                command.Add("plot 'emweak0_data.txt' using 1:2 title 'Error', \\");
                command.Add("     'emweak0_data.txt' using 1:1 title 'Slope = 1'");
                break;
            default:
                command.Add("plot 'emweak1_data.txt' using 1:2 title 'Error', \\");
                command.Add("     'emweak1_data.txt' using 1:1 title 'Slope = 1'");
                break;
        }

        command.Add("quit");

        File.WriteAllLines(command_filename, command);
        Console.WriteLine("  EMWEAK plot commands stored in \"" + command_filename + "\".");

    }
}