using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.Types;

namespace Burkardt.StochasticDifferentialEquations
{
    public static class EulerMaruyamaStrong
    {
        public static void emstrong(ref int seed, int m, int n, int p_max, ref double[] dtvals,
                ref double[] xerr)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    EMSTRONG tests the strong convergence of the EM method.
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
            //    The Euler-Maruyama method uses 5 different timesteps: 
            //      16*dt, 8*dt, 4*dt, 2*dt, dt.
            //
            //    We are interested in examining strong convergence at T=1,
            //    that is
            //
            //      E | X_L - X(T) |.
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
            //    Input/output, int &SEED, a seed for the random 
            //    number generator.
            //
            //    Input, int M, the number of simulations to perform.
            //    A typical value is M = 1000.
            //
            //    Input, int N, the number of time steps to take.
            //    A typical value is N = 512.
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
            double[] a;
            double dt;
            double dt2;
            double[] dw;
            double e;
            int i;
            int j;
            int k;
            int l;
            double lambda;
            double mu;
            int p;
            int r;
            double resid;
            double[] rhs;
            int s;
            double[] sol;
            double tmax;
            double[] w;
            double winc;
            double xtemp;
            double xtrue;
            double xzero;
            //
            //  Set problem parameters.
            //
            lambda = 2.0;
            mu = 1.0;
            xzero = 1.0;
            //
            //  Set stepping parameters.
            //
            tmax = 1.0;
            dt = tmax / (double) (n);

            for (p = 0; p < p_max; p++)
            {
                dtvals[p] = dt * Math.Pow(2.0, p);
            }

            //
            //  Sample over discrete Brownian paths.
            //
            for (p = 0; p < p_max; p++)
            {
                xerr[p] = 0.0;
            }

            for (s = 0; s < m; s++)
            {
                //
                //  Define the increments dW.
                //
                dw = typeMethods.r8vec_normal_01_new(n, ref seed);

                for (j = 0; j < n; j++)
                {
                    dw[j] = Math.Sqrt(dt) * dw[j];
                }

                //
                //  Sum the increments to get the Brownian path.
                //
                w = new double[n + 1];
                w[0] = 0.0;
                for (j = 1; j <= n; j++)
                {
                    w[j] = w[j - 1] + dw[j - 1];
                }

                //
                //  Determine the true solution.
                //
                xtrue = xzero * Math.Exp((lambda - 0.5 * mu * mu) + mu * w[n]);
                //
                //  Use the Euler-Maruyama method with 5 different time steps dt2 = r * dt
                //  to estimate the solution value at time TMAX.
                //
                for (p = 0; p < p_max; p++)
                {
                    dt2 = dtvals[p];
                    r = (int) Math.Pow(2, p);
                    l = n / r;
                    xtemp = xzero;
                    for (j = 0; j < l; j++)
                    {
                        winc = 0.0;
                        for (k = r * j; k < r * (j + 1); k++)
                        {
                            winc = winc + dw[k];
                        }

                        xtemp = xtemp + dt2 * lambda * xtemp + mu * xtemp * winc;
                    }

                    xerr[p] = xerr[p] + Math.Abs(xtemp - xtrue);
                }
            }

            for (p = 0; p < p_max; p++)
            {
                xerr[p] = xerr[p] / (double) (m);
            }

            //
            //  Least squares fit of error = c * dt^q.
            //
            a = new double[p_max * 2];
            rhs = new double[p_max];

            for (i = 0; i < p_max; i++)
            {
                a[i + 0 * p_max] = 1.0;
                a[i + 1 * p_max] = Math.Log(dtvals[i]);
                rhs[i] = Math.Log(xerr[i]);
            }

            sol = QRSolve.qr_solve(p_max, 2, a, rhs);

            Console.WriteLine("");
            Console.WriteLine("EMSTRONG:");
            Console.WriteLine("  Least squares solution to Error = c * dt ^ q");
            Console.WriteLine("  (Expecting Q to be about 1/2.)");
            Console.WriteLine("  Computed Q = " + sol[1] + "");

            resid = 0.0;
            for (i = 0; i < p_max; i++)
            {
                e = a[i + 0 * p_max] * sol[0] + a[i + 1 * p_max] * sol[1] - rhs[i];
                resid = resid + e * e;
            }

            resid = Math.Sqrt(resid);
            Console.WriteLine("  Residual is " + resid + "");

        }

        public static void emstrong_gnuplot(int p_max, double[] dtvals, double[] xerr)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    EMSTRONG_GNUPLOT writes a GNUPLOT input file to plot EMSTRONG data.
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
            string command_filename = "emstrong_commands.txt";
            List<string> command = new List<string>();
            string data_filename = "emstrong_data.txt";
            List<string> data = new List<string>();
            int i;

            for (i = 0; i < p_max; i++)
            {
                data.Add("  " + dtvals[i]
                              + "  " + xerr[i]
                              + "  " + Math.Sqrt(dtvals[i]) + "");
            }

            File.WriteAllLines(data_filename, data);
            Console.WriteLine("");
            Console.WriteLine("  EMSTRONG data stored in \"" + data_filename + "\".");

            command.Add("# emstrong_commands.txt");
            command.Add("# created by sde::emstrong_gnuplot.");
            command.Add("#");
            command.Add("# Usage:");
            command.Add("#  gnuplot < emstrong_commands.txt");
            command.Add("#");
            command.Add("set term png");
            command.Add("set output 'emstrong.png'");
            command.Add("set xlabel 'Log(dt)'");
            command.Add("set ylabel 'Log(Averaged Error at final T)'");
            command.Add("set logscale xy 10");
            command.Add("set title 'Euler-Maruyama Error as function of DT'");
            command.Add("set grid");
            command.Add("set style data linespoints");
            command.Add("plot 'emstrong_data.txt' using 1:2 title 'Error', \\");
            command.Add("     'emstrong_data.txt' using 1:3 title 'Slope = 1/2'");
            command.Add("quit");

            File.WriteAllLines(command_filename, command);

            Console.WriteLine("  EMSTRONG plot commands stored in \"" + command_filename + "\".");

        }
    }
}