using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.Types;

namespace Burkardt.StochasticDifferentialEquations
{
    public static class BrownianPath
    {
        public static double[] bpath(ref int seed, int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BPATH performs a Brownian path simulation.
            //
            //  Discussion:
            //
            //    This routine computes one simulation of discretized Brownian 
            //    motion over the time interval [0,1] using N time steps.
            //    The user specifies a random number seed.  Different values of
            //    the seed will result in different realizations of the path.
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
            //    Input/output, int &SEED, a seed for the random number 
            //    generator.
            //
            //    Input, int N, the number of steps.
            //
            //    Output, double BPATH[N+1], the Brownian path.
            //
        {
            double dt;
            double[] dw;
            int j;
            double tmax;
            double[] w;

            tmax = 1.0;
            dt = tmax / (double) (n);
            //
            //  Define the increments dW.
            //
            dw = typeMethods.r8vec_normal_01_new(n, ref seed);

            for (j = 0; j < n; j++)
            {
                dw[j] = Math.Sqrt(dt) * dw[j];
            }

            //
            //  W is the sum of the previous increments.
            //
            w = new double[n + 1];

            w[0] = 0.0;
            for (j = 1; j <= n; j++)
            {
                w[j] = w[j - 1] + dw[j];
            }

            return w;
        }

        public static void bpath_gnuplot(int n, double[] w)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BPATH_GNUPLOT writes a GNUPLOT input file to plot BPATH data.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 September 2012
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
            //    Input, double W[N+1], the Brownian path.
            //
        {
            string command_filename = "bpath_commands.txt";
            List<string> command = new List<string>();
            string data_filename = "bpath_data.txt";
            List<string> data = new List<string>();
            int i;
            double t;

            for (i = 0; i <= n; i++)
            {
                t = (double) (i) / (double) (n);
                data.Add("  " + t
                              + "  " + w[i] + "");
            }

            File.WriteAllLines(data_filename, data);

            Console.WriteLine("");
            Console.WriteLine("  BPATH data stored in \"" + data_filename + "\".");

            command.Add("# bpath_commands.txt");
            command.Add("# created by sde::bpath_gnuplot.");
            command.Add("#");
            command.Add("# Usage:");
            command.Add("#  gnuplot < bpath_commands.txt");
            command.Add("#");
            command.Add("set term png");
            command.Add("set output 'bpath.png'");
            command.Add("set xlabel 't'");
            command.Add("set ylabel 'W(t)'");
            command.Add("set title 'Brownian motion by BPATH'");
            command.Add("set grid");
            command.Add("set style data lines");
            command.Add("plot 'bpath_data.txt' using 1:2");
            command.Add("quit");

            File.WriteAllLines(command_filename, command);

            Console.WriteLine("  BPATH plot commands stored in \"" + command_filename + "\".");
        }

        public static void bpath_average(ref int seed, int m, int n, ref double[] u, ref double[] umean,
                ref double error)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BPATH_AVERAGE: displays the average of 1000 Brownian paths.
            //
            //  Discussion:
            //
            //    This routine computes M simulations of discretized Brownian 
            //    motion W(t) over the time interval [0,1] using N time steps.
            //    The user specifies a random number seed.  Different values of
            //    the seed will result in a different set of realizations of the path.
            //
            //    Actually, we are interested in a function u(W(t)):
            //
            //      u(W(t)) = exp ( t + W(t)/2 )
            //
            //    The routine plots 5 of the simulations, as well as the average
            //    of all the simulations.  
            //
            //    The plot of the average should be quite smooth.  Its expected
            //    value is exp ( 9 * t / 8 ), and we compute the 'error', that is,
            //    the difference between the averaged value and this expected
            //    value.  This 'error' should decrease as the number of simulation
            //    is increased.
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
            //    Input, int M, the number of simulations to compute 
            //    and average.  A typical value is 1000.
            //
            //    Input, int N, the number of steps.  A typical value
            //    is 500.
            //
            //    Output, double U[M*(N+1)], the M paths.
            //
            //    Output, double UMEAN[N+1], the averaged path.
            //
            //    Output, double &ERROR, the maximum difference between the
            //    averaged path and the exact expected value.
            //
        {
            double dt;
            double[] dw;
            int i;
            int j;
            double[] t;
            double tmax;
            double[] w;

            tmax = 1.0;
            dt = tmax / (double) (n);

            t = new double[n + 1];
            for (j = 0; j <= n; j++)
            {
                t[j] = (double) (j) * tmax / (double) (n);
            }

            w = new double[n + 1];

            for (i = 0; i < m; i++)
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
                //  W is the sum of the previous increments.
                //
                w[0] = 0.0;
                for (j = 1; j <= n; j++)
                {
                    w[j] = w[j - 1] + dw[j - 1];
                }

                for (j = 0; j <= n; j++)
                {
                    u[i + j * m] = Math.Exp(t[j] + 0.5 * w[j]);
                }
            }

            //
            //  Average the M estimates of the path.
            //
            for (j = 0; j <= n; j++)
            {
                umean[j] = 0.0;
                for (i = 0; i < m; i++)
                {
                    umean[j] = umean[j] + u[i + j * m];
                }

                umean[j] = umean[j] / (double) (m);
            }

            error = 0.0;
            for (j = 0; j <= n; j++)
            {
                error = Math.Max(error, Math.Abs(umean[j] - Math.Exp(9.0 * t[j] / 8.0)));
            }
        }

        public static void bpath_average_gnuplot(int m, int n, double[] u, double[] umean)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BPATH_AVERAGE_GNUPLOT writes a GNUPLOT input file to plot BPATH_AVERAGE data.
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
            //    Input, int M, the number of simulations.
            //
            //    Input, int N, the number of steps. 
            //
            //    Input, double U[M*(N+1)], the M paths.
            //
            //    Input, double UMEAN[N+1], the averaged path.
            //
        {
            string command_filename = "bpath_average_commands.txt";
            List<string> command = new List<string>();
            string data_filename = "bpath_average_data.txt";
            List<string> data = new List<string>();
            int i;
            int j;
            double t;

            for (i = 0; i <= n; i++)
            {
                t = (double) (i) / (double) (n);
                string tmp = "  " + t;
                for (j = 0; j < 5; j++)
                {
                    tmp += "  " + u[j + i * m];
                }

                data.Add(tmp + "  " + umean[i] + "");
            }

            File.WriteAllLines(data_filename, data);

            Console.WriteLine("");
            Console.WriteLine("  BPATH_AVERAGE data stored in \"" + data_filename + "\".");

            command.Add("# bpath_average_commands.txt");
            command.Add("# created by sde::bpath_average_gnuplot.");
            command.Add("#");
            command.Add("# Usage:");
            command.Add("#  gnuplot < bpath_average_commands.txt");
            command.Add("#");
            command.Add("set term png");
            command.Add("set output 'bpath_average.png'");
            command.Add("set xlabel 't'");
            command.Add("set ylabel 'W(t)'");
            command.Add("set title 'Averaged Brownian paths'");
            command.Add("set grid");
            command.Add("set style data lines");
            command.Add("plot 'bpath_average_data.txt' using 1:2 title 'sample 1', \\");
            command.Add("     'bpath_average_data.txt' using 1:3 title 'sample 2', \\");
            command.Add("     'bpath_average_data.txt' using 1:4 title 'sample 3', \\");
            command.Add("     'bpath_average_data.txt' using 1:5 title 'sample 4', \\");
            command.Add("     'bpath_average_data.txt' using 1:6 title 'sample 5', \\");
            command.Add("     'bpath_average_data.txt' using 1:7 title 'average' lw 3");
            command.Add("quit");

            File.WriteAllLines(command_filename, command);

            Console.WriteLine("  BPATH_AVERAGE plot commands stored in \"" + command_filename + "\".");

        }
    }
}