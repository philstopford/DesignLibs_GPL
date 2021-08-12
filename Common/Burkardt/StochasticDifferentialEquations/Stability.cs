using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.IO;
using Burkardt.Types;

namespace Burkardt.StochasticDifferentialEquations
{
    public static class Stability
    {
        public static void stab_asymptotic(ref typeMethods.r8vecNormalData vdata, ref typeMethods.r8NormalData data, ref int seed, int n, int p_max)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    STAB_ASYMPTOTIC examines asymptotic stability.
            //
            //  Discussion:
            //
            //    The function tests the asymptotic stability
            //    of the Euler-Maruyama method applied to a stochastic differential
            //    equation (SDE).
            //
            //    The SDE is
            //
            //      dX = lambda*X dt + mu*X dW,
            //      X(0) = Xzero,
            //
            //    where 
            //
            //      lambda is a constant,
            //      mu is a constant,
            //      Xzero = 1.
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
            //    Input, int N, the number of time steps for the
            //    first solution.
            //
            //    Input, int P_MAX, the number of time step sizes.
            //
        {
            string command_filename = "stab_asymptotic_commands.txt";
            List<string> command = new List<string>();
            string data_filename;
            string data_filename0 = "stab_asymptotic0_data.txt";
            List<string> out_data = new List<string>();
            double dt;
            double[] dtvals;
            int i;
            int j;
            double lambda;
            double mu;
            int nval;
            int p;
            double t;
            double test;
            double tmax;
            double[] u;
            double winc;
            double[] xemabs;
            double xmin;
            double xtemp;
            double xzero;

            Console.WriteLine("");
            Console.WriteLine("STAB_ASYMPTOTIC:");
            Console.WriteLine("  Investigate asymptotic stability of Euler-Maruyama");
            Console.WriteLine("  solution with stepsize DT and MU.");
            Console.WriteLine("");
            Console.WriteLine("  SDE is asymptotically stable if");
            Console.WriteLine("    Real ( lambda - 1/2 mu^2 ) < 0.");
            Console.WriteLine("");
            Console.WriteLine("  EM with DT is asymptotically stable if");
            Console.WriteLine("    E log ( | 1 + lambda dt - sqrt(dt) mu n(0,1) | ) < 0.");
            Console.WriteLine("  where n(0,1) is a normal random value.");
            //
            //  Problem parameters.
            //
            lambda = 0.5;
            mu = Math.Sqrt(6.0);
            xzero = 1.0;
            //
            //  Test the SDE.
            //
            Console.WriteLine("");
            Console.WriteLine("  Lambda = " + lambda + "");
            Console.WriteLine("  Mu =     " + mu + "");
            test = lambda - 0.5 * mu * mu;
            Console.WriteLine("  SDE asymptotic stability test = " + test + "");
            //
            //  Step parameters.
            //
            tmax = 500.0;
            //
            //  For each stepsize, compute the Euler-Maruyama solution.
            //
            data_filename = data_filename0;
            dtvals = new double[p_max];

            for (p = 0; p < p_max; p++)
            {
                nval = n * (int) Math.Pow(2, p);
                dt = tmax / (double) (nval);
                dtvals[p] = dt;
                //
                //  Test the EM for this DT.
                //
                Console.WriteLine("");
                Console.WriteLine("  dt = " + dt + "");
                u = typeMethods.r8vec_normal_01_new(1000, ref vdata, ref seed);
                for (i = 0; i < 1000; i++)
                {
                    u[i] = Math.Log(Math.Abs(1.0 + lambda * dt - Math.Sqrt(dt) * mu * u[i]));
                }

                test = typeMethods.r8vec_mean(1000, u);
                Console.WriteLine("  EM asymptotic test = " + test + "");

                xtemp = xzero;
                xemabs = new double[nval + 1];
                xemabs[0] = xtemp;

                for (j = 1; j <= nval; j++)
                {
                    winc = Math.Sqrt(dt) * typeMethods.r8_normal_01(ref data, ref seed);
                    xtemp = xtemp + dt * lambda * xtemp + mu * xtemp * winc;
                    xemabs[j] = Math.Abs(xtemp);
                }

                //
                //  Write this data to a file.
                //
                Files.filename_inc(ref data_filename);

                //
                //  We have to impose a tiny lower bound on the values because we
                //  will end up plotting their logs.
                //
                xmin = Math.Exp(-200.0);
                for (i = 0; i <= nval; i++)
                {
                    t = tmax * (double) (i) / (double) (nval);
                    out_data.Add("  " + t
                                  + "  " + Math.Max(xemabs[i], xmin) + "");
                }

                File.WriteAllLines(data_filename, out_data);

                Console.WriteLine("");
                Console.WriteLine("  Data for DT = " + dt + " stored in \"" + data_filename + "\"");
            }

            command.Add("# stab_asymptotic_commands.txt");
            command.Add("# created by sde::stab_asymptotic.");
            command.Add("#");
            command.Add("# Usage:");
            command.Add("#  gnuplot < stab_asymptotic_commands.txt");
            command.Add("#");
            command.Add("set term png");
            command.Add("set output 'stab_asymptotic.png'");
            command.Add("set xlabel 't'");
            command.Add("set ylabel '|X(t)|'");
            command.Add("set title 'Absolute value of EM Solution'");
            command.Add("set grid");
            command.Add("set logscale y 10");
            command.Add("set style data lines");

            data_filename = data_filename0;

            Files.filename_inc(ref data_filename);
            command.Add("plot '" + data_filename + "' using 1:2, \\");

            for (p = 1; p < p_max - 1; p++)
            {
                Files.filename_inc(ref data_filename);
                command.Add("     '" + data_filename + "' using 1:2, \\");
            }

            Files.filename_inc(ref data_filename);
            command.Add("     '" + data_filename + "' using 1:2");

            command.Add("quit");

            File.WriteAllLines(command_filename, command);

            Console.WriteLine("  STAB_ASYMPTOTIC plot stored in \"" + command_filename + "\".");
        }

        public static void stab_meansquare(ref typeMethods.r8vecNormalData vdata, ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    STAB_MEANSQUARE examines mean-square stability.
            //
            //  Discussion:
            //
            //    The function tests the mean-square stability
            //    of the Euler-Maruyama method applied to a stochastic differential
            //    equation (SDE).
            //
            //    The SDE is
            //
            //      dX = lambda*X dt + mu*X dW,
            //      X(0) = Xzero,
            //
            //    where 
            //
            //      lambda is a constant,
            //      mu is a constant,
            //      Xzero = 1.
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
            //    In the reference, this value is set to 100.
            //
        {
            string command_filename = "stab_meansquare_commands.txt";
            List<string> command = new List<string>();
            string data_filename0 = "stab_meansquare0_data.txt";
            string data_filename;
            List<string> data = new List<string>();
            double dt;
            int i;
            int j;
            int k;
            double lambda;
            int m;
            double mu;
            int n;
            double t;
            double test;
            double tmax;
            double[] winc;
            double[] xms;
            double[] xtemp;
            double xzero;

            Console.WriteLine("");
            Console.WriteLine("STAB_MEANSQUARE:");
            Console.WriteLine("  Investigate mean square stability of Euler-Maruyama");
            Console.WriteLine("  solution with stepsize DT and MU.");
            Console.WriteLine("");
            Console.WriteLine("  SDE is mean square stable if");
            Console.WriteLine("    Real ( lambda + 1/2 |mu|^2 ) < 0.");
            Console.WriteLine("");
            Console.WriteLine("  EM with DT is mean square stable if");
            Console.WriteLine("    |1+dt^2| + dt * |mu|^2 - 1.0 < 0.");
            //
            //  Set problem parameters.
            //
            tmax = 20.0;
            m = 50000;
            xzero = 1.0;
            //
            //  Problem parameters.
            //
            lambda = -3.0;
            mu = Math.Sqrt(3.0);
            //
            //  Test the SDE.
            //
            Console.WriteLine("");
            Console.WriteLine("  Lambda = " + lambda + "");
            Console.WriteLine("  Mu =     " + mu + "");
            test = lambda + 0.5 * mu * mu;
            Console.WriteLine("  SDE mean square stability test = " + test + "");
            //
            //  XMS is the mean square estimate of M paths.
            //
            data_filename = data_filename0;

            for (k = 0; k < 3; k++)
            {
                dt = Math.Pow(2.0, -k);
                n = 20 * (int) Math.Pow(2, k);
                //
                //  Test the EM for this DT.
                //
                Console.WriteLine("");
                Console.WriteLine("  dt = " + dt + "");
                test = Math.Pow(1.0 + dt * lambda, 2) + dt * mu * mu - 1.0;
                Console.WriteLine("  EM mean square stability test = " + test + "");

                xms = new double[n + 1];
                xtemp = new double[m];
                for (i = 0; i < m; i++)
                {
                    xtemp[i] = xzero;
                }

                xms[0] = xzero;

                for (j = 0; j <= n; j++)
                {
                    winc = typeMethods.r8vec_normal_01_new(m, ref vdata, ref seed);
                    for (i = 0; i < m; i++)
                    {
                        winc[i] = Math.Sqrt(dt) * winc[i];
                    }

                    for (i = 0; i < m; i++)
                    {
                        xtemp[i] = xtemp[i]
                                   + dt * lambda * xtemp[i]
                                   + mu * xtemp[i] * winc[i];
                    }

                    xms[j] = 0.0;
                    for (i = 0; i < m; i++)
                    {
                        xms[j] = xms[j] + xtemp[i] * xtemp[i];
                    }

                    xms[j] = xms[j] / (double) (m);
                }

                //
                //  Write this data to a file.
                //
                Files.filename_inc(ref data_filename);

                for (j = 0; j <= n; j++)
                {
                    t = tmax * (double) (j) / (double) (n);
                    data.Add("  " + t
                                  + "  " + xms[j] + "");
                }

                File.WriteAllLines(data_filename, data);

                Console.WriteLine("");
                Console.WriteLine("  Data for DT = " + dt + " stored in \"" + data_filename + "\".");


                command.Add("# stab_meansquare_commands.txt");
                command.Add("# created by sde::stab_meansquare.");
                command.Add("#");
                command.Add("# Usage:");
                command.Add("#  gnuplot < stab_meansquare_commands.txt");
                command.Add("#");
                command.Add("set term png");
                command.Add("set output 'stab_meansquare.png'");
                command.Add("set xlabel 't'");
                command.Add("set ylabel 'E|X^2(t)|'");
                command.Add("set title 'Mean Square of EM Solution'");
                command.Add("set grid");
                command.Add("set logscale y 10");
                command.Add("set style data lines");

                data_filename = data_filename0;

                Files.filename_inc(ref data_filename);
                command.Add("plot '" + data_filename + "' using 1:2, \\");

                for (k = 1; k <= 1; k++)
                {
                    Files.filename_inc(ref data_filename);
                    command.Add("     '" + data_filename + "' using 1:2, \\");
                }

                Files.filename_inc(ref data_filename);
                command.Add("     '" + data_filename + "' using 1:2");

                command.Add("quit");

                File.WriteAllLines(command_filename, command);

                Console.WriteLine("  STAB_MEANSQUARE plot commands stored in \"" + command_filename + "\".");

            }
        }
    }
}