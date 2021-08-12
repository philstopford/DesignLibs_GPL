using System;
using System.Collections.Generic;
using System.IO;
using Burkardt;
using Burkardt.ODE;
using Burkardt.Types;

namespace MidpointExplicitTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    midpoint_explicit_test() tests midpoint_explicit.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 April 2021
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int n;
            double[] tspan = new double[2];
            double[] y0 = new double[2];

            Console.WriteLine("");
            Console.WriteLine("midpoint_explicit_test:");
            Console.WriteLine("  Test midpoint_explicit() on several ODE's.");

            tspan[0] = 0.0;
            tspan[1] = 5.0;
            y0[0] = 5000.0;
            y0[1] = 100.0;
            n = 200;
            predator_prey_midpoint_explicit_test(tspan, y0, n);

            tspan[0] = 0.0;
            tspan[1] = 1.0;
            y0[0] = 0.0;
            n = 27;
            stiff_midpoint_explicit_test(tspan, y0, n);
            //
            //  Terminate.
            //
            Console.WriteLine("");
            Console.WriteLine("midpoint_explicit_test:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static double[] predator_prey_deriv(double t, double[] rf, int rfIndex)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    predator_prey_deriv() evaluates the right hand side of the system.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    28 April 2020
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    George Lindfield, John Penny,
            //    Numerical Methods Using MATLAB,
            //    Second Edition,
            //    Prentice Hall, 1999,
            //    ISBN: 0-13-012641-1,
            //    LC: QA297.P45.
            //
            //  Input:
            //
            //    double T: the current time.
            //
            //    double RF[2]: the current solution variables, rabbits and foxes.
            //
            //  Output:
            //
            //    double PREDATOR_PREY_DERIV[2]: the right hand side of the 2 ODE's.
            //
        {
            double[] drfdt;

            drfdt = new double[2];

            drfdt[0] = 2.0 * rf[rfIndex + 0] - 0.001 * rf[rfIndex + 0] * rf[rfIndex + 1];
            drfdt[1] = -10.0 * rf[rfIndex + 1] + 0.002 * rf[rfIndex + 0] * rf[rfIndex + 1];

            return drfdt;
        }

        static void predator_prey_midpoint_explicit_test(double[] tspan, double[] p0, int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    predator_prey_midpoint_explicit_test() tests predator_prey_midpoint_explicit().
            //
            //  Discussion:
            //
            //    The physical system under consideration is a pair of animal populations.
            //
            //    The PREY reproduce rapidly for each animal alive at the beginning of the
            //    year, two more will be born by the end of the year.  The prey do not have
            //    a natural death rate instead, they only die by being eaten by the predator.
            //    Every prey animal has 1 chance in 1000 of being eaten in a given year by
            //    a given predator.
            //
            //    The PREDATORS only die of starvation, but this happens very quickly.
            //    If unfed, a predator will tend to starve in about 1/10 of a year.
            //    On the other hand, the predator reproduction rate is dependent on
            //    eating prey, and the chances of this depend on the number of available prey.
            //
            //    The resulting differential equations can be written:
            //
            //      PREY(0) = 5000         
            //      PRED(0) =  100
            //
            //      d PREY / dT =    2 * PREY(T) - 0.001 * PREY(T) * PRED(T)
            //      d PRED / dT = - 10 * PRED(T) + 0.002 * PREY(T) * PRED(T)
            //
            //    Here, the initial values (5000,100) are a somewhat arbitrary starting point.
            //
            //    The pair of ordinary differential equations that result have an interesting
            //    behavior.  For certain choices of the interaction coefficients (such as
            //    those given here), the populations of predator and prey will tend to 
            //    a periodic oscillation.  The two populations will be out of phase the number
            //    of prey will rise, then after a delay, the predators will rise as the prey
            //    begins to fall, causing the predator population to crash again.
            //
            //    There is a conserved quantity, which here would be:
            //      E(r,f) = 0.002 r + 0.001 f - 10 ln(r) - 2 ln(f)
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 April 2021
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    George Lindfield, John Penny,
            //    Numerical Methods Using MATLAB,
            //    Second Edition,
            //    Prentice Hall, 1999,
            //    ISBN: 0-13-012641-1,
            //    LC: QA297.P45.
            //
            //  Input:
            //
            //    double TSPAN = [ T0, TMAX ], the initial and final times.
            //    A reasonable value might be [ 0, 5 ].
            //
            //    double P0 = [ PREY, PRED ], the initial number of prey and predators.
            //    A reasonable value might be [ 5000, 100 ].
            //
            //    int N: the number of time steps.
            //
        {
            string command_filename;
            List<string> command = new List<string>();
            string data_filename;
            List<string> data = new List<string>();
            string header = "predator_prey_midpoint_explicit";
            int i;
            const int m = 2;
            double[] pout;
            double[] t;

            Console.WriteLine("");
            Console.WriteLine("predator_prey_midpoint_test");
            Console.WriteLine("  A pair of ordinary differential equations for a population");
            Console.WriteLine("  of predators and prey are solved using midpoint_explicit().");

            t = new double[n + 1];
            pout = new double[(n + 1) * m];

            MidpointExplicit.midpoint_explicit(predator_prey_deriv, tspan, p0, n, m, ref t, ref pout);
            //
            //  Create the data file.
            //
            data_filename = header + "_data.txt";

            for (i = 0; i < n; i++)
            {
                data.Add("  " + t[i]
                              + "  " + pout[0 + i * m]
                              + "  " + pout[1 + i * m] + "");
            }

            File.WriteAllLines(data_filename, data);

            Console.WriteLine("");
            Console.WriteLine("  predator_prey_midpoint_explicit_test: data stored in '" + data_filename + "'.");
            //
            //  Create the command file.
            //
            command_filename = header + "_commands.txt";

            command.Add("# " + command_filename + "");
            command.Add("#");
            command.Add("# Usage:");
            command.Add("#  gnuplot < " + command_filename + "");
            command.Add("#");
            command.Add("set term png");
            command.Add("set output '" + header + ".png'");
            command.Add("set xlabel '<-- PREDATOR -->'");
            command.Add("set ylabel '<-- PREY -->'");
            command.Add("set title 'Predator prey: midpoint explicit'");
            command.Add("set grid");
            command.Add("set style data lines");
            command.Add("plot '" + data_filename + "' using 2:3 with lines lw 3");
            command.Add("quit");

            File.WriteAllLines(command_filename, command);

            Console.WriteLine("  predator_prey_midpoint_explicit_test: plot commands stored in '" + command_filename +
                              "'.");
        }

        static double[] stiff_deriv(double t, double[] y, int yindex)

            //****************************************************************************80
            //
            //  stiff_deriv() evaluates the right hand side of the stiff ODE.
            //
            //  Discussion:
            //
            //    y' = 50 * ( cos(t) - y )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 April 2021
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Input:
            //
            //    double T, Y[1]: the time and solution value.
            //
            //  Output:
            //
            //    double DYDT[1]: the derivative value.
            //
        {
            double[] dydt;

            dydt = new double[1];

            dydt[0] = 50.0 * (Math.Cos(t) - y[yindex + 0]);

            return dydt;
        }

        static double[] stiff_exact(int n, double[] t)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    stiff_exact() evaluates the exact solution of the stiff ODE.
            //
            //  Discussion:
            //
            //    y' = 50 * ( cos(t) - y )
            //    y(0) = 0
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    29 April 2020
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Input:
            //
            //    int N: the number of values.
            //
            //    double T[N]: the evaluation times.
            //
            //  Output:
            //
            //    double STIFF_EXACT[N]: the exact solution values.
            //
        {
            int i;
            double[] y;

            y = new double[n];

            for (i = 0; i < n; i++)
            {
                y[i] = 50.0 * (Math.Sin(t[i]) + 50.0 * Math.Cos(t[i])
                               - 50.0 * Math.Exp(-50.0 * t[i])) / 2501.0;
            }

            return y;
        }

        static void stiff_midpoint_explicit_test(double[] tspan, double[] y0, int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    stiff_midpoint_explicit_test() tests stiff_midpoint_explicit().
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 April 2021
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Input:
            //
            //    double TSPAN(2): the first and last times.
            //
            //    double Y0: the initial condition.
            //
            //    int N: the number of steps to take.
            //
        {
            int m = 1;
            const int n2 = 101;
            double[] t1;
            double[] t2;
            double[] y1;
            double[] y2;

            Console.WriteLine("");
            Console.WriteLine("stiff_midpoint_explicit_test");
            Console.WriteLine("  Solve stiff ODE using the midpoint_explicit method.");

            t1 = new double[n + 1];
            y1 = new double[n + 1];
            MidpointExplicit.midpoint_explicit(stiff_deriv, tspan, y0, n, m, ref t1, ref y1);

            t2 = typeMethods.r8vec_linspace_new(n2, tspan[0], tspan[1]);
            y2 = stiff_exact(n2, t2);

            plot2(n + 1, t1, y1, n2, t2, y2, "stiff_midpoint_explicit",
                "Stiff ODE: midpoint explicit method");
        }

        static void plot2(int n1, double[] t1, double[] y1, int n2, double[] t2,
                double[] y2, string header, string title)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    plot2() plots two curves together.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    29 April 2020
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Input:
            //
            //    int N1: the size of the first data set.
            //
            //    double T1(N1), Y1(N1), the first dataset.
            //
            //    int N2: the size of the second data set.
            //
            //    double T2(N2), Y2(N2), the secod dataset.
            //
            //    string HEADER: an identifier for the data.
            //
            //    string TITLE: a title to appear in the plot.
            //
        {
            string command_filename;
            List<string> command = new List<string>();
            string data1_filename;
            string data2_filename;
            List<string> data = new List<string>();
            int i;
            //
            //  Create the data files.
            //
            data1_filename = header + "_data1.txt";
            for (i = 0; i < n1; i++)
            {
                data.Add("  " + t1[i] + "  " + y1[i] + "");
            }

            File.WriteAllLines(data1_filename, data);
            data.Clear();

            data2_filename = header + "_data2.txt";
            for (i = 0; i < n2; i++)
            {
                data.Add("  " + t2[i] + "  " + y2[i] + "");
            }

            File.WriteAllLines(data2_filename, data);
            data.Clear();

            Console.WriteLine("");
            Console.WriteLine("  plot2: data stored in '" + data1_filename
                                                          + "' and '" + data2_filename + "'");
            //
            //  Create the command file.
            //
            command_filename = header + "_commands.txt";

            command.Add("# " + command_filename + "");
            command.Add("#");
            command.Add("# Usage:");
            command.Add("#  gnuplot < " + command_filename + "");
            command.Add("#");
            command.Add("set term png");
            command.Add("set output '" + header + ".png'");
            command.Add("set xlabel '<-- T -->'");
            command.Add("set ylabel '<-- Y(T) -->'");
            command.Add("set title '" + title + "'");
            command.Add("set grid");
            command.Add("set style data lines");
            command.Add("plot '" + data1_filename + "' using 1:2 with lines lw 3 lt rgb 'red',\\");
            command.Add("     '" + data2_filename + "' using 1:2 with lines lw 3 lt rgb 'blue'");
            command.Add("quit");

            File.WriteAllLines(command_filename, command);

            Console.WriteLine("  plot2: plot commands stored in '" + command_filename + "'.");
        }
    }
}