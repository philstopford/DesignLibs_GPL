using System;
using System.Collections.Generic;
using System.IO;

namespace Burkardt.ODENS
{
    public static class PredatorPrey
    {
        public static double[] predator_prey_conserved(int n, double[] rf)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    predator_prey_conserved evaluates a conserved quantity.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    29 October 2020
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
            //    int N: the number of sets of variables.
            //
            //    double RF[N*2]: the current solution variables, rabbits and foxes.
            //
            //  Output:
            //
            //    double PREDATOR_PREY_CONSERVED[N]: the value of the conserved quantity.
            //
        {
            double alpha = 0;
            double beta = 0;
            double delta = 0;
            double gamma = 0;
            double[] h;
            int i;

            predator_prey_parameters(ref alpha, ref beta, ref gamma, ref delta);

            h = new double[n];

            for (i = 0; i < n; i++)
            {
                h[i] = delta * rf[0 + i * 2] - gamma * Math.Log(rf[0 + i * 2])
                    + beta * rf[1 + i * 2] - alpha * Math.Log(rf[1 + i * 2]);
            }

            return h;
        }

        public static double[] predator_prey_deriv(double t, double[] rf, int rfIndex)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    predator_prey_deriv evaluates the right hand side of the system.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    29 October 2020
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
            double alpha = 0;
            double beta = 0;
            double delta = 0;
            double[] drfdt;
            double gamma = 0;

            drfdt = new double[2];

            predator_prey_parameters(ref alpha, ref beta, ref gamma, ref delta);

            drfdt[0] = alpha * rf[rfIndex + 0] - beta * rf[rfIndex + 0] * rf[rfIndex + 1];
            drfdt[1] = -gamma * rf[rfIndex + 1] + delta * rf[rfIndex + 0] * rf[rfIndex + 1];

            return drfdt;
        }

        public static void predator_prey_parameters(ref double alpha, ref double beta, ref double gamma,
                ref double delta)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    predator_prey_parameters returns the problem parameters.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    29 October 2020
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Output:
            //
            //    double &ALPHA, &BETA, &GAMMA, &DELTA, the coefficient values.
            //
        {
            alpha = 2.0;
            beta = 0.001;
            gamma = 10.0;
            delta = 0.002;
        }

        public static void predator_prey_euler(double[] tspan, double[] p0, int n )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    predator_prey_euler solves the predator-prey system using euler().
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
        //    double TSPAN[2]: the initial and final times.
        //
        //    double P0[2]: the initial number of prey and predators.
        //
        //    int N: the number of time steps.
        //
        {
            string command_filename;
            List<string> command = new List<string>();
            string data_filename;
            List<string> data = new List<string>();
            string header = "predator_prey_euler";
            int i;
            const int m = 2;
            double[] pout;
            double[] t;

            Console.WriteLine("");
            Console.WriteLine("predator_prey_euler");
            Console.WriteLine("  A pair of ordinary differential equations for a population");
            Console.WriteLine("  of predators and prey are solved using euler().");

            t = new double[n + 1];
            pout = new double[(n + 1) * m];

            Euler.euler(predator_prey_deriv, tspan, p0, n, m, t, pout);
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
            Console.WriteLine("  predator_prey_euler: data stored in '" + data_filename + "'.");
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
            command.Add("set title 'Predator prey: euler'");
            command.Add("set grid");
            command.Add("set style data lines");
            command.Add("plot '" + data_filename + "' using 2:3 with lines lw 3");
            command.Add("quit");

            File.WriteAllLines(command_filename, command);

            Console.WriteLine("  predator_prey_euler: plot commands stored in '" + command_filename + "'.");
        }

        public static void predator_prey_midpoint(double[] tspan, double[] p0, int n )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    predator_prey_midpoint solves the predator-prey system using midpoint_fixed().
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
            string header = "predator_prey_midpoint";
            int i;
            const int m = 2;
            double[] pout;
            double[] t;

            Console.WriteLine("");
            Console.WriteLine("predator_prey_midpoint");
            Console.WriteLine("  A pair of ordinary differential equations for a population");
            Console.WriteLine("  of predators and prey are solved using midpoint_fixed().");

            t = new double[n + 1];
            pout = new double[(n + 1) * m];

            MidpointFixed.midpoint_fixed(predator_prey_deriv, tspan, p0, n, m, 0, ref t, ref pout);
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
            Console.WriteLine("  predator_prey_midpoint: data stored in '" + data_filename + "'.");
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
            command.Add("set title 'midpoint: predator prey'");
            command.Add("set grid");
            command.Add("set style data lines");
            command.Add("plot '" + data_filename + "' using 2:3 with lines lw 3");
            command.Add("quit");

            File.WriteAllLines(command_filename, command);

            Console.WriteLine("  predator_prey_midpoint: plot commands stored in '" + command_filename + "'.");
        }
    }
}