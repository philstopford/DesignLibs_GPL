using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.Simulation;
using Burkardt.Types;

namespace PoissonSimulationTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for POISSON_SIMULATION_TEST.
        //
        //  Discussion:
        //
        //    POISSON_SIMULATION_TEST tests the POISSON_SIMULATION library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 September 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("POISSON_SIMULATION_TEST");
        Console.WriteLine("  Test the POISSON_SIMULATION library.");

        test01();
        test02();

        Console.WriteLine("");
        Console.WriteLine("POISSON_SIMULATION_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 simulates waiting for a given number of events.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 September 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int bin_num = 30;
        List<string> command = new();
        List<string> data = new();
        const int event_num = 1000;
        int i;

        Console.WriteLine("");
        Console.WriteLine("TEST01:");
        Console.WriteLine("  POISSON_FIXED_EVENTS simulates a Poisson process");
        Console.WriteLine("  until a given number of events have occurred.");
        Console.WriteLine("");
        Console.WriteLine("  Simulate a Poisson process, for which, on average,");
        Console.WriteLine("  LAMBDA events occur per unit time.");
        Console.WriteLine("  Run until you have observed EVENT_NUM events.");

        double lambda = 0.5;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("  LAMBDA = " + lambda + "");
        Console.WriteLine("  EVENT_NUM = " + event_num + "");

        double[] t = new double[event_num + 1];
        double[] w = new double[event_num + 1];
        Poisson.poisson_fixed_events(lambda, event_num, ref seed, ref t, ref w);

        double w_min = typeMethods.r8vec_min(event_num + 1, w);
        double w_max = typeMethods.r8vec_max(event_num + 1, w);
        double w_ave = typeMethods.r8vec_mean(event_num + 1, w);

        Console.WriteLine("");
        Console.WriteLine("  Minimum wait = " + w_min + "");
        Console.WriteLine("  Average wait = " + w_ave + "");
        Console.WriteLine("  Maximum wait = " + w_max + "");

        Console.WriteLine("");
        Console.WriteLine(" Count            Time            Wait");
        Console.WriteLine("");
        for (i = 0; i <= 5; i++)
        {
            Console.WriteLine("  " + i
                                   + "  " + t[i]
                                   + "  " + w[i] + "");
        }

        Console.WriteLine("  ....  ..............  ..............");
        for (i = event_num - 5; i <= event_num; i++)
        {
            Console.WriteLine("  " + i
                                   + "  " + t[i]
                                   + "  " + w[i] + "");
        }

        //
        //  Create the data file.
        //
        string data_filename = "poisson_timeline_data.txt";

        for (i = 0; i <= event_num; i++)
        {
            data.Add("  " + t[i]
                          + "  " + i + "");
        }

        File.WriteAllLines(data_filename, data);
        Console.WriteLine(" ");
        Console.WriteLine("  Data stored in \"" + data_filename + "\".");
        //
        //  Create the command file.
        //
        string command_filename = "poisson_timeline_commands.txt";

        command.Add("# poisson_timeline_commands.txt");
        command.Add("#");
        command.Add("# Usage:");
        command.Add("#  gnuplot < poisson_timeline_commands.txt");
        command.Add("#");
        command.Add("set term png");
        command.Add("set output 'poisson_timeline.png'");
        command.Add("set style data lines");
        command.Add("set xlabel 'Time'");
        command.Add("set ylabel 'Number of Events'");
        command.Add("set title 'Observation of Fixed Number of Poisson Events'");
        command.Add("set grid");
        command.Add("plot 'poisson_timeline_data.txt' using 1:2 lw 2");
        command.Add("quit");

        File.WriteAllLines(command_filename, command);

        Console.WriteLine("  Plot commands stored in \"" + command_filename + "\".");
        //
        //  Determine bin information.
        //
        w_min = typeMethods.r8vec_min(event_num + 1, w);
        w_max = typeMethods.r8vec_max(event_num + 1, w);

        double[] w_bin = typeMethods.r8vec_midspace_new(bin_num, w_min, w_max);
        int[] f_bin = new int[bin_num];

        for (i = 0; i < bin_num; i++)
        {
            f_bin[i] = 0;
        }

        for (i = 0; i < event_num; i++)
        {
            int j = 1 + (int)(bin_num * (w[i] - w_min) / (w_max - w_min));
            j = Math.Min(j, bin_num-1);
            f_bin[j] += 1;
        }

        //
        //  Create the data file.
        //
        data_filename = "poisson_times_data.txt";

        data.Clear();

        for (i = 0; i < bin_num; i++)
        {
            data.Add("  " + w_bin[i]
                          + "  " + f_bin[i] + "");
        }

        File.WriteAllLines(data_filename, data);

        Console.WriteLine(" ");
        Console.WriteLine("  Data stored in \"" + data_filename + "\".");
        //
        //  Create the command file.
        //
        command_filename = "poisson_times_commands.txt";

        command.Clear();
        command.Add("# poisson_times_commands.txt");
        command.Add("#");
        command.Add("# Usage:");
        command.Add("#  gnuplot < poisson_times_commands.txt");
        command.Add("#");
        command.Add("set term png");
        command.Add("set output 'poisson_times.png'");
        command.Add("set xlabel 'Waiting Time'");
        command.Add("set ylabel 'Frequency'");
        command.Add("set title 'Waiting Times Observed Over Fixed Time'");
        command.Add("set grid");
        command.Add("set style fill solid");
        double width = 0.85 * (w_max - w_min) / bin_num;
        command.Add("plot 'poisson_times_data.txt' using 1:2:(" + width + ") with boxes");
        command.Add("quit");

        File.WriteAllLines(command_filename, command);

        Console.WriteLine("  Plot commands stored in \"" + command_filename + "\".");
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 simulates waiting for a given length of time.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 September 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int bin_num = 30;
        List<string> command = new();
        List<string> data = new();
        int i;
        int test;
        const int test_num = 20000;

        const double lambda = 0.5;
        const double t = 1000.0;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("TEST02:");
        Console.WriteLine("  POISSON_FIXED_EVENTS simulates a Poisson process");
        Console.WriteLine("  counting the number of events that occur during");
        Console.WriteLine("  a given time.");
        Console.WriteLine("");
        Console.WriteLine("  Simulate a Poisson process, for which, on average,");
        Console.WriteLine("  LAMBDA events occur per unit time.");
        Console.WriteLine("  Run for a total of " + t + " time units.");
        Console.WriteLine("  LAMBDA = " + lambda + "");

        int[] n = new int[test_num];

        for (test = 0; test < test_num; test++)
        {
            n[test] = Poisson.poisson_fixed_time(lambda, t, ref seed);
        }

        double n_mean = typeMethods.i4vec_mean(test_num, n);
        double n_var = typeMethods.i4vec_variance(test_num, n);
        Console.WriteLine("");
        Console.WriteLine("  Mean number of events = " + n_mean + "");
        Console.WriteLine("  Variance = " + n_var + "");
        Console.WriteLine("  STD = " + Math.Sqrt(n_var) + "");

        double n_min = typeMethods.i4vec_min(test_num, n);
        double n_max = typeMethods.i4vec_max(test_num, n);

        double[] n_bin = typeMethods.r8vec_midspace_new(bin_num, n_min, n_max);

        int[] f_bin = new int[bin_num];
        for (i = 0; i < bin_num; i++)
        {
            f_bin[i] = 0;
        }

        for (test = 0; test < test_num; test++)
        {
            i = 1 + (int)(bin_num * (n[test] - n_min)
                          / (n_max - n_min));
            i = Math.Min(i, bin_num-1);
            f_bin[i] += 1;
        }

        //
        //  Create the data file.
        //
        const string data_filename = "poisson_events_data.txt";

        data.Clear();

        for (i = 0; i < bin_num; i++)
        {
            data.Add("  " + n_bin[i]
                          + "  " + f_bin[i] + "");
        }

        File.WriteAllLines(data_filename, data);
        Console.WriteLine(" ");
        Console.WriteLine("  Data stored in \"" + data_filename + "\".");
        //
        //  Create the command file.
        //
        const string command_filename = "poisson_events_commands.txt";

        command.Clear();
        command.Add("# poisson_events_commands.txt");
        command.Add("#");
        command.Add("# Usage:");
        command.Add("#  gnuplot < poisson_events_commands.txt");
        command.Add("#");
        command.Add("set term png");
        command.Add("set output 'poisson_events.png'");
        command.Add("set xlabel 'Number of Events'");
        command.Add("set ylabel 'Frequency'");
        command.Add("set title 'Number of Poisson Events Over Fixed Time'");
        command.Add("set grid");
        command.Add("set style fill solid");
        double w = 0.85 * (n_max - n_min) / bin_num;
        command.Add("plot 'poisson_events_data.txt' using 1:2:(" + w + ") with boxes");
        command.Add("quit");

        File.WriteAllLines(command_filename, command);

        Console.WriteLine("  Plot commands stored in \"" + command_filename + "\".");

    }
}