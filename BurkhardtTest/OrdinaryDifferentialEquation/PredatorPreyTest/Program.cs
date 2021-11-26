using System;
using Burkardt.ODENS;

namespace PredatorPreyTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    predator_prey_ode_test tests predator_prey_ode.
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
    {
        double alpha = 0;
        double beta = 0;
        double delta = 0;
        double gamma = 0;
        double[] p0 = new double[2];
        double[] tspan = new double[2];

        Console.WriteLine("");
        Console.WriteLine("predator_prey_ode_test:");
        Console.WriteLine("  Test predator_prey_ode using euler, midpoint.");

        PredatorPrey.predator_prey_parameters ( ref alpha, ref beta, ref gamma, ref delta );

        Console.WriteLine("");
        Console.WriteLine("  Parameters:");
        Console.WriteLine("    alpha = " + alpha + "");
        Console.WriteLine("    beta  = " + beta + "");
        Console.WriteLine("    gamma = " + gamma + "");
        Console.WriteLine("    delta = " + delta + "");

        tspan[0] = 0.0;
        tspan[1] = 5.0;
        p0[0] = 5000.0;
        p0[1] = 100.0;
        int n = 200;

        PredatorPrey.predator_prey_euler ( tspan, p0, n );
        PredatorPrey.predator_prey_midpoint ( tspan, p0, n );

        Console.WriteLine("");
        Console.WriteLine("predator_prey_ode_test:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }    }