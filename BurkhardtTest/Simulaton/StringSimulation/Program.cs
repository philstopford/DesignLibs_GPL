using System;
using System.Collections.Generic;
using System.IO;

namespace StringSimulation;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for STRING_SIMULATION.
        //
        //  Discussion:
        //
        //    This program solves the 1D wave equation of the form:
        //
        //      Utt = c^2 Uxx
        //
        //    over the spatial interval [X1,X2] and time interval [T1,T2],
        //    with initial conditions:
        //
        //      U(T1,X)  = U_T1(X),
        //      Ut(T1,X) = UT_T1(X),
        //
        //    and boundary conditions of Dirichlet type:
        //
        //      U(T,X1) = U_X1(T),
        //      U(T,X2) = U_X2(T).
        //
        //    The value C represents the propagation speed of waves.
        //
        //    The program uses the finite difference method, and marches
        //    forward in time, solving for all the values of U at the next
        //    time step by using the values known at the previous two time steps.
        //
        //    Central differences may be used to approximate both the time
        //    and space derivatives in the original differential equation.
        //
        //    Thus, assuming we have available the approximated values of U
        //    at the current and previous times, we may write a discretized
        //    version of the wave equation as follows:
        //
        //      Uxx(T,X) = ( U(T,   X+dX) - 2 U(T,X) + U(T,   X-dX) ) / dX^2
        //      Utt(T,X) = ( U(T+dt,X   ) - 2 U(T,X) + U(T-dt,X   ) ) / dT^2
        //
        //    If we multiply the first term by C^2 and solve for the single
        //    unknown value U(T+dt,X), we have:
        //
        //      U(T+dT,X) =        (     C^2 * dT^2 / dX^2 ) * U(T,   X+dX)
        //                  +  2 * ( 1 - C^2 * dT^2 / dX^2 ) * U(T,   X   )
        //                  +      (     C^2 * dT^2 / dX^2 ) * U(T,   X-dX)
        //                  -                                  U(T-dT,X   )
        //
        //    (Equation to advance from time T to time T+dT, except for FIRST step)
        //
        //    However, on the very first step, we only have the values of U
        //    for the initial time, but not for the previous time step.
        //    In that case, we use the initial condition information for dUdT
        //    which can be approximated by a central difference that involves
        //    U(T+dT,X) and U(T-dT,X):
        //
        //      dU/dT(T,X) = ( U(T+dT,X) - U(T-dT,X) ) / ( 2 * dT )
        //
        //    and so we can estimate U(T-dT,X) as
        //
        //      U(T-dT,X) = U(T+dT,X) - 2 * dT * dU/dT(T,X)
        //
        //    If we replace the "missing" value of U(T-dT,X) by the known values
        //    on the right hand side, we now have U(T+dT,X) on both sides of the
        //    equation, so we have to rearrange to get the formula we use
        //    for just the first time step:
        //
        //      U(T+dT,X) =   1/2 * (     C^2 * dT^2 / dX^2 ) * U(T,   X+dX)
        //                  +       ( 1 - C^2 * dT^2 / dX^2 ) * U(T,   X   )
        //                  + 1/2 * (     C^2 * dT^2 / dX^2 ) * U(T,   X-dX)
        //                  +  dT *                         dU/dT(T,   X   )
        //
        //    (Equation to advance from time T to time T+dT for FIRST step.)
        //
        //    It should be clear now that the quantity ALPHA = C * DT / DX will affect
        //    the stability of the calculation.  If it is greater than 1, then
        //    the middle coefficient 1-C^2 DT^2 / DX^2 is negative, and the
        //    sum of the magnitudes of the three coefficients becomes unbounded.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 December 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Local Parameters:
        //
        //    Local, double ALPHA, the CFL stability parameter.
        //
        //    Local, double C, the wave speed.
        //
        //    Local, double DT, the time step.
        //
        //    Local, double DX, the spatial step.
        //
        //    Local, int M, the number of time steps.
        //
        //    Local, int N, the number of spatial intervals.
        //
        //    Local, double T1, T2, the initial and final times.
        //
        //    Local, double U[M+1,N+1], the computed solution.
        //
        //    Local, double X1, X2, the left and right spatial endpoints.
        //
    {
        const int m = 30;
        const int n = 40;

        const double c = 0.25;
        List<string> command_unit = new();
        List<string> data_unit = new();
        int i;
        int j;
        const double t1 = 0.0;
        const double t2 = 3.0;
        double[,] u = new double[m + 1, n + 1];
        double x;
        const double x1 = 0.0;
        const double x2 = 1.0;

        Console.WriteLine("");
        Console.WriteLine("STRING_SIMULATION:");
        Console.WriteLine("  Simulate the behavior of a vibrating string.");

        const double dx = (x2 - x1) / n;
        const double dt = (t2 - t1) / m;
        double alpha = Math.Pow(c * dt / dx, 2);
        Console.WriteLine("  ALPHA = ( C * dT / dX )^2 = " + alpha + "");
        switch (alpha)
        {
            //
            //  Warn the user if ALPHA will cause an unstable computation.
            //
            case > 1.0:
                Console.WriteLine("");
                Console.WriteLine("  Warning!");
                Console.WriteLine("  ALPHA is greater than 1.");
                Console.WriteLine("  The computation is unstable.");
                break;
        }

        //
        //  Time step 0: 
        //  Use the initial condition for U.
        //
        u[0, 0] = 0.0;
        for (j = 1; j < n; j++)
        {
            x = j * dx;
            u[0, j] = f(x);
        }

        u[0, n] = 0.0;
        //
        //  Time step 1:
        //  Use the initial condition for dUdT.
        //
        u[1, 0] = 0.0;
        for (j = 1; j < n; j++)
        {
            x = j * dx;
            u[1, j] =
                alpha / 2.0 * u[0, j - 1]
                + (1.0 - alpha) * u[0, j]
                + alpha / 2.0 * u[0, j + 1]
                + dt * g(x);
        }

        u[1, n] = 0.0;
        //
        //  Time steps 2 through M:
        //
        for (i = 2; i <= m; i++)
        {
            u[i, 0] = 0.0;
            for (j = 1; j < n; j++)
            {
                u[i, j] =
                    alpha * u[i - 1, j - 1]
                    + 2.0 * (1.0 - alpha) * u[i - 1, j]
                    + alpha * u[i - 1, j + 1]
                    - u[i - 2, j];
            }

            u[i, n] = 0.0;
        }

        for (i = 0; i <= m; i++)
        {
            double t = i * dt;
            for (j = 0; j <= n; j++)
            {
                x = j * dx;
                data_unit.Add("  " + x
                                   + "  " + t
                                   + "  " + u[i, j] + "");
            }

            data_unit.Add("");
        }

        File.WriteAllLines("string_data.txt", data_unit);

        Console.WriteLine("");
        Console.WriteLine("  Plot data written to the file \"string_data.txt\".");

        command_unit.Add("set term png");
        command_unit.Add("set output \"string.png\"");
        command_unit.Add("set grid");
        command_unit.Add("set style data lines");
        command_unit.Add("unset key");
        command_unit.Add("set xlabel '<---X--->'");
        command_unit.Add("set ylabel '<---Time--->'");
        command_unit.Add("splot \"string_data.txt\" using 1:2:3 with lines");
        command_unit.Add("quit");

        File.WriteAllLines("string_commands.txt", command_unit);


        Console.WriteLine("  Gnuplot command data written to the file \"string_commands.txt\".");

        Console.WriteLine("");
        Console.WriteLine("STRING_SIMULATION:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static double f(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F supplies the initial condition.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 December 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the location.
        //
        //    Output, double F, the value of the solution at time 0 and location X.
        //
    {
        double value = x switch
        {
            >= 0.25 and <= 0.50 => (x - 0.25) * (0.50 - x),
            _ => 0.0
        };

        return value;
    }

    private static double g(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    G supplies the initial derivative.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 December 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the location.
        //
        //    Output, double G, the value of the time derivative of the solution 
        //    at time 0 and location X.
        //
    {
        const double value = 0.0;

        return value;
    }
}