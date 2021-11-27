﻿using System;
using Burkardt.FDM;
using Burkardt.Types;

namespace FD1DHeatSteadyTest;

internal static class Program
{
    private static void Main()
    {
        problem1();
        problem2();
        problem3();
        problem4();
    }

    private static void problem1()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for problem 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 May 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int n = 11;

        Console.WriteLine("");
        Console.WriteLine("PROBLEM1:");
        Console.WriteLine("  A test problem for FD1D_HEAT_STEADY.");

        const double a = 0.0;
        const double b = 1.0;

        double[] x = typeMethods.r8vec_even(n, a, b);

        const double ua = 0.0;
        const double ub = 1.0;

        double[] u = FD1D_Heat_Steady.fd1d_heat_steady(n, a, b, ua, ub, k1, f1, x);

        const string x_file = "problem1_nodes.txt";
        typeMethods.r8mat_write(x_file, 1, n, x);

        Console.WriteLine("");
        Console.WriteLine("  X data written to \"" + x_file + "\".");

        const string u_file = "problem1_values.txt";
        typeMethods.r8mat_write(u_file, 1, n, u);

        Console.WriteLine("  U data written to \"" + u_file + "\".");
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("PROBLEM1:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static double k1(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    K1 evaluates the heat transfer coefficient K(X).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 May 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the position.
        //
        //    Output, double K1, the value of K(X).
        //
    {
        const double value = 1.0;

        return value;
    }

    private static double f1(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F1 evaluates the right hand side of the steady state heat equation.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 May 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the position.
        //
        //    Output, double F1, the value of F(X).
        //
    {
        const double value = 0.0;

        return value;
    }

    private static void problem2()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for problem 2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 May 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int n = 11;

        Console.WriteLine("");
        Console.WriteLine("PROBLEM2:");
            
        Console.WriteLine("  A test problem for FD1D_HEAT_STEADY.");
        Console.WriteLine("  Low K, then high K, then moderate K.");

        const double a = 0.0;
        const double b = 1.0;

        double[] x = typeMethods.r8vec_even(n, a, b);

        const double ua = 0.0;
        const double ub = 1.0;

        double[] u = FD1D_Heat_Steady.fd1d_heat_steady(n, a, b, ua, ub, k2, f2, x);

        const string x_file = "problem2_nodes.txt";
        typeMethods.r8mat_write(x_file, 1, n, x);

        Console.WriteLine("");
        Console.WriteLine("  X data written to \"" + x_file + "\".");

        const string u_file = "problem2_values.txt";
        typeMethods.r8mat_write(u_file, 1, n, u);

        Console.WriteLine("  U data written to \"" + u_file + "\".");
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("PROBLEM2");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static double k2(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    K2 evaluates the heat transfer coefficient K(X).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 May 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the position.
        //
        //    Output, double K2, the value of K(X).
        //
    {
        double value = x switch
        {
            < 0.5 => 0.25,
            < 0.75 => 4.0,
            _ => 1.0
        };

        return value;
    }

    private static double f2(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F2 evaluates the right hand side of the steady state heat equation.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 May 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the position.
        //
        //    Output, double F2, the value of F(X).
        //
    {
        const double value = 0.0;

        return value;
    }

    private static void problem3()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for problem 3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 May 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int n = 21;

        Console.WriteLine("");
        Console.WriteLine("PROBLEM3:");
            
        Console.WriteLine("  A test problem for FD1D_HEAT_STEADY.");
        Console.WriteLine("  Interior source term.");

        const double a = 0.0;
        const double b = 1.0;

        double[] x = typeMethods.r8vec_even(n, a, b);

        const double ua = 0.0;
        const double ub = 100.0;

        double[] u = FD1D_Heat_Steady.fd1d_heat_steady(n, a, b, ua, ub, k3, f3, x);

        const string x_file = "problem3_nodes.txt";
        typeMethods.r8mat_write(x_file, 1, n, x);

        Console.WriteLine("");
        Console.WriteLine("  X data written to \"" + x_file + "\".");

        const string u_file = "problem3_values.txt";
        typeMethods.r8mat_write(u_file, 1, n, u);

        Console.WriteLine("  U data written to \"" + u_file + "\".");
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("PROBLEM3:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static double k3(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    K3 evaluates the heat transfer coefficient K(X).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 May 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the position.
        //
        //    Output, double K3, the value of K(X).
        //
    {
        const double value = 1.0;

        return value;
    }

    private static double f3(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F1 evaluates the right hand side of the steady state heat equation.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 May 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the position.
        //
        //    Output, double F3, the value of F(X).
        //
    {
        double value = x switch
        {
            < 0.15 => 0.0,
            < 0.45 => 200.0,
            _ => 0.0
        };

        return value;
    }

    private static void problem4()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for problem 4.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 May 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int n = 21;

        Console.WriteLine("");
        Console.WriteLine("PROBLEM4:");
        Console.WriteLine("  A test problem for FD1D_HEAT_STEADY.");
        Console.WriteLine("  A heat source and a heat sink.");

        const double a = 0.0;
        const double b = 1.0;

        double[] x = typeMethods.r8vec_even(n, a, b);

        const double ua = 0.0;
        const double ub = 1.0;

        double[] u = FD1D_Heat_Steady.fd1d_heat_steady(n, a, b, ua, ub, k4, f4, x);

        const string x_file = "problem4_nodes.txt";
        typeMethods.r8mat_write(x_file, 1, n, x);

        Console.WriteLine("");
        Console.WriteLine("  X data written to \"" + x_file + "\".");

        const string u_file = "problem4_values.txt";
        typeMethods.r8mat_write(u_file, 1, n, u);

        Console.WriteLine("  U data written to \"" + u_file + "\".");
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("PROBLEM4:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static double k4(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    K4 evaluates the heat transfer coefficient K(X).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 May 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the position.
        //
        //    Output, double K4, the value of K(X).
        //
    {
        const double value = 1.0;

        return value;
    }

    private static double f4(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F4 evaluates the right hand side of the steady state heat equation.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 May 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the position.
        //
        //    Output, double F4, the value of F(X).
        //
    {
        double value = x switch
        {
            < 0.15 => 0.0,
            < 0.35 => 1.0,
            < 0.75 => 0.0,
            < 0.85 => -2.0,
            _ => 0.0
        };

        return value;
    }
}