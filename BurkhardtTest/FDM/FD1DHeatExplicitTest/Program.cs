﻿using System;
using Burkardt.FDM;
using Burkardt.Types;

namespace FD1DHeatExplicitTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for FD1D_HEAT_EXPLICIT_TEST.
        //
        //  Discussion:
        //
        //    FD1D_HEAT_EXPLICIT_TEST tests the FD1D_HEAT_EXPLICIT library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    26 January 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("FD1D_HEAT_EXPLICIT_TEST:");
            
        Console.WriteLine("  Test the FD1D_HEAT_EXPLICIT library.");

        fd1d_heat_explicit_test01();

        Console.WriteLine("");
        Console.WriteLine("FD1D_HEAT_EXPLICIT_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void fd1d_heat_explicit_test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FD1D_HEAT_EXPLICIT_TEST01 does a simple test problem
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    26 January 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("FD1D_HEAT_EXPLICIT_TEST01:");
        Console.WriteLine("  Compute an approximate solution to the time-dependent");
        Console.WriteLine("  one dimensional heat equation:");
        Console.WriteLine("");
        Console.WriteLine("    dH/dt - K * d2H/dx2 = f(x,t)");
        Console.WriteLine("");
        Console.WriteLine("  Run a simple test case.");
        //
        //  Heat coefficient.
        //
        const double k = 0.002;
        //
        //  X_NUM is the number of equally spaced nodes to use between 0 and 1.
        //
        const int x_num = 21;
        const double x_min = 0.0;
        const double x_max = 1.0;
        double[] x = typeMethods.r8vec_linspace_new(x_num, x_min, x_max);
        //
        //  T_NUM is the number of equally spaced time points between 0 and 10.0.
        //
        const int t_num = 201;
        const double t_min = 0.0;
        const double t_max = 80.0;
        double dt = (t_max - t_min) / (t_num - 1);
        double[] t = typeMethods.r8vec_linspace_new(t_num, t_min, t_max);
        //
        //  Get the CFL coefficient.
        //
        double cfl = FD1D_Heat_Explicit.fd1d_heat_explicit_cfl(k, t_num, t_min, t_max, x_num, x_min, x_max);
        //
        //  Running the code produces an array H of temperatures H(t,x),
        //  and vectors x and t.
        //
        double[] h = ic_test01(x_num, x, t[0]);
        bc_test01(x_num, x, t[0], h);

        double[] hmat = new double[x_num * t_num];

        int j = 0;
        for (i = 0; i < x_num; i++)
        {
            hmat[i + j * x_num] = h[i];
        }

        for (j = 1; j < t_num; j++)
        {
            double[] h_new = FD1D_Heat_Explicit.fd1d_heat_explicit(x_num, x, t[j - 1], dt, cfl, rhs_test01, bc_test01, h);

            for (i = 0; i < x_num; i++)
            {
                hmat[i + j * x_num] = h_new[i];
                h[i] = h_new[i];
            }
        }

        //
        //  Write the data to files.
        //
        typeMethods.r8mat_write("h_test01.txt", x_num, t_num, hmat);
        typeMethods.r8vec_write("t_test01.txt", t_num, t);
        typeMethods.r8vec_write("x_test01.txt", x_num, x);
    }

    private static double[] bc_test01(int x_num, double[] x, double t, double[] h )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BC_TEST01 evaluates the boundary conditions for problem 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    26 January 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int X_NUM, the number of nodes.
        //
        //    Input, double X[X_NUM], the node coordinates.
        //
        //    Input, double T, the current time.
        //
        //    Input, double H[X_NUM], the current heat values.
        //
        //    Output, double H[X_NUM], the current heat values, after boundary
        //    conditions have been imposed.
        //
    {
        h[0] = 90.0;
        h[x_num - 1] = 70.0;

        return h;
    }

    private static double[] ic_test01(int x_num, double[] x, double t )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    IC_TEST01 evaluates the initial condition for problem 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    26 January 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int X_NUM, the number of nodes.
        //
        //    Input, double X[X_NUM], the node coordinates.
        //
        //    Input, double T, the initial time.
        //
        //    Output, double H[X_NUM], the heat values at the initial time.
        //
    {
        int j;

        double[] h = new double[x_num];

        for (j = 0; j < x_num; j++)
        {
            h[j] = 50.0;
        }

        return h;
    }

    private static double[] rhs_test01(int x_num, double[] x, double t )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RHS_TEST01 evaluates the right hand side for problem 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 January 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int X_NUM, the number of nodes.
        //
        //    Input, double X[X_NUM], the node coordinates.
        //
        //    Input, double T, the current time.
        //
        //    Output, double RHS_TEST01[X_NUM], the source term.
        //
    {
        int i;

        double[] value = new double[x_num];

        for (i = 0; i < x_num; i++)
        {
            value[i] = 0.0;
        }

        return value;
    }
}