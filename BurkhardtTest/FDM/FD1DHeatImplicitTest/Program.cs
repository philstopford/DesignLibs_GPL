using System;
using Burkardt;
using Burkardt.IO;
using Burkardt.Types;

namespace FD1DHeatImplicitTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for FD1D_HEAT_IMPLICIT.
        //
        //  Discussion:
        //
        //    FD1D_HEAT_IMPLICIT solves the 1D heat equation with an implicit method.
        //
        //    This program solves
        //
        //      dUdT - k * d2UdX2 = F(X,T)
        //
        //    over the interval [A,B] with boundary conditions
        //
        //      U(A,T) = UA(T),
        //      U(B,T) = UB(T),
        //
        //    over the time interval [T0,T1] with initial conditions
        //
        //      U(X,T0) = U0(X)
        //
        //    The code uses the finite difference method to approximate the
        //    second derivative in space, and an implicit backward Euler approximation
        //    to the first derivative in time.
        //
        //    The finite difference form can be written as
        //
        //      U(X,T+dt) - U(X,T)                  ( U(X-dx,T) - 2 U(X,T) + U(X+dx,T) )
        //      ------------------  = F(X,T) + k *  ------------------------------------
        //               dt                                   dx * dx
        //
        //    so that we have the following linear system for the values of U at time T+dt:
        //
        //            -     k * dt / dx / dx   * U(X-dt,T+dt)
        //      + ( 1 + 2 * k * dt / dx / dx ) * U(X,   T+dt)
        //            -     k * dt / dx / dx   * U(X+dt,T+dt)
        //      =               dt             * F(X,   T+dt)
        //      +                                U(X,   T)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 May 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] a;
        double[] b;
        double[] fvec;
        bool header;
        int i;
        int j;
        int job;
        double k;
        double[] t;
        double t_delt;
        string t_file;
        double t_max;
        double t_min;
        int t_num;
        double[] u;
        string u_file;
        double w;
        double[] x;
        double x_delt;
        string x_file;
        double x_max;
        double x_min;
        int x_num;

        Console.WriteLine("");
        Console.WriteLine("FD1D_HEAT_IMPLICIT");
        Console.WriteLine("");
        Console.WriteLine("  Finite difference solution of");
        Console.WriteLine("  the time dependent 1D heat equation");
        Console.WriteLine("");
        Console.WriteLine("    Ut - k * Uxx = F(x,t)");
        Console.WriteLine("");
        Console.WriteLine("  for space interval A <= X <= B with boundary conditions");
        Console.WriteLine("");
        Console.WriteLine("    U(A,t) = UA(t)");
        Console.WriteLine("    U(B,t) = UB(t)");
        Console.WriteLine("");
        Console.WriteLine("  and time interval T0 <= T <= T1 with initial condition");
        Console.WriteLine("");
        Console.WriteLine("    U(X,T0) = U0(X).");
        Console.WriteLine("");
        Console.WriteLine("  A second order difference approximation is used for Uxx.");
        Console.WriteLine("");
        Console.WriteLine("  A first order backward Euler difference approximation");
        Console.WriteLine("  is used for Ut.");

        k = 5.0E-07;
        //
        //  Set X values.
        //
        x_min = 0.0;
        x_max = 0.3;
        x_num = 11;
        x_delt = (x_max - x_min) / (x_num - 1);

        x = new double[x_num];

        for (i = 0; i < x_num; i++)
        {
            x[i] = ((x_num - i - 1) * x_min
                    + i * x_max)
                   / (x_num - 1);
        }

        // 
        //  Set T values.
        //
        t_min = 0.0;
        t_max = 22000.0;
        t_num = 51;
        t_delt = (t_max - t_min) / (t_num - 1);

        t = new double[t_num];

        for (j = 0; j < t_num; j++)
        {
            t[j] = ((t_num - j - 1) * t_min
                    + j * t_max)
                   / (t_num - 1);
        }

        //
        //  Set the initial data, for time T_MIN.
        //
        u = new double[x_num * t_num];

        u0(x_min, x_max, t_min, x_num, x, ref u);
        //
        //  The matrix A does not change with time.  We can set it once,
        //  factor it once, and solve repeatedly.
        //
        w = k * t_delt / x_delt / x_delt;

        a = new double[3 * x_num];

        a[0 + 0 * 3] = 0.0;

        a[1 + 0 * 3] = 1.0;
        a[0 + 1 * 3] = 0.0;

        for (i = 1; i < x_num - 1; i++)
        {
            a[2 + (i - 1) * 3] = -w;
            a[1 + i * 3] = 1.0 + 2.0 * w;
            a[0 + (i + 1) * 3] = -w;
        }

        a[2 + (x_num - 2) * 3] = 0.0;
        a[1 + (x_num - 1) * 3] = 1.0;

        a[2 + (x_num - 1) * 3] = 0.0;
        //
        //  Factor the matrix.
        //
        typeMethods.r83_np_fa(x_num, ref a);

        b = new double[x_num];
        fvec = new double[x_num];

        for (j = 1; j < t_num; j++)
        {
            //
            //  Set the right hand side B.
            //
            b[0] = ua(x_min, x_max, t_min, t[j]);

            f(x_min, x_max, t_min, t[j - 1], x_num, x, ref fvec);

            for (i = 1; i < x_num - 1; i++)
            {
                b[i] = u[i + (j - 1) * x_num] + t_delt * fvec[i];
            }

            b[x_num - 1] = ub(x_min, x_max, t_min, t[j]);

            job = 0;
            fvec = typeMethods.r83_np_sl(x_num, a, b, job);

            for (i = 0; i < x_num; i++)
            {
                u[i + j * x_num] = fvec[i];
            }
        }

        x_file = "x.txt";
        header = false;
        DTable.dtable_write(x_file, 1, x_num, x, header);

        Console.WriteLine("");
        Console.WriteLine("  X data written to \"" + x_file + "\".");

        t_file = "t.txt";
        header = false;
        DTable.dtable_write(t_file, 1, t_num, t, header);

        Console.WriteLine("  T data written to \"" + t_file + "\".");

        u_file = "u.txt";
        header = false;
        DTable.dtable_write(u_file, x_num, t_num, u, header);

        Console.WriteLine("  U data written to \"" + u_file + "\".");

        Console.WriteLine("");
        Console.WriteLine("FD1D_HEAT_IMPLICIT");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void f(double a, double b, double t0, double t, int n, double[] x,
            ref double[] value)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F returns the right hand side of the heat equation.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 May 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, the left and right endpoints.
        //
        //    Input, double T0, the initial time.
        //
        //    Input, double T, the current time.
        //
        //    Input, int N, the number of points.
        //
        //    Input, double X[N], the current spatial positions.
        //
        //    Output, double VALUE[N], the prescribed value of U(X(:),T0).
        //
    {
        int i;

        for (i = 0; i < n; i++)
        {
            value[i] = 0.0;
        }
    }

    private static void u0(double a, double b, double t0, int n, double[] x, ref double[] value )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    U0 returns the initial condition at the starting time.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 May 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, the left and right endpoints
        //
        //    Input, double T0, the initial time.
        //
        //    Input, double T, the current time.
        //
        //    Input, int N, the number of points where initial data is needed.
        //
        //    Input, double X[N], the positions where initial data is needed.
        //
        //    Output, double VALUE[N], the prescribed value of U(X,T0).
        //
    {
        int i;

        for (i = 0; i < n; i++)
        {
            value[i] = 100.0;
        }
    }

    private static double ua(double a, double b, double t0, double t)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UA returns the Dirichlet boundary condition at the left endpoint.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 May 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, the left and right endpoints
        //
        //    Input, double T0, the initial time.
        //
        //    Input, double T, the current time.
        //
        //    Output, double UA, the prescribed value of U(A,T).
        //
    {
        double value = 20.0;

        return value;
    }

    private static double ub(double a, double b, double t0, double t)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UB returns the Dirichlet boundary condition at the right endpoint.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 May 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, the left and right endpoints
        //
        //    Input, double T0, the initial time.
        //
        //    Input, double T, the current time.
        //
        //    Output, double UB, the prescribed value of U(B,T).
        //
    {
        double value = 20.0;

        return value;
    }
}