﻿using System;
using Burkardt.Types;

namespace FD1DBurgersLeapfrogTest
{
    class Program
    {
        static void Main(string[] args)
            //
            //  Purpose:
            //
            //    FD1D_BURGERS_LEAP solves the nonviscous Burgers equation using leapfrogging.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    18 August 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    None
            //
        {
            double a;
            double b;
            double dt;
            double dx;
            int i;
            int ihi;
            int ilo;
            int n;
            int step;
            int step_num;
            double t;
            double t_init;
            double t_last;
            double[] uc;
            double[] un;
            double[] uo;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("FD1D_BURGERS_LEAP:");
            Console.WriteLine("  Solve the non-viscous time-dependent Burgers equation,");
            Console.WriteLine("  using the leap-frog method.");
            Console.WriteLine("");
            Console.WriteLine("  Equation to be solved:");
            Console.WriteLine("");
            Console.WriteLine("    du/dt + u * du/dx = 0");
            Console.WriteLine("");
            Console.WriteLine("  for x in [ a, b ], for t in [t_init, t_last]");
            Console.WriteLine("");
            Console.WriteLine("  with initial conditions:");
            Console.WriteLine("");
            Console.WriteLine("    u(x,o) = u_init");
            Console.WriteLine("");
            Console.WriteLine("  and boundary conditions:");
            Console.WriteLine("");
            Console.WriteLine("    u(a,t) = u_a(t), u(b,t) = u_b(t)");
            //
            //  Set and report the problem parameters.
            //
            n = 21;
            a = -1.0;
            b = +1.0;
            dx = (b - a) / (double) (n - 1);
            step_num = 30;
            t_init = 0.0;
            t_last = 3.0;
            dt = (t_last - t_init) / (double) (step_num);

            Console.WriteLine("");
            Console.WriteLine("  " + a + " <= X <= " + b + "");
            Console.WriteLine("  Number of nodes = " + n + "");
            Console.WriteLine("  DX = " + dx + "");
            Console.WriteLine("");
            Console.WriteLine("  " + t_init + " <= T <= " + t_last + "");
            Console.WriteLine("  Number of time steps = " + step_num + "");
            Console.WriteLine("  DT = " + dt + "");

            uc = new double[n];
            un = new double[n];
            uo = new double[n];

            x = typeMethods.r8vec_even(n, a, b);

            Console.WriteLine("");
            Console.WriteLine("  X:");
            Console.WriteLine("");
            for (ilo = 0; ilo < n; ilo = ilo + 5)
            {
                ihi = Math.Min(ilo + 5, n - 1);
                string cout = "";
                for (i = ilo; i <= ihi; i++)
                {
                    cout += "  " + x[i].ToString().PadLeft(14);
                }

                Console.WriteLine(cout);
            }

            //
            //  Set the initial condition,
            //  and apply boundary conditions to first and last entries.
            //
            step = 0;
            t = t_init;
            u_init(n, x, t, ref un);
            un[0] = u_a(x[0], t);
            un[n - 1] = u_b(x[n - 1], t);

            report(step, step_num, n, x, t, un);
            //
            //  Use Euler's method to get the first step.
            //
            step = 1;
            t = ((double) (step_num - step) * t_init
                 + (double) (step) * t_last)
                / (double) (step_num);

            for (i = 0; i < n; i++)
            {
                uc[i] = un[i];
            }

            for (i = 1; i < n - 1; i++)
            {
                un[i] = uc[i] - dt * uc[i] * (uc[i + 1] - uc[i - 1]) / 2.0 / dx;
            }

            un[0] = u_a(x[0], t);
            un[n - 1] = u_b(x[n - 1], t);

            report(step, step_num, n, x, t, un);
            //
            //  Subsequent steps use the leapfrog method.
            //
            for (step = 2; step <= step_num; step++)
            {
                t = ((double) (step_num - step) * t_init
                     + (double) (step) * t_last)
                    / (double) (step_num);

                for (i = 0; i < n; i++)
                {
                    uo[i] = uc[i];
                    uc[i] = un[i];
                }

                for (i = 1; i < n - 1; i++)
                {
                    un[i] = uo[i] - dt * uc[i] * (uc[i + 1] - uc[i - 1]) / dx;
                }

                un[0] = u_a(x[0], t);
                un[n - 1] = u_b(x[n - 1], t);

                report(step, step_num, n, x, t, un);
            }

            Console.WriteLine("");
            Console.WriteLine("FD1D_BURGERS_LEAP:");
            Console.WriteLine("  Normal end of execution.");

            Console.WriteLine("");
        }

        static void report(int step, int step_num, int n, double[] x, double t, double[] u)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    REPORT prints or plots or saves the data at the current time step.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    18 August 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int STEP, the index of the current step,
            //    between 0 and STEP_NUM.
            //
            //    Input, int STEP_NUM, the number of steps to take.
            //
            //    Input, int N, the number of nodes.
            //
            //    Input, double X[N], the coordinates of the nodes.
            //
            //    Input, double T, the current time.
            //
            //    Input, double U[N], the initial values U(X,T).
            //
        {
            int i;
            int ihi;
            int ilo;

            Console.WriteLine("");
            Console.WriteLine("  STEP = " + step + "");
            Console.WriteLine("  TIME = " + t + "");
            Console.WriteLine("");
            for (ilo = 0; ilo < n; ilo = ilo + 5)
            {
                ihi = Math.Min(ilo + 4, n - 1);
                string cout = "";
                for (i = ilo; i <= ihi; i++)
                {
                    cout += "  " + u[i].ToString().PadLeft(14);
                }

                Console.WriteLine(cout);
            }
        }

        static double u_a(double x, double t)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    U_A sets the boundary condition for U at A.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    18 August 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double X, T, the position and time.
            //
            //    Output, double U_A, the prescribed value of U(X,T).
            //
        {
            double ua = +0.5;

            return ua;
        }

        static double u_b(double x, double t)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    U_B sets the boundary condition for U at B.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    18 August 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double X, T, the position and time.
            //
            //    Output, double U_B, the prescribed value of U(X,T).
            //
        {
            double ub;

            ub = -0.5;

            return ub;
        }

        static void u_init(int n, double[] x, double t, ref double[] u)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    U_INIT sets the initial condition for U.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    18 August 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of nodes.
            //
            //    Input, double X[N], the coordinates of the nodes.
            //
            //    Input, double T, the current time.
            //
            //    Output, double U[N], the initial values U(X,T).
            //
        {
            int i;
            double q;
            double r;
            double s;
            double ua;
            double ub;

            ua = u_a(x[0], t);
            ub = u_b(x[n - 1], t);

            q = 2.0 * (ua - ub) / Math.PI;
            r = (ua + ub) / 2.0;
            //
            //  S can be varied.  It is the slope of the initial condition at the midpoint.
            //
            s = 1.0;

            for (i = 0; i < n; i++)
            {
                u[i] = -q * Math.Atan(s * (2.0 * x[i] - x[0] - x[n - 1])
                                 / (x[n - 1] - x[0])) + r;
            }
        }

    }
}