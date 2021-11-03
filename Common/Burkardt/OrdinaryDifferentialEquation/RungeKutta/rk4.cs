using System;

namespace Burkardt.ODENS
{
    public static partial class RungeKutta
    {
        public static double[] rk4vec(double t0, int m, double[] u0, double dt,
                Func<double, int, double[], int, double[]> f, int index = 0)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    RK4VEC takes one Runge-Kutta step for a vector ODE.
            //
            //  Discussion:
            //
            //    It is assumed that an initial value problem, of the form
            //
            //      du/dt = f ( t, u )
            //      u(t0) = u0
            //
            //    is being solved.
            //
            //    If the user can supply current values of t, u, a stepsize dt, and a
            //    function to evaluate the derivative, this function can compute the
            //    fourth-order Runge Kutta estimate to the solution at time t+dt.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    09 October 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double T0, the current time.
            //
            //    Input, int M, the spatial dimension.
            //
            //    Input, double U0[M], the solution estimate at the current time.
            //
            //    Input, double DT, the time step.
            //
            //    Input, double[] F ( double T, int M, double U[] ), a function which evaluates
            //    the derivative, or right hand side of the problem.
            //
            //    Output, double RK4VEC[M], the fourth-order Runge-Kutta solution estimate
            //    at time T0+DT.
            //
        {
            double[] f0;
            double[] f1;
            double[] f2;
            double[] f3;
            int i;
            double t1;
            double t2;
            double t3;
            double[] u;
            double[] u1;
            double[] u2;
            double[] u3;
            //
            //  Get four sample values of the derivative.
            //
            f0 = f(t0, m, u0, index);

            t1 = t0 + dt / 2.0;
            u1 = new double[m];
            for (i = 0; i < m; i++)
            {
                u1[i] = u0[index + i] + dt * f0[i] / 2.0;
            }

            f1 = f(t1, m, u1, 0);

            t2 = t0 + dt / 2.0;
            u2 = new double[m];
            for (i = 0; i < m; i++)
            {
                u2[i] = u0[index + i] + dt * f1[i] / 2.0;
            }

            f2 = f(t2, m, u2, 0);

            t3 = t0 + dt;
            u3 = new double[m];
            for (i = 0; i < m; i++)
            {
                u3[i] = u0[index + i] + dt * f2[i];
            }

            f3 = f(t3, m, u3, 0);
            //
            //  Combine them to estimate the solution.
            //
            u = new double[m];
            for (i = 0; i < m; i++)
            {
                u[i] = u0[index + i] + dt * (f0[i] + 2.0 * f1[i] + 2.0 * f2[i] + f3[i]) / 6.0;
            }

            return u;
        }
    }
}