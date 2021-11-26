using System;
using Burkardt.ODENS.RungeKuttaFehlberg;
using Burkardt.Types;

namespace ThreeBodySimulation;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for SIMPLE_RKF45.
        //
        //  Discussion:
        //
        //    SIMPLE_RKF45 uses RKF45 as an integrator for the simple version
        //    of the three-body problem.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {

        Console.WriteLine("");
        Console.WriteLine("SIMPLE_RKF45");
        Console.WriteLine("  Simulate the behavior of three bodies which are");
        Console.WriteLine("  constrained to lie in a plane, moving under the");
        Console.WriteLine("  influence of gravity.");
        Console.WriteLine(" ");
        Console.WriteLine("  Use RKF45 for the ODE integrator.");

        simple_rkf45_run();
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("SIMPLE_RKF45");
        Console.WriteLine("  Normal end of execution.");

        Console.WriteLine("");
    }

    private static void simple_rkf45_run()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SIMPLE_RKF45_RUN runs the simple three body ODE system.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        const int neqn = 12;
        int step;
        const int step_num = 630;
        const string t_filename = "simple_rkf45_t.txt";
        const string y_filename = "simple_rkf45_y.txt";
        RungeKuttaFehlberg.r8RKFData data = new();

        double[] ts = new double[step_num + 1];
        double[] y = new double[neqn];
        double[] yp = new double[neqn];
        double[] ys = new double[neqn * (step_num + 1)];

        Console.WriteLine("");
        Console.WriteLine("SIMPLE_RKF45_RUN");
        Console.WriteLine("  Simulate the planar three-body problem as an ODE system");
        Console.WriteLine("  using RKF45 for the ODE integration.");

        const double abserr = 1.0E-10;
        double relerr = 1.0E-10;

        int flag = 1;

        const double t_start = 0.0;
        const double t_stop = 63.0;

        double t = 0.0;

        y[0] = 1.0;
        y[1] = -1.0;
        y[2] = 0.0;
        y[3] = 0.0;
        y[4] = 1.0;
        y[5] = 3.0;
        y[6] = 0.0;
        y[7] = 0.0;
        y[8] = -2.0;
        y[9] = -1.0;
        y[10] = 0.0;
        y[11] = 0.0;

        simple_f(t, y, yp);

        for (i = 0; i < neqn; i++)
        {
            ys[i + 0 * neqn] = y[i];
        }

        ts[0] = t;

        for (step = 1; step <= step_num; step++)
        {
            t = ((step_num - step + 1) * t_start
                 + (step - 1) * t_stop)
                / step_num;

            double t_out = ((step_num - step) * t_start
                            + step * t_stop)
                           / step_num;

            flag = RungeKuttaFehlberg.r8_rkf45(ref data, simple_f, neqn, ref y, ref yp, ref t, t_out, ref relerr,
                abserr, flag);

            if (Math.Abs(flag) != 2)
            {
                Console.WriteLine("");
                Console.WriteLine("SIMPLE_RKF45_RUN - Warning");
                Console.WriteLine("  Output value of FLAG = " + flag
                                                              + " at output time T_OUT = " + t_out + "");
            }

            flag = flag switch
            {
                7 => 2,
                _ => flag
            };

            for (i = 0; i < neqn; i++)
            {
                ys[i + step * neqn] = y[i];
            }

            ts[step] = t_out;
        }

        typeMethods.r8mat_write(t_filename, 1, step_num + 1, ts);
        typeMethods.r8mat_write(y_filename, neqn, step_num + 1, ys);

        Console.WriteLine("");
        Console.WriteLine("SIMPLE_RKF45_RUN:");
        Console.WriteLine("  Time data written to \"" + t_filename + "\".");
        Console.WriteLine("  Solution data written to \"" + y_filename + "\".");

    }

    private static double[] simple_f(double t, double[] y, double[] yp)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SIMPLE_F returns the right hand side of the three body ODE system.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double T, the value of the independent variable.
        //
        //    Input, double Y[NEQN], the value of the dependent variable.
        //
        //    Output, double YP[NEQN], the value of the derivatives.
        //
    {
        const double m0 = 5.0;
        const double m1 = 3.0;
        const double m2 = 4.0;

        double x0 = y[0];
        double y0 = y[1];

        double x1 = y[4];
        double y1 = y[5];

        double x2 = y[8];
        double y2 = y[9];

        double n0 = Math.Sqrt(Math.Pow(Math.Pow(x2 - x1, 2) + Math.Pow(y2 - y1, 2), 3));
        double n1 = Math.Sqrt(Math.Pow(Math.Pow(x0 - x2, 2) + Math.Pow(y0 - y2, 2), 3));
        double n2 = Math.Sqrt(Math.Pow(Math.Pow(x1 - x0, 2) + Math.Pow(y1 - y0, 2), 3));

        yp[0] = y[2];
        yp[1] = y[3];
        yp[2] = -m1 * (x0 - x1) / n2 - m2 * (x0 - x2) / n1;
        yp[3] = -m1 * (y0 - y1) / n2 - m2 * (y0 - y2) / n1;
        yp[4] = y[6];
        yp[5] = y[7];
        yp[6] = -m2 * (x1 - x0) / n0 - m0 * (x1 - x2) / n2;
        yp[7] = -m2 * (y1 - y0) / n0 - m0 * (y1 - y2) / n2;
        yp[8] = y[10];
        yp[9] = y[11];
        yp[10] = -m0 * (x2 - x0) / n1 - m1 * (x2 - x1) / n0;
        yp[11] = -m0 * (y2 - y0) / n1 - m1 * (y2 - y1) / n0;

        return yp;
    }
}