using System;
using Burkardt.ODENS;

namespace StochasticRungeKuttaTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for STOCHASTIC_RK_TEST.
            //
            //  Discussion:
            //
            //    STOCHASTIC_RK_TEST tests the STOCHASTIC_RK library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 July 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("STOCHASTIC_RK_TEST");
            Console.WriteLine("  Test the STOCHASTIC_RK library.");

            test01();

            Console.WriteLine("");
            Console.WriteLine("STOCHASTIC_RK_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void test01()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST01 tests RK1_TI_STEP.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 July 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double h;
            int i;
            int n;
            double q;
            int seed;
            double t;
            double t0 = 0.0;
            double tn = 1.0;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST01");
            Console.WriteLine("  RK1_TI_STEP uses a first order RK method");
            Console.WriteLine("  for a problem whose right hand side does not");
            Console.WriteLine("  depend explicitly on time.");

            n = 10;
            x = new double[n + 1];
            h = (tn - t0) / (double) (n);
            q = 1.0;
            seed = 123456789;

            i = 0;
            t = t0;
            x[i] = 0.0;

            Console.WriteLine("");
            Console.WriteLine("         I           T             X");
            Console.WriteLine("");
            Console.WriteLine("  " + i.ToString().PadLeft(8)
                                   + "  " + t.ToString().PadLeft(14)
                                   + "  " + x[i].ToString().PadLeft(14) + "");

            for (i = 1; i <= n; i++)
            {
                t = ((double) (n - i) * t0
                     + (double) (i) * tn)
                    / (double) (n);

                x[i] = RungeKutta.rk1_ti_step(x[i - 1], t, h, q, fi, gi, ref seed);

                Console.WriteLine("  " + i.ToString().PadLeft(8)
                                       + "  " + t.ToString().PadLeft(14)
                                       + "  " + x[i].ToString().PadLeft(14) + "");
            }
        }

        static double fi(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FI is a time invariant deterministic right hand side.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 July 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double X, the argument.
            //
            //    Output, double FI, the value.
            //
        {
            double value;

            value = 1.0;

            return value;
        }

        static double gi(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GI is a time invariant stochastic right hand side.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 July 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double X, the argument.
            //
            //    Output, double GI, the value.
            //
        {
            double value;

            value = 1.0;

            return value;
        }
    }
}