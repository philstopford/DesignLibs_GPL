using System;
using Burkardt.TOMSNS;

namespace TOMS178Test
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for TOMS178_TEST.
            //
            //  Discussion:
            //
            //    TOMS178_TEST tests the TOMS178 library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    12 February 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("TOMS178_TEST:");
            Console.WriteLine("  Test the TOMS178 library.");

            test01();
            test02();

            Console.WriteLine("");
            Console.WriteLine("TOMS178_TEST:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void test01()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST01 tests HOOKE with the Rosenbrock function.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    12 February 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] endpt;
            double eps;
            int i;
            int it;
            int itermax;
            int nvars = 2;
            double rho;
            double[] startpt;
            double value;

            endpt = new double[nvars];
            startpt = new double[nvars];

            Console.WriteLine("");
            Console.WriteLine("TEST01");
            Console.WriteLine("  HOOKE seeks a minimizer of F(X).");
            Console.WriteLine("  Here we use the Rosenbrock function.");
            //
            //  Starting guess for Rosenbrock.
            //
            startpt[0] = -1.2;
            startpt[1] = 1.0;

            Console.WriteLine("");
            Console.WriteLine("  Initial estimate X =");
            Console.WriteLine("");
            for (i = 0; i < nvars; i++)
            {
                Console.WriteLine("  " + (i + 1).ToString().PadLeft(8)
                                       + "  " + startpt[i].ToString().PadLeft(14) + "");
            }

            value = rosenbrock(startpt, nvars);

            Console.WriteLine("");
            Console.WriteLine("  F(X) = " + value + "");
            //
            //  Call HOOKE.
            //
            itermax = 5000;
            rho = 0.5;
            eps = 1.0E-06;

            it = TOMS.hooke(nvars, startpt, ref endpt, rho, eps, itermax, rosenbrock);
            //
            //  Results.
            //
            Console.WriteLine("");
            Console.WriteLine("  Number of iterations taken = " + it + "");
            Console.WriteLine("");
            Console.WriteLine("  X* = ");
            Console.WriteLine("");
            for (i = 0; i < nvars; i++)
            {
                Console.WriteLine("  " + (i + 1).ToString().PadLeft(8)
                                       + "  " + endpt[i].ToString().PadLeft(14) + "");
            }

            value = rosenbrock(endpt, nvars);

            Console.WriteLine("");
            Console.WriteLine("  F(X*) = " + value + "");
        }

        static void test02()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST02 tests HOOKE with the WOODS function.
            //
            //  Discussion:
            //
            //    The Hooke and Jeeves algorithm works well when RHO = 0.5, but
            //    does poorly when RHO = 0.6, and better when RHO = 0.8
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    12 February 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] endpt;
            double eps;
            int i;
            int it;
            int itermax;
            int nvars = 4;
            double rho;
            double[] startpt;
            double value;

            endpt = new double[nvars];
            startpt = new double[nvars];

            Console.WriteLine("");
            Console.WriteLine("TEST02");
            Console.WriteLine("  HOOKE seeks a minimizer of F(X).");
            Console.WriteLine("  Here we use the Rosenbrock function.");
            Console.WriteLine("  Here we use the Woods function.");
            //
            //  Starting guess.
            //
            startpt[0] = -3.0;
            startpt[1] = -1.0;
            startpt[2] = -3.0;
            startpt[3] = -1.0;

            Console.WriteLine("");
            Console.WriteLine("  Initial estimate X =");
            Console.WriteLine("");
            for (i = 0; i < nvars; i++)
            {
                Console.WriteLine("  " + (i + 1).ToString().PadLeft(8)
                                       + "  " + startpt[i].ToString().PadLeft(14) + "");
            }

            value = woods(startpt, nvars);

            Console.WriteLine("");
            Console.WriteLine("  F(X) = " + value + "");
            //
            //  Call HOOKE.
            //
            itermax = 5000;
            rho = 0.5;
            eps = 1.0E-06;

            it = TOMS.hooke(nvars, startpt, ref endpt, rho, eps, itermax, woods);
            //
            //  Results.
            //
            Console.WriteLine("");
            Console.WriteLine("  Number of iterations taken = " + it + "");
            Console.WriteLine("");
            Console.WriteLine("  X* = ");
            Console.WriteLine("");
            for (i = 0; i < nvars; i++)
            {
                Console.WriteLine("  " + (i + 1).ToString().PadLeft(8)
                                       + "  " + endpt[i].ToString().PadLeft(14) + "");
            }

            value = woods(endpt, nvars);

            Console.WriteLine("");
            Console.WriteLine("  F(X*) = " + value + "");

        }

        static double rosenbrock(double[] x, int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ROSENBROCK evaluates the Rosenbrock function.
            //
            //  Discussion:
            //
            //    The Hooke and Jeeves algorithm works reasonably well on
            //    Rosenbrock's test function, depending on the value of RHO chosen.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    12 February 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, real X(N), the argument of the function.
            //
            //    Input, integer N, the spatial dimension.
            //
            //    Output, real ROSENBROCK, the value of the function.
            //
        {
            double value;

            value = 100.0 * Math.Pow(x[1] - x[0] * x[0], 2)
                    + Math.Pow(1.0 - x[0], 2);

            return value;
        }

        static double woods(double[] x, int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    WOODS evaluates the Woods function.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    12 February 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, real X(N), the argument of the function.
            //
            //    Input, int N, the spatial dimension.
            //
            //    Output, real WOODS, the value of the function.
            //
        {
            double s1;
            double s2;
            double s3;
            double t1;
            double t2;
            double t3;
            double t4;
            double t5;
            double value;

            s1 = x[1] - x[0] * x[0];
            s2 = 1.0 - x[0];
            s3 = x[1] - 1.0;
            t1 = x[3] - x[2] * x[2];
            t2 = 1.0 - x[2];
            t3 = x[3] - 1.0;
            t4 = s3 + t3;
            t5 = s3 - t3;

            value = 100.0 * s1 * s1
                    + s2 * s2
                    + 90.0 * t1 * t1
                    + t2 * t2
                    + 10.0 * t4 * t4
                    + 0.1 * t5 * t5;

            return value;
        }
    }
}