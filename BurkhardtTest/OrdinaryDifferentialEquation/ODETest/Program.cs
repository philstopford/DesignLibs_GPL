using System;
using Burkardt.ODENS;

namespace ODETest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for ODE_TEST.
            //
            //  Discussion:
            //
            //    ODE_TEST tests the ODE library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    15 January 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("ODE_TEST");
            
            Console.WriteLine("  Test the ODE library.");

            test01();
            test02();

            Console.WriteLine("");
            Console.WriteLine("ODE_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void test01()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST01 tests ODE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    15 November 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double abserr;
            int i;
            int iflag;
            int[] iwork = new int[5];
            int neqn = 2;
            double pi = 3.141592653589793;
            double relerr;
            int step_num = 12;
            double t;
            double tout;
            double[] work;
            double[] y;

            Console.WriteLine("");
            Console.WriteLine("TEST01");
            Console.WriteLine("  ODE solves a system of ordinary differential");
            Console.WriteLine("  equations.");
            Console.WriteLine("");
            Console.WriteLine("      T           Y(1)         Y(2)");
            Console.WriteLine("");

            abserr = 0.00001;
            relerr = 0.00001;

            iflag = 1;

            t = 0.0;
            y = new double[neqn];
            y[0] = 1.0;
            y[1] = 0.0;

            Console.WriteLine("  " + t.ToString().PadLeft(8)
                                   + "  " + y[0].ToString().PadLeft(14)
                                   + "  " + y[1].ToString().PadLeft(14) + "");

            work = new double[100 + 21 * neqn];

            for (i = 1; i <= step_num; i++)
            {
                tout = (double) (i) * 2.0 * pi / (double) (step_num);

                ODE.ode(f01, neqn, y, ref t, tout, relerr, abserr, ref iflag, ref work, ref iwork);

                if (iflag != 2)
                {
                    Console.WriteLine("");
                    Console.WriteLine("TEST01 - Fatal error!");
                    Console.WriteLine("  ODE returned IFLAG = " + iflag + "");
                    break;
                }

                Console.WriteLine("  " + t.ToString().PadLeft(8)
                                       + "  " + y[0].ToString().PadLeft(14)
                                       + "  " + y[1].ToString().PadLeft(14) + "");
            }
        }

        static void test02()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST02 tests ODE by integrating in the NEGATIVE time direction.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    02 February 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double abserr;
            int i;
            int iflag;
            int[] iwork = new int[5];
            int neqn = 2;
            double pi = 3.141592653589793;
            double relerr;
            int step_num = 12;
            double t;
            double tout;
            double[] work;
            double[] y;

            Console.WriteLine("");
            Console.WriteLine("TEST02");
            Console.WriteLine("  ODE solves a system of ordinary differential");
            Console.WriteLine("  equations.");
            Console.WriteLine("");
            Console.WriteLine("  In this example, we integrate in the negative");
            Console.WriteLine("  time direction.");
            Console.WriteLine("");
            Console.WriteLine("      T           Y(1)         Y(2)");
            Console.WriteLine("");

            abserr = 0.00001;
            relerr = 0.00001;

            iflag = 1;

            t = 0.0;
            y = new double[neqn];
            y[0] = 1.0;
            y[1] = 0.0;

            Console.WriteLine("  " + t.ToString().PadLeft(8)
                                   + "  " + y[0].ToString().PadLeft(14)
                                   + "  " + y[1].ToString().PadLeft(14) + "");

            work = new double[100 + 21 * neqn];

            for (i = 1; i <= step_num; i++)
            {
                tout = -(double) (i) * 2.0 * pi / (double) (step_num);

                ODE.ode(f01, neqn, y, ref t, tout, relerr, abserr, ref iflag, ref work, ref iwork);

                if (iflag != 2)
                {
                    Console.WriteLine("");
                    Console.WriteLine("TEST02 - Fatal error!");
                    Console.WriteLine("  ODE returned IFLAG = " + iflag + "");
                    break;
                }

                Console.WriteLine("  " + t.ToString().PadLeft(8)
                                       + "  " + y[0].ToString().PadLeft(14)
                                       + "  " + y[1].ToString().PadLeft(14) + "");
            }

        }

        static double[] f01(double t, double[] y, int yIndex, double[] yp, int ypIndex)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    F01 supplies the right hand side of the ODE for problem 1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    15 January 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double T, the time.
            //
            //    Input, double Y[], the dependent variable.
            //
            //    Output, double YP[], the value of the derivative.
            //
        {
            yp[(0 + ypIndex) % yp.Length] = y[(1 + yIndex) % y.Length];
            yp[(1 + ypIndex) % yp.Length] = -y[(0 + yIndex) % y.Length];

            return yp;
        }
    }
}