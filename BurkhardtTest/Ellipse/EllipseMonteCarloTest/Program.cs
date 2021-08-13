using System;
using Burkardt;
using Burkardt.MonomialNS;
using Burkardt.Types;

namespace EllipseMonteCarloTest
{
    using MonteCarlo = Burkardt.Ellipse.MonteCarlo;
    
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for ELLIPSE_MONTE_CARLO_TEST.
            //
            //  Discussion:
            //
            //    ELLIPSE_MONTE_CARLO_TEST tests the ELLIPSE_MONTE_CARLO library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    08 November 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("ELLIPSE_MONTE_CARLO_TEST");
            Console.WriteLine("  Test the ELLIPSE_MONTE_CARLO library.");

            ellipse_area1_test();
            ellipse_area2_test();
            test01();

            Console.WriteLine("");
            Console.WriteLine("ELLIPSE_MONTE_CARLO_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void ellipse_area1_test()

            //*****************************************************************************/
            //
            //  Purpose:
            //
            //    ELLIPSE_AREA1_TEST tests ELLIPSE_AREA1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    08 November 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] a =  {
                5.0, 1.0, 1.0, 2.0
            }
            ;
            double area;
            double r;

            Console.WriteLine("");
            Console.WriteLine("ELLIPSE_AREA1_TEST");
            Console.WriteLine("  C++ version");
            Console.WriteLine("  ELLIPSE_AREA1 computes the area of an ellipse.");

            r = 10.0;

            area = MonteCarlo.ellipse_area1(a, r);

            Console.WriteLine("");
            Console.WriteLine("  R = " + r + "");
            typeMethods.r8mat_print(2, 2, a, "  Matrix A in ellipse definition x*A*x=r^2");
            Console.WriteLine("  Area = " + area + "");
            //
            //  Terminate.
            //
            Console.WriteLine("");
            Console.WriteLine("ELLIPSE_AREA1_TEST");
            Console.WriteLine("  Normal end of execution.");
        }

        static void ellipse_area2_test()

            //*****************************************************************************/
            //
            //  Purpose:
            //
            //    ELLIPSE_AREA2_TEST tests ELLIPSE_AREA2.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    08 November 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double a;
            double area;
            double b;
            double c;
            double d;

            Console.WriteLine("");
            Console.WriteLine("ELLIPSE_AREA2_TEST");
            Console.WriteLine("  C++ version");
            Console.WriteLine("  ELLIPSE_AREA2 computes the area of an ellipse.");

            a = 5.0;
            b = 2.0;
            c = 2.0;
            d = 10.0;

            area = MonteCarlo.ellipse_area2(a, b, c, d);

            Console.WriteLine("");
            Console.WriteLine("  Ellipse: " + a
                + " * x^2 + " + b
                + " * xy + " + c
                + " * y^2 = " + d + "");
            Console.WriteLine("  Area = " + area + "");
            //
            //  Terminate.
            //
            Console.WriteLine("");
            Console.WriteLine("ELLIPSE_AREA2_TEST");
            Console.WriteLine("  Normal end of execution.");
        }

        static void test01()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST01 uses ELLIPSE01_SAMPLE with an increasing number of points.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    20 April 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] a =  {
                9.0, 1.0, 1.0, 4.0
            }
            ;
            int[] e = new int[2];
            int[] e_test =  {
                0, 0,
                1, 0,
                0, 1,
                2, 0,
                1, 1,
                0, 2,
                3, 0
            }
            ;
            int j;
            int n;
            double r = 2.0;
            double result;
            int seed;
            double[] value;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST01");
            Console.WriteLine("  Use ELLIPSE01_SAMPLE to estimate integrals");
            Console.WriteLine("  in the ellipse x' * A * x <= r^2.");

            seed = 123456789;
            typeMethods.r8vecNormalData data = new typeMethods.r8vecNormalData();

            Console.WriteLine("");
            Console.WriteLine("         N        1              X               Y  "
                              + "             X^2               XY             Y^2             X^3");
            Console.WriteLine("");

            n = 1;

            while (n <= 65536)
            {

                x = MonteCarlo.ellipse_sample(n, a, r, ref data, ref seed);

                string cout = "  " + n.ToString().PadLeft(8);

                for (j = 0; j < 7; j++)
                {
                    e[0] = e_test[0 + j * 2];
                    e[1] = e_test[1 + j * 2];

                    value = Monomial.monomial_value(2, n, e, x);

                    result = MonteCarlo.ellipse_area1(a, r) * typeMethods.r8vec_sum(n, value)
                             / (double) (n);

                    cout += "  " + result.ToString().PadLeft(14);
                }

                Console.WriteLine("");
                n = 2 * n;
            }
        }
    }
}