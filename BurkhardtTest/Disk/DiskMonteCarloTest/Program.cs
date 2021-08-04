using System;
using Burkardt;
using Burkardt.Types;
using Burkardt.Uniform;

namespace DiskMonteCarloTest
{
    using MonteCarlo = Burkardt.Disk.MonteCarlo;

    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for DISK_MONTE_CARLO_TEST.
            //
            //  Discussion:
            //
            //    DISK_MONTE_CARLO_TEST tests the DISK_MONTE_CARLO library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 July 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] center = new double[2];
            double r;

            Console.WriteLine("");
            Console.WriteLine("DISK_MONTE_CARLO_TEST");
            Console.WriteLine("  Test the DISK_MONTE_CARLO library.");

            disk_area_test();

            center[0] = 0.0;
            center[1] = 0.0;
            r = 1.0;
            disk_sample_test(center, r);

            center[0] = 1.0;
            center[1] = 0.0;
            r = 1.0;
            disk_sample_test(center, r);

            center[0] = 1.0;
            center[1] = 2.0;
            r = 3.0;
            disk_sample_test(center, r);

            Console.WriteLine("");
            Console.WriteLine("DISK_MONTE_CARLO_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void disk_area_test()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DISK_AREA_TEST test DISK_AREA.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 July 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double area;
            double[] center = new double[2];
            double[] data;
            int i;
            double r;
            int seed;

            Console.WriteLine("");
            Console.WriteLine("DISK_AREA_TEST");
            Console.WriteLine("  DISK_AREA computes the area of a disk with");
            Console.WriteLine("  center = (CX,CY) and radius R.");

            seed = 123456789;

            Console.WriteLine("");
            Console.WriteLine("  (   CX        CY     )    R          Area");
            Console.WriteLine("");

            for (i = 1; i <= 10; i++)
            {
                data = UniformRNG.r8vec_uniform_01_new(3, ref seed);
                center[0] = 10.0 * data[0] - 5.0;
                center[1] = 10.0 * data[1] - 5.0;
                r = data[2];
                area = MonteCarlo.disk_area(center, r);
                Console.WriteLine("  (" + center[0].ToString().PadLeft(9)
                                        + ", " + center[1].ToString().PadLeft(9)
                                        + ")  " + r.ToString().PadLeft(9)
                                        + "  " + area + "");
            }
        }

        static void disk_sample_test(double[] center, double r)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DISK_SAMPLE_TEST uses DISK_SAMPLE to estimate integrals.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 July 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int[] e = new int[2];
            int[] e_test =
                {
                    0, 0,
                    2, 0,
                    0, 2,
                    4, 0,
                    2, 2,
                    0, 4,
                    6, 0
                }
                ;
            int i;
            int j;
            int n;
            double result;
            int seed;
            double[] value;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("DISK_SAMPLE_TEST");
            Console.WriteLine("  Use DISK_SAMPLE to estimate integrals in the unit disk");
            Console.WriteLine("  with center (" + center[0]
                                                + "," + center[1]
                                                + ") and radius " + r + "");

            seed = 123456789;
            typeMethods.r8vecNormalData data = new typeMethods.r8vecNormalData();

            Console.WriteLine("");
            Console.WriteLine("         N        1              X^2             Y^2" 
                              + "             X^4             X^2Y^2           Y^4             X^6");
            Console.WriteLine("");

            n = 1;

            while (n <= 65536)
            {
                x = MonteCarlo.disk_sample(center, r, n, ref data, ref seed);

                string cout = "  " + n.ToString().PadLeft(8);
                for (j = 0; j < 7; j++)
                {
                    for (i = 0; i < 2; i++)
                    {
                        e[i] = e_test[i + j * 2];
                    }

                    value = Monomial.monomial_value(2, n, e, x);

                    result = MonteCarlo.disk_area(center, r) * typeMethods.r8vec_sum(n, value) / (double) (n);
                    cout += "  " + result.ToString().PadLeft(14);
                }

                Console.WriteLine(cout);

                n = 2 * n;
            }

            if (
                center[0] == 0.0 &&
                center[1] == 0.0 &&
                r == 1.0)
            {
                Console.WriteLine("");
                string cout = "     Exact";
                for (j = 0; j < 7; j++)
                {
                    for (i = 0; i < 2; i++)
                    {
                        e[i] = e_test[i + j * 2];
                    }

                    result = MonteCarlo.disk01_monomial_integral(e);
                    cout += "  " + result.ToString().PadLeft(14);
                }

                Console.WriteLine(cout);
            }
        }
    }
}