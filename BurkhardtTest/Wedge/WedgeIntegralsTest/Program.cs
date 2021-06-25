using System;
using Burkardt;
using Burkardt.Types;
using Integrals = Burkardt.Wedge.Integrals;

namespace WedgeIntegralsTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for WEDGE_INTEGRALS_TEST.
            //
            //  Discussion:
            //
            //    WEDGE_INTEGRALS_TEST tests the WEDGE_INTEGRALS library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    21 August 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("WEDGE_INTEGRALS_TEST");
            Console.WriteLine("  C++ version");
            Console.WriteLine("  Test the WEDGE_INTEGRALS library.");

            test01();

            Console.WriteLine("");
            Console.WriteLine("WEDGE_INTEGRALS_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void test01()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST01 compares exact and estimated monomial integrals.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    21 August 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int e_max = 6;
            int e1;
            int e2;
            int e3;
            double error;
            double exact;
            int[] expon = new int[3];
            int m = 3;
            int n = 500000;
            double q;
            int seed;
            double[] value;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST01");
            Console.WriteLine("  Compare exact and estimated integrals ");
            Console.WriteLine("  over the unit wedge in 3D.");
            //
            //  Get sample points.
            //
            seed = 123456789;
            x = Integrals.wedge01_sample(n, ref seed);

            Console.WriteLine("");
            Console.WriteLine("  Number of sample points used is " + n + "");
            Console.WriteLine("");
            Console.WriteLine("   E1  E2  E3     MC-Estimate      Exact           Error");
            Console.WriteLine("");
            //
            //  Check all monomials up to total degree E_MAX.
            //
            for (e3 = 0; e3 <= e_max; e3++)
            {
                expon[2] = e3;
                for (e2 = 0; e2 <= e_max - e3; e2++)
                {
                    expon[1] = e2;
                    for (e1 = 0; e1 <= e_max - e3 - e2; e1++)
                    {
                        expon[0] = e1;

                        value = Monomial.monomial_value(m, n, expon, x);

                        q = Integrals.wedge01_volume() * typeMethods.r8vec_sum(n, value) / (double) (n);
                        exact = Integrals.wedge01_integral(expon);
                        error = Math.Abs(q - exact);

                        Console.WriteLine(expon[0].ToString().PadLeft(4) + "  "
                            + expon[1].ToString().PadLeft(2) + "  "
                            + expon[2].ToString().PadLeft(2) + "  "
                            + q.ToString().PadLeft(14) + "  "
                            + exact.ToString().PadLeft(14) + "  "
                            + error.ToString().PadLeft(14) + "");
                    }
                }
            }
        }
    }
}