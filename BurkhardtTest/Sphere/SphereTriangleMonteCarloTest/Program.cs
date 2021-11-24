using System;
using System.Globalization;
using Burkardt.MonomialNS;
using Burkardt.SphereNS;
using Burkardt.Types;

namespace SphereTriangleMonteCarloTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for SPHERE_TRIANGLE_MONTE_CARLO_TEST.
        //
        //  Discussion:
        //
        //    SPHERE_TRIANGLE_MONTE_CARLO_TEST tests the SPHERE_TRIANGLE_MONTE_CARLO 
        //    library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 April 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("SPHERE_TRIANGLE_MONTE_CARLO_TEST");
        Console.WriteLine("  Test the SPHERE_TRIANGLE_MONTE_CARLO library.");

        test01();
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("SPHERE_TRIANGLE_MONTE_CARLO_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //*****************************************************************************/
        //
        //  Purpose:
        //
        //    TEST01 uses SPHERE_TRIANGLE_SAMPLE_01 with an increasing number of points.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 April 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int M = 3;

        int[] e = new int[M];
        int[] e_test =
        {
            0, 0, 0,
            2, 0, 0,
            0, 2, 0,
            0, 0, 2,
            4, 0, 0,
            2, 2, 0,
            0, 0, 4
        };
        int i;
        int k;
        double[] v1 = new double[M];
        double[] v2 = new double[M];
        double[] v3 = new double[M];
        double[] wc = new double[M];

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Estimate monomial integrals over a sphere triangle");
        Console.WriteLine("  using the Monte Carlo method.");

        int seed = 123456789;
        //
        //  Choose three points at random to define a spherical triangle.
        //
        double[] w1 = MonteCarlo.sphere01_sample(1, ref seed);
        double[] w2 = MonteCarlo.sphere01_sample(1, ref seed);
        double[] w3 = MonteCarlo.sphere01_sample(1, ref seed);

        for (i = 0; i < M; i++)
        {
            wc[i] = (w1[i] + w2[i] + w3[i]) / 3.0;
        }

        typeMethods.r8vec_normalize(M, ref wc);
        //
        //  Shrink triangle by factor F.
        //
        double shrink = 2.0;

        for (k = 1; k <= 3; k++)
        {
            shrink /= 2.0;

            for (i = 0; i < M; i++)
            {
                v1[i] = wc[i] + shrink * (w1[i] - wc[i]);
                v2[i] = wc[i] + shrink * (w2[i] - wc[i]);
                v3[i] = wc[i] + shrink * (w3[i] - wc[i]);
            }

            typeMethods.r8vec_normalize(M, ref v1);
            typeMethods.r8vec_normalize(M, ref v2);
            typeMethods.r8vec_normalize(M, ref v3);

            double area = Triangle.sphere01_triangle_vertices_to_area(v1, v2, v3);

            Console.WriteLine("");
            Console.WriteLine("  Vertices of random spherical triangle");
            Console.WriteLine("  with shrink factor = " + shrink + "");
            Console.WriteLine("  and area = " + area + "");
            Console.WriteLine("");
            typeMethods.r8vec_transpose_print(M, v1, "  V1:");
            typeMethods.r8vec_transpose_print(M, v2, "  V2:");
            typeMethods.r8vec_transpose_print(M, v3, "  V3:");
            //
            //  Estimate integrals.
            //
            Console.WriteLine("");
            Console.WriteLine("         N        1              X^2             Y^2" +
                              "             Z^2             X^4           X^2Y^2           Z^4");
            Console.WriteLine("");

            int n = 1;

            while (n <= 4 * 65536)
            {
                double[] x = Triangle.sphere01_triangle_sample(n, v1, v2, v3, ref seed);

                string cout = "  " + n.ToString().PadLeft(8);
                int j;
                for (j = 0; j < 7; j++)
                {
                    for (i = 0; i < M; i++)
                    {
                        e[i] = e_test[i + j * M];
                    }

                    double[] value = Monomial.monomial_value(M, n, e, x);

                    double result = area * typeMethods.r8vec_sum(n, value) / n;
                    cout += "  " + result.ToString(CultureInfo.InvariantCulture).PadLeft(14);
                }

                Console.WriteLine(cout);

                n = 2 * n;
            }
        }
    }
}