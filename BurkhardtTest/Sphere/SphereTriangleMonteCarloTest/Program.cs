using System;
using Burkardt.MonomialNS;
using Burkardt.SphereNS;
using Burkardt.Types;

namespace SphereTriangleMonteCarloTest;

internal class Program
{
    private static void Main(string[] args)
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
        int M = 3;

        double area;
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
        int j;
        int k;
        int m = M;
        int n;
        double result;
        int seed;
        double shrink;
        double[] v1 = new double[M];
        double[] v2 = new double[M];
        double[] v3 = new double[M];
        double[] wc = new double[M];
        double[] w1;
        double[] w2;
        double[] w3;
        double[] value;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Estimate monomial integrals over a sphere triangle");
        Console.WriteLine("  using the Monte Carlo method.");

        seed = 123456789;
        //
        //  Choose three points at random to define a spherical triangle.
        //
        w1 = MonteCarlo.sphere01_sample(1, ref seed);
        w2 = MonteCarlo.sphere01_sample(1, ref seed);
        w3 = MonteCarlo.sphere01_sample(1, ref seed);

        for (i = 0; i < m; i++)
        {
            wc[i] = (w1[i] + w2[i] + w3[i]) / 3.0;
        }

        typeMethods.r8vec_normalize(m, ref wc);
        //
        //  Shrink triangle by factor F.
        //
        shrink = 2.0;

        for (k = 1; k <= 3; k++)
        {
            shrink /= 2.0;

            for (i = 0; i < m; i++)
            {
                v1[i] = wc[i] + shrink * (w1[i] - wc[i]);
                v2[i] = wc[i] + shrink * (w2[i] - wc[i]);
                v3[i] = wc[i] + shrink * (w3[i] - wc[i]);
            }

            typeMethods.r8vec_normalize(m, ref v1);
            typeMethods.r8vec_normalize(m, ref v2);
            typeMethods.r8vec_normalize(m, ref v3);

            area = Triangle.sphere01_triangle_vertices_to_area(v1, v2, v3);

            Console.WriteLine("");
            Console.WriteLine("  Vertices of random spherical triangle");
            Console.WriteLine("  with shrink factor = " + shrink + "");
            Console.WriteLine("  and area = " + area + "");
            Console.WriteLine("");
            typeMethods.r8vec_transpose_print(m, v1, "  V1:");
            typeMethods.r8vec_transpose_print(m, v2, "  V2:");
            typeMethods.r8vec_transpose_print(m, v3, "  V3:");
            //
            //  Estimate integrals.
            //
            Console.WriteLine("");
            Console.WriteLine("         N        1              X^2             Y^2" +
                              "             Z^2             X^4           X^2Y^2           Z^4");
            Console.WriteLine("");

            n = 1;

            while (n <= 4 * 65536)
            {
                x = Triangle.sphere01_triangle_sample(n, v1, v2, v3, ref seed);

                string cout = "  " + n.ToString(CultureInfo.InvariantCulture).PadLeft(8);
                for (j = 0; j < 7; j++)
                {
                    for (i = 0; i < m; i++)
                    {
                        e[i] = e_test[i + j * m];
                    }

                    value = Monomial.monomial_value(m, n, e, x);

                    result = area * typeMethods.r8vec_sum(n, value) / n;
                    cout += "  " + result.ToString(CultureInfo.InvariantCulture).PadLeft(14);
                }

                Console.WriteLine(cout);

                n = 2 * n;
            }
        }
    }
}