using System;
using Burkardt.TriangleNS;

namespace TriangleInterpolateTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_INTERPOLATE_TEST tests the TRIANGLE_INTERPOLATE library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 January 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("TRIANGLE_INTERPOLATE_TEST");
        Console.WriteLine("  Test the TRIANGLE_INTERPOLATE library.");

        triangle_interpolate_linear_test();

        Console.WriteLine("");
        Console.WriteLine("TRIANGLE_INTERPOLATE_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void triangle_interpolate_linear_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_INTERPOLATE_LINEAR_TEST tests TRIANGLE_INTERPOLATE_LINEAR.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 January 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int j;
        int m = 3;
        int n = 10;
        double[] p;
        double[] p1 = {0.0, 1.0};
        double[] p2 = {5.0, 0.0};
        double[] p3 = {4.0, 4.0};
        int seed;
        double[] v;
        double[] v1 = {1.0, 0.0, 0.0};
        double[] v2 = {0.0, 1.0, 0.0};
        double[] v3 = {0.0, 0.0, 1.0};

        Console.WriteLine("");
        Console.WriteLine("TRIANGLE_INTERPOLATE_LINEAR_TEST");
        Console.WriteLine("  TRIANGLE_INTERPOLATE_LINEAR uses linear interpolation");
        Console.WriteLine("  on vertex data to estimate values in the interior.");
        //
        //  Get N sample points inside the triangle.
        //
        seed = 123456789;
        p = Burkardt.Uniform.Triangle.uniform_in_triangle_map1(p1, p2, p3, n, ref seed);
        //
        //  Request an intepolated value for R, G and B at each point.
        //
        v = Interpolate.triangle_interpolate_linear(m, n, p1, p2, p3, p, v1, v2, v3);
        //
        //  Report the data.
        //
        Console.WriteLine("");
        Console.WriteLine("       X               Y               V(1)            V(2)            V(3)");
        Console.WriteLine("");
        Console.WriteLine("  " + p1[0].ToString().PadLeft(14)
                               + "  " + p1[1].ToString().PadLeft(14)
                               + "  " + v1[0].ToString().PadLeft(14)
                               + "  " + v1[1].ToString().PadLeft(14)
                               + "  " + v1[2].ToString().PadLeft(14) + "");
        Console.WriteLine("  " + p2[0].ToString().PadLeft(14)
                               + "  " + p2[1].ToString().PadLeft(14)
                               + "  " + v2[0].ToString().PadLeft(14)
                               + "  " + v2[1].ToString().PadLeft(14)
                               + "  " + v2[2].ToString().PadLeft(14) + "");
        Console.WriteLine("  " + p3[0].ToString().PadLeft(14)
                               + "  " + p3[1].ToString().PadLeft(14)
                               + "  " + v3[0].ToString().PadLeft(14)
                               + "  " + v3[1].ToString().PadLeft(14)
                               + "  " + v3[2].ToString().PadLeft(14) + "");
        Console.WriteLine("");
        for (j = 0; j < n; j++)
        {
            Console.WriteLine("  " + p[(0 + j * m) % p.Length].ToString().PadLeft(14)
                                   + "  " + p[(1 + j * m) % p.Length].ToString().PadLeft(14)
                                   + "  " + v[(0 + j * m) % v.Length].ToString().PadLeft(14)
                                   + "  " + v[(1 + j * m) % v.Length].ToString().PadLeft(14)
                                   + "  " + v[(2 + j * m) % v.Length].ToString().PadLeft(14) + "");
        }
    }
}