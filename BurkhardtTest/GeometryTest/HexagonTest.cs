using System;
using System.Globalization;
using Burkardt.Hexagon;
using Burkardt.Types;

namespace GeometryTest;

public static class HexagonTest
{
    public static void test0315()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST0315 tests HEXAGON_CONTAINS_POINT_2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 July 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 40;
        int DIM_NUM = 2;

        double[] h =
        {
            0.2, 0.4,
            0.4, 0.2,
            0.8, 0.0,
            1.0, 0.6,
            0.4, 1.0,
            0.2, 0.8
        };
        int i;
        double[] p = new double[DIM_NUM];
        const double xhi = 1.0;
        const double xlo = 0.0;
        const double yhi = 1.0;
        const double ylo = 0.0;

        Console.WriteLine("");
        Console.WriteLine("TEST0315");
        Console.WriteLine("  HEXAGON_CONTAINS_POINT_2D reports if a hexagon");
        Console.WriteLine("  contains a point.");
        Console.WriteLine("");
        Console.WriteLine("  We will call the function repeatedly, and draw");
        Console.WriteLine("  a sketch of an irregular hexagon in the unit square.");
        Console.WriteLine("");

        for (i = 1; i <= N; i++)
        {
            p[1] = ((N - i) * yhi
                    + (i - 1) * ylo)
                   / (N - 1);

            string cout = "  ";
            int j;
            for (j = 1; j <= N; j++)
            {
                p[0] = ((N - j) * xlo
                        + (j - 1) * xhi)
                       / (N - 1);

                if (Geometry.hexagon_contains_point_2d(h, p))
                {
                    cout += '*';
                }
                else
                {
                    cout += '-';
                }
            }

            Console.WriteLine(cout);
        }
    }

    public static void test032()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST032 tests HEXAGON_SHAPE_2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 July 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int DIM_NUM = 2;

        int i;
        double[] p = new double[DIM_NUM];

        Console.WriteLine("");
        Console.WriteLine("TEST032");
        Console.WriteLine("  HEXAGON_SHAPE_2D: points on a unit hexagon.");
        Console.WriteLine("");
        Console.WriteLine("  Angle    X    Y");
        Console.WriteLine("");

        for (i = -10; i <= 370; i += 10)
        {
            double angle = i;
            Geometry.hexagon_shape_2d(angle, ref p);
            Console.WriteLine("  " + angle.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + p[0].ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + p[1].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    public static void test0321()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST0321 tests HEXAGON_VERTICES_2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 July 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;

        double[] p = new double[DIM_NUM * 6];

        Console.WriteLine("");
        Console.WriteLine("TEST0321");
        Console.WriteLine("  HEXAGON_VERTICES_2D: the vertices of the unit hexagon.");

        Geometry.hexagon_vertices_2d(ref p);

        typeMethods.r8mat_transpose_print(DIM_NUM, 6, p, "  Vertices:");

    }

}