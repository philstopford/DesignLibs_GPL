using System;
using System.Globalization;
using Burkardt.Polygon;
using Burkardt.Types;

namespace GeometryTest;

public static class PolyloopTest
{
    public static void test0845()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST0845 tests POLYLOOP_ARCLENGTH_ND.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int DIM_NUM = 2;
        const int N = 4;

        int j;
        double[] p =
        {
            0.0, 0.0,
            1.0, 1.0,
            2.0, 0.0,
            0.0, 0.0
        };

        Console.WriteLine("");
        Console.WriteLine("TEST0845");
        Console.WriteLine("  POLYLOOP_ARCLENGTH_ND computes the arclength");
        Console.WriteLine("  of the nodes of a polyloop.");

        double[] s = Geometry.polyloop_arclength_nd(DIM_NUM, N, p);

        Console.WriteLine("");
        Console.WriteLine("           P            Arclength(P)");
        Console.WriteLine("");

        for (j = 0; j <= N; j++)
        {
            string cout = "";
            int j2 = typeMethods.i4_wrap(j, 0, N - 1);
            int i;
            for (i = 0; i < DIM_NUM; i++)
            {
                cout += "  " + p[i + j2 * DIM_NUM].ToString(CultureInfo.InvariantCulture).PadLeft(12);
            }

            Console.WriteLine(cout + "  " + s[j2] + "");
        }
    }

    public static void test0846()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST0846 tests POLYLOOP_POINTS_ND.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int DIM_NUM = 2;
        const int NK = 4;
        const int NT = 12;

        double[] pk =
        {
            0.0, 2.0,
            0.0, 0.0,
            1.0, 0.0,
            1.0, 2.0
        };

        Console.WriteLine("");
        Console.WriteLine("TEST0846");
        Console.WriteLine("  POLYLOOP_POINTS_ND computes points on a polyloop.");

        typeMethods.r8mat_transpose_print(DIM_NUM, NK, pk, "  The defining points:");

        double[] pt = Geometry.polyloop_points_nd(DIM_NUM, NK, pk, NT);

        typeMethods.r8mat_transpose_print(DIM_NUM, NT, pt, "  The computed points:");

    }

}