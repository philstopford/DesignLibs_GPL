using System;
using Burkardt.Polygon;
using Burkardt.Types;

namespace GeometryTest;

public static class PolylineTest
{

    public static void test084()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST084 tests POLYLINE_INDEX_POINT_ND and POLYLINE_ARCLENGTH_ND.
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
        int DIM_NUM = 2;
        int N = 4;

        int i;
        int j;
        double[] p =
        {
            0.0, 0.0,
            1.0, 1.0,
            2.0, 0.0,
            0.0, 0.0
        };
        double[] pt;
        double[] s;
        double t;

        t = 2.0;

        string cout = "";

        Console.WriteLine("");
        Console.WriteLine("TEST084");
        Console.WriteLine("  POLYLINE_INDEX_POINT_ND finds a point on a ");
        Console.WriteLine("  polyline with given arclength.");
        Console.WriteLine("  POLYLINE_ARCLENGTH_ND computes the arclength ");
        Console.WriteLine("  of the polyline, and its nodes.");
        Console.WriteLine("");
        Console.WriteLine("  The line we examine is defined by these points:");
        //
        //  The call to POLYLINE_ARCLENGTH_ND is just to help us believe the final result.
        //
        s = Geometry.polyline_arclength_nd(DIM_NUM, N, p);

        Console.WriteLine("");
        Console.WriteLine("           P              Arclength(X,Y)");
        Console.WriteLine("");
        for (j = 0; j < N; j++)
        {
            cout = "";
            for (i = 0; i < DIM_NUM; i++)
            {
                cout += "  " + p[i + j * DIM_NUM].ToString().PadLeft(12);
            }

            Console.WriteLine(cout + "  " + s[j].ToString().PadLeft(12) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  We search for the point with coordinate " + t + "");

        pt = Geometry.polyline_index_point_nd(DIM_NUM, N, p, t);

        typeMethods.r8vec_print(DIM_NUM, pt, "  The computed point:");

    }

    public static void polyline_points_nd_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYLINE_POINTS_ND_TEST tests POLYLINE_POINTS_ND.
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
        int DIM_NUM = 2;
        int NK = 4;
        int NT = 13;

        double[] pk =
        {
            0.0, 1.0,
            0.0, 0.0,
            1.0, 0.0,
            1.0, 2.0
        };
        double[] pt;

        Console.WriteLine("");
        Console.WriteLine("POLYLINE_POINTS_ND_TEST");
        Console.WriteLine("  POLYLINE_POINTS_ND computes points on a polyline.");

        typeMethods.r8mat_transpose_print(DIM_NUM, NK, pk, "  The defining points:");

        pt = Geometry.polyline_points_nd(DIM_NUM, NK, pk, NT);

        typeMethods.r8mat_transpose_print(DIM_NUM, NT, pt, "  The computed points:");

    }


}