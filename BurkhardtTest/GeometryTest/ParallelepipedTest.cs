using System;
using System.Globalization;
using System.Linq;
using Burkardt.Parallelepiped;

namespace GeometryTest;

public static class ParallelepipedTest
{
    public static void test0495()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST0495 tests PARALLELEPIPED_POINT_DIST_3D.
        //
        //  Discussion:
        //
        //    The points tested are:
        //
        //    1: Center of box.
        //    2: The middle of a face.
        //    3: The middle of an edge.
        //    4: A corner.
        //    5: Close to a face.
        //    6: Close to an edge.
        //    7: Close to a corner.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;
        int TEST_NUM = 7;

        double dist;
        int i;
        double[] p;
        double[] p_test =
        {
            1.0, 4.0, 0.5,
            1.0, 0.0, 0.5,
            0.0, 4.0, 1.0,
            2.0, 8.0, 1.0,
            -0.5, 4.0, 0.5,
            1.0, -1.0, -1.0,
            3.0, 9.0, 2.0
        };
        double[] p1 = {0.0, 0.0, 0.0};
        double[] p2 = {2.0, 0.0, 0.0};
        double[] p3 = {0.0, 8.0, 0.0};
        double[] p4 = {0.0, 0.0, 1.0};
        int test;

        Console.WriteLine("");
        Console.WriteLine("TEST0495");
        Console.WriteLine("  PARALLELEPIPED_POINT_DIST_3D computes the distance");
        Console.WriteLine("  from a point to a box (parallelipiped) in 3D.");
        Console.WriteLine("");
        Console.WriteLine("  The 4 box corners that are specified:");
        Console.WriteLine("");
        string cout = "";
        for (i = 0; i < DIM_NUM; i++)
        {
            cout += "  " + p1[i].ToString(CultureInfo.InvariantCulture).PadLeft(12);
        }

        Console.WriteLine(cout);
        cout = "";
        for (i = 0; i < DIM_NUM; i++)
        {
            cout += "  " + p2[i].ToString(CultureInfo.InvariantCulture).PadLeft(12);
        }

        Console.WriteLine(cout);
        cout = "";
        for (i = 0; i < DIM_NUM; i++)
        {
            cout += "  " + p3[i].ToString(CultureInfo.InvariantCulture).PadLeft(12);
        }

        Console.WriteLine(cout);
        cout = "";
        for (i = 0; i < DIM_NUM; i++)
        {
            cout += "  " + p4[i].ToString(CultureInfo.InvariantCulture).PadLeft(12);
        }

        Console.WriteLine(cout);
        Console.WriteLine("");
        Console.WriteLine(" TEST          P                  Distance to box");
        Console.WriteLine("");

        for (test = 0; test < TEST_NUM; test++)
        {
            p = p_test.Skip(+test * DIM_NUM).ToArray();

            dist = Geometry.parallelepiped_point_dist_3d(p1, p2, p3, p4, p);

            cout = "  " + test.ToString(CultureInfo.InvariantCulture).PadLeft(3);
            for (i = 0; i < DIM_NUM; i++)
            {
                cout += "  " + p[i].ToString(CultureInfo.InvariantCulture).PadLeft(8);
            }

            Console.WriteLine(cout + "  " + dist.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");

        }
    }

}