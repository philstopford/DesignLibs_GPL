using System;
using System.Linq;
using Burkardt.Parallelogram;
using Burkardt.Types;

namespace GeometryTest;

public static class ParallelogramTest
{

    public static void test0477()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST0477 tests PARALLELOGRAM_AREA_2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 May 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double area;
        double[] p =
        {
            2.0, 7.0,
            5.0, 7.0,
            6.0, 9.0,
            3.0, 9.0
        };

        Console.WriteLine("");
        Console.WriteLine("TEST0477");
        Console.WriteLine("  PARALLELOGRAM_AREA_2D finds the area of a");
        Console.WriteLine("  parallelogram in 2D.");

        typeMethods.r8mat_transpose_print(2, 4, p, "  Vertices:");

        area = Geometry.parallelogram_area_2d(p);

        Console.WriteLine("");
        Console.WriteLine("  AREA = " + area + "");
    }

    public static void test0478()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST0478 tests PARALLELOGRAM_AREA_3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 May 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double area;
        double[] p =
        {
            1.0, 2.0, 3.0,
            2.4142137, 3.4142137, 3.0,
            1.7071068, 2.7071068, 4.0,
            0.2928931, 0.2928931, 4.0
        };

        Console.WriteLine("");
        Console.WriteLine("TEST0478");
        Console.WriteLine("  PARALLELOGRAM_AREA_3D finds the area of a");
        Console.WriteLine("  parallelogram in 3D.");

        typeMethods.r8mat_transpose_print(3, 4, p, "  Vertices:");

        area = Geometry.parallelogram_area_3d(p);

        Console.WriteLine("");
        Console.WriteLine("  AREA = " + area + "");
    }

    public static void parallelogram_contains_point_2d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PARALLELOGRAM_CONTAINS_POINT_2D_TEST tests PARALLELOGRAM_CONTAINS_POINT_2D.
        //
        //  Discussion:
        //
        //    The four points are In, Out, Out, and Out.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 July 2005
        //
    {
        int DIM_NUM = 2;
        int TEST_NUM = 4;

        int i;
        bool inside;
        double[] p;
        double[] p_test =
        {
            1.0, 0.5,
            2.0, 0.0,
            0.5, -0.1,
            0.1, 0.5
        };
        double[] p1 = {0.0, 0.0};
        double[] p2 = {1.0, 0.0};
        double[] p3 = {1.0, 1.0};
        int test;

        Console.WriteLine("");
        Console.WriteLine("PARALLELOGRAM_CONTAINS_POINT_2D_TEST");
        Console.WriteLine("  PARALLELOGRAM_CONTAINS_POINT_2D determines if a point ");
        Console.WriteLine("  is within a parallelogram in 2D.");
        Console.WriteLine("");
        Console.WriteLine("       P     Inside?");
        Console.WriteLine("");

        for (test = 0; test < TEST_NUM; test++)
        {
            p = p_test.Skip(+test * DIM_NUM).ToArray();

            inside = Geometry.parallelogram_contains_point_2d(p1, p2, p3, p);

            string cout = "";
            for (i = 0; i < DIM_NUM; i++)
            {
                cout += "  " + p[i].ToString().PadLeft(12);
            }

            Console.WriteLine(cout + "  " + inside + "");
        }
    }

    public static void parallelogram_contains_point_2d_test2()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PARALLELOGRAM_CONTAINS_POINT_2D_TEST2 tests PARALLELOGRAM_CONTAINS_POINT_2D
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int N = 51;

        int i;
        int j;
        double[] p = new double[DIM_NUM];
        double[] p1 = {0.2, 0.0};
        double[] p2 = {0.4, 0.6};
        double[] p3 = {0.6, 0.4};
        double xhi = 1.0;
        double xlo = 0.0;
        double yhi = 1.0;
        double ylo = 0.0;
        string cout = "";

        Console.WriteLine("");
        Console.WriteLine("PARALLELOGRAM_CONTAINS_POINT_2D_TEST2");
        Console.WriteLine("  PARALLELOGRAM_CONTAINS_POINT_2D reports if a parallelogram");
        Console.WriteLine("  contains a point.");
        Console.WriteLine("");
        Console.WriteLine("  We will call the function repeatedly, and draw");
        Console.WriteLine("  a sketch of the unit square.");
        Console.WriteLine("");

        for (i = 1; i <= N; i++)
        {
            p[1] = ((N - i) * yhi
                    + (i - 1) * ylo)
                   / (N - 1);

            cout = "  ";
            for (j = 1; j <= N; j++)
            {
                p[0] = ((N - j) * xlo
                        + (j - 1) * xhi)
                       / (N - 1);

                if (Geometry.parallelogram_contains_point_2d(p1, p2, p3, p))
                {
                    cout += "*";
                }
                else
                {
                    cout += "-";
                }
            }

            Console.WriteLine(cout);
        }
    }

    public static void parallelogram_contains_point_3d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PARALLELOGRAM_CONTAINS_POINT_3D_TEST tests PARALLELOGRAM_CONTAINS_POINT_3D.
        //
        //  Discussion:
        //
        //    The points are In, Out, Out, Out, Out
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;
        int TEST_NUM = 5;

        int i;
        bool inside;
        double[] p;
        double[] p_test =
        {
            1.0, 1.0, 0.5,
            3.0, 3.0, 0.0,
            0.5, 0.5, -0.1,
            0.1, 0.1, 0.5,
            1.5, 1.6, 0.5
        };
        double[] p1 = {0.0, 0.0, 0.0};
        double[] p2 = {2.0, 2.0, 0.0};
        double[] p3 = {1.0, 1.0, 1.0};
        int test;

        Console.WriteLine("");
        Console.WriteLine("PARALLELOGRAM_CONTAINS_POINT_3D_TEST");
        Console.WriteLine("  PARALLELOGRAM_CONTAINS_POINT_3D determines if a point");
        Console.WriteLine("  is within a parallelogram in 3D.");
        Console.WriteLine("");
        Console.WriteLine("           P           Inside?");
        Console.WriteLine("");

        for (test = 0; test < TEST_NUM; test++)
        {
            p = p_test.Skip(+test * DIM_NUM).ToArray();

            inside = Geometry.parallelogram_contains_point_3d(p1, p2, p3, p);

            string cout = "";
            for (i = 0; i < DIM_NUM; i++)
            {
                cout += "  " + p[i].ToString().PadLeft(12);
            }

            Console.WriteLine(cout + "  " + inside + "");
        }
    }

}