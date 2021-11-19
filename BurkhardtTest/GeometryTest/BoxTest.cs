using System;
using Burkardt.Cube;
using Burkardt.Types;

namespace GeometryTest;

public static class BoxTest
{
    public static void box_contains_point_2d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BOX_CONTAINS_POINT_2D_TEST tests BOX_CONTAINS_POINT_2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 46;

        int i;
        int j;
        double[] p = new double[2];
        double[] p1 = {-0.1, 0.3};
        double[] p2 = {1.1, 0.9};
        double xhi = 1.2;
        double xlo = -0.3;
        double yhi = 1.4;
        double ylo = -0.1;
        string cout = "";

        Console.WriteLine("");
        Console.WriteLine("BOX_CONTAINS_POINT_2D_TEST");
        Console.WriteLine("  BOX_CONTAINS_POINT_2D reports if a box");
        Console.WriteLine("  contains a point.");
        Console.WriteLine("");
        Console.WriteLine("  We will call the function repeatedly, and draw");
        Console.WriteLine("  a sketch of the box.");
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

                if (Geometry.box_contains_point_2d(p1, p2, p))
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

    public static void box_segment_clip_2d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BOX_SEGMENT_CLIP_2D_TEST tests BOX_SEGMENT_CLIP_2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int TEST_NUM = 5;

        int i;
        int ival;
        double[] p1 = {-10.0, 10.0};
        double[] p2 = {10.0, 20.0};
        double[] pa = new double[DIM_NUM];
        double[] pb = new double[DIM_NUM];
        double[] qa = new double[DIM_NUM];
        double[] qb = new double[DIM_NUM];
        double[] p1_test =
        {
            1.0, 2.0,
            -3.0, 12.0,
            -20.0, 20.0,
            -20.0, 40.0,
            10.0, 40.0
        };
        double[] p2_test =
        {
            8.0, 16.0,
            5.0, 12.0,
            7.0, 20.0,
            0.0, 0.0,
            20.0, 30.0
        };
        int test;
        string cout = "";

        Console.WriteLine("");
        Console.WriteLine("BOX_SEGMENT_CLIP_2D_TEST");
        Console.WriteLine("  BOX_SEGMENT_CLIP_2D clips a line with respect to a box.");
        Console.WriteLine("");
        Console.WriteLine("  The lower left box corner is:");
        Console.WriteLine("");
        for (i = 0; i < DIM_NUM; i++)
        {
            cout += "  " + p1[i].ToString().PadLeft(8);
        }

        Console.WriteLine("");
        Console.WriteLine("");
        Console.WriteLine("  The upper right box corner is:");
        Console.WriteLine("");
        for (i = 0; i < DIM_NUM; i++)
        {
            cout += "  " + p2[i].ToString().PadLeft(8);
        }

        Console.WriteLine(cout);
        cout = "";
        Console.WriteLine("");
        Console.WriteLine("  We list the points PA and PB, and then");
        Console.WriteLine("  the clipped values.");
        Console.WriteLine("");

        for (test = 0; test < TEST_NUM; test++)
        {
            cout = "";
            typeMethods.r8vec_copy(DIM_NUM, p1_test, ref pa, a1index: + test * DIM_NUM);
            typeMethods.r8vec_copy(DIM_NUM, p2_test, ref pb, a1index: + test * DIM_NUM);
            typeMethods.r8vec_copy(DIM_NUM, pa, ref qa);
            typeMethods.r8vec_copy(DIM_NUM, pb, ref qb);

            ival = Geometry.box_segment_clip_2d(p1, p2, qa, qb);

            Console.WriteLine("");
            for (i = 0; i < DIM_NUM; i++)
            {
                cout += "  " + pa[i].ToString().PadLeft(8);
            }

            for (i = 0; i < DIM_NUM; i++)
            {
                cout += "  " + pb[i].ToString().PadLeft(8);
            }

            Console.WriteLine(cout);
            cout = "";
            Console.WriteLine("");

            switch (ival)
            {
                case -1:
                    Console.WriteLine("  Line is outside the box.");
                    break;
                case 0:
                    Console.WriteLine("  Line is inside the box.");
                    break;
                case 1:
                {
                    for (i = 0; i < DIM_NUM; i++)
                    {
                        cout += "  " + qa[i].ToString().PadLeft(8);
                    }

                    break;
                }
                case 2:
                {
                    for (i = 0; i < DIM_NUM; i++)
                    {
                        cout += "  " + "        ";
                    }

                    for (i = 0; i < DIM_NUM; i++)
                    {
                        cout += "  " + qb[i].ToString().PadLeft(8);
                    }

                    Console.WriteLine(cout);
                    cout = "";
                    break;
                }
                case 3:
                {
                    for (i = 0; i < DIM_NUM; i++)
                    {
                        cout += "  " + qa[i].ToString().PadLeft(8);
                    }

                    for (i = 0; i < DIM_NUM; i++)
                    {
                        cout += "  " + qb[i].ToString().PadLeft(8);
                    }

                    Console.WriteLine(cout);
                    cout = "";
                    break;
                }
            }

        }
    }

    public static void box_ray_int_2d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BOX_RAY_INT_2D_TEST tests BOX_RAY_INT_2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int TEST_NUM = 3;

        double[] p1 = {0.0, 0.0};
        double[] p2 = {5.0, 3.0};
        double[] pa = new double[DIM_NUM];
        double[] pa_test =
        {
            3.0, 1.0,
            4.0, 1.0,
            3.0, 1.0
        };
        double[] pb = new double[DIM_NUM];
        double[] pb_test =
        {
            5.0, 5.0,
            3.0, 1.0,
            4.0, 2.0
        };
        double[] pc_test =
        {
            4.0, 3.0,
            0.0, 1.0,
            5.0, 3.0
        };
        double[] pint = new double[DIM_NUM];
        int test;

        Console.WriteLine("");
        Console.WriteLine("BOX_RAY_INT_2D_TEST");
        Console.WriteLine("  BOX_RAY_INT_2D computes the intersection of");
        Console.WriteLine("  a box with coordinate line sides");
        Console.WriteLine("  and a ray whose origin is within the box.");
        Console.WriteLine("");
        Console.WriteLine("  Lower left box corner:");
        Console.WriteLine("");
        Console.WriteLine("  " + p1[0].ToString().PadLeft(12)
                               + "  " + p1[1].ToString().PadLeft(12) + "");
        Console.WriteLine("");
        Console.WriteLine("  Upper right box corner:");
        Console.WriteLine("");
        Console.WriteLine("  " + p2[0].ToString().PadLeft(12)
                               + "  " + p2[1].ToString().PadLeft(12) + "");
        Console.WriteLine("");

        for (test = 0; test < TEST_NUM; test++)
        {
            typeMethods.r8vec_copy(DIM_NUM, pa_test, ref pa, a1index: + test * DIM_NUM);
            typeMethods.r8vec_copy(DIM_NUM, pb_test, ref pb, a1index: + test * DIM_NUM);

            Geometry.box_ray_int_2d(p1, p2, pa, pb, pint);

            Console.WriteLine("");
            Console.WriteLine("  Origin:       " + pa[0].ToString().PadLeft(12)
                                                 + "  " + pa[1].ToString().PadLeft(12) + "");
            Console.WriteLine("  Point 2:      " + pb[0].ToString().PadLeft(12)
                                                 + "  " + pb[1].ToString().PadLeft(12) + "");
            Console.WriteLine("  Intersection: " + pint[0].ToString().PadLeft(12)
                                                 + "  " + pint[1].ToString().PadLeft(12) + "");
            Console.WriteLine("  Correct:      " + pc_test[0 + test * DIM_NUM]
                                                 + "  " + pc_test[1 + test * DIM_NUM].ToString().PadLeft(12) + "");
        }
    }

    public static void box01_contains_point_2d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BOX01_CONTAINS_POINT_2D_TEST tests BOX01_CONTAINS_POINT_2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 46;

        int i;
        int j;
        double[] p = new double[2];
        double xhi = 1.2;
        double xlo = -0.3;
        double yhi = 1.4;
        double ylo = -0.1;
        string cout = "";

        Console.WriteLine("");
        Console.WriteLine("BOX01_CONTAINS_POINT_2D_TEST");
        Console.WriteLine("  BOX01_CONTAINS_POINT_2D reports if the unit box");
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

                if (Geometry.box01_contains_point_2d(p))
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

}