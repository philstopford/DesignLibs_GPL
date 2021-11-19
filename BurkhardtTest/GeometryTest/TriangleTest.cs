using System;
using System.Globalization;
using System.Linq;
using Burkardt;
using Burkardt.TriangleNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace GeometryTest;

public static class TriangleTest
{
    public static void triangle_angles_2d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    Geometry.triangle_angles_2d_TEST tests Geometry.triangle_angles_2d;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;

        double[] angle = new double[3];
        int i;
        double[] t =
        {
            0.0, 1.0,
            0.0, 0.0,
            1.0, 0.0
        };

        Console.WriteLine("");
        Console.WriteLine("Geometry.triangle_angles_2d_TEST");
        Console.WriteLine("  For a triangle in 2D,");
        Console.WriteLine("  Geometry.triangle_angles_2d computes the angles;");

        typeMethods.r8mat_transpose_print(DIM_NUM, 3, t, "  Triangle vertices:");

        Geometry.triangle_angles_2d(t, ref angle);

        Console.WriteLine("");
        Console.WriteLine("      Radians      Degrees");
        Console.WriteLine("");
        for (i = 0; i < 3; i++)
        {
            Console.WriteLine("  " + angle[i].ToString().PadLeft(12)
                                   + "  " + Helpers.radians_to_degrees(angle[i]).ToString().PadLeft(12) + "");
        }

    }

    public static void test20605()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST20605 tests Geometry.triangle_angles_3d;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;

        double[] angle = new double[3];
        int i;
        double[] t =
        {
            1.0, 2.0, 3.0,
            2.4142137, 3.4142137, 3.0,
            1.7071068, 2.7071068, 4.0
        };

        Console.WriteLine("");
        Console.WriteLine("TEST206");
        Console.WriteLine("  For a triangle in 3D,");
        Console.WriteLine("  Geometry.triangle_angles_3d computes the angles;");

        typeMethods.r8mat_transpose_print(DIM_NUM, 3, t, "  Triangle vertices:");

        Geometry.triangle_angles_3d(t, ref angle);

        Console.WriteLine("");
        Console.WriteLine("      Radians      Degrees");
        Console.WriteLine("");
        for (i = 0; i < 3; i++)
        {
            Console.WriteLine("  " + angle[i].ToString().PadLeft(12)
                                   + "  " + Helpers.radians_to_degrees(angle[i]).ToString().PadLeft(12) + "");
        }

    }

    public static void test2061()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST2061 tests Geometry.triangle_area_2d;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;

        double area;
        double[] t =
        {
            0.0, 1.0,
            0.0, 0.0,
            1.0, 0.0
        };

        Console.WriteLine("");
        Console.WriteLine("TEST2061");
        Console.WriteLine("  For a triangle in 2D,");
        Console.WriteLine("  Geometry.triangle_area_2d computes the area;");

        typeMethods.r8mat_transpose_print(DIM_NUM, 3, t, "  Triangle vertices:");

        area = Geometry.triangle_area_2d(t);

        Console.WriteLine("  Area is " + area + "");

    }

    public static void test2062()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST2062 tests Geometry.triangle_area_heron;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;

        double area;
        int i;
        int j;
        int jp1;
        double[] s = new double[3];
        double[] t =
        {
            0.0, 1.0,
            0.0, 0.0,
            1.0, 0.0
        };

        Console.WriteLine("");
        Console.WriteLine("TEST2062");
        Console.WriteLine("  For a triangle in 2D,");
        Console.WriteLine("  Geometry.triangle_area_heron computes the area;");

        typeMethods.r8mat_transpose_print(DIM_NUM, 3, t, "  Triangle vertices:");

        for (j = 0; j < 3; j++)
        {
            s[j] = 0.0;

            jp1 = (j + 1) % 3;

            for (i = 0; i < DIM_NUM; i++)
            {
                s[j] += Math.Pow(t[i + j * DIM_NUM] - t[i + jp1 * DIM_NUM], 2);
            }

            s[j] = Math.Sqrt(s[j]);
        }

        typeMethods.r8vec_print(3, s, "  Side lengths:");

        area = Geometry.triangle_area_heron(s);

        Console.WriteLine("  Area is " + area + "");

    }

    public static void test209()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST209 tests Geometry.triangle_area_3d, Geometry.triangle_area_3d_2 and Geometry.triangle_area_3d_3;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;

        double area;
        double[] t =
        {
            1.0, 2.0, 3.0,
            2.4142137, 3.4142137, 3.0,
            1.7071068, 2.7071068, 4.0
        };

        Console.WriteLine("");
        Console.WriteLine("TEST209");
        Console.WriteLine("  For a triangle in 3D:");
        Console.WriteLine("  Geometry.triangle_area_3d   computes the area;");
        Console.WriteLine("  Geometry.triangle_area_3d_2 computes the area;");
        Console.WriteLine("  Geometry.triangle_area_3d_3 computes the area;");

        typeMethods.r8mat_print(DIM_NUM, 3, t, "  Triangle (vertices are columns)");

        area = Geometry.triangle_area_3d(t);

        Console.WriteLine("");
        Console.WriteLine("  Area #1 " + area + "");

        area = Geometry.triangle_area_3d_2(t);

        Console.WriteLine("  Area #2 " + area + "");

        area = Geometry.triangle_area_3d_3(t);

        Console.WriteLine("  Area #3 " + area + "");

    }

    public static void test20655()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST20655 tests TRIANGLE_BARYCENTRIC_2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int TEST_NUM = 7;

        double[] p;
        double[] p_test =
        {
            0.25, 0.25,
            0.75, 0.25,
            1.00, 1.00,
            11.00, 0.50,
            0.00, 1.00,
            0.50, -10.00,
            0.60, 0.60
        };
        double[] t =
        {
            0.0, 1.0,
            0.0, 0.0,
            1.0, 0.0
        };
        int test;
        double[] xsi;

        Console.WriteLine("");
        Console.WriteLine("TEST20655");
        Console.WriteLine("  For a triangle in 2D,");
        Console.WriteLine("  TRIANGLE_BARYCENTRIC_2D converts XY coordinates");
        Console.WriteLine("  to barycentric XSI coordinates;");

        typeMethods.r8mat_transpose_print(DIM_NUM, 3, t, "  Triangle vertices:");

        Console.WriteLine("");
        Console.WriteLine("          P       XSI");
        Console.WriteLine("");

        for (test = 0; test < TEST_NUM; test++)
        {
            p = p_test.Skip(+test * DIM_NUM).ToArray();

            xsi = Geometry.triangle_barycentric_2d(t, p);

            Console.WriteLine("  " + p[0].ToString().PadLeft(8)
                                   + "  " + p[1].ToString().PadLeft(8)
                                   + "  "
                                   + "  " + xsi[0].ToString().PadLeft(8)
                                   + "  " + xsi[1].ToString().PadLeft(8)
                                   + "  " + xsi[2].ToString().PadLeft(8) + "");

        }

    }

    public static void test2066()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST2066 tests TRIANGLE_CENTROID_2D;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int TEST_NUM = 4;

        double[] centroid;
        double[] t;
        double[] t_test =
        {
            0.0, 0.0,
            1.0, 0.0,
            0.0, 1.0,
            0.0, 0.0,
            1.0, 0.0,
            0.5, 0.86602539,
            0.0, 0.0,
            1.0, 0.0,
            0.5, 10.0,
            0.0, 0.0,
            1.0, 0.0,
            10.0, 2.0
        };
        int test;

        Console.WriteLine("");
        Console.WriteLine("TEST2066");
        Console.WriteLine("  For a triangle in 2D:");
        Console.WriteLine("  TRIANGLE_CENTROID_2D computes the centroid.");

        for (test = 0; test < TEST_NUM; test++)
        {
            t = t_test.Skip(+test * DIM_NUM * 3).ToArray();

            typeMethods.r8mat_transpose_print(DIM_NUM, 3, t, "  Triangle vertices:");

            centroid = Geometry.triangle_centroid_2d(t);

            typeMethods.r8vec_print(DIM_NUM, centroid, "  Centroid:");

        }

    }

    public static void test2094()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST2094 tests TRIANGLE_CENTROID_3D;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;

        double[] centroid;
        double[] t =
        {
            1.0, 2.0, 3.0,
            2.4142137, 3.4142137, 3.0,
            1.7071068, 2.7071068, 4.0
        };

        Console.WriteLine("");
        Console.WriteLine("TEST2094");
        Console.WriteLine("  For a triangle in 3D:");
        Console.WriteLine("  TRIANGLE_CENTROID_3D computes the centroid.");

        typeMethods.r8mat_transpose_print(DIM_NUM, 3, t, "  Triangle vertices:");

        centroid = Geometry.triangle_centroid_3d(t);

        typeMethods.r8vec_print(DIM_NUM, centroid, "  Centroid:");

    }

    public static void test2101()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST2101 tests tests TRIANGLE_CIRCUMCENTER_2D and others.
        //
        //  Discussion:
        //
        //    The functions tested include
        //    * TRIANGLE_CIRCUMCENTER_2D;
        //    * TRIANGLE_CIRCUMCENTER_2D_2;
        //    * TRIANGLE_CIRCUMCENTER.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 October 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int M = 2;
        int TEST_NUM = 4;

        int m = M;
        double[] pc;
        double[] t;
        double[] t_test =
        {
            10.0, 5.0,
            11.0, 5.0,
            10.0, 6.0,
            10.0, 5.0,
            11.0, 5.0,
            10.5, 5.86602539,
            10.0, 5.0,
            11.0, 5.0,
            10.5, 15.0,
            10.0, 5.0,
            11.0, 5.0,
            20.0, 7.0
        };
        int test;
        int test_num = TEST_NUM;

        Console.WriteLine("");
        Console.WriteLine("TEST2101");
        Console.WriteLine("  For a triangle in 2D, the circumcenter can be computed by:");
        Console.WriteLine("  TRIANGLE_CIRCUMCENTER_2D;");
        Console.WriteLine("  TRIANGLE_CIRCUMCENTER_2D_2;");
        Console.WriteLine("  TRIANGLE_CIRCUMCENTER (any dimension);");

        for (test = 0; test < test_num; test++)
        {
            t = t_test.Skip(+test * m * 3).ToArray();

            typeMethods.r8mat_transpose_print(m, 3, t, "  Triangle vertices:");

            pc = Geometry.triangle_circumcenter_2d(t);
            typeMethods.r8vec_print(m, pc, "  Circumcenter by TRIANGLE_CIRCUMCENTER_2D:");

            pc = Geometry.triangle_circumcenter_2d_2(t);
            typeMethods.r8vec_print(m, pc, "  Circumcenter by TRIANGLE_CIRCUMCENTER_2D_2:");

            pc = Geometry.triangle_circumcenter(m, t);
            typeMethods.r8vec_print(m, pc, "  Circumcenter by TRIANGLE_CIRCUMCENTER:");
        }

    }

    public static void test21011()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST21011 tests TRIANGLE_CIRCUMCENTER.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 October 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int M1 = 2;
        int TEST_NUM = 4;

        double[] a12;
        int i = 0;
        int j;
        int k;
        int m2;
        double[] o1;
        double[] o2;
        int m1 = M1;
        double[] pc2;
        int seed;
        double[] t1 = new double[M1 * 3];
        double[] t2;
        double[] t_test =
        {
            10.0, 5.0,
            11.0, 5.0,
            10.0, 6.0,
            10.0, 5.0,
            11.0, 5.0,
            10.5, 5.86602539,
            10.0, 5.0,
            11.0, 5.0,
            10.5, 15.0,
            10.0, 5.0,
            11.0, 5.0,
            20.0, 7.0
        };
        int test;
        int test_num = TEST_NUM;

        Console.WriteLine("");
        Console.WriteLine("TEST21011");
        Console.WriteLine("  For a triangle in M dimensions, the circumenter can be computed by:");
        Console.WriteLine("  TRIANGLE_CIRCUMCENTER;");
        //
        //  Vary the dimension.
        //
        for (m2 = 2; m2 <= 5; m2++)
        {
            seed = 123456789;

            Console.WriteLine("");
            Console.WriteLine("  M2 = " + m2 + "");

            t2 = new double[m2 * 3];
            //
            //  Randomly choose a mapping P2 = O2 + A12 * ( P1 - O1 )
            //
            a12 = UniformRNG.r8mat_uniform_01_new(m2, m1, ref seed);
            o1 = UniformRNG.r8vec_uniform_01_new(m1, ref seed);
            o2 = UniformRNG.r8vec_uniform_01_new(m2, ref seed);
            //
            //  Map each M1-dimensional triangle into M2 space.
            //
            for (test = 0; test < test_num; test++)
            {
                for (j = 0; j < 3; j++)
                {
                    for (i = 0; i < m1; i++)
                    {
                        t1[(i + j * m1) % t1.Length] = t_test[(i + j * m1 + test * m1 * 3) % t_test.Length];
                    }
                }

                for (j = 0; j < 3; j++)
                {
                    t1[(i + j * m1) % t1.Length] -= o1[i % o1.Length];
                }

                for (j = 0; j < 3; j++)
                {
                    for (i = 0; i < m2; i++)
                    {
                        t2[(i + j * m2) % t2.Length] = 0.0;
                        for (k = 0; k < m1; k++)
                        {
                            t2[(i + j * m2) % t2.Length] += a12[(i + k * m2) % a12.Length] * t1[(k + j * m1) % t1.Length];
                        }
                    }
                }

                for (j = 0; j < 3; j++)
                {
                    for (i = 0; i < m2; i++)
                    {
                        t2[(i + j * m2) % t2.Length] += o2[i];
                    }
                }

                pc2 = Geometry.triangle_circumcenter(m2, t2);

                typeMethods.r8vec_print(m2, pc2, "  Circumcenter by TRIANGLE_CIRCUMCENTER:");
                Console.WriteLine("");
                Console.WriteLine("  Distances from circumcenter to vertices:");
                Console.WriteLine("");
                for (j = 0; j < 3; j++)
                {
                    Console.WriteLine("  " + typeMethods.r8vec_norm_affine(m2, pc2, t2, v1Index: +j * m2) + "");
                }
            }
        }
    }

    public static void test2067()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST2067 tests TRIANGLE_CIRCUMCIRCLE_2D and TRIANGLE_CIRCUMCIRCLE_2D_2;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int TEST_NUM = 4;

        double[] pc = new double[DIM_NUM];
        double r = 0;
        double[] t;
        double[] t_test =
        {
            0.0, 0.0,
            1.0, 0.0,
            0.0, 1.0,
            0.0, 0.0,
            1.0, 0.0,
            0.5, 0.86602539,
            0.0, 0.0,
            1.0, 0.0,
            0.5, 10.0,
            0.0, 0.0,
            1.0, 0.0,
            10.0, 2.0
        };
        int test;

        Console.WriteLine("");
        Console.WriteLine("TEST2067");
        Console.WriteLine("  For a triangle in 2D:");
        Console.WriteLine("  TRIANGLE_CIRCUMCIRCLE_2D computes the circumcenter.");
        Console.WriteLine("  TRIANGLE_CIRCUMCIRCLE_2D_2 computes the circumcenter.");

        for (test = 0; test < TEST_NUM; test++)
        {
            t = t_test.Skip(+test * DIM_NUM * 3).ToArray();

            typeMethods.r8mat_transpose_print(DIM_NUM, 3, t, "  Triangle vertices:");

            Geometry.triangle_circumcircle_2d(t, ref r, ref pc);

            typeMethods.r8vec_print(DIM_NUM, pc, "  Circumcenter");
            Console.WriteLine("  Circumradius: " + r + "");

            Geometry.triangle_circumcircle_2d_2(t, ref r, ref pc);

            typeMethods.r8vec_print(DIM_NUM, pc, "  Circumcenter2");
            Console.WriteLine("  Circumradius2: " + r + "");
        }

    }

    public static void test21015()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST21015 tests TRIANGLE_CIRCUMRADIUS_2D;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int TEST_NUM = 4;

        double r;
        double[] t;
        double[] t_test =
        {
            0.0, 0.0,
            1.0, 0.0,
            0.0, 1.0,
            0.0, 0.0,
            1.0, 0.0,
            0.5, 0.86602539,
            0.0, 0.0,
            1.0, 0.0,
            0.5, 10.0,
            0.0, 0.0,
            1.0, 0.0,
            10.0, 2.0
        };
        int test;

        Console.WriteLine("");
        Console.WriteLine("TEST21015");
        Console.WriteLine("  For a triangle in 2D:");
        Console.WriteLine("  TRIANGLE_CIRCUMRADIUS_2D computes the circumradius.");

        for (test = 0; test < TEST_NUM; test++)
        {
            t = t_test.Skip(+test * DIM_NUM * 3).ToArray();

            typeMethods.r8mat_transpose_print(DIM_NUM, 3, t, "  Triangle vertices:");

            r = Geometry.triangle_circumradius_2d(t);

            Console.WriteLine("  Circumradius: " + r + "");
        }

    }

    public static void triangle_contains_line_exp_3d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_CONTAINS_LINE_EXP_3D_TEST tests TRIANGLE_CONTAINS_LINE_EXP_3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;

        bool inside = false;
        double[] p1 = {3.0, 0.0, -7.0};
        double[] p2 = {5.0, 1.0, -2.0};
        double[] pint = new double[DIM_NUM];
        double[] t =
        {
            8.0, 4.0, 2.0,
            9.0, 0.0, 5.0,
            2.0, 1.0, 2.0
        };

        Console.WriteLine("");
        Console.WriteLine("TRIANGLE_CONTAINS_LINE_EXP_3D_TEST");
        Console.WriteLine("  TRIANGLE_CONTAINS_LINE_EXP_3D determines whether ");
        Console.WriteLine("  a triangle contains an explicit line in 3D.");

        typeMethods.r8mat_transpose_print(DIM_NUM, 3, t, "  Triangle vertices:");

        typeMethods.r8vec_print(DIM_NUM, p1, "  Line point P1:");
        typeMethods.r8vec_print(DIM_NUM, p2, "  Line point P2:");

        Geometry.triangle_contains_line_exp_3d(t, p1, p2, ref inside, ref pint);

        switch (inside)
        {
            case true:
                Console.WriteLine("");
                Console.WriteLine("  The triangle contains the line.");
                typeMethods.r8vec_print(DIM_NUM, pint, "  Intersection point:");
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("  The triangle does not contain the line.");
                typeMethods.r8vec_print(DIM_NUM, pint, "  The intersection point:");
                break;
        }

        Console.WriteLine("");
        Console.WriteLine("  Expected answer:");
        Console.WriteLine("");
        Console.WriteLine("    The triangle contains the line, and");
        Console.WriteLine("    the intersection point is at:");
        Console.WriteLine("      7, 2, 3.");

    }

    public static void triangle_contains_line_par_3d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_CONTAINS_LINE_PAR_3D_TEST tests TRIANGLE_CONTAINS_LINE_PAR_3D.
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
        int DIM_NUM = 3;

        int dim;
        bool inside = false;
        double norm;
        double[] p0 =
        {
            3.0, 0.0, -7.0
        };
        double[] pd =
        {
            2.0, 1.0, 5.0
        };
        double[] pint = new double[DIM_NUM];
        double[] t =
        {
            8.0, 4.0, 2.0,
            9.0, 0.0, 5.0,
            2.0, 1.0, 2.0
        };

        Console.WriteLine("");
        Console.WriteLine("TRIANGLE_CONTAINS_LINE_PAR_3D_TEST");
        Console.WriteLine("  TRIANGLE_CONTAINS_LINE_PAR_3D determines whether ");
        Console.WriteLine("  a triangle 'contains' a parametric line in 3D.");

        typeMethods.r8mat_transpose_print(DIM_NUM, 3, t, "  Triangle vertices:");

        norm = 0.0;
        for (dim = 0; dim < DIM_NUM; dim++)
        {
            norm += Math.Pow(pd[dim], 2);
        }

        norm = Math.Sqrt(norm);

        for (dim = 0; dim < DIM_NUM; dim++)
        {
            pd[dim] /= norm;
        }

        typeMethods.r8vec_print(DIM_NUM, p0, "  Parametric base point P0:");
        typeMethods.r8vec_print(DIM_NUM, pd, "  Parametric direction PD:");

        Geometry.triangle_contains_line_par_3d(t, p0, pd, ref inside, ref pint);

        switch (inside)
        {
            case true:
                Console.WriteLine("");
                Console.WriteLine("  The triangle contains the line.");
                typeMethods.r8vec_print(DIM_NUM, pint, "  Intersection point:");
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("  The triangle does not contain the line.");
                typeMethods.r8vec_print(DIM_NUM, pint, "  The intersection point:");
                break;
        }

        Console.WriteLine("");
        Console.WriteLine("  Expected answer:");
        Console.WriteLine("");
        Console.WriteLine("    The triangle contains the line, and");
        Console.WriteLine("    the intersection point is at:");
        Console.WriteLine("      ( 7, 2, 3 ).");

    }

    public static void test207()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST207 tests TRIANGLE_CONTAINS_POINT_2D_1, TRIANGLE_CONTAINS_POINT_2D_2, TRIANGLE_CONTAINS_POINT_2D_3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 June 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int TEST_NUM = 7;

        bool inside1;
        bool inside2;
        bool inside3;
        int j;
        double[] p;
        double[] p_test =
        {
            0.25, 0.25,
            0.75, 0.25,
            1.00, 1.00,
            11.00, 0.50,
            0.00, 1.00,
            0.50, -10.00,
            0.60, 0.60
        };
        double[] t =
        {
            0.0, 1.0,
            0.0, 0.0,
            1.0, 0.0
        };
        double[] t2 = new double[DIM_NUM * 3];
        int test;

        Console.WriteLine("");
        Console.WriteLine("TEST207");
        Console.WriteLine("  For a triangle in 2D,");
        Console.WriteLine("  TRIANGLE_CONTAINS_POINT_2D_1 reports if a point");
        Console.WriteLine("  is inside a triangle (and doesn't care about");
        Console.WriteLine("  the ordering of the vertices);");
        Console.WriteLine("  TRIANGLE_CONTAINS_POINT_2D_2 reports if a point ");
        Console.WriteLine("  is inside a triangle (and DOES care about");
        Console.WriteLine("  the ordering of the vertices);");
        Console.WriteLine("  TRIANGLE_CONTAINS_POINT_2D_3 reports if a point");
        Console.WriteLine("  is inside a triangle (and doesn't care about");
        Console.WriteLine("  the ordering of the vertices);");

        typeMethods.r8mat_transpose_print(DIM_NUM, 3, t, "  Triangle vertices:");

        Console.WriteLine("");
        Console.WriteLine("         X         Y  In1  In2  In3");
        Console.WriteLine("");

        for (test = 0; test < TEST_NUM; test++)
        {
            p = p_test.Skip(+test * DIM_NUM).ToArray();

            inside1 = Geometry.triangle_contains_point_2d_1(t, p);
            inside2 = Geometry.triangle_contains_point_2d_2(t, p);
            inside3 = Geometry.triangle_contains_point_2d_3(t, p);

            Console.WriteLine("  " + p[0].ToString().PadLeft(8)
                                   + "  " + p[1].ToString().PadLeft(8)
                                   + "    " + inside1.ToString().PadLeft(1)
                                   + "    " + inside2.ToString().PadLeft(1)
                                   + "    " + inside3.ToString().PadLeft(1) + "");
        }

        //
        //  Make a copy of the triangle with vertices in reverse order.
        //
        Console.WriteLine("");
        Console.WriteLine("  Repeat the test, but reverse the triangle vertex");
        Console.WriteLine("  ordering.");

        for (j = 0; j < 3; j++)
        {
            t2[0 + j * 2] = t[0 + (2 - j) * 2];
            t2[1 + j * 2] = t[1 + (2 - j) * 2];
        }

        typeMethods.r8mat_transpose_print(DIM_NUM, 3, t2, "  Triangle vertices (reversed):");

        Console.WriteLine("");
        Console.WriteLine("         X         Y  In1  In2  In3");
        Console.WriteLine("");

        for (test = 0; test < TEST_NUM; test++)
        {
            p = p_test.Skip(+test * DIM_NUM).ToArray();

            inside1 = Geometry.triangle_contains_point_2d_1(t2, p);
            inside2 = Geometry.triangle_contains_point_2d_2(t2, p);
            inside3 = Geometry.triangle_contains_point_2d_3(t2, p);

            Console.WriteLine("  " + p[0].ToString().PadLeft(8)
                                   + "  " + p[1].ToString().PadLeft(8)
                                   + "    " + inside1.ToString().PadLeft(1)
                                   + "    " + inside2.ToString().PadLeft(1)
                                   + "    " + inside3.ToString().PadLeft(1) + "");
        }

    }

    public static void test2075()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST2075 tests Geometry.triangle_diameter_2d.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int TEST_NUM = 3;

        double diameter;
        double[] t;
        double[] t_test =
        {
            4.0, 2.0,
            1.0, 5.0,
            -2.0, 2.0,
            4.0, 2.0,
            5.0, 4.0,
            6.0, 6.0,
            4.0, 2.0,
            1.0, 5.0,
            4.0, 2.0
        };
        int test;

        Console.WriteLine("");
        Console.WriteLine("TEST2075");
        Console.WriteLine("  Geometry.triangle_diameter_2d computes the diameter of ");
        Console.WriteLine("  the SMALLEST circle around the triangle.");

        for (test = 0; test < TEST_NUM; test++)
        {
            t = t_test.Skip(+test * DIM_NUM * 3).ToArray();

            typeMethods.r8mat_transpose_print(DIM_NUM, 3, t, "  Triangle vertices:");

            diameter = Geometry.triangle_diameter_2d(t);

            Console.WriteLine("  Diameter =       " + diameter + "");
        }

    }

    public static void test208()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST208 tests TRIANGLE_GRIDPOINTS_2D;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int GRID_MAX = 50;

        double[] g = new double[DIM_NUM * GRID_MAX];
        int grid_num = 0;
        int sub_num;
        double[] t =
        {
            0.0, 1.0,
            0.0, 0.0,
            1.0, 0.0
        };

        sub_num = 3;

        Console.WriteLine("");
        Console.WriteLine("TEST208");
        Console.WriteLine("  For a triangle in 2D,");
        Console.WriteLine("  TRIANGLE_GRIDPOINTS_2D produces a set of");
        Console.WriteLine("  gridpoints in or on the triangle.");

        typeMethods.r8mat_transpose_print(DIM_NUM, 3, t, "  Triangle vertices:");

        Geometry.triangle_gridpoints_2d(t, sub_num, GRID_MAX, ref grid_num, ref g);

        Console.WriteLine("");
        Console.WriteLine("  Number of grid points is " + grid_num + "");

        typeMethods.r8mat_print(DIM_NUM, grid_num, g, "  Grid points: ");

    }

    public static void test2102()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST2102 tests TRIANGLE_INCENTER_2D;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int TEST_NUM = 4;

        double[] pc = new double[DIM_NUM];
        double[] t;
        double[] t_test =
        {
            0.0, 0.0,
            1.0, 0.0,
            0.0, 1.0,
            0.0, 0.0,
            1.0, 0.0,
            0.5, 0.86602539,
            0.0, 0.0,
            1.0, 0.0,
            0.5, 10.0,
            0.0, 0.0,
            1.0, 0.0,
            10.0, 2.0
        };
        int test;

        Console.WriteLine("");
        Console.WriteLine("TEST2102");
        Console.WriteLine("  For a triangle in 2D:");
        Console.WriteLine("  TRIANGLE_INCENTER_2D computes the incenter.");

        for (test = 0; test < TEST_NUM; test++)
        {
            t = t_test.Skip(+test * DIM_NUM * 3).ToArray();

            typeMethods.r8mat_transpose_print(DIM_NUM, 3, t, "  Triangle vertices:");

            Geometry.triangle_incenter_2d(t, ref pc);

            typeMethods.r8vec_print(DIM_NUM, pc, "  Incenter");
        }

    }

    public static void test2070()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST2070 tests TRIANGLE_INCIRCLE_2D;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    22 July 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;

        double[] pc = new double[DIM_NUM];
        double r = 0;
        double[] t =
        {
            0.0, 1.0,
            0.0, 0.0,
            1.0, 0.0
        };

        Console.WriteLine("");
        Console.WriteLine("TEST2070");
        Console.WriteLine("  For a triangle in 2D,");
        Console.WriteLine("  TRIANGLE_INCIRCLE_2D computes the incircle;");

        typeMethods.r8mat_transpose_print(DIM_NUM, 3, t, "  Triangle vertices:");

        Geometry.triangle_incircle_2d(t, ref pc, ref r);

        Console.WriteLine("  Incircle center is (" + pc[0] + "  " + pc[1] + ").");
        Console.WriteLine("  Incircle radius is " + r + "");

    }

    public static void test20701()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST20701 tests TRIANGLE_INRADIUS_2D;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    22 July 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;

        double r;
        double[] t =
        {
            0.0, 1.0,
            0.0, 0.0,
            1.0, 0.0
        };

        Console.WriteLine("");
        Console.WriteLine("TEST20701");
        Console.WriteLine("  For a triangle in 2D,");
        Console.WriteLine("  TRIANGLE_INRADIUS_2D computes the inradius;");

        typeMethods.r8mat_transpose_print(DIM_NUM, 3, t, "  Triangle vertices:");

        r = Geometry.triangle_inradius_2d(t);

        Console.WriteLine("  Incircle radius is " + r + "");

    }

    public static void test2104()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST2104 tests TRIANGLE_LATTICE_LAYER_POINT_NEXT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 July 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] c = new int[3];
        int i;
        int j;
        int layer;
        int n = 2;
        bool more;
        int[] v = new int[2];
        string cout = "";

        Console.WriteLine("");
        Console.WriteLine("TEST2104");
        Console.WriteLine("  TRIANGLE_LATTICE_LAYER_POINT_NEXT returns the next");
        Console.WriteLine("  point in a triangle lattice layer defined by:");
        Console.WriteLine("");
        Console.WriteLine("    C[2] - 1 < X[0]/C[0] + X[1]/C[1] <= C[2].");

        c[0] = 2;
        c[1] = 3;
        v[0] = 0;
        v[1] = 0;

        Console.WriteLine("");
        Console.WriteLine("  N = " + n + "");
        cout = "  C =       ";
        for (i = 0; i < n; i++)
        {
            cout += "  " + c[i].ToString().PadLeft(4);
        }

        Console.WriteLine(cout);

        for (layer = 0; layer <= 4; layer++)
        {
            Console.WriteLine("");
            Console.WriteLine("  Layer " + layer + "");
            Console.WriteLine("");

            c[2] = layer;
            more = false;
            i = 0;

            for (;;)
            {
                Geometry.triangle_lattice_layer_point_next(c, ref v, ref more);
                if (!more)
                {
                    Console.WriteLine("  No more.");
                    break;
                }

                i += 1;
                cout += "  " + i.ToString().PadLeft(4);
                for (j = 0; j < n; j++)
                {
                    cout += "  " + v[j].ToString().PadLeft(4);
                }

                Console.WriteLine(cout);
            }
        }
    }


    public static void test2105()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST2105 tests TRIANGLE_LATTICE_POINT_NEXT.
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
        int N = 2;

        int[] c = new int[N + 1];
        int i;
        int j;
        bool more;
        int n = N;
        int[] v = new int[N];
        string cout = "";

        Console.WriteLine("");
        Console.WriteLine("TEST2105");
        Console.WriteLine("  TRIANGLE_LATTICE_POINT_NEXT returns the next lattice");
        Console.WriteLine("  point in a triangle defined by:");
        Console.WriteLine("");
        Console.WriteLine("    0 <= X(1)/C(1) + X(2)/C(2) <= C(3).");

        for (i = 0; i < n + 1; i++)
        {
            c[i] = n + 1 - i;
        }

        for (i = 0; i < n; i++)
        {
            v[i] = 0;
        }

        more = false;

        Console.WriteLine("");
        Console.WriteLine("  N = " + n + "");
        cout = "  C = ";
        for (i = 0; i < n + 1; i++)
        {
            cout += "  " + c[i].ToString().PadLeft(4);
        }

        Console.WriteLine(cout);
        Console.WriteLine("");

        i = 0;

        for (;;)
        {
            Geometry.triangle_lattice_point_next(c, ref v, ref more);
            if (!more)
            {
                Console.WriteLine("  No more.");
                break;
            }

            i += 1;
            cout += "  " + i.ToString().PadLeft(4);
            for (j = 0; j < n; j++)
            {
                cout += "  " + v[j].ToString().PadLeft(4);
            }

            Console.WriteLine(cout);
        }
    }

    public static void test211()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST211 tests TRIANGLE_ORIENTATION_2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int TEST_NUM = 4;

        int i;
        double[] t;
        double[] t_test =
        {
            4.0, 2.0,
            1.0, 5.0,
            -2.0, 2.0,
            1.0, 5.0,
            4.0, 2.0,
            1.0, -1.0,
            1.0, 5.0,
            2.0, 7.0,
            3.0, 9.0,
            1.0, 5.0,
            4.0, 2.0,
            1.0, 5.0
        };
        int test;

        Console.WriteLine("");
        Console.WriteLine("TEST211");
        Console.WriteLine("  TRIANGLE_ORIENTATION_2D determines orientation");
        Console.WriteLine("  of a triangle.");

        for (test = 0; test < TEST_NUM; test++)
        {
            t = t_test.Skip(+test * DIM_NUM * 3).ToArray();

            i = Geometry.triangle_orientation_2d(t);

            typeMethods.r8mat_transpose_print(DIM_NUM, 3, t, "  Triangle vertices:");

            switch (i)
            {
                case 0:
                    Console.WriteLine("  The points are counterclockwise.");
                    break;
                case 1:
                    Console.WriteLine("  The points are clockwise.");
                    break;
                case 2:
                    Console.WriteLine("  The points are colinear.");
                    break;
                case 3:
                    Console.WriteLine("  The points are not distinct.");
                    break;
                default:
                    Console.WriteLine("  The return value makes no sense.");
                    break;
            }
        }

    }

    public static void test2103()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST2103 tests TRIANGLE_ORTHOCENTER_2D;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 July 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int TEST_NUM = 4;

        bool flag = false;
        double[] pc = new double[DIM_NUM];
        double[] t;
        double[] t_test =
        {
            0.0, 0.0,
            1.0, 0.0,
            0.0, 1.0,
            0.0, 0.0,
            1.0, 0.0,
            0.5, 0.86602539,
            0.0, 0.0,
            1.0, 0.0,
            0.5, 10.0,
            0.0, 0.0,
            1.0, 0.0,
            10.0, 2.0
        };
        int test;

        Console.WriteLine("");
        Console.WriteLine("TEST2103");
        Console.WriteLine("  For a triangle in 2D:");
        Console.WriteLine("  TRIANGLE_ORTHOCENTER_2D computes the orthocenter.");

        for (test = 0; test < TEST_NUM; test++)
        {
            t = t_test.Skip(+test * DIM_NUM * 3).ToArray();

            typeMethods.r8mat_transpose_print(DIM_NUM, 3, t, "  Triangle vertices:");

            Geometry.triangle_orthocenter_2d(t, ref pc, ref flag);

            typeMethods.r8vec_print(DIM_NUM, pc, "  Orthocenter");
        }

    }

    public static void test2071()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST2071 tests TRIANGLE_POINT_DIST_2D;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int TEST_NUM = 7;

        double dist;
        double[] p;
        double[] p_test =
        {
            0.25, 0.25,
            0.75, 0.25,
            1.00, 1.00,
            11.00, 0.50,
            0.00, 1.00,
            0.50, -10.00,
            0.60, 0.60
        };
        double[] t =
        {
            0.0, 1.0,
            0.0, 0.0,
            1.0, 0.0
        };
        int test;

        Console.WriteLine("");
        Console.WriteLine("TEST2071");
        Console.WriteLine("  For a triangle in 2D,");
        Console.WriteLine("  TRIANGLE_POINT_DIST_2D computes the distance");
        Console.WriteLine("  to a point;");

        typeMethods.r8mat_transpose_print(DIM_NUM, 3, t, "  Triangle vertices:");

        Console.WriteLine("");
        Console.WriteLine("           P            DIST");
        Console.WriteLine("");

        for (test = 0; test < TEST_NUM; test++)
        {
            p = p_test.Skip(+test * DIM_NUM).ToArray();

            dist = Geometry.triangle_point_dist_2d(t, p);

            Console.WriteLine("  " + p[0].ToString().PadLeft(8)
                                   + "  " + p[1].ToString().PadLeft(8)
                                   + "  " + dist.ToString().PadLeft(8) + "");
        }

    }

    public static void test20715()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST20715 tests TRIANGLE_POINT_DIST_SIGNED_2D;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int TEST_NUM = 7;

        double dist_signed;
        double[] p;
        double[] p_test =
        {
            0.25, 0.25,
            0.75, 0.25,
            1.00, 1.00,
            11.00, 0.50,
            0.00, 1.00,
            0.50, -10.00,
            0.60, 0.60
        };
        double[] t =
        {
            0.0, 1.0,
            0.0, 0.0,
            1.0, 0.0
        };
        int test;

        Console.WriteLine("");
        Console.WriteLine("TEST20715");
        Console.WriteLine("  For a triangle in 2D,");
        Console.WriteLine("  TRIANGLE_POINT_DIST_SIGNED_2D computes signed");
        Console.WriteLine("  distance to a point;");

        typeMethods.r8mat_transpose_print(DIM_NUM, 3, t, "  Triangle vertices:");

        Console.WriteLine("");
        Console.WriteLine("           P       DIST_SIGNED");
        Console.WriteLine("");

        for (test = 0; test < TEST_NUM; test++)
        {
            p = p_test.Skip(+test * DIM_NUM).ToArray();

            dist_signed = Geometry.triangle_point_dist_signed_2d(t, p);

            Console.WriteLine("  " + p[0].ToString().PadLeft(8)
                                   + "  " + p[1].ToString().PadLeft(8)
                                   + "  " + dist_signed.ToString().PadLeft(8) + "");
        }

    }

    public static void triangle_point_dist_3d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_POINT_DIST_3D_TEST tests TRIANGLE_POINT_DIST_3D;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;
        int TEST_NUM = 3;

        double dist;
        double[] p;
        double[] p_test =
        {
            1.0, 2.0, 3.0,
            1.3535534, 2.3535534, 3.0,
            0.0, 0.0, 0.0
        };
        double[] t =
        {
            1.0, 2.0, 3.0,
            2.4142137, 3.4142137, 3.0,
            1.7071068, 2.7071068, 4.0
        };
        int test;

        Console.WriteLine("");
        Console.WriteLine("TRIANGLE_POINT_DIST_3D_TEST");
        Console.WriteLine("  For a triangle in 3D:");
        Console.WriteLine("  TRIANGLE_POINT_DIST_3D computes the distance");
        Console.WriteLine("  to a point;");

        typeMethods.r8mat_transpose_print(DIM_NUM, 3, t, "  Triangle vertices:");

        Console.WriteLine("");
        Console.WriteLine("                   P                          DIST");
        Console.WriteLine("");

        for (test = 0; test < TEST_NUM; test++)
        {
            p = p_test.Skip(+test * DIM_NUM).ToArray();

            dist = Geometry.triangle_point_dist_3d(t, p);

            Console.WriteLine("  " + p[0].ToString().PadLeft(10)
                                   + "  " + p[1].ToString().PadLeft(10)
                                   + "  " + p[2].ToString().PadLeft(10)
                                   + "  " + dist.ToString().PadLeft(12) + "");
        }

    }

    public static void triangle_point_near_2d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_POINT_NEAR_2D_TEST tests TRIANGLE_POINT_NEAR_2D;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int TEST_NUM = 7;

        double dist = 0;
        double[] p;
        double[] p_test =
        {
            0.25, 0.25,
            0.75, 0.25,
            1.00, 1.00,
            11.00, 0.50,
            0.00, 1.00,
            0.50, -10.00,
            0.60, 0.60
        };
        double[] pn = new double[DIM_NUM];
        double[] t =
        {
            0.0, 1.0,
            0.0, 0.0,
            1.0, 0.0
        };
        int test;

        Console.WriteLine("");
        Console.WriteLine("TRIANGLE_POINT_NEAR_2D_TEST");
        Console.WriteLine("  For a triangle in 2D,");
        Console.WriteLine("  TRIANGLE_POINT_NEAR_2D computes the nearest");
        Console.WriteLine("  point to a point.");

        typeMethods.r8mat_transpose_print(DIM_NUM, 3, t, "  Triangle vertices:");

        Console.WriteLine("");
        Console.WriteLine("           P                PN");
        Console.WriteLine("");

        for (test = 0; test < TEST_NUM; test++)
        {
            p = p_test.Skip(+test * DIM_NUM).ToArray();

            Geometry.triangle_point_near_2d(t, p, ref pn, ref dist);

            Console.WriteLine("  " + p[0].ToString().PadLeft(10)
                                   + "  " + p[1].ToString().PadLeft(10)
                                   + "  " + pn[0].ToString().PadLeft(10)
                                   + "  " + pn[1].ToString().PadLeft(10) + "");
        }

    }

    public static void triangle_quality_2d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_QUALITY_2D_TEST tests TRIANGLE_QUALITY_2D;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int TEST_NUM = 4;

        double quality;
        double[] t;
        double[] t_test =
        {
            0.0, 0.0,
            1.0, 0.0,
            0.0, 1.0,
            0.0, 0.0,
            1.0, 0.0,
            0.5, 0.86602539,
            0.0, 0.0,
            1.0, 0.0,
            0.5, 10.0,
            0.0, 0.0,
            1.0, 0.0,
            10.0, 2.0
        };
        int test;

        Console.WriteLine("");
        Console.WriteLine("TRIANGLE_QUALITY_2D_TEST");
        Console.WriteLine("  For a triangle in 2D:");
        Console.WriteLine("  TRIANGLE_QUALITY_2D computes the quality.");

        for (test = 0; test < TEST_NUM; test++)
        {
            t = t_test.Skip(+test * DIM_NUM * 3).ToArray();

            typeMethods.r8mat_transpose_print(DIM_NUM, 3, t, "  Triangle vertices:");

            quality = Geometry.triangle_quality_2d(t);

            Console.WriteLine("  Quality = " + quality + "");
        }

    }

    public static void test212()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST212 tests Geometry.triangle_sample, Geometry.triangle_xy_to_xsi_2d.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int TEST_NUM = 10;

        double[] p = new double[DIM_NUM];
        int seed = 123456789;
        double[] t =
        {
            4.0, 2.0,
            1.0, 5.0,
            -2.0, 2.0
        };
        int test;
        double[] xsi = new double[DIM_NUM + 1];

        Console.WriteLine("");
        Console.WriteLine("TEST212");
        Console.WriteLine("  Geometry.triangle_sample samples a triangle.");
        Console.WriteLine("  Geometry.triangle_xy_to_xsi_2d converts XY to XSI coordinates.");
        Console.WriteLine("");
        Console.WriteLine("  We are computing the XSI coordinates just to verify");
        Console.WriteLine("  that the points are inside the triangle.");

        typeMethods.r8mat_transpose_print(DIM_NUM, 3, t, "  Triangle vertices:");

        Console.WriteLine("");
        Console.WriteLine("  Sample points (X,Y) and (XSI1,XSI2,XSI3) coordinates:");
        Console.WriteLine("");

        for (test = 0; test < TEST_NUM; test++)
        {
            Geometry.triangle_sample(t, 1, ref seed, ref p);
            Geometry.triangle_xy_to_xsi_2d(t, p, xsi);

            Console.WriteLine("  " + p[0].ToString().PadLeft(8)
                                   + "  " + p[1].ToString().PadLeft(8)
                                   + "  " + xsi[0].ToString().PadLeft(8).ToString().PadLeft(8)
                                   + "  " + xsi[1].ToString().PadLeft(8)
                                   + "  " + xsi[2].ToString().PadLeft(8) + "");
        }

    }

    public static void test213()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST213 tests Geometry.triangle_sample, Geometry.triangle_xy_to_xsi_2d, Geometry.triangle_xsi_to_xy_2d.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 April 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int TEST_NUM = 10;

        int i;
        int j;
        double[] p = new double[DIM_NUM];
        double[] p2 = new double[DIM_NUM];
        int seed = 123456789;
        double[] t =
        {
            4.0, 2.0,
            1.0, 5.0,
            -2.0, 2.0
        };
        int test;
        double[] xsi = new double[DIM_NUM + 1];

        Console.WriteLine("");
        Console.WriteLine("TEST213");
        Console.WriteLine("  Geometry.triangle_sample samples a triangle.");
        Console.WriteLine("  Geometry.triangle_xy_to_xsi_2d converts XY to XSI coordinates.");
        Console.WriteLine("  Geometry.triangle_xsi_to_xy_2d converts XSI to XY coordinates.");
        Console.WriteLine("");
        Console.WriteLine("  We verify that (X,Y) -> (XSI1,XSI2,XSI3) -> (X,Y)");
        Console.WriteLine("  works properly.");

        typeMethods.r8mat_transpose_print(DIM_NUM, 3, t, "  Triangle vertices:");

        Console.WriteLine("");
        Console.WriteLine("  Sample points:");
        Console.WriteLine("");

        for (test = 0; test < TEST_NUM; test++)
        {
            switch (test)
            {
                case 1:
                {
                    for (i = 0; i < DIM_NUM; i++)
                    {
                        p[i] = 0.0;
                        for (j = 0; j < 3; j++)
                        {
                            p[i] += t[i + j * DIM_NUM];
                        }

                        p[i] /= 3.0;
                    }

                    break;
                }
                case 2:
                    p[0] = 3.0;
                    p[1] = 0.0;
                    break;
                default:
                    Geometry.triangle_sample(t, 1, ref seed, ref p);
                    break;
            }

            Geometry.triangle_xy_to_xsi_2d(t, p, xsi);

            Geometry.triangle_xsi_to_xy_2d(t, xsi, ref p2);

            Console.WriteLine("");
            Console.WriteLine("  " + p[0].ToString().PadLeft(8)
                                   + "  " + p[1].ToString().PadLeft(8)
                                   + "  " + xsi[0].ToString().PadLeft(8)
                                   + "  " + xsi[1].ToString().PadLeft(8)
                                   + "  " + xsi[2].ToString().PadLeft(8) + "");
            Console.WriteLine("  " + p2[0].ToString().PadLeft(8)
                                   + "  " + p2[1].ToString().PadLeft(8) + "");
        }

    }

}