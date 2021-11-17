using System;
using System.Linq;
using Burkardt;
using Burkardt.Polygon;
using Burkardt.Types;

namespace GeometryTest;

public static class PolygonTest
{
    public static void test0755()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST0755 tests POLYGON_1_2D, POLYGON_X_2D, POLYGON_Y_2D, POLYGON_XX_2D etc.
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

        double result;
        double[] v =
        {
            0.0, 0.0,
            1.0, 0.0,
            1.0, 1.0,
            0.0, 1.0
        };

        Console.WriteLine("");
        Console.WriteLine("TEST0755");
        Console.WriteLine("  For a polygon in 2D:");
        Console.WriteLine("  POLYGON_1_2D integrates 1");
        Console.WriteLine("  POLYGON_X_2D integrates X");
        Console.WriteLine("  POLYGON_Y_2D integrates Y");
        Console.WriteLine("  POLYGON_XX_2D integrates X*X");
        Console.WriteLine("  POLYGON_XY_2D integrates X*Y");
        Console.WriteLine("  POLYGON_YY_2D integrates Y*Y");

        typeMethods.r8mat_transpose_print(DIM_NUM, N, v, "  The polygon vertices:");

        Console.WriteLine("");
        Console.WriteLine("  F(X,Y)    Integral");
        Console.WriteLine("");

        result = Geometry.polygon_1_2d(N, v);
        Console.WriteLine("    1  " + "  " + result + "");

        result = Geometry.polygon_x_2d(N, v);
        Console.WriteLine("    X  " + "  " + result + "");

        result = Geometry.polygon_y_2d(N, v);
        Console.WriteLine("    Y  " + "  " + result + "");

        result = Geometry.polygon_xx_2d(N, v);
        Console.WriteLine("  X*X  " + "  " + result + "");

        result = Geometry.polygon_xy_2d(N, v);
        Console.WriteLine("  X*Y  " + "  " + result + "");

        result = Geometry.polygon_yy_2d(N, v);
        Console.WriteLine("  Y*Y  " + "  " + result + "");

    }

    public static void polygon_angles_2d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_ANGLES_2D_TEST tests POLYGON_ANGLES_2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 March 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int N = 6;

        double[] angle;
        int i;
        double[] v =
        {
            0.0, 0.0,
            1.0, 0.0,
            2.0, 1.0,
            3.0, 0.0,
            3.0, 2.0,
            1.0, 2.0
        };

        Console.WriteLine("");
        Console.WriteLine("POLYGON_ANGLES_2D_TEST");
        Console.WriteLine("  For a polygon in 2D:");
        Console.WriteLine("  POLYGON_ANGLES_2D computes the angles.");

        Console.WriteLine("");
        Console.WriteLine("  Number of polygonal vertices = " + N + "");

        typeMethods.r8mat_transpose_print(DIM_NUM, N, v, "  The polygon vertices:");

        angle = Geometry.polygon_angles_2d(N, v);

        Console.WriteLine("");
        Console.WriteLine("  Polygonal angles in degrees:");
        Console.WriteLine("");

        for (i = 0; i < N; i++)
        {
            Console.WriteLine("  " + i.ToString().PadLeft(6)
                                   + "  " + Helpers.radians_to_degrees(angle[i]).ToString().PadLeft(14) + "");
        }
    }

    public static void test076()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST076 tests POLYGON_AREA_2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int TEST_NUM = 2;

        double area;
        double area_exact;
        double[] area_exact_test = {2.0, 6.0};
        int n;
        int[] n_test = {4, 8};
        int test;
        double[] v;
        double[] v1 =
        {
            1.0, 0.0,
            2.0, 1.0,
            1.0, 2.0,
            0.0, 1.0
        };
        double[] v2 =
        {
            0.0, 0.0,
            3.0, 0.0,
            3.0, 3.0,
            2.0, 3.0,
            2.0, 1.0,
            1.0, 1.0,
            1.0, 2.0,
            0.0, 2.0
        };

        Console.WriteLine("");
        Console.WriteLine("TEST076");
        Console.WriteLine("  For a polygon in 2D:");
        Console.WriteLine("  POLYGON_AREA_2D computes the area.");

        for (test = 0; test < TEST_NUM; test++)
        {
            n = n_test[test];
            area_exact = area_exact_test[test];
            v = new double[DIM_NUM * n];

            switch (test)
            {
                case 0:
                    typeMethods.r8mat_copy(DIM_NUM, n, v1, ref v);
                    break;
                case 1:
                    typeMethods.r8mat_copy(DIM_NUM, n, v2, ref v);
                    break;
            }

            Console.WriteLine("");
            Console.WriteLine("  Number of polygonal vertices = " + n + "");

            typeMethods.r8mat_transpose_print(DIM_NUM, n, v, "  The polygon vertices:");

            area = Geometry.polygon_area_2d(n, v);

            Console.WriteLine("");
            Console.WriteLine("  Exact area is        " + area_exact + "");
            Console.WriteLine("  The computed area is " + area + "");

        }

    }

    public static void test0765()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST0765 tests POLYGON_AREA_2D_2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int TEST_NUM = 2;

        double area;
        double area_exact;
        double[] area_exact_test = {2.0, 6.0};
        int n;
        int[] n_test = {4, 8};
        int test;
        double[] v;
        double[] v1 =
        {
            1.0, 0.0,
            2.0, 1.0,
            1.0, 2.0,
            0.0, 1.0
        };
        double[] v2 =
        {
            0.0, 0.0,
            3.0, 0.0,
            3.0, 3.0,
            2.0, 3.0,
            2.0, 1.0,
            1.0, 1.0,
            1.0, 2.0,
            0.0, 2.0
        };

        Console.WriteLine("");
        Console.WriteLine("TEST0765");
        Console.WriteLine("  For a polygon in 2D:");
        Console.WriteLine("  POLYGON_AREA_2D_2 computes the area.");

        for (test = 0; test < TEST_NUM; test++)
        {
            n = n_test[test];
            area_exact = area_exact_test[test];
            v = new double[DIM_NUM * n];

            switch (test)
            {
                case 0:
                    typeMethods.r8mat_copy(DIM_NUM, n, v1, ref v);
                    break;
                case 1:
                    typeMethods.r8mat_copy(DIM_NUM, n, v2, ref v);
                    break;
            }

            Console.WriteLine("");
            Console.WriteLine("  Number of polygonal vertices = " + n + "");

            typeMethods.r8mat_transpose_print(DIM_NUM, n, v, "  The polygon vertices:");

            area = Geometry.polygon_area_2d_2(n, v);

            Console.WriteLine("");
            Console.WriteLine("  Exact area is        " + area_exact + "");
            Console.WriteLine("  The computed area is " + area + "");

        }

    }

    public static void test078()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST078 tests POLYGON_AREA_3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;
        int TEST_NUM = 2;

        double area;
        double area_exact;
        double[] area_exact_test = {2.4494898, 6.0};
        int n;
        int[] n_test = {4, 8};
        double[] normal = new double[DIM_NUM];
        int test;
        double[] v;
        double[] v1 =
        {
            1.0, 0.0, 0.0,
            2.0, 1.0, 1.0,
            1.0, 2.0, 1.0,
            0.0, 1.0, 0.0
        };
        double[] v2 =
        {
            0.00000, 0.00000, 0.00000,
            2.62679, 1.26009, -0.715657,
            1.48153, 3.97300, -0.142512,
            0.605932, 3.55297, 0.0960401,
            1.36944, 1.74437, -0.286056,
            0.493842, 1.32433, -0.0475041,
            0.112090, 2.22864, 0.143544,
            -0.763505, 1.80861, 0.382097
        };

        Console.WriteLine("");
        Console.WriteLine("TEST078");
        Console.WriteLine("  For a polygon in 3D:");
        Console.WriteLine("  POLYGON_AREA_3D computes the area;");

        for (test = 0; test < TEST_NUM; test++)
        {
            area_exact = area_exact_test[test];
            n = n_test[test];

            v = new double[DIM_NUM * n];

            switch (test)
            {
                case 0:
                    typeMethods.r8mat_copy(DIM_NUM, n, v1, ref v);
                    break;
                case 1:
                    typeMethods.r8mat_copy(DIM_NUM, n, v2, ref v);
                    break;
            }

            typeMethods.r8mat_transpose_print(DIM_NUM, n, v, "  The polygon vertices:");

            area = Geometry.polygon_area_3d(n, v, ref normal);

            Console.WriteLine("");
            Console.WriteLine("  Exact area is        " + area_exact + "");
            Console.WriteLine("  The computed area is " + area + "");
        }

    }

    public static void test0782()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST0782 tests POLYGON_AREA_3D_2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;
        int TEST_NUM = 2;

        double area;
        double area_exact;
        double[] area_exact_test = {2.4494898, 6.0};
        int n;
        int[] n_test = {4, 8};
        int test;
        double[] v;
        double[] v1 =
        {
            1.0, 0.0, 0.0,
            2.0, 1.0, 1.0,
            1.0, 2.0, 1.0,
            0.0, 1.0, 0.0
        };
        double[] v2 =
        {
            0.00000, 0.00000, 0.00000,
            2.62679, 1.26009, -0.715657,
            1.48153, 3.97300, -0.142512,
            0.605932, 3.55297, 0.0960401,
            1.36944, 1.74437, -0.286056,
            0.493842, 1.32433, -0.0475041,
            0.112090, 2.22864, 0.143544,
            -0.763505, 1.80861, 0.382097
        };

        Console.WriteLine("");
        Console.WriteLine("TEST0782");
        Console.WriteLine("  For a polygon in 3D:");
        Console.WriteLine("  POLYGON_AREA_3D_2 computes the area;");

        for (test = 0; test < TEST_NUM; test++)
        {
            area_exact = area_exact_test[test];
            n = n_test[test];

            v = new double[DIM_NUM * n];

            switch (test)
            {
                case 0:
                    typeMethods.r8mat_copy(DIM_NUM, n, v1, ref v);
                    break;
                case 1:
                    typeMethods.r8mat_copy(DIM_NUM, n, v2, ref v);
                    break;
            }

            typeMethods.r8mat_transpose_print(DIM_NUM, n, v, "  The polygon vertices:");

            area = Geometry.polygon_area_3d_2(n, v);

            Console.WriteLine("");
            Console.WriteLine("  Exact area is        " + area_exact + "");
            Console.WriteLine("  The computed area is " + area + "");

        }

    }

    public static void test0784()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST0784 tests POLYGON_CENTROID_2D and POLYGON_CENTROID_2D_2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int N = 4;

        double[] centroid;
        double[] v =
        {
            1.0, 0.0,
            2.0, 1.0,
            1.0, 2.0,
            0.0, 1.0
        };

        Console.WriteLine("");
        Console.WriteLine("TEST0784");
        Console.WriteLine("  For a polygon in 2D:");
        Console.WriteLine("  POLYGON_CENTROID_2D computes the centroid.");
        Console.WriteLine("  POLYGON_CENTROID_2D_2 computes the centroid.");

        typeMethods.r8mat_transpose_print(DIM_NUM, N, v, "  The polygon vertices:");

        centroid = Geometry.polygon_centroid_2d(N, v);

        typeMethods.r8vec_print(DIM_NUM, centroid, "  POLYGON_CENTROID_2D:");

        centroid = Geometry.polygon_centroid_2d_2(N, v);

        typeMethods.r8vec_print(DIM_NUM, centroid, "  POLYGON_CENTROID_2D_2:");

    }

    public static void polygon_centroid_3d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_CENTROID_3D_TEST tests POLYGON_CENTROID_3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 May 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;
        int N = 4;

        double[] centroid;
        double[] v =
        {
            1.0, 0.0, 0.0,
            2.0, 1.0, 1.0,
            1.0, 2.0, 1.0,
            0.0, 1.0, 0.0
        };

        Console.WriteLine("");
        Console.WriteLine("POLYGON_CENTROID_3D_TEST");
        Console.WriteLine("  For a polygon in 3D:");
        Console.WriteLine("  POLYGON_CENTROID_3D computes the centroid.");

        typeMethods.r8mat_transpose_print(DIM_NUM, N, v, "  The polygon vertices:");

        centroid = Geometry.polygon_centroid_3d(N, v);

        typeMethods.r8vec_print(DIM_NUM, centroid, "  The centroid:");

    }

    public static void polygon_contains_point_2d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_CONTAINS_POINT_2D_TEST tests POLYGON_CONTAINS_POINT_2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 November 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int dim_num = 2;
        bool inside;
        int n = 5;
        double[] p;
        double[] p_test =
        {
            1.0, 1.0,
            3.0, 4.0,
            0.0, 2.0,
            0.5, -0.25
        };
        int test;
        int test_num = 4;
        double[] v =
        {
            0.0, 0.0,
            1.0, 0.0,
            2.0, 1.0,
            1.0, 2.0,
            0.0, 2.0
        };

        Console.WriteLine("");
        Console.WriteLine("POLYGON_CONTAINS_POINT_2D_TEST");
        Console.WriteLine("  POLYGON_CONTAINS_POINT_2D determines if ");
        Console.WriteLine("  a point is in a polygon.");

        typeMethods.r8mat_transpose_print(dim_num, n, v, "  The polygon vertices:");

        Console.WriteLine("");
        Console.WriteLine("          P          Inside");
        Console.WriteLine("");

        for (test = 0; test < test_num; test++)
        {
            p = p_test.Skip(+dim_num * test).ToArray();

            inside = Geometry.polygon_contains_point_2d(n, v, p);

            Console.WriteLine("  " + p[0].ToString().PadLeft(12)
                                   + "  " + p[1].ToString().PadLeft(12)
                                   + "  " + inside + "");
        }

    }

    public static void polygon_contains_point_2d_2_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_CONTAINS_POINT_2D_2_TEST tests POLYGON_CONTAINS_POINT_2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 November 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int dim_num = 2;
        bool inside;
        int n = 5;
        double[] p;
        double[] p_test =
        {
            1.0, 1.0,
            3.0, 4.0,
            0.0, 2.0,
            0.5, -0.25
        };
        int test;
        int test_num = 4;
        double[] v =
        {
            0.0, 0.0,
            1.0, 0.0,
            2.0, 1.0,
            1.0, 2.0,
            0.0, 2.0
        };

        Console.WriteLine("");
        Console.WriteLine("POLYGON_CONTAINS_POINT_2D_2_TEST");
        Console.WriteLine("  POLYGON_CONTAINS_POINT_2D determines if ");
        Console.WriteLine("  a point is in a polygon.");

        typeMethods.r8mat_transpose_print(dim_num, n, v, "  The polygon vertices:");

        Console.WriteLine("");
        Console.WriteLine("          P          Inside");
        Console.WriteLine("");

        for (test = 0; test < test_num; test++)
        {
            p = p_test.Skip(+dim_num * test).ToArray();

            string cout = "  " + p[0].ToString().PadLeft(12)
                               + "  " + p[1].ToString().PadLeft(12);

            inside = Geometry.polygon_contains_point_2d_2(n, v, p);

            Console.WriteLine(cout + "  " + inside + "");
        }
    }

    public static void polygon_contains_point_2d_3_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_CONTAINS_POINT_2D_3_TEST tests POLYGON_CONTAINS_POINT_2D_3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 November 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int dim_num = 2;
        bool inside;
        int n = 5;
        double[] p;
        double[] p_test =
        {
            1.0, 1.0,
            3.0, 4.0,
            0.0, 2.0,
            0.5, -0.25
        };
        int test;
        int test_num = 4;
        double[] v =
        {
            0.0, 0.0,
            1.0, 0.0,
            2.0, 1.0,
            1.0, 2.0,
            0.0, 2.0
        };

        Console.WriteLine("");
        Console.WriteLine("POLYGON_CONTAINS_POINT_2D_3_TEST");
        Console.WriteLine("  POLYGON_CONTAINS_POINT_2D_3 determines if ");
        Console.WriteLine("  a point is in a polygon.");

        typeMethods.r8mat_transpose_print(dim_num, n, v, "  The polygon vertices:");

        Console.WriteLine("");
        Console.WriteLine("          P          Inside");
        Console.WriteLine("");

        for (test = 0; test < test_num; test++)
        {
            p = p_test.Skip(+dim_num * test).ToArray();

            string cout = "  " + p[0].ToString().PadLeft(12)
                               + "  " + p[1].ToString().PadLeft(12);

            inside = Geometry.polygon_contains_point_2d_3(n, v, p);

            Console.WriteLine(cout + "  " + inside + "");
        }
    }

    public static void test080()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST080 tests POLYGON_DIAMETER_2D;
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

        double diameter;
        double diameter_exact = 2.0;
        double[] v =
        {
            1.0, 0.0,
            2.0, 1.0,
            1.0, 2.0,
            0.0, 1.0
        };

        Console.WriteLine("");
        Console.WriteLine("TEST080");
        Console.WriteLine("  For a polygon in 2D:");
        Console.WriteLine("  POLYGON_DIAMETER_2D computes the diameter;");

        typeMethods.r8mat_transpose_print(DIM_NUM, N, v, "  The polygon vertices:");

        diameter = Geometry.polygon_diameter_2d(N, v);

        Console.WriteLine("");
        Console.WriteLine("  Diameter ( computed ) " + diameter + "");
        Console.WriteLine("  Diameter ( exact )    " + diameter_exact + "");

    }

    public static void test0801()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST0801 tests POLYGON_EXPAND_2D;
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

        double h;
        double[] v =
        {
            1.0, 1.0,
            5.0, 1.0,
            2.0, 4.0,
            1.0, 3.0
        };
        double[] w;

        Console.WriteLine("");
        Console.WriteLine("TEST0801");
        Console.WriteLine("  For a polygon in 2D:");
        Console.WriteLine("  POLYGON_EXPAND_2D expands it by an amount H.");

        h = 0.5;

        typeMethods.r8mat_transpose_print(DIM_NUM, N, v, "  The polygon vertices:");

        Console.WriteLine("");
        Console.WriteLine("  The expansion amount H = " + h + "");

        w = Geometry.polygon_expand_2d(N, v, h);

        typeMethods.r8mat_transpose_print(DIM_NUM, N, w, "  The expanded polygon:");

    }

    public static void test0803()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST0803 tests POLYGON_INRAD_DATA_2D, POLYGON_OUTRAD_DATA_2D, and POLYGON_SIDE_DATA_2D;
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
        double area = 0;
        int n;
        double radin = 0;
        double radout = 0;
        double side = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST0803");
        Console.WriteLine("  For a REGULAR polygon in 2D:");
        Console.WriteLine("  the inradius, outradius and side are related.");
        Console.WriteLine("  POLYGON_INRAD_DATA_2D uses the inradius;");
        Console.WriteLine("  POLYGON_OUTRAD_DATA_2D uses the inradius;");
        Console.WriteLine("  POLYGON_SIDE_DATA_2D uses the inradius;");

        for (n = 3; n <= 5; n++)
        {
            Console.WriteLine("");
            Console.WriteLine("  Number of polygonal sides = " + n + "");

            side = 1.0;

            Console.WriteLine("");
            Console.WriteLine("  Assuming SIDE = " + side + "");
            Geometry.polygon_side_data_2d(n, side, ref area, ref radin, ref radout);
            Console.WriteLine("    AREA =   " + area + "");
            Console.WriteLine("    RADIN =  " + radin + "");
            Console.WriteLine("    RADOUT = " + radout + "");

            Console.WriteLine("");
            Console.WriteLine("  Assuming RADIN = " + radin + "");
            Geometry.polygon_inrad_data_2d(n, radin, ref area, ref radout, ref side);
            Console.WriteLine("    AREA =   " + area + "");
            Console.WriteLine("    RADOUT = " + radout + "");
            Console.WriteLine("    SIDE =   " + side + "");

            Console.WriteLine("");
            Console.WriteLine("  Assuming RADOUT = " + radout + "");
            Geometry.polygon_outrad_data_2d(n, radout, ref area, ref radin, ref side);
            Console.WriteLine("    AREA =   " + area + "");
            Console.WriteLine("    RADIN =  " + radin + "");
            Console.WriteLine("    SIDE =   " + side + "");
        }
    }

    public static void test0805()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST0805 tests POLYGON_DIAMETER_2D;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;
        int N = 4;

        double diameter;
        double diameter_exact = 2.0;
        double[] v =
        {
            1.0, 0.0,
            2.0, 1.0,
            1.0, 2.0,
            0.0, 1.0
        };

        Console.WriteLine("");
        Console.WriteLine("TEST0805");
        Console.WriteLine("  For a polygon in 2D:");
        Console.WriteLine("  POLYGON_DIAMETER_2D computes the diameter;");

        typeMethods.r8mat_transpose_print(DIM_NUM, N, v, "  The polygonal vertices:");

        diameter = Geometry.polygon_diameter_2d(N, v);

        Console.WriteLine("");
        Console.WriteLine("  Diameter ( computed ) " + diameter + "");
        Console.WriteLine("  Diameter ( exact )    " + diameter_exact + "");
    }

    public static void polygon_solid_angle_3d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_SOLID_ANGLE_3D_TEST tests POLYGON_SOLID_ANGLE_3D.
        //
        //  Discussion:
        //
        //    The polygonal vertices seemed to be stored by rows instead of
        //    by columns.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 May 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n = 0;
        double[] p = new double[3];
        double solid_angle;
        int test;
        int test_num = 4;
        double[] v = null;

        Console.WriteLine("");
        Console.WriteLine("POLYGON_SOLID_ANGLE_3D_TEST");
        Console.WriteLine("  POLYGON_SOLID_ANGLE_3D computes the solid angle");
        Console.WriteLine("  subtended by a planar polygon in 3D as viewed from");
        Console.WriteLine("  a point P.");

        for (test = 1; test <= test_num; test++)
        {
            switch (test)
            {
                //
                //  One eighth of sphere surface, on the unit sphere surface.
                //
                case 1:
                    n = 3;

                    v = new double[3 * n];

                    v[0 + 0 * 3] = 1.0;
                    v[0 + 1 * 3] = 0.0;
                    v[0 + 2 * 3] = 0.0;

                    v[1 + 0 * 3] = 0.0;
                    v[1 + 1 * 3] = 1.0;
                    v[1 + 2 * 3] = 0.0;

                    v[2 + 0 * 3] = 0.0;
                    v[2 + 1 * 3] = 0.0;
                    v[2 + 2 * 3] = 1.0;

                    p[0] = 0.0;
                    p[1] = 0.0;
                    p[2] = 0.0;
                    break;
                //
                //  Reverse order of vertices.
                //
                case 2:
                    n = 3;

                    v = new double[3 * n];

                    v[0 + 0 * 3] = 1.0;
                    v[0 + 1 * 3] = 0.0;
                    v[0 + 2 * 3] = 0.0;

                    v[1 + 0 * 3] = 0.0;
                    v[1 + 1 * 3] = 0.0;
                    v[1 + 2 * 3] = 1.0;

                    v[2 + 0 * 3] = 0.0;
                    v[2 + 1 * 3] = 1.0;
                    v[2 + 2 * 3] = 0.0;

                    p[0] = 0.0;
                    p[1] = 0.0;
                    p[2] = 0.0;
                    break;
                //
                //  One eighth of sphere surface, on the unit sphere surface, 
                //  translated by (1,2,3).
                //
                case 3:
                    n = 3;

                    v = new double[3 * n];

                    v[0 + 0 * 3] = 2.0;
                    v[0 + 1 * 3] = 1.0;
                    v[0 + 2 * 3] = 1.0;

                    v[1 + 0 * 3] = 2.0;
                    v[1 + 1 * 3] = 3.0;
                    v[1 + 2 * 3] = 2.0;

                    v[2 + 0 * 3] = 3.0;
                    v[2 + 1 * 3] = 3.0;
                    v[2 + 2 * 3] = 4.0;

                    p[0] = 1.0;
                    p[1] = 2.0;
                    p[2] = 3.0;
                    break;
                //
                //  One eighth of sphere surface, but on sphere of radius 2.
                //
                case 4:
                    n = 3;

                    v = new double[3 * n];

                    v[0 + 0 * 3] = 2.0;
                    v[0 + 1 * 3] = 0.0;
                    v[0 + 2 * 3] = 0.0;

                    v[1 + 0 * 3] = 0.0;
                    v[1 + 1 * 3] = 2.0;
                    v[1 + 2 * 3] = 0.0;

                    v[2 + 0 * 3] = 0.0;
                    v[2 + 1 * 3] = 0.0;
                    v[2 + 2 * 3] = 2.0;

                    p[0] = 0.0;
                    p[1] = 0.0;
                    p[2] = 0.0;
                    break;
            }

            Console.WriteLine("");
            Console.WriteLine("  TEST #" + test + "");
            Console.WriteLine("");

            typeMethods.r8vec_print(3, p, "  The viewing point P:");

            typeMethods.r8mat_transpose_print(3, n, v, "  The polygon vertices V:");

            solid_angle = Geometry.polygon_solid_angle_3d(n, v, p);

            Console.WriteLine("");
            Console.WriteLine("  Solid angle subtended: " + solid_angle + "");

        }
    }

}