using System;
using Burkardt.Polygon;
using Burkardt.Types;
using Burkardt.Uniform;

namespace PolygonPropertiesTest;

using Properties = Properties;

internal class Program
{
    private static void Main(string[] args)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for POLYGON_PROPERTIES_TEST.
        //
        //  Discussion:
        //
        //    POLYGON_PROPERTIES_TEST tests the POLYGON_PROPERTIES library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 October 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("POLYGON_PROPERTIES_TEST");
        Console.WriteLine("  Test the POLYGON_PROPERTIES library.");

        polygon_angles_test();
        polygon_area_test();
        polygon_area_2_test();
        polygon_centroid_test();
        polygon_centroid_2_test();
        polygon_contains_point_test();
        polygon_contains_point_2_test();
        polygon_diameter_test();
        polygon_expand_test();
        polygon_inrad_data_test();
        polygon_integral_1_test();
        polygon_integral_x_test();
        polygon_integral_xx_test();
        polygon_integral_xy_test();
        polygon_integral_y_test();
        polygon_integral_yy_test();
        polygon_is_convex_test();
        polygon_lattice_area_test();
        polygon_outrad_data_test();
        polygon_perimeter_test();
        polygon_perimeter_quad_test();
        polygon_point_dist_test();
        polygon_point_near_test();
        polygon_sample_test();
        polygon_side_data_test();
        polygon_triangulate_test();
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("POLYGON_PROPERTIES_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void polygon_angles_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_ANGLES_TEST tests POLYGON_ANGLES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 October 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] angle;
        int i;
        int n = 6;
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
        Console.WriteLine("POLYGON_ANGLES_TEST");
        Console.WriteLine("  POLYGON_ANGLES computes the angles of a polygon.");

        Console.WriteLine("");
        Console.WriteLine("  Number of polygonal vertices = " + n + "");

        typeMethods.r8mat_transpose_print(2, n, v, "  The polygon vertices:");

        angle = Properties.polygon_angles(n, v);

        Console.WriteLine("");
        Console.WriteLine("  Polygonal angles in degrees:");
        Console.WriteLine("");

        for (i = 0; i < n; i++)
        {
            Console.WriteLine(i.ToString().PadLeft(8) + "  "
                                                      + typeMethods.r8_degrees(angle[i]).ToString().PadLeft(14) +
                                                      "");
        }
    }

    private static void polygon_area_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_AREA_TEST tests POLYGON_AREA.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 October 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double area;
        double area_exact1 = 2.0;
        double area_exact2 = 6.0;
        int n1 = 4;
        int n2 = 8;
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
        Console.WriteLine("POLYGON_AREA_TEST");
        Console.WriteLine("  POLYGON_AREA computes the area of a polygon.");

        Console.WriteLine("");
        Console.WriteLine("  Number of polygonal vertices = " + n1 + "");
        typeMethods.r8mat_transpose_print(2, n1, v1, "  The polygon vertices:");
        area = Properties.polygon_area(n1, v1);
        Console.WriteLine("");
        Console.WriteLine("  Exact area is        " + area_exact1 + "");
        Console.WriteLine("  The computed area is " + area + "");

        Console.WriteLine("");
        Console.WriteLine("  Number of polygonal vertices = " + n2 + "");
        typeMethods.r8mat_transpose_print(2, n2, v2, "  The polygon vertices:");
        area = Properties.polygon_area(n2, v2);
        Console.WriteLine("");
        Console.WriteLine("  Exact area is        " + area_exact2 + "");
        Console.WriteLine("  The computed area is " + area + "");

    }

    private static void polygon_area_2_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_AREA_2_TEST tests POLYGON_AREA_2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 October 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double area;
        double area_exact1 = 2.0;
        double area_exact2 = 6.0;
        int n1 = 4;
        int n2 = 8;
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
        Console.WriteLine("POLYGON_AREA_2_TEST");
        Console.WriteLine("  POLYGON_AREA_2 computes the area of a polygon.");

        Console.WriteLine("");
        Console.WriteLine("  Number of polygonal vertices = " + n1 + "");
        typeMethods.r8mat_transpose_print(2, n1, v1, "  The polygon vertices:");
        area = Properties.polygon_area_2(n1, v1);
        Console.WriteLine("");
        Console.WriteLine("  Exact area is        " + area_exact1 + "");
        Console.WriteLine("  The computed area is " + area + "");

        Console.WriteLine("");
        Console.WriteLine("  Number of polygonal vertices = " + n2 + "");
        typeMethods.r8mat_transpose_print(2, n2, v2, "  The polygon vertices:");
        area = Properties.polygon_area_2(n2, v2);
        Console.WriteLine("");
        Console.WriteLine("  Exact area is        " + area_exact2 + "");
        Console.WriteLine("  The computed area is " + area + "");

    }

    private static void polygon_centroid_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_CENTROID_TEST tests POLYGON_CENTROID.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 October 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] centroid;
        int n = 4;
        double[] v =
        {
            1.0, 0.0,
            2.0, 1.0,
            1.0, 2.0,
            0.0, 1.0
        };

        Console.WriteLine("");
        Console.WriteLine("POLYGON_CENTROID_TEST");
        Console.WriteLine("  POLYGON_CENTROID computes the centroid of a polygon.");

        typeMethods.r8mat_transpose_print(2, n, v, "  The polygon vertices:");

        centroid = Properties.polygon_centroid(n, v);
        typeMethods.r8vec_print(2, centroid, "  POLYGON_CENTROID:");
    }

    private static void polygon_centroid_2_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_CENTROID_2_TEST tests POLYGON_CENTROID_2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 October 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] centroid;
        int n = 4;
        double[] v =
        {
            1.0, 0.0,
            2.0, 1.0,
            1.0, 2.0,
            0.0, 1.0
        };

        Console.WriteLine("");
        Console.WriteLine("POLYGON_CENTROID_2_TEST");
        Console.WriteLine("  POLYGON_CENTROID_2 computes the centroid of a polygon.");

        typeMethods.r8mat_transpose_print(2, n, v, "  The polygon vertices:");

        centroid = Properties.polygon_centroid_2(n, v);
        typeMethods.r8vec_print(2, centroid, "  POLYGON_CENTROID_2:");
    }

    private static void polygon_contains_point_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_CONTAINS_POINT_TEST tests POLYGON_CONTAINS_POINT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 October 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        bool inside;
        int n = 5;
        double[] p = new double[2];
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
        Console.WriteLine("POLYGON_CONTAINS_POINT_TEST");
        Console.WriteLine("  POLYGON_CONTAINS_POINT determines if a point is in a polygon.");

        typeMethods.r8mat_transpose_print(2, n, v, "  The polygon vertices:");

        Console.WriteLine("");
        Console.WriteLine("          P          Inside");
        Console.WriteLine("");

        for (test = 0; test < test_num; test++)
        {
            p[0] = p_test[0 + test * 2];
            p[1] = p_test[1 + test * 2];

            inside = Properties.polygon_contains_point(n, v, p);

            Console.WriteLine(p[0].ToString().PadLeft(14) + "  "
                                                          + p[1].ToString().PadLeft(14) + "  "
                                                          + inside + "");
        }

    }

    private static void polygon_contains_point_2_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_CONTAINS_POINT_2_TEST tests POLYGON_CONTAINS_POINT_2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 October 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        bool inside;
        int n = 5;
        double[] p = new double[2];
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
        Console.WriteLine("POLYGON_CONTAINS_POINT_2_TEST");
        Console.WriteLine("  POLYGON_CONTAINS_POINT_2 determines if a point is in a polygon.");

        typeMethods.r8mat_transpose_print(2, n, v, "  The polygon vertices:");

        Console.WriteLine("");
        Console.WriteLine("          P          Inside");
        Console.WriteLine("");

        for (test = 0; test < test_num; test++)
        {
            p[0] = p_test[0 + test * 2];
            p[1] = p_test[1 + test * 2];

            inside = Properties.polygon_contains_point_2(n, v, p);

            Console.WriteLine(p[0].ToString().PadLeft(14) + "  "
                                                          + p[1].ToString().PadLeft(14) + "  "
                                                          + inside + "");
        }
    }

    private static void polygon_diameter_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_DIAMETER_TEST tests POLYGON_DIAMETER;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 October 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double diameter;
        double diameter_exact = 2.0;
        int n = 4;
        double[] v =
        {
            1.0, 0.0,
            2.0, 1.0,
            1.0, 2.0,
            0.0, 1.0
        };

        Console.WriteLine("");
        Console.WriteLine("POLYGON_DIAMETER_TEST");
        Console.WriteLine("  POLYGON_DIAMETER computes the diameter of a polygon.");

        typeMethods.r8mat_transpose_print(2, n, v, "  The polygon vertices:");

        diameter = Properties.polygon_diameter(n, v);

        Console.WriteLine("");
        Console.WriteLine("  Diameter ( computed ) " + diameter + "");
        Console.WriteLine("  Diameter ( exact )    " + diameter_exact + "");

    }

    private static void polygon_expand_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_EXPAND_TEST tests POLYGON_EXPAND;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 October 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double h;
        int n = 4;
        double[] v =
        {
            1.0, 1.0,
            5.0, 1.0,
            2.0, 4.0,
            1.0, 3.0
        };
        double[] w;

        Console.WriteLine("");
        Console.WriteLine("POLGON_EXPAND_TEST");
        Console.WriteLine("  POLYGON_EXPAND expands a polygon by an amount H.");

        h = 0.5;

        typeMethods.r8mat_transpose_print(2, n, v, "  The polygon vertices:");

        Console.WriteLine("");
        Console.WriteLine("  The expansion amount H = " + h + "");

        w = Properties.polygon_expand(n, v, h);

        typeMethods.r8mat_transpose_print(2, n, w, "  The expanded polygon:");

    }

    private static void polygon_inrad_data_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_INRAD_DATA_TEST tests POLYGON_INRAD_DATA.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 October 2015
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
        Console.WriteLine("POLYGON_INRAD_DATA_TEST");
        Console.WriteLine("  POLYGON_INRAD_DATA uses the inradius of a regular polygon");
        Console.WriteLine("  to compute the area, outradius, and side length.");

        for (n = 3; n <= 5; n++)
        {
            Console.WriteLine("");
            Console.WriteLine("  Number of polygonal sides = " + n + "");
            radin = 1.0;
            Console.WriteLine("");
            Console.WriteLine("  Assuming RADIN = " + radin + "");
            Properties.polygon_inrad_data(n, radin, ref area, ref radout, ref side);
            Console.WriteLine("    AREA =   " + area + "");
            Console.WriteLine("    RADOUT = " + radout + "");
            Console.WriteLine("    SIDE =   " + side + "");
        }
    }

    private static void polygon_integral_1_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_INTEGRAL_1_TEST tests POLYGON_INTEGRAL_1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 October 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n1 = 4;
        int n2 = 3;
        double result;
        double[] v1 =
        {
            0.0, 0.0,
            1.0, 0.0,
            1.0, 1.0,
            0.0, 1.0
        };
        double[] v2 =
        {
            1.0, 1.0,
            4.0, 3.0,
            2.0, 5.0
        };

        Console.WriteLine("");
        Console.WriteLine("POLYGON_INTEGRAL_1_TEST");
        Console.WriteLine("  POLYGON_INTEGRAL_1 integrates 1 over a polygon.");

        typeMethods.r8mat_transpose_print(2, n1, v1, "  The polygon vertices:");

        result = Properties.polygon_integral_1(n1, v1);
        Console.WriteLine("");
        Console.WriteLine("    1    " + result + "");

        typeMethods.r8mat_transpose_print(2, n2, v2, "  The polygon vertices:");

        result = Properties.polygon_integral_1(n2, v2);
        Console.WriteLine("");
        Console.WriteLine("    1    " + result + "");

    }

    private static void polygon_integral_x_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_INTEGRAL_X_TEST tests POLYGON_INTEGRAL_X.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 October 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n1 = 4;
        int n2 = 3;
        double result;
        double[] v1 =
        {
            0.0, 0.0,
            1.0, 0.0,
            1.0, 1.0,
            0.0, 1.0
        };
        double[] v2 =
        {
            1.0, 1.0,
            4.0, 3.0,
            2.0, 5.0
        };

        Console.WriteLine("");
        Console.WriteLine("POLYGON_INTEGRAL_X_TEST");
        Console.WriteLine("  POLYGON_INTEGRAL_X integrates x over a polygon.");

        typeMethods.r8mat_transpose_print(2, n1, v1, "  The polygon vertices:");

        result = Properties.polygon_integral_x(n1, v1);
        Console.WriteLine("");
        Console.WriteLine("    x    " + result + "");

        typeMethods.r8mat_transpose_print(2, n2, v2, "  The polygon vertices:");

        result = Properties.polygon_integral_x(n2, v2);
        Console.WriteLine("");
        Console.WriteLine("    x    " + result + "");

    }

    private static void polygon_integral_xx_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_INTEGRAL_XX_TEST tests POLYGON_INTEGRAL_XX.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 October 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n1 = 4;
        int n2 = 3;
        double result;
        double[] v1 =
        {
            0.0, 0.0,
            1.0, 0.0,
            1.0, 1.0,
            0.0, 1.0
        };
        double[] v2 =
        {
            1.0, 1.0,
            4.0, 3.0,
            2.0, 5.0
        };

        Console.WriteLine("");
        Console.WriteLine("POLYGON_INTEGRAL_XX_TEST");
        Console.WriteLine("  POLYGON_INTEGRAL_XX integrates x^2 over a polygon.");

        typeMethods.r8mat_transpose_print(2, n1, v1, "  The polygon vertices:");

        result = Properties.polygon_integral_xx(n1, v1);
        Console.WriteLine("");
        Console.WriteLine("    x^2  " + result + "");

        typeMethods.r8mat_transpose_print(2, n2, v2, "  The polygon vertices:");

        result = Properties.polygon_integral_xx(n2, v2);
        Console.WriteLine("");
        Console.WriteLine("    x^2  " + result + "");

    }

    private static void polygon_integral_xy_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_INTEGRAL_XY_TEST tests POLYGON_INTEGRAL_XY.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 October 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n1 = 4;
        int n2 = 3;
        double result;
        double[] v1 =
        {
            0.0, 0.0,
            1.0, 0.0,
            1.0, 1.0,
            0.0, 1.0
        };
        double[] v2 =
        {
            1.0, 1.0,
            4.0, 3.0,
            2.0, 5.0
        };

        Console.WriteLine("");
        Console.WriteLine("POLYGON_INTEGRAL_XY_TEST");
        Console.WriteLine("  POLYGON_INTEGRAL_XY integrates xy over a polygon.");

        typeMethods.r8mat_transpose_print(2, n1, v1, "  The polygon vertices:");

        result = Properties.polygon_integral_xy(n1, v1);
        Console.WriteLine("");
        Console.WriteLine("    x*y  " + result + "");

        typeMethods.r8mat_transpose_print(2, n2, v2, "  The polygon vertices:");

        result = Properties.polygon_integral_xy(n2, v2);
        Console.WriteLine("");
        Console.WriteLine("    x*y  " + result + "");

    }

    private static void polygon_integral_y_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_INTEGRAL_Y_TEST tests POLYGON_INTEGRAL_Y.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 October 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n1 = 4;
        int n2 = 3;
        double result;
        double[] v1 =
        {
            0.0, 0.0,
            1.0, 0.0,
            1.0, 1.0,
            0.0, 1.0
        };
        double[] v2 =
        {
            1.0, 1.0,
            4.0, 3.0,
            2.0, 5.0
        };

        Console.WriteLine("");
        Console.WriteLine("POLYGON_INTEGRAL_Y_TEST");
        Console.WriteLine("  POLYGON_INTEGRAL_Y integrates y over a polygon.");

        typeMethods.r8mat_transpose_print(2, n1, v1, "  The polygon vertices:");

        result = Properties.polygon_integral_y(n1, v1);
        Console.WriteLine("");
        Console.WriteLine("    y    " + result + "");

        typeMethods.r8mat_transpose_print(2, n2, v2, "  The polygon vertices:");

        result = Properties.polygon_integral_y(n2, v2);
        Console.WriteLine("");
        Console.WriteLine("    y    " + result + "");

    }

    private static void polygon_integral_yy_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_INTEGRAL_YY_TEST tests POLYGON_INTEGRAL_YY.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 October 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n1 = 4;
        int n2 = 3;
        double result;
        double[] v1 =
        {
            0.0, 0.0,
            1.0, 0.0,
            1.0, 1.0,
            0.0, 1.0
        };
        double[] v2 =
        {
            1.0, 1.0,
            4.0, 3.0,
            2.0, 5.0
        };

        Console.WriteLine("");
        Console.WriteLine("POLYGON_INTEGRAL_YY_TEST");
        Console.WriteLine("  POLYGON_INTEGRAL_YY integrates y^2 over a polygon.");

        typeMethods.r8mat_transpose_print(2, n1, v1, "  The polygon vertices:");

        result = Properties.polygon_integral_yy(n1, v1);
        Console.WriteLine("");
        Console.WriteLine("    y^2  " + result + "");

        typeMethods.r8mat_transpose_print(2, n2, v2, "  The polygon vertices:");

        result = Properties.polygon_integral_yy(n2, v2);
        Console.WriteLine("");
        Console.WriteLine("    y^2  " + result + "");
    }

    private static void polygon_is_convex_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_IS_CONVEX_TEST tests POLYGON_IS_CONVEX.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 October 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double angle = 0;
        int i = 0;
        int n = 0;
        int n01 = 1;
        int n02 = 2;
        int n03 = 3;
        int n04 = 3;
        int n05 = 3;
        int n06 = 4;
        int n07 = 5;
        int n08 = 5;
        int n09 = 6;
        int n10 = 6;
        int n11 = 8;
        const double r8_pi = 3.141592653589793;
        int result = 0;
        int test;
        int test_num = 11;
        string title = "";
        double[] v = null;
        double[] v01 =
        {
            0.0, 0.0
        };
        double[] v02 =
        {
            0.0, 0.0,
            1.0, 2.0
        };
        double[] v03 =
        {
            0.0, 0.0,
            2.0, 0.0,
            1.0, 0.0
        };
        double[] v04 =
        {
            0.0, 0.0,
            1.0, 0.0,
            0.0, 2.0
        };
        double[] v05 =
        {
            0.0, 0.0,
            0.0, 2.0,
            1.0, 0.0
        };
        double[] v06 =
        {
            1.0, 0.0,
            2.0, 0.0,
            3.0, 1.0,
            0.0, 1.0
        };
        double[] v07 =
        {
            0.0, 0.0,
            0.5, 0.5,
            1.0, 0.0,
            1.0, 1.0,
            0.0, 1.0
        };
        double[] v08;
        double[] v09;
        double[] v10 =
        {
            0.0, 0.0,
            2.0, 0.0,
            1.0, 1.0,
            0.0, 0.0,
            2.0, 0.0,
            1.0, 1.0
        };
        double[] v11 =
        {
            1.0, 0.0,
            3.0, 0.0,
            3.0, 3.0,
            0.0, 3.0,
            0.0, 1.0,
            2.0, 1.0,
            2.0, 2.0,
            1.0, 2.0
        };

        Console.WriteLine("");
        Console.WriteLine("POLYGON_IS_CONVEX_TEST");
        Console.WriteLine("  POLYGON_IS_CONVEX determines if a polygon is convex.");

        for (test = 1; test <= test_num; test++)
        {
            switch (test)
            {
                case 1:
                    n = n01;
                    v = v01;
                    title = "  A point:";
                    break;
                case 2:
                    n = n02;
                    v = v02;
                    title = "  A line:";
                    break;
                case 3:
                    n = n03;
                    v = v03;
                    title = "  A triangle:";
                    break;
                case 4:
                    n = n04;
                    v = v04;
                    title = "  A CCW triangle:";
                    break;
                case 5:
                    n = n05;
                    v = v05;
                    title = "  A CW triangle:";
                    break;
                case 6:
                    n = n06;
                    v = v06;
                    title = "  Polygon with large angle:";
                    break;
                case 7:
                    n = n07;
                    v = v07;
                    title = "  Polygon with huge angle:";
                    break;
                case 8:
                {
                    n = n08;
                    v08 = new double[2 * n];
                    for (i = 0; i < n; i++)
                    {
                        angle = i * 4.0 * r8_pi / n;
                        v08[0 + i * 2] = Math.Cos(angle);
                        v08[1 + i * 2] = Math.Sin(angle);
                    }

                    v = v08;
                    title = "  A five-pointed star:";
                    break;
                }
                case 9:
                {
                    n = n09;
                    v09 = new double[2 * n];
                    for (i = 0; i < n; i++)
                    {
                        angle = i * 2.0 * r8_pi / n;
                        v09[0 + i * 2] = Math.Cos(angle);
                        v09[1 + i * 2] = Math.Sin(angle);
                    }

                    v = v09;
                    title = "  A hexagon:";
                    break;
                }
                case 10:
                    n = n10;
                    v = v10;
                    title = "  A triangle twice:";
                    break;
                case 11:
                    n = n11;
                    v = v11;
                    title = "  Square knot:";
                    break;
            }

            typeMethods.r8mat_transpose_print(2, n, v, title);
            result = Properties.polygon_is_convex(n, v);

            switch (result)
            {
                case -1:
                    Console.WriteLine("  The polygon is not convex.");
                    break;
                case 0:
                    Console.WriteLine("  The polygon is degenerate and convex.");
                    break;
                case 1:
                    Console.WriteLine("  The polygon is convex and counterclockwise.");
                    break;
                case 2:
                    Console.WriteLine("  The polygon is convex and clockwise.");
                    break;
            }
        }
    }

    private static void polygon_lattice_area_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_LATTICE_AREA_TEST tests POLYGON_LATTICE_AREA.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 October 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double area;
        int b;
        int i;

        Console.WriteLine("");
        Console.WriteLine("POLYGON_LATTICE_AREA_TEST");
        Console.WriteLine("  POLYGON_LATTICE_AREA returns the area");
        Console.WriteLine("  of a polygon, measured in lattice points.");

        i = 5;
        b = 6;

        Console.WriteLine("");
        Console.WriteLine("  Number of interior lattice points = " + i + "");
        Console.WriteLine("  Number of boundary lattice points = " + b + "");

        area = Properties.polygon_lattice_area(i, b);

        Console.WriteLine("  Area of polygon is " + area + "");
    }

    private static void polygon_outrad_data_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_OUTRAD_DATA_TEST tests POLYGON_OUTRAD_DATA.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 October 2015
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
        Console.WriteLine("POLYGON_OUTRAD_DATA_TEST");
        Console.WriteLine("  POLYGON_OUTRAD_DATA uses the outradius of a regular polygon");
        Console.WriteLine("  to determine the area, inradius, and side length.");

        for (n = 3; n <= 5; n++)
        {
            Console.WriteLine("");
            Console.WriteLine("  Number of polygonal sides = " + n + "");
            radout = 1.0;
            Console.WriteLine("");
            Console.WriteLine("  Assuming RADOUT = " + radout + "");
            Properties.polygon_outrad_data(n, radout, ref area, ref radin, ref side);
            Console.WriteLine("    AREA =   " + area + "");
            Console.WriteLine("    RADIN =  " + radin + "");
            Console.WriteLine("    SIDE =   " + side + "");
        }
    }

    private static void polygon_perimeter_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_PERIMETER_TEST tests POLYGON_PERIMETER.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 October 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n1 = 4;
        int n2 = 3;
        double result;
        double[] v1 =
        {
            0.0, 0.0,
            1.0, 0.0,
            1.0, 1.0,
            0.0, 1.0
        };
        double[] v2 =
        {
            1.0, 1.0,
            4.0, 3.0,
            2.0, 5.0
        };

        Console.WriteLine("");
        Console.WriteLine("POLYGON_PERIMETER_TEST");
        Console.WriteLine("  POLYGON_PERIMETER computes the perimeter of a polygon");

        typeMethods.r8mat_transpose_print(2, n1, v1, "  Vertices of polygon V1:");

        result = Properties.polygon_perimeter(n1, v1);
        Console.WriteLine("");
        Console.WriteLine("    Perimeter = " + result + "");

        typeMethods.r8mat_transpose_print(2, n2, v2, "  Vertices of polygon V2:");

        result = Properties.polygon_perimeter(n2, v2);
        Console.WriteLine("");
        Console.WriteLine("    Perimeter = " + result + "");

    }

    private static void polygon_perimeter_quad_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_PERIMETER_QUAD_TEST tests POLYGON_PERIMETER_QUAD.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 October 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double hmax;
        int i;
        int n1 = 4;
        int n2 = 3;
        double[] v1 =
        {
            0.0, 0.0,
            1.0, 0.0,
            1.0, 1.0,
            0.0, 1.0
        };
        double[] v2 =
        {
            1.0, 1.0,
            4.0, 3.0,
            2.0, 5.0
        };
        double value = 0;

        Console.WriteLine("");
        Console.WriteLine("POLYGON_PERIMETER_QUAD_TEST");
        Console.WriteLine("  POLYGON_PERIMETER_QUAD estimates the integral of");
        Console.WriteLine("  a function over the perimeter of a polygon using");
        Console.WriteLine("  the composite midpoint rule over each side.");

        typeMethods.r8mat_transpose_print(2, n1, v1, "  Vertices of polygon V1:");

        hmax = 0.5;
        value = Properties.polygon_perimeter_quad(n1, v1, hmax, f1);
        Console.WriteLine("");
        Console.WriteLine("  Using HMAX = " + hmax + ", estimated integral of 1 over perimeter = " + value + "");

        Console.WriteLine("");
        hmax = 2.0;
        for (i = 1; i <= 3; i++)
        {
            hmax /= 2.0;
            value = Properties.polygon_perimeter_quad(n1, v1, hmax, fx2);
            Console.WriteLine("  Using HMAX = " + hmax + ", estimated integral of x^2 over perimeter = " + value +
                              "");
        }

        typeMethods.r8mat_transpose_print(2, n2, v2, "  Vertices of polygon V2:");

        hmax = 0.5;
        value = Properties.polygon_perimeter_quad(n2, v2, hmax, f1);
        Console.WriteLine("");
        Console.WriteLine("  Using HMAX = " + hmax + ", estimated integral of 1 over perimeter = " + value + "");

        Console.WriteLine("");
        hmax = 2.0;
        for (i = 1; i <= 3; i++)
        {
            hmax /= 2.0;
            value = Properties.polygon_perimeter_quad(n2, v2, hmax, fx2);
            Console.WriteLine("  Using HMAX = " + hmax + ", estimated integral of x^2 over perimeter = " + value +
                              "");
        }

    }

    private static double f1(double x, double y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F1 evaluates f(x,y) = 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 October 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double value = 0;

        value = 1.0;

        return value;
    }

    private static double fx2(double x, double y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FX2 evaluates f(x,y) = x^2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 October 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double value = 0;

        value = x * x;

        return value;
    }

    private static void polygon_point_dist_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_POINT_DIST_TEST tests POLYGON_POINT_DIST.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 October 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double dist;
        int n = 3;
        double[] p = new double[2];
        double[] p_test =
        {
            4.0, 5.0,
            2.0, 3.0,
            -2.0, -1.0
        };
        double[] v =
        {
            1.0, 1.0,
            4.0, 3.0,
            2.0, 5.0
        };
        int test;

        Console.WriteLine("");
        Console.WriteLine("POLYGON_POINT_DIST_TEST");
        Console.WriteLine("  POLYGON_POINT_DIST computes polygon-point distance.");

        typeMethods.r8mat_transpose_print(2, n, v, "  Vertices of polygon:");

        Console.WriteLine("");
        Console.WriteLine("       X             Y             DIST");
        Console.WriteLine("");

        for (test = 0; test < 3; test++)
        {
            p[0] = p_test[0 + test * 2];
            p[1] = p_test[1 + test * 2];
            dist = Properties.polygon_point_dist(n, v, p);
            Console.WriteLine("  " + p[0].ToString().PadLeft(14)
                                   + "  " + p[1].ToString().PadLeft(14)
                                   + "  " + dist.ToString().PadLeft(14) + "");
        }
    }

    private static void polygon_point_near_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_POINT_NEAR_TEST tests POLYGON_POINT_NEAR.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 October 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n = 3;
        double[] p = new double[2];
        double[] p_test =
        {
            4.0, 5.0,
            2.0, 3.0,
            -2.0, -1.0
        };
        double[] pn;
        double[] v =
        {
            1.0, 1.0,
            4.0, 3.0,
            2.0, 5.0
        };
        int test;

        Console.WriteLine("");
        Console.WriteLine("POLYGON_POINT_NEAR_TEST");
        Console.WriteLine("  POLYGON_POINT_NEAR computes nearest point on polygon.");

        typeMethods.r8mat_transpose_print(2, n, v, "  Vertices of polygon:");

        Console.WriteLine("");
        Console.WriteLine("       X             Y             XN             YN");
        Console.WriteLine("");

        for (test = 0; test < 3; test++)
        {
            p[0] = p_test[0 + test * 2];
            p[1] = p_test[1 + test * 2];
            pn = Properties.polygon_point_near(n, v, p);
            Console.WriteLine("  " + p[0].ToString().PadLeft(14)
                                   + "  " + p[1].ToString().PadLeft(14)
                                   + "  " + pn[0].ToString().PadLeft(14)
                                   + "  " + pn[1].ToString().PadLeft(14) + "");
        }

    }

    private static void polygon_sample_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_SAMPLE_TEST tests POLYGON_SAMPLE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 October 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n = 20;
        int nv = 6;
        int seed;
        double[] v =
        {
            0.0, 0.0,
            2.0, 0.0,
            2.0, 1.0,
            1.0, 1.0,
            1.0, 2.0,
            0.0, 1.0
        };
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("POLYGON_SAMPLE_TEST");
        Console.WriteLine("  POLYGON_SAMPLE samples a polygon.");

        seed = 123456789;

        x = Properties.polygon_sample(nv, v, n, ref seed);

        typeMethods.r8mat_transpose_print(2, n, x, "  Sample points:");

    }

    private static void polygon_side_data_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_SIDE_DATA_TEST tests POLYGON_SIDE_DATA;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 October 2015
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
        Console.WriteLine("POLYGON_SIDE_DATA_TEST");
        Console.WriteLine("  POLYGON_SIDE_DATA uses the side length of a regular polygon");
        Console.WriteLine("  to determine the area, inradius and outradius.");

        for (n = 3; n <= 5; n++)
        {
            Console.WriteLine("");
            Console.WriteLine("  Number of polygonal sides = " + n + "");
            side = 1.0;
            Console.WriteLine("");
            Console.WriteLine("  Assuming SIDE = " + side + "");
            Properties.polygon_side_data(n, side, ref area, ref radin, ref radout);
            Console.WriteLine("    AREA =   " + area + "");
            Console.WriteLine("    RADIN =  " + radin + "");
            Console.WriteLine("    RADOUT = " + radout + "");
        }
    }

    private static void polygon_triangulate_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_TRIANGULATE_TEST tests POLYGON_TRIANGULATE.
        //
        //  Discussion:
        //
        //    There are N-3 triangles in the triangulation.
        //
        //    For the first N-2 triangles, the first edge listed is always an
        //    internal diagonal.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 October 2015
        //
    {
        int j;
        int n = 10;
        int[] triangles;
        double[] x =
        {
            8.0, 7.0, 6.0, 5.0, 4.0,
            3.0, 2.0, 1.0, 0.0, 4.0
        };
        double[] y =
        {
            0.0, 10.0, 0.0, 10.0, 0.0,
            10.0, 0.0, 10.0, 0.0, -2.0
        };

        Console.WriteLine("");
        Console.WriteLine("POLYGON_TRIANGULATE_TEST");
        Console.WriteLine("  POLYGON_TRIANGULATE triangulates a polygon.");
        Console.WriteLine("  Here, we triangulate the comb_10 polygon.");

        triangles = Triangulate.polygon_triangulate(n, x, y);

        Console.WriteLine("");
        Console.WriteLine("  Triangles:");
        Console.WriteLine("");

        for (j = 0; j < n - 2; j++)
        {
            Console.WriteLine("  " + j.ToString().PadLeft(2) + ":  "
                              + "  " + triangles[0 + j * 3].ToString().PadLeft(2)
                              + "  " + triangles[1 + j * 3].ToString().PadLeft(2)
                              + "  " + triangles[2 + j * 3].ToString().PadLeft(2) + "");
        }
    }
}