using System;
using System.Linq;
using Burkardt;
using Burkardt.SphereNS;
using Burkardt.Types;

namespace GeometryTest;

public static class SphereTest
{
    public static void test068()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST068 tests the SPHERE_DISTANCE routines.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 October 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int TEST_NUM = 6;

        double dist1;
        double dist2;
        double dist3;
        string[] name =
        {
            "Atlanta, Georgia  ",
            "North Pole        ",
            "South Pole        ",
            "Timbuktu          ",
            "San Antonio, Texas",
            "Savannah, Georgia "
        };
        int[] lat_d = {33, 90, -90, 16, 29, 32};
        int[] lat_m = {11, 0, 0, 49, 25, 5};
        int[] long_d = {82, 0, 0, 3, 98, 81};
        int[] long_m = {34, 0, 0, 0, 30, 6};
        double lat1;
        double lat2;
        double long1;
        double long2;
        double radius = 3957.0;
        int test1;
        int test2;

        Console.WriteLine("");
        Console.WriteLine("TEST068");
        Console.WriteLine("  SPHERE_DISTANCE1, SPHERE_DISTANCE2 and SPHERE_DISTANCE3");
        Console.WriteLine("  measure the distance between two points on a sphere.");
        Console.WriteLine("");
        Console.WriteLine("  All tests uses RADIUS = " + radius + "");
        Console.WriteLine("  which is the radius of the earth in miles.");
        Console.WriteLine("");

        for (test1 = 0; test1 < TEST_NUM - 1; test1++)
        {
            lat1 = Helpers.dms_to_radians(lat_d[test1], lat_m[test1], 0);
            long1 = Helpers.dms_to_radians(long_d[test1], long_m[test1], 0);

            Console.WriteLine("");
            Console.WriteLine("  Distance from " + name[test1] + "");

            for (test2 = test1 + 1; test2 < TEST_NUM; test2++)
            {
                lat2 = Helpers.dms_to_radians(lat_d[test2], lat_m[test2], 0);
                long2 = Helpers.dms_to_radians(long_d[test2], long_m[test2], 0);

                dist1 = Geometry.sphere_distance1(lat1, long1, lat2, long2, radius);
                dist2 = Geometry.sphere_distance2(lat1, long1, lat2, long2, radius);
                dist3 = Geometry.sphere_distance3(lat1, long1, lat2, long2, radius);

                Console.WriteLine("             to " + name[test2]
                                                     + "  " + dist1
                                                     + "  " + dist2
                                                     + "  " + dist3 + "");
            }
        }
    }

    public static void test0125()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST0125 tests SPHERE_CAP_VOLUME_2D and SPHERE_CAP_VOLUME_ND.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;

        double h;
        double haver_sine;
        double pi = 3.141592653589793;
        double r = 1.0;
        int test;
        int test_num = 12;
        double theta1;
        double theta2;
        double volume1;
        double volume2;

        Console.WriteLine("");
        Console.WriteLine("TEST0125");
        Console.WriteLine("  SPHERE_CAP_VOLUME_2D computes the volume (area) of a");
        Console.WriteLine("  spherical cap, defined by a plane that cuts the");
        Console.WriteLine("  sphere to a thickness of H units.");
        Console.WriteLine("  SPHERE_CAP_VOLUME_ND does the same operation,");
        Console.WriteLine("  but in N dimensions.");
        Console.WriteLine("");
        Console.WriteLine("  The two routines should get the same results");
        Console.WriteLine("  if THETA1, THETA2 and H correspond.");
        Console.WriteLine("");
        Console.WriteLine("  Using a radius R = " + r + "");

        Console.WriteLine("");
        Console.WriteLine("        Theta1      Theta2      H           Cap        Cap");
        Console.WriteLine("                                            vol_3d     vol_nd");
        Console.WriteLine("");

        for (test = 1; test <= test_num; test++)
        {
            h = 2.0 * r * test / test_num;

            haver_sine = Math.Sqrt(r * r - (r - h) * (r - h));

            if (h <= r)
            {
                theta2 = typeMethods.r8_asin(haver_sine / r);
            }
            else
            {
                theta2 = pi - typeMethods.r8_asin(haver_sine / r);
            }

            theta1 = -theta2;

            volume1 = Geometry.sphere_cap_volume_2d(r, h);

            volume2 = Geometry.sphere_cap_volume_nd(DIM_NUM, r, h);

            Console.WriteLine("  " + theta1.ToString().PadLeft(10)
                                   + "  " + theta2.ToString().PadLeft(10)
                                   + "  " + h.ToString().PadLeft(10)
                                   + "  " + volume1.ToString().PadLeft(10)
                                   + "  " + volume2.ToString().PadLeft(10) + "");
        }

    }

    public static void test0126()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST0126 tests SPHERE_CAP_VOLUME_3D and SPHERE_CAP_VOLUME_ND.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;

        double h;
        int test;
        int test_num = 12;
        double r = 1.0;
        double volume1;
        double volume2;

        Console.WriteLine("");
        Console.WriteLine("TEST0126");
        Console.WriteLine("  SPHERE_CAP_VOLUME_3D computes the volume of a");
        Console.WriteLine("  spherical cap, defined by a plane that cuts the");
        Console.WriteLine("  sphere to a thickness of H units.");
        Console.WriteLine("  SPHERE_CAP_VOLUME_ND does the same operation,");
        Console.WriteLine("  but in N dimensions.");
        Console.WriteLine("");
        Console.WriteLine("  Using a radius R = " + r + "");

        Console.WriteLine("");
        Console.WriteLine("        H           Cap        Cap");
        Console.WriteLine("               volume_3d  volume_nd");
        Console.WriteLine("");

        for (test = 0; test <= test_num; test++)
        {
            h = 2.0 * r * test / test_num;

            volume1 = Geometry.sphere_cap_volume_3d(r, h);

            volume2 = Geometry.sphere_cap_volume_nd(DIM_NUM, r, h);

            Console.WriteLine("  " + h.ToString().PadLeft(12)
                                   + "  " + volume1.ToString().PadLeft(12)
                                   + "  " + volume2.ToString().PadLeft(12) + "");
        }

    }

    public static void test0127()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST0127 tests SPHERE_CAP_AREA_3D and SPHERE_CAP_AREA_ND.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;

        double area1;
        double area2;
        double h;
        int test;
        int test_num = 12;
        double r = 1.0;

        Console.WriteLine("");
        Console.WriteLine("TEST0127");
        Console.WriteLine("  SPHERE_CAP_AREA_3D computes the volume of a");
        Console.WriteLine("  3D spherical cap, defined by a plane that cuts the");
        Console.WriteLine("  sphere to a thickness of H units.");
        Console.WriteLine("  SPHERE_CAP_AREA_ND computes the volume of an");
        Console.WriteLine("  ND spherical cap, defined by a plane that cuts the");
        Console.WriteLine("  sphere to a thickness of H units.");
        Console.WriteLine("");
        Console.WriteLine("        R           H           Cap         Cap");
        Console.WriteLine("                                area_3d     area_nd");
        Console.WriteLine("");

        for (test = 0; test <= test_num; test++)
        {
            h = 2.0 * r * test / test_num;

            area1 = Geometry.sphere_cap_area_3d(r, h);

            area2 = Geometry.sphere_cap_area_nd(DIM_NUM, r, h);

            Console.WriteLine("  " + r.ToString().PadLeft(12)
                                   + "  " + h.ToString().PadLeft(12)
                                   + "  " + area1.ToString().PadLeft(12)
                                   + "  " + area2.ToString().PadLeft(12) + "");
        }

    }

    public static void sphere_dia2imp_3d_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_DIA2IMP_3D_TEST tests SPHERE_DIA2IMP_3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;

        double[] p1 = {-1.0, -1.0, 4.0};
        double[] p2 = {5.0, 7.0, 4.0};
        double[] pc = new double[DIM_NUM];
        double r = 0;

        Console.WriteLine("");
        Console.WriteLine("SPHERE_DIA2IMP_3D_TEST");
        Console.WriteLine("  SPHERE_DIA2IMP_3D converts a sphere from");
        Console.WriteLine("  diameter to implicit form.");

        typeMethods.r8vec_print(DIM_NUM, p1, "  Point P1:");
        typeMethods.r8vec_print(DIM_NUM, p2, "  Point P2:");

        Geometry.sphere_dia2imp_3d(p1, p2, ref r, ref pc);

        Console.WriteLine("");
        Console.WriteLine("    Radius: " + r + "");

        typeMethods.r8vec_print(DIM_NUM, pc, "  The center:");

    }

    public static void test182()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST182 tests SPHERE_EXP_CONTAINS_POINT_3D, SPHERE_IMP_CONTAINS_POINT_3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;
        int TEST_NUM = 4;

        int i;
        bool inside;
        double[] p;
        double[] p_test =
        {
            1.0, 2.0, 3.0,
            7.0, 2.0, 3.0,
            1.0, 5.0, 3.0,
            2.5, 3.5, 4.5
        };
        double[] pc = {1.0, 2.0, 3.0};
        double[] p1 = {4.0, 2.0, 3.0};
        double[] p2 = {1.0, 5.0, 3.0};
        double[] p3 = {1.0, 2.0, 6.0};
        double[] p4 = {-2.0, 2.0, 3.0};
        double r = 3.0;
        int test;
        string cout = "";
        Console.WriteLine("");
        Console.WriteLine("TEST182");
        Console.WriteLine("  SPHERE_EXP_CONTAINS_POINT_3D determines if a");
        Console.WriteLine("  point is within an explicit sphere;");
        Console.WriteLine("  SPHERE_IMP_CONTAINS_POINT_3D determines if a");
        Console.WriteLine("  point is within an implicit sphere;");
        Console.WriteLine("");
        Console.WriteLine("  SPHERE_EXP_CONTAINS_POINT_3D:");
        Console.WriteLine("    Inside      P");
        Console.WriteLine("");

        for (test = 0; test < TEST_NUM; test++)
        {
            p = p_test.Skip(+test * DIM_NUM).ToArray();

            inside = Geometry.sphere_exp_contains_point_3d(p1, p2, p3, p4, p);

            cout = "  " + inside.ToString().PadLeft(1);
            for (i = 0; i < DIM_NUM; i++)
            {
                cout += "  " + p[i].ToString().PadLeft(12);
            }

            Console.WriteLine(cout);
        }

        Console.WriteLine("");
        Console.WriteLine("  SPHERE_IMP_CONTAINS_POINT_3D:");
        Console.WriteLine("    Inside      P");
        Console.WriteLine("");

        for (test = 0; test < TEST_NUM; test++)
        {
            p = p_test.Skip(+test * DIM_NUM).ToArray();

            inside = Geometry.sphere_imp_contains_point_3d(r, pc, p);

            cout = "  " + inside.ToString().PadLeft(1);
            for (i = 0; i < DIM_NUM; i++)
            {
                cout += "  " + p[i].ToString().PadLeft(12);
            }

            Console.WriteLine(cout);
        }

    }

    public static void test183()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST183 tests SPHERE_EXP_POINT_NEAR_3D and SPHERE_IMP_POINT_NEAR_3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;
        int TEST_NUM = 4;

        int i;
        double[] p;
        double[] p_test =
        {
            1.0, 2.0, 3.0,
            7.0, 2.0, 3.0,
            1.0, 5.0, 3.0,
            2.5, 3.5, 4.5
        };
        double[] p1 = {4.0, 2.0, 3.0};
        double[] p2 = {1.0, 5.0, 3.0};
        double[] p3 = {1.0, 2.0, 6.0};
        double[] p4 = {-2.0, 2.0, 3.0};
        double[] pc = {1.0, 2.0, 3.0};
        double[] pn = new double[DIM_NUM];
        double r = 3.0;
        int test;
        string cout = "";

        Console.WriteLine("");
        Console.WriteLine("TEST183");
        Console.WriteLine("  SPHERE_EXP_POINT_NEAR_3D determines if a");
        Console.WriteLine("  point is within an explicit sphere;");
        Console.WriteLine("  SPHERE_IMP_POINT_NEAR_3D determines if a");
        Console.WriteLine("  point is within an implicit sphere;");
        Console.WriteLine("");
        Console.WriteLine("  Sphere radius " + r + "");

        typeMethods.r8vec_print(DIM_NUM, pc, "  Sphere center:");

        Console.WriteLine("");
        Console.WriteLine("  SPHERE_EXP_POINT_NEAR_3D:");
        Console.WriteLine("    P          PN");
        Console.WriteLine("");

        for (test = 0; test < TEST_NUM; test++)
        {
            p = p_test.Skip(+test * DIM_NUM).ToArray();

            Geometry.sphere_exp_point_near_3d(p1, p2, p3, p4, p, ref pn);

            for (i = 0; i < DIM_NUM; i++)
            {
                cout = "  " + p[i].ToString().PadLeft(12);
            }

            for (i = 0; i < DIM_NUM; i++)
            {
                cout += "  " + pn[i].ToString().PadLeft(12);
            }

            Console.WriteLine(cout);
        }

        Console.WriteLine("");
        Console.WriteLine("  SPHERE_IMP_POINT_NEAR_3D:");
        Console.WriteLine("    P          PN");
        Console.WriteLine("");

        for (test = 0; test < TEST_NUM; test++)
        {
            p = p_test.Skip(+test * DIM_NUM).ToArray();

            Geometry.sphere_imp_point_near_3d(r, pc, p, ref pn);

            cout = "";
            for (i = 0; i < DIM_NUM; i++)
            {
                cout += "  " + p[i].ToString().PadLeft(12);
            }

            for (i = 0; i < DIM_NUM; i++)
            {
                cout += "  " + pn[i].ToString().PadLeft(12);
            }

            Console.WriteLine(cout);
        }

    }

    public static void test1835()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST1835 tests SPHERE_EXP2IMP_3D and SPHERE_IMP2EXP_3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;

        int i;
        double[] pc = {1.0, 2.0, 3.0};
        double[] p1 = {4.0, 2.0, 3.0};
        double[] p2 = {1.0, 5.0, 3.0};
        double[] p3 = {1.0, 2.0, 6.0};
        double[] p4 = {-2.0, 2.0, 3.0};
        double r = 3.0;
        string cout = "";

        Console.WriteLine("");
        Console.WriteLine("TEST1835");
        Console.WriteLine("  SPHERE_EXP2IMP_3D: explicit sphere => implicit form;");
        Console.WriteLine("  SPHERE_IMP2EXP_3D: implicit sphere => explicit form.");

        Console.WriteLine("");
        Console.WriteLine("  Initial form of explicit sphere:");
        Console.WriteLine("");
        cout = "  P1:";
        for (i = 0; i < DIM_NUM; i++)
        {
            cout += "  " + p1[i].ToString().PadLeft(12);
        }

        Console.WriteLine(cout);
        cout = "  P2:";
        for (i = 0; i < DIM_NUM; i++)
        {
            cout += "  " + p2[i].ToString().PadLeft(12);
        }

        Console.WriteLine(cout);
        cout = "  P3:";
        for (i = 0; i < DIM_NUM; i++)
        {
            cout += "  " + p3[i].ToString().PadLeft(12);
        }

        Console.WriteLine(cout);
        cout = "  P4:";
        for (i = 0; i < DIM_NUM; i++)
        {
            cout += "  " + p4[i].ToString().PadLeft(12);
        }

        Console.WriteLine(cout);

        Geometry.sphere_exp2imp_3d(p1, p2, p3, p4, ref r, ref pc);

        Console.WriteLine("");
        Console.WriteLine("  Computed form of implicit sphere:");
        Console.WriteLine("");
        Console.WriteLine("  Imputed radius = " + r + "");

        typeMethods.r8vec_print(DIM_NUM, pc, "  Imputed center:");

        Geometry.sphere_imp2exp_3d(r, pc, ref p1, ref p2, ref p3, ref p4);

        Console.WriteLine("");
        Console.WriteLine("  Computed form of explicit sphere:");
        Console.WriteLine("");
        cout = "  P1:";
        for (i = 0; i < DIM_NUM; i++)
        {
            cout += "  " + p1[i].ToString().PadLeft(12);
        }

        Console.WriteLine(cout);
        cout = "  P2:";
        for (i = 0; i < DIM_NUM; i++)
        {
            cout += "  " + p2[i].ToString().PadLeft(12);
        }

        Console.WriteLine(cout);
        cout = "  P3:";
        for (i = 0; i < DIM_NUM; i++)
        {
            cout += "  " + p3[i].ToString().PadLeft(12);
        }

        Console.WriteLine(cout);
        cout = "  P4:";
        for (i = 0; i < DIM_NUM; i++)
        {
            cout += "  " + p4[i].ToString().PadLeft(12);
        }

        Console.WriteLine(cout);

    }

    public static void test1836()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST1836 tests SPHERE_EXP2IMP_ND.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    26 July 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 3;

        int n = N;
        double[] p =
        {
            4.0, 2.0, 3.0,
            1.0, 5.0, 3.0,
            1.0, 2.0, 6.0,
            -2.0, 2.0, 3.0
        };
        double[] pc = new double[N];
        double[] pc_true = {1.0, 2.0, 3.0};
        double r = 0;
        double r_true = 3.0;

        Console.WriteLine("");
        Console.WriteLine("TEST1836");
        Console.WriteLine("  SPHERE_EXP2IMP_ND: explicit sphere => implicit form;");

        typeMethods.r8mat_transpose_print(n, n + 1, p, "  Initial form of explicit sphere:");

        Geometry.sphere_exp2imp_nd(n, p, ref r, ref pc);

        Console.WriteLine("");
        Console.WriteLine("  Computed form of implicit sphere:");
        Console.WriteLine("");
        Console.WriteLine("  Imputed radius = " + r + "");
        Console.WriteLine("  True radius =    " + r_true + "");

        typeMethods.r8vec_print(n, pc, "  Imputed center");

        typeMethods.r8vec_print(n, pc_true, "  True center");

    }

    public static void test187()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST187 tests SPHERE_IMP_GRIDFACES_3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;
        int TRI_MAX = 1000;

        int nlat = 3;
        int nlong = 4;
        int triangle_num = 0;
        int[] tri = new int[DIM_NUM * TRI_MAX];

        Console.WriteLine("");
        Console.WriteLine("TEST187");
        Console.WriteLine("  SPHERE_IMP_GRIDFACES_3D computes gridfaces");
        Console.WriteLine("  on a sphere in 3D.");
        Console.WriteLine("");
        Console.WriteLine("  Number of intermediate latitudes is " + nlat + "");
        Console.WriteLine("  Number of longitudes is " + nlong + "");

        Geometry.sphere_imp_gridfaces_3d(TRI_MAX, nlat, nlong, ref triangle_num, ref tri);

        Console.WriteLine("");
        Console.WriteLine("  The number of triangles is " + triangle_num + "");

        typeMethods.i4mat_transpose_print(DIM_NUM, triangle_num, tri, "  Triangle vertices:");

    }

    public static void test188()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST188 tests SPHERE_IMP_POINT_PROJECT_3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;
        int TEST_NUM = 4;

        int i;
        double[] p_test =
        {
            2.0, 0.0, 0.0,
            0.0, 4.0, 0.0,
            2.0, 4.0, 10.0,
            3.0, 5.0, 0.0
        };
        double[] p1;
        double[] p2 = new double[DIM_NUM];
        double[] pc = {2.0, 4.0, 0.0};
        double r = 2.0;
        int test;

        Console.WriteLine("");
        Console.WriteLine("TEST188");
        Console.WriteLine("  SPHERE_IMP_POINT_PROJECT_3D projects a 3D point");
        Console.WriteLine("  onto a sphere.");
        Console.WriteLine("");
        Console.WriteLine("        P1       projection P2");
        Console.WriteLine("");

        for (test = 0; test < TEST_NUM; test++)
        {
            p1 = p_test.Skip(+test * DIM_NUM).ToArray();

            Geometry.sphere_imp_point_project_3d(r, pc, p1, ref p2);

            Console.WriteLine("");
            string cout = "";
            for (i = 0; i < DIM_NUM; i++)
            {
                cout += "  " + p1[i].ToString().PadLeft(12);
            }

            for (i = 0; i < DIM_NUM; i++)
            {
                cout += "  " + p2[i].ToString().PadLeft(12);
            }

            Console.WriteLine(cout);
        }

    }

    public static void test189()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST189 tests SPHERE_IMP_AREA_ND and SPHERE_IMP_VOLUME_ND.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double area;
        int dim_num;
        double r = 1.0;
        double volume;

        Console.WriteLine("");
        Console.WriteLine("TEST189");
        Console.WriteLine("  For a implicit sphere in N dimensions:");
        Console.WriteLine("  SPHERE_IMP_AREA_ND computes the area;");
        Console.WriteLine("  SPHERE_IMP_VOLUME_ND computes the volume.");
        Console.WriteLine("");
        Console.WriteLine("  We use a radius of R = " + r + "");
        Console.WriteLine("");
        Console.WriteLine("       DIM_NUM        Area            Volume");
        Console.WriteLine("");

        for (dim_num = 2; dim_num <= 10; dim_num++)
        {
            area = Geometry.sphere_imp_area_nd(dim_num, r);
            volume = Geometry.sphere_imp_volume_nd(dim_num, r);
            Console.WriteLine("  " + dim_num.ToString().PadLeft(6)
                                   + "  " + area.ToString().PadLeft(14)
                                   + "  " + volume.ToString().PadLeft(14) + "");
        }

    }

    public static void test1892()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST1892 tests SPHERE01_POLYGON_AREA and SPHERE01_POLYGON_AREA_KARNEY.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 October 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double area;
        double area1;
        double area2;
        double[] area_test =
        {
            0.001,
            0.01,
            0.1,
            1.0,
            10.0,
            100.0,
            1000.0,
            1.0E+04,
            1.0E+05,
            1.0E+06,
            1.0E+07,
            1.0E+08,
            1.0E+09,
            1.0E+10,
            1.0E+11
        };
        int i;
        int i_test;
        double[] lat;
        double[] lat_test =
        {
            -80.0, -79.999999567493945456, -79.999999783746965784,
            -70.0, -69.999998632295765827, -69.999999316147849276,
            -60.0, -60.000003745612234573, -59.999999999999717257,
            -50.0, -49.999999999998054558, -49.999988155333397067,
            -40.0, -39.999999999986302375, -39.999962543873523006,
            -30.0, -30.000118446637603681, -29.999999999905752108,
            -20.0, -20.000374561081985306, -19.999999999405844361,
            -10.0, -9.999999997121604111, -9.998815532668760469,
            0.0, 0.002162530270410797, 0.004325060543902236,
            10.0, 9.993161263027051324, 10.00683830531098623,
            20.0, 19.999994058480001852, 20.037454636874350123,
            30.0, 29.931544224229581478, 30.06831450158650766,
            40.0, 39.567495485769567709, 39.78272828426239574,
            50.0, 49.301837019419162475, 50.669079015386089932,
            60.0, 55.676479026604850271, 57.644171397676860716
        };
        double[] lon;
        double[] lon_test =
        {
            0.0, 0.0, -0.000002157012112311,
            0.0, 0.0, -0.000003463148577445,
            0.0, 0.000004325061035168, 0.000008650121090885,
            0.0, 0.000021277700652005, 0.000010638847704917,
            0.0, 0.000056459655628232, 0.000028229812328744,
            0.0, 0.000078964535025114, 0.000157928881554143,
            0.0, 0.000230132212430763, 0.000460263329704454,
            0.0, 0.001388803276510279, 0.000694399107048577,
            0.0, 0.003745612305703688, -0.00000000000355722,
            0.0, 0.012027136076063851, 0.012027642262720488,
            0.0, 0.046026330173423388, 0.023018642559700887,
            0.0, 0.136676192075098179, 0.136864622375766145,
            0.0, 0.0, 0.487409388137601056,
            0.0, 1.816527281102587718, 1.868836724548604625,
            0.0, 0.0, 7.01051280511881837
        };
        int n;
        int[] n_test =
        {
            3, 3, 3, 3, 3,
            3, 3, 3, 3, 3,
            3, 3, 3, 3, 3
        };
        double pi = 3.141592653589793;
        double r;
        int test;
        int test_num = 15;

        Console.WriteLine("");
        Console.WriteLine("TEST1892");
        Console.WriteLine("  For a polygon on the surface of a unit sphere in 3D,");
        Console.WriteLine("  SPHERE01_POLYGON_AREA and");
        Console.WriteLine("  SPHERE01_POLYGON_AREA_KARNEY compute the area.");
        Console.WriteLine("");
        Console.WriteLine("     I     N  AREA           AREA1          ERROR1" +
                          "         AREA2          ERROR2");
        Console.WriteLine("");

        i_test = 0;
        r = 20000000.0 / pi;

        for (test = 0; test < test_num; test++)
        {
            n = n_test[test];
            area = area_test[test] / r / r;
            lat = new double[n];
            lon = new double[n];
            for (i = 0; i < n; i++)
            {
                lat[i] = lat_test[i_test] * pi / 180.0;
                lon[i] = lon_test[i_test] * pi / 180.0;
                i_test += 1;
            }

            area1 = Geometry.sphere01_polygon_area(n, lat, lon);
            area2 = Geometry.sphere01_polygon_area_karney(n, lat, lon);
            Console.WriteLine("  " + test.ToString().PadLeft(4)
                                   + "  " + n.ToString().PadLeft(4)
                                   + "  " + area.ToString().PadLeft(13)
                                   + "  " + area1.ToString().PadLeft(13)
                                   + "  " + Math.Abs(area - area1).ToString().PadLeft(13)
                                   + "  " + area2.ToString().PadLeft(13)
                                   + "  " + Math.Abs(area - area2).ToString().PadLeft(13) + "");

        }

    }

    public static void test1895()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST1895 tests SPHERE_UNIT_AREA_ND and SPHERE_UNIT_AREA_VALUES.
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
        double area2 = 0;
        int dim_num = 0;
        int n_data = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST1895:");
        Console.WriteLine("  SPHERE_UNIT_AREA_ND evaluates the area of the unit");
        Console.WriteLine("  sphere in N dimensions.");
        Console.WriteLine("  SPHERE_UNIT_AREA_VALUES returns some test values.");
        Console.WriteLine("");
        Console.WriteLine("    DIM_NUM    Exact          Computed");
        Console.WriteLine("               Area           Area");
        Console.WriteLine("");

        n_data = 0;

        for (;;)
        {
            Geometry.sphere_unit_area_values(ref n_data, ref dim_num, ref area);

            if (n_data == 0)
            {
                break;
            }

            area2 = Geometry.sphere_unit_area_nd(dim_num);

            Console.WriteLine("  " + dim_num.ToString().PadLeft(6)
                                   + "  " + area.ToString().PadLeft(14)
                                   + "  " + area2.ToString().PadLeft(14) + "");
        }

    }

    public static void test190()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST190 tests SPHERE_UNIT_SAMPLE_2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 2;

        double[] average = new double[DIM_NUM];
        double dot_average;
        int i;
        int j;
        int sample_num = 1000;
        int seed = 123456789;
        double[] v;
        double[] x;
        string cout = "";

        Console.WriteLine("");
        Console.WriteLine("TEST190");
        Console.WriteLine("  For the unit sphere in 2 dimensions (the circle):");
        Console.WriteLine("  SPHERE_UNIT_SAMPLE_2D samples;");

        Console.WriteLine("");
        Console.WriteLine("  A few sample values:");
        Console.WriteLine("");

        for (j = 1; j <= 5; j++)
        {
            cout = "";
            x = Geometry.sphere_unit_sample_2d(ref seed);
            for (i = 0; i < DIM_NUM; i++)
            {
                cout += "  " + x[i].ToString().PadLeft(12);
            }

            Console.WriteLine(cout);
        }

        Console.WriteLine("");
        Console.WriteLine("  Number of sample points = " + sample_num + "");

        typeMethods.r8vec_zero(DIM_NUM, ref average);

        for (j = 1; j <= sample_num; j++)
        {
            x = Geometry.sphere_unit_sample_2d(ref seed);
            for (i = 0; i < DIM_NUM; i++)
            {
                average[i] += x[i];
            }
        }

        for (i = 0; i < DIM_NUM; i++)
        {
            average[i] /= sample_num;
        }

        Console.WriteLine("");
        Console.WriteLine("  Now average the points, which should get a value");
        Console.WriteLine("  close to zero, and closer as sample_num increases.");
        Console.WriteLine("");
        cout = "  Average:        ";
        for (i = 0; i < DIM_NUM; i++)
        {
            cout += "  " + average[i].ToString().PadLeft(12);
        }

        Console.WriteLine(cout);

        Console.WriteLine("");
        Console.WriteLine("  Now choose a random direction, sample the same");
        Console.WriteLine("  number of points, and compute the dot product with");
        Console.WriteLine("  the direction.");
        Console.WriteLine("  Take the absolute value of each dot product ");
        Console.WriteLine("  and sum and average.");
        Console.WriteLine("");
        Console.WriteLine("  We expect a value near 2 / PI = 0.6366...");

        for (j = 1; j <= 5; j++)
        {
            v = Geometry.sphere_unit_sample_2d(ref seed);

            dot_average = 0.0;

            for (i = 1; i <= sample_num; i++)
            {
                x = Geometry.sphere_unit_sample_2d(ref seed);
                dot_average += Math.Abs(typeMethods.r8vec_dot_product(DIM_NUM, x, v));
            }

            dot_average /= sample_num;

            Console.WriteLine("");
            cout = "  V:                ";
            for (i = 0; i < DIM_NUM; i++)
            {
                cout += "  " + v[i].ToString().PadLeft(12);
            }

            Console.WriteLine(cout);

            Console.WriteLine("  Average |(XdotV)| " + dot_average + "");
        }

    }

    public static void test191()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST191 tests SPHERE_UNIT_SAMPLE_3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;

        double[] average = new double[DIM_NUM];
        double dot_average;
        int i;
        int j;
        int k;
        int sample_num = 1000;
        int seed = 123456789;
        double[] v;
        double[] x;
        string cout = "";

        Console.WriteLine("");
        Console.WriteLine("TEST191");
        Console.WriteLine("  For the unit sphere in 3 dimensions:");
        Console.WriteLine("  SPHERE_UNIT_SAMPLE_3D samples;");

        Console.WriteLine("");
        Console.WriteLine("  A few sample values:");
        Console.WriteLine("");

        for (j = 1; j <= 5; j++)
        {
            x = Geometry.sphere_unit_sample_3d(ref seed);
            for (i = 0; i < DIM_NUM; i++)
            {
                cout += "  " + x[i].ToString().PadLeft(8);
            }
        }

        Console.WriteLine(cout);

        typeMethods.r8vec_zero(DIM_NUM, ref average);

        for (j = 1; j <= sample_num; j++)
        {
            x = Geometry.sphere_unit_sample_3d(ref seed);
            for (i = 0; i < DIM_NUM; i++)
            {
                average[i] += x[i];
            }
        }

        for (i = 0; i < DIM_NUM; i++)
        {
            average[i] /= sample_num;
        }

        Console.WriteLine("");
        Console.WriteLine("  Now average the points, which should get a value");
        Console.WriteLine("  close to zero, and closer as sample_num increases.");
        Console.WriteLine("");
        cout = "  Average:        ";
        for (i = 0; i < DIM_NUM; i++)
        {
            cout += "  " + average[i].ToString().PadLeft(8);
        }

        Console.WriteLine(cout);

        Console.WriteLine("");
        Console.WriteLine("  Now choose a random direction, sample the same");
        Console.WriteLine("  number of points, and compute the dot product with");
        Console.WriteLine("  the direction.");
        Console.WriteLine("  Take the absolute value of each dot product ");
        Console.WriteLine("  and sum and average.");

        for (k = 1; k <= 5; k++)
        {
            v = Geometry.sphere_unit_sample_3d(ref seed);

            dot_average = 0.0;

            for (j = 1; j <= sample_num; j++)
            {
                x = Geometry.sphere_unit_sample_3d(ref seed);
                dot_average += Math.Abs(typeMethods.r8vec_dot_product(DIM_NUM, x, v));
            }

            dot_average /= sample_num;

            Console.WriteLine("");
            cout = "  V:                ";
            for (i = 0; i < DIM_NUM; i++)
            {
                cout += "  " + v[i].ToString().PadLeft(8);
            }

            Console.WriteLine(cout);
            Console.WriteLine("  Average |(XdotV)| " + dot_average + "");
        }

    }

    public static void test192()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST192 tests SPHERE_UNIT_SAMPLE_3D_2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;

        double[] average = new double[DIM_NUM];
        double dot_average;
        int i;
        int j;
        int k;
        int sample_num = 1000;
        int seed = 123456789;
        double[] v;
        double[] x;
        string cout = "";

        Console.WriteLine("");
        Console.WriteLine("TEST192");
        Console.WriteLine("  For the unit sphere in 3 dimensions:");
        Console.WriteLine("  SPHERE_UNIT_SAMPLE_3D_2 samples;");
        Console.WriteLine("");
        Console.WriteLine("  Warning: SPHERE_UNIT_SAMPLE_3D_2 is NOT a good code!");
        Console.WriteLine("  I only implemented it for comparison.");

        Console.WriteLine("");
        Console.WriteLine("  A few sample values:");
        Console.WriteLine("");

        for (j = 1; j <= 5; j++)
        {
            x = Geometry.sphere_unit_sample_3d_2(ref seed);
            for (i = 0; i < DIM_NUM; i++)
            {
                cout += "  " + x[i].ToString().PadLeft(8);
            }
        }

        Console.WriteLine(cout);

        typeMethods.r8vec_zero(DIM_NUM, ref average);

        for (j = 1; j <= sample_num; j++)
        {
            x = Geometry.sphere_unit_sample_3d_2(ref seed);
            for (i = 0; i < DIM_NUM; i++)
            {
                average[i] += x[i];
            }
        }

        for (i = 0; i < DIM_NUM; i++)
        {
            average[i] /= sample_num;
        }

        Console.WriteLine("");
        Console.WriteLine("  Now average the points, which should get a value");
        Console.WriteLine("  close to zero, and closer as sample_num increases.");
        Console.WriteLine("");
        cout = "  Average:        ";
        for (i = 0; i < DIM_NUM; i++)
        {
            cout += "  " + average[i].ToString().PadLeft(8);
        }

        Console.WriteLine(cout);

        Console.WriteLine("");
        Console.WriteLine("  Now choose a random direction, sample the same");
        Console.WriteLine("  number of points, and compute the dot product with");
        Console.WriteLine("  the direction.");
        Console.WriteLine("  Take the absolute value of each dot product ");
        Console.WriteLine("  and sum and average.");

        for (k = 1; k <= 5; k++)
        {
            v = Geometry.sphere_unit_sample_3d_2(ref seed);

            dot_average = 0.0;

            for (j = 1; j <= sample_num; j++)
            {
                x = Geometry.sphere_unit_sample_3d_2(ref seed);
                dot_average += Math.Abs(typeMethods.r8vec_dot_product(DIM_NUM, x, v));
            }

            dot_average /= sample_num;

            Console.WriteLine("");
            cout = "  V:                ";
            for (i = 0; i < DIM_NUM; i++)
            {
                cout += "  " + v[i].ToString().PadLeft(8);
            }

            Console.WriteLine(cout);
            Console.WriteLine("  Average |(XdotV)| " + dot_average + "");
        }

    }

    public static void test193()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST193 tests SPHERE_UNIT_SAMPLE_ND.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;

        double[] average = new double[DIM_NUM];
        double dot_average;
        int i;
        int j;
        int k;
        int sample_num = 1000;
        int seed = 123456789;
        double[] v;
        double[] x;
        string cout = "";

        Console.WriteLine("");
        Console.WriteLine("TEST193");
        Console.WriteLine("  For the unit sphere in N dimensions:");
        Console.WriteLine("  SPHERE_UNIT_SAMPLE_ND samples;");

        Console.WriteLine("");
        Console.WriteLine("  A few sample values:");
        Console.WriteLine("");

        for (j = 1; j <= 5; j++)
        {
            x = Geometry.sphere_unit_sample_nd(DIM_NUM, ref seed);
            for (i = 0; i < DIM_NUM; i++)
            {
                cout += "  " + x[i].ToString().PadLeft(8);
            }

            Console.WriteLine(cout);
        }

        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension = " + DIM_NUM + "");
        Console.WriteLine("  Number of sample points = " + sample_num + "");

        typeMethods.r8vec_zero(DIM_NUM, ref average);

        for (j = 1; j <= sample_num; j++)
        {
            x = Geometry.sphere_unit_sample_nd(DIM_NUM, ref seed);
            for (i = 0; i < DIM_NUM; i++)
            {
                average[i] += x[i];
            }
        }

        for (i = 0; i < DIM_NUM; i++)
        {
            average[i] /= sample_num;
        }

        Console.WriteLine("");
        Console.WriteLine("  Now average the points, which should get a value");
        Console.WriteLine("  close to zero, and closer as N increases.");
        Console.WriteLine("");
        Console.WriteLine("  Average:        ");
        cout = "";
        for (i = 0; i < DIM_NUM; i++)
        {
            cout += "  " + average[i].ToString().PadLeft(8);
        }

        Console.WriteLine(cout);
        Console.WriteLine("");
        Console.WriteLine("  Now choose a random direction, sample the same");
        Console.WriteLine("  number of points, and compute the dot product with");
        Console.WriteLine("  the direction.");
        Console.WriteLine("  Take the absolute value of each dot product ");
        Console.WriteLine("  and sum and average.");

        for (k = 1; k <= 5; k++)
        {
            v = Geometry.sphere_unit_sample_nd(DIM_NUM, ref seed);

            dot_average = 0.0;

            for (j = 1; j <= sample_num; j++)
            {
                x = Geometry.sphere_unit_sample_nd(DIM_NUM, ref seed);
                dot_average += Math.Abs(typeMethods.r8vec_dot_product(DIM_NUM, x, v));
            }

            dot_average /= sample_num;

            Console.WriteLine("");
            cout = "  V:                ";
            for (i = 0; i < DIM_NUM; i++)
            {
                cout += "  " + v[i].ToString().PadLeft(8);
            }

            Console.WriteLine(cout);
            Console.WriteLine("  Average |(XdotV)| " + dot_average + "");

        }

    }

    public static void sphere_unit_sample_nd_2_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_UNIT_SAMPLE_ND_2_TEST tests SPHERE_UNIT_SAMPLE_ND_2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;

        double[] average = new double[DIM_NUM];
        double dot_average;
        int i;
        int j;
        int k;
        int sample_num = 1000;
        int seed = 123456789;
        double[] v;
        double[] x;
        string cout = "";
        Console.WriteLine("");
        Console.WriteLine("SPHERE_UNIT_SAMPLE_ND_2_TEST");
        Console.WriteLine("  For the unit sphere in N dimensions:");
        Console.WriteLine("  SPHERE_UNIT_SAMPLE_ND_2 samples;");

        Console.WriteLine("");
        Console.WriteLine("  A few sample values:");
        Console.WriteLine("");

        for (j = 1; j <= 5; j++)
        {
            x = Geometry.sphere_unit_sample_nd_2(DIM_NUM, ref seed);
            for (i = 0; i < DIM_NUM; i++)
            {
                cout += "  " + x[i].ToString().PadLeft(8);
            }

            Console.WriteLine(cout);
        }

        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension = " + DIM_NUM + "");
        Console.WriteLine("  Number of sample points = " + sample_num + "");

        typeMethods.r8vec_zero(DIM_NUM, ref average);

        for (j = 1; j <= sample_num; j++)
        {
            x = Geometry.sphere_unit_sample_nd_2(DIM_NUM, ref seed);
            for (i = 0; i < DIM_NUM; i++)
            {
                average[i] += x[i];
            }
        }

        for (i = 0; i < DIM_NUM; i++)
        {
            average[i] /= sample_num;
        }

        Console.WriteLine("");
        Console.WriteLine("  Now average the points, which should get a value");
        Console.WriteLine("  close to zero, and closer as N increases.");
        Console.WriteLine("");
        Console.WriteLine("  Average:        ");
        cout = "";
        for (i = 0; i < DIM_NUM; i++)
        {
            cout += "  " + average[i].ToString().PadLeft(8);
        }

        Console.WriteLine(cout);
        Console.WriteLine("");
        Console.WriteLine("  Now choose a random direction, sample the same");
        Console.WriteLine("  number of points, and compute the dot product with");
        Console.WriteLine("  the direction.");
        Console.WriteLine("  Take the absolute value of each dot product ");
        Console.WriteLine("  and sum and average.");

        for (k = 1; k <= 5; k++)
        {
            v = Geometry.sphere_unit_sample_nd_2(DIM_NUM, ref seed);

            dot_average = 0.0;

            for (j = 1; j <= sample_num; j++)
            {
                x = Geometry.sphere_unit_sample_nd_2(DIM_NUM, ref seed);
                dot_average += Math.Abs(typeMethods.r8vec_dot_product(DIM_NUM, x, v));
            }

            dot_average /= sample_num;

            Console.WriteLine("");
            cout = "  V:                ";
            for (i = 0; i < DIM_NUM; i++)
            {
                cout += "  " + v[i].ToString().PadLeft(8);
            }

            Console.WriteLine(cout);
            Console.WriteLine("  Average |(XdotV)| " + dot_average + "");

        }

    }

    public static void test195()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST195 tests SPHERE_UNIT_SAMPLE_ND_3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int DIM_NUM = 3;

        double[] average = new double[DIM_NUM];
        double dot_average;
        int i;
        int j;
        int k;
        int sample_num = 1000;
        int seed = 123456789;
        double[] v;
        double[] x;
        string cout = "";

        Console.WriteLine("");
        Console.WriteLine("TEST195");
        Console.WriteLine("  For the unit sphere in N dimensions:");
        Console.WriteLine("  SPHERE_UNIT_SAMPLE_ND_3 samples;");

        Console.WriteLine("");
        Console.WriteLine("  A few sample values:");
        Console.WriteLine("");

        for (j = 1; j <= 5; j++)
        {
            x = Geometry.sphere_unit_sample_nd_3(DIM_NUM, ref seed);
            for (i = 0; i < DIM_NUM; i++)
            {
                cout += "  " + x[i].ToString().PadLeft(8);
            }

            Console.WriteLine(cout);
        }

        Console.WriteLine("");
        Console.WriteLine("  Spatial dimension = " + DIM_NUM + "");
        Console.WriteLine("  Number of sample points = " + sample_num + "");

        typeMethods.r8vec_zero(DIM_NUM, ref average);

        for (j = 1; j <= sample_num; j++)
        {
            x = Geometry.sphere_unit_sample_nd_3(DIM_NUM, ref seed);
            for (i = 0; i < DIM_NUM; i++)
            {
                average[i] += x[i];
            }
        }

        for (i = 0; i < DIM_NUM; i++)
        {
            average[i] /= sample_num;
        }

        Console.WriteLine("");
        Console.WriteLine("  Now average the points, which should get a value");
        Console.WriteLine("  close to zero, and closer as N increases.");
        Console.WriteLine("");
        Console.WriteLine("  Average:        ");
        cout = "";
        for (i = 0; i < DIM_NUM; i++)
        {
            cout += "  " + average[i].ToString().PadLeft(8);
        }

        Console.WriteLine(cout);
        Console.WriteLine("");
        Console.WriteLine("  Now choose a random direction, sample the same");
        Console.WriteLine("  number of points, and compute the dot product with");
        Console.WriteLine("  the direction.");
        Console.WriteLine("  Take the absolute value of each dot product ");
        Console.WriteLine("  and sum and average.");

        for (k = 1; k <= 5; k++)
        {
            v = Geometry.sphere_unit_sample_nd_3(DIM_NUM, ref seed);

            dot_average = 0.0;

            for (j = 1; j <= sample_num; j++)
            {
                x = Geometry.sphere_unit_sample_nd_3(DIM_NUM, ref seed);
                dot_average += Math.Abs(typeMethods.r8vec_dot_product(DIM_NUM, x, v));
            }

            dot_average /= sample_num;

            Console.WriteLine("");
            cout = "  V:                ";
            for (i = 0; i < DIM_NUM; i++)
            {
                cout += "  " + v[i].ToString().PadLeft(8);
            }

            Console.WriteLine(cout);
            Console.WriteLine("  Average |(XdotV)| " + dot_average + "");

        }

    }

    public static void test1955()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST1955 tests SPHERE_UNIT_VOLUME_ND and SPHERE_UNIT_VOLUME_VALUES.
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
        int dim_num = 0;
        int n_data = 0;
        double volume = 0;
        double volume2 = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST1955:");
        Console.WriteLine("  SPHERE_UNIT_VOLUME_ND evaluates the area of the unit");
        Console.WriteLine("  sphere in N dimensions.");
        Console.WriteLine("  SPHERE_UNIT_VOLUME_VALUES returns some test values.");
        Console.WriteLine("");
        Console.WriteLine("     DIM_NUM    Exact          Computed");
        Console.WriteLine("                Volume         Volume");
        Console.WriteLine("");

        n_data = 0;

        for (;;)
        {
            Geometry.sphere_unit_volume_values(ref n_data, ref dim_num, ref volume);

            if (n_data == 0)
            {
                break;
            }

            volume2 = Geometry.sphere_unit_volume_nd(dim_num);

            Console.WriteLine("  " + dim_num.ToString().PadLeft(6)
                                   + "  " + volume.ToString().PadLeft(10)
                                   + "  " + volume2.ToString().PadLeft(10) + "");
        }
    }

    public static void test200()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST200 tests SPHERE_TRIANGLE_SIDES_TO_ANGLES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a = 0;
        double as_;
        double b = 0;
        double bs = 0;
        double c = 0;
        double cs = 0;
        double r = 10.0;

        Console.WriteLine("");
        Console.WriteLine("TEST200");
        Console.WriteLine("  SPHERE_TRIANGLE_SIDES_TO_ANGLES takes the sides of a");
        Console.WriteLine("  spherical triangle and determines the angles.");

        as_ = 121.0 + 15.4 / 60.0;
        bs = 104.0 + 54.7 / 60.0;
        cs = 65.0 + 42.5 / 60.0;

        as_ = Helpers.degrees_to_radians(as_);
        bs = Helpers.degrees_to_radians(bs);
        cs = Helpers.degrees_to_radians(cs);

        as_ = r * as_;
        bs = r * bs;
        cs = r * cs;
        //
        //  Get the spherical angles.
        //
        Geometry.sphere_triangle_sides_to_angles(r, as_, bs, cs, ref a, ref b, ref c);

        Console.WriteLine("");
        Console.WriteLine("  A       = " + a + " (radians)");
        a = Helpers.radians_to_degrees(a);
        Console.WriteLine("          = " + a + " ( degrees )");
        a = 117.0 + 58.0 / 60.0;
        Console.WriteLine("  Correct = " + a + " (degrees)");

        Console.WriteLine("");
        Console.WriteLine("  B       = " + b + " (radians)");
        b = Helpers.radians_to_degrees(b);
        Console.WriteLine("          = " + b + " ( degrees )");
        b = 93.0 + 13.8 / 60.0;
        Console.WriteLine("  Correct = " + b + " (degrees)");

        Console.WriteLine("");
        Console.WriteLine("  C       = " + c + " (radians)");
        c = Helpers.radians_to_degrees(c);
        Console.WriteLine("          = " + c + " ( degrees )");
        c = 70.0 + 20.6 / 60.0;
        Console.WriteLine("  Correct = " + c + " (degrees)");

    }

}