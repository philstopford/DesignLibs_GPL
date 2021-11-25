using System;
using System.Globalization;
using Burkardt.Stroud;

namespace StroudTest;

using Sphere = Sphere;

public static class tests4
{
    public static void test30()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST30 tests SPHERE_UNIT_**_3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //       
        //    01 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        string name = "";

        Console.WriteLine("");
        Console.WriteLine("TEST30");
        Console.WriteLine("  For integrals on the unit sphere in 3D:");
        Console.WriteLine("  SPHERE_UNIT_07_3D uses a formula of degree 7.");
        Console.WriteLine("  SPHERE_UNIT_11_3D uses a formula of degree 11.");
        Console.WriteLine("  SPHERE_UNIT_14_3D uses a formula of degree 14.");
        Console.WriteLine("  SPHERE_UNIT_15_3D uses a formula of degree 15.");
        Console.WriteLine("");
        Console.WriteLine("  Unit sphere area = " + Sphere.sphere_unit_area_nd(3) + "");
        Console.WriteLine("");
        Console.WriteLine("    F(X)    S3S07        S3S11         S3S14 "
                          + "         S3S15");
        Console.WriteLine("");

        int num = functions.function_3d_num();

        for (i = 1; i <= num; i++)
        {
            functions.function_3d_name(i, ref name);

            double result1 = Sphere.sphere_unit_07_3d(i, functions.function_3d);
            double result2 = Sphere.sphere_unit_11_3d(i, functions.function_3d);
            double result3 = Sphere.sphere_unit_14_3d(i, functions.function_3d);
            double result4 = Sphere.sphere_unit_15_3d(i, functions.function_3d);

            Console.WriteLine("  " + name
                                   + "  " + result1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + result2.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + result3.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + result4.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    public static void test31()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST31 tests SPHERE_UNIT_**_ND.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n;
        string name = "";

        Console.WriteLine("");
        Console.WriteLine("TEST31");
        Console.WriteLine("  For integrals on the unit sphere in ND:");
        Console.WriteLine("  SPHERE_UNIT_03_ND uses a formula of degree 3;");
        Console.WriteLine("  SPHERE_UNIT_04_ND uses a formula of degree 4;");
        Console.WriteLine("  SPHERE_UNIT_05_ND uses a formula of degree 5.");
        Console.WriteLine("  SPHERE_UNIT_07_1_ND uses a formula of degree 7.");
        Console.WriteLine("  SPHERE_UNIT_07_2_ND uses a formula of degree 7.");
        Console.WriteLine("  SPHERE_UNIT_11_ND uses a formula of degree 11.");
        Console.WriteLine("");

        for (n = 3; n <= 10; n++)
        {
            Console.WriteLine("");
            Console.WriteLine("  Spatial dimension N = " + n + "");
            Console.WriteLine("  Unit sphere area = " + Sphere.sphere_unit_area_nd(n) + "");
            Console.WriteLine("");
            Console.WriteLine("    Rule:     #3            #4            #5");
            Console.WriteLine("              #7.1          #7.2          #11");
            Console.WriteLine("    Function");
            Console.WriteLine("");

            int num = functions.function_nd_num();

            int i;
            for (i = 1; i <= num; i++)
            {
                int function_nd_index = i;
                functions.function_nd_name(function_nd_index, ref name);

                double result1 = Sphere.sphere_unit_03_nd(function_nd_index, functions.function_nd, n);
                double result2 = Sphere.sphere_unit_04_nd(function_nd_index, functions.function_nd, n);
                double result3 = Sphere.sphere_unit_05_nd(function_nd_index, functions.function_nd, n);
                double result4 = Sphere.sphere_unit_07_1_nd(function_nd_index, functions.function_nd, n);
                double result5 = Sphere.sphere_unit_07_2_nd(function_nd_index, functions.function_nd, n);
                double result6 = Sphere.sphere_unit_11_nd(function_nd_index, functions.function_nd, n);

                Console.WriteLine("  " + name
                                       + result1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + result2.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + result3.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
                Console.WriteLine("  " + "       "
                                       + result4.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + result5.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + result6.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }
        }
    }

    public static void test32()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST32 tests SPHERE_05_ND, SPHERE_07_1_ND.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n;
        string name = "";

        Console.WriteLine("");
        Console.WriteLine("TEST32");
        Console.WriteLine("  For integrals on a sphere in ND:");
        Console.WriteLine("  SPHERE_05_ND uses a formula of degree 5.");
        Console.WriteLine("  SPHERE_07_1_ND uses a formula of degree 7.");
        Console.WriteLine("");
        double r = 2.0;

        for (n = 2; n <= 4; n++)
        {
            double[] center = new double[n];
            int dim;
            for (dim = 0; dim < n; dim++)
            {
                center[dim] = 1.0;
            }

            Console.WriteLine("");
            Console.WriteLine("  Spatial dimension N = " + n + "");
            Console.WriteLine("  Sphere center = ");
            string cout = "";
            for (dim = 0; dim < n; dim++)
            {
                cout += center[dim].ToString(CultureInfo.InvariantCulture).PadLeft(14);
            }

            Console.WriteLine(cout);
            Console.WriteLine("  Sphere radius = " + r + "");
            Console.WriteLine("  Sphere area = " + Sphere.sphere_area_nd(n, r) + "");
            Console.WriteLine("");
            Console.WriteLine("    Rule:     #5           #7.1");
            Console.WriteLine("    Function");
            Console.WriteLine("");

            int num = functions.function_nd_num();

            int i;
            for (i = 1; i <= num; i++)
            {
                int function_nd_index = i;
                functions.function_nd_name(function_nd_index, ref name);

                double result1 = Sphere.sphere_05_nd(function_nd_index, functions.function_nd, n, center, r);
                double result2 = Sphere.sphere_07_1_nd(function_nd_index, functions.function_nd, n, center, r);

                Console.WriteLine("  " + name
                                       + result1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + result2.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }

        }
    }

    public static void test322()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST322 tests SPHERE_CAP_AREA_3D, SPHERE_CAP_AREA_ND.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int dim_num = 3;
        int i;
        const int ntest = 12;
        const double r = 1.0;

        Console.WriteLine("");
        Console.WriteLine("TEST322");
        Console.WriteLine("  SPHERE_CAP_AREA_3D computes the volume of a");
        Console.WriteLine("    3D spherical cap, defined by a plane that cuts the");
        Console.WriteLine("    sphere to a thickness of H units.");
        Console.WriteLine("  SPHERE_CAP_AREA_ND computes the volume of an");
        Console.WriteLine("    ND spherical cap, defined by a plane that cuts the");
        Console.WriteLine("    sphere to a thickness of H units.");

        double area1 = Sphere.sphere_area_3d(r);

        Console.WriteLine("");
        Console.WriteLine("  Area of the total sphere in 3D = " + area1 + "");

        Console.WriteLine("");
        Console.WriteLine("        R           H           Cap         Cap");
        Console.WriteLine("                                area_3d     area_nd");
        Console.WriteLine("");

        for (i = 0; i <= ntest + 1; i++)
        {
            double h = 2.0 * r * i / ntest;

            area1 = Sphere.sphere_cap_area_3d(r, h);

            double area2 = Sphere.sphere_cap_area_nd(dim_num, r, h);

            Console.WriteLine("  " + r.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + h.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + area1.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + area2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    public static void test324()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST324 tests SPHERE_CAP_VOLUME_2D, SPHERE_CAP_VOLUME_ND.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int dim_num = 2;
        int i;
        const int ntest = 12;
        const double r = 1.0;

        Console.WriteLine("");
        Console.WriteLine("TEST324");
        Console.WriteLine("  SPHERE_CAP_VOLUME_2D computes the volume (area) of a");
        Console.WriteLine("    spherical cap, defined by a plane that cuts the");
        Console.WriteLine("    sphere to a thickness of H units.");
        Console.WriteLine("  SPHERE_CAP_VOLUME_ND does the same operation,");
        Console.WriteLine("    but in N dimensions.");
        Console.WriteLine("");
        Console.WriteLine("  Using a radius R = " + r + "");

        double volume1 = Sphere.sphere_volume_2d(r);

        Console.WriteLine("");
        Console.WriteLine("  Volume of the total sphere in 2D = " + volume1 + "");

        Console.WriteLine("");
        Console.WriteLine("        H           Cap        Cap");
        Console.WriteLine("                    vol_2d     vol_nd");
        Console.WriteLine("");

        for (i = 0; i <= ntest + 1; i++)
        {
            double h = 2.0 * r * i / ntest;

            volume1 = Sphere.sphere_cap_volume_2d(r, h);

            double volume2 = Sphere.sphere_cap_volume_nd(dim_num, r, h);

            Console.WriteLine("  " + h.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + volume1.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + volume2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    public static void test326()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST326 tests SPHERE_CAP_VOLUME_3D, SPHERE_CAP_VOLUME_ND.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int dim_num = 3;
        int i;
        const int ntest = 12;
        const double r = 1.0;

        Console.WriteLine("");
        Console.WriteLine("TEST326");
        Console.WriteLine("  SPHERE_CAP_VOLUME_3D computes the volume of a");
        Console.WriteLine("    spherical cap, defined by a plane that cuts the");
        Console.WriteLine("    sphere to a thickness of H units.");
        Console.WriteLine("  SPHERE_CAP_VOLUME_ND does the same operation,");
        Console.WriteLine("    but in N dimensions.");
        Console.WriteLine("");
        Console.WriteLine("  Using a radius R = " + r + "");

        double volume1 = Sphere.sphere_volume_3d(r);

        Console.WriteLine("");
        Console.WriteLine("  Volume of the total sphere in 3D = " + volume1 + "");

        Console.WriteLine("");
        Console.WriteLine("        H           Cap        Cap");
        Console.WriteLine("                    volume_3d  volume_nd");
        Console.WriteLine("");

        for (i = 0; i <= ntest + 1; i++)
        {
            double h = 2.0 * r * i / ntest;

            volume1 = Sphere.sphere_cap_volume_3d(r, h);

            double volume2 = Sphere.sphere_cap_volume_nd(dim_num, r, h);

            Console.WriteLine("  " + h.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + volume1.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + volume2.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

    }

    public static void test33()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST33 tests SPHERE_CAP_AREA_ND, SPHERE_CAP_VOLUME_ND.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int n = 12;
        int dim_num;

        Console.WriteLine("");
        Console.WriteLine("TEST33");
        Console.WriteLine("  For a sphere in ND:");
        Console.WriteLine("  SPHERE_CAP_AREA_ND computes the area");
        Console.WriteLine("    of a spherical cap.");
        Console.WriteLine("  SPHERE_CAP_VOLUME_ND computes the volume");
        Console.WriteLine("    of a spherical cap.");
        Console.WriteLine("");

        double r = 1.0;

        for (dim_num = 2; dim_num <= 5; dim_num++)
        {
            Console.WriteLine("");
            Console.WriteLine("  Spatial dimension N = " + dim_num + "");
            Console.WriteLine("  Radius =       " + r + "");
            Console.WriteLine("  Area =         " + Sphere.sphere_area_nd(dim_num, r) + "");
            double volume = Sphere.sphere_volume_nd(dim_num, r);
            Console.WriteLine("  Volume =       " + volume + "");

            Console.WriteLine("");
            Console.WriteLine("                 Sphere         Sphere");
            Console.WriteLine("                 cap            cap");
            Console.WriteLine("      H          area           volume");
            Console.WriteLine("");

            int i;
            for (i = 0; i <= n + 1; i++)
            {
                double h = 2 * i * r / n;
                double area = Sphere.sphere_cap_area_nd(dim_num, r, h);
                volume = Sphere.sphere_cap_volume_nd(dim_num, r, h);
                Console.WriteLine("  " + h.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + area.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + volume.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }
        }

    }

    public static void test335()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST335 tests SPHERE_SHELL_03_ND.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] center = new double[3];
        int j;
        const int N_MAX = 3;
        string name = "";

        Console.WriteLine("");
        Console.WriteLine("TEST335");
        Console.WriteLine("  For integrals inside a spherical shell in ND:");
        Console.WriteLine("  SPHERE_SHELL_03_ND approximates the integral.");
        Console.WriteLine("");
        Console.WriteLine("  We compare these results with those computed by");
        Console.WriteLine("  from the difference of two ball integrals:");
        Console.WriteLine("");
        Console.WriteLine("  BALL_F1_ND approximates the integral;");
        Console.WriteLine("  BALL_F3_ND approximates the integral");
        Console.WriteLine("");

        for (j = 1; j <= 2; j++)
        {
            double r1;
            double r2;
            switch (j)
            {
                case 1:
                    r1 = 0.0;
                    r2 = 1.0;
                    center[0] = 0.0;
                    center[1] = 0.0;
                    center[2] = 0.0;
                    break;
                default:
                    r1 = 2.0;
                    r2 = 3.0;
                    center[0] = 1.0;
                    center[1] = -1.0;
                    center[2] = 2.0;
                    break;
            }

            int n;
            for (n = 2; n <= N_MAX; n++)
            {
                Console.WriteLine("");
                Console.WriteLine("  Spatial dimension N = " + n + "");
                Console.WriteLine("  Sphere center:");
                string cout = "";
                int i;
                for (i = 0; i < n; i++)
                {
                    cout += center[i].ToString(CultureInfo.InvariantCulture).PadLeft(10);
                }

                Console.WriteLine(cout);
                Console.WriteLine("  Inner sphere radius = " + r1 + "");
                Console.WriteLine("  Outer sphere radius = " + r2 + "");
                Console.WriteLine("  Spherical shell volume = "
                                  + Sphere.sphere_shell_volume_nd(n, r1, r2) + "");
                Console.WriteLine("");
                Console.WriteLine("");
                Console.WriteLine("    Rule:      #3       F1(R2)-F1(R1)  "
                                  + "F3(R2)-F3(R1)");
                Console.WriteLine("    F(X)");
                Console.WriteLine("");

                int num = functions.function_nd_num();

                for (i = 1; i <= num; i++)
                {
                    int function_nd_index = i;
                    functions.function_nd_name(function_nd_index, ref name);
                    cout = "  " + name;

                    double result1 = Sphere.sphere_shell_03_nd(function_nd_index, functions.function_nd, n, center, r1,
                        r2);

                    double result3 = Ball.ball_f1_nd(function_nd_index, functions.function_nd, n, center, r1);
                    double result4 = Ball.ball_f1_nd(function_nd_index, functions.function_nd, n, center, r2);

                    double result5 = Ball.ball_f3_nd(function_nd_index, functions.function_nd, n, center, r1);
                    double result6 = Ball.ball_f3_nd(function_nd_index, functions.function_nd, n, center, r2);

                    Console.WriteLine(result1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                      + (result4 - result3).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                      + (result6 - result5).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
                }
            }
        }
    }

    public static void test34()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST34 tests SPHERE_UNIT_AREA_ND, SPHERE_UNIT_AREA_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double area = 0;
        int dim_num = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST34:");
        Console.WriteLine("  SPHERE_UNIT_AREA_ND evaluates the area of the unit");
        Console.WriteLine("  sphere in N dimensions.");
        Console.WriteLine("  SPHERE_UNIT_AREA_VALUES returns some test values.");
        Console.WriteLine("");
        Console.WriteLine("     dim_num    Exact          Computed");
        Console.WriteLine("                Area           Area");
        Console.WriteLine("");

        int n_data = 0;

        for (;;)
        {
            Sphere.sphere_unit_area_values(ref n_data, ref dim_num, ref area);

            if (n_data == 0)
            {
                break;
            }

            double area2 = Sphere.sphere_unit_area_nd(dim_num);

            Console.WriteLine("  " + dim_num.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + area.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + area2.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }
    }

    public static void test345()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST345 tests SPHERE_UNIT_VOLUME_ND, SPHERE_UNIT_VOLUME_VALUES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int dim_num = 0;
        double volume = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST345:");
        Console.WriteLine("  SPHERE_UNIT_VOLUME_ND evaluates the area of the unit");
        Console.WriteLine("  sphere in N dimensions.");
        Console.WriteLine("  SPHERE_UNIT_VOLUME_VALUES returns some test values.");
        Console.WriteLine("");
        Console.WriteLine("     dim_num    Exact          Computed");
        Console.WriteLine("                Volume         Volume");
        Console.WriteLine("");

        int n_data = 0;

        for (;;)
        {
            Sphere.sphere_unit_volume_values(ref n_data, ref dim_num, ref volume);

            if (n_data == 0)
            {
                break;
            }

            double volume2 = Sphere.sphere_unit_volume_nd(dim_num);

            Console.WriteLine("  " + dim_num.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + volume.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + volume2.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }
    }

    public static void test35()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST35 tests SQUARE_UNIT_SET, RECTANGLE_SUB_2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        string name = "";
        int[] nsub = new int[2];
        double[] xval = { 1.0, 3.0 };
        double[] yval = { 2.0, 3.0 };

        Console.WriteLine("");
        Console.WriteLine("TEST35");
        Console.WriteLine("  SQUARE_UNIT_SET sets up a quadrature rule ");
        Console.WriteLine("    on a unit square.");
        Console.WriteLine("  RECTANGLE_SUB_2D applies it to subrectangles of an");
        Console.WriteLine("    arbitrary rectangle.");
        Console.WriteLine("");

        Console.WriteLine("");
        Console.WriteLine("  The corners of the rectangle are:");
        Console.WriteLine("");
        Console.WriteLine(xval[0].ToString(CultureInfo.InvariantCulture).PadLeft(14) + yval[0].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine(xval[1].ToString(CultureInfo.InvariantCulture).PadLeft(14) + yval[1].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        //
        //  Get the quadrature abscissas and weights for a unit square.
        //
        int rule = 2;
        int order = Square.square_unit_size(rule);

        double[] xtab = new double[order];
        double[] ytab = new double[order];
        double[] weight = new double[order];

        Square.square_unit_set(rule, order, xtab, ytab, weight);

        Console.WriteLine("");
        Console.WriteLine("  Using unit square integration rule number " + rule + "");
        Console.WriteLine("  Order of rule is " + order + "");
        //
        //  Set the function.
        //
        int num = functions.function_2d_num();

        for (i = 1; i <= num; i++)
        {
            int function_2d_index = i;
            functions.function_2d_name(function_2d_index, ref name);
            //
            //  Try an increasing number of subdivisions.
            //
            Console.WriteLine("");
            Console.WriteLine("    Function  Subdivisions  Integral");
            Console.WriteLine("");

            int j;
            for (j = 1; j <= 5; j++)
            {
                nsub[0] = j;
                nsub[1] = 2 * j;

                double result = Rectangle.rectangle_sub_2d(function_2d_index, functions.function_2d, xval, yval, nsub,
                    order, xtab,
                    ytab, weight);

                Console.WriteLine("  " + name
                                       + nsub[0].ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                       + nsub[1].ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                       + result.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }
        }

    }

    public static void test36()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST36 tests SQUARE_UNIT_SET and SQUARE_SUM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] center = { 2.0, 2.0 };
        int ilo;
        string name = "";
        const int rule_max = 6;
        const double r = 3.0;

        Console.WriteLine("");
        Console.WriteLine("TEST36");
        Console.WriteLine("  SQUARE_UNIT_SET sets up quadrature on the unit square;");
        Console.WriteLine("  SQUARE_SUM carries it out on an arbitrary square.");
        Console.WriteLine("");
        Console.WriteLine("  Square center:");
        Console.WriteLine("  CENTER = ( " + center[0] + ", " + center[1] + ")");
        Console.WriteLine("  Square radius is " + r + "");

        for (ilo = 1; ilo <= rule_max; ilo += 5)
        {
            int ihi = Math.Min(ilo + 4, rule_max);

            Console.WriteLine("");
            string cout = "  Rule:";
            int rule;
            for (rule = ilo; rule <= ihi; rule++)
            {
                cout += "       " + rule.ToString(CultureInfo.InvariantCulture).PadLeft(6);
            }

            Console.WriteLine(cout);
            Console.WriteLine("  Function ");
            Console.WriteLine("");

            int num = functions.function_2d_num();

            int i;
            for (i = 1; i <= num; i++)
            {
                int function_2d_index = i;
                functions.function_2d_name(function_2d_index, ref name);

                cout = "  " + name;

                for (rule = ilo; rule <= ihi; rule++)
                {
                    int order = Square.square_unit_size(rule);

                    double[] xtab = new double[order];
                    double[] ytab = new double[order];
                    double[] weight = new double[order];

                    Square.square_unit_set(rule, order, xtab, ytab, weight);

                    double result = Square.square_sum(function_2d_index, functions.function_2d, center, r, order, xtab,
                        ytab,
                        weight);

                    cout += result.ToString(CultureInfo.InvariantCulture).PadLeft(13);

                }

                Console.WriteLine(cout);
            }
        }
    }

    public static void test37()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST37 tests SQUARE_UNIT_SET and SQUARE_UNIT_SUM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int ilo;
        string name = "";
        const int rule_max = 6;

        Console.WriteLine("");
        Console.WriteLine("TEST37");
        Console.WriteLine("  SQUARE_UNIT_SET sets up quadrature on the unit square;");
        Console.WriteLine("  SQUARE_UNIT_SUM carries it out on the unit square.");
        Console.WriteLine("");

        for (ilo = 1; ilo <= rule_max; ilo += 5)
        {
            int ihi = Math.Min(ilo + 4, rule_max);

            Console.WriteLine("");
            string cout = "  Rule:   ";
            int rule;
            for (rule = ilo; rule <= ihi; rule++)
            {
                cout += rule.ToString(CultureInfo.InvariantCulture).PadLeft(6);
            }

            Console.WriteLine("");
            Console.WriteLine("  Function");
            Console.WriteLine("");

            int num = functions.function_2d_num();

            int i;
            for (i = 1; i <= num; i++)
            {
                int function_2d_index = i;
                functions.function_2d_name(function_2d_index, ref name);
                cout = "  " + name;

                for (rule = ilo; rule <= ihi; rule++)
                {
                    int order = Square.square_unit_size(rule);

                    double[] xtab = new double[order];
                    double[] ytab = new double[order];
                    double[] weight = new double[order];

                    Square.square_unit_set(rule, order, xtab, ytab, weight);

                    double result = Square.square_unit_sum(function_2d_index, functions.function_2d, order, xtab, ytab,
                        weight);

                    cout += result.ToString(CultureInfo.InvariantCulture).PadLeft(13);
                }

                Console.WriteLine(cout);
            }
        }
    }

    public static void test38()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST38 tests TETRA_07, TETRA_TPRODUCT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        string name = "";
        const int order_max = 9;
        double[] x = { 1.0, 4.0, 1.0, 1.0 };
        double[] y = { 2.0, 2.0, 3.0, 2.0 };
        double[] z = { 6.0, 6.0, 6.0, 8.0 };

        Console.WriteLine("");
        Console.WriteLine("TEST38");
        Console.WriteLine("  For integrals inside an arbitrary tetrahedron:");
        Console.WriteLine("  TETRA_07 uses a formula of degree 7;");
        Console.WriteLine("  TETRA_TPRODUCT uses a triangular product formula");
        Console.WriteLine("    of varying degree.");
        Console.WriteLine("");
        Console.WriteLine("  Tetrahedron vertices:");
        Console.WriteLine("");
        for (i = 0; i < 4; i++)
        {
            Console.WriteLine("  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + y[i].ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + z[i].ToString(CultureInfo.InvariantCulture).PadLeft(4) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  Tetrahedron unit volume = " + Tetrahedron.tetra_unit_volume() + "");
        Console.WriteLine("  Tetrahedron Volume = " + Tetrahedron.tetra_volume(x, y, z) + "");
        Console.WriteLine("");
        Console.WriteLine("");
        Console.WriteLine("  F(X)    TETRA_07");
        Console.WriteLine("          TETRA_TPRODUCT(1:4)");
        Console.WriteLine("          TETRA_TPRODUCT(5:8)");
        Console.WriteLine("          TETRA_TPRODUCT(9)");
        Console.WriteLine("");

        int num = functions.function_3d_num();

        for (i = 1; i <= num; i++)
        {
            int function_3d_index = i;
            functions.function_3d_name(function_3d_index, ref name);
            string cout = "  " + name;
            double result = Tetrahedron.tetra_07(function_3d_index, functions.function_3d, x, y, z);
            Console.WriteLine(cout + result.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

            int order_lo;
            for (order_lo = 1; order_lo <= order_max; order_lo += 4)
            {
                int order_hi = Math.Min(order_lo + 3, order_max);
                cout = "  " + "       ";
                int order;
                for (order = order_lo; order <= order_hi; order++)
                {
                    result = Tetrahedron.tetra_tproduct(function_3d_index, functions.function_3d, order, x, y, z);
                    cout += result.ToString(CultureInfo.InvariantCulture).PadLeft(16);
                }

                Console.WriteLine(cout);
            }

            Console.WriteLine("");
        }
    }

    public static void test39()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST39 tests TETRA_UNIT_SET and TETRA_UNIT_SUM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int ilo;
        string name = "";
        const int rule_max = 8;

        Console.WriteLine("");
        Console.WriteLine("TEST39");
        Console.WriteLine("  TETRA_UNIT_SET sets quadrature rules");
        Console.WriteLine("    for the unit tetrahedron;");
        Console.WriteLine("  TETRA_UNIT_SUM applies them to the unit tetrahedron.");
        Console.WriteLine("");

        for (ilo = 1; ilo <= rule_max; ilo += 5)
        {
            int ihi = Math.Min(ilo + 4, rule_max);

            Console.WriteLine("");
            string cout = "  Rule:   ";
            int rule;
            for (rule = ilo; rule <= ihi; rule++)
            {
                cout += rule.ToString(CultureInfo.InvariantCulture).PadLeft(6);
            }

            Console.WriteLine(cout);
            Console.WriteLine("Function");
            Console.WriteLine("");

            int num = functions.function_3d_num();

            int i;
            for (i = 1; i <= num; i++)
            {
                int function_3d_index = i;
                functions.function_3d_name(function_3d_index, ref name);
                cout = "  " + name;

                for (rule = ilo; rule <= ihi; rule++)
                {
                    int order = Tetrahedron.tetra_unit_size(rule);

                    double[] weight = new double[order];
                    double[] xtab = new double[order];
                    double[] ytab = new double[order];
                    double[] ztab = new double[order];

                    Tetrahedron.tetra_unit_set(rule, order, ref xtab, ref ytab, ref ztab, ref weight);

                    double result = Tetrahedron.tetra_unit_sum(function_3d_index, functions.function_3d, order, xtab, ytab,
                        ztab,
                        weight);

                    cout += result.ToString(CultureInfo.InvariantCulture).PadLeft(14);

                }

                Console.WriteLine(cout);
            }
        }
    }

}