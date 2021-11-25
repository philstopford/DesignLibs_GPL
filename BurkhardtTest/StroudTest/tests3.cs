using System;
using System.Globalization;
using Burkardt.IntegralNS;
using Burkardt.Stroud;
using Burkardt.Types;
using Burkardt.Uniform;

namespace StroudTest;

using Polygon = Burkardt.Stroud.Polygon;
using Simplex = Burkardt.Stroud.Simplex;

public static class tests3
{
    public static void test20()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST20 tests CUBE_UNIT_ND.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i_test;
        int[] k_test = { 10, 5 };
        const int max_test = 2;
        int[] n_test = { 2, 3 };
        string name = "";
        double[] qa = new double[10];
        double[] qb = new double[10];

        Console.WriteLine("");
        Console.WriteLine("TEST20");
        Console.WriteLine("  CUBE_UNIT_ND approximates integrals inside ");
        Console.WriteLine("    the unit cube in ND.");

        for (i_test = 0; i_test < max_test; i_test++)
        {
            int n = n_test[i_test];
            int k = k_test[i_test];
            Console.WriteLine("");
            Console.WriteLine("");
            Console.WriteLine("  Spatial dimension N = " + n + "");
            Console.WriteLine("  Value of K = " + k + "");
            Console.WriteLine("");
            Console.WriteLine("    F(X)    CUBE_UNIT_ND");
            Console.WriteLine("");

            int num = functions.function_nd_num();

            int i;
            for (i = 1; i <= num; i++)
            {
                int function_nd_index = i;
                functions.function_nd_name(function_nd_index, ref name);

                Burkardt.Stroud.Cube.cube_unit_nd(function_nd_index, functions.function_nd, ref qa, ref qb, n, k);

                int klo;
                string cout = "";
                int khi;
                int k2;
                for (klo = 0; klo <= k - 1; klo += 5)
                {
                    khi = Math.Min(klo + 4, k - 1);
                    cout = klo switch
                    {
                        0 => "  " + name,
                        _ => "  " + "       "
                    };

                    for (k2 = klo; k2 <= khi; k2++)
                    {
                        cout += qa[k2].ToString(CultureInfo.InvariantCulture).PadLeft(14);
                    }

                    Console.WriteLine(cout);
                }

                for (klo = 0; klo <= k - 1; klo += 5)
                {
                    khi = Math.Min(klo + 4, k - 1);
                    cout = "  " + "       ";
                    for (k2 = klo; k2 <= khi; k2++)
                    {
                        cout += qb[k2].ToString(CultureInfo.InvariantCulture).PadLeft(14);
                    }

                    Console.WriteLine(cout);
                }
            }
        }
    }

    public static void test205()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST205 tests ELLIPSE_AREA_2D, ELLIPSE_CIRCUMFERENCE_2D, ELLIPSE_ECCENTRICITY_2D.
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
        int i;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("TEST205");
        Console.WriteLine("  ELLIPSE_AREA_2D returns the area of an ellipse.");
        Console.WriteLine("  ELLIPSE_ECCENTRICITY_2D returns the");
        Console.WriteLine("    eccentricity of an ellipse.");
        Console.WriteLine("  ELLIPSE_CIRCUMFERENCE_2D returns the ");
        Console.WriteLine("    circumference of an ellipse.");
        Console.WriteLine("");
        Console.WriteLine("        R1        R2        E         Circum    Area");
        Console.WriteLine("");

        for (i = 1; i <= 5; i++)
        {
            double r1;
            double r2;
            switch (i)
            {
                case 1:
                    r1 = 25.0;
                    r2 = 20.0;
                    break;
                default:
                    r1 = UniformRNG.r8_uniform_01(ref seed);
                    r2 = UniformRNG.r8_uniform_01(ref seed);
                    break;
            }

            double e = Ellipse.ellipse_eccentricity_2d(r1, r2);
            double p = Ellipse.ellipse_circumference_2d(r1, r2);
            double area = Ellipse.ellipse_area_2d(r1, r2);

            Console.WriteLine("  " + r1.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + r2.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + e.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + p.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + area.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  (For the first example, ");
        Console.WriteLine("  the eccentricity should be 0.6,");
        Console.WriteLine("  the circumference should be about 141.8).");

    }

    public static void test207()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST207 tests the Stroud EN_R2 rules on monomials.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 April 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n;

        Console.WriteLine("");
        Console.WriteLine("TEST207");
        Console.WriteLine("  Demonstrate the use of Stroud rules for the region");
        Console.WriteLine("  EN_R2, that is, all of N-dimensional space, with the");
        Console.WriteLine("  weight function W(X) = exp ( - X1^2 - X2^2 ... -XN^2 )");
        Console.WriteLine("");
        Console.WriteLine("  We use the formulas to integrate various monomials of");
        Console.WriteLine("  the form X1^EXPON1 * X2^EXPON2 * ... XN^EXPONN");
        Console.WriteLine("  and compare to the exact integral.");
        Console.WriteLine("");
        Console.WriteLine("  The precision of each formula is known, and we only use");
        Console.WriteLine("  a formula if its precision indicates it should be able to");
        Console.WriteLine("  produce an exact result.");

        for (n = 1; n <= 7; n++)
        {
            int expon_dim = Math.Max(n, 2);
            int[] expon = new int[expon_dim];

            int i;
            for (i = 0; i < n; i++)
            {
                expon[i] = 0;
            }

            en.en_r2_test(n, expon);

            for (i = 0; i < n; i++)
            {
                expon[i] = 0;
            }

            expon[0] = 2;
            en.en_r2_test(n, expon);

            for (i = 0; i < n; i++)
            {
                expon[i] = 0;
            }

            expon[1] = 4;
            en.en_r2_test(n, expon);

            for (i = 0; i < n; i++)
            {
                expon[i] = 0;
            }

            i = 3 % n;
            expon[i] = 6;
            en.en_r2_test(n, expon);

            for (i = 0; i < n; i++)
            {
                expon[i] = 0;
            }

            expon[0] = 2;
            expon[1] = 4;
            en.en_r2_test(n, expon);

            for (i = 0; i < n; i++)
            {
                expon[i] = 0;
            }

            i = 4 % n;
            expon[i] = 8;
            en.en_r2_test(n, expon);

            for (i = 0; i < n; i++)
            {
                expon[i] = 0;
            }

            i = 5 % n;
            expon[i] = 10;
            en.en_r2_test(n, expon);

            for (i = 0; i < n; i++)
            {
                expon[i] = i + 1;
            }

            en.en_r2_test(n, expon);

            for (i = 0; i < n; i++)
            {
                expon[i] = 2;
            }

            en.en_r2_test(n, expon);

        }

    }

    public static void test2075()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST2075 tests the rules for EPN with GLG on monomials.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 January 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int TEST_NUM = 5;

        double[] alpha_test = { -0.5, 0.0, 0.5, 1.0, 2.0 };
        int n;

        Console.WriteLine("");
        Console.WriteLine("TEST2075");
        Console.WriteLine("  Demonstrate the use of quadrature rules for the region");
        Console.WriteLine("  EPN_GLG, that is, the positive half space [0,+oo)^N, with the");
        Console.WriteLine("  weight W(ALPHA;X) = product ( 1 <= I <= N ) X(I)^ALPHA exp ( -X(I) )");
        Console.WriteLine("");
        Console.WriteLine("  We use the formulas to integrate various monomials of");
        Console.WriteLine("  the form X(1)^E(1) * X(2)^E(2) * ... X(N)^E(N)");
        Console.WriteLine("  and compare to the exact integral.");
        Console.WriteLine("");
        Console.WriteLine("  The precision of each formula is known, and we only use");
        Console.WriteLine("  a formula if its precision indicates it should be able to");
        Console.WriteLine("  produce an exact result.");

        for (n = 1; n <= 6; n++)
        {
            int[] expon = new int[n];

            typeMethods.i4vec_zero(n, ref expon);
            double alpha;
            int test;
            for (test = 0; test < TEST_NUM; test++)
            {
                alpha = alpha_test[test];
                epn.epn_glg_test(n, expon, alpha);
            }

            typeMethods.i4vec_zero(n, ref expon);
            expon[n - 1] = 1;
            for (test = 0; test < TEST_NUM; test++)
            {
                alpha = alpha_test[test];
                epn.epn_glg_test(n, expon, alpha);
            }

            switch (n)
            {
                case >= 2:
                {
                    typeMethods.i4vec_zero(n, ref expon);
                    expon[0] = 1;
                    expon[1] = 1;
                    for (test = 0; test < TEST_NUM; test++)
                    {
                        alpha = alpha_test[test];
                        epn.epn_glg_test(n, expon, alpha);
                    }

                    break;
                }
            }

            typeMethods.i4vec_zero(n, ref expon);
            expon[0] = 2;
            for (test = 0; test < TEST_NUM; test++)
            {
                alpha = alpha_test[test];
                epn.epn_glg_test(n, expon, alpha);
            }
        }
    }

    public static void test208()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST208 tests the rules for EPN with Laguerre weight on monomials.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 January 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n;

        Console.WriteLine("");
        Console.WriteLine("TEST208");
        Console.WriteLine("  Demonstrate the use of quadrature rules for the region");
        Console.WriteLine("  EPN_LAG, that is, the positive half space [0,+oo)^N, with the");
        Console.WriteLine("  weight W(X) = product ( 1 <= I <= N ) exp ( -X(I) )");
        Console.WriteLine("");
        Console.WriteLine("  We use the formulas to integrate various monomials of");
        Console.WriteLine("  the form X(1)^E(1) * X(2)^E(2) * ... X(N)^E(N)");
        Console.WriteLine("  and compare to the exact integral.");
        Console.WriteLine("");
        Console.WriteLine("  The precision of each formula is known, and we only use");
        Console.WriteLine("  a formula if its precision indicates it should be able to");
        Console.WriteLine("  produce an exact result.");

        for (n = 1; n <= 6; n++)
        {
            int[] expon = new int[n];

            typeMethods.i4vec_zero(n, ref expon);
            epn.epn_lag_test(n, expon);

            typeMethods.i4vec_zero(n, ref expon);
            expon[n - 1] = 1;
            epn.epn_lag_test(n, expon);

            switch (n)
            {
                case >= 2:
                    typeMethods.i4vec_zero(n, ref expon);
                    expon[0] = 1;
                    expon[1] = 1;
                    epn.epn_lag_test(n, expon);
                    break;
            }

            typeMethods.i4vec_zero(n, ref expon);
            expon[0] = 2;
            epn.epn_lag_test(n, expon);

        }

    }

    public static void test21()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //   TEST21 tests HEXAGON_UNIT_SET and HEXAGON_SUM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] center = { 0.0, 0.0 };
        int ilo;
        string name = "";
        const int rule_max = 4;

        const double rad = 2.0;

        Console.WriteLine("");
        Console.WriteLine("TEST21");
        Console.WriteLine("  HEXAGON_UNIT_SET sets a quadrature rule for the");
        Console.WriteLine("    unit hexagon.");
        Console.WriteLine("  HEXAGON_SUM evaluates the quadrature rule");
        Console.WriteLine("    in an arbitrary hexagon.");
        Console.WriteLine("");
        Console.WriteLine("  We use a radius " + rad + "");
        Console.WriteLine("  and center:");
        Console.WriteLine("  CENTER = (" + center[0] + ", " + center[1] + ")");
        Console.WriteLine("");

        for (ilo = 1; ilo <= rule_max; ilo += 5)
        {
            int ihi = Math.Min(ilo + 4, rule_max);

            Console.WriteLine("");
            string cout = "  Rule:   ";
            int rule;
            for (rule = ilo; rule <= ihi; rule++)
            {
                cout += rule.ToString().PadLeft(6);
            }

            Console.WriteLine(cout);
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
                    int order = Hexagon.hexagon_unit_size(rule);

                    double[] xtab = new double[order];
                    double[] ytab = new double[order];
                    double[] weight = new double[order];

                    Hexagon.hexagon_unit_set(rule, order, ref xtab, ref ytab, ref weight);

                    double result = Hexagon.hexagon_sum(function_2d_index, functions.function_2d, center, rad, order, xtab,
                        ytab,
                        weight);

                    cout += result.ToString(CultureInfo.InvariantCulture).PadLeft(14);

                }

                Console.WriteLine(cout);
            }
        }
    }

    public static void test215()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST215 tests LENS_HALF_2D.
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
        double area;
        double[] center = { 0.0, 0.0 };
        int i;
        const int i_max = 8;
        int order;
        double theta1;
        double theta2;
        double value;

        const double r = 1.0;

        Console.WriteLine("");
        Console.WriteLine("TEST215");
        Console.WriteLine("  LENS_HALF_2D approximates an integral within a");
        Console.WriteLine("    circular half lens, defined by joining the endpoints");
        Console.WriteLine("    of a circular arc.");
        Console.WriteLine("");
        Console.WriteLine("  Integrate F(X,Y) = 1");
        Console.WriteLine("");
        Console.WriteLine("      R            Theta1      Theta2        "
                          + "Area        Order Integral");
        Console.WriteLine("");

        for (i = 0; i <= i_max; i++)
        {
            theta1 = 0.0;
            theta2 = i * 2.0 * Math.PI / i_max;

            area = Lens.lens_half_area_2d(r, theta1, theta2);

            Console.WriteLine("");

            for (order = 2; order <= 16; order += 2)
            {
                value = Lens.lens_half_2d(functions.f_1_2d, center, r, theta1, theta2, order);
                Console.WriteLine(r.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                  + theta1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                  + theta2.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                  + area.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                  + order.ToString().PadLeft(8)
                                  + value.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }
        }

        Console.WriteLine("");
        Console.WriteLine("");
        Console.WriteLine("  Integrate F(X,Y) = X");
        Console.WriteLine("");
        Console.WriteLine("      R            Theta1      Theta2        "
                          + "Area        Order Integral");
        Console.WriteLine("");

        for (i = 0; i <= i_max; i++)
        {
            theta1 = 0.0;
            theta2 = i * 2.0 * Math.PI / i_max;

            area = Lens.lens_half_area_2d(r, theta1, theta2);

            Console.WriteLine("");

            for (order = 2; order <= 16; order += 2)
            {
                value = Lens.lens_half_2d(functions.f_x_2d, center, r, theta1, theta2, order);
                Console.WriteLine(r.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                  + theta1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                  + theta2.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                  + area.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                  + order.ToString().PadLeft(8)
                                  + value.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }
        }

        Console.WriteLine("");
        Console.WriteLine("");
        Console.WriteLine("  Integrate F(X,Y) = R");
        Console.WriteLine("");
        Console.WriteLine("      R            Theta1      Theta2        "
                          + "Area        Order Integral");
        Console.WriteLine("");

        for (i = 0; i <= i_max; i++)
        {
            theta1 = 0.0;
            theta2 = i * 2.0 * Math.PI / i_max;

            area = Lens.lens_half_area_2d(r, theta1, theta2);

            Console.WriteLine("");

            for (order = 2; order <= 16; order += 2)
            {
                value = Lens.lens_half_2d(functions.f_r_2d, center, r, theta1, theta2, order);
                Console.WriteLine(r.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                  + theta1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                  + theta2.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                  + area.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                  + order.ToString().PadLeft(8)
                                  + value.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }
        }

    }

    public static void test22()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST22 tests OCTAHEDRON_UNIT_ND.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        const int N_MAX = 3;
        string name = "";

        Console.WriteLine("");
        Console.WriteLine("TEST22");
        Console.WriteLine("  OCTAHEDRON_UNIT_ND approximates integrals in a unit");
        Console.WriteLine("    octahedron in N dimensions.");
        Console.WriteLine("");
        Console.WriteLine("");
        Console.WriteLine("    F(X)    N = 1    N = 2   N = 3");
        Console.WriteLine("");

        int num = functions.function_nd_num();

        for (i = 1; i <= num; i++)
        {
            int function_nd_index = i;
            functions.function_nd_name(function_nd_index, ref name);
            string cout = "  " + name;

            int n;
            for (n = 1; n <= N_MAX; n++)
            {
                double result = Octahedron.octahedron_unit_nd(function_nd_index, functions.function_nd, n);
                cout += result.ToString(CultureInfo.InvariantCulture).PadLeft(14);
            }

            Console.WriteLine(cout);
        }
    }

    public static void test23()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST23 tests PARALLELIPIPED_VOLUME_ND.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n;

        Console.WriteLine("");
        Console.WriteLine("TEST23");
        Console.WriteLine("  PARALLELIPIPED_VOLUME_ND computes the volume of a");
        Console.WriteLine("    parallelipiped in N dimensions.");
        Console.WriteLine("");

        for (n = 2; n <= 4; n++)
        {
            Console.WriteLine("");
            Console.WriteLine("  Spatial dimension N = " + n + "");
            //
            //  Set the values of the parallelipiped.
            //
            double[] v = misc.setsim(n);

            Console.WriteLine("");
            Console.WriteLine("  Parallelipiped vertices:");
            Console.WriteLine("");

            int i;
            for (i = 0; i < n; i++)
            {
                string cout = "";
                int j;
                for (j = 0; j < n + 1; j++)
                {
                    cout += v[i + j * n].ToString(CultureInfo.InvariantCulture).PadLeft(14);
                }

                Console.WriteLine(cout);
            }

            double volume = Parallelipiped.parallelipiped_volume_nd(n, v);

            Console.WriteLine("");
            Console.WriteLine("  Volume is " + volume + "");

        }
    }

    public static void test24()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST24 tests POLYGON_**_2D;
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int npts = 4;
        double[] x = { 0.0, 1.0, 1.0, 0.0 };
        double[] y = { 0.0, 0.0, 1.0, 1.0 };

        Console.WriteLine("");
        Console.WriteLine("TEST24");
        Console.WriteLine("  For a polygon in 2D:");
        Console.WriteLine("  POLYGON_1_2D integrates 1");
        Console.WriteLine("  POLYGON_X_2D integrates X");
        Console.WriteLine("  POLYGON_Y_2D integrates Y");
        Console.WriteLine("  POLYGON_XX_2D integrates X*X");
        Console.WriteLine("  POLYGON_XY_2D integrates X*Y");
        Console.WriteLine("  POLYGON_YY_2D integrates Y*Y");
        Console.WriteLine("");
        Console.WriteLine("  F(X,Y)    Integral");
        Console.WriteLine("");

        double result = Polygon.polygon_1_2d(npts, x, y);
        Console.WriteLine("     1    " + result + "");

        result = Polygon.polygon_x_2d(npts, x, y);
        Console.WriteLine("     X    " + result + "");

        result = Polygon.polygon_y_2d(npts, x, y);
        Console.WriteLine("     Y    " + result + "");

        result = Polygon.polygon_xx_2d(npts, x, y);
        Console.WriteLine("   X*X    " + result + "");

        result = Polygon.polygon_xy_2d(npts, x, y);
        Console.WriteLine("   X*Y    " + result + "");

        result = Polygon.polygon_yy_2d(npts, x, y);
        Console.WriteLine("   Y*Y    " + result + "");

    }

    public static void test25()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST25 tests PYRAMID_UNIT_O**_3D, PYRAMID_VOLUME_3D.
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
        int jlo;
        string name = "";
        int[] order = { 1, 5, 6, 8, 8, 9, 13, 18, 27, 48 };

        Console.WriteLine("");
        Console.WriteLine("TEST25");
        Console.WriteLine("  For the unit pyramid, we approximate integrals with:");
        Console.WriteLine("  PYRAMID_UNIT_O01_3D, a 1 point rule.");
        Console.WriteLine("  PYRAMID_UNIT_O05_3D, a 5 point rule.");
        Console.WriteLine("  PYRAMID_UNIT_O06_3D, a 6 point rule.");
        Console.WriteLine("  PYRAMID_UNIT_O08_3D, an 8 point rule.");
        Console.WriteLine("  PYRAMID_UNIT_O08b_3D, an 8 point rule.");
        Console.WriteLine("  PYRAMID_UNIT_O09_3D, a 9 point rule.");
        Console.WriteLine("  PYRAMID_UNIT_O13_3D, a 13 point rule.");
        Console.WriteLine("  PYRAMID_UNIT_O18_3D, a 18 point rule.");
        Console.WriteLine("  PYRAMID_UNIT_O27_3D, a 27 point rule.");
        Console.WriteLine("  PYRAMID_UNIT_O48_3D, a 48 point rule.");
        Console.WriteLine("");
        Console.WriteLine("  PYRAMID_UNIT_VOLUME_3D computes the volume of a unit pyramid.");
        Console.WriteLine("");
        Console.WriteLine("  Volume = " + Pyramid.pyramid_unit_volume_3d() + "");
        Console.WriteLine("");

        int num = functions.function_3d_num();

        for (jlo = 1; jlo <= num; jlo += 5)
        {
            int jhi = Math.Min(jlo + 4, num);
            Console.WriteLine("");
            string cout = "  Order   ";
            int j;
            int function_3d_index;
            for (j = jlo; j <= jhi; j++)
            {
                function_3d_index = j;
                functions.function_3d_name(function_3d_index, ref name);
                cout += "   " + name + "    ";
            }

            Console.WriteLine(cout);
            Console.WriteLine("");
            cout = "  " + order[0].ToString().PadLeft(5);
            for (j = jlo; j <= jhi; j++)
            {
                function_3d_index = j;
                cout += Pyramid.pyramid_unit_o01_3d(function_3d_index, functions.function_3d).ToString(CultureInfo.InvariantCulture)
                    .PadLeft(14);
            }

            Console.WriteLine(cout);
            cout = order[1].ToString().PadLeft(5);
            for (j = jlo; j <= jhi; j++)
            {
                function_3d_index = j;
                cout += Pyramid.pyramid_unit_o05_3d(function_3d_index, functions.function_3d).ToString(CultureInfo.InvariantCulture)
                    .PadLeft(14);
            }

            Console.WriteLine(cout);
            cout = order[2].ToString().PadLeft(5);
            for (j = jlo; j <= jhi; j++)
            {
                function_3d_index = j;
                cout += Pyramid.pyramid_unit_o06_3d(function_3d_index, functions.function_3d).ToString(CultureInfo.InvariantCulture)
                    .PadLeft(14);
            }

            Console.WriteLine(cout);
            cout = order[3].ToString().PadLeft(5);
            for (j = jlo; j <= jhi; j++)
            {
                function_3d_index = j;
                cout += Pyramid.pyramid_unit_o08_3d(function_3d_index, functions.function_3d).ToString(CultureInfo.InvariantCulture)
                    .PadLeft(14);
            }

            Console.WriteLine(cout);
            cout = order[4].ToString().PadLeft(5);
            for (j = jlo; j <= jhi; j++)
            {
                function_3d_index = j;
                cout += Pyramid.pyramid_unit_o08b_3d(function_3d_index, functions.function_3d).ToString(CultureInfo.InvariantCulture)
                    .PadLeft(14);
            }

            Console.WriteLine(cout);
            cout = order[5].ToString().PadLeft(5);
            for (j = jlo; j <= jhi; j++)
            {
                function_3d_index = j;
                cout += Pyramid.pyramid_unit_o09_3d(function_3d_index, functions.function_3d).ToString(CultureInfo.InvariantCulture)
                    .PadLeft(14);
            }

            Console.WriteLine(cout);
            cout = order[6].ToString().PadLeft(5);
            for (j = jlo; j <= jhi; j++)
            {
                function_3d_index = j;
                cout += Pyramid.pyramid_unit_o13_3d(function_3d_index, functions.function_3d).ToString(CultureInfo.InvariantCulture)
                    .PadLeft(14);
            }

            Console.WriteLine(cout);
            cout = order[7].ToString().PadLeft(5);
            for (j = jlo; j <= jhi; j++)
            {
                function_3d_index = j;
                cout += Pyramid.pyramid_unit_o18_3d(function_3d_index, functions.function_3d).ToString(CultureInfo.InvariantCulture)
                    .PadLeft(14);
            }

            Console.WriteLine(cout);
            cout = order[8].ToString().PadLeft(5);
            for (j = jlo; j <= jhi; j++)
            {
                function_3d_index = j;
                cout += Pyramid.pyramid_unit_o27_3d(function_3d_index, functions.function_3d).ToString(CultureInfo.InvariantCulture)
                    .PadLeft(14);
            }

            Console.WriteLine(cout);
            cout = order[9].ToString().PadLeft(5);
            for (j = jlo; j <= jhi; j++)
            {
                function_3d_index = j;
                cout += Pyramid.pyramid_unit_o48_3d(function_3d_index, functions.function_3d).ToString(CultureInfo.InvariantCulture)
                    .PadLeft(14);
            }

            Console.WriteLine(cout);
        }
    }

    public static void test255()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST255 tests PYRAMID_UNIT_MONOMIAL_3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 March 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int alpha;
        const int degree_max = 4;

        Console.WriteLine("");
        Console.WriteLine("TEST255");
        Console.WriteLine("  For the unit pyramid,");
        Console.WriteLine("  PYRAMID_UNIT_MONOMIAL_3D returns the exact value of the");
        Console.WriteLine(" integral of X^ALPHA Y^BETA Z^GAMMA");
        Console.WriteLine("");
        Console.WriteLine("  Volume = " + Pyramid.pyramid_unit_volume_3d() + "");
        Console.WriteLine("");
        Console.WriteLine("     ALPHA      BETA     GAMMA      INTEGRAL");
        Console.WriteLine("");

        for (alpha = 0; alpha <= degree_max; alpha++)
        {
            int beta;
            for (beta = 0; beta <= degree_max - alpha; beta++)
            {
                int gamma;
                for (gamma = 0; gamma <= degree_max - alpha - beta; gamma++)
                {
                    double value = Pyramid.pyramid_unit_monomial_3d(alpha, beta, gamma);

                    Console.WriteLine("  " + alpha.ToString().PadLeft(8)
                                           + "  " + beta.ToString().PadLeft(8)
                                           + "  " + gamma.ToString().PadLeft(8)
                                           + "  " + value.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
                }
            }
        }
    }

    public static void test26()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST26 tests QMULT_1D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        string name = "";

        const double a = -1.0;
        const double b = 1.0;

        Console.WriteLine("");
        Console.WriteLine("TEST26");
        Console.WriteLine("  QMULT_1D approximates an integral on a");
        Console.WriteLine("    one-dimensional interval.");
        Console.WriteLine("");
        Console.WriteLine("  We use the interval:");
        Console.WriteLine("  A = " + a + "");
        Console.WriteLine("  B = " + b + "");
        Console.WriteLine("");
        Console.WriteLine("    F(X)     QMULT_1D");
        Console.WriteLine("");

        int num = functions.function_1d_num();

        for (i = 1; i <= num; i++)
        {
            int function_1d_index = i;
            functions.function_1d_name(function_1d_index, ref name);

            double result = QMult.qmult_1d(function_1d_index, functions.function_1d, a, b);
            Console.WriteLine("  " + name
                                   + "  " + result.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

    }

    public static void test27()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST27 tests SIMPLEX_ND.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 April 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n;
        string name = "";

        Console.WriteLine("");
        Console.WriteLine("TEST27");
        Console.WriteLine("  SIMPLEX_ND approximates integrals inside an");
        Console.WriteLine("    arbitrary simplex in ND.");
        Console.WriteLine("");

        for (n = 2; n <= 4; n++)
        {
            Console.WriteLine("");
            Console.WriteLine("  Spatial dimension N = " + n + "");
            //
            //  Restore values of simplex.
            //
            double[] v = misc.setsim(n);

            Console.WriteLine("");
            Console.WriteLine("  Simplex vertices:");
            Console.WriteLine("");

            int i;
            for (i = 0; i < n; i++)
            {
                string cout = "";
                int j;
                for (j = 0; j < n + 1; j++)
                {
                    cout += v[i + j * n].ToString(CultureInfo.InvariantCulture).PadLeft(4);
                }

                Console.WriteLine(cout);
            }

            Console.WriteLine("");
            Console.WriteLine("  F(X)    SIMPLEX_ND");
            Console.WriteLine("");

            int num = functions.function_nd_num();

            for (i = 1; i <= num; i++)
            {
                v = misc.setsim(n);
                int function_nd_index = i;
                functions.function_nd_name(function_nd_index, ref name);

                double result = Simplex.simplex_nd(function_nd_index, functions.function_nd, n, ref v);
                Console.WriteLine("  " + name
                                       + "  " + result.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }
        }
    }

    public static void test28()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST28 tests SIMPLEX_VOLUME_ND.
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
        int n;

        Console.WriteLine("");
        Console.WriteLine("TEST28");
        Console.WriteLine("  SIMPLEX_VOLUME_ND computes the volume of a simplex");
        Console.WriteLine("    in N dimensions.");
        Console.WriteLine("");

        for (n = 2; n <= 4; n++)
        {
            Console.WriteLine("");
            Console.WriteLine("  Spatial dimension N = " + n + "");
            //
            //  Set the values of the simplex.
            //
            double[] v = misc.setsim(n);

            Console.WriteLine("");
            Console.WriteLine("  Simplex vertices:");
            Console.WriteLine("");

            int i;
            for (i = 0; i < n; i++)
            {
                string cout = "";
                int j;
                for (j = 0; j < n + 1; j++)
                {
                    cout += v[i + j * n];
                }

                Console.WriteLine(cout);
            }

            double volume = Simplex.simplex_volume_nd(n, v);

            Console.WriteLine("");
            Console.WriteLine("  Volume is " + volume + "");

        }
    }

    public static void test29()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST29 tests SIMPLEX_UNIT_**_ND.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        string name = "";

        Console.WriteLine("");
        Console.WriteLine("TEST29");
        Console.WriteLine("  For integrals in the unit simplex in ND,");
        Console.WriteLine("  SIMPLEX_UNIT_01_ND uses a formula of degree 1.");
        Console.WriteLine("  SIMPLEX_UNIT_03_ND uses a formula of degree 3.");
        Console.WriteLine("  SIMPLEX_UNIT_05_ND uses a formula of degree 5.");
        Console.WriteLine("  SIMPLEX_UNIT_05_2_ND uses a formula of degree 5.");

        for (i = 1; i <= 6; i++)
        {
            int function_nd_index = i;
            functions.function_nd_name(function_nd_index, ref name);

            Console.WriteLine("");
            Console.WriteLine("  Check the integral of " + name + "");
            Console.WriteLine("");
            Console.WriteLine("  N     Volume         #1              #3"
                              + "              #5              #5.2");
            Console.WriteLine("");

            int n;
            for (n = 2; n <= 16; n++)
            {
                double result1 = Simplex.simplex_unit_01_nd(function_nd_index, functions.function_nd, n);
                double result2 = Simplex.simplex_unit_03_nd(function_nd_index, functions.function_nd, n);
                double result3 = Simplex.simplex_unit_05_nd(function_nd_index, functions.function_nd, n);
                double result4 = Simplex.simplex_unit_05_2_nd(function_nd_index, functions.function_nd, n);

                double volume = Simplex.simplex_unit_volume_nd(n);

                Console.WriteLine("  " + n.ToString().PadLeft(2)
                                       + "  " + volume.ToString(CultureInfo.InvariantCulture).PadLeft(13)
                                       + "  " + result1.ToString(CultureInfo.InvariantCulture).PadLeft(13)
                                       + "  " + result2.ToString(CultureInfo.InvariantCulture).PadLeft(13)
                                       + "  " + result3.ToString(CultureInfo.InvariantCulture).PadLeft(13)
                                       + "  " + result4.ToString(CultureInfo.InvariantCulture).PadLeft(13) + "");
            }
        }
    }
}