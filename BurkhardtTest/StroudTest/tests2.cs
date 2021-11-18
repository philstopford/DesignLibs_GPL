using System;
using Burkardt.IntegralNS;
using Burkardt.Stroud;
using Burkardt.Types;

namespace StroudTest;

public static class tests2
{
    public static void test10()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST10 tests CIRCLE_CUM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] center =
            {
                0.0, 0.0
            }
            ;
        int i;
        int j;
        string name = "";
        int num;
        int order;
        double r = 3.0;

        Console.WriteLine("");
        Console.WriteLine("TEST10");
        Console.WriteLine("  CIRCLE_CUM approximates an integral over a circle.");
        Console.WriteLine("");
        Console.WriteLine("  We use radius R = " + r + "");
        Console.WriteLine("  and center:");
        Console.WriteLine("  CENTER = ( " + center[0] + ", " + center[1] + ").");
        Console.WriteLine("");

        Console.WriteLine("");
        Console.WriteLine("    Order:      2             4              8            16");
        Console.WriteLine("  F(X)");
        Console.WriteLine("");

        num = functions.function_2d_num();

        for (i = 1; i <= num; i++)
        {
            int function_2d_index = i;
            functions.function_2d_name(function_2d_index, ref name);

            string cout = "  " + name;

            for (j = 1; j <= 4; j++)
            {
                order = (int)Math.Pow(2, j);

                cout += Circle.circle_cum(function_2d_index, functions.function_2d, center, r, order).ToString(CultureInfo.InvariantCulture)
                    .PadLeft(14);
            }

            Console.WriteLine(cout);
        }

    }

    public static void test11()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST11 tests LENS_HALF_AREA_2D, CIRCLE_SECTOR_AREA_2D, CIRCLE_TRIANGLE_AREA_2D.
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
        double area1;
        double area2;
        double area3;
        int i;
        double pi = 3.141592653589793;
        double r;
        double theta1;
        double theta2;

        r = 1.0;

        Console.WriteLine("");
        Console.WriteLine("TEST11");
        Console.WriteLine("  LENS_HALF_AREA_2D computes the area of a");
        Console.WriteLine("    circular half lens, defined by joining the endpoints");
        Console.WriteLine("    of a circular arc.");
        Console.WriteLine("  CIRCLE_SECTOR_AREA_2D computes the area of a");
        Console.WriteLine("    circular sector, defined by joining the endpoints");
        Console.WriteLine("    of a circular arc to the center.");
        Console.WriteLine("  CIRCLE_TRIANGLE_AREA_2D computes the signed area of a");
        Console.WriteLine("    triangle, defined by joining the endpoints");
        Console.WriteLine("    of a circular arc and the center.");
        Console.WriteLine("");
        Console.WriteLine("      R       Theta1 Theta2        "
                          + "Sector       Triangle     Half Lens");
        Console.WriteLine("");

        for (i = 0; i <= 12; i++)
        {
            theta1 = 0.0;
            theta2 = i * 2.0 * pi / 12.0;

            area1 = Circle.circle_sector_area_2d(r, theta1, theta2);

            area2 = Circle.circle_triangle_area_2d(r, theta1, theta2);

            area3 = Lens.lens_half_area_2d(r, theta1, theta2);

            Console.WriteLine("  " + r.ToString(CultureInfo.InvariantCulture).PadLeft(9)
                                   + "  " + theta1.ToString(CultureInfo.InvariantCulture).PadLeft(9)
                                   + "  " + theta2.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + area1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + area2.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + area3.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

    }

    public static void test12()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST12 tests LENS_HALF_AREA_2D, LENS_HALF_H_AREA_2D, LENS_HALF_W_AREA_2D.
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
        double area1;
        double area2;
        double area3;
        double pi = 3.141592653589793;
        double h;
        int i;
        double r;
        double theta1;
        double theta2;
        double w;

        r = 50.0;

        Console.WriteLine("");
        Console.WriteLine("TEST12");
        Console.WriteLine("  For the area of a circular half lens,");
        Console.WriteLine("  LENS_HALF_AREA_2D uses two angles;");
        Console.WriteLine("  LENS_HALF_H_AREA_2D works from the height;");
        Console.WriteLine("  LENS_HALF_W_AREA_2D works from the width.");
        Console.WriteLine("");
        Console.WriteLine("  The circle has radius R = " + r + "");
        Console.WriteLine("");
        Console.WriteLine("  THETA1 THETA2  H     W  Area(THETA) Area(H)  Area(W)");
        Console.WriteLine("");

        for (i = 0; i <= 12; i++)
        {
            theta1 = 0.0;
            theta2 = i * 2.0 * pi / 12.0;
            w = 2.0 * r * Math.Sin(0.5 * (theta2 - theta1));
            h = r * (1.0 - Math.Cos(0.5 * (theta2 - theta1)));

            area1 = Lens.lens_half_area_2d(r, theta1, theta2);

            area2 = Lens.lens_half_h_area_2d(r, h);

            area3 = Lens.lens_half_w_area_2d(r, w);

            Console.WriteLine("  " + theta1.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + theta2.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + h.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + w.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + area1.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + area2.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + area3.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

    }

    public static void test13()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST13 tests CIRCLE_SECTOR, CIRCLE_SECTOR_AREA_2D.
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
        double area;
        double[] center = new double[2];
        double[] center_test =
            {
                0.0, 0.0,
                0.0, 0.0,
                0.0, 0.0,
                0.0, 0.0
            }
            ;
        int dim;
        int dim_num = 2;
        int i;
        int j;
        string name = "";
        int num;
        int nr;
        int nrhi = 5;
        int nrlo = 1;
        double pi = 3.141592653589793;
        double radius;
        double[] radius_test =
            {
                1.0, 2.0, 4.0, 8.0
            }
            ;
        double result;
        int test_num = 4;
        double theta1;
        double[] theta1_test =
            {
                0.0, 0.0, 0.0, 0.0
            }
            ;
        double theta2;
        double[] theta2_test =
            {
                2.0, 1.0, 0.5, 0.25
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("TEST13");
        Console.WriteLine("  CIRCLE_SECTOR_AREA_2D computes the area");
        Console.WriteLine("    of a circular sector.");
        Console.WriteLine("  CIRCLE_SECTOR estimates an integral ");
        Console.WriteLine("    in a circular sector.");
        Console.WriteLine("");
        Console.WriteLine("  The user can specify NR, the number of radial values");
        Console.WriteLine("  used to approximated the integral.");
        Console.WriteLine("");
        Console.WriteLine("  In this test, computations will use values of NR");
        Console.WriteLine("  from " + nrlo + "");
        Console.WriteLine("  to   " + nrhi + "");
        Console.WriteLine("");

        for (i = 0; i < test_num; i++)
        {
            for (dim = 0; dim < dim_num; dim++)
            {
                center[dim] = center_test[dim + i * dim_num];
            }

            radius = radius_test[i];
            theta1 = theta1_test[i] * pi;
            theta2 = theta2_test[i] * pi;

            area = Circle.circle_sector_area_2d(radius, theta1, theta2);

            Console.WriteLine("");
            Console.WriteLine("       CENTER      RADIUS  THETA1  THETA2  Area");
            Console.WriteLine("");
            Console.WriteLine(center[0].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                              + center[1].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                              + radius.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                              + theta1.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                              + theta2.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                              + area.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");
            Console.WriteLine("");
            string cout = "       F   ";
            for (nr = nrlo; nr <= nrhi; nr++)
            {
                cout += "      " + nr.ToString(CultureInfo.InvariantCulture).PadLeft(2) + "      ";
            }

            Console.WriteLine(cout);
            Console.WriteLine("");

            num = functions.function_2d_num();

            for (j = 1; j <= num; j++)
            {
                int function_2d_index = j;

                functions.function_2d_name(function_2d_index, ref name);

                cout = "  " + name;

                for (nr = nrlo; nr <= nrhi; nr++)
                {
                    result = Circle.circle_sector(function_2d_index, functions.function_2d, center, radius, theta1,
                        theta2,
                        nr);
                    cout += result.ToString(CultureInfo.InvariantCulture).PadLeft(14);
                }

                Console.WriteLine(cout);
            }
        }

    }

    public static void test14()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST14 tests CIRCLE_SECTOR, CIRCLE_RT_SET, CIRCLE_RT_SUM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double area1;
        double area2;
        double area3;
        double[] center = new double[2];
        double[] center_test =
            {
                0.0, 0.0,
                0.0, 0.0,
                0.0, 0.0,
                0.0, 0.0
            }
            ;
        int dim = 0;
        int dim_num = 2;
        int i = 0;
        int j = 0;
        string name = "";
        int nc = 0;
        int num = 0;
        int nr = 0;
        int nr2 = 0;
        int nt = 0;
        double pi = 3.141592653589793;
        double[] ra = new double[5];
        double radius;
        double[] radius_test =
            {
                1.0, 2.0, 4.0, 8.0
            }
            ;
        double result1 = 0;
        double result2 = 0;
        double resulta = 0;
        double resultb = 0;
        int rule = 0;
        double[] rw = new double[5];
        double[] ta = new double[20];
        int test_num = 4;
        double theta1 = 0;
        double[] theta1_test =
            {
                0.0, 0.0, 0.0, 0.0
            }
            ;
        double theta2 = 0;
        double[] theta2_test =
            {
                2.0, 1.0, 0.5, 0.25
            }
            ;
        double theta3 = 0;
        double[] tw = new double[20];
        double zw = 0;

        nr = 5;

        rule = 9;
        Circle.circle_rt_size(rule, ref nr2, ref nt, ref nc);
        Circle.circle_rt_set(rule, nr2, nt, nc, ref ra, ref rw, ref ta, ref tw, ref zw);

        Console.WriteLine("");
        Console.WriteLine("TEST14");
        Console.WriteLine("  CIRCLE_SECTOR estimates integrals in a circular sector.");
        Console.WriteLine("  CIRCLE_RT_SET sets an integration rule in a circle.");
        Console.WriteLine("  CIRCLE_RT_SUM uses an integration rule in a circle.");
        Console.WriteLine("");
        Console.WriteLine("  To test CIRCLE_SECTOR, we estimate an integral over");
        Console.WriteLine("  a sector, and over its complement and add the results");
        Console.WriteLine("  to get RESULT1.");
        Console.WriteLine("");
        Console.WriteLine("  We also estimate the integral over the whole circle");
        Console.WriteLine("  using CIRCLE_RT_SET and CIRCLE_RT_SUM to get RESULT2.");
        Console.WriteLine("");
        Console.WriteLine("  We will then compare RESULT1 and RESULT2.");
        Console.WriteLine("");
        Console.WriteLine("  CIRCLE_SECTOR computations will use NR = " + nr + "");
        Console.WriteLine("  CIRCLE_RT_SET/CIRCLE_RT_SUM will use rule " + rule + "");
        Console.WriteLine("");
        Console.WriteLine("  'Sector1' and 'Sector2' are the CIRCLE_SECTOR");
        Console.WriteLine("  computations");
        Console.WriteLine("  for the sector and its complement.");
        Console.WriteLine("  'Sum' is the sum of Sector1 and Sector2.");
        Console.WriteLine("  'Circle' is the computation for ");
        Console.WriteLine("  CIRCLE_RT_SET + CIRCLE_RT_SUM.");
        Console.WriteLine("");

        for (i = 0; i < test_num; i++)
        {
            for (dim = 0; dim < dim_num; dim++)
            {
                center[dim] = center_test[dim + i * dim_num];
            }

            radius = radius_test[i];

            theta1 = theta1_test[i] * pi;
            theta2 = theta2_test[i] * pi;
            theta3 = theta2 + 2.0 * pi - (theta2 - theta1);

            area1 = Circle.circle_sector_area_2d(radius, theta1, theta2);
            area2 = Circle.circle_sector_area_2d(radius, theta2, theta3);
            area3 = Circle.circle_area_2d(radius);

            Console.WriteLine("");
            Console.WriteLine("      CENTER       RADIUS   THETA1   THETA2   Area1   Area2  Circle");
            Console.WriteLine("");
            Console.WriteLine(center[0].ToString(CultureInfo.InvariantCulture).PadLeft(9)
                              + center[1].ToString(CultureInfo.InvariantCulture).PadLeft(9)
                              + radius.ToString(CultureInfo.InvariantCulture).PadLeft(9)
                              + theta1.ToString(CultureInfo.InvariantCulture).PadLeft(9)
                              + theta2.ToString(CultureInfo.InvariantCulture).PadLeft(9)
                              + area1.ToString(CultureInfo.InvariantCulture).PadLeft(9)
                              + area2.ToString(CultureInfo.InvariantCulture).PadLeft(9)
                              + area3.ToString(CultureInfo.InvariantCulture).PadLeft(9) + "");
            Console.WriteLine("");
            Console.WriteLine("       F   Sector1       Sector2         Sum         Circle");
            Console.WriteLine("");

            num = functions.function_2d_num();

            for (j = 1; j <= num; j++)
            {
                int function_2d_index = j;
                functions.function_2d_name(function_2d_index, ref name);

                resulta = Circle.circle_sector(function_2d_index, functions.function_2d, center, radius, theta1,
                    theta2, nr);

                resultb = Circle.circle_sector(function_2d_index, functions.function_2d, center, radius, theta2,
                    theta3, nr);

                result1 = resulta + resultb;

                result2 = Circle.circle_rt_sum(function_2d_index, functions.function_2d, center, radius, nr2, ra,
                    rw,
                    nt, ta, tw, zw);

                Console.WriteLine("  " + name
                                       + resulta.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + resultb.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + result1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + result2.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }
        }
    }

    public static void test15()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST15 tests CIRCLE_RT_SET and CIRCLE_RT_SUM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] center =
            {
                1.0, 1.0
            }
            ;
        int i = 0;
        int ihi = 0;
        int ilo = 0;
        string name = "";
        int nc = 0;
        int num = 0;
        int nr = 0;
        int nt = 0;
        double r = 1.0;
        double[] ra;
        double result = 0;
        double[] rw;
        int rule = 0;
        int rule_max = 9;
        double[] ta;
        double[] tw;
        double zw = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST15");
        Console.WriteLine("  For R, Theta product rules on the unit circle,");
        Console.WriteLine("  CIRCLE_RT_SET sets a rule.");
        Console.WriteLine("  CIRCLE_RT_SUM uses the rule in an arbitrary circle.");
        Console.WriteLine("");
        Console.WriteLine("  We use a radius " + r + "");
        Console.WriteLine("  and center:");
        Console.WriteLine("  CENTER = " + center[0].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                        + center[1].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine("");

        for (ilo = 1; ilo <= rule_max; ilo += 5)
        {
            ihi = Math.Min(ilo + 4, rule_max);

            Console.WriteLine("");
            string cout = "  Rule:  ";
            for (rule = ilo; rule <= ihi; rule++)
            {
                cout += "       " + rule.ToString(CultureInfo.InvariantCulture).PadLeft(7);
            }

            Console.WriteLine(cout);
            Console.WriteLine("Function");

            num = functions.function_2d_num();

            for (i = 1; i <= num; i++)
            {
                int function_2d_index = i;
                functions.function_2d_name(function_2d_index, ref name);
                cout = "  " + name;

                for (rule = ilo; rule <= ihi; rule++)
                {
                    Circle.circle_rt_size(rule, ref nr, ref nt, ref nc);

                    ra = new double[nr];
                    rw = new double[nr];
                    ta = new double[nt];
                    tw = new double[nt];

                    Circle.circle_rt_set(rule, nr, nt, nc, ref ra, ref rw, ref ta, ref tw, ref zw);

                    result = Circle.circle_rt_sum(function_2d_index, functions.function_2d, center, r, nr, ra, rw,
                        nt, ta,
                        tw, zw);
                    cout += result.ToString(CultureInfo.InvariantCulture).PadLeft(14);
                }

                Console.WriteLine(cout);
            }
        }
    }

    public static void test16()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST16 tests CIRCLE_XY_SET and CIRCLE_XY_SUM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 March 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] center =
            {
                1.0, 1.0
            }
            ;
        int i;
        int ihi;
        int ilo;
        string name = "";
        int num;
        int order;
        double r;
        double result;
        int rule;
        int rule_max = 13;
        double[] weight;
        double[] xtab;
        double[] ytab;

        r = 1.0;
        Console.WriteLine("");
        Console.WriteLine("TEST16");
        Console.WriteLine("  CIRCLE_XY_SET sets a quadrature rule");
        Console.WriteLine("    for the unit circle.");
        Console.WriteLine("  CIRCLE_XY_SUM evaluates the quadrature rule");
        Console.WriteLine("    in an arbitrary circle.");
        Console.WriteLine("");
        Console.WriteLine("  We use a radius " + r + "");
        Console.WriteLine("  and center:");
        Console.WriteLine("  CENTER = (" + center[0] + ", " + center[1] + ").");
        Console.WriteLine("");

        for (ilo = 1; ilo <= rule_max; ilo += 5)
        {
            ihi = Math.Min(ilo + 4, rule_max);

            Console.WriteLine("");
            string cout = "  Rule:  ";
            for (rule = ilo; rule <= ihi; rule++)
            {
                cout += "       " + rule.ToString(CultureInfo.InvariantCulture).PadLeft(7);
            }

            Console.WriteLine(cout);
            Console.WriteLine(" 'Function");

            num = functions.function_2d_num();

            for (i = 1; i <= num; i++)
            {
                int function_2d_index = i;
                functions.function_2d_name(function_2d_index, ref name);

                cout = "  " + name;

                for (rule = ilo; rule <= ihi; rule++)
                {
                    order = Circle.circle_xy_size(rule);

                    xtab = new double[order];
                    ytab = new double[order];
                    weight = new double[order];

                    Circle.circle_xy_set(rule, order, xtab, ytab, weight);

                    result = Circle.circle_xy_sum(function_2d_index, functions.function_2d, center, r, order, xtab,
                        ytab,
                        weight);

                    cout += result.ToString(CultureInfo.InvariantCulture).PadLeft(14);

                }

                Console.WriteLine(cout);
            }
        }
    }

    public static void test163()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST163 tests the rules for CN with Gegenbauer weight on monomials.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 January 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int TEST_NUM = 5;

        double alpha;
        double[] alpha_test =
            {
                -0.5, 0.0, 0.5, 1.0, 1.5
            }
            ;
        int[] expon;
        int n;
        int test;

        Console.WriteLine("");
        Console.WriteLine("TEST163");
        Console.WriteLine("  Demonstrate the use of quadrature rules for the region");
        Console.WriteLine("  CN_GEG, that is, the hypercube [-1,+1]^N, with the");
        Console.WriteLine("  weight W(ALPHA;X) = product ( 1 <= I <= N )");
        Console.WriteLine("    (1-X(I)^2)^ALPHA");
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
            expon = new int[n];

            for (test = 0; test < TEST_NUM; test++)
            {
                alpha = alpha_test[test];

                typeMethods.i4vec_zero(n, ref expon);
                cn.cn_geg_test(n, alpha, expon);
            }

            for (test = 0; test < TEST_NUM; test++)
            {
                alpha = alpha_test[test];

                typeMethods.i4vec_zero(n, ref expon);
                expon[n - 1] = 1;
                cn.cn_geg_test(n, alpha, expon);
            }

            switch (n)
            {
                case >= 2:
                {
                    for (test = 0; test < TEST_NUM; test++)
                    {
                        alpha = alpha_test[test];

                        typeMethods.i4vec_zero(n, ref expon);
                        expon[0] = 1;
                        expon[1] = 1;
                        cn.cn_geg_test(n, alpha, expon);
                    }

                    break;
                }
            }

            for (test = 0; test < TEST_NUM; test++)
            {
                alpha = alpha_test[test];

                typeMethods.i4vec_zero(n, ref expon);
                expon[0] = 2;
                cn.cn_geg_test(n, alpha, expon);

            }
        }

    }

    public static void test165()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST165 tests the rules for CN with Jacobi weight on monomials.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 January 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int TEST_NUM = 4;

        double alpha;
        double[] alpha_test =
            {
                0.0, 1.0, 0.0, 0.5
            }
            ;
        double beta;
        double[] beta_test =
            {
                0.0, 0.0, 2.0, 1.5
            }
            ;
        int[] expon;
        int n;
        int test;

        Console.WriteLine("");
        Console.WriteLine("TEST165");
        Console.WriteLine("  Demonstrate the use of quadrature rules for the region");
        Console.WriteLine("  CN_JAC, that is, the hypercube [-1,+1]^N, with the");
        Console.WriteLine("  weight W(ALPHA,BETA;X) = product ( 1 <= I <= N )");
        Console.WriteLine("    (1-X(I))^ALPHA (1+X(I))^BETA");
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
            expon = new int[n];

            for (test = 0; test < TEST_NUM; test++)
            {
                alpha = alpha_test[test];
                beta = beta_test[test];

                typeMethods.i4vec_zero(n, ref expon);
                cn.cn_jac_test(n, alpha, beta, expon);
            }

            for (test = 0; test < TEST_NUM; test++)
            {
                alpha = alpha_test[test];
                beta = beta_test[test];

                typeMethods.i4vec_zero(n, ref expon);
                expon[n - 1] = 1;
                cn.cn_jac_test(n, alpha, beta, expon);
            }

            switch (n)
            {
                case >= 2:
                {
                    for (test = 0; test < TEST_NUM; test++)
                    {
                        alpha = alpha_test[test];
                        beta = beta_test[test];

                        typeMethods.i4vec_zero(n, ref expon);
                        expon[0] = 1;
                        expon[1] = 1;
                        cn.cn_jac_test(n, alpha, beta, expon);
                    }

                    break;
                }
            }

            for (test = 0; test < TEST_NUM; test++)
            {
                alpha = alpha_test[test];
                beta = beta_test[test];

                typeMethods.i4vec_zero(n, ref expon);
                expon[0] = 2;
                cn.cn_jac_test(n, alpha, beta, expon);

            }
        }
    }

    public static void test167()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST167 tests the rules for CN with Legendre weight on monomials.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 March 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int[] expon;
        int n;

        Console.WriteLine("");
        Console.WriteLine("TEST167");
        Console.WriteLine("  Demonstrate the use of quadrature rules for the region");
        Console.WriteLine("  CN_LEG, that is, the hypercube [-1,+1]^N, with the");
        Console.WriteLine("  Legendre weight W(X) = 1.");
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
            expon = new int[n];

            typeMethods.i4vec_zero(n, ref expon);
            cn.cn_leg_test(n, expon);

            typeMethods.i4vec_zero(n, ref expon);
            expon[n - 1] = 1;
            cn.cn_leg_test(n, expon);

            switch (n)
            {
                case >= 2:
                    typeMethods.i4vec_zero(n, ref expon);
                    expon[0] = 1;
                    expon[1] = 1;
                    cn.cn_leg_test(n, expon);
                    break;
            }

            typeMethods.i4vec_zero(n, ref expon);
            expon[0] = 2;
            cn.cn_leg_test(n, expon);

            typeMethods.i4vec_zero(n, ref expon);
            expon[0] = 3;
            cn.cn_leg_test(n, expon);

            typeMethods.i4vec_zero(n, ref expon);
            expon[n - 1] = 4;
            cn.cn_leg_test(n, expon);

            switch (n)
            {
                case >= 2:
                    typeMethods.i4vec_zero(n, ref expon);
                    expon[0] = 3;
                    expon[1] = 2;
                    cn.cn_leg_test(n, expon);
                    break;
            }
        }

    }

    public static void test17()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST17 tests CONE_UNIT_3D, CONE_VOLUME_3D.
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
        double h;
        int i;
        string name = "";
        int num;
        double r = 1.0;
        double result;

        h = 1.0;

        Console.WriteLine("");
        Console.WriteLine("TEST17");
        Console.WriteLine("  CONE_UNIT_3D approximates integrals in a unit cone.");
        Console.WriteLine("");
        Console.WriteLine("  Volume = " + Cone.cone_volume_3d(r, h) + "");
        Console.WriteLine("");
        Console.WriteLine("    F(X)    CONE_3D");
        Console.WriteLine("");

        num = functions.function_3d_num();

        for (i = 1; i <= num; i++)
        {
            int function_3d_index = i;
            functions.function_3d_name(function_3d_index, ref name);

            result = Cone.cone_unit_3d(function_3d_index, functions.function_3d);

            Console.WriteLine("  " + name
                                   + "  " + result.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    public static void test18()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST18 tests CUBE_SHELL_ND.
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
        int n;
        const int N_MAX = 4;
        string name = "";
        int num;
        double r1;
        double[] r1_test = { 0.0, 1.0 };
        double r2;
        double[] r2_test = { 1.0, 2.0 };
        double result;
        int test;
        int test_num = 2;

        Console.WriteLine("");
        Console.WriteLine("TEST18");
        Console.WriteLine("  CUBE_SHELL_ND approximates integrals in a");
        Console.WriteLine("    cubical shell in ND.");

        for (test = 0; test < test_num; test++)
        {
            r1 = r1_test[test];
            r2 = r2_test[test];

            Console.WriteLine("");
            Console.WriteLine("  Inner radius = " + r1 + "");
            Console.WriteLine("  Outer radius = " + r2 + "");
            Console.WriteLine("");

            for (n = 2; n <= n_max; n++)
            {
                Console.WriteLine("");
                Console.WriteLine("  Spatial dimension N = " + n + "");
                Console.WriteLine("  Volume = " + Cube.cube_shell_volume_nd(n, r1, r2) + "");
                Console.WriteLine("");
                Console.WriteLine("    F(X)      CUBE_SHELL_ND");
                Console.WriteLine("");

                num = functions.function_nd_num();

                for (i = 1; i <= num; i++)
                {
                    int function_nd_index = i;
                    functions.function_nd_name(function_nd_index, ref name);

                    result = Cube.cube_shell_nd(function_nd_index, functions.function_nd, n, r1, r2);

                    Console.WriteLine("  " + name
                                           + result.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
                }
            }
        }
    }

    public static void test19()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST19 tests CUBE_UNIT_3D, QMULT_3D, RECTANGLE_3D.
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
        double a1;
        double[] a = new double[3];
        double b1;
        double[] b = new double[3];
        int i;
        string name = "";
        int num;
        double result1 = 0;
        double result2 = 0;
        double result3;

        a1 = -1.0;
        b1 = +1.0;

        a[0] = -1.0;
        a[1] = -1.0;
        a[2] = -1.0;
        b[0] = 1.0;
        b[1] = 1.0;
        b[2] = 1.0;

        Console.WriteLine("");
        Console.WriteLine("TEST19");
        Console.WriteLine("  CUBE_UNIT_3D approximates integrals");
        Console.WriteLine("    in the unit cube in 3D.");
        Console.WriteLine("  QMULT_3D approximates triple integrals.");
        Console.WriteLine("  RECTANGLE_3D approximates integrals");
        Console.WriteLine("    in a rectangular block.");
        Console.WriteLine("");
        Console.WriteLine("    F(X)      CUBE_UNIT_3D    QMULT_3D        RECTANGLE_3D");
        Console.WriteLine("");

        num = functions.function_3d_num();

        for (i = 1; i <= num; i++)
        {
            int function_3d_index = i;
            functions.function_3d_name(function_3d_index, ref name);

            result1 = Cube.cube_unit_3d(function_3d_index, functions.function_3d);
            result2 = QMult.qmult_3d(function_3d_index, functions.function_3d, a1, b1, integrationlimits.fu18,
                integrationlimits.fl18, integrationlimits.fu28, integrationlimits.fl28);
            result3 = Rectangle.rectangle_3d(function_3d_index, functions.function_3d, a, b);

            Console.WriteLine("  " + name
                                   + "  " + result1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + result2.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + result3.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }


}