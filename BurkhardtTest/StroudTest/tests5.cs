﻿using System;
using System.Globalization;
using Burkardt.Stroud;
using Burkardt.Types;

namespace StroudTest;

public static class tests5
{
    public static void test40()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST40 tests TETRA_UNIT_SET and TETRA_SUM.
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
        int i;
        int ilo;
        string name = "";
        const int rule_max = 8;
        double[] x = { 1.0, 4.0, 1.0, 1.0 };
        double[] y = { 2.0, 2.0, 3.0, 2.0 };
        double[] z = { 6.0, 6.0, 6.0, 8.0 };

        Console.WriteLine("");
        Console.WriteLine("TEST40");
        Console.WriteLine("  TETRA_UNIT_SET sets quadrature rules");
        Console.WriteLine("    for the unit tetrahedron;");
        Console.WriteLine("  TETRA_SUM applies them to an arbitrary tetrahedron.");
        Console.WriteLine("");
        Console.WriteLine("  Tetrahedron vertices:");
        Console.WriteLine("");
        for (i = 0; i < 4; i++)
        {
            Console.WriteLine("  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + y[i].ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + z[i].ToString(CultureInfo.InvariantCulture).PadLeft(6) + "");
        }


        for (ilo = 1; ilo <= rule_max; ilo += 5)
        {
            int ihi = Math.Min(ilo + 4, rule_max);

            Console.WriteLine("");
            string cout = "  Rule:   ";
            int rule;
            for (rule = ilo; rule <= ihi; rule++)
            {
                cout += "       " + rule.ToString().PadLeft(7);
            }

            Console.WriteLine(cout);
            Console.WriteLine("  Function");
            Console.WriteLine("");

            int num = functions.function_3d_num();

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

                    double result = Tetrahedron.tetra_sum(function_3d_index, functions.function_3d, x, y, z, order, xtab,
                        ytab, ztab,
                        weight);

                    cout += result.ToString(CultureInfo.InvariantCulture).PadLeft(14);
                }

                Console.WriteLine(cout);
            }
        }
    }

    public static void test41()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST41 tests TRIANGLE_UNIT_SET, TRIANGLE_SUB.
        //
        //  Discussion:
        //
        //    Break up the triangle into NSUB*NSUB equal subtriangles.  Approximate 
        //    the integral over the triangle by the sum of the integrals over each
        //    subtriangle.
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
        int i;
        string name = "";
        double[] xval = { 0.0, 0.0, 1.0 };
        double[] yval = { 0.0, 1.0, 0.0 };

        Console.WriteLine("");
        Console.WriteLine("TEST41");
        Console.WriteLine("  TRIANGLE_UNIT_SET sets up a quadrature rule");
        Console.WriteLine("    on a triangle.");
        Console.WriteLine("  TRIANGLE_SUB applies it to subtriangles of an");
        Console.WriteLine("    arbitrary triangle.");
        Console.WriteLine("");
        Console.WriteLine("  Triangle vertices:");
        Console.WriteLine("");
        Console.WriteLine(xval[0].ToString(CultureInfo.InvariantCulture).PadLeft(14) + yval[0].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine(xval[1].ToString(CultureInfo.InvariantCulture).PadLeft(14) + yval[1].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        Console.WriteLine(xval[2].ToString(CultureInfo.InvariantCulture).PadLeft(14) + yval[2].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        //
        //  Get the quadrature abscissas and weights for a unit triangle.
        //
        int rule = 3;
        int order = Triangle.triangle_unit_size(rule);

        double[] xtab = new double[order];
        double[] ytab = new double[order];
        double[] weight = new double[order];

        Triangle.triangle_unit_set(rule, order, ref xtab, ref ytab, ref weight);

        Console.WriteLine("");
        Console.WriteLine("  Using unit triangle quadrature rule " + rule + "");
        Console.WriteLine("  Rule order = " + order + "");
        Console.WriteLine("");
        Console.WriteLine("  Function Nsub  Result");
        Console.WriteLine("");
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
            int nsub;
            for (nsub = 1; nsub <= 5; nsub++)
            {
                double result = Triangle.triangle_sub(function_2d_index, functions.function_2d, xval, yval, nsub, order,
                    xtab,
                    ytab, weight);
                Console.WriteLine("  " + name + nsub.ToString().PadLeft(4) + result.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }
        }

    }

    public static void test42()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST42 tests TRIANGLE_UNIT_SET and TRIANGLE_UNIT_SUM.
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
        const int rule_max = 20;

        Console.WriteLine("");
        Console.WriteLine("TEST42");
        Console.WriteLine("  TRIANGLE_UNIT_SET sets up a quadrature");
        Console.WriteLine("    in the unit triangle,");
        Console.WriteLine("  TRIANGLE_UNIT_SUM applies it.");
        Console.WriteLine("");

        for (ilo = 1; ilo <= rule_max; ilo += 5)
        {
            int ihi = Math.Min(ilo + 4, rule_max);

            Console.WriteLine("");
            string cout = "  Rule:   ";
            int rule;
            for (rule = ilo; rule <= ihi; rule++)
            {
                cout += "       " + rule;
            }

            Console.WriteLine(cout);
            Console.WriteLine("Function");
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
                    int order = Triangle.triangle_unit_size(rule);

                    double[] xtab = new double[order];
                    double[] ytab = new double[order];
                    double[] weight = new double[order];

                    Triangle.triangle_unit_set(rule, order, ref xtab, ref ytab, ref weight);

                    double result = Triangle.triangle_unit_sum(function_2d_index, functions.function_2d, order, xtab, ytab,
                        weight);

                    cout += result.ToString(CultureInfo.InvariantCulture).PadLeft(14);
                }

                Console.WriteLine(cout);
            }
        }
    }

    public static void test425()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST425 tests TRIANGLE_UNIT_SET.
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
        int a;
        const int rule_max = 20;

        Console.WriteLine("");
        Console.WriteLine("TEST425");
        Console.WriteLine("  TRIANGLE_UNIT_SET sets up a quadrature");
        Console.WriteLine("    in the unit triangle,");
        Console.WriteLine("");
        Console.WriteLine("  Estimate integral of X^A * Y^B.");

        for (a = 0; a <= 10; a++)
        {
            int b;
            for (b = 0; b <= 10 - a; b++)
            {
                double coef = (a + b + 2) * (double)(a + b + 1);
                int i;
                for (i = 1; i <= b; i++)
                {
                    coef = coef * (a + i) / i;
                }

                Console.WriteLine("");
                Console.WriteLine("  A = " + a
                                           + "  B = " + b + "");
                Console.WriteLine("");
                Console.WriteLine("  Rule       QUAD           ERROR");
                Console.WriteLine("");

                int rule;
                for (rule = 1; rule <= rule_max; rule++)
                {
                    int order = Triangle.triangle_unit_size(rule);

                    double[] xtab = new double[order];
                    double[] ytab = new double[order];
                    double[] weight = new double[order];

                    Triangle.triangle_unit_set(rule, order, ref xtab, ref ytab, ref weight);

                    double quad = 0.0;

                    for (i = 0; i < order; i++)
                    {
                        double value = 0;
                        switch (a)
                        {
                            case 0 when b == 0:
                                value = coef;
                                break;
                            case 0 when true:
                                value = coef * Math.Pow(ytab[i], b);
                                break;
                            default:
                            {
                                if (b == 0)
                                {
                                    value = coef * Math.Pow(xtab[i], a);
                                }
                                else
                                {
                                    value = coef * Math.Pow(xtab[i], a) * Math.Pow(ytab[i], b);
                                }

                                break;
                            }
                        }

                        quad += 0.5 * weight[i] * value;
                    }

                    double exact = 1.0;
                    double err = Math.Abs(exact - quad);
                    Console.WriteLine("  " + rule.ToString().PadLeft(4)
                                           + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                           + "  " + err.ToString(CultureInfo.InvariantCulture).PadLeft(11) + "");
                }
            }
        }
    }

    public static void test43()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST43 tests TRIANGLE_UNIT_PRODUCT_SET and TRIANGLE_UNIT_SUM.
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
        int ilo;
        string name = "";
        const int rule_max = 8;

        Console.WriteLine("");
        Console.WriteLine("TEST43");
        Console.WriteLine("  TRIANGLE_UNIT_PRODUCT_SET sets up a product quadrature");
        Console.WriteLine("    rule in the unit triangle,");
        Console.WriteLine("  TRIANGLE_UNIT_SUM applies it.");
        Console.WriteLine("");

        for (ilo = 1; ilo <= rule_max; ilo += 5)
        {
            int ihi = Math.Min(ilo + 4, rule_max);

            Console.WriteLine("");
            string cout = "  Rule Order: ";
            int rule;
            for (rule = ilo; rule <= ihi; rule++)
            {
                cout += rule.ToString().PadLeft(6);
            }

            Console.WriteLine(cout);
            Console.WriteLine("Function");
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
                    int order = Triangle.triangle_unit_product_size(rule);

                    double[] xtab = new double[order];
                    double[] ytab = new double[order];
                    double[] weight = new double[order];

                    Triangle.triangle_unit_product_set(rule, order, ref xtab, ref ytab, ref weight);

                    double result = Triangle.triangle_unit_sum(function_2d_index, functions.function_2d, order, xtab, ytab,
                        weight);

                    cout += result.ToString(CultureInfo.InvariantCulture).PadLeft(14);

                }

                Console.WriteLine(cout);
            }
        }
    }

    public static void test44()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST44 tests TRIANGLE_UNIT_SET and TRIANGLE_SUM.
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
        int ilo;
        string name = "";
        const int rule_max = 20;
        double[] xval = { 1.0, 3.0, 1.0 };
        double[] yval = { 1.0, 1.0, 4.0 };

        Console.WriteLine("");
        Console.WriteLine("TEST44");
        Console.WriteLine("  TRIANGLE_UNIT_SET sets up quadrature");
        Console.WriteLine("    in the unit triangle,");
        Console.WriteLine("  TRIANGLE_SUM applies it to an arbitrary triangle.");
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

            Console.WriteLine("Function");
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
                    int order = Triangle.triangle_unit_size(rule);

                    double[] xtab = new double[order];
                    double[] ytab = new double[order];
                    double[] weight = new double[order];

                    Triangle.triangle_unit_set(rule, order, ref xtab, ref ytab, ref weight);

                    double result = Triangle.triangle_sum(function_2d_index, functions.function_2d, xval, yval, order,
                        xtab, ytab,
                        weight);

                    cout += result.ToString(CultureInfo.InvariantCulture).PadLeft(14);

                }

                Console.WriteLine(cout);
            }
        }
    }

    public static void test45()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST45 tests TORUS_1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int j;
        int j2;
        string name = "";

        const double r1 = 0.5;
        const double r2 = 1.0;

        Console.WriteLine("");
        Console.WriteLine("TEST45");
        Console.WriteLine("  TORUS_1 approximates integrals on a torus.");
        Console.WriteLine("");
        Console.WriteLine("  The order N will be varied.");
        Console.WriteLine("");
        Console.WriteLine("  Inner radius = " + r1 + "");
        Console.WriteLine("  Outer radius = " + r2 + "");
        Console.WriteLine("  Area = " + Torus.torus_area_3d(r1, r2) + "");
        Console.WriteLine("");
        string cout = "  " + "  F(X)  ";
        for (j = 1; j <= 5; j++)
        {
            j2 = 2 * (j - 1);
            cout += ((int)Math.Pow(2, j2)).ToString().PadLeft(14);
        }

        Console.WriteLine(cout);
        Console.WriteLine("");

        int num = functions.function_3d_num();

        for (i = 1; i <= num; i++)
        {
            int function_3d_index = i;
            functions.function_3d_name(function_3d_index, ref name);

            cout = "  " + name;

            for (j = 1; j <= 5; j++)
            {
                j2 = 2 * (j - 1);
                int n = (int)Math.Pow(2, j2);
                double result = Torus.torus_1(function_3d_index, functions.function_3d, r1, r2, n);
                cout += result.ToString(CultureInfo.InvariantCulture).PadLeft(14);
            }

            Console.WriteLine(cout);
        }
    }

    public static void test46()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST46 tests TORUS_5S2, TORUS_6S2 and TORUS_14S.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        string name = "";

        const double r1 = 0.5;
        const double r2 = 1.0;

        Console.WriteLine("");
        Console.WriteLine("TEST46");
        Console.WriteLine("  For the interior of a torus,");
        Console.WriteLine("  TORUS_5S2,");
        Console.WriteLine("  TORUS_6S2, and");
        Console.WriteLine("  TORUS_5S2 approximate integrals.");
        Console.WriteLine("");
        Console.WriteLine("  Inner radius = " + r1 + "");
        Console.WriteLine("  Outer radius = " + r2 + "");
        Console.WriteLine("  Volume = " + Torus.torus_volume_3d(r1, r2) + "");
        Console.WriteLine("");
        Console.WriteLine("    Rule:        #5S2          #6S2          #14S");
        Console.WriteLine("    F(X)");
        Console.WriteLine("");

        int num = functions.function_3d_num();

        for (i = 1; i <= num; i++)
        {
            int function_3d_index = i;
            functions.function_3d_name(function_3d_index, ref name);

            double result1 = Torus.torus_5s2(function_3d_index, functions.function_3d, r1, r2);
            double result2 = Torus.torus_6s2(function_3d_index, functions.function_3d, r1, r2);
            double result3 = Torus.torus_14s(function_3d_index, functions.function_3d, r1, r2);

            Console.WriteLine("  " + name
                                   + "  " + result1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + result2.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + result3.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    public static void test47()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST47 tests TORUS_SQUARE_5C2 and TORUS_SQUARE_14C.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        string name = "";

        const double r1 = 1.0;
        const double r2 = 0.125;

        Console.WriteLine("");
        Console.WriteLine("TEST47");
        Console.WriteLine("  For integrals inside a torus with square cross-section:");
        Console.WriteLine("  TORUS_SQUARE_5C2 approximates the integral;");
        Console.WriteLine("  TORUS_SQUARE_14C approximates the integral.");
        Console.WriteLine("");
        Console.WriteLine("  Inner radius = " + r1 + "");
        Console.WriteLine("  Outer radius = " + r2 + "");
        Console.WriteLine("  Volume = " + Torus.torus_square_volume_3d(r1, r2) + "");
        Console.WriteLine("");
        Console.WriteLine("    F(X)    5C2           14C");
        Console.WriteLine("");

        int num = functions.function_3d_num();

        for (i = 1; i <= num; i++)
        {
            int function_3d_index = i;
            functions.function_3d_name(function_3d_index, ref name);

            double result1 = Torus.torus_square_5c2(function_3d_index, functions.function_3d, r1, r2);
            double result2 = Torus.torus_square_14c(function_3d_index, functions.function_3d, r1, r2);

            Console.WriteLine("  " + name
                                   + "  " + result1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + result2.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    public static void test48()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST48 tests TVEC_EVEN, TVEC_EVEN2 and TVEC_EVEN3.
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
        Console.WriteLine("");
        Console.WriteLine("TEST48");
        Console.WriteLine("  For evenly spaced angles between 0 and 2*PI:");
        Console.WriteLine("  TVEC_EVEN");
        Console.WriteLine("  TVEC_EVEN2");
        Console.WriteLine("  TVEC_EVEN3");

        int nt = 4;
        double[] t = typeMethods.tvec_even(nt);
        typeMethods.r8vec_print(nt, t, "  TVEC_EVEN:");

        nt = 4;
        t = typeMethods.tvec_even2(nt);
        typeMethods.r8vec_print(nt, t, "  TVEC_EVEN2:");

        nt = 4;
        t = typeMethods.tvec_even3(nt);
        typeMethods.r8vec_print(nt, t, "  TVEC_EVEN3:");

    }

    public static void test49()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST49 tests TVEC_EVEN_BRACKET, TVEC_EVEN_BRACKET2 and TVEC_EVEN_BRACKET3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("TEST49");
        Console.WriteLine("  For evenly spaced angles between THETA1 and THETA2:");
        Console.WriteLine("  TVEC_EVEN_BRACKET");
        Console.WriteLine("  TVEC_EVEN_BRACKET2.");
        Console.WriteLine("  TVEC_EVEN_BRACKET3.");

        const double theta1 = 30.0;
        const double theta2 = 90.0;

        Console.WriteLine("");
        Console.WriteLine("  THETA1 = " + theta1 + "");
        Console.WriteLine("  THETA2 = " + theta2 + "");

        int nt = 4;
        double[] t = typeMethods.tvec_even_bracket(nt, theta1, theta2);
        typeMethods.r8vec_print(nt, t, "  TVEC_EVEN_BRACKET");

        nt = 5;
        t = typeMethods.tvec_even_bracket2(nt, theta1, theta2);
        typeMethods.r8vec_print(nt, t, "  TVEC_EVEN_BRACKET2");

        nt = 3;
        t = typeMethods.tvec_even_bracket3(nt, theta1, theta2);
        typeMethods.r8vec_print(nt, t, "  TVEC_EVEN_BRACKET3");

    }
}