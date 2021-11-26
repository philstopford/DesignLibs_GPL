﻿using System;
using System.Globalization;
using Burkardt.Composition;
using Burkardt.TriangleNS;
using Burkardt.Types;
using Monomial = Burkardt.MonomialNS.Monomial;

namespace TriangleFelippaRuleTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for TRIANGLE_FELIPPA_RULE_TEST.
        //
        //  Discussion:
        //
        //    TRIANGLE_FELIPPA_RULE_TEST tests the TRIANGLE_FELIPPA_RULE library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    25 August 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("TRIANGLE_FELIPPA_RULE_TEST");
        Console.WriteLine("  C++ version");
        Console.WriteLine("  Test the TRIANGLE_FELIPPA_RULE library.");

        int degree_max = 4;
        triangle_unit_monomial_test(degree_max);

        degree_max = 7;
        triangle_unit_quad_test(degree_max);
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("TRIANGLE_FELIPPA_RULE_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");

    }

    private static void triangle_unit_monomial_test(int degree_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_UNIT_MONOMIAL_TEST tests TRIANGLE_UNIT_MONOMIAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 April 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DEGREE_MAX, the maximum total degree of the
        //    monomials to check.
        //
    {
        int alpha;
        int beta;
        int[] expon = new int[2];

        Console.WriteLine("");
        Console.WriteLine("TRIANGLE_UNIT_MONOMIAL_TEST");
        Console.WriteLine("  For the unit triangle,");
        Console.WriteLine("  TRIANGLE_UNIT_MONOMIAL returns the exact value of the");
        Console.WriteLine("  integral of X^ALPHA Y^BETA");
        Console.WriteLine("");
        Console.WriteLine("  Volume = " + QuadratureRule.triangle_unit_volume() + "");
        Console.WriteLine("");
        Console.WriteLine("     ALPHA      BETA      INTEGRAL");
        Console.WriteLine("");

        for (alpha = 0; alpha <= degree_max; alpha++)
        {
            expon[0] = alpha;
            for (beta = 0; beta <= degree_max - alpha; beta++)
            {
                expon[1] = beta;

                double value = QuadratureRule.triangle_unit_monomial(expon);

                Console.WriteLine("  " + expon[0].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + expon[1].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + value.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }
        }

    }

    private static void triangle_unit_quad_test(int degree_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_UNIT_QUAD_TEST tests the rules for the unit triangle.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 April 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DEGREE_MAX, the maximum total degree of the
        //    monomials to check.
        //
    {
        const int DIM_NUM = 2;

        int[] expon = new int[DIM_NUM];
        int h = 0;
        int t = 0;
        SubCompData data = new();

        Console.WriteLine("");
        Console.WriteLine("TRIANGLE_UNIT_QUAD_TEST");
        Console.WriteLine("  For the unit triangle,");
        Console.WriteLine("  we approximate monomial integrals with:");
        Console.WriteLine("  QuadratureRule.triangle_unit_o01,");
        Console.WriteLine("  QuadratureRule.triangle_unit_o03,");
        Console.WriteLine("  QuadratureRule.triangle_unit_o03b,");
        Console.WriteLine("  QuadratureRule.triangle_unit_o06,");
        Console.WriteLine("  QuadratureRule.triangle_unit_o06b,");
        Console.WriteLine("  QuadratureRule.triangle_unit_o07,");
        Console.WriteLine("  QuadratureRule.triangle_unit_o12,");

        bool more = false;

        for (;;)
        {
            SubComp.subcomp_next(ref data, degree_max, DIM_NUM, ref expon, ref more, ref h, ref t);

            Console.WriteLine("");
            string cout = "  Monomial exponents: ";
            int dim;
            for (dim = 0; dim < DIM_NUM; dim++)
            {
                cout += "  " + expon[dim].ToString(CultureInfo.InvariantCulture).PadLeft(2);
            }

            Console.WriteLine(cout);
            Console.WriteLine("");

            int order = 1;
            double[] w = new double[order];
            double[] xy = new double[DIM_NUM * order];
            QuadratureRule.triangle_unit_o01(ref w, ref xy);
            double[] v = Monomial.monomial_value(DIM_NUM, order, expon, xy);
            double quad = QuadratureRule.triangle_unit_volume() * typeMethods.r8vec_dot_product(order, w, v);
            Console.WriteLine("  " + order.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

            order = 3;
            w = new double[order];
            xy = new double[DIM_NUM * order];
            QuadratureRule.triangle_unit_o03(ref w, ref xy);
            v = Monomial.monomial_value(DIM_NUM, order, expon, xy);
            quad = QuadratureRule.triangle_unit_volume() * typeMethods.r8vec_dot_product(order, w, v);
            Console.WriteLine("  " + order.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

            order = 3;
            w = new double[order];
            xy = new double[DIM_NUM * order];
            QuadratureRule.triangle_unit_o03b(ref w, ref xy);
            v = Monomial.monomial_value(DIM_NUM, order, expon, xy);
            quad = QuadratureRule.triangle_unit_volume() * typeMethods.r8vec_dot_product(order, w, v);
            Console.WriteLine("  " + order.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

            order = 6;
            w = new double[order];
            xy = new double[DIM_NUM * order];
            QuadratureRule.triangle_unit_o06(ref w, ref xy);
            v = Monomial.monomial_value(DIM_NUM, order, expon, xy);
            quad = QuadratureRule.triangle_unit_volume() * typeMethods.r8vec_dot_product(order, w, v);
            Console.WriteLine("  " + order.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

            order = 6;
            w = new double[order];
            xy = new double[DIM_NUM * order];
            QuadratureRule.triangle_unit_o06b(ref w, ref xy);
            v = Monomial.monomial_value(DIM_NUM, order, expon, xy);
            quad = QuadratureRule.triangle_unit_volume() * typeMethods.r8vec_dot_product(order, w, v);
            Console.WriteLine("  " + order.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

            order = 7;
            w = new double[order];
            xy = new double[DIM_NUM * order];
            QuadratureRule.triangle_unit_o07(ref w, ref xy);
            v = Monomial.monomial_value(DIM_NUM, order, expon, xy);
            quad = QuadratureRule.triangle_unit_volume() * typeMethods.r8vec_dot_product(order, w, v);
            Console.WriteLine("  " + order.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

            order = 12;
            w = new double[order];
            xy = new double[DIM_NUM * order];
            QuadratureRule.triangle_unit_o12(ref w, ref xy);
            v = Monomial.monomial_value(DIM_NUM, order, expon, xy);
            quad = QuadratureRule.triangle_unit_volume() * typeMethods.r8vec_dot_product(order, w, v);
            Console.WriteLine("  " + order.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

            Console.WriteLine("");
            quad = QuadratureRule.triangle_unit_monomial(expon);
            Console.WriteLine("  " + " Exact"
                                   + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

            if (!more)
            {
                break;
            }
        }

    }
}