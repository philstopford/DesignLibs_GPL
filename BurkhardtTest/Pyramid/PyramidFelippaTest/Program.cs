using System;
using System.Globalization;
using Burkardt.Composition;
using Burkardt.MonomialNS;
using Burkardt.PyramidNS;
using Burkardt.Types;

namespace PyramidFelippaTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for PYRAMID_FELIPPA_RULE_TEST.
        //
        //  Discussion:
        //
        //    PYRAMID_FELIPPA_RULE_TEST tests the PYRAMID_FELIPPA_RULE library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 August 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("PYRAMID_FELIPPA_RULE_TEST");
        Console.WriteLine("  Test the PYRAMID_FELIPPA_RULE library.");

        int degree_max = 4;
        pyramid_unit_monomial_test(degree_max);

        degree_max = 5;
        pyramid_unit_quad_test(degree_max);

        Console.WriteLine("");
        Console.WriteLine("PYRAMID_FELIPPA_RULE_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void pyramid_unit_monomial_test(int degree_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PYRAMID__UNIT_MONOMIAL_TEST tests PYRAMID__UNIT_MONOMIAL.
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
        int[] expon = new int[3];

        Console.WriteLine("");
        Console.WriteLine("PYRAMID__UNIT_MONOMIAL_TEST");
        Console.WriteLine("  For the unit pyramid,");
        Console.WriteLine("  PYRAMID__UNIT_MONOMIAL returns the exact value of the");
        Console.WriteLine("  integral of X^ALPHA Y^BETA Z^GAMMA");
        Console.WriteLine("");
        Console.WriteLine("  Volume = " + FelippaRule.pyramid_unit_volume() + "");
        Console.WriteLine("");
        Console.WriteLine("     ALPHA      BETA     GAMMA      INTEGRAL");
        Console.WriteLine("");

        for (alpha = 0; alpha <= degree_max; alpha++)
        {
            expon[0] = alpha;
            int beta;
            for (beta = 0; beta <= degree_max - alpha; beta++)
            {
                expon[1] = beta;
                int gamma;
                for (gamma = 0; gamma <= degree_max - alpha - beta; gamma++)
                {
                    expon[2] = gamma;

                    double value = FelippaRule.pyramid_unit_monomial(expon);

                    Console.WriteLine("  " + expon[0].ToString().PadLeft(8)
                                           + "  " + expon[1].ToString().PadLeft(8)
                                           + "  " + expon[2].ToString().PadLeft(8)
                                           + "  " + value.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
                }
            }
        }
    }

    private static void pyramid_unit_quad_test(int degree_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PYRAMID__UNIT_QUAD_TEST tests the rules for the unit pyramid.
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
        //  Parameters:
        //
        //    Input, int DEGREE_MAX, the maximum total degree of the
        //    monomials to check.
        //
    {
        const int DIM_NUM = 3;

        int[] expon = new int[DIM_NUM];
        int h = 0;
        int t = 0;

        Console.WriteLine("");
        Console.WriteLine("PYRAMID__UNIT_QUAD_TEST");
        Console.WriteLine("  For the unit pyramid,");
        Console.WriteLine("  we approximate monomial integrals with:");
        Console.WriteLine("  PYRAMID__UNIT_O01,");
        Console.WriteLine("  PYRAMID__UNIT_O05,");
        Console.WriteLine("  PYRAMID__UNIT_O06,");
        Console.WriteLine("  PYRAMID__UNIT_O08,");
        Console.WriteLine("  PYRAMID__UNIT_O08b,");
        Console.WriteLine("  PYRAMID__UNIT_O09,");
        Console.WriteLine("  PYRAMID__UNIT_O13,");
        Console.WriteLine("  PYRAMID__UNIT_O18,");
        Console.WriteLine("  PYRAMID__UNIT_O27,");
        Console.WriteLine("  PYRAMID__UNIT_O48,");

        bool more = false;

        SubCompData data = new();

        for (;;)
        {
            SubComp.subcomp_next(ref data, degree_max, DIM_NUM, ref expon, ref more, ref h, ref t);

            if (expon[0] % 2 == 1 || expon[1] % 2 == 1)
            {
                continue;
            }

            Console.WriteLine("");
            string cout = "  Monomial exponents: ";
            int dim;
            for (dim = 0; dim < DIM_NUM; dim++)
            {
                cout += "  " + expon[dim].ToString().PadLeft(2);
            }

            Console.WriteLine(cout);
            Console.WriteLine("");

            int order = 1;
            double[] w = new double[order];
            double[] xyz = new double[DIM_NUM * order];
            FelippaRule.pyramid_unit_o01(ref w, ref xyz);
            double[] v = Monomial.monomial_value(DIM_NUM, order, expon, xyz);
            double quad = FelippaRule.pyramid_unit_volume() * typeMethods.r8vec_dot_product(order, w, v);
            Console.WriteLine("  " + order.ToString().PadLeft(6)
                                   + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

            order = 5;
            w = new double[order];
            xyz = new double[DIM_NUM * order];
            FelippaRule.pyramid_unit_o05(ref w, ref xyz);
            v = Monomial.monomial_value(DIM_NUM, order, expon, xyz);
            quad = FelippaRule.pyramid_unit_volume() * typeMethods.r8vec_dot_product(order, w, v);
            Console.WriteLine("  " + order.ToString().PadLeft(6)
                                   + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

            order = 6;
            w = new double[order];
            xyz = new double[DIM_NUM * order];
            FelippaRule.pyramid_unit_o06(ref w, ref xyz);
            v = Monomial.monomial_value(DIM_NUM, order, expon, xyz);
            quad = FelippaRule.pyramid_unit_volume() * typeMethods.r8vec_dot_product(order, w, v);
            Console.WriteLine("  " + order.ToString().PadLeft(6)
                                   + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

            order = 8;
            w = new double[order];
            xyz = new double[DIM_NUM * order];
            FelippaRule.pyramid_unit_o08(ref w, ref xyz);
            v = Monomial.monomial_value(DIM_NUM, order, expon, xyz);
            quad = FelippaRule.pyramid_unit_volume() * typeMethods.r8vec_dot_product(order, w, v);
            Console.WriteLine("  " + order.ToString().PadLeft(6)
                                   + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

            order = 8;
            w = new double[order];
            xyz = new double[DIM_NUM * order];
            FelippaRule.pyramid_unit_o08b(ref w, ref xyz);
            v = Monomial.monomial_value(DIM_NUM, order, expon, xyz);
            quad = FelippaRule.pyramid_unit_volume() * typeMethods.r8vec_dot_product(order, w, v);
            Console.WriteLine("  " + order.ToString().PadLeft(6)
                                   + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

            order = 9;
            w = new double[order];
            xyz = new double[DIM_NUM * order];
            FelippaRule.pyramid_unit_o09(ref w, ref xyz);
            v = Monomial.monomial_value(DIM_NUM, order, expon, xyz);
            quad = FelippaRule.pyramid_unit_volume() * typeMethods.r8vec_dot_product(order, w, v);
            Console.WriteLine("  " + order.ToString().PadLeft(6)
                                   + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

            order = 13;
            w = new double[order];
            xyz = new double[DIM_NUM * order];
            FelippaRule.pyramid_unit_o13(ref w, ref xyz);
            v = Monomial.monomial_value(DIM_NUM, order, expon, xyz);
            quad = FelippaRule.pyramid_unit_volume() * typeMethods.r8vec_dot_product(order, w, v);
            Console.WriteLine("  " + order.ToString().PadLeft(6)
                                   + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

            order = 18;
            w = new double[order];
            xyz = new double[DIM_NUM * order];
            FelippaRule.pyramid_unit_o18(ref w, ref xyz);
            v = Monomial.monomial_value(DIM_NUM, order, expon, xyz);
            quad = FelippaRule.pyramid_unit_volume() * typeMethods.r8vec_dot_product(order, w, v);
            Console.WriteLine("  " + order.ToString().PadLeft(6)
                                   + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

            order = 27;
            w = new double[order];
            xyz = new double[DIM_NUM * order];
            FelippaRule.pyramid_unit_o27(ref w, ref xyz);
            v = Monomial.monomial_value(DIM_NUM, order, expon, xyz);
            quad = FelippaRule.pyramid_unit_volume() * typeMethods.r8vec_dot_product(order, w, v);
            Console.WriteLine("  " + order.ToString().PadLeft(6)
                                   + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

            order = 48;
            w = new double[order];
            xyz = new double[DIM_NUM * order];
            FelippaRule.pyramid_unit_o48(ref w, ref xyz);
            v = Monomial.monomial_value(DIM_NUM, order, expon, xyz);
            quad = FelippaRule.pyramid_unit_volume() * typeMethods.r8vec_dot_product(order, w, v);
            Console.WriteLine("  " + order.ToString().PadLeft(6)
                                   + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

            Console.WriteLine("");
            quad = FelippaRule.pyramid_unit_monomial(expon);
            Console.WriteLine("  " + " Exact"
                                   + "  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

            if (!more)
            {
                break;
            }
        }
    }
}