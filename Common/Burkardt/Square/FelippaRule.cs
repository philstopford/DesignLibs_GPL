using System;
using Burkardt.Composition;
using Burkardt.MonomialNS;
using Burkardt.Types;

namespace Burkardt.Square;

public static class FelippaRule
{

    public static double square_monomial(double[] a, double[] b, int[] expon)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SQUARE_MONOMIAL integrates a monomial over a square in 2D.
        //
        //  Discussion:
        //
        //    This routine integrates a monomial of the form
        //
        //      product ( 1 <= dim <= 2 ) x(dim)^expon(dim)
        //
        //    where the exponents are nonnegative integers.  Note that
        //    if the combination 0^0 is encountered, it should be treated
        //    as 1.
        //
        //    The integration region is:
        //      A(1) <= X <= B(1)
        //      A(2) <= Y <= B(2)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 September 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A[2], B[2], the lower and upper limits.
        //
        //    Input, int EXPON[2], the exponents.
        //
        //    Output, double SQUARE_MONOMIAL, the integral of the monomial.
        //
    {
        int i;
        double value = 0;

        for (i = 0; i < 2; i++)
        {
            switch (expon[i])
            {
                case -1:
                    Console.WriteLine("");
                    Console.WriteLine("SQUARE_MONOMIAL - Fatal error!");
                    Console.WriteLine("  Exponent of -1 encountered.");
                    return 1;
            }
        }

        value = 1.0;

        for (i = 0; i < 2; i++)
        {
            value = value * (Math.Pow(b[i], expon[i] + 1) - Math.Pow(a[i], expon[i] + 1))
                    / (expon[i] + 1);
        }

        return value;
    }

    public static void square_monomial_test(int degree_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SQUARE_MONOMIAL_TEST tests SQUARE_MONOMIAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 September 2014
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
        double[] a = { -1.0, -1.0 };
        int alpha;
        double[] b = { +1.0, +1.0 };
        int beta;
        int[] expon = new int [2];
        double value = 0;

        Console.WriteLine("");
        Console.WriteLine("SQUARE_MONOMIAL_TEST");
        Console.WriteLine("  For a square in 2D,");
        Console.WriteLine("  SQUARE_MONOMIAL returns the exact value of the");
        Console.WriteLine("  integral of X^ALPHA Y^BETA");
        Console.WriteLine("");
        Console.WriteLine("  Volume = " + square_volume(a, b) + "");
        Console.WriteLine("");
        Console.WriteLine("     ALPHA      BETA      INTEGRAL");
        Console.WriteLine("");

        for (alpha = 0; alpha <= degree_max; alpha++)
        {
            expon[0] = alpha;
            for (beta = 0; beta <= degree_max - alpha; beta++)
            {
                expon[1] = beta;

                value = square_monomial(a, b, expon);

                Console.WriteLine("  " + expon[0].ToString().PadLeft(8)
                                       + "  " + expon[1].ToString().PadLeft(8)
                                       + "  " + value.ToString().PadLeft(14) + "");
            }
        }
    }

    public static void square_quad_test(int degree_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SQUARE_QUAD_TEST tests the rules for a square in 2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 September 2014
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
        double[] a = { -1.0, -1.0 };
        double[] b = { +1.0, +1.0 };
        int[] expon = new int[2];
        int h = 0;
        int i;
        int k;
        bool more;
        int order;
        int[] order_1d = new int[2];
        double quad;
        int t = 0;
        double[] v;
        double[] w;
        double[] xy;
        SubCompData data = new();

        Console.WriteLine("");
        Console.WriteLine("SQUARE_QUAD_TEST");
        Console.WriteLine("  For a square in 2D,");
        Console.WriteLine("  we approximate monomial integrals with:");
        Console.WriteLine("  SQUARE_RULE, which returns M by N point rules..");

        more = false;

        for (;;)
        {
            SubComp.subcomp_next(ref data, degree_max, 2, ref expon, ref more, ref h, ref t);

            Console.WriteLine("");
            string cout = "  Monomial exponents: ";
            for (i = 0; i < 2; i++)
            {
                cout += "  " + expon[i].ToString().PadLeft(2);
            }

            Console.WriteLine(cout);
            Console.WriteLine("");

            for (k = 1; k <= 5; k++)
            {
                for (i = 0; i < 2; i++)
                {
                    order_1d[i] = k;
                }

                order = order_1d[0] * order_1d[1];
                w = new double[order];
                xy = new double[2 * order];
                square_rule(a, b, order_1d, ref w, ref xy);
                v = Monomial.monomial_value(2, order, expon, xy);
                quad = typeMethods.r8vec_dot_product(order, w, v);
                Console.WriteLine("  " + order_1d[0].ToString().PadLeft(6)
                                       + "  " + order_1d[1].ToString().PadLeft(6)
                                       + "  " + quad.ToString().PadLeft(14) + "");
            }

            //
            //  Try a rule of mixed orders.
            //
            order_1d[0] = 3;
            order_1d[1] = 5;
            order = order_1d[0] * order_1d[1];
            w = new double[order];
            xy = new double[2 * order];
            square_rule(a, b, order_1d, ref w, ref xy);
            v = Monomial.monomial_value(2, order, expon, xy);
            quad = typeMethods.r8vec_dot_product(order, w, v);
            Console.WriteLine("  " + order_1d[0].ToString().PadLeft(6)
                                   + "  " + order_1d[1].ToString().PadLeft(6)
                                   + "  " + quad.ToString().PadLeft(14) + "");

            Console.WriteLine("");
            quad = square_monomial(a, b, expon);
            Console.WriteLine("  " + " Exact"
                                   + "  " + "      "
                                   + "  " + quad.ToString().PadLeft(14) + "");

            if (!more)
            {
                break;
            }
        }
    }

    public static void square_rule(double[] a, double[] b, int[] order_1d, ref double[] w,
            ref double[] xy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SQUARE_RULE returns a quadrature rule for a square in 2D.
        //
        //  Discussion:
        //
        //    The integration region is:
        //      A(1) <= X <= B(1)
        //      A(2) <= Y <= B(2)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 September 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Carlos Felippa,
        //    A compendium of FEM integration formulas for symbolic work,
        //    Engineering Computation,
        //    Volume 21, Number 8, 2004, pages 867-890.
        //
        //  Parameters:
        //
        //    Input, double A[2], B[2], the lower and upper limits.
        //
        //    Input, int ORDER_1D[2], the order of the rule in 
        //    each dimension.  1 <= ORDER_1D(I) <= 5.
        //
        //    Output, double W[ORDER_1D[0]*ORDER_1D[1]], the weights.
        //
        //    Output, double XY[2*ORDER_1D[0]*ORDER_1D[1]], the abscissas.
        //
    {
        int i;
        int j;
        int o;
        int order;
        double[] w_1d;
        double[] x_1d;
        typeMethods.r8vecDPData data1 = new();
        typeMethods.r8vecDPData data2 = new();

        order = order_1d[0] * order_1d[1];

        for (i = 0; i < 2; i++)
        {
            o = order_1d[i];

            w_1d = new double[o];
            x_1d = new double[o];

            switch (o)
            {
                case 1:
                    LineNS.Felippa.line_unit_o01(ref w_1d, ref x_1d);
                    break;
                case 2:
                    LineNS.Felippa.line_unit_o02(ref w_1d, ref x_1d);
                    break;
                case 3:
                    LineNS.Felippa.line_unit_o03(ref w_1d, ref x_1d);
                    break;
                case 4:
                    LineNS.Felippa.line_unit_o04(ref w_1d, ref x_1d);
                    break;
                case 5:
                    LineNS.Felippa.line_unit_o05(ref w_1d, ref x_1d);
                    break;
                default:
                    Console.WriteLine("");
                    Console.WriteLine("SQUARE_RULE - Fatal error!");
                    Console.WriteLine("  Illegal value of ORDER_1D[*].");
                    return;
            }

            //
            //  Transform from [-1,+1] to [Ai,Bi]
            //
            for (j = 0; j < o; j++)
            {
                w_1d[j] = w_1d[j] * (b[i] - a[i]) / 2.0;
                x_1d[j] = ((1.0 - x_1d[j]) * a[i]
                           + (1.0 + x_1d[j]) * b[i])
                          / 2.0;
            }

            //
            //  Add this information to the rule.
            //
            typeMethods.r8vec_direct_product(ref data1, i, o, x_1d, 2, order, ref xy);

            typeMethods.r8vec_direct_product2(ref data2, i, o, w_1d, 2, order, ref w);
        }

    }

    public static double square_volume(double[] a, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SQUARE_VOLUME: volume of a unit quadrilateral.
        //
        //  Discussion:
        //
        //    The integration region is:
        //      A(1) <= X <= B(1)
        //      A(2) <= Y <= B(2)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 September 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A[2], B[2], the lower and upper limits.
        //
        //    Output, double SQUARE_VOLUME, the volume.
        //
    {
        double value = 0;

        value = (b[0] - a[0]) * (b[1] - a[1]);

        return value;
    }
}