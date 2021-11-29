using System;
using System.Globalization;
using Burkardt.Types;

namespace Burkardt.LineNS;

public static class Felippa
{

    public static double line_monomial(double a, double b, int expon)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_MONOMIAL: monomial integral over a line segment in 1D.
        //
        //  Discussion:
        //
        //    This function returns the integral of X^EXPON.
        //
        //    The integration region is:
        //    A <= X <= B
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 September 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Arthur Stroud,
        //    Approximate Calculation of Multiple Integrals,
        //    Prentice Hall, 1971,
        //    ISBN: 0130438936,
        //    LC: QA311.S85.
        //
        //  Parameters:
        //
        //    Input, double A, B, the lower and upper limits.
        //
        //    Input, int EXPON, the exponent of X.  The exponent must not be -1.
        //
        //    Output, double LINE_MONOMIAL, the integral of X^EXPON.
        //
    {
        switch (expon)
        {
            case -1:
                Console.WriteLine("");
                Console.WriteLine("LINE_MONOMIAL - Fatal error!");
                Console.WriteLine("  Exponent = -1 is not a legal input.");
                return 1;
            default:
                double value = (Math.Pow(b, expon + 1) - Math.Pow(a, expon + 1))
                               / (expon + 1);

                return value;
        }
    }

    public static void line_monomial_test(int degree_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_MONOMIAL_TEST tests LINE_MONOMIAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 September 2014
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
        const double a = 0.0;
        const double b = 1.0;
        int expon;

        Console.WriteLine("");
        Console.WriteLine("LINE_MONOMIAL_TEST");
        Console.WriteLine("  For a line segment in 1D,");
        Console.WriteLine("  LINE_MONOMIAL returns the exact value of the");
        Console.WriteLine("  integral of X^EXPON");
        Console.WriteLine("");
        Console.WriteLine("  Volume = " + line_volume(a, b) + "");
        Console.WriteLine("");
        Console.WriteLine("     EXPON      INTEGRAL");
        Console.WriteLine("");

        for (expon = 0; expon <= degree_max; expon++)
        {
            double value = line_monomial(a, b, expon);

            Console.WriteLine("  " + expon.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + value.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    public static void line_quad_test(int degree_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_QUAD_TEST tests the rules for a line segment in 1D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 September 2014
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
        const double a = 0.0;
        const double b = 1.0;
        int expon;

        Console.WriteLine("");
        Console.WriteLine("LINE_QUAD_TEST");
        Console.WriteLine("  For a line segment in 1D,");
        Console.WriteLine("  we approximate monomial integrals with:");
        Console.WriteLine("  LINE_UNIT_O01, a 1 point rule.");
        Console.WriteLine("  LINE_UNIT_O02, a 2 point rule.");
        Console.WriteLine("  LINE_UNIT_O03, a 3 point rule.");
        Console.WriteLine("  LINE_UNIT_O04, a 4 point rule.");
        Console.WriteLine("  LINE_UNIT_O05, a 5 point rule.");

        for (expon = 0; expon <= degree_max; expon++)
        {
            Console.WriteLine("");
            Console.WriteLine("  Monomial exponent:   " + expon + "");
            Console.WriteLine("");

            int order;
            double quad;
            for (order = 1; order <= 5; order++)
            {
                double[] v = new double[order];
                double[] w = new double[order];
                double[] x = new double[order];

                line_rule(a, b, order, ref w, ref x);
                int j;
                for (j = 0; j < order; j++)
                {
                    v[j] = Math.Pow(x[j], expon);
                }

                quad = typeMethods.r8vec_dot_product(order, w, v);
                Console.WriteLine(order.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "  "
                                                              + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }

            Console.WriteLine("");
            quad = line_monomial(a, b, expon);
            Console.WriteLine("   Exact  " + quad.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    public static void line_rule(double a, double b, int order, ref double[] w, ref double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_RULE returns a quadrature rule for a line segment in 1D.
        //
        //  Discussion:
        //
        //    The integration region is:
        //      A <= X <= B
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 September 2014
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
        //    Input, double A, B, the lower and upper limits.
        //
        //    Input, int ORDER, the order of the rule.
        //
        //    Output, double W[ORDER], the weights.
        //
        //    Output, double X[ORDER], the abscissas.
        //
    {
        int j;

        switch (order)
        {
            case 1:
                line_unit_o01(ref w, ref x);
                break;
            case 2:
                line_unit_o02(ref w, ref x);
                break;
            case 3:
                line_unit_o03(ref w, ref x);
                break;
            case 4:
                line_unit_o04(ref w, ref x);
                break;
            case 5:
                line_unit_o05(ref w, ref x);
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("LINE_RULE - Fatal error!");
                Console.WriteLine("  Illegal value of ORDER.");
                return;
        }

        //
        //  Transform from [-1,+1] to [A,B]
        //
        for (j = 0; j < order; j++)
        {
            w[j] = w[j] * (b - a) / 2.0;
            x[j] = ((1.0 - x[j]) * a
                    + (1.0 + x[j]) * b)
                   / 2.0;
        }
    }

    public static void line_unit_o01(ref double[] w, ref double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_UNIT_O01 returns a 1 point quadrature rule for the unit line.
        //
        //  Discussion:
        //
        //    The integration region is:
        //
        //    - 1.0 <= X <= 1.0
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 April 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Arthur Stroud,
        //    Approximate Calculation of Multiple Integrals,
        //    Prentice Hall, 1971,
        //    ISBN: 0130438936,
        //    LC: QA311.S85.
        //
        //  Parameters:
        //
        //    Output, double W[1], the weights.
        //
        //    Output, double X[1], the abscissas.
        //
    {
        const int order = 1;
        double[] w_save = {2.0};
        double[] x_save = {0.0};

        typeMethods.r8vec_copy(order, w_save, ref w);
        typeMethods.r8vec_copy(order, x_save, ref x);
    }

    public static void line_unit_o02(ref double[] w, ref double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_UNIT_O02 returns a 2 point quadrature rule for the unit line.
        //
        //  Discussion:
        //
        //    The integration region is:
        //
        //    - 1.0 <= X <= 1.0
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 April 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Arthur Stroud,
        //    Approximate Calculation of Multiple Integrals,
        //    Prentice Hall, 1971,
        //    ISBN: 0130438936,
        //    LC: QA311.S85.
        //
        //  Parameters:
        //
        //    Output, double W[2], the weights.
        //
        //    Output, double X[2], the abscissas.
        //
    {
        const int order = 2;
        double[] w_save =
        {
            1.0000000000000000000,
            1.0000000000000000000
        };
        double[] x_save =
        {
            -0.57735026918962576451,
            0.57735026918962576451
        };

        typeMethods.r8vec_copy(order, w_save, ref w);
        typeMethods.r8vec_copy(order, x_save, ref x);
    }

    public static void line_unit_o03(ref double[] w, ref double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_UNIT_O03 returns a 3 point quadrature rule for the unit line.
        //
        //  Discussion:
        //
        //    The integration region is:
        //
        //    - 1.0 <= X <= 1.0
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 April 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Arthur Stroud,
        //    Approximate Calculation of Multiple Integrals,
        //    Prentice Hall, 1971,
        //    ISBN: 0130438936,
        //    LC: QA311.S85.
        //
        //  Parameters:
        //
        //    Output, double W[3], the weights.
        //
        //    Output, double X[3], the abscissas.
        //
    {
        const int order = 3;
        double[] w_save =
        {
            0.55555555555555555556,
            0.88888888888888888889,
            0.55555555555555555556
        };
        double[] x_save =
        {
            -0.77459666924148337704,
            0.00000000000000000000,
            0.77459666924148337704
        };

        typeMethods.r8vec_copy(order, w_save, ref w);
        typeMethods.r8vec_copy(order, x_save, ref x);
    }

    public static void line_unit_o04(ref double[] w, ref double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_UNIT_O04 returns a 4 point quadrature rule for the unit line.
        //
        //  Discussion:
        //
        //    The integration region is:
        //
        //    - 1.0 <= X <= 1.0
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 April 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Arthur Stroud,
        //    Approximate Calculation of Multiple Integrals,
        //    Prentice Hall, 1971,
        //    ISBN: 0130438936,
        //    LC: QA311.S85.
        //
        //  Parameters:
        //
        //    Output, double W[4], the weights.
        //
        //    Output, double X[4], the abscissas.
        //
    {
        const int order = 4;
        double[] w_save =
        {
            0.34785484513745385737,
            0.65214515486254614263,
            0.65214515486254614263,
            0.34785484513745385737
        };
        double[] x_save =
        {
            -0.86113631159405257522,
            -0.33998104358485626480,
            0.33998104358485626480,
            0.86113631159405257522
        };

        typeMethods.r8vec_copy(order, w_save, ref w);
        typeMethods.r8vec_copy(order, x_save, ref x);
    }

    public static void line_unit_o05(ref double[] w, ref double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_UNIT_O05 returns a 5 point quadrature rule for the unit line.
        //
        //  Discussion:
        //
        //    The integration region is:
        //
        //    - 1.0 <= X <= 1.0
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 April 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Arthur Stroud,
        //    Approximate Calculation of Multiple Integrals,
        //    Prentice Hall, 1971,
        //    ISBN: 0130438936,
        //    LC: QA311.S85.
        //
        //  Parameters:
        //
        //    Output, double W[5], the weights.
        //
        //    Output, double X[5], the abscissas.
        //
    {
        const int order = 5;
        double[] w_save =
        {
            0.23692688505618908751,
            0.47862867049936646804,
            0.56888888888888888889,
            0.47862867049936646804,
            0.23692688505618908751
        };
        double[] x_save =
        {
            -0.90617984593866399280,
            -0.53846931010568309104,
            0.00000000000000000000,
            0.53846931010568309104,
            0.90617984593866399280
        };

        typeMethods.r8vec_copy(order, w_save, ref w);
        typeMethods.r8vec_copy(order, x_save, ref x);
    }

    public static double line_volume(double a, double b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_VOLUME: volume of a line segment in 1D.
        //
        //  Discussion:
        //
        //    The integration region is:
        //    A <= X <= B
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 September 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, the lower and upper limits.
        //
        //    Output, double LINE_VOLUME, the volume of the line.
        //
    {
        double volume = b - a;

        return volume;
    }
}