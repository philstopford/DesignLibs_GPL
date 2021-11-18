using System;
using Burkardt.Graph;
using Burkardt.PolynomialNS;
using Burkardt.Types;

using Integrals = Burkardt.TriangleNS.Integrals;

namespace TriangleIntegralsTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for TRIANGLE_INTEGRALS_TEST.
        //
        //  Discussion:
        //
        //    TRIANGLE_INTEGRALS_TEST tests the TRIANGLE_INTEGRALS library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 April 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("TRIANGLE_INTEGRALS_TEST:");
            
        Console.WriteLine("  Test the TRIANGLE_INTEGRALS library.");

        i4_to_pascal_test();
        i4_to_pascal_degree_test();
        pascal_to_i4_test();
        r8mat_print_test();
        r8mat_print_some_test();
        trinomial_test();

        rs_to_xy_map_test();
        xy_to_rs_map_test();

        poly_print_test();
        poly_power_linear_test();
        poly_power_test();
        poly_product_test();

        triangle01_monomial_integral_test();
        triangle01_poly_integral_test();
        triangle_area_test();
        triangle_xy_integral_test();
        triangle_monomial_integral_test();
        triangle_poly_integral_test();

        Console.WriteLine("");
        Console.WriteLine("TRIANGLE_INTEGRALS_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void i4_to_pascal_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_TO_PASCAL_TEST tests I4_TO_PASCAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 April 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i = 0;
        int j = 0;
        int k;

        Console.WriteLine("");
        Console.WriteLine("I4_TO_PASCAL_TEST");
        Console.WriteLine("  I4_TO_PASCAL converts a linear index to");
        Console.WriteLine("  Pascal triangle indices.");
        Console.WriteLine("");
        Console.WriteLine("     K  =>   I     J");
        Console.WriteLine("");

        for (k = 1; k <= 20; k++)
        {
            typeMethods.i4_to_pascal(k, ref i, ref j);
            Console.WriteLine("  " + k.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "    " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + j.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "");
        }
    }

    private static void i4_to_pascal_degree_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    I4_TO_PASCAL_DEGREE_TEST tests I4_TO_PASCAL_DEGREE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 April 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int d;
        int k;

        Console.WriteLine("");
        Console.WriteLine("I4_TO_PASCAL_DEGREE_TEST");
        Console.WriteLine("  I4_TO_PASCAL_DEGREE converts a linear index to");
        Console.WriteLine("  the degree of the corresponding Pascal triangle indices.");
        Console.WriteLine("");
        Console.WriteLine("     K  =>   D");
        Console.WriteLine("");

        for (k = 1; k <= 20; k++)
        {
            d = typeMethods.i4_to_pascal_degree(k);
            Console.WriteLine("  " + k.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "    " + d.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "");
        }
    }

    private static void pascal_to_i4_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PASCAL_TO_I4_TEST tests PASCAL_TO_I4.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 April 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int d;
        int i;
        int j;
        int k;

        Console.WriteLine("");
        Console.WriteLine("PASCAL_TO_I4_TEST");
        Console.WriteLine("  PASCAL_TO_I4 converts Pascal triangle indices to a");
        Console.WriteLine("  linear index.");
        Console.WriteLine("");
        Console.WriteLine("     I     J =>    K");
        Console.WriteLine("");

        for (d = 0; d <= 4; d++)
        {
            for (i = d; 0 <= i; i--)
            {
                j = d - i;
                k = typeMethods.pascal_to_i4(i, j);
                Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                       + "  " + j.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                       + "    " + k.ToString(CultureInfo.InvariantCulture).PadLeft(4) + "");
            }

            Console.WriteLine("");
        }
    }

    private static void poly_power_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLY_POWER_TEST tests POLY_POWER.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 April 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n1 = 2;
        int n4 = 3;

        int d1 = 1;
        int d2;
        int d3 = 2;
        int d4 = 2;
        int d5;
        int d6 = 6;

        double[] p1 =  {
                1.0, 2.0, 3.0
            }
            ;
        double[] p2;
        double[] p3 =  {
                1.0, 4.0, 6.0, 4.0, 12.0, 9.0
            }
            ;
        double[] p4 =  {
                1.0, -2.0, 3.0, -4.0, +5.0, -6.0
            }
            ;
        double[] p5;
        double[] p6 =  {
                1.0,
                -6.0, 9.0,
                0.0, -21.0, 9.0,
                40.0, -96.0, 108.0, -81.0,
                0.0, 84.0, -141.0, 171.0, -54.0,
                -96.0, 384.0, -798.0, 1017.0, -756.0, 324.0,
                -64.0, 240.0, -588.0, 845.0, -882.0, 540.0, -216.0
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("POLY_POWER_TEST:");
        Console.WriteLine("  POLY_POWER computes the N-th power of an X,Y polynomial.");
        //
        //  P1 = ( 1 + 2 x + 3 y )
        //  P2 = P1^2 = 1 + 4x + 6y + 4x^2 + 12xy + 9y^2 
        //  P3 = correct value
        //
        Console.WriteLine("");
        typeMethods.poly_print(d1, ref p1, "  p1(x,y)");

        d2 = n1 * d1;
        p2 = typeMethods.poly_power(d1, p1, n1);
        Console.WriteLine("");
        typeMethods.poly_print(d2, ref p2, "  p2(x,y) = p1(x,y)^2");

        Console.WriteLine("");
        typeMethods.poly_print(d3, ref p3, "  p3(x,y)=correct answer");
        //
        //  P4 = ( 1 - 2 x + 3 y - 4 x^2 + 5 xy - 6 y^2 )
        //  P5 = P4^3 =
        //    1
        //    -6x +9y
        //    +0x^2 - 21xy + 9y^2
        //    +40x^3 - 96x^2y  + 108x^y2 - 81y^3
        //    +0x^4 + 84x^3y - 141 x^2y^2 +171xy^3 - 54y^4
        //    -96x^5 + 384x^4y -798x^3y^2 + 1017 x^2y^3 - 756 xy^4 + 324 y^5
        //    -64x^6 + 240x^5y - 588x^4y^2 + 845 x^3y^3 - 882 x^2y^4 +540 xy^5 - 216y^6
        //
        Console.WriteLine("");
        typeMethods.poly_print(d4, ref p4, "  p4(x,y)");

        d5 = n4 * d4;
        p5 = typeMethods.poly_power(d4, p4, n4);
        Console.WriteLine("");
        typeMethods.poly_print(d5, ref p5, "  p5(x,y) = p1(x,y)^3");

        Console.WriteLine("");
        typeMethods.poly_print(d6, ref p6, "  p6(x,y)=correct answer");
    }

    private static void poly_power_linear_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLY_POWER_LINEAR_TEST tests POLY_POWER_LINEAR.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 April 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int n1 = 2;
        int n4 = 3;

        int d1 = 1;
        int d2 = d1 * n1;
        int d3 = 2;
        int d4 = 1;
        int d5 = d4 * n4;
        int d6 = 3;

        double[] p1 =  {
                1.0, 2.0, 3.0
            }
            ;
        double[] p2;
        double[] p3 =  {
                1.0, 4.0, 6.0, 4.0, 12.0, 9.0
            }
            ;
        double[] p4 =  {
                2.0, -1.0, 3.0
            }
            ;
        double[] p5;
        double[] p6 =  {
                8.0, -12.0, 36.0, 6.0, -36.0, 54.0, -1.0, 9.0, -27.0, 27.0
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("POLY_POWER_LINEAR_TEST:");
        Console.WriteLine("  POLY_POWER_LINEAR computes the N-th power of");
        Console.WriteLine("  a linear polynomial in X and Y.");
        //
        //  P1 = ( 1 + 2 x + 3 y )
        //  P2 = P1^2
        //  P3 = correct value
        //
        Console.WriteLine("");
        typeMethods.poly_print(d1, ref p1, "  p1(x,y)");

        p2 = typeMethods.poly_power_linear(d1, p1, n1);
        Console.WriteLine("");
        typeMethods.poly_print(d2, ref p2, "  p2(x,y) = p1(x,y)^n");

        Console.WriteLine("");
        typeMethods.poly_print(d3, ref p3, "  Correct answer");

        //
        //  P4 = ( 2 - x + 3 y )
        //  P5 = P4^3
        //  P6 = correct value
        //
        Console.WriteLine("");
        typeMethods.poly_print(d4, ref p4, "  p4(x,y)");

        p5 = typeMethods.poly_power_linear(d4, p4, n4);
        Console.WriteLine("");
        typeMethods.poly_print(d5, ref p5, "  p5(x,y) = p4(x,y)^3");

        Console.WriteLine("");
        typeMethods.poly_print(d6, ref p6, "  Correct answer");
    }

    private static void poly_print_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLY_PRINT_TEST tests POLY_PRINT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 April 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int d1 = 0;
        int d2 = 1;
        int d3 = 2;
        int d4 = 3;

        double[] p1 =  {
                12.34
            }
            ;
        double[] p2 =  {
                1.0, 2.0, 3.0
            }
            ;
        double[] p3 =  {
                0.0, 0.0, 0.0, 0.0, 1.0, 0.0
            }
            ;
        double[] p4 =  {
                1.0, -2.1, +3.2, -4.3, +5.4,
                -6.5, +7.6, -8.7, +9.8, -10.9
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("POLY_PRINT_TEST:");
        Console.WriteLine("  POLY_PRINT can print a D-degree polynomial in X and Y.");
        //
        //  P1 = 12.34
        //
        Console.WriteLine("");
        Console.WriteLine("  P1(x,y) = 12.34");
        typeMethods.poly_print(d1, ref p1, "  p1(x,y)");
        //
        //  P2 = 1.0 + 2.0 * x + 3.0 * Y
        //
        Console.WriteLine("");
        Console.WriteLine("  P2(x,y) = 1 + 2 * x + 3 * Y");
        typeMethods.poly_print(d2, ref p2, "  p2(x,y)");
        //
        //  P3 = XY
        //
        Console.WriteLine("");
        Console.WriteLine("  P3(x,y) = xy");
        typeMethods.poly_print(d3, ref p3, "  p3(x,y) = xy");
        //
        //  P4 = 1 - 2.1 * x + 3.2 * y - 4.3 * x^2 + 5.4 * xy - 6.5 * y^2
        //    + 7.6 * x^3 - 8.7 * x^2y + 9.8 * xy^2 - 10.9 * y^3.
        //
        Console.WriteLine("");
        Console.WriteLine("  P4(x,y) = 1.0 - 2.1 * x + 3.2 * y - 4.3 * x^2 ");
        Console.WriteLine("          + 5.4 * xy - 6.5 * y^2 + 7.6 * x^3 ");
        Console.WriteLine("          - 8.7 * x^2y + 9.8 * xy^2 - 10.9 * y^3.");
        typeMethods.poly_print(d4, ref p4, "  p4(x,y)");
    }

    private static void poly_product_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLY_PRODUCT_TEST tests POLY_PRODUCT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 April 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int d1 = 1;
        int d2 = 1;
        int d3;
        int d4 = d1 + d2;

        int d5 = 2;
        int d6 = 2;
        int d7;
        int d8 = d5 + d6;

        double[] p1 =  {
                1.0, 2.0, 3.0
            }
            ;
        double[] p2 =  {
                4.0, 5.0, 0.0
            }
            ;
        double[] p3;
        double[] p4 =  {
                4.0, 13.0, 12.0, 10.0, 15.0, 0.0
            }
            ;
        double[] p5 =  {
                1.0, -2.0, 3.0, -4.0, +5.0, -6.0
            }
            ;
        double[] p6 =  {
                7.0, 0.0, 0.0, 3.0, 0.0, 0.0
            }
            ;
        double[] p7;
        double[] p8 =  {
                7.0,
                -14.0, 21.0,
                -25.0, +35.0, -42.0,
                -6.0, 9.0, 0.0, 0.0,
                -12.0, +15.0, -18.0, 0.0, 0.0
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("POLY_PRODUCT_TEST:");
        Console.WriteLine("  POLY_PRODUCT computes the product of two X,Y polynomials.");
        //
        //  P1 = ( 1 + 2 x + 3 y )
        //  P2 = ( 4 + 5 x )
        //  P3 = P1 * P2
        //  P4 = 4 + 13x + 12y + 10x^2 + 15xy + 0y^2 
        //
        Console.WriteLine("");
        typeMethods.poly_print(d1, ref p1, "  p1(x,y)");

        Console.WriteLine("");
        typeMethods.poly_print(d2, ref p2, "  p2(x,y)");

        d3 = d1 + d2;
        p3 = typeMethods.poly_product(d1, p1, d2, p2);
        Console.WriteLine("");
        typeMethods.poly_print(d3, ref p3, "  p3(x,y) = p1(x,y) * p2(x,y)");

        Console.WriteLine("");
        typeMethods.poly_print(d4, ref p4, "  p4(x,y) = correct answer");
        //
        //  P5 = ( 1 - 2 x + 3 y - 4x^2 + 5xy - 6y^2)
        //  P6 = ( 7 + 3x^2 )
        //  P7 = P5 * P6
        //  P8 =    7 
        //       - 14x   + 21   y 
        //       - 25x^2 + 35x  y - 42   y^2 
        //       -  6x^3 +  9x^2y +  0x  y^2 + 0  y^3
        //       - 12x^4 + 15x^3y - 18x^2y^2 + 0 xy^3 + 0y^4
        //
        Console.WriteLine("");
        typeMethods.poly_print(d5, ref p5, "  p5(x,y)");

        Console.WriteLine("");
        typeMethods.poly_print(d6, ref p6, "  p6(x,y)");

        d7 = d5 + d6;
        p7 = typeMethods.poly_product(d5, p5, d6, p6);
        Console.WriteLine("");
        typeMethods.poly_print(d7, ref p7, "  p7(x,y) = p5(x,y) * p6(x,y)");

        Console.WriteLine("");
        typeMethods.poly_print(d8, ref p8, "  p8(x,y) = Correct answer");
    }

    private static void r8mat_print_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_PRINT_TEST tests R8MAT_PRINT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 August 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int M = 6;
        int N = 4;

        double[] a = new double[M * N];
        int i;
        int j;
        int m = M;
        int n = N;

        Console.WriteLine("");
        Console.WriteLine("R8MAT_PRINT_TEST");
        Console.WriteLine("  R8MAT_PRINT prints an R8MAT.");

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < m; i++)
            {
                a[i + j * m] = (i + 1) * 10 + j + 1;
            }
        }

        typeMethods.r8mat_print(m, n, a, "  The R8MAT:");
    }

    private static void r8mat_print_some_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_PRINT_SOME_TEST tests R8MAT_PRINT_SOME.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 August 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int M = 6;
        int N = 4;

        double[] a = new double[M * N];
        int i;
        int j;
        int m = M;
        int n = N;

        Console.WriteLine("");
        Console.WriteLine("R8MAT_PRINT_SOME_TEST");
        Console.WriteLine("  R8MAT_PRINT_SOME prints some of an R8MAT.");

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < m; i++)
            {
                a[i + j * m] = (i + 1) * 10 + j + 1;
            }
        }

        typeMethods.r8mat_print_some(m, n, a, 2, 1, 4, 2, "  The R8MAT, rows 2:4, cols 1:2:");
    }

    private static void rs_to_xy_map_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RS_TO_XY_MAP_TEST tests RS_TO_XY_MAP.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 April 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a = 0;
        double b = 0;
        double c = 0;
        double d = 0;
        double e = 0;
        double f = 0;
        int j;
        double[] t =  {
                2.0, 0.0,
                3.0, 4.0,
                0.0, 3.0
            }
            ;
        double[] tr =  {
                0.0, 0.0,
                1.0, 0.0,
                0.0, 1.0
            }
            ;
        double x;
        double y;

        Console.WriteLine("");
        Console.WriteLine("RS_TO_XY_MAP_TEST:");
        Console.WriteLine("  RS_TO_XY_MAP determines the coefficients of");
        Console.WriteLine("  the linear map from a the reference in RS coordinates");
        Console.WriteLine("  to the physical triangle in XY coordinates:");
        Console.WriteLine("    X = a + b * R + c * S");
        Console.WriteLine("    Y = d + e * R + f * S");

        typeMethods.r8mat_print(2, 3, t, "  XY triangle vertices:");

        Map.rs_to_xy_map(t, ref a, ref b, ref c, ref d, ref e, ref f);

        Console.WriteLine("");
        Console.WriteLine("  Mapping coefficients are:");
        Console.WriteLine("");
        Console.WriteLine("    X = " + a + " + " + b + " * R + " + c + " * S");
        Console.WriteLine("    Y = " + d + " + " + e + " * R + " + f + " * S");

        Console.WriteLine("");
        Console.WriteLine("  Apply map to RS triangle vertices.");
        Console.WriteLine("  Recover XY vertices (2,0), (3,4) and (0,3).");
        Console.WriteLine("");
        for (j = 0; j < 3; j++)
        {
            x = a + b * tr[0 + j * 2] + c * tr[1 + j * 2];
            y = d + e * tr[0 + j * 2] + f * tr[1 + j * 2];
            Console.WriteLine("  V(" + j + ") = (" + x + "," + y + ")");
        }
    }

    private static void triangle_area_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_AREA_TEST tests TRIANGLE_AREA_MAP.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 April 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double angled;
        double angler;
        double area;
        int i;
        double r;
        const double r8_pi = 3.141592653589793;
        double[] t =  {
                0.0, 0.0,
                2.0, 0.0,
                0.0, 1.0
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("TRIANGLE_AREA_TEST:");
        Console.WriteLine("  TRIANGLE_AREA determines the (signed) area of a triangle.");

        Console.WriteLine("");
        Console.WriteLine("  Triangle vertices are:");
        Console.WriteLine("    (X1,Y1) = (0,0)");
        Console.WriteLine("    (X2,Y2) = 2*(cos(angle),sin(angle))");
        Console.WriteLine("    (X3,Y3) = (0,1)");
        Console.WriteLine("  where angle will sweep from 0 to 360 degrees.");

        r = 2.0;

        Console.WriteLine("");
        Console.WriteLine("   I      Angle         X2          Y2          Area");
        Console.WriteLine("        (degrees)");
        Console.WriteLine("");
        for (i = 0; i <= 24; i++)
        {
            angled = i * 180.0 / 12.0;
            angler = i * r8_pi / 12.0;
            t[0 + 1 * 2] = r * Math.Cos(angler);
            t[1 + 1 * 2] = r * Math.Sin(angler);
            area = Integrals.triangle_area(t);
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                   + "  " + angled.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + t[0 + 1 * 2].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + t[1 + 1 * 2].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + area.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }
    }

    private static void triangle_monomial_integral_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_MONOMIAL_INTEGRAL_TEST estimates integrals over a triangle.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 April 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int j;
        double q;
        double q2;
        double[] t1 =  {
                0.0, 0.0,
                1.0, 0.0,
                0.0, 1.0
            }
            ;
        double[] t2 =  {
                0.0, 0.0,
                1.0, 0.0,
                1.0, 2.0
            }
            ;
        double[] t3 =  {
                -3.0, 0.0,
                6.0, 0.0,
                0.0, 3.0
            }
            ;
        double[] t4 =  {
                0.0, 0.0,
                4.0, 0.0,
                0.0, 1.0
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("TRIANGLE_MONOMIAL_INTEGRAL_TEST");
        Console.WriteLine("  TRIANGLE_MONOMIAL_INTEGRAL returns the integral Q of");
        Console.WriteLine("  a monomial X^I Y^J over the interior of a triangle.");
        //
        //  Test 1:
        //
        i = 1;
        j = 0;

        Console.WriteLine("");
        Console.WriteLine("  Triangle vertices:");
        Console.WriteLine("    (" + t1[0 + 0 * 2] + "," + t1[1 + 0 * 2] + ")");
        Console.WriteLine("    (" + t1[0 + 1 * 2] + "," + t1[1 + 1 * 2] + ")");
        Console.WriteLine("    (" + t1[0 + 2 * 2] + "," + t1[1 + 2 * 2] + ")");
        Console.WriteLine("  Integrand = x^" + i + " * y^" + j + "");

        q = Integrals.triangle_monomial_integral(i, j, t1);
        q2 = 1.0 / 6.0;

        Console.WriteLine("  Computed Q = " + q + "");
        Console.WriteLine("  Exact Q =    " + q2 + "");
        //
        //  Test 2:
        //
        i = 1;
        j = 1;

        Console.WriteLine("");
        Console.WriteLine("  Triangle vertices:");
        Console.WriteLine("    (" + t2[0 + 0 * 2] + "," + t2[1 + 0 * 2] + ")");
        Console.WriteLine("    (" + t2[0 + 1 * 2] + "," + t2[1 + 1 * 2] + ")");
        Console.WriteLine("    (" + t2[0 + 2 * 2] + "," + t2[1 + 2 * 2] + ")");
        Console.WriteLine("  Integrand = x^" + i + " * y^" + j + "");

        q = Integrals.triangle_monomial_integral(i, j, t2);
        q2 = 0.5;

        Console.WriteLine("  Computed Q = " + q + "");
        Console.WriteLine("  Exact Q =    " + q2 + "");
        //
        //  Test 3:
        //
        i = 1;
        j = 0;

        Console.WriteLine("");
        Console.WriteLine("  Triangle vertices:");
        Console.WriteLine("    (" + t3[0 + 0 * 2] + "," + t3[1 + 0 * 2] + ")");
        Console.WriteLine("    (" + t3[0 + 1 * 2] + "," + t3[1 + 1 * 2] + ")");
        Console.WriteLine("    (" + t3[0 + 2 * 2] + "," + t3[1 + 2 * 2] + ")");
        Console.WriteLine("  Integrand = x^" + i + " * y^" + j + "");

        q = Integrals.triangle_monomial_integral(i, j, t3);
        q2 = 13.5;

        Console.WriteLine("  Computed Q = " + q + "");
        Console.WriteLine("  Exact Q =    " + q2 + "");
        //
        //  Test 4:
        //
        i = 1;
        j = 1;

        Console.WriteLine("");
        Console.WriteLine("  Triangle vertices:");
        Console.WriteLine("    (" + t4[0 + 0 * 2] + "," + t4[1 + 0 * 2] + ")");
        Console.WriteLine("    (" + t4[0 + 1 * 2] + "," + t4[1 + 1 * 2] + ")");
        Console.WriteLine("    (" + t4[0 + 2 * 2] + "," + t4[1 + 2 * 2] + ")");
        Console.WriteLine("  Integrand = x^" + i + " * y^" + j + "");

        q = Integrals.triangle_monomial_integral(i, j, t4);
        q2 = 2.0 / 3.0;

        Console.WriteLine("  Computed Q = " + q + "");
        Console.WriteLine("  Exact Q =    " + q2 + "");
    }

    private static void triangle_poly_integral_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_POLY_INTEGRAL_TEST estimates integrals over a triangle.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 April 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int d1 = 1;
        int d2 = 2;
        int d3 = 2;
        int d4 = 2;

        double[] p1 =  {
                0.0, 1.0, 0.0
            }
            ;
        double[] p2 =  {
                0.0, 0.0, 0.0, 0.0, 1.0, 0.0
            }
            ;
        double[] p3 =  {
                2.0, -3.0, 0.0, 0.0, 1.0, 0.0
            }
            ;
        double[] p4 =  {
                0.0, 0.0,-40.0, 6.0, 0.0, 0.0
            }
            ;
        double q;
        double q2;
        double[] t1 =  {
                0.0, 0.0,
                1.0, 0.0,
                0.0, 1.0
            }
            ;
        double[] t2 =  {
                0.0, 0.0,
                1.0, 0.0,
                1.0, 2.0
            }
            ;
        double[] t3 =  {
                0.0, 0.0,
                1.0, 0.0,
                1.0, 3.0
            }
            ;
        double[] t4 =  {
                0.0, 3.0,
                1.0, 1.0,
                5.0, 3.0
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("TRIANGLE_POLY_INTEGRAL_TEST");
        Console.WriteLine("  TRIANGLE_POLY_INTEGRAL returns the integral Q of");
        Console.WriteLine("  a polynomial over the interior of a triangle.");
        //
        //  Test 1:
        //  Integrate x over reference triangle.
        //
        Console.WriteLine("");
        Console.WriteLine("  Triangle vertices:");
        Console.WriteLine("    (" + t1[0 + 0 * 2] + "," + t1[1 + 0 * 2] + ")");
        Console.WriteLine("    (" + t1[0 + 1 * 2] + "," + t1[1 + 1 * 2] + ")");
        Console.WriteLine("    (" + t1[0 + 2 * 2] + "," + t1[1 + 2 * 2] + ")");

        typeMethods.poly_print(d1, ref p1, "  Integrand p1(x,y)");

        q = Integrals.triangle_poly_integral(d1, p1, t1);
        q2 = 1.0 / 6.0;

        Console.WriteLine("  Computed Q = " + q + "");
        Console.WriteLine("  Exact Q =    " + q2 + "");
        //
        //  Test 2:
        //  Integrate xy over a general triangle.
        //
        Console.WriteLine("");
        Console.WriteLine("  Triangle vertices:");
        Console.WriteLine("    (" + t2[0 + 0 * 2] + "," + t2[1 + 0 * 2] + ")");
        Console.WriteLine("    (" + t2[0 + 1 * 2] + "," + t2[1 + 1 * 2] + ")");
        Console.WriteLine("    (" + t2[0 + 2 * 2] + "," + t2[1 + 2 * 2] + ")");

        typeMethods.poly_print(d2, ref p2, "  Integrand p2(x,y)");

        q = Integrals.triangle_poly_integral(d2, p2, t2);
        q2 = 0.5;

        Console.WriteLine("  Computed Q = " + q + "");
        Console.WriteLine("  Exact Q =    " + q2 + "");
        //
        //  Test 3:
        //  Integrate 2-3x+xy over a general triangle.
        //
        Console.WriteLine("");
        Console.WriteLine("  Triangle vertices:");
        Console.WriteLine("    (" + t3[0 + 0 * 2] + "," + t3[1 + 0 * 2] + ")");
        Console.WriteLine("    (" + t3[0 + 1 * 2] + "," + t3[1 + 1 * 2] + ")");
        Console.WriteLine("    (" + t3[0 + 2 * 2] + "," + t3[1 + 2 * 2] + ")");

        typeMethods.poly_print(d3, ref p3, "  Integrand p3(x,y)");

        q = Integrals.triangle_poly_integral(d3, p3, t3);
        q2 = 9.0 / 8.0;

        Console.WriteLine("  Computed Q = " + q + "");
        Console.WriteLine("  Exact Q =    " + q2 + "");
        //
        //  Test 4:
        //  Integrate -40y + 6x^2 over a general triangle.
        //
        Console.WriteLine("");
        Console.WriteLine("  Triangle vertices:");
        Console.WriteLine("    (" + t4[0 + 0 * 2] + "," + t4[1 + 0 * 2] + ")");
        Console.WriteLine("    (" + t4[0 + 1 * 2] + "," + t4[1 + 1 * 2] + ")");
        Console.WriteLine("    (" + t4[0 + 2 * 2] + "," + t4[1 + 2 * 2] + ")");

        typeMethods.poly_print(d4, ref p4, "  Integrand p4(x,y)");

        q = Integrals.triangle_poly_integral(d4, p4, t4);
        q2 = -935.0 / 3.0;

        Console.WriteLine("  Computed Q = " + q + "");
        Console.WriteLine("  Exact Q =    " + q2 + "");
    }
    //****************************************************************************80

    private static void triangle01_monomial_integral_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE01_MONOMIAL_INTEGRAL_TEST estimates integrals over the unit triangle.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 April 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int d;
        int i;
        int j;
        double q;

        Console.WriteLine("");
        Console.WriteLine("TRIANGLE01_MONOMIAL_INTEGRAL_TEST");
        Console.WriteLine("  TRIANGLE01_MONOMIAL_INTEGRAL returns the integral Q of");
        Console.WriteLine("  a monomial X^I Y^J over the interior of the unit triangle.");

        Console.WriteLine("");
        Console.WriteLine("   I   J         Q(I,J)");

        for (d = 0; d <= 5; d++)
        {
            Console.WriteLine("");
            for (i = 0; i <= d; i++)
            {
                j = d - i;
                q = Integrals.triangle01_monomial_integral(i, j);
                Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                       + "  " + j.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                       + "  " + q + "");
            }
        }
    }

    private static void triangle01_poly_integral_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE01_POLY_INTEGRAL_TEST: polynomial integrals over the unit triangle.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 April 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int d_max = 6;
        int d1 = 1;
        int d2 = 2;
        int d3 = 2;

        int i = 0;
        int j = 0;
        int k;
        int km1;
        int m_max = (d_max + 1) * (d_max + 2) / 2;
        int m1 = (d1 + 1) * (d1 + 2) / 2;
        int m2 = (d2 + 1) * (d2 + 2) / 2;
        int m3 = (d3 + 1) * (d3 + 2) / 2;
        double[] p1 =  {
                1.0, 2.0, 3.0
            }
            ;
        double[] p2 =  {
                0.0, 0.0, 0.0, 0.0, 1.0, 0.0
            }
            ;
        double[] p3 =  {
                1.0, -2.0, 3.0, -4.0, 5.0, -6.0
            }
            ;
        double q;
        double q2;
        double[] qm = new double[28];

        for (k = 1; k <= m_max; k++)
        {
            typeMethods.i4_to_pascal(k, ref i, ref j);
            km1 = k - 1;
            qm[km1] = Integrals.triangle01_monomial_integral(i, j);
        }

        Console.WriteLine("");
        Console.WriteLine("TRIANGLE01_POLY_INTEGRAL_TEST");
        Console.WriteLine("  TRIANGLE01_POLY_INTEGRAL returns the integral Q of");
        Console.WriteLine("  a polynomial P(X,Y) over the interior of the unit triangle.");

        Console.WriteLine("");
        typeMethods.poly_print(d1, ref p1, "  p(x,y)");
        q = Integrals.triangle01_poly_integral(d1, p1);
        Console.WriteLine("");
        Console.WriteLine("  Q =         " + q + "");
        q2 = typeMethods.r8vec_dot_product(m1, p1, qm);
        Console.WriteLine("  Q (exact) = " + q2 + "");

        Console.WriteLine("");
        typeMethods.poly_print(d2, ref p2, "  p(x,y)");
        q = Integrals.triangle01_poly_integral(d2, p2);
        Console.WriteLine("");
        Console.WriteLine("  Q =         " + q + "");
        q2 = typeMethods.r8vec_dot_product(m2, p2, qm);
        Console.WriteLine("  Q (exact) = " + q2 + "");

        Console.WriteLine("");
        typeMethods.poly_print(d3,ref  p3, "  p(x,y)");
        q = Integrals.triangle01_poly_integral(d3, p3);
        Console.WriteLine("");
        Console.WriteLine("  Q =         " + q + "");
        q2 = typeMethods.r8vec_dot_product(m3, p3, qm);
        Console.WriteLine("  Q (exact) = " + q2 + "");
    }

    private static void triangle_xy_integral_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_XY_INTEGRAL_TEST tests TRIANGLE_XY_INTEGRAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 April 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double q;
        double x1;
        double x2;
        double x3;
        double y1;
        double y2;
        double y3;

        Console.WriteLine("");
        Console.WriteLine("TRIANGLE_XY_INTEGRAL_TEST");
        Console.WriteLine("  TRIANGLE_XY_INTEGRAL determines Q, the integral of the");
        Console.WriteLine("  monomial X*Y over a triangle (X1,Y1), (X2,Y2), (X3,Y3).");

        x1 = 0.0;
        y1 = 0.0;

        x2 = 1.0;
        y2 = 0.0;

        x3 = 1.0;
        y3 = 2.0;

        q = Integrals.triangle_xy_integral(x1, y1, x2, y2, x3, y3);

        Console.WriteLine("");
        Console.WriteLine("  (X1,Y1) = (" + x1 + "," + y1 + ")");
        Console.WriteLine("  (X2,Y2) = (" + x2 + "," + y2 + ")");
        Console.WriteLine("  (X3,Y3) = (" + x3 + "," + y3 + ")");
        Console.WriteLine("  Q = " + q + "");
        Console.WriteLine("  (Expecting answer 1/2.");

        x1 = 0.0;
        y1 = 0.0;

        x2 = 4.0;
        y2 = 0.0;

        x3 = 0.0;
        y3 = 1.0;

        q = Integrals.triangle_xy_integral(x1, y1, x2, y2, x3, y3);

        Console.WriteLine("");
        Console.WriteLine("  (X1,Y1) = (" + x1 + "," + y1 + ")");
        Console.WriteLine("  (X2,Y2) = (" + x2 + "," + y2 + ")");
        Console.WriteLine("  (X3,Y3) = (" + x3 + "," + y3 + ")");
        Console.WriteLine("  Q = " + q + "");
        Console.WriteLine("  (Expecting answer 2/3.");
    }

    private static void trinomial_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRINOMIAL_TEST tests TRINOMIAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 April 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int j;
        int k;
        int t;

        Console.WriteLine("");
        Console.WriteLine("TRINOMIAL_TEST");
        Console.WriteLine("  TRINOMIAL evaluates the trinomial coefficient:");
        Console.WriteLine("");
        Console.WriteLine("  T(I,J,K) = (I+J+K)! / I! / J! / K!");
        Console.WriteLine("");
        Console.WriteLine("     I     J     K    T(I,J,K)");
        Console.WriteLine("");

        for (k = 0; k <= 4; k++)
        {
            for (j = 0; j <= 4; j++)
            {
                for (i = 0; i <= 4; i++)
                {
                    t = Trinomial.trinomial(i, j, k);
                    Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                           + "  " + j.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                           + "  " + k.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                           + "  " + t.ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");
                }
            }
        }
    }

    private static void xy_to_rs_map_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    XY_TO_RS_MAP_TEST tests XY_TO_RS_MAP.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 April 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a = 0;
        double b = 0;
        double c = 0;
        double d = 0;
        double e = 0;
        double f = 0;
        int j;
        double r;
        double s;
        double[] t =  {
                2.0, 0.0,
                3.0, 4.0,
                0.0, 3.0
            }
            ;

        Console.WriteLine("");
        Console.WriteLine("XY_TO_RS_MAP_TEST:");
        Console.WriteLine("  XY_TO_RS_MAP determines the coefficients of the linear");
        Console.WriteLine("  map from a general triangle in XY coordinates");
        Console.WriteLine("  to the reference triangle in RS coordinates:");
        Console.WriteLine("    R = a + b * X + c * Y");
        Console.WriteLine("    S = d + e * X + f * Y");

        typeMethods.r8mat_print(2, 3, t, "  XY triangle vertices:");

        Map.xy_to_rs_map(t, ref a, ref b, ref c, ref d, ref e, ref f);

        Console.WriteLine("");
        Console.WriteLine("  Mapping coefficients are:");
        Console.WriteLine("");
        Console.WriteLine("    R = " + a + " + " + b + " * X + " + c + " * Y");
        Console.WriteLine("    S = " + d + " + " + e + " * X + " + f + " * Y");

        Console.WriteLine("");
        Console.WriteLine("  Apply map to XY triangle vertices.");
        Console.WriteLine("  Recover RS vertices (0,0), (1,0) and (0,1).");
        Console.WriteLine("");
        for (j = 0; j < 3; j++)
        {
            r = a + b * t[0 + j * 2] + c * t[1 + j * 2];
            s = d + e * t[0 + j * 2] + f * t[1 + j * 2];
            Console.WriteLine("  V[" + j + "] = (" + r + "," + s + ")");
        }
    }
}