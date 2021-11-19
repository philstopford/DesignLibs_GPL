using System;
using Burkardt.PolynomialNS;

namespace PolynomialTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for POLYNOMIAL_TEST.
        //
        //  Discussion:
        //
        //    POLYNOMIAL_TEST tests the POLYNOMIAL library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 October 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("POLYNOMIAL_TEST");
        Console.WriteLine("  Test the POLYNOMIAL library.");

        polynomial_add_test();
        polynomial_axpy_test();
        polynomial_compress_test();
        polynomial_dif_test();
        polynomial_mul_test();
        polynomial_print_test();
        polynomial_scale_test();
        polynomial_sort_test();
        polynomial_value_test();

        Console.WriteLine("");
        Console.WriteLine("POLYNOMIAL_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void polynomial_add_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYNOMIAL_ADD_TEST tests POLYNOMIAL_ADD.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 December 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] c = new double[11];
        double[] c1 = { 7.0, -5.0, 9.0, 11.0, 0.0, -13.0 };
        double[] c2 = { 2.0, 3.0, -8.0, 4.0, 9.0 };
        int m = 3;
        int[] e = new int[11];
        int[] e1 = { 1, 2, 4, 5, 12, 33 };
        int[] e2 = { 1, 3, 4, 30, 33 };
        int o = 0;
        int o1 = 6;
        int o2 = 5;
        string title = "  P1(X) + P2(X) =";
        string title1 = "  P1(X) =";
        string title2 = "  P2(X) =";

        Console.WriteLine("");
        Console.WriteLine("POLYNOMIAL_ADD_TEST");
        Console.WriteLine("  POLYNOMIAL_ADD adds two polynomials.");

        Console.WriteLine("");
        Polynomial.polynomial_print(m, o1, c1, e1, title1);

        Console.WriteLine("");
        Polynomial.polynomial_print(m, o2, c2, e2, title2);

        Polynomial.polynomial_add(o1, c1, e1, o2, c2, e2, ref o, ref c, ref e);
        Console.WriteLine("");
        Polynomial.polynomial_print(m, o, c, e, title);
    }

    private static void polynomial_axpy_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYNOMIAL_AXPY_TEST tests POLYNOMIAL_ADD.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] c = new double[11];
        double[] c1 = { 7.0, -5.0, 9.0, 11.0, 0.0, -13.0 };
        double[] c2 = { 2.0, 3.0, -8.0, 4.0, 9.0 };
        int m = 3;
        int[] e = new int[11];
        int[] e1 = { 1, 2, 4, 5, 12, 33 };
        int[] e2 = { 1, 3, 4, 30, 33 };
        int o = 0;
        int o1 = 6;
        int o2 = 5;
        double s = 10.0;
        string title = "  S * P1(X) + P2(X) =";
        string title1 = "  P1(X) =";
        string title2 = "  P2(X) =";

        Console.WriteLine("");
        Console.WriteLine("POLYNOMIAL_AXPY_TEST");
        Console.WriteLine("  POLYNOMIAL_AXPY adds a multiple of one polynomial to another.");

        Console.WriteLine("");
        Polynomial.polynomial_print(m, o1, c1, e1, title1);

        Console.WriteLine("");
        Polynomial.polynomial_print(m, o2, c2, e2, title2);

        Polynomial.polynomial_axpy(s, o1, c1, e1, o2, c2, e2, ref o, ref c, ref e);
        Console.WriteLine("");
        Polynomial.polynomial_print(m, o, c, e, title);
    }

    private static void polynomial_compress_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYNOMIAL_COMPRESS_TEST tests POLYNOMIAL_COMPRESS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 October 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] c = { 7.0, -5.0, 5.0, 9.0, 11.0, 3.0, 6.0, 0.0, -13.0, 1.0E-20 };
        double[] c2 = new double[10];
        int m = 3;
        int[] e = { 1, 2, 2, 4, 5, 5, 5, 12, 33, 35 };
        int[] e2 = new int[10];
        int o = 10;
        int o2 = 0;
        string title;

        Console.WriteLine("");
        Console.WriteLine("POLYNOMIAL_COMPRESS_TEST");
        Console.WriteLine("  POLYNOMIAL_COMPRESS compresses a polynomial.");

        Console.WriteLine("");
        title = "  Uncompressed P(X) = ";
        Polynomial.polynomial_print(m, o, c, e, title);

        Polynomial.polynomial_compress(o, c, e, ref o2, ref c2, ref e2);

        Console.WriteLine("");
        title = "  Compressed P(X) = ";
        Polynomial.polynomial_print(m, o2, c2, e2, title);
    }

    private static void polynomial_dif_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYNOMIAL_DIF_TEST tests POLYNOMIAL_DIF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 December 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] c = new double[4];
        double[] c1 = { 2.0, 3.0, 4.0, 5.0 };
        int m = 2;
        int[] dif = { 2, 1 };
        int[] e = new int[4];
        int[] e1 = { 1, 10, 12, 32 };
        int o = 0;
        int o1 = 4;
        string title = "  d3 P(X) dx1 dx1 dx2 =";
        string title1 = "  P(X) =";

        Console.WriteLine("");
        Console.WriteLine("POLYNOMIAL_DIF_TEST");
        Console.WriteLine("  POLYNOMIAL_DIF computes derivatives of a polynomial.");

        Console.WriteLine("");
        Polynomial.polynomial_print(m, o1, c1, e1, title1);

        Polynomial.polynomial_dif(m, o1, c1, e1, dif, ref o, ref c, ref e);
        Console.WriteLine("");
        Polynomial.polynomial_print(m, o, c, e, title);
    }

    private static void polynomial_mul_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYNOMIAL_MUL_TEST tests POLYNOMIAL_MUL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 December 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] c = new double[8];
        double[] c1 = { 2.0, 3.0, 4.0, 5.0 };
        double[] c2 = { 6.0, 7.0 };
        int m = 3;
        int[] e = new int[8];
        int[] e1 = { 1, 3, 4, 6 };
        int[] e2 = { 2, 5 };
        int o = 0;
        int o1 = 4;
        int o2 = 2;
        string title = "  P1(X) * P2(X) =";
        string title1 = "  P1(X) =";
        string title2 = "  P2(X) =";

        Console.WriteLine("");
        Console.WriteLine("POLYNOMIAL_MUL_TEST");
        Console.WriteLine("  POLYNOMIAL_MUL multiplies two polynomials.");

        Console.WriteLine("");
        Polynomial.polynomial_print(m, o1, c1, e1, title1);

        Console.WriteLine("");
        Polynomial.polynomial_print(m, o2, c2, e2, title2);

        Polynomial.polynomial_mul(m, o1, c1, e1, o2, c2, e2, ref o, ref c, ref e);
        Console.WriteLine("");
        Polynomial.polynomial_print(m, o, c, e, title);
    }

    private static void polynomial_print_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYNOMIAL_PRINT_TEST tests POLYNOMIAL_PRINT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 October 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] c = { 7.0, -5.0, 9.0, 11.0, 0.0, -13.0 };
        int m = 3;
        int[] e = { 1, 2, 4, 5, 12, 33 };
        int o = 6;
        string title = "  P1(X) =";

        Console.WriteLine("");
        Console.WriteLine("POLYNOMIAL_PRINT_TEST");
        Console.WriteLine("  POLYNOMIAL_PRINT prints a polynomial.");

        Console.WriteLine("");
        Polynomial.polynomial_print(m, o, c, e, title);
    }

    private static void polynomial_scale_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYNOMIAL_SCALE_TEST tests POLYNOMIAL_SCALE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 October 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] c = { 7.0, -5.0, 9.0, 11.0, 0.0, -13.0 };
        int m = 3;
        int[] e = { 1, 2, 4, 5, 12, 33 };
        int o = 6;
        double s;
        string title;

        Console.WriteLine("");
        Console.WriteLine("POLYNOMIAL_PRINT_TEST");
        Console.WriteLine("  POLYNOMIAL_PRINT prints a polynomial.");

        Console.WriteLine("");
        title = "  Original P(X):";
        Polynomial.polynomial_print(m, o, c, e, title);

        s = -0.5;
        Console.WriteLine("");
        Console.WriteLine("  Apply scale factor S = " + s + "");
        Polynomial.polynomial_scale(s, m, o, ref c, e);

        Console.WriteLine("");
        title = "  S * P(X):";
        Polynomial.polynomial_print(m, o, c, e, title);
    }

    private static void polynomial_sort_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYNOMIAL_SORT_TEST tests POLYNOMIAL_SORT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 October 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] c = { 0.0, 9.0, -5.0, -13.0, 7.0, 11.0 };
        int m = 3;
        int[] e = { 12, 4, 2, 33, 1, 5 };
        int o = 6;
        string title;

        Console.WriteLine("");
        Console.WriteLine("POLYNOMIAL_SORT_TEST");
        Console.WriteLine("  POLYNOMIAL_SORT sorts a polynomial by exponent index.");

        Console.WriteLine("");
        title = "  Unsorted polynomial:";
        Polynomial.polynomial_print(m, o, c, e, title);

        Polynomial.polynomial_sort(o, ref c, ref e);

        Console.WriteLine("");
        title = "  Sorted polynomial:";
        Polynomial.polynomial_print(m, o, c, e, title);
    }

    private static void polynomial_value_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYNOMIAL_VALUE_TEST tests POLYNOMIAL_VALUE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 October 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] c = { 7.0, -5.0, 9.0, 11.0, 0.0, -13.0 };
        int m = 3;
        int[] e = { 1, 2, 4, 5, 12, 33 };
        int j;
        int nx = 2;
        int o = 6;
        double[] p;
        string title = "  P(X) =";
        double[] x =
        {
            1.0, 2.0, 3.0,
            -2.0, 4.0, 1.0
        };

        Console.WriteLine("");
        Console.WriteLine("POLYNOMIAL_VALUE_TEST");
        Console.WriteLine("  POLYNOMIAL_VALUE evaluates a polynomial.");

        Console.WriteLine("");
        Polynomial.polynomial_print(m, o, c, e, title);

        p = Polynomial.polynomial_value(m, o, c, e, nx, x);

        Console.WriteLine("");
        for (j = 0; j < nx; j++)
        {
            Console.WriteLine("  P(" + x[0 + j * m]
                                     + "," + x[1 + j * m]
                                     + "," + x[2 + j * m]
                                     + ") = " + p[j] + "");
        }
    }
}