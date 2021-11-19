using System;
using System.Globalization;
using Burkardt.ChebyshevNS;

namespace ChebyshevTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for CHEBYSHEV_TEST.
        //
        //  Discussion:
        //
        //    CHEBYSHEV_TEST tests the CHEBYSHEV library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 September 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("CHEBYSHEV_TEST");
        Console.WriteLine("  Test the CHEBYSHEV library.");

        test01();

        Console.WriteLine("");
        Console.WriteLine("CHEBYSHEV_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests CHEBYSHEV_COEFFICIENTS and CHEBYSHEV_INTERPOLANT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 September 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("CHEBYSHEV_TEST01");
        Console.WriteLine("  CHEBYSHEV_COEFFICIENTS computes the coefficients of the");
        Console.WriteLine("  Chebyshev interpolant.");
        Console.WriteLine("  CHEBYSHEV_INTERPOLANT evaluates the interpolant.");

        int n = 5;
        double a = -1.0;
        double b = +1.0;

        double[] c = Chebyshev.chebyshev_coefficients(a, b, n, f1);

        double[] x = Chebyshev.chebyshev_zeros(n);
        for (i = 0; i < n; i++)
        {
            x[i] = 0.5 * (a + b) + x[i] * 0.5 * (b - a);
        }

        int m = n;
        double[] fc = Chebyshev.chebyshev_interpolant(a, b, n, c, m, x);

        Console.WriteLine("");
        Console.WriteLine("  F(X) is a trig function:");
        Console.WriteLine("");
        Console.WriteLine("          X               C(I)            F(X)           C(F)(X)");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + c[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + f1(x[i]).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + fc[i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        //
        //  Try a variant interval.
        //
        n = 5;
        a = 0.0;
        b = +3.0;

        c = Chebyshev.chebyshev_coefficients(a, b, n, f1);

        x = Chebyshev.chebyshev_zeros(n);
        for (i = 0; i < n; i++)
        {
            x[i] = 0.5 * (a + b) + x[i] * 0.5 * (b - a);
        }

        m = n;
        fc = Chebyshev.chebyshev_interpolant(a, b, n, c, m, x);

        Console.WriteLine("");
        Console.WriteLine("  Consider the same F(X), but now over [0,3]:");
        Console.WriteLine("");
        Console.WriteLine("          X               C(I)            F(X)           C(F)(X)");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + c[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + f1(x[i]).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + fc[i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        //
        //  Try a higher order.
        //
        n = 10;
        a = -1.0;
        b = +1.0;

        c = Chebyshev.chebyshev_coefficients(a, b, n, f1);

        x = Chebyshev.chebyshev_zeros(n);
        for (i = 0; i < n; i++)
        {
            x[i] = 0.5 * (a + b) + x[i] * 0.5 * (b - a);
        }

        m = n;
        fc = Chebyshev.chebyshev_interpolant(a, b, n, c, m, x);

        Console.WriteLine("");
        Console.WriteLine("  Consider the same F(X), but now with higher order:");
        Console.WriteLine("");
        Console.WriteLine("          X               C(I)            F(X)           C(F)(X)");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + c[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + f1(x[i]).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + fc[i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        //
        //  Try a polynomial.
        //
        n = 10;
        a = -1.0;
        b = +1.0;

        c = Chebyshev.chebyshev_coefficients(a, b, n, f3);

        x = Chebyshev.chebyshev_zeros(n);
        for (i = 0; i < n; i++)
        {
            x[i] = 0.5 * (a + b) + x[i] * 0.5 * (b - a);
        }

        m = n;
        fc = Chebyshev.chebyshev_interpolant(a, b, n, c, m, x);

        Console.WriteLine("");
        Console.WriteLine("  F(X) is a degree 4 polynomial:");
        Console.WriteLine("");
        Console.WriteLine("          X               C(I)            F(X)           C(F)(X)");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + c[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + f3(x[i]).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + fc[i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        //
        //  Try a function with decaying behavior.
        //
        n = 10;
        a = -1.0;
        b = +1.0;

        c = Chebyshev.chebyshev_coefficients(a, b, n, f2);

        x = Chebyshev.chebyshev_zeros(n);
        for (i = 0; i < n; i++)
        {
            x[i] = 0.5 * (a + b) + x[i] * 0.5 * (b - a);
        }

        m = n;
        fc = Chebyshev.chebyshev_interpolant(a, b, n, c, m, x);

        Console.WriteLine("");
        Console.WriteLine("  The polynomial approximation to F(X) decays:");
        Console.WriteLine("");
        Console.WriteLine("          X               C(I)            F(X)           C(F)(X)");
        Console.WriteLine("");
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + c[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + f2(x[i]).ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + fc[i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

    }

    private static double f1(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F1 evaluates a function that can be used for Chebyshev interpolation.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 September 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, a point where the function is to be evaluated.
        //
        //    Output, double F1, the function value.
        //
    {

        double value = Math.Cos(2.0 * Math.PI * x) * Math.Sin(3.0 * Math.PI * x);

        return value;
    }

    private static double f2(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F2 evaluates a function that can be used for Chebyshev interpolation.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 September 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input,  double X, a point where the function is to be evaluated.
        //
        //    Output, double F2, the function value.
        //
    {
        double value = Math.Exp(x);

        return value;
    }

    private static double f3(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F3 evaluates a function that can be used for Chebyshev interpolation.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 September 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, a point where the function is to be evaluated.
        //
        //    Output, double F3, the function values.
        //
    {
        double value = (x - 3.0) * (x - 1.0) * (x - 1.0) * (x + 2.0);

        return value;
    }
}