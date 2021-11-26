using System;
using System.Globalization;
using Burkardt.IntegralNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace HermiteCubicTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for HERMITE_CUBIC_TEST.
        //
        //  Discussion:
        //
        //    HERMITE_CUBIC_TEST tests the HERMITE_CUBIC library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 March 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("HERMITE_CUBIC_TEST");
        Console.WriteLine("  Test the HERMITE_CUBIC library.");

        test01();
        test02();
        test03();
        test04();
        test05();
        test06();
        test07();
        test08();
        test09();
        test10();
        test11();
        test12();
        test13();
        test14();
        test15();

        Console.WriteLine("");
        Console.WriteLine("HERMITE_CUBIC_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests HERMITE_CUBIC_VALUE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 February 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] d = new double[1];
        double[] f = new double[1];
        const int n = 1;
        double[] s = new double[1];
        double[] t = new double[1];
        double[] x = new double[1];
        int x_interval;

        Console.WriteLine("");
        Console.WriteLine("TEST01:");
        Console.WriteLine("  HERMITE_CUBIC_VALUE evaluates a Hermite cubic polynomial.");
        Console.WriteLine("  Try out four sets of data:");
        Console.WriteLine("  (F1,D1,F2,D2) = (1,0,0,0), (0,1,0,0), (0,0,1,0), (0,0,0,1)");
        Console.WriteLine("  on [0,1] and [1.0,-2.0] (interval reversed)");

        for (x_interval = 1; x_interval <= 2; x_interval++)
        {
            double x1;
            double x2;
            switch (x_interval)
            {
                case 1:
                    x1 = 0.0;
                    x2 = 1.0;
                    break;
                default:
                    x1 = 1.0;
                    x2 = -2.0;
                    break;
            }

            int i;
            for (i = 1; i <= 4; i++)
            {
                double f1 = 0.0;
                double d1 = 0.0;
                double f2 = 0.0;
                double d2 = 0.0;

                switch (i)
                {
                    case 1:
                        f1 = 1.0;
                        break;
                    case 2:
                        d1 = 1.0;
                        break;
                    case 3:
                        f2 = 1.0;
                        break;
                    case 4:
                        d2 = 1.0;
                        break;
                }

                Console.WriteLine("");
                Console.WriteLine("    J      X           F           D");
                Console.WriteLine("");

                int j;
                for (j = -3; j <= 12; j++)
                {
                    x[0] = ((10 - j) * x1
                            + j * x2)
                           / 10;

                    HermiteCubic.hermite_cubic_value(x1, f1, d1, x2, f2, d2, n, x, ref f, ref d, ref s, ref t);

                    switch (j)
                    {
                        case 0:
                            Console.WriteLine("*Data"
                                              + "  " + x1.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                              + "  " + f1.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                              + "  " + d1.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
                            break;
                    }

                    Console.WriteLine("  " + j.ToString(CultureInfo.InvariantCulture).PadLeft(3)
                                           + "  " + x[0].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                           + "  " + f[0].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                           + "  " + d[0].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
                    switch (j)
                    {
                        case 10:
                            Console.WriteLine("*Data"
                                              + "  " + x2.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                              + "  " + f2.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                              + "  " + d2.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
                            break;
                    }
                }
            }
        }
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests HERMITE_CUBIC_VALUE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 February 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] d = new double[1];
        double dc = 0;
        double d1 = 0;
        double d2 = 0;
        double[] f = new double[1];
        double fc = 0;
        double f1 = 0;
        double f2 = 0;
        int n = 1;
        double[] s = new double[1];
        double s1 = 0;
        double s2 = 0;
        double sc = 0;
        double[] t = new double[1];
        double t1 = 0;
        double t2 = 0;
        double tc = 0;
        double[] x = new double[1];
        int x_interval;
        double x1 = 0;
        double x2 = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST02:");
        Console.WriteLine("  HERMITE_CUBIC_VALUE evaluates a Hermite cubic polynomial.");
        Console.WriteLine("  Try out data from a cubic function:");
        Console.WriteLine("  on [0,10] and [-1.0,1.0] and [0.5,0.75]");

        for (x_interval = 1; x_interval <= 3; x_interval++)
        {
            switch (x_interval)
            {
                case 1:
                    x1 = 0.0;
                    x2 = 10.0;
                    break;
                case 2:
                    x1 = -1.0;
                    x2 = +1.0;
                    break;
                case 3:
                    x1 = 0.5;
                    x2 = 0.75;
                    break;
            }

            cubic_value(x1, ref f1, ref d1, ref s1, ref t1);
            cubic_value(x2, ref f2, ref d2, ref s2, ref t2);

            Console.WriteLine("");
            Console.WriteLine("    J      X           F           D           S           T");
            Console.WriteLine("");

            int j;
            for (j = -3; j <= 12; j++)
            {
                x[0] = ((10 - j) * x1
                        + j * x2)
                       / 10;

                HermiteCubic.hermite_cubic_value(x1, f1, d1, x2, f2, d2, n, x, ref f, ref d, ref s, ref t);
                cubic_value(x[0], ref fc, ref dc, ref sc, ref tc);

                switch (j)
                {
                    case 0:
                        Console.WriteLine("*Data"
                                          + "  " + x1.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                          + "  " + f1.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                          + "  " + d1.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
                        break;
                }

                Console.WriteLine("Exact"
                                  + "  " + x[0].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                  + "  " + fc.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                  + "  " + dc.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                  + "  " + sc.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                  + "  " + tc.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
                Console.WriteLine("  " + j.ToString(CultureInfo.InvariantCulture).PadLeft(3)
                                       + "  " + x[0].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + f[0].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + d[0].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + s[0].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + t[0].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
                switch (j)
                {
                    case 10:
                        Console.WriteLine("*Data"
                                          + "  " + x2.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                          + "  " + f2.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                          + "  " + d2.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
                        break;
                }
            }
        }
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests HERMITE_CUBIC_INTEGRATE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 February 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double d1 = 0;
        double d2 = 0;
        double f1 = 0;
        double f2 = 0;
        double s1 = 0;
        double s2 = 0;
        double t1 = 0;
        double t2 = 0;
        int x_interval;
        double x1 = 0;
        double x2 = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST03:");
        Console.WriteLine("  HERMITE_CUBIC_INTEGRATE integrates a Hermite cubic");
        Console.WriteLine("  polynomial from A to B.");

        for (x_interval = 1; x_interval <= 3; x_interval++)
        {
            switch (x_interval)
            {
                case 1:
                    x1 = 0.0;
                    x2 = 10.0;
                    break;
                case 2:
                    x1 = -1.0;
                    x2 = +1.0;
                    break;
                case 3:
                    x1 = 0.5;
                    x2 = 0.75;
                    break;
            }

            cubic_value(x1, ref f1, ref d1, ref s1, ref t1);
            cubic_value(x2, ref f2, ref d2, ref s2, ref t2);

            Console.WriteLine("");
            Console.WriteLine("                                     Exact           Computed");
            Console.WriteLine("    J          A           B         Integral        Integral");
            Console.WriteLine("");

            double a = x1 - 1.0;

            int j;
            for (j = -3; j <= 12; j++)
            {
                double b = ((10 - j) * x1
                            + j * x2)
                           / 10;

                double q_exact = cubic_integrate(a, b);

                double q_computed = HermiteCubic.hermite_cubic_integrate(x1, f1, d1, x2, f2, d2, a, b);

                Console.WriteLine("  " + j.ToString(CultureInfo.InvariantCulture).PadLeft(3)
                                       + "  " + a.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + b.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + q_exact.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + q_computed.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }
        }
    }

    private static void test04()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 tests HERMITE_CUBIC_SPLINE_VALUE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 February 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        const int n = 51;
        const int nn = 11;
        double u;
        double v;

        Console.WriteLine("");
        Console.WriteLine("TEST04:");
        Console.WriteLine("  HERMITE_CUBIC_SPLINE_VALUE evaluates a Hermite cubic spline.");

        const double x1 = 0.0;
        const double x2 = 10.0;

        double[] xn = typeMethods.r8vec_even_new(nn, x1, x2);
        double[] fn = new double[nn];
        double[] dn = new double[nn];

        for (i = 0; i < nn; i++)
        {
            fn[i] = Math.Sin(xn[i]);
            dn[i] = Math.Cos(xn[i]);
        }

        double[] x = typeMethods.r8vec_even_new(n, x1, x2);
        double[] f = new double[n];
        double[] d = new double[n];
        double[] s = new double[n];
        double[] t = new double[n];

        HermiteCubic.hermite_cubic_spline_value(nn, xn, fn, dn, n, x, ref f, ref d, ref s, ref t);

        Console.WriteLine("");
        Console.WriteLine("     I      X       F computed     F exact      Error");
        Console.WriteLine("");

        for (i = 0; i < n; i++)
        {
            u = Math.Sin(x[i]);
            v = Math.Abs(f[i] - u);
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + f[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + u.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + v.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("     I      X       D computed     D exact      Error");
        Console.WriteLine("");

        for (i = 0; i < n; i++)
        {
            u = Math.Cos(x[i]);
            v = Math.Abs(d[i] - u);
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + d[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + u.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + v.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    private static void test05()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 tests HERMITE_CUBIC_TO_POWER_CUBIC
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 February 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double c0 = 0;
        double c1 = 0;
        double c2 = 0;
        double c3 = 0;
        double[] d = new double[1];
        double d1 = 0;
        double d1r = 0;
        double d2 = 0;
        double d2r = 0;
        double[] f = new double[1];
        double f1 = 0;
        double f1r = 0;
        double f2 = 0;
        double f2r = 0;
        int j;
        const int n = 1;
        double[] s = new double[1];
        double s1 = 0;
        double s2 = 0;
        double[] t = new double[1];
        double t1 = 0;
        double t2 = 0;
        double[] x = new double[1];

        Console.WriteLine("");
        Console.WriteLine("TEST05:");
        Console.WriteLine("  HERMITE_CUBIC_TO_POWER_CUBIC converts the Hermite data");
        Console.WriteLine("  to the coefficients of the power form of the polynomial");
        Console.WriteLine("  POWER_CUBIC_TO_HERMITE_CUBIC converts the power form");
        Console.WriteLine("  to Hermite form");

        const double x1 = -1.0;
        const double x2 = +1.0;

        cubic_value(x1, ref f1, ref d1, ref s1, ref t1);
        cubic_value(x2, ref f2, ref d2, ref s2, ref t2);

        Console.WriteLine("");
        Console.WriteLine("  Hermite data:");
        Console.WriteLine("");
        Console.WriteLine("  X1, F1, D1:" + x1.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                          + "  " + f1.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                          + "  " + d1.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        Console.WriteLine("  X2, F2, D2:" + x2.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                          + "  " + f2.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                          + "  " + d2.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");

        HermiteCubic.hermite_cubic_to_power_cubic(x1, f1, d1, x2, f2, d2, ref c0, ref c1, ref c2, ref c3);

        Console.WriteLine("");
        Console.WriteLine("  Power form:");
        Console.WriteLine("  p(x) = " + c0 + " + " + c1 + " * x + "
                          + c2 + " * x^2 + " + c3 + " * x^3");
        Console.WriteLine("");
        Console.WriteLine("      X       F (Hermite)  F (power)");
        Console.WriteLine("");

        for (j = -3; j <= 12; j++)
        {
            x[0] = ((10 - j) * x1
                    + j * x2)
                   / 10;

            HermiteCubic.hermite_cubic_value(x1, f1, d1, x2, f2, d2, n, x, ref f, ref d, ref s, ref t);

            double fp = c0 + x[0] * (c1 + x[0] * (c2 + x[0] * c3));

            Console.WriteLine("  " + x[0].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + f[0].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + fp.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

        HermiteCubic.power_cubic_to_hermite_cubic(c0, c1, c2, c3, x1, x2, ref f1r, ref d1r,
            ref f2r, ref d2r);

        Console.WriteLine("");
        Console.WriteLine("  Use POWER_CUBIC_TO_HERMITE_CUBIC to recover the");
        Console.WriteLine("  original Hermite data:");
        Console.WriteLine("");
        Console.WriteLine("         Original   Recovered");
        Console.WriteLine("");
        Console.WriteLine("  F1:  " + "  " + f1.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                          + "  " + f1r.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        Console.WriteLine("  D1:  " + "  " + d1.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                          + "  " + d1r.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        Console.WriteLine("  F2:  " + "  " + f2.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                          + "  " + f2r.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        Console.WriteLine("  D2:  " + "  " + d2.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                          + "  " + d2r.ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
    }

    private static void test06()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST06 tests HERMITE_CUBIC_INTEGRATE using vectors.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 February 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double d1 = 0;
        double d2 = 0;
        double f1 = 0;
        double f2 = 0;
        int i;
        double s1 = 0;
        double s2 = 0;
        double t1 = 0;
        double t2 = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST06:");
        Console.WriteLine("  HERMITE_CUBIC_INTEGRATE integrates a Hermite cubic");
        Console.WriteLine("  polynomial from A to B.");
        Console.WriteLine("  Use A, B vectors for the calculation.");

        const double x1 = 0.0;
        const double x2 = 10.0;

        cubic_value(x1, ref f1, ref d1, ref s1, ref t1);
        cubic_value(x2, ref f2, ref d2, ref s2, ref t2);

        Console.WriteLine("");
        Console.WriteLine("                                 Exact       Computed");
        Console.WriteLine("    J      A           B         Integral    Integral");
        Console.WriteLine("");

        for (i = -3; i <= 12; i++)
        {
            double a = x1 - 1.0;
            double b = ((10 - i) * x1
                        + i * x2)
                       / 10;

            double q_exact = cubic_integrate(a, b);

            double q_computed = HermiteCubic.hermite_cubic_integrate(x1, f1, d1, x2, f2, d2, a, b);

            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(3)
                                   + "  " + a.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + b.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + q_exact.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + q_computed.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    private static void test07()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST07 tests HERMITE_CUBIC_INTEGRAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 February 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double d1 = 0;
        double d2 = 0;
        double f1 = 0;
        double f2 = 0;
        double s1 = 0;
        double s2 = 0;
        double t1 = 0;
        double t2 = 0;
        int x_interval;
        double x1 = 0;
        double x2 = 0;

        Console.WriteLine("");
        Console.WriteLine("TEST07:");
        Console.WriteLine("  HERMITE_CUBIC_INTEGRAL integrates a Hermite cubic");
        Console.WriteLine("  polynomial over the definition interval [X1,X2].");
        Console.WriteLine("");
        Console.WriteLine("                            Exact       Computed");
        Console.WriteLine("     X1          X2         Integral    Integral");
        Console.WriteLine("");

        for (x_interval = 1; x_interval <= 3; x_interval++)
        {
            switch (x_interval)
            {
                case 1:
                    x1 = 0.0;
                    x2 = 10.0;
                    break;
                case 2:
                    x1 = -1.0;
                    x2 = +1.0;
                    break;
                case 3:
                    x1 = 0.5;
                    x2 = 0.75;
                    break;
            }

            cubic_value(x1, ref f1, ref d1, ref s1, ref t1);
            cubic_value(x2, ref f2, ref d2, ref s2, ref t2);

            double q_exact = cubic_integrate(x1, x2);

            double q_computed = HermiteCubic.hermite_cubic_integral(x1, f1, d1, x2, f2, d2);

            Console.WriteLine("  " + x1.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + x2.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + q_exact.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + q_computed.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    private static void test08()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST08 tests HERMITE_CUBIC_SPLINE_INTEGRAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 February 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int nn = 11;
        double q_exact = 0;
        int test;
        double[] xn = new double[1];

        Console.WriteLine("");
        Console.WriteLine("TEST08:");
        Console.WriteLine("  HERMITE_CUBIC_SPLINE_INTEGRAL integrates a Hermite");
        Console.WriteLine("  cubic spline over the definition interval [X1,XNN].");
        Console.WriteLine("");
        Console.WriteLine("                            Exact       Computed");
        Console.WriteLine("     X1          XNN        Integral    Integral");
        Console.WriteLine("");

        double[] fn = new double[nn];
        double[] dn = new double[nn];

        for (test = 1; test <= 3; test++)
        {
            double b;
            double a;
            int i;
            switch (test)
            {
                case 1:
                {
                    a = 0.0;
                    b = 1.0;

                    xn = typeMethods.r8vec_even_new(nn, a, b);
                    for (i = 0; i < nn; i++)
                    {
                        fn[i] = xn[i] * (4.0 * xn[i] - 1.0) * (xn[i] - 1.0);
                        dn[i] = 1.0 + xn[i] * (-10.0 + xn[i] * 12.0);
                    }

                    q_exact =
                        xn[nn - 1] * xn[nn - 1] * (0.5 + xn[nn - 1] * (-(5.0 / 3.0) + xn[nn - 1]))
                        - xn[0] * xn[0] * (0.5 + xn[0] * (-(5.0 / 3.0) + xn[0]));
                    break;
                }
                //
                //  Use variable spacing.
                //
                case 2:
                {
                    a = 0.0;
                    b = 1.0;

                    xn = typeMethods.r8vec_even_new(nn, a, b);
                    for (i = 0; i < nn; i++)
                    {
                        xn[i] = Math.Sqrt(xn[i]);
                        fn[i] = xn[i] * (4.0 * xn[i] - 1.0) * (xn[i] - 1.0);
                        dn[i] = 1.0 + xn[i] * (-10.0 + xn[i] * 12.0);
                    }

                    q_exact =
                        xn[nn - 1] * xn[nn - 1] * (0.5 + xn[nn - 1] * (-(5.0 / 3.0) + xn[nn - 1]))
                        - xn[0] * xn[0] * (0.5 + xn[0] * (-(5.0 / 3.0) + xn[0]));
                    break;
                }
                //
                //  Try a non-cubic.
                //
                case 3:
                {
                    a = 0.0;
                    b = Math.PI;

                    xn = typeMethods.r8vec_even_new(nn, a, b);
                    for (i = 0; i < nn; i++)
                    {
                        fn[i] = Math.Sin(xn[i]);
                        dn[i] = Math.Cos(xn[i]);
                    }

                    q_exact = -Math.Cos(xn[nn - 1]) + Math.Cos(xn[0]);
                    break;
                }
            }

            double q_computed = HermiteCubic.hermite_cubic_spline_integral(nn, xn, fn, dn);

            Console.WriteLine("  " + xn[0].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + xn[nn - 1].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + q_exact.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + q_computed.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    private static void test09()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST09 tests HERMITE_CUBIC_SPLINE_INTEGRATE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 February 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        const int n = 25;
        const int nn = 11;

        Console.WriteLine("");
        Console.WriteLine("TEST09:");
        Console.WriteLine("  HERMITE_CUBIC_SPLINE_INTEGRATE integrates a Hermite");
        Console.WriteLine("  cubic spline from A to B.");
        //
        //  Define the cubic spline.
        //
        const double x1 = 0.0;
        const double x2 = 10.0;

        double[] xn = typeMethods.r8vec_even_new(nn, x1, x2);
        double[] fn = new double[nn];
        double[] dn = new double[nn];
        double[] sn = new double[nn];
        double[] tn = new double[nn];

        for (i = 0; i < nn; i++)
        {
            cubic_value(xn[i], ref fn[i], ref dn[i], ref sn[i], ref tn[i]);
        }

        double[] a = new double[n];
        for (i = 0; i < n; i++)
        {
            a[i] = 2.5;
        }

        double[] b = typeMethods.r8vec_even_new(n, x1 - 1.0, x2 + 1.0);

        double[] q = HermiteCubic.hermite_cubic_spline_integrate(nn, xn, fn, dn, n, a, b);

        Console.WriteLine("");
        Console.WriteLine("                                 Exact       Computed");
        Console.WriteLine("    I      A           B         Integral    Integral");
        Console.WriteLine("");

        for (i = 0; i < n; i++)
        {
            double q_exact = cubic_integrate(a[i], b[i]);

            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(3)
                                   + "  " + a[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + b[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + q_exact.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + q[i].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }
    }

    private static void test10()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST10 tests HERMITE_CUBIC_SPLINE_INTEGRAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 February 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        string comment = "";
        const int nn = 11;
        double q_exact = 0;
        int test;
        double[] xn = new double[1];

        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("TEST10:");
        Console.WriteLine("  HERMITE_CUBIC_SPLINE_INTEGRAL integrates a Hermite");
        Console.WriteLine("  cubic spline over the definition interval [X1,XNN].");
        Console.WriteLine("");
        Console.WriteLine("  If the subintervals are equally spaced, the derivative");
        Console.WriteLine("  information has no effect on the result, except for");
        Console.WriteLine("  the first and last values, DN(1) and DN(NN).");
        Console.WriteLine("");
        Console.WriteLine("                            Exact       Computed");
        Console.WriteLine("     X1          XNN        Integral    Integral  Comment");
        Console.WriteLine("");

        double[] fn = new double[nn];
        double[] dn = new double[nn];

        for (test = 1; test <= 5; test++)
        {
            int i;
            switch (test)
            {
                //
                //  Equal spacing.
                //
                case 1:
                {
                    xn = typeMethods.r8vec_even_new(nn, 0.0, Math.PI);
                    for (i = 0; i < nn; i++)
                    {
                        fn[i] = Math.Sin(xn[i]);
                        dn[i] = Math.Cos(xn[i]);
                    }

                    q_exact = -Math.Cos(xn[nn - 1]) + Math.Cos(xn[0]);
                    comment = "Equal spacing, correct DN";
                    break;
                }
                //
                //  Equal spacing, reset DN(2:NN-1) to random numbers.
                //
                case 2:
                {
                    xn = typeMethods.r8vec_even_new(nn, 0.0, Math.PI);
                    for (i = 0; i < nn; i++)
                    {
                        fn[i] = Math.Sin(xn[i]);
                        if (i == 0 || i == nn - 1)
                        {
                            dn[i] = Math.Cos(xn[i]);
                        }
                        else
                        {
                            dn[i] = 1000.0 * UniformRNG.r8_uniform_01(ref seed);
                        }
                    }

                    q_exact = -Math.Cos(xn[nn - 1]) + Math.Cos(xn[0]);
                    comment = "Equal spacing, DN(2:N-1) random";
                    break;
                }
                //
                //  Equal spacing, now reset all of DN to random numbers.
                //
                case 3:
                {
                    xn = typeMethods.r8vec_even_new(nn, 0.0, Math.PI);
                    for (i = 0; i < nn; i++)
                    {
                        fn[i] = Math.Sin(xn[i]);
                        dn[i] = 1000.0 * UniformRNG.r8_uniform_01(ref seed);
                    }

                    q_exact = -Math.Cos(xn[nn - 1]) + Math.Cos(xn[0]);
                    comment = "Equal spacing, DN(1:N) random";
                    break;
                }
                //
                //  Variable spacing, correct data.
                //
                case 4:
                {
                    xn = typeMethods.r8vec_even_new(nn, 0.0, Math.PI * Math.PI);
                    for (i = 0; i < nn; i++)
                    {
                        xn[i] = Math.Sqrt(xn[i]);
                        fn[i] = Math.Sin(xn[i]);
                        dn[i] = Math.Cos(xn[i]);
                    }

                    q_exact = -Math.Cos(xn[nn - 1]) + Math.Cos(xn[0]);
                    comment = "Variable spacing, correct DN";
                    break;
                }
                //
                //  Variable spacing, change one entry in DN.
                //
                case 5:
                {
                    xn = typeMethods.r8vec_even_new(nn, 0.0, Math.PI * Math.PI);
                    for (i = 0; i < nn; i++)
                    {
                        xn[i] = Math.Sqrt(xn[i]);
                        fn[i] = Math.Sin(xn[i]);
                        dn[i] = Math.Cos(xn[i]);
                    }

                    dn[(nn - 1) / 2] = 1000.0 * UniformRNG.r8_uniform_01(ref seed);
                    q_exact = -Math.Cos(xn[nn - 1]) + Math.Cos(xn[0]);
                    comment = "Variable spacing, a single internal DN randomized.";
                    break;
                }
            }

            double q_computed = HermiteCubic.hermite_cubic_spline_integral(nn, xn, fn, dn);

            Console.WriteLine("  " + xn[0].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + xn[nn - 1].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + q_exact.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + q_computed.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + comment + "");
        }
    }

    private static void test11()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST11 tests HERMITE_CUBIC_LAGRANGE_VALUE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 February 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int j;
        const int n = 11;

        Console.WriteLine("");
        Console.WriteLine("TEST11:");
        Console.WriteLine("  HERMITE_CUBIC_LAGRANGE_VALUE evaluates the four");
        Console.WriteLine("  Lagrange basis functions associated with F1, D1,");
        Console.WriteLine("  F2 and D2 such that");
        Console.WriteLine("");
        Console.WriteLine("  P(X) = F1 * LF1(X) + D1 * LD1(X)");
        Console.WriteLine("       + F2 * LF2(X) + D2 * LD2(X).");
        Console.WriteLine("");
        Console.WriteLine("  The first, second and third derivatives of these four");
        Console.WriteLine("  Lagrange basis functions are also computed.");

        const double x1 = 1.0;
        const double x2 = 2.0;
        double[] x = typeMethods.r8vec_even_new(n, 0.0, 2.5);

        double[] f = new double[4 * n];
        double[] d = new double[4 * n];
        double[] s = new double[4 * n];
        double[] t = new double[4 * n];

        HermiteCubic.hermite_cubic_lagrange_value(x1, x2, n, x, ref f, ref d, ref s, ref t);

        Console.WriteLine("");
        Console.WriteLine("  The Lagrange basis functions:");
        Console.WriteLine("");
        Console.WriteLine("     I        X           LF1         LD1         LF2         LD2");
        Console.WriteLine("");
        for (j = 0; j < n; j++)
        {
            Console.WriteLine("  " + j.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + x[j].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + f[0 + j * 4].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + f[1 + j * 4].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + f[2 + j * 4].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + f[3 + j * 4].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("  The derivative of the Lagrange basis functions:");
        Console.WriteLine("");
        Console.WriteLine("     I        X           LF1         LD1         LF2         LD2");
        Console.WriteLine("");
        for (j = 0; j < n; j++)
        {
            Console.WriteLine("  " + j.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + x[j].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + d[0 + j * 4].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + d[1 + j * 4].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + d[2 + j * 4].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + d[3 + j * 4].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }
    }

    private static void test12()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST12 tests HERMITE_CUBIC_LAGRANGE_INTEGRAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 February 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("TEST12:");
        Console.WriteLine("  HERMITE_CUBIC_LAGRANGE_INTEGRAL returns the integrals");
        Console.WriteLine("  of the four Lagrange basis functions associated");
        Console.WriteLine("  with F1, D1, F2 and D2 such that");
        Console.WriteLine("");
        Console.WriteLine("  P(X) = F1 * LF1(X) + D1 * LD1(X)");
        Console.WriteLine("       + F2 * LF2(X) + D2 * LD2(X).");
        Console.WriteLine("");
        Console.WriteLine("  The Lagrange basis function integrals:");
        Console.WriteLine("");
        Console.WriteLine("        X1          X2          LF1         LD1         LF2         LD2");
        Console.WriteLine("");

        const double x2 = 1.0;

        for (i = -6; i <= 2; i++)
        {
            double x1 = i;
            double[] q = HermiteCubic.hermite_cubic_lagrange_integral(x1, x2);
            Console.WriteLine("  " + x1.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + x2.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + q[0].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + q[1].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + q[2].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + q[3].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }
    }

    private static void test13()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST13 tests HERMITE_CUBIC_LAGRANGE_INTEGRATE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 February 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int j;
        double[] p = new double[4];

        Console.WriteLine("");
        Console.WriteLine("TEST13:");
        Console.WriteLine("  HERMITE_CUBIC_LAGRANGE_INTEGRATE integrates a Hermite cubic");
        Console.WriteLine("  Lagrange polynomial from A to B.");
        Console.WriteLine("");
        Console.WriteLine("  Compute each result TWICE:");
        Console.WriteLine("  First row computed using HERMITE_CUBIC_INTEGRATE.");
        Console.WriteLine("  Second row computed using HERMITE_CUBIC_LAGRANGE_INTEGRATE.");

        const double x1 = 0.0;
        const double x2 = 10.0;

        Console.WriteLine("");
        Console.WriteLine("        A           B           LF1         LD1         LF2         LD2");
        Console.WriteLine("");

        const double a = x1 - 1.0;

        for (j = -3; j <= 12; j++)
        {
            double b = ((10 - j) * x1
                        + j * x2)
                       / 10;

            double f1 = 1.0;
            double d1 = 0.0;
            double f2 = 0.0;
            double d2 = 0.0;
            p[0] = HermiteCubic.hermite_cubic_integrate(x1, f1, d1, x2, f2, d2, a, b);

            f1 = 0.0;
            d1 = 1.0;
            f2 = 0.0;
            d2 = 0.0;
            p[1] = HermiteCubic.hermite_cubic_integrate(x1, f1, d1, x2, f2, d2, a, b);

            f1 = 0.0;
            d1 = 0.0;
            f2 = 1.0;
            d2 = 0.0;
            p[2] = HermiteCubic.hermite_cubic_integrate(x1, f1, d1, x2, f2, d2, a, b);

            f1 = 0.0;
            d1 = 0.0;
            f2 = 0.0;
            d2 = 1.0;
            p[3] = HermiteCubic.hermite_cubic_integrate(x1, f1, d1, x2, f2, d2, a, b);

            double[] q = HermiteCubic.hermite_cubic_lagrange_integrate(x1, x2, a, b);

            Console.WriteLine("  " + a.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + b.ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + p[0].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + p[1].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + p[2].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + p[3].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
            Console.WriteLine("  " + "          "
                                   + "  " + "          "
                                   + "  " + q[0].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + q[1].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + q[2].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                   + "  " + q[3].ToString(CultureInfo.InvariantCulture).PadLeft(10) + "");
        }
    }

    private static void test14()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST14 tests HERMITE_CUBIC_SPLINE_QUAD_RULE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 February 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 11;

        double[] dn = new double[N];
        double[] fn = new double[N];
        int k;
        double[] x = new double[N];

        Console.WriteLine("");
        Console.WriteLine("TEST14:");
        Console.WriteLine("  HERMITE_CUBIC_SPLINE_QUAD_RULE returns a quadrature rule");
        Console.WriteLine("  for Hermite cubic splines.");

        int seed = 123456789;

        for (k = 1; k <= 2; k++)
        {
            Console.WriteLine("");
            int i;
            switch (k)
            {
                case 1:
                {
                    Console.WriteLine("  Case 1: Random spacing");
                    double[] r = UniformRNG.r8vec_uniform_01_new(N, ref seed);

                    x[0] = r[0];
                    for (i = 1; i < N; i++)
                    {
                        x[i] = x[i - 1] + r[i];
                    }

                    break;
                }
                case 2:
                {
                    Console.WriteLine("  Case 2: Uniform spacing");
                    Console.WriteLine("  F(2:N-1) have equal weight.");
                    Console.WriteLine("  D(2:N-1) have zero weight.");
                    for (i = 0; i < N; i++)
                    {
                        x[i] = (10 + i) / 20.0;
                    }

                    break;
                }
            }

            double[] w = HermiteCubic.hermite_cubic_spline_quad_rule(N, x);

            Console.WriteLine("");
            Console.WriteLine("   I   J        X         W                Q");
            Console.WriteLine("");

            for (i = 0; i <= 1; i++)
            {
                int j;
                for (j = 0; j < N; j++)
                {
                    int l;
                    for (l = 0; l < N; l++)
                    {
                        fn[l] = 0.0;
                        dn[l] = 0.0;
                    }

                    switch (i)
                    {
                        case 0:
                            fn[j] = 1.0;
                            break;
                        default:
                            dn[j] = 1.0;
                            break;
                    }

                    double q = HermiteCubic.hermite_cubic_spline_integral(N, x, fn, dn);

                    Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                           + "  " + j.ToString(CultureInfo.InvariantCulture).PadLeft(2)
                                           + "  " + x[j].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                           + "  " + q.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                           + "  " + w[i + j * 2].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
                }
            }
        }
    }

    private static void test15()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST15 tests HERMITE_CUBIC_SPLINE_QUAD_RULE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 March 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int N = 11;

        double dn = 0;
        double fn = 0;
        int j;
        double s = 0;
        double t = 0;
        double[] x = new double[N];

        Console.WriteLine("");
        Console.WriteLine("TEST15:");
        Console.WriteLine("  HERMITE_CUBIC_SPLINE_QUAD_RULE returns a quadrature rule");
        Console.WriteLine("  for Hermite cubic splines.");

        int seed = 123456789;

        double[] r = UniformRNG.r8vec_uniform_01_new(N, ref seed);

        x[0] = r[0];
        for (j = 1; j < N; j++)
        {
            x[j] = x[j - 1] + r[j];
        }
            
        Console.WriteLine("");
        Console.WriteLine("  Random spacing");
        Console.WriteLine("  Number of points N = " + N);
        Console.WriteLine("  Interval = [" + x[0] + ", " + x[N - 1] + "]");

        double[] w = HermiteCubic.hermite_cubic_spline_quad_rule(N, x);

        double q = 0.0;

        for (j = 0; j < N; j++)
        {
            cubic_value(x[j], ref fn, ref dn, ref s, ref t);
            q = q + w[0 + j * 2] * fn + w[1 + j * 2] * dn;
        }

        double q_exact = cubic_integrate(x[0], x[N - 1]);

        Console.WriteLine("");
        Console.WriteLine("  Q         = " + q + "");
        Console.WriteLine("  Q (exact) = " + q_exact + "");
    }

    private static double cubic_antiderivative(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CUBIC_ANTIDERIVATIVE evaluates the antiderivative function of a cubic.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 January 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument.
        //
        //    Output, double CUBIC_ANTIDERIVATIVE, the value.
        //
    {
        double value = x * x * (5.0 + x * (-7.0 / 3.0 + x * 1.0 / 4.0));

        return value;
    }

    private static double cubic_integrate(double a, double b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CUBIC_INTEGRATE integrates the cubic from A to B.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 February 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, the integration interval.
        //
        //    Output, double Q, the integral from A to B.
        //
    {
        double q = cubic_antiderivative(b) - cubic_antiderivative(a);

        return q;
    }

    private static void cubic_value(double x, ref double f, ref double d, ref double s, ref double t)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CUBIC_VALUE evaluates a cubic function.
        //
        //  Discussion:
        //
        //    f(x) =   x^3 -  7 x^2 + 10 x
        //    d(x) = 3 x^2 - 14 x   + 10
        //    s(x) = 6 x   - 14
        //    t(x) = 6
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 February 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument.
        //
        //    Output, double F, D, S, T, the value and first three
        //    derivatives of the cubic function.
        //
    {
        f = 0.0 + x * (10.0 + x * (-7.0 + x * 1.0));
        d = 10.0 + x * (-14.0 + x * 3.0);
        s = -14.0 + x * 6.0;
        t = 6.0;
    }
}