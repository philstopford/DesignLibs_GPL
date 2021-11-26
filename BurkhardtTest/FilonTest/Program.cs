﻿using System;
using System.Globalization;
using Burkardt.IntegralNS;

namespace FilonTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for FILON_TEST.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("FILON_TEST");
        Console.WriteLine("  Test the FILON library.");

        test01();
        test02();
        test03();
        test04();
        test05();
        test06();

        Console.WriteLine("");
        Console.WriteLine("FILON_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests FILON_TAB_COS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double exact = 0;
        double[] ftab = new double[1];
        int i;
        int k;
        double t = 0;

        const double a = 0.0;
        const double b = 2.0 * Math.PI;

        int n = 11;
        //
        //  Set the X values.
        //
        double[] x = new double[n];
        for (i = 0; i < n; i++)
        {
            x[i] = ((n - i - 1) * a
                    + i * b)
                   / (n - 1);
        }

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  FILON_TAB_COS estimates the integral of.");
        Console.WriteLine("  F(X) * COS ( T * X )");
        Console.WriteLine("  Use integrands F(X)=1, X, X^2.");
        Console.WriteLine("");
        Console.WriteLine("  A = " + a + "");
        Console.WriteLine("  B = " + b + "");
        Console.WriteLine("  N = " + n + "");
        Console.WriteLine("");
        Console.WriteLine("       T                      Approximate             Exact");

        for (k = 1; k <= 3; k++)
        {
            t = k switch
            {
                1 => 1.0,
                2 => 2.0,
                3 => 10.0,
                _ => t
            };

            Console.WriteLine("");

            for (i = 1; i <= 3; i++)
            {
                ftab = i switch
                {
                    1 => zero_integrand(n, x),
                    2 => one_integrand(n, x),
                    3 => two_integrand(n, x),
                    _ => ftab
                };

                double result = Filon.filon_tab_cos(n, ftab, a, b, t);

                exact = i switch
                {
                    1 => (Math.Sin(t * b) - Math.Sin(t * a)) / t,
                    2 => (Math.Cos(t * b) + t * b * Math.Sin(t * b) - (Math.Cos(t * a) + t * a * Math.Sin(t * a))) / t /
                         t,
                    3 => (2.0 * t * b * Math.Cos(t * b) + (t * t * b * b - 2.0) * Math.Sin(t * b) -
                          (2.0 * t * a * Math.Cos(t * a) + (t * t * a * a - 2.0) * Math.Sin(t * a))) / t / t / t,
                    _ => exact
                };

                Console.WriteLine(t.ToString(CultureInfo.InvariantCulture).PadLeft(24) + "  "
                                                           + result.ToString(CultureInfo.InvariantCulture).PadLeft(24) + "  "
                                                           + exact.ToString(CultureInfo.InvariantCulture).PadLeft(24) + "");
            }
        }
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests FILON_TAB_COS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int j;
        //
        //  Example suggested by James Roedder.
        //
        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  Integrate F(X) = log(1+X)*Math.Cos(T*X):");
        Console.WriteLine("  Supply integrand as a table.");
        Console.WriteLine("  T = 10, and N increases");
        Console.WriteLine("");
        Console.WriteLine("       N    Approximate             Exact                   Error");
        Console.WriteLine("");

        const double a = 0.0;
        const double b = 2.0 * Math.PI;

        for (j = 1; j <= 6; j++)
        {
            int n = (int)Math.Pow(2, j) * 10 + 1;
            //
            //  Set the X values.
            //
            double[] x = new double[n];
            int i;
            for (i = 0; i < n; i++)
            {
                x[i] = ((n - i - 1) * a
                        + i * b)
                       / (n - 1);
            }

            double[] ftab = log_integrand(n, x);

            double t = 10.0;

            double result = Filon.filon_tab_cos(n, ftab, a, b, t);

            double exact = -0.008446594405;
            double error = result - exact;

            Console.WriteLine(n.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                                                      + result.ToString(CultureInfo.InvariantCulture).PadLeft(24) + "  "
                                                      + exact.ToString(CultureInfo.InvariantCulture).PadLeft(24) + "  "
                                                      + error.ToString(CultureInfo.InvariantCulture).PadLeft(24) + "");

        }
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests FILON_FUN_COS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int j;
        //
        //  Example suggested by James Roedder.
        //
        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  Integrate F(X)=log(1+X)*Math.Cos(T*X):");
        Console.WriteLine("  Supply integrand as a function.");
        Console.WriteLine("  T = 10, and N increases");
        Console.WriteLine("");
        Console.WriteLine("       N    Approximate             Exact                   Error");
        Console.WriteLine("");

        double a = 0.0;
        double b = 2.0 * Math.PI;

        for (j = 1; j <= 6; j++)
        {
            int n = (int)Math.Pow(2, j) * 10 + 1;

            double t = 10.0;

            double result = Filon.filon_fun_cos(n, log_integrand, a, b, t);

            double exact = -0.008446594405;
            double error = result - exact;

            Console.WriteLine(n.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                                                      + result.ToString(CultureInfo.InvariantCulture).PadLeft(24) + "  "
                                                      + exact.ToString(CultureInfo.InvariantCulture).PadLeft(24) + "  "
                                                      + error.ToString(CultureInfo.InvariantCulture).PadLeft(24) + "");
        }
    }

    private static void test04()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 tests FILON_TAB_SIN.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double exact = 0;
        double[] ftab = new double[1];
        int i;
        int k;
        double t = 0;

        const double a = 0.0;
        const double b = 2.0 * Math.PI;
        const int n = 11;
        //
        //  Set the X values.
        //
        double[] x = new double[n];
        for (i = 0; i < n; i++)
        {
            x[i] = ((n - i - 1) * a
                    + i * b)
                   / (n - 1);
        }

        Console.WriteLine("");
        Console.WriteLine("TEST04");
        Console.WriteLine("  FILON_TAB_SIN estimates the integral of.");
        Console.WriteLine("  F(X) * SIN ( T * X )");
        Console.WriteLine("  Use integrands 1, X, X^2.");
        Console.WriteLine("");
        Console.WriteLine("  A = " + a + "");
        Console.WriteLine("  B = " + b + "");
        Console.WriteLine("  N = " + n + "");
        Console.WriteLine("");
        Console.WriteLine("       T                      Approximate             Exact");
        Console.WriteLine("");

        for (k = 1; k <= 3; k++)
        {
            t = k switch
            {
                1 => 1.0,
                2 => 2.0,
                3 => 10.0,
                _ => t
            };

            Console.WriteLine("");

            for (i = 1; i <= 3; i++)
            {
                ftab = i switch
                {
                    1 => zero_integrand(n, x),
                    2 => one_integrand(n, x),
                    3 => two_integrand(n, x),
                    _ => ftab
                };

                double result = Filon.filon_tab_sin(n, ftab, a, b, t);

                exact = i switch
                {
                    1 => (-Math.Cos(t * b) + Math.Cos(t * a)) / t,
                    2 => (Math.Sin(t * b) - t * b * Math.Cos(t * b) - (Math.Sin(t * a) - t * a * Math.Cos(t * a))) / t /
                         t,
                    3 => (2.0 * t * b * Math.Sin(t * b) + (2.0 - t * t * b * b) * Math.Cos(t * b) -
                          (2.0 * t * a * Math.Sin(t * a) + (2.0 - t * t * a * a) * Math.Cos(t * a))) / t / t / t,
                    _ => exact
                };

                Console.WriteLine(t.ToString(CultureInfo.InvariantCulture).PadLeft(24) + "  "
                                                           + result.ToString(CultureInfo.InvariantCulture).PadLeft(24) + "  "
                                                           + exact.ToString(CultureInfo.InvariantCulture).PadLeft(24) + "");
            }
        }
    }

    private static void test05()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 tests FILON_TAB_COS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int j;
        //
        //  Example suggested by James Roedder.
        //
        Console.WriteLine("");
        Console.WriteLine("TEST05");
        Console.WriteLine("  Integrate F(X)=log(1+X)*Math.Sin(T*X):");
        Console.WriteLine("  Supply integrand as a table.");
        Console.WriteLine("  T = 10, and N increases");
        Console.WriteLine("");
        Console.WriteLine("       N    Approximate             Exact                   Error");
        Console.WriteLine("");

        const double a = 0.0;
        const double b = 2.0 * Math.PI;

        for (j = 1; j <= 6; j++)
        {
            int n = (int)Math.Pow(2, j) * 10 + 1;
            //
            //  Set the X values.
            //
            double[] x = new double[n];
            int i;
            for (i = 0; i < n; i++)
            {
                x[i] = ((n - i - 1) * a
                        + i * b)
                       / (n - 1);
            }

            double[] ftab = log_integrand(n, x);

            double t = 10.0;

            double result = Filon.filon_tab_sin(n, ftab, a, b, t);

            double exact = -0.19762680771872;
            double error = result - exact;

            Console.WriteLine(n.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                                                      + result.ToString(CultureInfo.InvariantCulture).PadLeft(24) + "  "
                                                      + exact.ToString(CultureInfo.InvariantCulture).PadLeft(24) + "  "
                                                      + error.ToString(CultureInfo.InvariantCulture).PadLeft(24) + "");

        }
    }

    private static void test06()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST06 tests FILON_FUN_COS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int j;
        //
        //  Example suggested by James Roedder.
        //
        Console.WriteLine("");
        Console.WriteLine("TEST06");
        Console.WriteLine("  Integrate F(X)=log(1+X)*Math.Sin(T*X):");
        Console.WriteLine("  Supply integrand as a function.");
        Console.WriteLine("  T = 10, and N increases");
        Console.WriteLine("");
        Console.WriteLine("       N    Approximate             Exact                   Error");
        Console.WriteLine("");

        double a = 0.0;
        double b = 2.0 * Math.PI;

        for (j = 1; j <= 6; j++)
        {
            int n = (int)Math.Pow(2, j) * 10 + 1;

            double t = 10.0;

            double result = Filon.filon_fun_sin(n, log_integrand, a, b, t);

            double exact = -0.19762680771872;
            double error = result - exact;

            Console.WriteLine(n.ToString(CultureInfo.InvariantCulture).PadLeft(6) + "  "
                                                      + result.ToString(CultureInfo.InvariantCulture).PadLeft(24) + "  "
                                                      + exact.ToString(CultureInfo.InvariantCulture).PadLeft(24) + "  "
                                                      + error.ToString(CultureInfo.InvariantCulture).PadLeft(24) + "");
        }
    }

    private static double[] zero_integrand(int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ZERO_INTEGRAND evaluates the integrand x^0.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of points.
        //
        //    Input, double X[N], the evaluation points.
        //
        //    Output, double ZERO_INTEGRAND[N], the function values.
        //
    {
        int i;

        double[] fx = new double[n];

        for (i = 0; i < n; i++)
        {
            fx[i] = 1.0;
        }

        return fx;
    }

    private static double[] one_integrand(int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ONE_INTEGRAND evaluates the integrand X.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of points.
        //
        //    Input, double X[N], the evaluation points.
        //
        //    Output, double ONE_INTEGRAND[N], the function values.
        //
    {
        int i;

        double[] fx = new double[n];

        for (i = 0; i < n; i++)
        {
            fx[i] = x[i];
        }

        return fx;
    }

    private static double[] two_integrand(int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TWO_INTEGRAND evaluates the integrand X^2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of points.
        //
        //    Input, double X[N], the evaluation points.
        //
        //    Output, double TWO_INTEGRAND[N], the function values.
        //
    {
        int i;

        double[] fx = new double[n];

        for (i = 0; i < n; i++)
        {
            fx[i] = x[i] * x[i];
        }

        return fx;
    }

    private static double[] log_integrand(int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LOG_INTEGRAND evaluates the logarithmic integrand.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of points.
        //
        //    Input, double X[N], the evaluation points.
        //
        //    Output, double LOG_INTEGRAND[N], the function values.
        //
    {
        int i;

        double[] fx = new double[n];

        for (i = 0; i < n; i++)
        {
            fx[i] = Math.Log(1.0 + x[i]);
        }

        return fx;
    }
}