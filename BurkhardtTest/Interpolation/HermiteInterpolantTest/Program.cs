using System;
using Burkardt.PolynomialNS;
using Burkardt.Types;
using Burkardt.Uniform;
using Hermite = Burkardt.Interpolation.Hermite;

namespace HermiteInterpolantTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for hermite_interpolant_test.
        //
        //  Discussion:
        //
        //    hermite_interpolant_test tests hermite_interpolant.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 November 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("hermite_interpolant_test");
        Console.WriteLine("  Test hermite_interpolant.");

        test01();
        test02();
        test03();
        test04();
        test05();
        test06();
        test07();
        test08();

        Console.WriteLine("");
        Console.WriteLine("hermite_interpolant_test");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 uses f(x) = 1 + 2x + 3x^2 at x = 0, 1, 2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 May 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 3;

        int n = N;
        double[] x = {0.0, 1.0, 2.0};
        double[] y = {1.0, 6.0, 17.0};
        double[] yp = {2.0, 8.0, 14.0};

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  HERMITE computes the Hermite interpolant to data.");
        Console.WriteLine("  Here, f(x) = 1 + 2x + 3x^2.");

        Hermite.hermite_demo(n, x, y, yp);

    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 uses f(x) = 6 + 5x + 4x^2 + 3x^3 + 2x^4 + x^5 at x = 0, 1, 2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 May 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 3;

        int i;
        int n = N;
        double[] x;
        double[] y;
        double[] yp;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  HERMITE computes the Hermite interpolant to data.");
        Console.WriteLine("  Here, f(x) = 6 + 5x + 4x^2 + 3x^3 + 2x^4 + x^5.");

        x = new double[n];
        y = new double[n];
        yp = new double[n];

        for (i = 0; i < n; i++)
        {
            x[i] = i;

            y[i] = 6.0 + x[i] * (
                5.0 + x[i] * (
                    4.0 + x[i] * (
                        3.0 + x[i] * (
                            2.0 + x[i]))));

            yp[i] = 5.0 + x[i] * (
                8.0 + x[i] * (
                    9.0 + x[i] * (
                        8.0 + x[i] *
                        5.0)));
        }

        Hermite.hermite_demo(n, x, y, yp);

    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 uses f(x) = r1 + r2x + r3x^2 + r4x^3 + r5x^4 + r6x^5 at x = r7 r8 r9
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 May 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int N = 3;

        double[] c;
        int i;
        int n = N;
        int seed;
        double[] x;
        double[] y;
        double[] yp;

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  HERMITE computes the Hermite interpolant to data.");
        Console.WriteLine("  Here, f(x) is a fifth order polynomial with random");
        Console.WriteLine("  coefficients, and the abscissas are random.");

        c = new double[2 * n];
        x = new double[n];
        y = new double[n];
        yp = new double[n];

        seed = 123456789;

        UniformRNG.r8vec_uniform_01(n, ref seed, ref x);
        typeMethods.r8vec_print(n, x, "  Random abscissas");

        UniformRNG.r8vec_uniform_01(2 * n, ref seed, ref c);
        typeMethods.r8vec_print(2 * n, c, "  Random polynomial coefficients.");

        for (i = 0; i < n; i++)
        {
            y[i] = c[0] + x[i] * (
                c[1] + x[i] * (
                    c[2] + x[i] * (
                        c[3] + x[i] * (
                            c[4] + x[i] * c[5]))));

            yp[i] = c[1] + x[i] * (
                c[2] * 2.0 + x[i] * (
                    c[3] * 3.0 + x[i] * (
                        c[4] * 4.0 + x[i] *
                        c[5] * 5.0)));
        }

        Hermite.hermite_demo(n, x, y, yp);
    }

    private static void test04()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 interpolates the Runge function using equally spaced data.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        double max_dif;
        int n;
        int nd;
        int ndp;
        int ns;
        double[] x;
        double[] xd;
        double[] xdp;
        double[] xs;
        double xhi;
        double xlo;
        double xt;
        double[] y;
        double[] yd;
        double[] ydp;
        double[] yp;
        double[] ys;
        double yt;

        Console.WriteLine("");
        Console.WriteLine("TEST04");
        Console.WriteLine("  HERMITE computes the Hermite interpolant to data.");
        Console.WriteLine("  Here, f(x) is the Runge function");
        Console.WriteLine("  and the data is evaluated at equally spaced points.");
        Console.WriteLine("  As N increases, the maximum error grows.");
        Console.WriteLine("");
        Console.WriteLine("     N     Max | F(X) - H(F(X)) |");
        Console.WriteLine("");

        for (n = 3; n <= 15; n += 2)
        {
            y = new double[n];
            yp = new double[n];

            nd = 2 * n;
            xd = new double[nd];
            yd = new double[nd];

            ndp = 2 * n - 1;
            xdp = new double[ndp];
            ydp = new double[ndp];

            ns = 10 * (n - 1) + 1;

            xlo = -5.0;
            xhi = +5.0;
            x = typeMethods.r8vec_linspace_new(n, xlo, xhi);

            for (i = 0; i < n; i++)
            {
                y[i] = 1.0 / (1.0 + x[i] * x[i]);
                yp[i] = -2.0 * x[i] / (1.0 + x[i] * x[i]) / (1.0 + x[i] * x[i]);
            }

            Hermite.hermite_interpolant(n, x, y, yp, ref xd, ref yd, ref xdp, ref ydp);
            //
            //  Compare exact and interpolant at sample points.
            //
            xs = typeMethods.r8vec_linspace_new(ns, xlo, xhi);

            ys = Dif.dif_vals(nd, xd, yd, ns, xs);

            max_dif = 0.0;
            for (i = 0; i < ns; i++)
            {
                xt = xs[i];
                yt = 1.0 / (1.0 + xt * xt);
                max_dif = Math.Max(max_dif, Math.Abs(ys[i] - yt));
            }

            Console.WriteLine("  " + n.ToString().PadLeft(4)
                                   + "  " + max_dif.ToString().PadLeft(14) + "");

        }

    }

    private static void test05()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 interpolates the Runge function using Chebyshev spaced data.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        double max_dif;
        int n;
        int nd;
        int ndp;
        int ns;
        double[] x;
        double[] xd;
        double[] xdp;
        double[] xs;
        double xhi;
        double xlo;
        double xt;
        double[] y;
        double[] yd;
        double[] ydp;
        double[] yp;
        double[] ys;
        double yt;

        Console.WriteLine("");
        Console.WriteLine("TEST05");
        Console.WriteLine("  HERMITE computes the Hermite interpolant to data.");
        Console.WriteLine("  Here, f(x) is the Runge function");
        Console.WriteLine("  and the data is evaluated at Chebyshev spaced points.");
        Console.WriteLine("  As N increases, the maximum error decreases.");
        Console.WriteLine("");
        Console.WriteLine("     N     Max | F(X) - H(F(X)) |");
        Console.WriteLine("");

        for (n = 3; n <= 15; n += 2)
        {
            y = new double[n];
            yp = new double[n];

            nd = 2 * n;
            xd = new double[nd];
            yd = new double[nd];

            ndp = 2 * n - 1;
            xdp = new double[ndp];
            ydp = new double[ndp];

            ns = 10 * (n - 1) + 1;

            xlo = -5.0;
            xhi = +5.0;
            x = typeMethods.r8vec_chebyshev_new(n, xlo, xhi);

            for (i = 0; i < n; i++)
            {
                y[i] = 1.0 / (1.0 + x[i] * x[i]);
                yp[i] = -2.0 * x[i] / (1.0 + x[i] * x[i]) / (1.0 + x[i] * x[i]);
            }

            Hermite.hermite_interpolant(n, x, y, yp, ref xd, ref yd, ref xdp, ref ydp);
            //
            //  Compare exact and interpolant at sample points.
            //
            xs = typeMethods.r8vec_linspace_new(ns, xlo, xhi);

            ys = Dif.dif_vals(nd, xd, yd, ns, xs);

            max_dif = 0.0;
            for (i = 0; i < ns; i++)
            {
                xt = xs[i];
                yt = 1.0 / (1.0 + xt * xt);
                max_dif = Math.Max(max_dif, Math.Abs(ys[i] - yt));
            }

            Console.WriteLine("  " + n.ToString().PadLeft(4)
                                   + "  " + max_dif.ToString().PadLeft(14) + "");

        }
    }

    private static void test06()

        //***************************************************************************80
        //
        //  Purpose:
        //
        //    TEST06 tests HERMITE_BASIS_0 and HERMITE_BASIS_1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 May 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int ND = 2;

        double f01;
        double f02;
        double f11;
        double f12;
        int j;
        int nd = ND;
        double[] xd = {0.0, 10.0};
        double xv;
        double[] yd = new double[ND];
        double yh;
        double[] ypd = new double[ND];
        double yv;

        Console.WriteLine("");
        Console.WriteLine("TEST06:");
        Console.WriteLine("  HERMITE_BASIS_0 and HERMITE_BASIS_1 evaluate the");
        Console.WriteLine("  Hermite global polynomial basis functions");
        Console.WriteLine("  of type 0: associated with function values, and");
        Console.WriteLine("  of type 1: associated with derivative values.");
        //
        //  Let y = x^3 + x^2 + x + 1,
        //  and compute the Hermite global polynomial interpolant based on two 
        //  abscissas:
        //
        for (j = 0; j < nd; j++)
        {
            yd[j] = Math.Pow(xd[j], 3) + Math.Pow(xd[j], 2) + xd[j] + 1.0;
            ypd[j] = 3.0 * Math.Pow(xd[j], 2) + 2.0 * xd[j] + 1.0;
        }

        Console.WriteLine("");
        Console.WriteLine("  Interpolate y = x^3 + x^2 + x + 1.");
        Console.WriteLine("");
        Console.WriteLine("     XD         Y(XD)      Y'(XD)");
        Console.WriteLine("");
        for (j = 0; j < nd; j++)
        {
            Console.WriteLine("  " + xd[j].ToString().PadLeft(10)
                                   + "  " + yd[j].ToString().PadLeft(10)
                                   + "  " + ypd[j].ToString().PadLeft(10) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("     XV         Y(XV)      H(XV)");
        Console.WriteLine("");

        for (j = 0; j <= 10; j++)
        {
            xv = j;

            yv = Math.Pow(xv, 3) + Math.Pow(xv, 2) + xv + 1.0;

            f01 = Hermite.hermite_basis_0(2, xd, 0, xv);
            f11 = Hermite.hermite_basis_1(2, xd, 0, xv);
            f02 = Hermite.hermite_basis_0(2, xd, 1, xv);
            f12 = Hermite.hermite_basis_1(2, xd, 1, xv);

            yh = yd[0] * f01 + ypd[0] * f11 + yd[1] * f02 + ypd[1] * f12;

            Console.WriteLine("  " + xv.ToString().PadLeft(10)
                                   + "  " + yv.ToString().PadLeft(10)
                                   + "  " + yh.ToString().PadLeft(10) + "");
        }
    }

    private static void test07()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST07 tests HERMITE_INTERPOLANT_RULE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 June 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a;
        double b;
        int i;
        int k;
        int n;
        double q;
        double[] w;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TEST07:");
        Console.WriteLine("  HERMITE_INTERPOLANT_RULE");
        Console.WriteLine("  is given a set of N abscissas for a Hermite interpolant");
        Console.WriteLine("  and returns N pairs of quadrature weights");
        Console.WriteLine("  for function and derivative values at the abscissas.");

        n = 3;
        a = 0.0;
        b = 10.0;
        x = typeMethods.r8vec_linspace_new(n, a, b);
        w = Hermite.hermite_interpolant_rule(n, a, b, x);

        Console.WriteLine("");
        Console.WriteLine("     I       X               W(F(X))        W(F'(X))");
        Console.WriteLine("");
        k = 0;
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString().PadLeft(4)
                                   + "  " + x[i].ToString().PadLeft(14)
                                   + "  " + w[k].ToString().PadLeft(14)
                                   + "  " + w[k + 1].ToString().PadLeft(14) + "");
            k += 2;
        }

        Console.WriteLine("");
        Console.WriteLine("  Use the quadrature rule over interval " + a + " to " + b + "");
        Console.WriteLine("");

        q = 0.0;
        k = 0;
        for (i = 0; i < n; i++)
        {
            q = q + w[k] * 1 + w[k + 1] * 0.0;
            k += 2;
        }

        Console.WriteLine("  Estimate integral of 1 = " + q + "");

        q = 0.0;
        k = 0;
        for (i = 0; i < n; i++)
        {
            q = q + w[k] * x[i] + w[k + 1] * 1.0;
            k += 2;
        }

        Console.WriteLine("  Estimate integral of X = " + q + "");

        q = 0.0;
        k = 0;
        for (i = 0; i < n; i++)
        {
            q = q + w[k] * x[i] * x[i] + w[k + 1] * 2.0 * x[i];
            k += 2;
        }

        Console.WriteLine("  Estimate integral of X^2 = " + q + "");

        q = 0.0;
        k = 0;
        for (i = 0; i < n; i++)
        {
            q = q + w[k] / (1.0 + x[i] * x[i])
                - w[k + 1] * 2.0 * x[i] / Math.Pow(1.0 + x[i] * x[i], 2);
            k += 2;
        }

        Console.WriteLine("  Estimate integral of 1/(1+x^2) = " + q + "");

        n = 3;
        a = 0.0;
        b = 1.0;
        x = typeMethods.r8vec_linspace_new(n, a, b);
        w = Hermite.hermite_interpolant_rule(n, a, b, x);

        Console.WriteLine("");
        Console.WriteLine("     I       X               W(F(X))        W(F'(X))");
        Console.WriteLine("");
        k = 0;
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString().PadLeft(4)
                                   + "  " + x[i].ToString().PadLeft(14)
                                   + "  " + w[k].ToString().PadLeft(14)
                                   + "  " + w[k + 1].ToString().PadLeft(14) + "");
            k += 2;
        }

        Console.WriteLine("");
        Console.WriteLine("  Use the quadrature rule over interval " + a + " to " + b + "");
        Console.WriteLine("");

        q = 0.0;
        k = 0;
        for (i = 0; i < n; i++)
        {
            q = q + w[k] * 1 + w[k + 1] * 0.0;
            k += 2;
        }

        Console.WriteLine("  Estimate integral of 1 = " + q + "");

        q = 0.0;
        k = 0;
        for (i = 0; i < n; i++)
        {
            q = q + w[k] * x[i] + w[k + 1] * 1.0;
            k += 2;
        }

        Console.WriteLine("  Estimate integral of X = " + q + "");

        q = 0.0;
        k = 0;
        for (i = 0; i < n; i++)
        {
            q = q + w[k] * x[i] * x[i] + w[k + 1] * 2.0 * x[i];
            k += 2;
        }

        Console.WriteLine("  Estimate integral of X^2 = " + q + "");

        q = 0.0;
        k = 0;
        for (i = 0; i < n; i++)
        {
            q = q + w[k] / (1.0 + x[i] * x[i])
                - w[k + 1] * 2.0 * x[i] / Math.Pow(1.0 + x[i] * x[i], 2);
            k += 2;
        }

        Console.WriteLine("  Estimate integral of 1/(1+x^2) = " + q + "");

        n = 11;
        a = 0.0;
        b = 10.0;
        x = typeMethods.r8vec_linspace_new(n, a, b);
        w = Hermite.hermite_interpolant_rule(n, a, b, x);

        Console.WriteLine("");
        Console.WriteLine("     I       X               W(F(X))        W(F'(X))");
        Console.WriteLine("");
        k = 0;
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString().PadLeft(4)
                                   + "  " + x[i].ToString().PadLeft(14)
                                   + "  " + w[k].ToString().PadLeft(14)
                                   + "  " + w[k + 1].ToString().PadLeft(14) + "");
            k += 2;
        }

        Console.WriteLine("");
        Console.WriteLine("  Use the quadrature rule over interval " + a + " to " + b + "");
        Console.WriteLine("");

        q = 0.0;
        k = 0;
        for (i = 0; i < n; i++)
        {
            q = q + w[k] * 1 + w[k + 1] * 0.0;
            k += 2;
        }

        Console.WriteLine("  Estimate integral of 1 = " + q + "");

        q = 0.0;
        k = 0;
        for (i = 0; i < n; i++)
        {
            q = q + w[k] * x[i] + w[k + 1] * 1.0;
            k += 2;
        }

        Console.WriteLine("  Estimate integral of X = " + q + "");

        q = 0.0;
        k = 0;
        for (i = 0; i < n; i++)
        {
            q = q + w[k] * x[i] * x[i] + w[k + 1] * 2.0 * x[i];
            k += 2;
        }

        Console.WriteLine("  Estimate integral of X^2 = " + q + "");

        q = 0.0;
        k = 0;
        for (i = 0; i < n; i++)
        {
            q = q + w[k] / (1.0 + x[i] * x[i])
                - w[k + 1] * 2.0 * x[i] / Math.Pow(1.0 + x[i] * x[i], 2);
            k += 2;
        }

        Console.WriteLine("  Estimate integral of 1/(1+x^2) = " + q + "");

        n = 11;
        a = 0.0;
        b = 1.0;
        x = typeMethods.r8vec_linspace_new(n, a, b);
        w = Hermite.hermite_interpolant_rule(n, a, b, x);

        Console.WriteLine("");
        Console.WriteLine("     I       X               W(F(X))        W(F'(X))");
        Console.WriteLine("");
        k = 0;
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString().PadLeft(4)
                                   + "  " + x[i].ToString().PadLeft(14)
                                   + "  " + w[k].ToString().PadLeft(14)
                                   + "  " + w[k + 1].ToString().PadLeft(14) + "");
            k += 2;
        }

        Console.WriteLine("");
        Console.WriteLine("  Use the quadrature rule over interval " + a + " to " + b + "");
        Console.WriteLine("");

        q = 0.0;
        k = 0;
        for (i = 0; i < n; i++)
        {
            q = q + w[k] * 1 + w[k + 1] * 0.0;
            k += 2;
        }

        Console.WriteLine("  Estimate integral of 1 = " + q + "");

        q = 0.0;
        k = 0;
        for (i = 0; i < n; i++)
        {
            q = q + w[k] * x[i] + w[k + 1] * 1.0;
            k += 2;
        }

        Console.WriteLine("  Estimate integral of X = " + q + "");

        q = 0.0;
        k = 0;
        for (i = 0; i < n; i++)
        {
            q = q + w[k] * x[i] * x[i] + w[k + 1] * 2.0 * x[i];
            k += 2;
        }

        Console.WriteLine("  Estimate integral of X^2 = " + q + "");

        q = 0.0;
        k = 0;
        for (i = 0; i < n; i++)
        {
            q = q + w[k] / (1.0 + x[i] * x[i])
                - w[k + 1] * 2.0 * x[i] / Math.Pow(1.0 + x[i] * x[i], 2);
            k += 2;
        }

        Console.WriteLine("  Estimate integral of 1/(1+x^2) = " + q + "");

        n = 11;
        a = 0.0;
        b = 1.0;
        x = typeMethods.r8vec_chebyshev_new(n, a, b);
        w = Hermite.hermite_interpolant_rule(n, a, b, x);

        Console.WriteLine("");
        Console.WriteLine("     I       X               W(F(X))        W(F'(X))");
        Console.WriteLine("");
        k = 0;
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString().PadLeft(4)
                                   + "  " + x[i].ToString().PadLeft(14)
                                   + "  " + w[k].ToString().PadLeft(14)
                                   + "  " + w[k + 1].ToString().PadLeft(14) + "");
            k += 2;
        }

        Console.WriteLine("");
        Console.WriteLine("  Use the quadrature rule over interval " + a + " to " + b + "");
        Console.WriteLine("");

        q = 0.0;
        k = 0;
        for (i = 0; i < n; i++)
        {
            q = q + w[k] * 1 + w[k + 1] * 0.0;
            k += 2;
        }

        Console.WriteLine("  Estimate integral of 1 = " + q + "");

        q = 0.0;
        k = 0;
        for (i = 0; i < n; i++)
        {
            q = q + w[k] * x[i] + w[k + 1] * 1.0;
            k += 2;
        }

        Console.WriteLine("  Estimate integral of X = " + q + "");

        q = 0.0;
        k = 0;
        for (i = 0; i < n; i++)
        {
            q = q + w[k] * x[i] * x[i] + w[k + 1] * 2.0 * x[i];
            k += 2;
        }

        Console.WriteLine("  Estimate integral of X^2 = " + q + "");

        q = 0.0;
        k = 0;
        for (i = 0; i < n; i++)
        {
            q = q + w[k] * Math.Sin(x[i]) + w[k + 1] * Math.Cos(x[i]);
            k += 2;
        }

        Console.WriteLine("  Estimate integral of SIN(X) = " + q + "");
    }

    private static void test08()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST08 tabulates the interpolant and its derivative. 
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 November 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int n;
        int nd;
        int ndp;
        int ns;
        double[] x;
        double[] xd;
        double[] xdp;
        double[] xs;
        double[] y;
        double[] yd;
        double[] ydp;
        double[] yp;
        double[] ys;
        double[] ysp;

        Console.WriteLine("");
        Console.WriteLine("TEST08");
        Console.WriteLine("  HERMITE_INTERPOLANT sets up the Hermite interpolant.");
        Console.WriteLine("  HERMITE_INTERPOLANT_VALUE evaluates it.");
        Console.WriteLine("  Consider data for y=sin(x) at x=0,1,2,3,4.");

        n = 5;
        y = new double[n];
        yp = new double[n];

        nd = 2 * n;
        xd = new double[nd];
        yd = new double[nd];

        ndp = 2 * n - 1;
        xdp = new double[ndp];
        ydp = new double[ndp];

        x = typeMethods.r8vec_linspace_new(n, 0.0, 4.0);
        for (i = 0; i < n; i++)
        {
            y[i] = Math.Sin(x[i]);
            yp[i] = Math.Cos(x[i]);
        }

        Hermite.hermite_interpolant(n, x, y, yp, ref xd, ref yd, ref xdp, ref ydp);
        /*
        Now sample the interpolant at NS points, which include data values.
        */
        ns = 4 * (n - 1) + 1;
        ys = new double[ns];
        ysp = new double[ns];

        xs = typeMethods.r8vec_linspace_new(ns, 0.0, 4.0);

        Hermite.hermite_interpolant_value(nd, xd, yd, xdp, ydp, ns, xs, ref ys, ref ysp);

        Console.WriteLine("");
        Console.WriteLine("  In the following table, there should be perfect");
        Console.WriteLine("  agreement between F and H, and F' and H'");
        Console.WriteLine("  at the data points X = 0, 1, 2, 3, and 4.");
        Console.WriteLine("");
        Console.WriteLine("  In between, H and H' approximate F and F'.");
        Console.WriteLine("");
        Console.WriteLine("     I       X(I)          F(X(I))         H(X(I)) " +
                          "        F'(X(I))        H'(X(I))");
        Console.WriteLine("");
        for (i = 0; i < ns; i++)
        {
            Console.WriteLine("  " + i.ToString().PadLeft(4)
                                   + "  " + xs[i].ToString().PadLeft(14)
                                   + "  " + Math.Sin(xs[i]).ToString().PadLeft(14)
                                   + "  " + ys[i].ToString().PadLeft(14)
                                   + "  " + Math.Cos(xs[i]).ToString().PadLeft(14)
                                   + "  " + ysp[i].ToString().PadLeft(14) + "");
        }

    }
}