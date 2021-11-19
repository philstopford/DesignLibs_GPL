using System;
using Burkardt;
using Burkardt.DifferNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace DifferTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for DIFFER_TEST.
        //
        //  Discussion:
        //
        //    DIFFER_TEST tests the DIFFER library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 November 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("DIFFER_TEST:");
        Console.WriteLine("  Test the DIFFER library.");

        test01();
        test02();
        test03();
        test04();
        test05();

        Console.WriteLine("");
        Console.WriteLine("DIFFER_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests DIFFER_MATRIX.
        //
        //  Discussion:
        //
        //    DIFFER_MATRIX computes a modified Vandermonde matrix A1.
        //
        //    The solution of a system A1 * X1 = B is related to the solution
        //    of the system A2 * X2 = B, where A2 is the standard Vandermonde
        //    matrix, simply by X2(I) = X1(I) * A(I,1).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 November 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] a;
        double[] b;
        int i;
        int info = 0;
        int job;
        int n = 4;
        double[] stencil =  {
                2.5, 3.3, -1.3, 0.5
            }
            ;
        double[] x =  {
                1.0, 2.0, 3.0, 4.0
            }
            ;
        double[] x1;
        double[] x2;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Demonstrate that the DIFFER matrix is 'really'");
        Console.WriteLine("  a Vandermonde matrix.");

        a = Differ.differ_matrix(n, stencil);
        typeMethods.r8mat_print(n, n, a, "  Stencil matrix:");
        b = typeMethods.r8mat_mv_new(n, n, a, x);
        //
        //  Set up and solve the DIFFER system.
        //
        a = Differ.differ_matrix(n, stencil);
        x1 = typeMethods.r8mat_fs_new(n, a, b);

        typeMethods.r8vec_print(n, x1, "  Solution of DIFFER system:");
        //
        //  R8VM_SL solves the related Vandermonde system.
        //  A simple transformation gives us the solution to the DIFFER system.
        //
        job = 0;
        x2 = typeMethods.r8vm_sl_new(n, stencil, b, job, ref info);

        if (info != 0)
        {
            Console.WriteLine("");
            Console.WriteLine("TEST01 - Warning!");
            Console.WriteLine("  VANDERMONDE system is singular.");
            return;
        }

        typeMethods.r8vec_print(n, x2, "  Solution of VANDERMONDE system:");

        for (i = 0; i < n; i++)
        {
            x2[i] /= stencil[i];
        }

        typeMethods.r8vec_print(n, x2, "  Transformed solution of VANDERMONDE system:");
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests DIFFER_INVERSE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 November 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] a;
        double[] b;
        double err;
        int n;
        const int N_MAX = 8;
        int seed;
        int test;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  DIFFER_INVERSE returns the inverse of a DIFFER matrix;");
        Console.WriteLine("");
        Console.WriteLine("   N    Inverse error");

        seed = 123456789;

        for (n = 2; n <= N_MAX; n++)
        {
            Console.WriteLine("");

            for (test = 1; test <= 5; test++)
            {
                x = UniformRNG.r8vec_uniform_01_new(n, ref seed);
                a = Differ.differ_matrix(n, x);
                b = Differ.differ_inverse(n, x);
                err = Helpers.inverse_error(n, a, b);
                Console.WriteLine("  " + n.ToString().PadLeft(2)
                                       + "  " + err.ToString().PadLeft(14) + "");
            }
        }
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests DIFFER_MATRIX.
        //
        //  Discussion:
        //
        //    Reproduce a specific example.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 October 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] a;
        double[] b;
        double[] c;
        double df;
        double dfdx;
        double dx;
        int i;
        int n = 4;
        int order;
        double[] stencil =  {
                -3.0, -2.0, -1.0, 1.0
            }
            ;
        double x0;

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  Reproduce a specific example.");
        //
        //  Compute the coefficients for a specific stencil.
        //
        b = new double[n];
        for (i = 0; i < n; i++)
        {
            b[i] = 0.0;
        }

        order = 1;
        b[order - 1] = 1.0;
        a = Differ.differ_matrix(n, stencil);
        c = typeMethods.r8mat_fs_new(n, a, b);

        typeMethods.r8vec_print(n, c, "  Solution of DIFFER system:");
        //
        //  Use the coefficients C to estimate the first derivative of EXP(X)
        //  at X0, using a spacing of DX = 0.1.
        //
        x0 = 1.3;
        dx = 0.1;
        df = 0.0;
        for (i = 0; i < n; i++)
        {
            df += c[i] * (Math.Exp(x0 + stencil[i] * dx) - Math.Exp(x0));
        }

        dfdx = df / dx;

        Console.WriteLine("");
        Console.WriteLine("  DFDX =         " + dfdx + "");
        Console.WriteLine("  d Math.Exp(x) /dx = " + Math.Exp(x0) + "");
    }

    private static void test04()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 tests DIFFER_FORWARD, DIFFER_BACKWARD, DIFFER_CENTRAL.
        //
        //  Discussion:
        //
        //    Evaluate the coefficients for uniformly spaced finite difference
        //    approximations of derivatives.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 November 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] c;
        double h;
        string label;
        int n;
        int o;
        int p;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TEST04");
        Console.WriteLine("  DIFFER_FORWARD,");
        Console.WriteLine("  DIFFER_BACKWARD, and");
        Console.WriteLine("  DIFFER_CENTRAL produce coefficients for difference");
        Console.WriteLine("  approximations of the O-th derivative,");
        Console.WriteLine("  with error of order H^P, for a uniform spacing of H.");

        h = 1.0;
        Console.WriteLine("");
        Console.WriteLine("  Use a spacing of H = " + h + " for all examples.");
        //
        //  Forward difference approximation to the third derivative with error of O(h).
        //
        o = 3;
        p = 1;
        n = o + p;
        c = new double[n];
        x = new double[n];
        Differ.differ_forward(h, o, p, ref c, ref x);
        label = "  Forward difference coefficients, O = " + o.ToString()
                                                          + ", P = " + p.ToString();
        typeMethods.r8vec2_print(n, x, c, label);
        //
        //  Backward difference approximation to the third derivative with error of O(h).
        //
        o = 3;
        p = 1;
        n = o + p;
        c = new double[n];
        x = new double[n];
        Differ.differ_backward(h, o, p, ref c, ref x);
        label = "  Backward difference coefficients, O = " + o.ToString()
                                                           + ", P = " + p.ToString();
        typeMethods.r8vec2_print(n, x, c, label);
        //
        //  Central difference approximation to the third derivative with error of O(h^2).
        //
        o = 3;
        p = 2;
        n = o + p;
        c = new double[n];
        x = new double[n];
        Differ.differ_central(h, o, p, ref c, ref x);
        label = "  Central difference coefficients, O = " + o.ToString()
                                                          + ", P = " + p.ToString();
        typeMethods.r8vec2_print(n, x, c, label);
        //
        //  Central difference approximation to the third derivative with error of O(h^4).
        //
        o = 3;
        p = 4;
        n = o + p;
        c = new double[n];
        x = new double[n];
        Differ.differ_central(h, o, p, ref c, ref x);
        label = "  Central difference coefficients, O = " + o.ToString()
                                                          + ", P = " + p.ToString();
        typeMethods.r8vec2_print(n, x, c, label);
        //
        //  Forward difference approximation to the fourth derivative with error of O(h).
        //
        o = 4;
        p = 1;
        n = o + p;
        c = new double[n];
        x = new double[n];
        Differ.differ_forward(h, o, p, ref c, ref x);
        label = "  Forward difference coefficients, O = " + o.ToString()
                                                          + ", P = " + p.ToString();
        typeMethods.r8vec2_print(n, x, c, label);
        //
        //  Backward difference approximation to the fourth derivative with error of O(h).
        //
        o = 4;
        p = 1;
        n = o + p;
        c = new double[n];
        x = new double[n];
        Differ.differ_backward(h, o, p, ref c, ref x);
        label = "  Backward difference coefficients, O = " + o.ToString()
                                                           + ", P = " + p.ToString();
        typeMethods.r8vec2_print(n, x, c, label);
        //
        //   Central difference approximation to the fourth derivative with error of O(h^3).
        //
        o = 4;
        p = 3;
        n = o + p;
        c = new double[n];
        x = new double[n];
        Differ.differ_central(h, o, p, ref c, ref x);
        label = "  Central difference coefficients, O = " + o.ToString()
                                                          + ", P = " + p.ToString();
        typeMethods.r8vec2_print(n, x, c, label);
    }

    private static void test05()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 tests DIFFER_STENCIL.
        //
        //  Discussion:
        //
        //    Evaluate the coefficients for uniformly spaced finite difference
        //    approximations of derivatives.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 November 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double[] c;
        double h;
        int i;
        string label;
        int n;
        int o;
        int p;
        double[] x;
        double x0;

        Console.WriteLine("");
        Console.WriteLine("TEST05");
        Console.WriteLine("  DIFFER_STENCIL produces coefficients for difference");
        Console.WriteLine("  approximations of the O-th derivative,");
        Console.WriteLine("  using arbitrarily spaced data, with maximum spacing H");
        Console.WriteLine("  with error of order H^P.");
        //
        //  Let X0 = 1.0.
        //
        x0 = 0.0;
        h = 1.0;
        Console.WriteLine("");
        Console.WriteLine("  Use a spacing of H = " + h + " for all examples.");
        //
        //  Forward difference approximation to the third derivative with error of O(h).
        //
        o = 3;
        p = 1;
        n = o + p;
        c = new double[n];
        x = new double[n];
        for (i = 0; i < n; i++)
        {
            x[i] = i * h;
        }

        Differ.differ_stencil(x0, o, p, x, ref c);
        label = "  Forward difference coefficients, O = " + o.ToString()
                                                          + ", P = " + p.ToString();
        typeMethods.r8vec2_print(n, x, c, label);
        //
        //  Backward difference approximation to the third derivative with error of O(h).
        //
        o = 3;
        p = 1;
        n = o + p;
        c = new double[n];
        x = new double[n];
        for (i = 0; i < n; i++)
        {
            x[i] = (i + 1 - n) * h;
        }

        Differ.differ_stencil(x0, o, p, x, ref c);
        label = "  Backward difference coefficients, O = " + o.ToString()
                                                           + ", P = " + p.ToString();
        typeMethods.r8vec2_print(n, x, c, label);
        //
        //  Central difference approximation to the third derivative with error of O(h^2).
        //
        o = 3;
        p = 2;
        n = o + p;
        c = new double[n];
        x = new double[n];
        for (i = 0; i < n; i++)
        {
            x[i] = (-n + 1 + 2 * i) * h / 2.0;
        }

        Differ.differ_stencil(x0, o, p, x, ref c);
        label = "  Central difference coefficients, O = " + o.ToString()
                                                          + ", P = " + p.ToString();
        typeMethods.r8vec2_print(n, x, c, label);
        //
        //  Central difference approximation to the third derivative with error of O(h^4).
        //
        o = 3;
        p = 4;
        n = o + p;
        c = new double[n];
        x = new double[n];
        for (i = 0; i < n; i++)
        {
            x[i] = (-n + 1 + 2 * i) * h / 2.0;
        }

        Differ.differ_stencil(x0, o, p, x, ref c);
        label = "  Central difference coefficients, O = " + o.ToString()
                                                          + ", P = " + p.ToString();
        typeMethods.r8vec2_print(n, x, c, label);
        //
        //  Forward difference approximation to the fourth derivative with error of O(h).
        //
        o = 4;
        p = 1;
        n = o + p;
        c = new double[n];
        x = new double[n];
        for (i = 0; i < n; i++)
        {
            x[i] = i * h;
        }

        Differ.differ_stencil(x0, o, p, x, ref c);
        label = "  Forward difference coefficients, O = " + o.ToString()
                                                          + ", P = " + p.ToString();
        typeMethods.r8vec2_print(n, x, c, label);
        //
        //  Backward difference approximation to the fourth derivative with error of O(h).
        //
        o = 4;
        p = 1;
        n = o + p;
        c = new double[n];
        x = new double[n];
        for (i = 0; i < n; i++)
        {
            x[i] = (i + 1 - n) * h;
        }

        Differ.differ_stencil(x0, o, p, x, ref c);
        label = "  Backward difference coefficients, O = " + o.ToString()
                                                           + ", P = " + p.ToString();
        typeMethods.r8vec2_print(n, x, c, label);
        //
        //   Central difference approximation to the fourth derivative with error of O(h^3).
        //
        o = 4;
        p = 3;
        n = o + p;
        c = new double[n];
        x = new double[n];
        for (i = 0; i < n; i++)
        {
            x[i] = (-n + 1 + 2 * i) * h / 2.0;
        }

        Differ.differ_stencil(x0, o, p, x, ref c);
        label = "  Central difference coefficients, O = " + o.ToString()
                                                          + ", P = " + p.ToString();
        typeMethods.r8vec2_print(n, x, c, label);
    }
}