﻿using System;
using System.Globalization;
using Burkardt.Grid;
using Burkardt.Types;

namespace MultiGridPoisson1DTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for MULTIGRID_POISSON_1D_TEST.
        //
        //  Discussion:
        //
        //    MULTIGRID_POISSON_1D_TEST tests the MULTIGRID_POISSON_1D library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 December 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("MULTIGRID_POISSON_1D:");
        Console.WriteLine("  Test the MULTIGRID_POISSON_1D library.");

        test01_mono();
        test01_multi();
        test02_mono();
        test02_multi();

        Console.WriteLine("");
        Console.WriteLine("MULTIGRID_POISSON_1D:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01_mono()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01_MONO tests MONOGRID_POISSON_1D on test case 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 December 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int it_num = 0;
        int k;

        Console.WriteLine("");
        Console.WriteLine("TEST01_MONO");
        Console.WriteLine("  MONOGRID_POISSON_1D solves a 1D Poisson BVP");
        Console.WriteLine("  using the Gauss-Seidel method.");

        const double a = 0.0;
        const double b = 1.0;
        const double ua = 0.0;
        const double ub = 0.0;

        Console.WriteLine("");
        Console.WriteLine("  -u''(x) = 1, for 0 < x < 1");
        Console.WriteLine("  u(0) = u(1) = 0.");
        Console.WriteLine("  Solution is u(x) = ( -x^2 + x ) / 2");

        for (k = 5; k <= 5; k++)
        {
            int n = (int) Math.Pow(2, k);

            double[] u = new double[n + 1];
            double[] x = typeMethods.r8vec_linspace_new(n + 1, a, b);

            Console.WriteLine("");
            Console.WriteLine("  Mesh index K = " + k + "");
            Console.WriteLine("  Number of intervals N=2^K = " + n + "");
            Console.WriteLine("  Number of nodes = 2^K+1 =   " + (n + 1) + "");

            Poisson.monogrid_poisson_1d(n, a, b, ua, ub, force1, exact1, ref it_num, ref u);

            Console.WriteLine("");
            Console.WriteLine("     I        X(I)      U(I)         U Exact(X(I))");
            Console.WriteLine("");
            int i;
            for (i = 0; i < n + 1; i++)
            {
                Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                       + "  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + exact1(x[i]).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }

            Console.WriteLine("");

            double difmax = 0.0;
            for (i = 0; i < n + 1; i++)
            {
                difmax = Math.Max(difmax, Math.Abs(u[i] - exact1(x[i])));
            }

            Console.WriteLine("  Maximum error = " + difmax + "");
            Console.WriteLine("  Number of iterations = " + it_num + "");
        }
    }

    private static void test01_multi()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01_MULTI tests MULTIGRID_POISSON_1D on test case 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int it_num = 0;
        int k;

        Console.WriteLine("");
        Console.WriteLine("TEST01_MULTI");
        Console.WriteLine("  MULTIGRID_POISSON_1D solves a 1D Poisson BVP");
        Console.WriteLine("  using the multigrid method.");

        const double a = 0.0;
        const double b = 1.0;
        const double ua = 0.0;
        const double ub = 0.0;

        Console.WriteLine("");
        Console.WriteLine("  -u''(x) = 1, for 0 < x < 1");
        Console.WriteLine("  u(0) = u(1) = 0.");
        Console.WriteLine("  Solution is u(x) = ( -x^2 + x ) / 2");

        for (k = 5; k <= 5; k++)
        {
            int n = (int) Math.Pow(2, k);

            double[] u = new double[n + 1];
            double[] x = typeMethods.r8vec_linspace_new(n + 1, a, b);

            Console.WriteLine("");
            Console.WriteLine("  Mesh index K = " + k + "");
            Console.WriteLine("  Number of intervals N=2^K = " + n + "");
            Console.WriteLine("  Number of nodes = 2^K+1 =   " + (n + 1) + "");

            Poisson.multigrid_poisson_1d(n, a, b, ua, ub, force1, exact1, ref it_num, ref u);

            Console.WriteLine("");
            Console.WriteLine("     I        X(I)      U(I)         U Exact(X(I))");
            Console.WriteLine("");
            int i;
            for (i = 0; i < n + 1; i++)
            {
                Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                       + "  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + exact1(x[i]).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }

            Console.WriteLine("");

            double difmax = 0.0;
            for (i = 0; i < n + 1; i++)
            {
                difmax = Math.Max(difmax, Math.Abs(u[i] - exact1(x[i])));
            }

            Console.WriteLine("  Maximum error = " + difmax + "");
            Console.WriteLine("  Number of iterations = " + it_num + "");
        }
    }

    private static double exact1(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT1 evaluates the exact solution.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 December 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    William Hager,
        //    Applied Numerical Linear Algebra,
        //    Prentice-Hall, 1988,
        //    ISBN13: 978-0130412942,
        //    LC: QA184.H33.
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double EXACT1, the value of the exact solution at X.
        //
    {
        double value = 0;

        value = 0.5 * (-x * x + x);

        return value;
    }

    private static double force1(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FORCE1 evaluates the forcing function.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 December 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    William Hager,
        //    Applied Numerical Linear Algebra,
        //    Prentice-Hall, 1988,
        //    ISBN13: 978-0130412942,
        //    LC: QA184.H33.
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double FORCE1, the value of the forcing function at X.
        //
    {
        double value = 0;

        value = 1.0;

        return value;
    }

    private static void test02_mono()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02_MONO tests MONOGRID_POISSON_1D on test case 2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int it_num = 0;
        int k;

        Console.WriteLine("");
        Console.WriteLine("TEST02_MONO");
        Console.WriteLine("  MONOGRID_POISSON_1D solves a 1D Poisson BVP");
        Console.WriteLine("  using the Gauss-Seidel method.");

        const double a = 0.0;
        const double b = 1.0;
        const double ua = 0.0;
        const double ub = 0.0;

        Console.WriteLine("");
        Console.WriteLine("  -u''(x) = - x * (x+3) * exp(x), for 0 < x < 1");
        Console.WriteLine("  u(0) = u(1) = 0.");
        Console.WriteLine("  Solution is u(x) = x * (x-1) * exp(x)");

        for (k = 5; k <= 5; k++)
        {
            int n = (int) Math.Pow(2, k);

            double[] u = new double[n + 1];
            double[] x = typeMethods.r8vec_linspace_new(n + 1, a, b);

            Console.WriteLine("");
            Console.WriteLine("  Mesh index K = " + k + "");
            Console.WriteLine("  Number of intervals N=2^K = " + n + "");
            Console.WriteLine("  Number of nodes = 2^K+1 =   " + (n + 1) + "");

            Poisson.monogrid_poisson_1d(n, a, b, ua, ub, force2, exact2, ref it_num, ref u);

            Console.WriteLine("");
            Console.WriteLine("     I        X(I)      U(I)         U Exact(X(I))");
            Console.WriteLine("");
            int i;
            for (i = 0; i < n + 1; i++)
            {
                Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                       + "  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + exact2(x[i]).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }

            Console.WriteLine("");

            double difmax = 0.0;
            for (i = 0; i < n + 1; i++)
            {
                difmax = Math.Max(difmax, Math.Abs(u[i] - exact2(x[i])));
            }

            Console.WriteLine("  Maximum error = " + difmax + "");
            Console.WriteLine("  Number of iterations = " + it_num + "");
        }
    }

    private static void test02_multi()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02_MULTI tests MULTIGRID_POISSON_1D on test case 2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int it_num = 0;
        int k;

        Console.WriteLine("");
        Console.WriteLine("TEST02_MULTI");
        Console.WriteLine("  MULTIGRID_POISSON_1D solves a 1D Poisson BVP");
        Console.WriteLine("  using the multigrid method.");

        const double a = 0.0;
        const double b = 1.0;
        const double ua = 0.0;
        const double ub = 0.0;

        Console.WriteLine("");
        Console.WriteLine("  -u''(x) = - x * (x+3) * exp(x), for 0 < x < 1");
        Console.WriteLine("  u(0) = u(1) = 0.");
        Console.WriteLine("  Solution is u(x) = x * (x-1) * exp(x)");

        for (k = 5; k <= 5; k++)
        {
            int n = (int) Math.Pow(2, k);

            double[] u = new double[n + 1];
            double[] x = typeMethods.r8vec_linspace_new(n + 1, a, b);

            Console.WriteLine("");
            Console.WriteLine("  Mesh index K = " + k + "");
            Console.WriteLine("  Number of intervals N=2^K = " + n + "");
            Console.WriteLine("  Number of nodes = 2^K+1 =   " + (n + 1) + "");

            Poisson.multigrid_poisson_1d(n, a, b, ua, ub, force2, exact2, ref it_num, ref u);

            Console.WriteLine("");
            Console.WriteLine("     I        X(I)      U(I)         U Exact(X(I))");
            Console.WriteLine("");
            int i;
            for (i = 0; i < n + 1; i++)
            {
                Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                       + "  " + x[i].ToString(CultureInfo.InvariantCulture).PadLeft(10)
                                       + "  " + u[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + exact2(x[i]).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }

            Console.WriteLine("");

            double difmax = 0.0;
            for (i = 0; i < n + 1; i++)
            {
                difmax = Math.Max(difmax, Math.Abs(u[i] - exact2(x[i])));
            }

            Console.WriteLine("  Maximum error = " + difmax + "");
            Console.WriteLine("  Number of iterations = " + it_num + "");
        }
    }

    private static double exact2(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT2 evaluates the exact solution.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 December 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    William Hager,
        //    Applied Numerical Linear Algebra,
        //    Prentice-Hall, 1988,
        //    ISBN13: 978-0130412942,
        //    LC: QA184.H33.
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double EXACT2, the value of the exact solution at X.
        //
    {
        double value = x * (x - 1.0) * Math.Exp(x);

        return value;
    }

    private static double force2(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FORCE2 evaluates the forcing function.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 December 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    William Hager,
        //    Applied Numerical Linear Algebra,
        //    Prentice-Hall, 1988,
        //    ISBN13: 978-0130412942,
        //    LC: QA184.H33.
        //
        //  Parameters:
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double FORCE2, the value of the forcing function at X.
        //
    {
        double value = -x * (x + 3.0) * Math.Exp(x);

        return value;
    }
}