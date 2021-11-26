using System;
using Burkardt.Tessellation;
using Burkardt.Types;
using Burkardt.Uniform;

namespace LloydCVT_LineTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_CVT_LLOYD_TEST tests the LINE_CVT_LLOYD library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("LINE_CVT_LLOYD_TEST");
        Console.WriteLine("  Test the LINE_CVT_LLOYD library.");

        test01();
        test02();
        //
        //  Repeat, using sorted initial points.
        //
        test03();
        test04();

        Console.WriteLine("");
        Console.WriteLine("LINE_CVT_LLOYD_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_CVT_LLOYD_TEST01 tests the unconstrained computation.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const string header = "test01";
        const int n = 25;

        Console.WriteLine("");
        Console.WriteLine("LINE_CVT_LLOYD_TEST01:");
        Console.WriteLine("  Test the unconstrained computation.");

        const double a = 0.0;
        const double b = 1.0;
        const int it_num = 200;
        int seed = 123456789;
        double[] x = UniformRNG.r8vec_uniform_ab_new(n, a, b, ref seed);

        Console.WriteLine("");
        Console.WriteLine("  Use " + n + " points in the interval [" + a + "," + b + "]");
        Console.WriteLine("  Number of iterations to take is " + it_num + "");
        Console.WriteLine("  Call this calculation '" + header + "'");
        const double h = (b - a) / (n - 1);
        Console.WriteLine("  Expect a uniform spacing of " + h + "");

        typeMethods.r8vec_print(n, x, "  Initial generators:");

        LloydCVT_Line.line_cvt_lloyd(n, a, b, it_num, header, x);

        typeMethods.r8vec_print(n, x, "  Final generators:");
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_CVT_LLOYD_TEST02 tests the constrained computation.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const string header = "test02";
        const int n = 25;

        Console.WriteLine("");
        Console.WriteLine("LINE_CVT_LLOYD_TEST02:");
        Console.WriteLine("  Test the constrained computation.");

        const double a = 0.0;
        const double b = 1.0;
        const int it_num = 200;
        int seed = 123456789;
        double[] x = UniformRNG.r8vec_uniform_ab_new(n, a, b, ref seed);

        Console.WriteLine("");
        Console.WriteLine("  Use " + n + " points in the interval [" + a + "," + b + "]");
        Console.WriteLine("  Number of iterations to take is " + it_num + "");
        Console.WriteLine("  Call this calculation '" + header + "'");
        const double h = (b - a) / n;
        Console.WriteLine("  Expect a uniform spacing of " + h + "");

        typeMethods.r8vec_print(n, x, "  Initial generators:");

        LloydCVT_Line.line_ccvt_lloyd(n, a, b, it_num, header, ref x);

        typeMethods.r8vec_print(n, x, "  Final generators:");
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_CVT_LLOYD_TEST03 tests the unconstrained computation.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const string header = "test03";
        const int n = 25;

        Console.WriteLine("");
        Console.WriteLine("LINE_CVT_LLOYD_TEST03:");
        Console.WriteLine("  Test the unconstrained computation.");
        Console.WriteLine("  SORT the random initial values before use.");

        const double a = 0.0;
        const double b = 1.0;
        const int it_num = 200;
        int seed = 123456789;
        double[] x = UniformRNG.r8vec_uniform_ab_new(n, a, b, ref seed);
        typeMethods.r8vec_sort_insert_a(n, ref x);

        Console.WriteLine("");
        Console.WriteLine("  Use " + n + " points in the interval [" + a + "," + b + "]");
        Console.WriteLine("  Number of iterations to take is " + it_num + "");
        Console.WriteLine("  Call this calculation '" + header + "'");
        double h = (b - a) / (n - 1);
        Console.WriteLine("  Expect a uniform spacing of " + h + "");

        typeMethods.r8vec_print(n, x, "  Initial generators:");

        LloydCVT_Line.line_cvt_lloyd(n, a, b, it_num, header, x);

        typeMethods.r8vec_print(n, x, "  Final generators:");
    }

    private static void test04()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LINE_CVT_LLOYD_TEST04 tests the constrained computation.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 July 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        string header = "test04";
        const int n = 25;

        Console.WriteLine("");
        Console.WriteLine("LINE_CVT_LLOYD_TEST04:");
        Console.WriteLine("  Test the constrained computation.");
        Console.WriteLine("  SORT the initial points before use.");

        const double a = 0.0;
        const double b = 1.0;
        const int it_num = 200;
        int seed = 123456789;
        double[] x = UniformRNG.r8vec_uniform_ab_new(n, a, b, ref seed);
        typeMethods.r8vec_sort_insert_a(n, ref x);

        Console.WriteLine("");
        Console.WriteLine("  Use " + n + " points in the interval [" + a + "," + b + "]");
        Console.WriteLine("  Number of iterations to take is " + it_num + "");
        Console.WriteLine("  Call this calculation '" + header + "'");
        const double h = (b - a) / n;
        Console.WriteLine("  Expect a uniform spacing of " + h + "");

        typeMethods.r8vec_print(n, x, "  Initial generators:");

        LloydCVT_Line.line_ccvt_lloyd(n, a, b, it_num, header, ref x);

        typeMethods.r8vec_print(n, x, "  Final generators:");
    }
}