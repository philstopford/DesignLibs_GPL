using System;

namespace HermiteIntegrandsTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for hermite_integrands_TEST.
        //
        //  Discussion:
        //
        //    hermite_integrands_TEST tests the hermite_integrands library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 May 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("hermite_integrands_TEST");
        Console.WriteLine("  Test the hermite_integrands library.");

        test01();
        test02();
        test03();
        test04();
        test05();
        test06();
        Console.WriteLine("");
        Console.WriteLine("hermite_integrands_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests P00_PROBLEM_NUM and P00_TITLE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    26 May 2009 
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int problem;
        int problem_num;
        string title;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  P00_PROBLEM_NUM returns the number of problems.");
        Console.WriteLine("  P00_TITLE returns the title of a problem.");

        problem_num = Problem00.p00_problem_num();

        Console.WriteLine("");
        Console.WriteLine("  P00_PROBLEM_NUM: number of problems is " + problem_num + "");
        Console.WriteLine("");
        Console.WriteLine("   Problem       Title");
        Console.WriteLine("");

        for (problem = 1; problem <= problem_num; problem++)
        {
            title = Problem00.p00_title(problem);

            Console.WriteLine("  " + problem.ToString().PadLeft(8)
                                   + "  \"" + title + "\".");
        }
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests P00_EXACT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    26 May 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double exact;
        int m;
        int problem;
        int problem_num;
        Problem00.p00Data data = new();

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  P00_EXACT returns the \"exact\" integral.");

        problem_num = Problem00.p00_problem_num();

        m = 4;
        Problem06.p06_param(ref data.p6data, 'S', 'M', ref m);

        Console.WriteLine("");
        Console.WriteLine("   Problem       EXACT");
        Console.WriteLine("");

        for (problem = 1; problem <= problem_num; problem++)
        {
            exact = Problem00.p00_exact(ref data, problem);

            Console.WriteLine("  " + problem.ToString().PadLeft(8)
                                   + "  " + exact.ToString("0.################").PadLeft(24) + "");
        }
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests P00_GAUSS_HERMITE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    26 May 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double error;
        double estimate;
        double exact;
        int m;
        int order;
        int order_log;
        int problem;
        int problem_num;

        Problem00.p00Data data = new();

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  P00_GAUSS_HERMITE applies a Gauss-Hermite rule");
        Console.WriteLine("  to estimate an integral on (-oo,+oo).");

        problem_num = Problem00.p00_problem_num();

        m = 4;
        Problem06.p06_param(ref data.p6data, 'S', 'M', ref  m);

        Console.WriteLine("");
        Console.WriteLine("   Problem     Order          Estimate        Exact          Error");

        for (problem = 1; problem <= problem_num; problem++)
        {
            exact = Problem00.p00_exact(ref data, problem);

            order = 1;

            Console.WriteLine("");

            for (order_log = 0; order_log <= 6; order_log++)
            {
                estimate = Problem00.p00_gauss_hermite(ref data, problem, order);

                error = Math.Abs(exact - estimate);

                Console.WriteLine("  " + problem.ToString().PadLeft(8)
                                       + "  " + order.ToString().PadLeft(8)
                                       + "  " + estimate.ToString("0.######").PadLeft(14)
                                       + "  " + exact.ToString("0.######").PadLeft(14)
                                       + "  " + error.ToString("0.######").PadLeft(14) + "");

                order *= 2;
            }
        }
    }

    private static void test04()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 tests P00_TURING.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    26 May 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double error;
        double estimate;
        double exact;
        double h;
        int m;
        int n = 0;
        int order_log;
        int problem;
        int problem_num;
        int test;
        double tol = 0;

        Problem00.p00Data data = new();

        Console.WriteLine("");
        Console.WriteLine("TEST04");
        Console.WriteLine("  P00_TURING applies a Turing procedure");
        Console.WriteLine("  to estimate an integral on (-oo,+oo).");

        problem_num = Problem00.p00_problem_num();

        m = 4;
        Problem06.p06_param(ref data.p6data, 'S', 'M', ref m);

        for (test = 1; test <= 2; test++)
        {
            tol = test switch
            {
                1 => 1.0E-4,
                2 => 1.0E-07,
                _ => tol
            };

            Console.WriteLine("");
            Console.WriteLine("  Using a tolerance of TOL = " + tol + "");
            Console.WriteLine("");
            Console.WriteLine(    "   Problem     Order          Estimate        Exact          Error");

            for (problem = 1; problem <= problem_num; problem++)
            {
                exact = Problem00.p00_exact(ref data, problem);

                h = 1.0;

                Console.WriteLine("");

                for (order_log = 0; order_log <= 6; order_log++)
                {
                    estimate = Problem00.p00_turing(ref data, problem, h, tol, ref n);

                    error = Math.Abs(exact - estimate);

                    Console.WriteLine("  " + problem.ToString().PadLeft(8)
                                           + "  " + h.ToString().PadLeft(10)
                                           + "  " + n.ToString().PadLeft(8)
                                           + "  " + estimate.ToString("0.######").PadLeft(14)
                                           + "  " + exact.ToString("0.######").PadLeft(14)
                                           + "  " + error.ToString("0.######").PadLeft(14) + "");

                    h /= 2.0;
                }
            }
        }
    }

    private static void test05()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 tests P00_GAUSS_HERMITE against the polynomials.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 May 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double error;
        double estimate;
        double exact;
        int m;
        int order;
        int problem;

        Problem00.p00Data data = new();

        Console.WriteLine("");
        Console.WriteLine("TEST05");
        Console.WriteLine("  P00_GAUSS_HERMITE applies a Gauss-Hermite rule to");
        Console.WriteLine("  estimate the integral x^m exp(-x*x) over (-oo,+oo).");

        problem = 6;

        Console.WriteLine("");
        Console.WriteLine("         M     Order      Estimate        Exact           Error");

        for (m = 0; m <= 6; m++)
        {
            Problem06.p06_param(ref data.p6data, 'S', 'M', ref  m);

            exact = Problem00.p00_exact(ref data, problem);

            Console.WriteLine("");

            for (order = 1; order <= 3 + m / 2; order++)
            {
                estimate = Problem00.p00_gauss_hermite(ref data, problem, order);

                error = Math.Abs(exact - estimate);

                Console.WriteLine("  " + m.ToString().PadLeft(8)
                                       + "  " + order.ToString().PadLeft(8)
                                       + "  " + estimate.ToString().PadLeft(14)
                                       + "  " + exact.ToString().PadLeft(14)
                                       + "  " + error.ToString().PadLeft(14) + "");
            }
        }
    }

    private static void test06()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST06 tests P00_MONTE_CARLO.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 May 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double error;
        double estimate;
        double exact;
        int m;
        int order;
        int order_log;
        int problem;
        int problem_num;

        Problem00.p00Data data = new();

        Console.WriteLine("");
        Console.WriteLine("TEST06");
        Console.WriteLine("  P00_MONTE_CARLO uses a weighted form of the Monte Carlo method");
        Console.WriteLine("  to estimate a Hermite integral on (-oo,+oo).");

        problem_num = Problem00.p00_problem_num();

        m = 4;
        Problem06.p06_param(ref data.p6data, 'S', 'M', ref  m);

        Console.WriteLine("");
        Console.WriteLine("   Problem     Order          Estimate        Exact          Error");

        for (problem = 1; problem <= problem_num; problem++)
        {
            exact = Problem00.p00_exact(ref data, problem);

            order = 128;

            Console.WriteLine("");

            for (order_log = 0; order_log <= 6; order_log++)
            {
                estimate = Problem00.p00_monte_carlo(ref data, problem, order);

                error = Math.Abs(exact - estimate);

                Console.WriteLine("  " + problem.ToString().PadLeft(8)
                                       + "  " + order.ToString().PadLeft(8)
                                       + "  " + estimate.ToString("0.######").PadLeft(14)
                                       + "  " + exact.ToString("0.######").PadLeft(14)
                                       + "  " + error.ToString("0.######").PadLeft(14) + "");

                order *= 4;
            }
        }
    }
}