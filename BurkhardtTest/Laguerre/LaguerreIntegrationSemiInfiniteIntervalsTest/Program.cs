using System;
using Burkardt.Laguerre;

namespace LaguerreIntegrationSemiInfiniteIntervalsTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for LAGUERRE_TEST_INT_TEST.
        //
        //  Discussion:
        //
        //    LAGUERRE_TEST_INT_TEST tests the LAGUERRE_TEST_INT library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 July 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("LAGUERRE_TEST_INT_TEST");
        Console.WriteLine("  Test the LAGUERRE_TEST_INT library.");

        test01();
        test02();
        test03();
        test04();
        test05();

        Console.WriteLine("");
        Console.WriteLine("LAGUERRE_TEST_INT_TEST");
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
        //    27 December 2011
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

        problem_num = Integrands.p00_problem_num();

        Console.WriteLine("");
        Console.WriteLine("  P00_PROBLEM_NUM: number of problems is " + problem_num + "");
        Console.WriteLine("");
        Console.WriteLine("   Problem       Title");
        Console.WriteLine("");

        for (problem = 1; problem <= problem_num; problem++)
        {
            title = Integrands.p00_title(problem);

            Console.WriteLine("  " + problem.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  \"" + title + "\".");
        }
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests P00_ALPHA and P00_EXACT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 July 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double alpha;
        double exact;
        int problem;
        int problem_num;
        Integrands.p00Data data = new();

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  P00_ALPHA returns the lower limit of integration.");
        Console.WriteLine("  P00_EXACT returns the \"exact\" integral.");

        problem_num = Integrands.p00_problem_num();

        Console.WriteLine("");
        Console.WriteLine("   Problem       ALPHA           EXACT");
        Console.WriteLine("");

        for (problem = 1; problem <= problem_num; problem++)
        {
            alpha = Integrands.p00_alpha(problem);

            exact = Integrands.p00_exact(ref data, problem);

            Console.WriteLine("  " + problem.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + alpha.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + exact.ToString("0.################").PadLeft(24) + "");
        }
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests P00_GAUSS_LAGUERRE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 December 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double error;
        double estimate;
        double exact;
        int order;
        int order_log;
        int problem;
        int problem_num;
        Integrands.p00Data data = new();

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  P00_GAUSS_LAGUERRE applies a Gauss-Laguerre rule");
        Console.WriteLine("  to estimate an integral on [ALPHA,+oo).");

        problem_num = Integrands.p00_problem_num();

        Console.WriteLine("");
        Console.WriteLine("                              Exact");
        Console.WriteLine("   Problem     Order          Estimate        Error");

        for (problem = 1; problem <= problem_num; problem++)
        {
            exact = Integrands.p00_exact(ref data, problem);

            order = 1;

            Console.WriteLine("");
            Console.WriteLine("  " + problem.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + "        "
                                   + "  " + exact.ToString("0.######").PadLeft(14) + "");

            for (order_log = 0; order_log <= 6; order_log++)
            {
                estimate = Integrands.p00_gauss_laguerre(ref data, problem, order);

                error = Math.Abs(exact - estimate);

                Console.WriteLine("  " + "        "
                                       + "  " + order.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + estimate.ToString("0.######").PadLeft(14)
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
        //    TEST04 tests P00_EXP_TRANSFORM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 December 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double error;
        double estimate;
        double exact;
        int order;
        int order_log;
        int problem;
        int problem_num;
        Integrands.p00Data data = new();

        Console.WriteLine("");
        Console.WriteLine("TEST04");
        Console.WriteLine("  P00_EXP_TRANSFORM applies an exponential transform");
        Console.WriteLine("  to estimate an integral on [ALPHA,+oo)");
        Console.WriteLine("  as a transformed integral on (0,exp(-ALPHA)],");
        Console.WriteLine("  and applying a Gauss-Legendre rule.");

        problem_num = Integrands.p00_problem_num();

        Console.WriteLine("");
        Console.WriteLine("                              Exact");
        Console.WriteLine("   Problem     Order          Estimate        Error");

        for (problem = 1; problem <= problem_num; problem++)
        {
            exact = Integrands.p00_exact(ref data, problem);

            order = 1;

            Console.WriteLine("");
            Console.WriteLine("  " + problem.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + "        "
                                   + "  " + exact.ToString("0.######").PadLeft(14) + "");

            for (order_log = 0; order_log <= 6; order_log++)
            {
                estimate = Integrands.p00_exp_transform(ref data, problem, order);

                error = Math.Abs(exact - estimate);

                Console.WriteLine("  " + "        "
                                       + "  " + order.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + estimate.ToString("0.######").PadLeft(14)
                                       + "  " + error.ToString("0.######").PadLeft(14) + "");

                order *= 2;
            }
        }
    }

    private static void test05()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 tests P00_RAT_TRANSFORM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 December 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double error;
        double estimate;
        double exact;
        int order;
        int order_log;
        int problem;
        int problem_num;
        Integrands.p00Data data = new();

        Console.WriteLine("");
        Console.WriteLine("TEST05");
        Console.WriteLine("  P00_RAT_TRANSFORM applies a rational transform");
        Console.WriteLine("  to estimate an integral on [ALPHA,+oo)");
        Console.WriteLine("  as a transformed integral on (0,1/(1+ALPHA)],");
        Console.WriteLine("  and applying a Gauss-Legendre rule.");

        problem_num = Integrands.p00_problem_num();

        Console.WriteLine("");
        Console.WriteLine("                              Exact");
        Console.WriteLine("   Problem     Order          Estimate        Error");

        for (problem = 1; problem <= problem_num; problem++)
        {
            exact = Integrands.p00_exact(ref data, problem);

            order = 1;

            Console.WriteLine("");
            Console.WriteLine("  " + problem.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + "        "
                                   + "  " + exact.ToString("0.######").PadLeft(14) + "");

            for (order_log = 0; order_log <= 6; order_log++)
            {
                estimate = Integrands.p00_rat_transform(ref data, problem, order);

                error = Math.Abs(exact - estimate);

                Console.WriteLine("  " + "        "
                                       + "  " + order.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                       + "  " + estimate.ToString("0.######").PadLeft(14)
                                       + "  " + error.ToString("0.######").PadLeft(14) + "");

                order *= 2;
            }
        }
    }
}