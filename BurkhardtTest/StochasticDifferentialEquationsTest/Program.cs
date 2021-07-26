using System;
using Burkardt.StochasticDifferentialEquations;

namespace StochasticDifferentialEquationsTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************
            //
            //  Purpose:
            //
            //    MAIN is the main program for SDE_TEST.
            //
            //  Discussion:
            //
            //    SDE_TEST tests the SDE library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    28 September 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("SDE_TEST");
            Console.WriteLine("  Test the SDE library.");

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

            Console.WriteLine("");
            Console.WriteLine("SDE_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void test01()

            //****************************************************************************
            //
            //  Purpose:
            //
            //    TEST01 tests BPATH.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 September 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int n = 500;
            int seed;
            double[] w;

            Console.WriteLine("");
            Console.WriteLine("TEST01:");
            Console.WriteLine("  BPATH generates a sample Brownian motion path.");

            seed = 123456789;

            w = BrownianPath.bpath(ref seed, n);

            BrownianPath.bpath_gnuplot(n, w);
        }

        static void test02()

            //****************************************************************************
            //
            //  Purpose:
            //
            //    TEST02 tests BPATH_AVERAGE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 September 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double error = 0;
            int m = 1000;
            int n = 500;
            int seed;
            double[] u;
            double[] umean;

            Console.WriteLine("");
            Console.WriteLine("TEST02:");
            Console.WriteLine("  BPATH_AVERAGE generates many Brownian paths");
            Console.WriteLine("  and averages them.");

            seed = 123456789;
            u = new double[m * (n + 1)];
            umean = new double[n + 1];

            BrownianPath.bpath_average(ref seed, m, n, ref u, ref umean, ref error);

            BrownianPath.bpath_average_gnuplot(m, n, u, umean);
        }

        static void test03()

            //****************************************************************************
            //
            //  Purpose:
            //
            //    TEST03 tests CHAIN.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    28 September 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double diff = 0;
            int n = 200;
            int seed;
            double[] vem;
            double[] xem;

            Console.WriteLine("");
            Console.WriteLine("TEST03:");
            Console.WriteLine("  CHAIN solves a stochastic differential equation for");
            Console.WriteLine("  a function of a stochastic variable X.");
            Console.WriteLine("  We can solve for X(t), and then evaluate V(X(t)).");
            Console.WriteLine("  Or, we can apply the stochastic chain rule to derive an");
            Console.WriteLine("  an SDE for V, and solve that.");

            seed = 123456789;
            xem = new double[n + 1];
            vem = new double[n + 1];
            ChainRule.chain(ref seed, n, ref xem, ref vem, ref diff);

            Console.WriteLine("");
            Console.WriteLine("  Maximum | Sqrt(X) - V | = " + diff + "");

            ChainRule.chain_gnuplot(n, xem, vem);
        }

        static void test04()

            //****************************************************************************
            //
            //  Purpose:
            //
            //    TEST04 tests EM.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    28 September 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double diff = 0;
            int n = 256;
            int seed;
            double[] t;
            double[] t2;
            double[] xem;
            double[] xtrue;

            Console.WriteLine("");
            Console.WriteLine("TEST04:");
            Console.WriteLine("  EM solves a stochastic differential equation");
            Console.WriteLine("  using the Euler-Maruyama method.");

            seed = 123456789;

            t = new double[n + 1];
            t2 = new double[1 + n / 4];
            xem = new double[1 + n / 4];
            xtrue = new double[n + 1];

            EulerMaruyama.em(ref seed, n, ref t, ref xtrue, ref t2, ref xem, ref diff);

            Console.WriteLine("");
            Console.WriteLine("  | Exact X(T) - EM X(T) | = " + diff + "");

            EulerMaruyama.em_gnuplot(n, t, xtrue, t2, xem);
        }

        static void test05()

            //****************************************************************************
            //
            //  Purpose:
            //
            //    TEST05 tests EMSTRONG.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    28 September 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] dtvals;
            int m = 100;
            int n = 512;
            int p_max = 6;
            int seed;
            double[] xerr;

            dtvals = new double[p_max];
            xerr = new double[p_max];

            Console.WriteLine("");
            Console.WriteLine("TEST05:");
            Console.WriteLine("  EMSTRONG investigates the strong convergence");
            Console.WriteLine("  of the Euler-Maruyama method.");

            seed = 123456789;

            EulerMaruyamaStrong.emstrong(ref seed, m, n, p_max, ref dtvals, ref xerr);

            EulerMaruyamaStrong.emstrong_gnuplot(p_max, dtvals, xerr);
        }

        static void test06()

            //****************************************************************************
            //
            //  Purpose:
            //
            //    TEST06 tests EMWEAK.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    28 September 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] dtvals;
            int m = 50000;
            int method;
            int p_max = 5;
            int seed;
            double[] xerr;

            dtvals = new double[p_max];
            xerr = new double[p_max];

            Console.WriteLine("");
            Console.WriteLine("TEST06:");
            Console.WriteLine("  EMWEAK investigates the weak convergence");
            Console.WriteLine("  of the Euler-Maruyama method.");

            seed = 123456789;
            method = 0;

            EulerMaruyamaWeak.emweak(ref seed, method, m, p_max, ref dtvals, ref xerr);

            EulerMaruyamaWeak.emweak_gnuplot(p_max, dtvals, xerr, method);

            seed = 123456789;
            method = 1;

            EulerMaruyamaWeak.emweak(ref seed, method, m, p_max, ref dtvals, ref xerr);

            EulerMaruyamaWeak.emweak_gnuplot(p_max, dtvals, xerr, method);
        }

        static void test07()

            //****************************************************************************
            //
            //  Purpose:
            //
            //    TEST07 tests MILSTRONG.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    28 September 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] dtvals;
            int p_max = 4;
            int seed;
            double[] xerr;

            dtvals = new double[p_max];
            xerr = new double[p_max];

            Console.WriteLine("");
            Console.WriteLine("TEST07:");
            Console.WriteLine("  MILSTRONG investigates the strong convergence");
            Console.WriteLine("  of the Milstein method.");

            seed = 123456789;

            MilsteinStrong.milstrong(ref seed, p_max, ref dtvals, ref xerr);

            MilsteinStrong.milstrong_gnuplot(p_max, dtvals, xerr);
        }

        static void test08()

            //****************************************************************************
            //
            //  Purpose:
            //
            //    TEST08 tests STAB_ASYMPTOTIC.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    28 September 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int n = 1000;
            int p_max = 3;
            int seed;

            Console.WriteLine("");
            Console.WriteLine("TEST08:");
            Console.WriteLine("  STAB_ASYMPTOTIC investigates the asymptotic");
            Console.WriteLine("  stability of the Euler-Maruyama method.");
            Console.WriteLine("");
            Console.WriteLine("  For technical reasons, the plotting is done");
            Console.WriteLine("  in the same routine as the computations.");

            seed = 123456789;

            Stability.stab_asymptotic(ref seed, n, p_max);
        }

        static void test09()

            //****************************************************************************
            //
            //  Purpose:
            //
            //    TEST09 tests STAB_MEANSQUARE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    28 September 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int seed;

            Console.WriteLine("");
            Console.WriteLine("TEST09:");
            Console.WriteLine("  STAB_MEANSQUARE investigates the mean square");
            Console.WriteLine("  stability of the Euler-Maruyama method.");
            Console.WriteLine("");
            Console.WriteLine("  For technical reasons, the plotting is done");
            Console.WriteLine("  in the same routine as the computations.");

            seed = 123456789;

            Stability.stab_meansquare(seed);

        }

        static void test10()

            //****************************************************************************
            //
            //  Purpose:
            //
            //    TEST10 tests STOCHASTIC_INTEGRAL_ITO.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    28 September 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double error = 0;
            double estimate = 0;
            double exact = 0;
            int i;
            int n;
            int seed;

            Console.WriteLine("");
            Console.WriteLine("TEST10:");
            Console.WriteLine("  Estimate the Ito integral of W(t) dW over [0,1].");
            Console.WriteLine("");
            Console.WriteLine("                                                 Abs          Rel");
            Console.WriteLine("         N        Exact        Estimate          Error        Error");
            Console.WriteLine("");

            n = 100;
            seed = 123456789;

            for (i = 1; i <= 7; i++)
            {
                Integrals.stochastic_integral_ito(n, ref seed, ref estimate, ref exact, ref error);

                Console.WriteLine("  " + n.ToString().PadLeft(8)
                                       + "  " + exact.ToString().PadLeft(16)
                                       + "  " + estimate.ToString().PadLeft(16)
                                       + "  " + error.ToString().PadLeft(16)
                                       + "  " + (error / exact).ToString().PadLeft(16) + "");

                n = n * 4;
            }
        }

        static void test11()

            //****************************************************************************
            //
            //  Purpose:
            //
            //    TEST11 tests STOCHASTIC_INTEGRAL_STRAT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    28 September 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double error = 0;
            double estimate = 0;
            double exact = 0;
            int i;
            int n;
            int seed;

            Console.WriteLine("");
            Console.WriteLine("TEST11:");
            Console.WriteLine("  Estimate the Stratonovich integral of W(t) dW over [0,1].");
            Console.WriteLine("");
            Console.WriteLine("                                                 Abs          Rel");
            Console.WriteLine("         N        Exact        Estimate          Error        Error");
            Console.WriteLine("");

            n = 100;
            seed = 123456789;

            for (i = 1; i <= 7; i++)
            {
                Integrals.stochastic_integral_strat(n, ref seed, ref estimate, ref exact, ref error);

                Console.WriteLine("  " + n.ToString().PadLeft(8)
                                       + "  " + exact.ToString().PadLeft(16)
                                       + "  " + estimate.ToString().PadLeft(16)
                                       + "  " + error.ToString().PadLeft(16)
                                       + "  " + (error / exact).ToString().PadLeft(16) + "");

                n = n * 4;
            }
        }
    }
}