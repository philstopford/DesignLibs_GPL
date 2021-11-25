using System;
using System.Globalization;
using Burkardt;
using Burkardt.Types;
using Burkardt.Uniform;

namespace AlpertRuleTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ALPERT_RULE_TEST tests the ALPERT_RULE library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 December 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("ALPERT_RULE_TEST");
            
        Console.WriteLine("  Test the ALPERT_RULE library.");

        monte_carlo_regular_test();
        monte_carlo_log_test();
        monte_carlo_power_test();

        trapezoid_regular_test();
        trapezoid_log_test();
        trapezoid_power_test();

        alpert_regular_test();
        alpert_log_test();
        alpert_power_test();
        Console.WriteLine("");
        Console.WriteLine("ALPERT_RULE_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void alpert_log_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ALPERT_LOG_TEST tests the Alpert rule on the log integrand.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 December 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int rule;

        Console.WriteLine("");
        Console.WriteLine("ALPERT_LOG_TEST");
        Console.WriteLine("  Test the Alpert rule on the log integrand.");
        Console.WriteLine("");
        Console.WriteLine(
            "  Rule  Order   J   A        N     N+2J               H        Estimate           Error");
        Console.WriteLine("");

        double v2 = AlpertRule.integral_log();

        int num_l = AlpertRule.num_log();
        //
        //  For the righthand interval, use the regular rule of the same index.
        //
        for (rule = 1; rule <= num_l; rule++)
        {
            int a_l = AlpertRule.a_log(rule);
            int j_l = AlpertRule.j_log(rule);
            int order_l = AlpertRule.order_log(rule);
            double[] x_l = new double[j_l];
            double[] w_l = new double[j_l];
            AlpertRule.rule_log(rule, j_l, ref x_l, ref w_l);
            double[] x1 = new double[j_l];

            int a_r = AlpertRule.a_regular(rule);
            int j_r = AlpertRule.j_regular(rule);
            AlpertRule.order_regular(rule);
            double[] x_r = new double[j_r];
            double[] w_r = new double[j_r];
            AlpertRule.rule_regular(rule, j_r, ref x_r, ref w_r);
            double[] x3 = new double[j_r];

            int n = 8;

            int nlog;
            for (nlog = 4; nlog <= 7; nlog++)
            {
                n *= 2;
                double h = 1.0 / (n + a_l + a_r - 1);

                int i;
                for (i = 0; i < j_l; i++)
                {
                    x1[i] = h * x_l[i];
                }

                double[] f1 = AlpertRule.integrand_log(j_l, x1);
                double s1 = typeMethods.r8vec_dot_product(j_l, w_l, f1);
                double[] x2 = typeMethods.r8vec_linspace_new(n, a_l * h, (a_l + n - 1) * h);
                double[] f2 = AlpertRule.integrand_log(n, x2);
                double s2 = typeMethods.r8vec_sum(n, f2);
                for (i = 0; i < j_r; i++)
                {
                    x3[i] = 1.0 - h * x_r[i];
                }

                double[] f3 = AlpertRule.integrand_log(j_r, x3);
                double s3 = typeMethods.r8vec_dot_product(j_r, w_r, f3);

                double v1 = h * (s1 + s2 + s3);

                Console.WriteLine("  " + rule.ToString().PadLeft(2)
                                       + "  " + order_l.ToString().PadLeft(4)
                                       + "  " + j_l.ToString().PadLeft(2)
                                       + "  " + a_l.ToString().PadLeft(2)
                                       + "  " + n.ToString().PadLeft(7)
                                       + "  " + (n + j_l + j_r).ToString().PadLeft(7)
                                       + "  " + h.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + v1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + Math.Abs(v1 - v2).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

            }

            Console.WriteLine("");

        }

        Console.WriteLine("");
        Console.WriteLine("                                                Exact:"
                          + v2.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

    }

    private static void alpert_power_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ALPERT_POWER_TEST tests the Alpert rule on the power integrand.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 December 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int rule;

        Console.WriteLine("");
        Console.WriteLine("ALPERT_POWER_TEST");
        Console.WriteLine("  Test the Alpert rule on the power integrand.");
        Console.WriteLine("");
        Console.WriteLine(
            "  Rule  Order   J   A        N     N+2J               H        Estimate           Error");
        Console.WriteLine("");

        double v2 = AlpertRule.integral_power();

        int num_p = AlpertRule.num_power();
        //
        //  For the righthand interval, use the regular rule of the same index.
        //
        for (rule = 1; rule <= num_p; rule++)
        {
            int a_p = AlpertRule.a_power(rule);
            int j_p = AlpertRule.j_power(rule);
            double order_p = AlpertRule.order_power(rule);
            double[] x_p = new double[j_p];
            double[] w_p = new double[j_p];
            AlpertRule.rule_power(rule, j_p, ref x_p, ref w_p);

            double[] x1 = new double[j_p];

            int a_r = AlpertRule.a_regular(rule);
            int j_r = AlpertRule.j_regular(rule);
            AlpertRule.order_regular(rule);
            double[] x_r = new double[j_r];
            double[] w_r = new double[j_r];
            AlpertRule.rule_regular(rule, j_r, ref x_r, ref w_r);

            double[] x3 = new double[j_r];

            int n = 8;

            int nlog;
            for (nlog = 4; nlog <= 6; nlog++)
            {
                n *= 2;
                double h = 1.0 / (n + a_p + a_r - 1);

                int i;
                for (i = 0; i < j_p; i++)
                {
                    x1[i] = h * x_p[i];
                }

                double[] f1 = AlpertRule.integrand_power(j_p, x1);
                double s1 = typeMethods.r8vec_dot_product(j_p, w_p, f1);

                double[] x2 = typeMethods.r8vec_linspace_new(n, a_p * h, (a_p + n - 1) * h);
                double[] f2 = AlpertRule.integrand_power(n, x2);
                double s2 = typeMethods.r8vec_sum(n, f2);

                for (i = 0; i < j_r; i++)
                {
                    x3[i] = 1.0 - h * x_r[i];
                }

                double[] f3 = AlpertRule.integrand_power(j_r, x3);
                double s3 = typeMethods.r8vec_dot_product(j_r, w_r, f3);

                double v1 = h * (s1 + s2 + s3);

                Console.WriteLine("  " + rule.ToString().PadLeft(2)
                                       + "  " + order_p.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                       + "  " + j_p.ToString().PadLeft(2)
                                       + "  " + a_p.ToString().PadLeft(2)
                                       + "  " + n.ToString().PadLeft(7)
                                       + "  " + (n + j_p + j_r).ToString().PadLeft(7)
                                       + "  " + h.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + v1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + Math.Abs(v1 - v2).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }

            Console.WriteLine("");
        }

        Console.WriteLine("");
        Console.WriteLine("                                                Exact:"
                          + v2.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

    }

    private static void alpert_regular_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ALPERT_REGULAR_TEST tests the Alpert rule on the regular integrand.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 December 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int rule;

        Console.WriteLine("");
        Console.WriteLine("ALPERT_REGULAR_TEST");
        Console.WriteLine("  Test the Alpert rule on the regular integrand.");
        Console.WriteLine("");
        Console.WriteLine(
            "  Rule  Order   J   A        N     N+2J               H        Estimate           Error");
        Console.WriteLine("");

        double v2 = AlpertRule.integral_regular();

        int num = AlpertRule.num_regular();

        for (rule = 1; rule <= num; rule++)
        {
            int a = AlpertRule.a_regular(rule);
            int j = AlpertRule.j_regular(rule);
            int order = AlpertRule.order_regular(rule);
            double[] x = new double[j];
            double[] w = new double[j];
            AlpertRule.rule_regular(rule, j, ref x, ref w);

            double[] x1 = new double[j];
            double[] x3 = new double[j];

            int n = 8;

            int nlog;
            for (nlog = 4; nlog <= 6; nlog++)
            {
                n *= 2;
                double h = 1.0 / (n + 2 * a - 1);

                int i;
                for (i = 0; i < j; i++)
                {
                    x1[i] = h * x[i];
                }

                double[] f1 = AlpertRule.integrand_regular(j, x1);
                double s1 = typeMethods.r8vec_dot_product(j, w, f1);

                double[] x2 = typeMethods.r8vec_linspace_new(n, a * h, (a + n - 1) * h);
                double[] f2 = AlpertRule.integrand_regular(n, x2);
                double s2 = typeMethods.r8vec_sum(n, f2);

                for (i = 0; i < j; i++)
                {
                    x3[i] = 1.0 - h * x[i];
                }

                double[] f3 = AlpertRule.integrand_regular(j, x3);
                double s3 = typeMethods.r8vec_dot_product(j, w, f3);

                double v1 = h * (s1 + s2 + s3);

                Console.WriteLine("  " + rule.ToString().PadLeft(2)
                                       + "  " + order.ToString().PadLeft(4)
                                       + "  " + j.ToString().PadLeft(2)
                                       + "  " + a.ToString().PadLeft(2)
                                       + "  " + n.ToString().PadLeft(7)
                                       + "  " + (n + 2 * j).ToString().PadLeft(7)
                                       + "  " + h.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + v1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                       + "  " + Math.Abs(v1 - v2).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
            }

            Console.WriteLine("");
        }

        Console.WriteLine("");
        Console.WriteLine("                                                Exact:"
                          + v2.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

    }

    private static void monte_carlo_log_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MONTE_CARLO_LOG_TEST tests the Monte Carlo rule on the log singular integrand.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 December 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int nlog;

        Console.WriteLine("");
        Console.WriteLine("MONTE_CARLO_LOG_TEST");
        Console.WriteLine("  Test the Monte Carlo rule on the log singular integrand.");
        Console.WriteLine("");
        Console.WriteLine("          N        Estimate           Error");
        Console.WriteLine("");

        double v2 = AlpertRule.integral_log();

        int seed = 123456789;

        int n = 17;

        for (nlog = 5; nlog <= 20; nlog++)
        {
            n = (n - 1) * 2 + 1;
            double h = 1.0 / n;
            double[] x = UniformRNG.r8vec_uniform_01_new(n, ref seed);
            double[] f = AlpertRule.integrand_log(n, x);
            double v1 = h * typeMethods.r8vec_sum(n, f);
            Console.WriteLine("  " + n.ToString().PadLeft(9)
                                   + "  " + v1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + Math.Abs(v1 - v2).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("      Exact: " + v2.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

    }

    private static void monte_carlo_power_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MONTE_CARLO_POWER_TEST tests the Monte Carlo rule on the power singular integrand.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 December 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int nlog;

        Console.WriteLine("");
        Console.WriteLine("MONTE_CARLO_POWER_TEST");
        Console.WriteLine("  Test the Monte Carlo rule on the power singular integrand.");
        Console.WriteLine("");
        Console.WriteLine("          N        Estimate           Error");
        Console.WriteLine("");

        double v2 = AlpertRule.integral_power();

        int seed = 123456789;

        int n = 17;

        for (nlog = 5; nlog <= 20; nlog++)
        {
            n = (n - 1) * 2 + 1;
            double h = 1.0 / n;
            double[] x = UniformRNG.r8vec_uniform_01_new(n, ref seed);
            double[] f = AlpertRule.integrand_power(n, x);
            double v1 = h * typeMethods.r8vec_sum(n, f);
            Console.WriteLine("  " + n.ToString().PadLeft(9)
                                   + "  " + v1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + Math.Abs(v1 - v2).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("      Exact: " + v2.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

    }

    private static void monte_carlo_regular_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MONTE_CARLO_REGULAR_TEST tests the Monte Carlo rule on the regular integrand.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 December 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int nlog;

        Console.WriteLine("");
        Console.WriteLine("MONTE_CARLO_REGULAR_TEST");
        Console.WriteLine("  Test the Monte Carlo rule on the regular integrand.");
        Console.WriteLine("");
        Console.WriteLine("          N        Estimate           Error");
        Console.WriteLine("");

        double v2 = AlpertRule.integral_regular();

        int seed = 123456789;

        int n = 17;

        for (nlog = 5; nlog <= 20; nlog++)
        {
            n = (n - 1) * 2 + 1;
            double h = 1.0 / n;
            double[] x = UniformRNG.r8vec_uniform_01_new(n, ref seed);
            double[] f = AlpertRule.integrand_regular(n, x);
            double v1 = h * typeMethods.r8vec_sum(n, f);
            Console.WriteLine("  " + n.ToString().PadLeft(9)
                                   + "  " + v1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + Math.Abs(v1 - v2).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("      Exact: " + v2.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

    }

    private static void trapezoid_log_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRAPEZOID_LOG_TEST tests the trapezoid rule on the log-singular integrand.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 December 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int nlog;

        Console.WriteLine("");
        Console.WriteLine("TRAPEZOID_LOG_TEST");
        Console.WriteLine("  Test the trapezoidal rule on the log-singular integrand.");
        Console.WriteLine("");
        Console.WriteLine("        N        Estimate           Error");
        Console.WriteLine("");

        double v2 = AlpertRule.integral_log();

        int n = 17;

        for (nlog = 5; nlog <= 12; nlog++)
        {
            n = (n - 1) * 2 + 1;
            double h = 1.0 / (n - 1);
            double[] x = typeMethods.r8vec_linspace_new(n, 0.0, 1.0);
            x[0] = 0.5 * (x[0] + x[1]);
            double[] f = AlpertRule.integrand_log(n, x);
            double v1 = h * (typeMethods.r8vec_sum(n, f) - 0.5 * (f[0] + f[n - 1]));
            Console.WriteLine("  " + n.ToString().PadLeft(7)
                                   + "  " + v1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + Math.Abs(v1 - v2).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("    Exact: " + v2.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

    }

    private static void trapezoid_power_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRAPEZOID_POWER_TEST tests the trapezoid rule on the power-singular integrand.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 December 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int nlog;

        Console.WriteLine("");
        Console.WriteLine("TRAPEZOID_POWER_TEST");
        Console.WriteLine("  Test the trapezoidal rule on the power-singular integrand.");
        Console.WriteLine("");
        Console.WriteLine("        N        Estimate           Error");
        Console.WriteLine("");

        double v2 = AlpertRule.integral_power();

        int n = 17;

        for (nlog = 5; nlog <= 12; nlog++)
        {
            n = (n - 1) * 2 + 1;
            double h = 1.0 / (n - 1);
            double[] x = typeMethods.r8vec_linspace_new(n, 0.0, 1.0);
            x[0] = 0.5 * (x[0] + x[1]);
            double[] f = AlpertRule.integrand_power(n, x);
            double v1 = h * (typeMethods.r8vec_sum(n, f) - 0.5 * (f[0] + f[n - 1]));
            Console.WriteLine("  " + n.ToString().PadLeft(7)
                                   + "  " + v1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + Math.Abs(v1 - v2).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("    Exact: " + v2.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

    }

    private static void trapezoid_regular_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRAPEZOID_REGULAR_TEST tests the trapezoid rule on the regular integrand.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 December 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int nlog;

        Console.WriteLine("");
        Console.WriteLine("TRAPEZOID_REGULAR_TEST");
        Console.WriteLine("  Test the trapezoidal rule on the regular integrand.");
        Console.WriteLine("");
        Console.WriteLine("        N        Estimate           Error");
        Console.WriteLine("");

        double v2 = AlpertRule.integral_regular();

        int n = 17;

        for (nlog = 5; nlog <= 12; nlog++)
        {
            n = (n - 1) * 2 + 1;
            double h = 1.0 / (n - 1);
            double[] x = typeMethods.r8vec_linspace_new(n, 0.0, 1.0);
            double[] f = AlpertRule.integrand_regular(n, x);
            double v1 = h * (typeMethods.r8vec_sum(n, f) - 0.5 * (f[0] + f[n - 1]));
            Console.WriteLine("  " + n.ToString().PadLeft(7)
                                   + "  " + v1.ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + Math.Abs(v1 - v2).ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("    Exact: " + v2.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");

    }

}