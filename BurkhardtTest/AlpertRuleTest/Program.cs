using System;
using Burkardt;
using Burkardt.Types;
using Burkardt.Uniform;

namespace AlpertRuleTest
{
    class Program
    {
        static void Main(string[] args)
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

        static void alpert_log_test()

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
            int a_l;
            int a_r;
            double[] f1;
            double[] f2;
            double[] f3;
            double h;
            int i;
            int j_l;
            int j_r;
            int n;
            int nlog;
            int num_l;
            int num_r;
            int order_l;
            int order_r;
            int rule;
            double s1;
            double s2;
            double s3;
            double v1;
            double v2;
            double[] w_l;
            double[] w_r;
            double[] x_l;
            double[] x_r;
            double[] x1;
            double[] x2;
            double[] x3;

            Console.WriteLine("");
            Console.WriteLine("ALPERT_LOG_TEST");
            Console.WriteLine("  Test the Alpert rule on the log integrand.");
            Console.WriteLine("");
            Console.WriteLine(
                "  Rule  Order   J   A        N     N+2J               H        Estimate           Error");
            Console.WriteLine("");

            v2 = AlpertRule.integral_log();

            num_l = AlpertRule.num_log();
            //
            //  For the righthand interval, use the regular rule of the same index.
            //
            for (rule = 1; rule <= num_l; rule++)
            {
                a_l = AlpertRule.a_log(rule);
                j_l = AlpertRule.j_log(rule);
                order_l = AlpertRule.order_log(rule);
                x_l = new double[j_l];
                w_l = new double[j_l];
                AlpertRule.rule_log(rule, j_l, ref x_l, ref w_l);
                x1 = new double[j_l];

                a_r = AlpertRule.a_regular(rule);
                j_r = AlpertRule.j_regular(rule);
                order_r = AlpertRule.order_regular(rule);
                x_r = new double[j_r];
                w_r = new double[j_r];
                AlpertRule.rule_regular(rule, j_r, ref x_r, ref w_r);
                x3 = new double[j_r];

                n = 8;

                for (nlog = 4; nlog <= 7; nlog++)
                {
                    n = n * 2;
                    h = 1.0 / (double) (n + a_l + a_r - 1);

                    for (i = 0; i < j_l; i++)
                    {
                        x1[i] = h * x_l[i];
                    }

                    f1 = AlpertRule.integrand_log(j_l, x1);
                    s1 = typeMethods.r8vec_dot_product(j_l, w_l, f1);
                    x2 = typeMethods.r8vec_linspace_new(n, a_l * h, (a_l + n - 1) * h);
                    f2 = AlpertRule.integrand_log(n, x2);
                    s2 = typeMethods.r8vec_sum(n, f2);
                    for (i = 0; i < j_r; i++)
                    {
                        x3[i] = 1.0 - h * x_r[i];
                    }

                    f3 = AlpertRule.integrand_log(j_r, x3);
                    s3 = typeMethods.r8vec_dot_product(j_r, w_r, f3);

                    v1 = h * (s1 + s2 + s3);

                    Console.WriteLine("  " + rule.ToString().PadLeft(2)
                                           + "  " + order_l.ToString().PadLeft(4)
                                           + "  " + j_l.ToString().PadLeft(2)
                                           + "  " + a_l.ToString().PadLeft(2)
                                           + "  " + n.ToString().PadLeft(7)
                                           + "  " + (n + j_l + j_r).ToString().PadLeft(7)
                                           + "  " + h.ToString().PadLeft(14)
                                           + "  " + v1.ToString().PadLeft(14)
                                           + "  " + Math.Abs(v1 - v2).ToString().PadLeft(14) + "");

                }

                Console.WriteLine("");

            }

            Console.WriteLine("");
            Console.WriteLine("                                                Exact:"
                              + v2.ToString().PadLeft(14) + "");

        }

        static void alpert_power_test()

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
            int a_p;
            int a_r;
            double[] f1;
            double[] f2;
            double[] f3;
            double h;
            int i;
            int j_p;
            int j_r;
            int n;
            int nlog;
            int num_p;
            int num_r;
            double order_p;
            int order_r;
            int rule;
            double s1;
            double s2;
            double s3;
            double v1;
            double v2;
            double[] w_p;
            double[] w_r;
            double[] x_p;
            double[] x_r;
            double[] x1;
            double[] x2;
            double[] x3;

            Console.WriteLine("");
            Console.WriteLine("ALPERT_POWER_TEST");
            Console.WriteLine("  Test the Alpert rule on the power integrand.");
            Console.WriteLine("");
            Console.WriteLine(
                "  Rule  Order   J   A        N     N+2J               H        Estimate           Error");
            Console.WriteLine("");

            v2 = AlpertRule.integral_power();

            num_p = AlpertRule.num_power();
            //
            //  For the righthand interval, use the regular rule of the same index.
            //
            for (rule = 1; rule <= num_p; rule++)
            {
                a_p = AlpertRule.a_power(rule);
                j_p = AlpertRule.j_power(rule);
                order_p = AlpertRule.order_power(rule);
                x_p = new double[j_p];
                w_p = new double[j_p];
                AlpertRule.rule_power(rule, j_p, ref x_p, ref w_p);

                x1 = new double[j_p];

                a_r = AlpertRule.a_regular(rule);
                j_r = AlpertRule.j_regular(rule);
                order_r = AlpertRule.order_regular(rule);
                x_r = new double[j_r];
                w_r = new double[j_r];
                AlpertRule.rule_regular(rule, j_r, ref x_r, ref w_r);

                x3 = new double[j_r];

                n = 8;

                for (nlog = 4; nlog <= 6; nlog++)
                {
                    n = n * 2;
                    h = 1.0 / (double) (n + a_p + a_r - 1);

                    for (i = 0; i < j_p; i++)
                    {
                        x1[i] = h * x_p[i];
                    }

                    f1 = AlpertRule.integrand_power(j_p, x1);
                    s1 = typeMethods.r8vec_dot_product(j_p, w_p, f1);

                    x2 = typeMethods.r8vec_linspace_new(n, a_p * h, (a_p + n - 1) * h);
                    f2 = AlpertRule.integrand_power(n, x2);
                    s2 = typeMethods.r8vec_sum(n, f2);

                    for (i = 0; i < j_r; i++)
                    {
                        x3[i] = 1.0 - h * x_r[i];
                    }

                    f3 = AlpertRule.integrand_power(j_r, x3);
                    s3 = typeMethods.r8vec_dot_product(j_r, w_r, f3);

                    v1 = h * (s1 + s2 + s3);

                    Console.WriteLine("  " + rule.ToString().PadLeft(2)
                                           + "  " + order_p.ToString().PadLeft(4)
                                           + "  " + j_p.ToString().PadLeft(2)
                                           + "  " + a_p.ToString().PadLeft(2)
                                           + "  " + n.ToString().PadLeft(7)
                                           + "  " + (n + j_p + j_r).ToString().PadLeft(7)
                                           + "  " + h.ToString().PadLeft(14)
                                           + "  " + v1.ToString().PadLeft(14)
                                           + "  " + Math.Abs(v1 - v2).ToString().PadLeft(14) + "");
                }

                Console.WriteLine("");
            }

            Console.WriteLine("");
            Console.WriteLine("                                                Exact:"
                              + v2.ToString().PadLeft(14) + "");

        }

        static void alpert_regular_test()

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
            int a;
            double[] f1;
            double[] f2;
            double[] f3;
            double h;
            int i;
            int j;
            int n;
            int nlog;
            int num;
            int order;
            int rule;
            double s1;
            double s2;
            double s3;
            double v1;
            double v2;
            double[] w;
            double[] x;
            double[] x1;
            double[] x2;
            double[] x3;

            Console.WriteLine("");
            Console.WriteLine("ALPERT_REGULAR_TEST");
            Console.WriteLine("  Test the Alpert rule on the regular integrand.");
            Console.WriteLine("");
            Console.WriteLine(
                "  Rule  Order   J   A        N     N+2J               H        Estimate           Error");
            Console.WriteLine("");

            v2 = AlpertRule.integral_regular();

            num = AlpertRule.num_regular();

            for (rule = 1; rule <= num; rule++)
            {
                a = AlpertRule.a_regular(rule);
                j = AlpertRule.j_regular(rule);
                order = AlpertRule.order_regular(rule);
                x = new double[j];
                w = new double[j];
                AlpertRule.rule_regular(rule, j, ref x, ref w);

                x1 = new double[j];
                x3 = new double[j];

                n = 8;

                for (nlog = 4; nlog <= 6; nlog++)
                {
                    n = n * 2;
                    h = 1.0 / (double) (n + 2 * a - 1);

                    for (i = 0; i < j; i++)
                    {
                        x1[i] = h * x[i];
                    }

                    f1 = AlpertRule.integrand_regular(j, x1);
                    s1 = typeMethods.r8vec_dot_product(j, w, f1);

                    x2 = typeMethods.r8vec_linspace_new(n, a * h, (a + n - 1) * h);
                    f2 = AlpertRule.integrand_regular(n, x2);
                    s2 = typeMethods.r8vec_sum(n, f2);

                    for (i = 0; i < j; i++)
                    {
                        x3[i] = 1.0 - h * x[i];
                    }

                    f3 = AlpertRule.integrand_regular(j, x3);
                    s3 = typeMethods.r8vec_dot_product(j, w, f3);

                    v1 = h * (s1 + s2 + s3);

                    Console.WriteLine("  " + rule.ToString().PadLeft(2)
                                           + "  " + order.ToString().PadLeft(4)
                                           + "  " + j.ToString().PadLeft(2)
                                           + "  " + a.ToString().PadLeft(2)
                                           + "  " + n.ToString().PadLeft(7)
                                           + "  " + (n + 2 * j).ToString().PadLeft(7)
                                           + "  " + h.ToString().PadLeft(14)
                                           + "  " + v1.ToString().PadLeft(14)
                                           + "  " + Math.Abs(v1 - v2).ToString().PadLeft(14) + "");
                }

                Console.WriteLine("");
            }

            Console.WriteLine("");
            Console.WriteLine("                                                Exact:"
                              + v2.ToString().PadLeft(14) + "");

        }

        static void monte_carlo_log_test()

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
            double[] f;
            double h;
            int n;
            int nlog;
            int seed;
            double v1;
            double v2;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("MONTE_CARLO_LOG_TEST");
            Console.WriteLine("  Test the Monte Carlo rule on the log singular integrand.");
            Console.WriteLine("");
            Console.WriteLine("          N        Estimate           Error");
            Console.WriteLine("");

            v2 = AlpertRule.integral_log();

            seed = 123456789;

            n = 17;

            for (nlog = 5; nlog <= 20; nlog++)
            {
                n = (n - 1) * 2 + 1;
                h = 1.0 / (double) (n);
                x = UniformRNG.r8vec_uniform_01_new(n, ref seed);
                f = AlpertRule.integrand_log(n, x);
                v1 = h * typeMethods.r8vec_sum(n, f);
                Console.WriteLine("  " + n.ToString().PadLeft(9)
                                       + "  " + v1.ToString().PadLeft(14)
                                       + "  " + Math.Abs(v1 - v2).ToString().PadLeft(14) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("      Exact: " + v2.ToString().PadLeft(14) + "");

        }

        static void monte_carlo_power_test()

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
            double[] f;
            double h;
            int n;
            int nlog;
            int seed;
            double v1;
            double v2;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("MONTE_CARLO_POWER_TEST");
            Console.WriteLine("  Test the Monte Carlo rule on the power singular integrand.");
            Console.WriteLine("");
            Console.WriteLine("          N        Estimate           Error");
            Console.WriteLine("");

            v2 = AlpertRule.integral_power();

            seed = 123456789;

            n = 17;

            for (nlog = 5; nlog <= 20; nlog++)
            {
                n = (n - 1) * 2 + 1;
                h = 1.0 / (double) (n);
                x = UniformRNG.r8vec_uniform_01_new(n, ref seed);
                f = AlpertRule.integrand_power(n, x);
                v1 = h * typeMethods.r8vec_sum(n, f);
                Console.WriteLine("  " + n.ToString().PadLeft(9)
                                       + "  " + v1.ToString().PadLeft(14)
                                       + "  " + Math.Abs(v1 - v2).ToString().PadLeft(14) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("      Exact: " + v2.ToString().PadLeft(14) + "");

        }

        static void monte_carlo_regular_test()

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
            double[] f;
            double h;
            int n;
            int nlog;
            int seed;
            double v1;
            double v2;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("MONTE_CARLO_REGULAR_TEST");
            Console.WriteLine("  Test the Monte Carlo rule on the regular integrand.");
            Console.WriteLine("");
            Console.WriteLine("          N        Estimate           Error");
            Console.WriteLine("");

            v2 = AlpertRule.integral_regular();

            seed = 123456789;

            n = 17;

            for (nlog = 5; nlog <= 20; nlog++)
            {
                n = (n - 1) * 2 + 1;
                h = 1.0 / (double) (n);
                x = UniformRNG.r8vec_uniform_01_new(n, ref seed);
                f = AlpertRule.integrand_regular(n, x);
                v1 = h * typeMethods.r8vec_sum(n, f);
                Console.WriteLine("  " + n.ToString().PadLeft(9)
                                       + "  " + v1.ToString().PadLeft(14)
                                       + "  " + Math.Abs(v1 - v2).ToString().PadLeft(14) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("      Exact: " + v2.ToString().PadLeft(14) + "");

        }

        static void trapezoid_log_test()

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
            double[] f;
            double h;
            int n;
            int nlog;
            double v1;
            double v2;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TRAPEZOID_LOG_TEST");
            Console.WriteLine("  Test the trapezoidal rule on the log-singular integrand.");
            Console.WriteLine("");
            Console.WriteLine("        N        Estimate           Error");
            Console.WriteLine("");

            v2 = AlpertRule.integral_log();

            n = 17;

            for (nlog = 5; nlog <= 12; nlog++)
            {
                n = (n - 1) * 2 + 1;
                h = 1.0 / (double) (n - 1);
                x = typeMethods.r8vec_linspace_new(n, 0.0, 1.0);
                x[0] = 0.5 * (x[0] + x[1]);
                f = AlpertRule.integrand_log(n, x);
                v1 = h * (typeMethods.r8vec_sum(n, f) - 0.5 * (f[0] + f[n - 1]));
                Console.WriteLine("  " + n.ToString().PadLeft(7)
                                       + "  " + v1.ToString().PadLeft(14)
                                       + "  " + Math.Abs(v1 - v2).ToString().PadLeft(14) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("    Exact: " + v2.ToString().PadLeft(14) + "");

        }

        static void trapezoid_power_test()

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
            double[] f;
            double h;
            int n;
            int nlog;
            double v1;
            double v2;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TRAPEZOID_POWER_TEST");
            Console.WriteLine("  Test the trapezoidal rule on the power-singular integrand.");
            Console.WriteLine("");
            Console.WriteLine("        N        Estimate           Error");
            Console.WriteLine("");

            v2 = AlpertRule.integral_power();

            n = 17;

            for (nlog = 5; nlog <= 12; nlog++)
            {
                n = (n - 1) * 2 + 1;
                h = 1.0 / (double) (n - 1);
                x = typeMethods.r8vec_linspace_new(n, 0.0, 1.0);
                x[0] = 0.5 * (x[0] + x[1]);
                f = AlpertRule.integrand_power(n, x);
                v1 = h * (typeMethods.r8vec_sum(n, f) - 0.5 * (f[0] + f[n - 1]));
                Console.WriteLine("  " + n.ToString().PadLeft(7)
                                       + "  " + v1.ToString().PadLeft(14)
                                       + "  " + Math.Abs(v1 - v2).ToString().PadLeft(14) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("    Exact: " + v2.ToString().PadLeft(14) + "");

        }

        static void trapezoid_regular_test()

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
            double[] f;
            double h;
            int n;
            int nlog;
            double v1;
            double v2;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TRAPEZOID_REGULAR_TEST");
            Console.WriteLine("  Test the trapezoidal rule on the regular integrand.");
            Console.WriteLine("");
            Console.WriteLine("        N        Estimate           Error");
            Console.WriteLine("");

            v2 = AlpertRule.integral_regular();

            n = 17;

            for (nlog = 5; nlog <= 12; nlog++)
            {
                n = (n - 1) * 2 + 1;
                h = 1.0 / (double) (n - 1);
                x = typeMethods.r8vec_linspace_new(n, 0.0, 1.0);
                f = AlpertRule.integrand_regular(n, x);
                v1 = h * (typeMethods.r8vec_sum(n, f) - 0.5 * (f[0] + f[n - 1]));
                Console.WriteLine("  " + n.ToString().PadLeft(7)
                                       + "  " + v1.ToString().PadLeft(14)
                                       + "  " + Math.Abs(v1 - v2).ToString().PadLeft(14) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("    Exact: " + v2.ToString().PadLeft(14) + "");

        }

    }
}