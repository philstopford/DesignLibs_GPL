using System;
using Burkardt.Chebyshev;
using Burkardt.ClenshawCurtisNS;
using Burkardt.IntegralNS;
using Burkardt.Interpolation;
using Burkardt.PointsNS;
using Burkardt.Quadrature;
using Burkardt.Types;
using Burkardt.Uniform;

namespace SandiaRulesTest
{
    internal class Program
    {
        private static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for SANDIA_RULES_TEST..
            //
            //  Discussion:
            //
            //    SANDIA_RULES_TEST tests the SANDIA_RULES library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    31 December 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int r;

            Console.WriteLine("");
            Console.WriteLine("SANDIA_RULES_TEST");
            Console.WriteLine("  Test the SANDIA_RULES library.");

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
            test12();
            test13();
            test14();
            test15();
            test16();
            test17();
            test18();
            test19();

            test20();
            test21();
            test22();
            //
            //  Repeat tests, but now call "NP" versions of routines.
            //
            test01_np();
            test02_np();
            test03_np();
            test04_np();
            test05_np();
            test06_np();
            test07_np();
            test08_np();
            test09_np();

            test10_np();
            test11_np();
            test12_np();
            test13_np();
            test14_np();
            test15_np();
            test16_np();
            test17_np();
            test18_np();
            test19_np();

            test20_np();
            test21_np();
            test22_np();
            //
            //  TEST23 takes an input argument of R, a rule index.
            //
            r = 1;
            test23(r);
            r = 3;
            test23(r);
            r = 4;
            test23(r);
            r = 11;
            test23(r);

            test24();
            test25();
            test26();
            test27();
            test28();
            test285();
            test29();

            test30();
            test31();
            test32();
            test33();
            test34();
            test35();
            test36();
            test37();
            test38();
            test39();
            //
            //  Terminate.
            //
            Console.WriteLine("");
            Console.WriteLine("SANDIA_RULES_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        private static void test01()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST01 tests CHEBYSHEV1_COMPUTE against CHEBYSHEV1_INTEGRAL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double error;
            double estimate;
            double exact;
            double[] f;
            int i;
            int n;
            int order;
            int order_max = 10;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST01");
            Console.WriteLine("  CHEBYSHEV1_COMPUTE computes a Gauss-Chebyshev type 1 rule");
            Console.WriteLine("  which is appropriate for integrands of the form");
            Console.WriteLine("    Integral ( -1 <= x <= +1 ) f(x) / Math.Sqrt ( 1 - x^2 ) dx.");
            Console.WriteLine("");
            Console.WriteLine("  CHEBYSHEV1_INTEGRAL determines the exact value of");
            Console.WriteLine("  this integal when f(x) = x^n.");
            Console.WriteLine("");
            Console.WriteLine("  A rule of order ORDER should be exact for monomials X^N");
            Console.WriteLine("  up to N = 2*ORDER-1");
            Console.WriteLine("");
            Console.WriteLine("  In the following table, for each order, the LAST THREE estimates");
            Console.WriteLine("  are made on monomials that exceed the exactness limit for the rule.");
            Console.WriteLine("");
            Console.WriteLine("     Order         N       Estimate       Exact            Error");

            for (order = 1; order <= order_max; order++)
            {
                Console.WriteLine("");

                f = new double[order];
                w = new double[order];
                x = new double[order];

                Chebyshev1.chebyshev1_compute(order, ref x, ref w);

                for (n = 0; n <= 2 * order + 2; n += 1)
                {
                    exact = Integral.chebyshev1_integral(n);

                    switch (n)
                    {
                        case 0:
                        {
                            for (i = 0; i < order; i++)
                            {
                                f[i] = 1.0;
                            }

                            break;
                        }
                        default:
                        {
                            for (i = 0; i < order; i++)
                            {
                                f[i] = Math.Pow(x[i], n);
                            }

                            break;
                        }
                    }

                    estimate = 0.0;
                    for (i = 0; i < order; i++)
                    {
                        estimate += w[i] * f[i];
                    }

                    error = typeMethods.r8_abs(exact - estimate);

                    Console.WriteLine("  " + order.ToString().PadLeft(8)
                                           + "  " + n.ToString().PadLeft(8)
                                           + "  " + estimate.ToString("0.######").PadLeft(14)
                                           + "  " + exact.ToString("0.######").PadLeft(14)
                                           + "  " + error.ToString("0.######").PadLeft(14) + "");
                }
            }
        }

        private static void test02()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST02 tests CHEBYSHEV1_COMPUTE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int i;
            int order;
            int order_max = 10;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST02");
            Console.WriteLine("  CHEBYSHEV1_COMPUTE computes a Gauss-Chebyshev type 1 rule");
            Console.WriteLine("  which is appropriate for integrands of the form");
            Console.WriteLine("    Integral ( -1 <= x <= +1 ) f(x) / sqrt(1-x^2) dx.");

            for (order = 1; order <= order_max; order++)
            {
                Console.WriteLine("");
                Console.WriteLine("  Order = " + order + "");

                w = new double[order];
                x = new double[order];

                Chebyshev1.chebyshev1_compute(order, ref x, ref w);

                for (i = 0; i < order; i += 1)
                {
                    Console.WriteLine("  " + i.ToString().PadLeft(8)
                                           + "  " + x[i].ToString("0.################").PadLeft(24)
                                           + "  " + w[i].ToString("0.################").PadLeft(24) + "");
                }
            }
        }

        private static void test03()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST03 tests CHEBYSHEV2_COMPUTE against CHEBYSHEV2_INTEGRAL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double error;
            double estimate;
            double exact;
            double[] f;
            int i;
            int n;
            int order;
            int order_max = 10;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST03");
            Console.WriteLine("  CHEBYSHEV2_COMPUTE computes a Gauss-Chebyshev type 2 rule");
            Console.WriteLine("  which is appropriate for integrands of the form");
            Console.WriteLine("    Integral ( -1 <= x <= +1 ) f(x) * Math.Sqrt ( 1 - x^2 ) dx.");
            Console.WriteLine("");
            Console.WriteLine("  CHEBYSHEV2_INTEGRAL determines the exact value of");
            Console.WriteLine("  this integal when f(x) = x^n.");
            Console.WriteLine("");
            Console.WriteLine("  A rule of order ORDER should be exact for monomials X^N");
            Console.WriteLine("  up to N = 2*ORDER-1");
            Console.WriteLine("");
            Console.WriteLine("  In the following table, for each order, the LAST THREE estimates");
            Console.WriteLine("  are made on monomials that exceed the exactness limit for the rule.");
            Console.WriteLine("");
            Console.WriteLine("     Order         N       Estimate       Exact            Error");

            for (order = 1; order <= order_max; order++)
            {
                Console.WriteLine("");

                f = new double[order];
                w = new double[order];
                x = new double[order];

                Chebyshev2.chebyshev2_compute(order, ref x, ref w);

                for (n = 0; n <= 2 * order + 2; n += 1)
                {
                    exact = Integral.chebyshev2_integral(n);

                    switch (n)
                    {
                        case 0:
                        {
                            for (i = 0; i < order; i++)
                            {
                                f[i] = 1.0;
                            }

                            break;
                        }
                        default:
                        {
                            for (i = 0; i < order; i++)
                            {
                                f[i] = Math.Pow(x[i], n);
                            }

                            break;
                        }
                    }

                    estimate = 0.0;
                    for (i = 0; i < order; i++)
                    {
                        estimate += w[i] * f[i];
                    }

                    error = typeMethods.r8_abs(exact - estimate);

                    Console.WriteLine("  " + order.ToString().PadLeft(8)
                                           + "  " + n.ToString().PadLeft(8)
                                           + "  " + estimate.ToString("0.######").PadLeft(14)
                                           + "  " + exact.ToString("0.######").PadLeft(14)
                                           + "  " + error.ToString("0.######").PadLeft(14) + "");
                }
            }
        }

        private static void test04()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST04 tests CHEBYSHEV2_COMPUTE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int i;
            int order;
            int order_max = 10;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST04");
            Console.WriteLine("  CHEBYSHEV2_COMPUTE computes a Gauss-Chebyshev type 2 rule");
            Console.WriteLine("  which is appropriate for integrands of the form");
            Console.WriteLine("    Integral ( -1 <= x <= +1 ) f(x) * sqrt(1-x^2) dx.");

            for (order = 1; order <= order_max; order++)
            {
                Console.WriteLine("");
                Console.WriteLine("  Order = " + order + "");

                w = new double[order];
                x = new double[order];

                Chebyshev2.chebyshev2_compute(order, ref x, ref w);

                for (i = 0; i < order; i += 1)
                {
                    Console.WriteLine("  " + i.ToString().PadLeft(8)
                                           + "  " + x[i].ToString("0.################").PadLeft(24)
                                           + "  " + w[i].ToString("0.################").PadLeft(24) + "");
                }
            }
        }

        private static void test05()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST05 tests CLENSHAW_CURTIS_COMPUTE against LEGENDRE_INTEGRAL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double error;
            double estimate;
            double exact;
            double[] f;
            int i;
            int n;
            int n_hi;
            int order;
            int order_max = 10;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST05");
            Console.WriteLine("  CLENSHAW_CURTIS_COMPUTE computes a Clenshaw Curtis rule");
            Console.WriteLine("  which is appropriate for integrands of the form");
            Console.WriteLine("    Integral ( -1 <= x <= +1 ) f(x)  dx.");
            Console.WriteLine("");
            Console.WriteLine("  LEGENDRE_INTEGRAL determines the exact value of");
            Console.WriteLine("  this integal when f(x) = x^n.");
            Console.WriteLine("");
            Console.WriteLine("  A rule of order ORDER should be exact for monomials X^N up to");
            Console.WriteLine("    N = ORDER+1 if ORDER is odd, or");
            Console.WriteLine("    N = ORDER   if ORDER is even");
            Console.WriteLine("");
            Console.WriteLine("  In the following table, for each order, the LAST THREE estimates");
            Console.WriteLine("  are made on monomials that exceed the exactness limit for the rule.");
            Console.WriteLine("");
            Console.WriteLine("     Order         N       Estimate       Exact            Error");

            for (order = 1; order <= order_max; order++)
            {
                Console.WriteLine("");

                f = new double[order];
                w = new double[order];
                x = new double[order];

                ClenshawCurtis.clenshaw_curtis_compute(order, ref x, ref w);

                n_hi = (order % 2) switch
                {
                    0 => order + 2,
                    _ => order + 3
                };

                for (n = 0; n <= n_hi; n += 1)
                {
                    exact = Integral.legendre_integral(n);

                    switch (n)
                    {
                        case 0:
                        {
                            for (i = 0; i < order; i++)
                            {
                                f[i] = 1.0;
                            }

                            break;
                        }
                        default:
                        {
                            for (i = 0; i < order; i++)
                            {
                                f[i] = Math.Pow(x[i], n);
                            }

                            break;
                        }
                    }

                    estimate = 0.0;
                    for (i = 0; i < order; i++)
                    {
                        estimate += w[i] * f[i];
                    }

                    error = typeMethods.r8_abs(exact - estimate);

                    Console.WriteLine("  " + order.ToString().PadLeft(8)
                                           + "  " + n.ToString().PadLeft(8)
                                           + "  " + estimate.ToString("0.######").PadLeft(14)
                                           + "  " + exact.ToString("0.######").PadLeft(14)
                                           + "  " + error.ToString("0.######").PadLeft(14) + "");
                }
            }
        }

        private static void test06()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST06 tests CLENSHAW_CURTIS_COMPUTE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int i;
            int order;
            int order_max = 10;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST06");
            Console.WriteLine("  CLENSHAW_CURTIS_COMPUTE computes a Clenshaw Curtis rule");
            Console.WriteLine("  which is appropriate for integrands of the form");
            Console.WriteLine("    Integral ( -1 <= x <= +1 ) f(x) dx.");

            for (order = 1; order <= order_max; order++)
            {
                Console.WriteLine("");
                Console.WriteLine("  Order = " + order + "");

                w = new double[order];
                x = new double[order];

                ClenshawCurtis.clenshaw_curtis_compute(order, ref x, ref w);

                for (i = 0; i < order; i += 1)
                {
                    Console.WriteLine("  " + i.ToString().PadLeft((8))
                                           + "  " + x[i].ToString("0.################").PadLeft(24)
                                           + "  " + w[i].ToString("0.################").PadLeft(24) + "");
                }
            }
        }

        private static void test07()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST07 tests FEJER2_COMPUTE against LEGENDRE_INTEGRAL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double error;
            double estimate;
            double exact;
            double[] f;
            int i;
            int n;
            int n_hi;
            int order;
            int order_max = 10;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST07");
            Console.WriteLine("  FEJER2_COMPUTE computes a Fejer Type 2 rule");
            Console.WriteLine("  which is appropriate for integrands of the form");
            Console.WriteLine("    Integral ( -1 <= x <= +1 ) f(x)  dx.");
            Console.WriteLine("");
            Console.WriteLine("  LEGENDRE_INTEGRAL determines the exact value of");
            Console.WriteLine("  this integal when f(x) = x^n.");
            Console.WriteLine("");
            Console.WriteLine("  A rule of order ORDER should be exact for monomials X^N up to");
            Console.WriteLine("    N = ORDER+1 if ORDER is odd, or");
            Console.WriteLine("    N = ORDER   if ORDER is even");
            Console.WriteLine("");
            Console.WriteLine("  In the following table, for each order, the LAST THREE estimates");
            Console.WriteLine("  are made on monomials that exceed the exactness limit for the rule.");
            Console.WriteLine("");
            Console.WriteLine("     Order         N       Estimate       Exact            Error");

            for (order = 1; order <= order_max; order++)
            {
                Console.WriteLine("");

                f = new double[order];
                w = new double[order];
                x = new double[order];

                Fejer2.fejer2_compute(order, ref x, ref w);

                n_hi = (order % 2) switch
                {
                    0 => order + 2,
                    _ => order + 3
                };

                for (n = 0; n <= n_hi; n += 1)
                {
                    exact = Integral.legendre_integral(n);

                    switch (n)
                    {
                        case 0:
                        {
                            for (i = 0; i < order; i++)
                            {
                                f[i] = 1.0;
                            }

                            break;
                        }
                        default:
                        {
                            for (i = 0; i < order; i++)
                            {
                                f[i] = Math.Pow(x[i], n);
                            }

                            break;
                        }
                    }

                    estimate = 0.0;
                    for (i = 0; i < order; i++)
                    {
                        estimate += w[i] * f[i];
                    }

                    error = typeMethods.r8_abs(exact - estimate);

                    Console.WriteLine("  " + order.ToString().PadLeft(8)
                                           + "  " + n.ToString().PadLeft(8)
                                           + "  " + estimate.ToString("0.######").PadLeft(14)
                                           + "  " + exact.ToString("0.######").PadLeft(14)
                                           + "  " + error.ToString("0.######").PadLeft(14) + "");
                }
            }
        }

        private static void test08()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST08 tests FEJER2_COMPUTE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int i;
            int order;
            int order_max = 10;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST08");
            Console.WriteLine("  FEJER2_COMPUTE computes a Fejer Type 2 rule");
            Console.WriteLine("  which is appropriate for integrands of the form");
            Console.WriteLine("    Integral ( -1 <= x <= +1 ) f(x) dx.");

            for (order = 1; order <= order_max; order++)
            {
                Console.WriteLine("");
                Console.WriteLine("  Order = " + order + "");

                w = new double[order];
                x = new double[order];

                Fejer2.fejer2_compute(order, ref x, ref w);

                for (i = 0; i < order; i += 1)
                {
                    Console.WriteLine("  " + i.ToString().PadLeft(8)
                                           + "  " + x[i].ToString("0.################").PadLeft(24)
                                           + "  " + w[i].ToString("0.################").PadLeft(24) + "");
                }
            }
        }

        private static void test09()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST09 tests GEGENBAUER_COMPUTE against GEGENBAUER_INTEGRAL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int TEST_NUM = 3;

            double alpha;
            double[] alpha_test = {0.5, 1.0, 2.5};
            double error;
            double estimate;
            double exact;
            double[] f;
            int i;
            int n;
            int order;
            int order_max = 10;
            int test;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST09");
            Console.WriteLine("  GEGENBAUER_COMPUTE computes a generalized Gauss-Gegenbauer rule");
            Console.WriteLine("  which is appropriate for integrands of the form");
            Console.WriteLine("    Integral ( 0 <= x < +oo ) f(x) (1-x^2)^alpha dx.");
            Console.WriteLine("");
            Console.WriteLine("  GEGENBAUER_INTEGRAL determines the exact value of");
            Console.WriteLine("  this integal when f(x) = x^n.");
            Console.WriteLine("");
            Console.WriteLine("  A rule of order ORDER should be exact for monomials X^N");
            Console.WriteLine("  up to N = 2*ORDER-1");
            Console.WriteLine("");
            Console.WriteLine("  In the following table, for each order, the LAST THREE estimates");
            Console.WriteLine("  are made on monomials that exceed the exactness limit for the rule.");
            Console.WriteLine("");
            Console.WriteLine("     Order         N       Alpha           Estimate       Exact            Error");

            for (test = 0; test < TEST_NUM; test++)
            {
                alpha = alpha_test[test];

                for (order = 1; order <= order_max; order++)
                {
                    Console.WriteLine("");

                    f = new double[order];
                    w = new double[order];
                    x = new double[order];

                    GegenbauerQuadrature.gegenbauer_compute(order, alpha, ref x, ref w);

                    for (n = 0; n <= 2 * order + 2; n += 1)
                    {
                        exact = Integral.gegenbauer_integral(n, alpha);

                        switch (n)
                        {
                            case 0:
                            {
                                for (i = 0; i < order; i++)
                                {
                                    f[i] = 1.0;
                                }

                                break;
                            }
                            default:
                            {
                                for (i = 0; i < order; i++)
                                {
                                    f[i] = Math.Pow(x[i], n);
                                }

                                break;
                            }
                        }

                        estimate = 0.0;
                        for (i = 0; i < order; i++)
                        {
                            estimate += w[i] * f[i];
                        }

                        error = typeMethods.r8_abs(exact - estimate);

                        Console.WriteLine("  " + order.ToString().PadLeft(8)
                                               + "  " + n.ToString().PadLeft(8)
                                               + "  " + alpha.ToString().PadLeft(14)
                                               + "  " + estimate.ToString("0.######").PadLeft(14)
                                               + "  " + exact.ToString("0.######").PadLeft(14)
                                               + "  " + error.ToString("0.######").PadLeft(14) + "");
                    }
                }
            }
        }

        private static void test10()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST10 tests GEGENBAUER_COMPUTE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int TEST_NUM = 3;

            double alpha;
            double[] alpha_test = {0.5, 1.0, 2.5};
            int i;
            int order;
            int order_max = 10;
            int test;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST10");
            Console.WriteLine("  GEGENBAUER_COMPUTE computes a generalized Gauss-Gegenbauer rule");
            Console.WriteLine("  which is appropriate for integrands of the form");
            Console.WriteLine("    Integral ( -1 <= x <= +1 ) f(x) (1-x^2)^alpha dx.");

            for (test = 0; test < TEST_NUM; test++)
            {
                alpha = alpha_test[test];

                for (order = 1; order <= order_max; order++)
                {
                    Console.WriteLine("");
                    Console.WriteLine("  Order = " + order + "");
                    Console.WriteLine("  ALPHA = " + alpha + "");

                    w = new double[order];
                    x = new double[order];

                    GegenbauerQuadrature.gegenbauer_compute(order, alpha, ref x, ref w);

                    for (i = 0; i < order; i += 1)
                    {
                        Console.WriteLine("  " + i.ToString().PadLeft(8)
                                               + "  " + x[i].ToString("0.################").PadLeft(24)
                                               + "  " + w[i].ToString("0.################").PadLeft(24) + "");
                    }
                }
            }
        }

        private static void test11()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST11 tests GEN_HERMITE_COMPUTE against GEN_HERMITE_INTEGRAL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int TEST_NUM = 3;

            double alpha;
            double[] alpha_test = {0.5, 1.0, 2.5};
            double error;
            double estimate;
            double exact;
            double[] f;
            int i;
            int n;
            int order;
            int order_max = 10;
            int test;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST11");
            Console.WriteLine("  GEN_HERMITE_COMPUTE computes a generalized Gauss-Hermite rule");
            Console.WriteLine("  which is appropriate for integrands of the form");
            Console.WriteLine("    Integral ( -oo < x < +oo ) f(x) x^alpha exp(-x*x) dx.");
            Console.WriteLine("");
            Console.WriteLine("  GEN_HERMITE_INTEGRAL determines the exact value of");
            Console.WriteLine("  this integal when f(x) = x^n.");
            Console.WriteLine("");
            Console.WriteLine("  A rule of order ORDER should be exact for monomials X^N");
            Console.WriteLine("  up to N = 2*ORDER-1");
            Console.WriteLine("");
            Console.WriteLine("  In the following table, for each order, the LAST THREE estimates");
            Console.WriteLine("  are made on monomials that exceed the exactness limit for the rule.");
            Console.WriteLine("");
            Console.WriteLine("     Order         N       Alpha           Estimate       Exact            Error");

            for (test = 0; test < TEST_NUM; test++)
            {
                alpha = alpha_test[test];

                for (order = 1; order <= order_max; order++)
                {
                    Console.WriteLine("");

                    f = new double[order];
                    w = new double[order];
                    x = new double[order];

                    HermiteQuadrature.gen_hermite_compute(order, alpha, ref x, ref w);

                    for (n = 0; n <= 2 * order + 2; n += 1)
                    {
                        exact = Integral.gen_hermite_integral(n, alpha);

                        switch (n)
                        {
                            case 0:
                            {
                                for (i = 0; i < order; i++)
                                {
                                    f[i] = 1.0;
                                }

                                break;
                            }
                            default:
                            {
                                for (i = 0; i < order; i++)
                                {
                                    f[i] = Math.Pow(x[i], n);
                                }

                                break;
                            }
                        }

                        estimate = 0.0;
                        for (i = 0; i < order; i++)
                        {
                            estimate += w[i] * f[i];
                        }

                        error = typeMethods.r8_abs(exact - estimate);

                        Console.WriteLine("  " + order.ToString().PadLeft(8)
                                               + "  " + n.ToString().PadLeft(8)
                                               + "  " + alpha.ToString("0.######").PadLeft(14)
                                               + "  " + estimate.ToString("0.######").PadLeft(14)
                                               + "  " + exact.ToString("0.######").PadLeft(14)
                                               + "  " + error.ToString("0.######").PadLeft(14) + "");
                    }
                }
            }
        }

        private static void test12()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST12 tests GEN_HERMITE_COMPUTE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int TEST_NUM = 3;

            double alpha;
            double[] alpha_test = {0.5, 1.0, 2.5};
            int i;
            int order;
            int order_max = 10;
            int test;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST12");
            Console.WriteLine("  GEN_HERMITE_COMPUTE computes a generalized Gauss-Hermite rule");
            Console.WriteLine("  which is appropriate for integrands of the form");
            Console.WriteLine("    Integral ( -oo < x < +oo ) f(x) x^alpha exp(-x*x) dx.");

            for (test = 0; test < TEST_NUM; test++)
            {
                alpha = alpha_test[test];

                for (order = 1; order <= order_max; order++)
                {
                    Console.WriteLine("");
                    Console.WriteLine("  Order = " + order + "");
                    Console.WriteLine("  ALPHA = " + alpha + "");

                    w = new double[order];
                    x = new double[order];

                    HermiteQuadrature.gen_hermite_compute(order, alpha, ref x, ref w);

                    for (i = 0; i < order; i += 1)
                    {
                        Console.WriteLine("  " + i.ToString().PadLeft(8)
                                               + "  " + x[i].ToString("0.################").PadLeft(24)
                                               + "  " + w[i].ToString("0.################").PadLeft(24) + "");
                    }
                }
            }
        }

        private static void test13()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST13 tests GEN_LAGUERRE_COMPUTE against GEN_LAGUERRE_INTEGRAL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int TEST_NUM = 3;

            double alpha;
            double[] alpha_test = {0.5, 1.0, 2.5};
            double error;
            double estimate;
            double exact;
            double[] f;
            int i;
            int n;
            int order;
            int order_max = 10;
            int test;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST13");
            Console.WriteLine("  GEN_LAGUERRE_COMPUTE computes a generalized Gauss-Laguerre rule");
            Console.WriteLine("  which is appropriate for integrands of the form");
            Console.WriteLine("    Integral ( 0 <= x < +oo ) f(x) x^alpha exp(-x) dx.");
            Console.WriteLine("");
            Console.WriteLine("  GEN_LAGUERRE_INTEGRAL determines the exact value of");
            Console.WriteLine("  this integal when f(x) = x^n.");
            Console.WriteLine("");
            Console.WriteLine("  A rule of order ORDER should be exact for monomials X^N");
            Console.WriteLine("  up to N = 2*ORDER-1");
            Console.WriteLine("");
            Console.WriteLine("  In the following table, for each order, the LAST THREE estimates");
            Console.WriteLine("  are made on monomials that exceed the exactness limit for the rule.");
            Console.WriteLine("");
            Console.WriteLine("     Order         N       Alpha           Estimate       Exact            Error");

            for (test = 0; test < TEST_NUM; test++)
            {
                alpha = alpha_test[test];

                for (order = 1; order <= order_max; order++)
                {
                    Console.WriteLine("");

                    f = new double[order];
                    w = new double[order];
                    x = new double[order];

                    Burkardt.Laguerre.QuadratureRule.gen_laguerre_compute(order, alpha, ref x, ref w);

                    for (n = 0; n <= 2 * order + 2; n += 1)
                    {
                        exact = Integral.gen_laguerre_integral(n, alpha);

                        switch (n)
                        {
                            case 0:
                            {
                                for (i = 0; i < order; i++)
                                {
                                    f[i] = 1.0;
                                }

                                break;
                            }
                            default:
                            {
                                for (i = 0; i < order; i++)
                                {
                                    f[i] = Math.Pow(x[i], n);
                                }

                                break;
                            }
                        }

                        estimate = 0.0;
                        for (i = 0; i < order; i++)
                        {
                            estimate += w[i] * f[i];
                        }

                        error = typeMethods.r8_abs(exact - estimate);

                        Console.WriteLine("  " + order.ToString().PadLeft(8)
                                               + "  " + n.ToString().PadLeft(8)
                                               + "  " + alpha.ToString("0.######").PadLeft(14)
                                               + "  " + estimate.ToString("0.######").PadLeft(14)
                                               + "  " + exact.ToString("0.######").PadLeft(14)
                                               + "  " + error.ToString("0.######").PadLeft(14) + "");
                    }
                }
            }
        }

        private static void test14()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST14 tests GEN_LAGUERRE_COMPUTE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int TEST_NUM = 3;

            double alpha;
            double[] alpha_test = {0.5, 1.0, 2.5};
            int i;
            int order;
            int order_max = 10;
            int test;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST14");
            Console.WriteLine("  GEN_LAGUERRE_COMPUTE computes a generalized Gauss-Laguerre rule");
            Console.WriteLine("  which is appropriate for integrands of the form");
            Console.WriteLine("    Integral ( 0 <= x < +oo ) f(x) x^alpha exp(-x) dx.");

            for (test = 0; test < TEST_NUM; test++)
            {
                alpha = alpha_test[test];

                for (order = 1; order <= order_max; order++)
                {
                    Console.WriteLine("");
                    Console.WriteLine("  Order = " + order + "");
                    Console.WriteLine("  ALPHA = " + alpha + "");

                    w = new double[order];
                    x = new double[order];

                    Burkardt.Laguerre.QuadratureRule.gen_laguerre_compute(order, alpha, ref x, ref w);

                    for (i = 0; i < order; i += 1)
                    {
                        Console.WriteLine("  " + i.ToString().PadLeft(8)
                                               + "  " + x[i].ToString("0.################").PadLeft(24)
                                               + "  " + w[i].ToString("0.################").PadLeft(24) + "");
                    }
                }
            }
        }

        private static void test15()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST15 tests HERMITE_COMPUTE against HERMITE_INTEGRAL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    19 April 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int d;
            double error;
            double estimate;
            double exact;
            double[] f;
            int i;
            int n;
            const int N_MAX = 10;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST15");
            Console.WriteLine("  HERMITE_COMPUTE computes a Gauss-Hermite rule");
            Console.WriteLine("  which is appropriate for integrands of the form");
            Console.WriteLine("    Integral ( -oo < x < +oo ) f(x) exp(-x*x) dx.");
            Console.WriteLine("");
            Console.WriteLine("  HERMITE_INTEGRAL determines the exact value of");
            Console.WriteLine("  this integal when f(x) = x^d.");
            Console.WriteLine("");
            Console.WriteLine("  A rule of order N should be exact for monomials x^d");
            Console.WriteLine("  up to degree D = 2*N-1");
            Console.WriteLine("");
            Console.WriteLine("  In the following table, for each order, the LAST THREE estimates");
            Console.WriteLine("  are made on monomials that exceed the exactness limit for the rule.");
            Console.WriteLine("");
            Console.WriteLine("         N         D       Estimate       Exact            Error");

            for (n = 1; n <= n_max; n++)
            {
                Console.WriteLine("");

                f = new double[n];
                w = new double[n];
                x = new double[n];

                HermiteQuadrature.hermite_compute(n, ref x, ref w);

                for (d = 0; d <= 2 * n + 2; d += 1)
                {
                    exact = Integral.hermite_integral(d);

                    switch (d)
                    {
                        case 0:
                        {
                            for (i = 0; i < n; i++)
                            {
                                f[i] = 1.0;
                            }

                            break;
                        }
                        default:
                        {
                            for (i = 0; i < n; i++)
                            {
                                f[i] = Math.Pow(x[i], d);
                            }

                            break;
                        }
                    }

                    estimate = 0.0;
                    for (i = 0; i < n; i++)
                    {
                        estimate += w[i] * f[i];
                    }

                    error = typeMethods.r8_abs(exact - estimate);

                    Console.WriteLine("  " + n.ToString().PadLeft(8)
                                           + "  " + d.ToString().PadLeft(8)
                                           + "  " + estimate.ToString("0.######").PadLeft(14)
                                           + "  " + exact.ToString("0.######").PadLeft(14)
                                           + "  " + error.ToString("0.######").PadLeft(14) + "");
                }
            }
        }

        private static void test16()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST16 tests HERMITE_COMPUTE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int i;
            int order;
            int order_max = 10;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST16");
            Console.WriteLine("  HERMITE_COMPUTE computes a Gauss-Hermite rule");
            Console.WriteLine("  which is appropriate for integrands of the form");
            Console.WriteLine("    Integral ( -oo < x < +oo ) f(x) exp(-x*x) dx.");

            for (order = 1; order <= order_max; order++)
            {
                Console.WriteLine("");
                Console.WriteLine("  Order = " + order + "");

                w = new double[order];
                x = new double[order];

                HermiteQuadrature.hermite_compute(order, ref x, ref w);

                for (i = 0; i < order; i += 1)
                {
                    Console.WriteLine("  " + i.ToString().PadLeft(8)
                                           + "  " + x[i].ToString("0.################").PadLeft(24)
                                           + "  " + w[i].ToString("0.################").PadLeft(24) + "");
                }
            }
        }

        private static void test17()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST17 tests JACOBI_COMPUTE against Integral.jacobi_intEGRAL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int TEST_NUM = 3;

            double alpha;
            double[] alpha_test = {0.5, 1.0, 2.5};
            double beta;
            double[] beta_test = {0.5, 1.0, 2.5};
            double error;
            double estimate;
            double exact;
            double[] f;
            int i;
            int n;
            int order;
            int order_max = 10;
            int test1;
            int test2;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST17");
            Console.WriteLine("  JACOBI_COMPUTE computes a Gauss-Jacobi rule");
            Console.WriteLine("  which is appropriate for integrands of the form");
            Console.WriteLine("    Integral ( -1 <= x <= +1 ) f(x) (1-x)^alpha (1+x)^beta dx.");
            Console.WriteLine("");
            Console.WriteLine("  Integral.jacobi_intEGRAL determines the exact value of");
            Console.WriteLine("  this integal when f(x) = x^n.");
            Console.WriteLine("");
            Console.WriteLine("  A rule of order ORDER should be exact for monomials X^N");
            Console.WriteLine("  up to N = 2*ORDER-1");
            Console.WriteLine("");
            Console.WriteLine("  In the following table, for each order, the LAST THREE estimates");
            Console.WriteLine("  are made on monomials that exceed the exactness limit for the rule.");
            Console.WriteLine("");
            Console.WriteLine(
                "     Order         N       Alpha           Beta            Estimate       Exact            Error");

            for (test1 = 0; test1 < TEST_NUM; test1++)
            {
                alpha = alpha_test[test1];

                for (test2 = 0; test2 < TEST_NUM; test2++)
                {
                    beta = beta_test[test2];

                    for (order = 1; order <= order_max; order++)
                    {
                        Console.WriteLine("");

                        f = new double[order];
                        w = new double[order];
                        x = new double[order];

                        JacobiQuadrature.jacobi_compute(order, alpha, beta, ref x, ref w);

                        for (n = 0; n <= 2 * order + 2; n += 1)
                        {
                            exact = Integral.jacobi_integral(n, alpha, beta);

                            switch (n)
                            {
                                case 0:
                                {
                                    for (i = 0; i < order; i++)
                                    {
                                        f[i] = 1.0;
                                    }

                                    break;
                                }
                                default:
                                {
                                    for (i = 0; i < order; i++)
                                    {
                                        f[i] = Math.Pow(x[i], n);
                                    }

                                    break;
                                }
                            }

                            estimate = 0.0;
                            for (i = 0; i < order; i++)
                            {
                                estimate += w[i] * f[i];
                            }

                            error = typeMethods.r8_abs(exact - estimate);

                            Console.WriteLine("  " + order.ToString().PadLeft(8)
                                                   + "  " + n.ToString().PadLeft(8)
                                                   + "  " + alpha.ToString("0.######").PadLeft(14)
                                                   + "  " + beta.ToString("0.######").PadLeft(14)
                                                   + "  " + estimate.ToString("0.######").PadLeft(14)
                                                   + "  " + exact.ToString("0.######").PadLeft(14)
                                                   + "  " + error.ToString("0.######").PadLeft(14) + "");
                        }
                    }
                }
            }
        }

        private static void test18()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST18 tests JACOBI_COMPUTE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int TEST_NUM = 3;

            double alpha;
            double[] alpha_test = {0.5, 1.0, 2.5};
            double beta;
            double[] beta_test = {0.5, 1.0, 2.5};
            int i;
            int n;
            const int N_MAX = 10;
            int test1;
            int test2;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST18");
            Console.WriteLine("  JACOBI_COMPUTE computes a Gauss-Jacobi rule");
            Console.WriteLine("  which is appropriate for integrands of the form");
            Console.WriteLine("    Integral ( -1 <= x <= +1 ) f(x) (1-x)^alpha (1+x)^beta dx.");

            for (test1 = 0; test1 < TEST_NUM; test1++)
            {
                alpha = alpha_test[test1];

                for (test2 = 0; test2 < TEST_NUM; test2++)
                {
                    beta = beta_test[test2];

                    for (n = 1; n <= n_max; n++)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("  N = " + n + "");
                        Console.WriteLine("  ALPHA = " + alpha + "");
                        Console.WriteLine("  BETA = " + beta + "");

                        w = new double[n];
                        x = new double[n];

                        JacobiQuadrature.jacobi_compute(n, alpha, beta, ref x, ref w);

                        for (i = 0; i < n; i++)
                        {
                            Console.WriteLine("  " + i.ToString().PadLeft(8)
                                                   + "  " + x[i].ToString("0.################").PadLeft(24)
                                                   + "  " + w[i].ToString("0.################").PadLeft(24) + "");
                        }
                    }
                }
            }
        }

        private static void test19()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST19 tests LAGUERRE_COMPUTE against Integral.laguerre_integral.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double error;
            double estimate;
            double exact;
            double[] f;
            int i;
            int n;
            int order;
            int order_max = 10;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST19");
            Console.WriteLine("  LAGUERRE_COMPUTE computes a Gauss-Laguerre rule");
            Console.WriteLine("  which is appropriate for integrands of the form");
            Console.WriteLine("    Integral ( 0 <= x < +oo ) f(x) exp(-x) dx.");
            Console.WriteLine("");
            Console.WriteLine("  Integral.laguerre_integral determines the exact value of");
            Console.WriteLine("  this integal when f(x) = x^n.");
            Console.WriteLine("");
            Console.WriteLine("  A rule of order ORDER should be exact for monomials X^N");
            Console.WriteLine("  up to N = 2*ORDER-1");
            Console.WriteLine("");
            Console.WriteLine("  In the following table, for each order, the LAST THREE estimates");
            Console.WriteLine("  are made on monomials that exceed the exactness limit for the rule.");
            Console.WriteLine("");
            Console.WriteLine("     Order         N       Estimate       Exact            Error");

            for (order = 1; order <= order_max; order++)
            {
                Console.WriteLine("");

                f = new double[order];
                w = new double[order];
                x = new double[order];

                Burkardt.Laguerre.QuadratureRule.laguerre_compute(order, ref x, ref w);

                for (n = 0; n <= 2 * order + 2; n += 1)
                {
                    exact = Integral.laguerre_integral(n);

                    switch (n)
                    {
                        case 0:
                        {
                            for (i = 0; i < order; i++)
                            {
                                f[i] = 1.0;
                            }

                            break;
                        }
                        default:
                        {
                            for (i = 0; i < order; i++)
                            {
                                f[i] = Math.Pow(x[i], n);
                            }

                            break;
                        }
                    }

                    estimate = 0.0;
                    for (i = 0; i < order; i++)
                    {
                        estimate += w[i] * f[i];
                    }

                    error = typeMethods.r8_abs(exact - estimate);

                    Console.WriteLine("  " + order.ToString().PadLeft(8)
                                           + "  " + n.ToString().PadLeft(8)
                                           + "  " + estimate.ToString("0.######").PadLeft(14)
                                           + "  " + exact.ToString("0.######").PadLeft(14)
                                           + "  " + error.ToString("0.######").PadLeft(14) + "");
                }
            }
        }

        private static void test20()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST20 tests LAGUERRE_COMPUTE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int i;
            int order;
            int order_max = 10;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST20");
            Console.WriteLine("  LAGUERRE_COMPUTE computes a generalized Gauss-Laguerre rule");
            Console.WriteLine("  which is appropriate for integrands of the form");
            Console.WriteLine("    Integral ( 0 <= x < +oo ) f(x) exp(-x) dx.");

            for (order = 1; order <= order_max; order++)
            {
                Console.WriteLine("");
                Console.WriteLine("  Order = " + order + "");

                w = new double[order];
                x = new double[order];

                Burkardt.Laguerre.QuadratureRule.laguerre_compute(order, ref x, ref w);

                for (i = 0; i < order; i += 1)
                {
                    Console.WriteLine("  " + i.ToString().PadLeft(8)
                                           + "  " + x[i].ToString("0.################").PadLeft(24)
                                           + "  " + w[i].ToString("0.################").PadLeft(24) + "");
                }
            }
        }

        private static void test21()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST21 tests LEGENDRE_COMPUTE against LEGENDRE_INTEGRAL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double error;
            double estimate;
            double exact;
            double[] f;
            int i;
            int n;
            int order;
            int order_max = 10;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST21");
            Console.WriteLine("  LEGENDRE_COMPUTE computes a Gauss-Legendre rule");
            Console.WriteLine("  which is appropriate for integrands of the form");
            Console.WriteLine("    Integral ( -1 <= x <= +1 ) f(x)  dx.");
            Console.WriteLine("");
            Console.WriteLine("  LEGENDRE_INTEGRAL determines the exact value of");
            Console.WriteLine("  this integal when f(x) = x^n.");
            Console.WriteLine("");
            Console.WriteLine("  A rule of order ORDER should be exact for monomials X^N");
            Console.WriteLine("  up to N = 2*ORDER-1");
            Console.WriteLine("");
            Console.WriteLine("  In the following table, for each order, the LAST THREE estimates");
            Console.WriteLine("  are made on monomials that exceed the exactness limit for the rule.");
            Console.WriteLine("");
            Console.WriteLine("     Order         N       Estimate       Exact            Error");

            for (order = 1; order <= order_max; order++)
            {
                Console.WriteLine("");

                f = new double[order];
                w = new double[order];
                x = new double[order];

                Burkardt.Legendre.QuadratureRule.legendre_compute(order, ref x, ref w);

                for (n = 0; n <= 2 * order + 2; n += 1)
                {
                    exact = Integral.legendre_integral(n);

                    switch (n)
                    {
                        case 0:
                        {
                            for (i = 0; i < order; i++)
                            {
                                f[i] = 1.0;
                            }

                            break;
                        }
                        default:
                        {
                            for (i = 0; i < order; i++)
                            {
                                f[i] = Math.Pow(x[i], n);
                            }

                            break;
                        }
                    }

                    estimate = 0.0;
                    for (i = 0; i < order; i++)
                    {
                        estimate += w[i] * f[i];
                    }

                    error = typeMethods.r8_abs(exact - estimate);

                    Console.WriteLine("  " + order.ToString().PadLeft(8)
                                           + "  " + n.ToString().PadLeft(8)
                                           + "  " + estimate.ToString("0.######").PadLeft(14)
                                           + "  " + exact.ToString("0.######").PadLeft(14)
                                           + "  " + error.ToString("0.######").PadLeft(14) + "");
                }
            }
        }

        private static void test22()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST22 tests LEGENDRE_COMPUTE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int i;
            int order;
            int order_max = 10;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST22");
            Console.WriteLine("  LEGENDRE_COMPUTE computes a Gauss-Legendre rule");
            Console.WriteLine("  which is appropriate for integrands of the form");
            Console.WriteLine("    Integral ( -1 <= x <= +1 ) f(x) dx.");

            for (order = 1; order <= order_max; order++)
            {
                Console.WriteLine("");
                Console.WriteLine("  Order = " + order + "");

                w = new double[order];
                x = new double[order];

                Burkardt.Legendre.QuadratureRule.legendre_compute(order, ref x, ref w);

                for (i = 0; i < order; i++)
                {
                    Console.WriteLine("  " + i.ToString().PadLeft(8)
                                           + "  " + x[i].ToString("0.################").PadLeft(24)
                                           + "  " + w[i].ToString("0.################").PadLeft(24) + "");
                }
            }
        }

        private static void test01_np()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST01_NP tests CHEBYSHEV1_COMPUTE_NP against CHEBYSHEV1_INTEGRAL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double error;
            double estimate;
            double exact;
            double[] f;
            int i;
            int n;
            int np = 0;
            int order;
            int order_max = 10;
            double[] p;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST01");
            Console.WriteLine("  CHEBYSHEV1_COMPUTE_NP computes a Gauss-Chebyshev type 1 rule");
            Console.WriteLine("  which is appropriate for integrands of the form");
            Console.WriteLine("    Integral ( -1 <= x <= +1 ) f(x) / Math.Sqrt ( 1 - x^2 ) dx.");
            Console.WriteLine("");
            Console.WriteLine("  CHEBYSHEV1_INTEGRAL determines the exact value of");
            Console.WriteLine("  this integal when f(x) = x^n.");
            Console.WriteLine("");
            Console.WriteLine("  A rule of order ORDER should be exact for monomials X^N");
            Console.WriteLine("  up to N = 2*ORDER-1");
            Console.WriteLine("");
            Console.WriteLine("  In the following table, for each order, the LAST THREE estimates");
            Console.WriteLine("  are made on monomials that exceed the exactness limit for the rule.");
            Console.WriteLine("");
            Console.WriteLine("     Order         N       Estimate       Exact            Error");

            p = new double[np];

            for (order = 1; order <= order_max; order++)
            {
                Console.WriteLine("");

                f = new double[order];
                w = new double[order];
                x = new double[order];

                Chebyshev1.chebyshev1_compute_np(order, np, p, ref x, ref w);

                for (n = 0; n <= 2 * order + 2; n += 1)
                {
                    exact = Integral.chebyshev1_integral(n);

                    switch (n)
                    {
                        case 0:
                        {
                            for (i = 0; i < order; i++)
                            {
                                f[i] = 1.0;
                            }

                            break;
                        }
                        default:
                        {
                            for (i = 0; i < order; i++)
                            {
                                f[i] = Math.Pow(x[i], n);
                            }

                            break;
                        }
                    }

                    estimate = 0.0;
                    for (i = 0; i < order; i++)
                    {
                        estimate += w[i] * f[i];
                    }

                    error = typeMethods.r8_abs(exact - estimate);

                    Console.WriteLine("  " + order.ToString().PadLeft(8)
                                           + "  " + n.ToString().PadLeft(8)
                                           + "  " + estimate.ToString("0.######").PadLeft(14)
                                           + "  " + exact.ToString("0.######").PadLeft(14)
                                           + "  " + error.ToString("0.######").PadLeft(14) + "");
                }
            }
        }

        private static void test02_np()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST02_NP tests CHEBYSHEV1_COMPUTE_NP.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int i;
            int np = 0;
            int order;
            int order_max = 10;
            double[] p;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST02");
            Console.WriteLine("  CHEBYSHEV1_COMPUTE_NP computes a Gauss-Chebyshev type 1 rule");
            Console.WriteLine("  which is appropriate for integrands of the form");
            Console.WriteLine("    Integral ( -1 <= x <= +1 ) f(x) / sqrt(1-x^2) dx.");

            p = new double[np];

            for (order = 1; order <= order_max; order++)
            {
                Console.WriteLine("");
                Console.WriteLine("  Order = " + order + "");

                w = new double[order];
                x = new double[order];

                Chebyshev1.chebyshev1_compute_np(order, np, p, ref x, ref w);

                for (i = 0; i < order; i += 1)
                {
                    Console.WriteLine("  " + i.ToString().PadLeft(8)
                                           + "  " + x[i].ToString("0.################").PadLeft(24)
                                           + "  " + w[i].ToString("0.################").PadLeft(24) + "");
                }
            }
        }

        private static void test03_np()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST03_NP tests CHEBYSHEV2_COMPUTE_NP against CHEBYSHEV2_INTEGRAL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double error;
            double estimate;
            double exact;
            double[] f;
            int i;
            int n;
            int np = 0;
            int order;
            int order_max = 10;
            double[] p;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST03");
            Console.WriteLine("  CHEBYSHEV2_COMPUTE_NP computes a Gauss-Chebyshev type 2 rule");
            Console.WriteLine("  which is appropriate for integrands of the form");
            Console.WriteLine("    Integral ( -1 <= x <= +1 ) f(x) * Math.Sqrt ( 1 - x^2 ) dx.");
            Console.WriteLine("");
            Console.WriteLine("  CHEBYSHEV2_INTEGRAL determines the exact value of");
            Console.WriteLine("  this integal when f(x) = x^n.");
            Console.WriteLine("");
            Console.WriteLine("  A rule of order ORDER should be exact for monomials X^N");
            Console.WriteLine("  up to N = 2*ORDER-1");
            Console.WriteLine("");
            Console.WriteLine("  In the following table, for each order, the LAST THREE estimates");
            Console.WriteLine("  are made on monomials that exceed the exactness limit for the rule.");
            Console.WriteLine("");
            Console.WriteLine("     Order         N       Estimate       Exact            Error");

            p = new double[np];

            for (order = 1; order <= order_max; order++)
            {
                Console.WriteLine("");

                f = new double[order];
                w = new double[order];
                x = new double[order];

                Chebyshev2.chebyshev2_compute_np(order, np, p, ref x, ref w);

                for (n = 0; n <= 2 * order + 2; n += 1)
                {
                    exact = Integral.chebyshev2_integral(n);

                    switch (n)
                    {
                        case 0:
                        {
                            for (i = 0; i < order; i++)
                            {
                                f[i] = 1.0;
                            }

                            break;
                        }
                        default:
                        {
                            for (i = 0; i < order; i++)
                            {
                                f[i] = Math.Pow(x[i], n);
                            }

                            break;
                        }
                    }

                    estimate = 0.0;
                    for (i = 0; i < order; i++)
                    {
                        estimate += w[i] * f[i];
                    }

                    error = typeMethods.r8_abs(exact - estimate);

                    Console.WriteLine("  " + order.ToString().PadLeft(8)
                                           + "  " + n.ToString().PadLeft(8)
                                           + "  " + estimate.ToString("0.######").PadLeft(14)
                                           + "  " + exact.ToString("0.######").PadLeft(14)
                                           + "  " + error.ToString("0.######").PadLeft(14) + "");
                }
            }

        }

        private static void test04_np()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST04_NP tests CHEBYSHEV2_COMPUTE_NP.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int i;
            int np = 0;
            int order;
            int order_max = 10;
            double[] p;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST04");
            Console.WriteLine("  CHEBYSHEV2_COMPUTE_NP computes a Gauss-Chebyshev type 2 rule");
            Console.WriteLine("  which is appropriate for integrands of the form");
            Console.WriteLine("    Integral ( -1 <= x <= +1 ) f(x) * sqrt(1-x^2) dx.");

            p = new double[np];

            for (order = 1; order <= order_max; order++)
            {
                Console.WriteLine("");
                Console.WriteLine("  Order = " + order + "");

                w = new double[order];
                x = new double[order];

                Chebyshev2.chebyshev2_compute_np(order, np, p, ref x, ref w);

                for (i = 0; i < order; i += 1)
                {
                    Console.WriteLine("  " + i.ToString().PadLeft(8)
                                           + "  " + x[i].ToString("0.################").PadLeft(24)
                                           + "  " + w[i].ToString("0.################").PadLeft(24) + "");
                }
            }

        }

        private static void test05_np()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST05_NP tests CLENSHAW_CURTIS_COMPUTE_NP against LEGENDRE_INTEGRAL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double error;
            double estimate;
            double exact;
            double[] f;
            int i;
            int n;
            int n_hi;
            int np = 0;
            int order;
            int order_max = 10;
            double[] p;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST05");
            Console.WriteLine("  CLENSHAW_CURTIS_COMPUTE_NP computes a Clenshaw Curtis rule");
            Console.WriteLine("  which is appropriate for integrands of the form");
            Console.WriteLine("    Integral ( -1 <= x <= +1 ) f(x)  dx.");
            Console.WriteLine("");
            Console.WriteLine("  LEGENDRE_INTEGRAL determines the exact value of");
            Console.WriteLine("  this integal when f(x) = x^n.");
            Console.WriteLine("");
            Console.WriteLine("  A rule of order ORDER should be exact for monomials X^N up to");
            Console.WriteLine("    N = ORDER+1 if ORDER is odd, or");
            Console.WriteLine("    N = ORDER   if ORDER is even");
            Console.WriteLine("");
            Console.WriteLine("  In the following table, for each order, the LAST THREE estimates");
            Console.WriteLine("  are made on monomials that exceed the exactness limit for the rule.");
            Console.WriteLine("");
            Console.WriteLine("     Order         N       Estimate       Exact            Error");

            p = new double[np];

            for (order = 1; order <= order_max; order++)
            {
                Console.WriteLine("");

                f = new double[order];
                w = new double[order];
                x = new double[order];

                ClenshawCurtis.clenshaw_curtis_compute_np(order, np, p, ref x, ref w);

                n_hi = (order % 2) switch
                {
                    0 => order + 2,
                    _ => order + 3
                };

                for (n = 0; n <= n_hi; n += 1)
                {
                    exact = Integral.legendre_integral(n);

                    switch (n)
                    {
                        case 0:
                        {
                            for (i = 0; i < order; i++)
                            {
                                f[i] = 1.0;
                            }

                            break;
                        }
                        default:
                        {
                            for (i = 0; i < order; i++)
                            {
                                f[i] = Math.Pow(x[i], n);
                            }

                            break;
                        }
                    }

                    estimate = 0.0;
                    for (i = 0; i < order; i++)
                    {
                        estimate += w[i] * f[i];
                    }

                    error = typeMethods.r8_abs(exact - estimate);

                    Console.WriteLine("  " + order.ToString().PadLeft(8)
                                           + "  " + n.ToString().PadLeft(8)
                                           + "  " + estimate.ToString("0.######").PadLeft(14)
                                           + "  " + exact.ToString("0.######").PadLeft(14)
                                           + "  " + error.ToString("0.######").PadLeft(14) + "");
                }
            }

        }

        private static void test06_np()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST06_NP tests CLENSHAW_CURTIS_COMPUTE_NP.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int i;
            int np = 0;
            int order;
            int order_max = 10;
            double[] p;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST06");
            Console.WriteLine("  CLENSHAW_CURTIS_COMPUTE_NP computes a Clenshaw Curtis rule");
            Console.WriteLine("  which is appropriate for integrands of the form");
            Console.WriteLine("    Integral ( -1 <= x <= +1 ) f(x) dx.");

            p = new double[np];

            for (order = 1; order <= order_max; order++)
            {
                Console.WriteLine("");
                Console.WriteLine("  Order = " + order + "");

                w = new double[order];
                x = new double[order];

                ClenshawCurtis.clenshaw_curtis_compute_np(order, np, p, ref x, ref w);

                for (i = 0; i < order; i += 1)
                {
                    Console.WriteLine("  " + i.ToString().PadLeft(8)
                                           + "  " + x[i].ToString("0.################").PadLeft(24)
                                           + "  " + w[i].ToString("0.################").PadLeft(24) + "");
                }
            }

        }

        private static void test07_np()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST07_NP tests FEJER2_COMPUTE_NP against LEGENDRE_INTEGRAL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double error;
            double estimate;
            double exact;
            double[] f;
            int i;
            int n;
            int n_hi;
            int np = 0;
            int order;
            int order_max = 10;
            double[] p;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST07");
            Console.WriteLine("  FEJER2_COMPUTE_NP computes a Fejer Type 2 rule");
            Console.WriteLine("  which is appropriate for integrands of the form");
            Console.WriteLine("    Integral ( -1 <= x <= +1 ) f(x)  dx.");
            Console.WriteLine("");
            Console.WriteLine("  LEGENDRE_INTEGRAL determines the exact value of");
            Console.WriteLine("  this integal when f(x) = x^n.");
            Console.WriteLine("");
            Console.WriteLine("  A rule of order ORDER should be exact for monomials X^N up to");
            Console.WriteLine("    N = ORDER+1 if ORDER is odd, or");
            Console.WriteLine("    N = ORDER   if ORDER is even");
            Console.WriteLine("");
            Console.WriteLine("  In the following table, for each order, the LAST THREE estimates");
            Console.WriteLine("  are made on monomials that exceed the exactness limit for the rule.");
            Console.WriteLine("");
            Console.WriteLine("     Order         N       Estimate       Exact            Error");

            p = new double[np];

            for (order = 1; order <= order_max; order++)
            {
                Console.WriteLine("");

                f = new double[order];
                w = new double[order];
                x = new double[order];

                Fejer2.fejer2_compute_np(order, np, p, ref x, ref w);

                n_hi = (order % 2) switch
                {
                    0 => order + 2,
                    _ => order + 3
                };

                for (n = 0; n <= n_hi; n += 1)
                {
                    exact = Integral.legendre_integral(n);

                    switch (n)
                    {
                        case 0:
                        {
                            for (i = 0; i < order; i++)
                            {
                                f[i] = 1.0;
                            }

                            break;
                        }
                        default:
                        {
                            for (i = 0; i < order; i++)
                            {
                                f[i] = Math.Pow(x[i], n);
                            }

                            break;
                        }
                    }

                    estimate = 0.0;
                    for (i = 0; i < order; i++)
                    {
                        estimate += w[i] * f[i];
                    }

                    error = typeMethods.r8_abs(exact - estimate);

                    Console.WriteLine("  " + order.ToString().PadLeft(8)
                                           + "  " + n.ToString().PadLeft(8)
                                           + "  " + estimate.ToString("0.######").PadLeft(14)
                                           + "  " + exact.ToString("0.######").PadLeft(14)
                                           + "  " + error.ToString("0.######").PadLeft(14) + "");
                }
            }

        }

        private static void test08_np()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST08_NP tests FEJER2_COMPUTE_NP.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int i;
            int np = 0;
            int order;
            int order_max = 10;
            double[] p;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST08");
            Console.WriteLine("  FEJER2_COMPUTE_NP computes a Fejer Type 2 rule");
            Console.WriteLine("  which is appropriate for integrands of the form");
            Console.WriteLine("    Integral ( -1 <= x <= +1 ) f(x) dx.");

            p = new double[np];

            for (order = 1; order <= order_max; order++)
            {
                Console.WriteLine("");
                Console.WriteLine("  Order = " + order + "");

                w = new double[order];
                x = new double[order];

                Fejer2.fejer2_compute_np(order, np, p, ref x, ref w);

                for (i = 0; i < order; i += 1)
                {
                    Console.WriteLine("  " + i.ToString().PadLeft(8)
                                           + "  " + x[i].ToString("0.################").PadLeft(24)
                                           + "  " + w[i].ToString("0.################").PadLeft(24) + "");
                }
            }

        }

        private static void test09_np()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST09_NP tests GEGENBAUER_COMPUTE_NP against GEGENBAUER_INTEGRAL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int TEST_NUM = 3;

            double alpha;
            double[] alpha_test = {0.5, 1.0, 2.5};
            double error;
            double estimate;
            double exact;
            double[] f;
            int i;
            int n;
            int np = 1;
            int order;
            int order_max = 10;
            double[] p;
            int test;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST09");
            Console.WriteLine("  GEGENBAUER_COMPUTE_NP computes a generalized Gauss-Gegenbauer rule");
            Console.WriteLine("  which is appropriate for integrands of the form");
            Console.WriteLine("    Integral ( 0 <= x < +oo ) f(x) (1-x^2)^alpha dx.");
            Console.WriteLine("");
            Console.WriteLine("  GEGENBAUER_INTEGRAL determines the exact value of");
            Console.WriteLine("  this integal when f(x) = x^n.");
            Console.WriteLine("");
            Console.WriteLine("  A rule of order ORDER should be exact for monomials X^N");
            Console.WriteLine("  up to N = 2*ORDER-1");
            Console.WriteLine("");
            Console.WriteLine("  In the following table, for each order, the LAST THREE estimates");
            Console.WriteLine("  are made on monomials that exceed the exactness limit for the rule.");
            Console.WriteLine("");
            Console.WriteLine("     Order         N       Alpha           Estimate       Exact            Error");

            p = new double[np];

            for (test = 0; test < TEST_NUM; test++)
            {
                alpha = alpha_test[test];
                p[0] = alpha;

                for (order = 1; order <= order_max; order++)
                {
                    Console.WriteLine("");

                    f = new double[order];
                    w = new double[order];
                    x = new double[order];

                    GegenbauerQuadrature.gegenbauer_compute_np(order, np, p, ref x, ref w);

                    for (n = 0; n <= 2 * order + 2; n += 1)
                    {
                        exact = Integral.gegenbauer_integral(n, alpha);

                        switch (n)
                        {
                            case 0:
                            {
                                for (i = 0; i < order; i++)
                                {
                                    f[i] = 1.0;
                                }

                                break;
                            }
                            default:
                            {
                                for (i = 0; i < order; i++)
                                {
                                    f[i] = Math.Pow(x[i], n);
                                }

                                break;
                            }
                        }

                        estimate = 0.0;
                        for (i = 0; i < order; i++)
                        {
                            estimate += w[i] * f[i];
                        }

                        error = typeMethods.r8_abs(exact - estimate);

                        Console.WriteLine("  " + order.ToString().PadLeft(8)
                                               + "  " + n.ToString().PadLeft(8)
                                               + "  " + alpha.ToString().PadLeft(14)
                                               + "  " + estimate.ToString("0.######").PadLeft(14)
                                               + "  " + exact.ToString("0.######").PadLeft(14)
                                               + "  " + error.ToString("0.######").PadLeft(14) + "");
                    }
                }
            }
        }

        private static void test10_np()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST10_NP tests GEGENBAUER_COMPUTE_NP.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int TEST_NUM = 3;

            double alpha;
            double[] alpha_test = {0.5, 1.0, 2.5};
            int i;
            int np = 0;
            int order;
            int order_max = 10;
            double[] p;
            int test;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST10");
            Console.WriteLine("  GEGENBAUER_COMPUTE_NP computes a generalized Gauss-Gegenbauer rule");
            Console.WriteLine("  which is appropriate for integrands of the form");
            Console.WriteLine("    Integral ( -1 <= x <= +1 ) f(x) (1-x^2)^alpha dx.");

            p = new double[np];

            for (test = 0; test < TEST_NUM; test++)
            {
                alpha = alpha_test[test];
                p[0] = alpha;

                for (order = 1; order <= order_max; order++)
                {
                    Console.WriteLine("");
                    Console.WriteLine("  Order = " + order + "");
                    Console.WriteLine("  ALPHA = " + alpha + "");

                    w = new double[order];
                    x = new double[order];

                    GegenbauerQuadrature.gegenbauer_compute_np(order, np, p, ref x, ref w);

                    for (i = 0; i < order; i += 1)
                    {
                        Console.WriteLine("  " + i.ToString().PadLeft(8)
                                               + "  " + x[i].ToString("0.################").PadLeft(24)
                                               + "  " + w[i].ToString("0.################").PadLeft(24) + "");
                    }
                }
            }
        }

        private static void test11_np()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST11_NP tests GEN_HERMITE_COMPUTE_NP against GEN_HERMITE_INTEGRAL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int TEST_NUM = 3;

            double alpha;
            double[] alpha_test = {0.5, 1.0, 2.5};
            double error;
            double estimate;
            double exact;
            double[] f;
            int i;
            int n;
            int np = 1;
            int order;
            int order_max = 10;
            double[] p;
            int test;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST11");
            Console.WriteLine("  GEN_HERMITE_COMPUTE_NP computes a generalized Gauss-Hermite rule");
            Console.WriteLine("  which is appropriate for integrands of the form");
            Console.WriteLine("    Integral ( -oo < x < +oo ) f(x) x^alpha exp(-x*x) dx.");
            Console.WriteLine("");
            Console.WriteLine("  GEN_HERMITE_INTEGRAL determines the exact value of");
            Console.WriteLine("  this integal when f(x) = x^n.");
            Console.WriteLine("");
            Console.WriteLine("  A rule of order ORDER should be exact for monomials X^N");
            Console.WriteLine("  up to N = 2*ORDER-1");
            Console.WriteLine("");
            Console.WriteLine("  In the following table, for each order, the LAST THREE estimates");
            Console.WriteLine("  are made on monomials that exceed the exactness limit for the rule.");
            Console.WriteLine("");
            Console.WriteLine("     Order         N       Alpha           Estimate       Exact            Error");

            p = new double[np];

            for (test = 0; test < TEST_NUM; test++)
            {
                alpha = alpha_test[test];
                p[0] = alpha;

                for (order = 1; order <= order_max; order++)
                {
                    Console.WriteLine("");

                    f = new double[order];
                    w = new double[order];
                    x = new double[order];

                    HermiteQuadrature.gen_hermite_compute_np(order, np, p, ref x, ref w);

                    for (n = 0; n <= 2 * order + 2; n += 1)
                    {
                        exact = Integral.gen_hermite_integral(n, alpha);

                        switch (n)
                        {
                            case 0:
                            {
                                for (i = 0; i < order; i++)
                                {
                                    f[i] = 1.0;
                                }

                                break;
                            }
                            default:
                            {
                                for (i = 0; i < order; i++)
                                {
                                    f[i] = Math.Pow(x[i], n);
                                }

                                break;
                            }
                        }

                        estimate = 0.0;
                        for (i = 0; i < order; i++)
                        {
                            estimate += w[i] * f[i];
                        }

                        error = typeMethods.r8_abs(exact - estimate);

                        Console.WriteLine("  " + order.ToString().PadLeft(8)
                                               + "  " + n.ToString().PadLeft(8)
                                               + "  " + alpha.ToString("0.######").PadLeft(14)
                                               + "  " + estimate.ToString("0.######").PadLeft(14)
                                               + "  " + exact.ToString("0.######").PadLeft(14)
                                               + "  " + error.ToString("0.######").PadLeft(14) + "");
                    }
                }
            }
        }

        private static void test12_np()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST12_NP tests GEN_HERMITE_COMPUTE_NP.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int TEST_NUM = 3;

            double alpha;
            double[] alpha_test = {0.5, 1.0, 2.5};
            int i;
            int np = 1;
            int order;
            int order_max = 10;
            double[] p;
            int test;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST12");
            Console.WriteLine("  GEN_HERMITE_COMPUTE_NP computes a generalized Gauss-Hermite rule");
            Console.WriteLine("  which is appropriate for integrands of the form");
            Console.WriteLine("    Integral ( -oo < x < +oo ) f(x) x^alpha exp(-x*x) dx.");

            p = new double[np];

            for (test = 0; test < TEST_NUM; test++)
            {
                alpha = alpha_test[test];
                p[0] = alpha;

                for (order = 1; order <= order_max; order++)
                {
                    Console.WriteLine("");
                    Console.WriteLine("  Order = " + order + "");
                    Console.WriteLine("  ALPHA = " + alpha + "");

                    w = new double[order];
                    x = new double[order];

                    HermiteQuadrature.gen_hermite_compute_np(order, np, p, ref x, ref w);

                    for (i = 0; i < order; i += 1)
                    {
                        Console.WriteLine("  " + i.ToString().PadLeft(8)
                                               + "  " + x[i].ToString("0.################").PadLeft(24)
                                               + "  " + w[i].ToString("0.################").PadLeft(24) + "");
                    }
                }
            }
        }

        private static void test13_np()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST13_NP tests GEN_LAGUERRE_COMPUTE_NP against GEN_LAGUERRE_INTEGRAL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int TEST_NUM = 3;

            double alpha;
            double[] alpha_test = {0.5, 1.0, 2.5};
            double error;
            double estimate;
            double exact;
            double[] f;
            int i;
            int n;
            int np = 1;
            int order;
            int order_max = 10;
            double[] p;
            int test;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST13");
            Console.WriteLine("  GEN_LAGUERRE_COMPUTE_NP computes a generalized Gauss-Laguerre rule");
            Console.WriteLine("  which is appropriate for integrands of the form");
            Console.WriteLine("    Integral ( 0 <= x < +oo ) f(x) x^alpha exp(-x) dx.");
            Console.WriteLine("");
            Console.WriteLine("  GEN_LAGUERRE_INTEGRAL determines the exact value of");
            Console.WriteLine("  this integal when f(x) = x^n.");
            Console.WriteLine("");
            Console.WriteLine("  A rule of order ORDER should be exact for monomials X^N");
            Console.WriteLine("  up to N = 2*ORDER-1");
            Console.WriteLine("");
            Console.WriteLine("  In the following table, for each order, the LAST THREE estimates");
            Console.WriteLine("  are made on monomials that exceed the exactness limit for the rule.");
            Console.WriteLine("");
            Console.WriteLine("     Order         N       Alpha           Estimate       Exact            Error");

            p = new double[np];

            for (test = 0; test < TEST_NUM; test++)
            {
                alpha = alpha_test[test];
                p[0] = alpha;

                for (order = 1; order <= order_max; order++)
                {
                    Console.WriteLine("");

                    f = new double[order];
                    w = new double[order];
                    x = new double[order];

                    Burkardt.Laguerre.QuadratureRule.gen_laguerre_compute_np(order, np, p, ref x, ref w);

                    for (n = 0; n <= 2 * order + 2; n += 1)
                    {
                        exact = Integral.gen_laguerre_integral(n, alpha);

                        switch (n)
                        {
                            case 0:
                            {
                                for (i = 0; i < order; i++)
                                {
                                    f[i] = 1.0;
                                }

                                break;
                            }
                            default:
                            {
                                for (i = 0; i < order; i++)
                                {
                                    f[i] = Math.Pow(x[i], n);
                                }

                                break;
                            }
                        }

                        estimate = 0.0;
                        for (i = 0; i < order; i++)
                        {
                            estimate += w[i] * f[i];
                        }

                        error = typeMethods.r8_abs(exact - estimate);

                        Console.WriteLine("  " + order.ToString().PadLeft(8)
                                               + "  " + n.ToString().PadLeft(8)
                                               + "  " + alpha.ToString("0.######").PadLeft(14)
                                               + "  " + estimate.ToString("0.######").PadLeft(14)
                                               + "  " + exact.ToString("0.######").PadLeft(14)
                                               + "  " + error.ToString("0.######").PadLeft(14) + "");
                    }
                }
            }
        }

        private static void test14_np()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST14_NP tests GEN_LAGUERRE_COMPUTE_NP.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int TEST_NUM = 3;

            double alpha;
            double[] alpha_test = {0.5, 1.0, 2.5};
            int i;
            int np = 1;
            int order;
            int order_max = 10;
            double[] p;
            int test;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST14");
            Console.WriteLine("  GEN_LAGUERRE_COMPUTE_NP computes a generalized Gauss-Laguerre rule");
            Console.WriteLine("  which is appropriate for integrands of the form");
            Console.WriteLine("    Integral ( 0 <= x < +oo ) f(x) x^alpha exp(-x) dx.");

            p = new double[np];

            for (test = 0; test < TEST_NUM; test++)
            {
                alpha = alpha_test[test];
                p[0] = alpha;

                for (order = 1; order <= order_max; order++)
                {
                    Console.WriteLine("");
                    Console.WriteLine("  Order = " + order + "");
                    Console.WriteLine("  ALPHA = " + alpha + "");

                    w = new double[order];
                    x = new double[order];

                    Burkardt.Laguerre.QuadratureRule.gen_laguerre_compute_np(order, np, p, ref x, ref w);

                    for (i = 0; i < order; i += 1)
                    {
                        Console.WriteLine("  " + i.ToString().PadLeft(8)
                                               + "  " + x[i].ToString("0.################").PadLeft(24)
                                               + "  " + w[i].ToString("0.################").PadLeft(24) + "");
                    }
                }
            }
        }

        private static void test15_np()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST15_NP tests HERMITE_COMPUTE_NP against HERMITE_INTEGRAL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double error;
            double estimate;
            double exact;
            double[] f;
            int i;
            int n;
            int np = 0;
            int order;
            int order_max = 10;
            double[] p;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST15");
            Console.WriteLine("  HERMITE_COMPUTE_NP computes a Gauss-Hermite rule");
            Console.WriteLine("  which is appropriate for integrands of the form");
            Console.WriteLine("    Integral ( -oo < x < +oo ) f(x) exp(-x*x) dx.");
            Console.WriteLine("");
            Console.WriteLine("  HERMITE_INTEGRAL determines the exact value of");
            Console.WriteLine("  this integal when f(x) = x^n.");
            Console.WriteLine("");
            Console.WriteLine("  A rule of order ORDER should be exact for monomials X^N");
            Console.WriteLine("  up to N = 2*ORDER-1");
            Console.WriteLine("");
            Console.WriteLine("  In the following table, for each order, the LAST THREE estimates");
            Console.WriteLine("  are made on monomials that exceed the exactness limit for the rule.");
            Console.WriteLine("");
            Console.WriteLine("     Order         N       Estimate       Exact            Error");

            p = new double[np];

            for (order = 1; order <= order_max; order++)
            {
                Console.WriteLine("");

                f = new double[order];
                w = new double[order];
                x = new double[order];

                HermiteQuadrature.hermite_compute_np(order, np, p, ref x, ref w);

                for (n = 0; n <= 2 * order + 2; n += 1)
                {
                    exact = Integral.hermite_integral(n);

                    switch (n)
                    {
                        case 0:
                        {
                            for (i = 0; i < order; i++)
                            {
                                f[i] = 1.0;
                            }

                            break;
                        }
                        default:
                        {
                            for (i = 0; i < order; i++)
                            {
                                f[i] = Math.Pow(x[i], n);
                            }

                            break;
                        }
                    }

                    estimate = 0.0;
                    for (i = 0; i < order; i++)
                    {
                        estimate += w[i] * f[i];
                    }

                    error = typeMethods.r8_abs(exact - estimate);

                    Console.WriteLine("  " + order.ToString().PadLeft(8)
                                           + "  " + n.ToString().PadLeft(8)
                                           + "  " + estimate.ToString("0.######").PadLeft(14)
                                           + "  " + exact.ToString("0.######").PadLeft(14)
                                           + "  " + error.ToString("0.######").PadLeft(14) + "");
                }
            }

        }

        private static void test16_np()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST16_NP tests HERMITE_COMPUTE_NP.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int i;
            int np = 0;
            int order;
            int order_max = 10;
            double[] p;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST16");
            Console.WriteLine("  HERMITE_COMPUTE_NP computes a Gauss-Hermite rule");
            Console.WriteLine("  which is appropriate for integrands of the form");
            Console.WriteLine("    Integral ( -oo < x < +oo ) f(x) exp(-x*x) dx.");

            p = new double[np];

            for (order = 1; order <= order_max; order++)
            {
                Console.WriteLine("");
                Console.WriteLine("  Order = " + order + "");

                w = new double[order];
                x = new double[order];

                HermiteQuadrature.hermite_compute_np(order, np, p, ref x, ref w);

                for (i = 0; i < order; i += 1)
                {
                    Console.WriteLine("  " + i.ToString().PadLeft(8)
                                           + "  " + x[i].ToString("0.################").PadLeft(24)
                                           + "  " + w[i].ToString("0.################").PadLeft(24) + "");
                }
            }

        }

        private static void test17_np()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST17_NP tests JACOBI_COMPUTE_NP against Integral.jacobi_intEGRAL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int TEST_NUM = 3;

            double alpha;
            double[] alpha_test = {0.5, 1.0, 2.5};
            double beta;
            double[] beta_test = {0.5, 1.0, 2.5};
            double error;
            double estimate;
            double exact;
            double[] f;
            int i;
            int n;
            int np = 2;
            int order;
            int order_max = 10;
            double[] p;
            int test1;
            int test2;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST17");
            Console.WriteLine("  JACOBI_COMPUTE_NP computes a Gauss-Jacobi rule");
            Console.WriteLine("  which is appropriate for integrands of the form");
            Console.WriteLine("    Integral ( -1 <= x <= +1 ) f(x) (1-x)^alpha (1+x)^beta dx.");
            Console.WriteLine("");
            Console.WriteLine("  Integral.jacobi_intEGRAL determines the exact value of");
            Console.WriteLine("  this integal when f(x) = x^n.");
            Console.WriteLine("");
            Console.WriteLine("  A rule of order ORDER should be exact for monomials X^N");
            Console.WriteLine("  up to N = 2*ORDER-1");
            Console.WriteLine("");
            Console.WriteLine("  In the following table, for each order, the LAST THREE estimates");
            Console.WriteLine("  are made on monomials that exceed the exactness limit for the rule.");
            Console.WriteLine("");
            Console.WriteLine(
                "     Order         N       Alpha           Beta            Estimate       Exact            Error");

            p = new double[np];

            for (test1 = 0; test1 < TEST_NUM; test1++)
            {
                alpha = alpha_test[test1];
                p[0] = alpha;

                for (test2 = 0; test2 < TEST_NUM; test2++)
                {
                    beta = beta_test[test2];
                    p[1] = beta;

                    for (order = 1; order <= order_max; order++)
                    {
                        Console.WriteLine("");

                        f = new double[order];
                        w = new double[order];
                        x = new double[order];

                        JacobiQuadrature.jacobi_compute_np(order, np, p, ref x, ref w);

                        for (n = 0; n <= 2 * order + 2; n += 1)
                        {
                            exact = Integral.jacobi_integral(n, alpha, beta);

                            switch (n)
                            {
                                case 0:
                                {
                                    for (i = 0; i < order; i++)
                                    {
                                        f[i] = 1.0;
                                    }

                                    break;
                                }
                                default:
                                {
                                    for (i = 0; i < order; i++)
                                    {
                                        f[i] = Math.Pow(x[i], n);
                                    }

                                    break;
                                }
                            }

                            estimate = 0.0;
                            for (i = 0; i < order; i++)
                            {
                                estimate += w[i] * f[i];
                            }

                            error = typeMethods.r8_abs(exact - estimate);

                            Console.WriteLine("  " + order.ToString().PadLeft(8)
                                                   + "  " + n.ToString().PadLeft(8)
                                                   + "  " + alpha.ToString("0.######").PadLeft(14)
                                                   + "  " + beta.ToString("0.######").PadLeft(14)
                                                   + "  " + estimate.ToString("0.######").PadLeft(14)
                                                   + "  " + exact.ToString("0.######").PadLeft(14)
                                                   + "  " + error.ToString("0.######").PadLeft(14) + "");
                        }
                    }
                }
            }
        }

        private static void test18_np()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST18_NP tests JACOBI_COMPUTE_NP.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int TEST_NUM = 3;

            double alpha;
            double[] alpha_test = {0.5, 1.0, 2.5};
            double beta;
            double[] beta_test = {0.5, 1.0, 2.5};
            int i;
            int np = 2;
            int order;
            int order_max = 10;
            double[] p;
            int test1;
            int test2;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST18");
            Console.WriteLine("  JACOBI_COMPUTE_NP computes a Gauss-Jacobi rule");
            Console.WriteLine("  which is appropriate for integrands of the form");
            Console.WriteLine("    Integral ( -1 <= x <= +1 ) f(x) (1-x)^alpha (1+x)^beta dx.");

            p = new double[np];

            for (test1 = 0; test1 < TEST_NUM; test1++)
            {
                alpha = alpha_test[test1];
                p[0] = alpha;

                for (test2 = 0; test2 < TEST_NUM; test2++)
                {
                    beta = beta_test[test2];
                    p[1] = beta;

                    for (order = 1; order <= order_max; order++)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("  Order = " + order + "");
                        Console.WriteLine("  ALPHA = " + alpha + "");
                        Console.WriteLine("  BETA = " + beta + "");

                        w = new double[order];
                        x = new double[order];

                        JacobiQuadrature.jacobi_compute_np(order, np, p, ref x, ref w);

                        for (i = 0; i < order; i += 1)
                        {
                            Console.WriteLine("  " + i.ToString().PadLeft(8)
                                                   + "  " + x[i].ToString("0.################").PadLeft(24)
                                                   + "  " + w[i].ToString("0.################").PadLeft(24) + "");
                        }
                    }
                }
            }
        }

        private static void test19_np()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST19_NP tests LAGUERRE_COMPUTE_NP against Integral.laguerre_integral.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double error;
            double estimate;
            double exact;
            double[] f;
            int i;
            int n;
            int np = 0;
            int order;
            int order_max = 10;
            double[] p;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST19");
            Console.WriteLine("  LAGUERRE_COMPUTE_NP computes a Gauss-Laguerre rule");
            Console.WriteLine("  which is appropriate for integrands of the form");
            Console.WriteLine("    Integral ( 0 <= x < +oo ) f(x) exp(-x) dx.");
            Console.WriteLine("");
            Console.WriteLine("  Integral.laguerre_integral determines the exact value of");
            Console.WriteLine("  this integal when f(x) = x^n.");
            Console.WriteLine("");
            Console.WriteLine("  A rule of order ORDER should be exact for monomials X^N");
            Console.WriteLine("  up to N = 2*ORDER-1");
            Console.WriteLine("");
            Console.WriteLine("  In the following table, for each order, the LAST THREE estimates");
            Console.WriteLine("  are made on monomials that exceed the exactness limit for the rule.");
            Console.WriteLine("");
            Console.WriteLine("     Order         N       Estimate       Exact            Error");

            p = new double[np];

            for (order = 1; order <= order_max; order++)
            {
                Console.WriteLine("");

                f = new double[order];
                w = new double[order];
                x = new double[order];

                Burkardt.Laguerre.QuadratureRule.laguerre_compute_np(order, np, p, ref x, ref w);

                for (n = 0; n <= 2 * order + 2; n += 1)
                {
                    exact = Integral.laguerre_integral(n);

                    switch (n)
                    {
                        case 0:
                        {
                            for (i = 0; i < order; i++)
                            {
                                f[i] = 1.0;
                            }

                            break;
                        }
                        default:
                        {
                            for (i = 0; i < order; i++)
                            {
                                f[i] = Math.Pow(x[i], n);
                            }

                            break;
                        }
                    }

                    estimate = 0.0;
                    for (i = 0; i < order; i++)
                    {
                        estimate += w[i] * f[i];
                    }

                    error = typeMethods.r8_abs(exact - estimate);

                    Console.WriteLine("  " + order.ToString().PadLeft(8)
                                           + "  " + n.ToString().PadLeft(8)
                                           + "  " + estimate.ToString("0.######").PadLeft(14)
                                           + "  " + exact.ToString("0.######").PadLeft(14)
                                           + "  " + error.ToString("0.######").PadLeft(14) + "");
                }
            }

        }

        private static void test20_np()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST20_NP tests LAGUERRE_COMPUTE_NP.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int i;
            int np = 0;
            int order;
            int order_max = 10;
            double[] p;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST20");
            Console.WriteLine("  LAGUERRE_COMPUTE_NP computes a generalized Gauss-Laguerre rule");
            Console.WriteLine("  which is appropriate for integrands of the form");
            Console.WriteLine("    Integral ( 0 <= x < +oo ) f(x) exp(-x) dx.");

            p = new double[np];

            for (order = 1; order <= order_max; order++)
            {
                Console.WriteLine("");
                Console.WriteLine("  Order = " + order + "");

                w = new double[order];
                x = new double[order];

                Burkardt.Laguerre.QuadratureRule.laguerre_compute_np(order, np, p, ref x, ref w);

                for (i = 0; i < order; i += 1)
                {
                    Console.WriteLine("  " + i.ToString().PadLeft(8)
                                           + "  " + x[i].ToString("0.################").PadLeft(24)
                                           + "  " + w[i].ToString("0.################").PadLeft(24) + "");
                }
            }

        }

        private static void test21_np()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST21_NP tests LEGENDRE_COMPUTE_NP against LEGENDRE_INTEGRAL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double error;
            double estimate;
            double exact;
            double[] f;
            int i;
            int n;
            int np = 0;
            int order;
            int order_max = 10;
            double[] p;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST21");
            Console.WriteLine("  LEGENDRE_COMPUTE_NP computes a Gauss-Legendre rule");
            Console.WriteLine("  which is appropriate for integrands of the form");
            Console.WriteLine("    Integral ( -1 <= x <= +1 ) f(x)  dx.");
            Console.WriteLine("");
            Console.WriteLine("  LEGENDRE_INTEGRAL determines the exact value of");
            Console.WriteLine("  this integal when f(x) = x^n.");
            Console.WriteLine("");
            Console.WriteLine("  A rule of order ORDER should be exact for monomials X^N");
            Console.WriteLine("  up to N = 2*ORDER-1");
            Console.WriteLine("");
            Console.WriteLine("  In the following table, for each order, the LAST THREE estimates");
            Console.WriteLine("  are made on monomials that exceed the exactness limit for the rule.");
            Console.WriteLine("");
            Console.WriteLine("     Order         N       Estimate       Exact            Error");

            p = new double[np];

            for (order = 1; order <= order_max; order++)
            {
                Console.WriteLine("");

                f = new double[order];
                w = new double[order];
                x = new double[order];

                Burkardt.Legendre.QuadratureRule.legendre_compute_np(order, np, p, ref x, ref w);

                for (n = 0; n <= 2 * order + 2; n += 1)
                {
                    exact = Integral.legendre_integral(n);

                    switch (n)
                    {
                        case 0:
                        {
                            for (i = 0; i < order; i++)
                            {
                                f[i] = 1.0;
                            }

                            break;
                        }
                        default:
                        {
                            for (i = 0; i < order; i++)
                            {
                                f[i] = Math.Pow(x[i], n);
                            }

                            break;
                        }
                    }

                    estimate = 0.0;
                    for (i = 0; i < order; i++)
                    {
                        estimate += w[i] * f[i];
                    }

                    error = typeMethods.r8_abs(exact - estimate);

                    Console.WriteLine("  " + order.ToString().PadLeft(8)
                                           + "  " + n.ToString().PadLeft(8)
                                           + "  " + estimate.ToString("0.######").PadLeft(14)
                                           + "  " + exact.ToString("0.######").PadLeft(14)
                                           + "  " + error.ToString("0.######").PadLeft(14) + "");
                }
            }

        }

        private static void test22_np()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST22_NP tests LEGENDRE_COMPUTE_NP.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int i;
            int np = 0;
            int order;
            int order_max = 10;
            double[] p;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST22");
            Console.WriteLine("  LEGENDRE_COMPUTE_NP computes a Gauss-Legendre rule");
            Console.WriteLine("  which is appropriate for integrands of the form");
            Console.WriteLine("    Integral ( -1 <= x <= +1 ) f(x) dx.");

            p = new double[np];

            for (order = 1; order <= order_max; order++)
            {
                Console.WriteLine("");
                Console.WriteLine("  Order = " + order + "");

                w = new double[order];
                x = new double[order];

                Burkardt.Legendre.QuadratureRule.legendre_compute_np(order, np, p, ref x, ref w);

                for (i = 0; i < order; i += 1)
                {
                    Console.WriteLine("  " + i.ToString().PadLeft(8)
                                           + "  " + x[i].ToString("0.################").PadLeft(24)
                                           + "  " + w[i].ToString("0.################").PadLeft(24) + "");
                }
            }

        }

        private static void test23(int r)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST23 tests LEVEL_GROWTH_TO_ORDER.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    16 October 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int R, the index of the rule to be examined.
            //
        {
            int dim;
            int dim_num = 11;
            int g;
            int[] growth;
            int[] level;
            int[] order;
            int[] rule;

            Console.WriteLine("");
            Console.WriteLine("TEST23");
            Console.WriteLine("  LEVEL_GROWTH_TO_ORDER uses Level, Growth and Rule");
            Console.WriteLine("  to determine the orders of each entry of a vector of 1D rules.");
            Console.WriteLine("");
            Console.WriteLine("  Here we examine rule " + r + ".");
            Console.WriteLine("");
            Console.WriteLine("       LEVEL:0     1     2     3     4     5     6     7     8     9    10");
            Console.WriteLine("GROWTH");

            growth = new int[dim_num];
            level = new int[dim_num];
            order = new int[dim_num];
            rule = new int[dim_num];

            for (dim = 0; dim < dim_num; dim++)
            {
                rule[dim] = r;
            }

            for (g = 0; g <= 6; g++)
            {
                switch (r)
                {
                    case 3:
                    case 10:
                    {
                        switch (g)
                        {
                            case 1:
                            case 2:
                            case 3:
                                continue;
                        }

                        break;
                    }
                }

                for (dim = 0; dim < dim_num; dim++)
                {
                    level[dim] = dim;
                }

                for (dim = 0; dim < dim_num; dim++)
                {
                    growth[dim] = g;
                }

                LevelToOrder.level_growth_to_order(dim_num, level, rule, growth, ref order);

                string cout = "  " + g.ToString().PadLeft(4) + "  ";

                for (dim = 0; dim < dim_num; dim++)
                {
                    cout += "  " + order[dim].ToString().PadLeft(4);
                }

                Console.WriteLine(cout);
            }

        }

        private static void test24()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST24 tests LEVEL_TO_ORDER_DEFAULT.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 February 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int dim;
            int dim_num = 11;
            int[] level;
            int[] order;
            int r;
            int[] rule;

            Console.WriteLine("");
            Console.WriteLine("TEST24");
            Console.WriteLine("  LEVEL_TO_ORDER_DEFAULT uses a default rule to");
            Console.WriteLine("  determine the order of a rule from its level.");
            Console.WriteLine("");
            Console.WriteLine("RULE/LEVEL   0     1     2     3     4     5     6     7     8     9    10");
            Console.WriteLine("");

            level = new int[dim_num];
            order = new int[dim_num];
            rule = new int[dim_num];

            for (r = 1; r <= 16; r++)
            {
                for (dim = 0; dim < dim_num; dim++)
                {
                    rule[dim] = r;
                }

                for (dim = 0; dim < dim_num; dim++)
                {
                    level[dim] = dim;
                }

                LevelToOrder.level_to_order_default(dim_num, level, rule, ref order);

                string cout = "  " + r.ToString().PadLeft(4) + "  ";

                for (dim = 0; dim < dim_num; dim++)
                {
                    cout += "  " + order[dim].ToString().PadLeft(4);
                }

                Console.WriteLine(cout);
            }
        }

        private static void test25()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST25 tests LEVEL_TO_ORDER_EXPONENTIAL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    17 December 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int dim;
            int dim_num = 11;
            int[] level;
            int[] order;
            int r;
            int[] rule;

            Console.WriteLine("");
            Console.WriteLine("TEST25");
            Console.WriteLine("  LEVEL_TO_ORDER_EXPONENTIAL uses an exponential rule to");
            Console.WriteLine("  determine the order of a rule from its level.");
            Console.WriteLine("");
            Console.WriteLine("RULE/LEVEL   0     1     2     3     4     5     6     7     8     9    10");
            Console.WriteLine("");

            level = new int[dim_num];
            order = new int[dim_num];
            rule = new int[dim_num];

            for (r = 1; r <= 10; r++)
            {
                for (dim = 0; dim < dim_num; dim++)
                {
                    rule[dim] = r;
                }

                for (dim = 0; dim < dim_num; dim++)
                {
                    level[dim] = dim;
                }

                LevelToOrder.level_to_order_exponential(dim_num, level, rule, ref order);

                string cout = "  " + r.ToString().PadLeft(4) + "  ";

                for (dim = 0; dim < dim_num; dim++)
                {
                    cout += "  " + order[dim].ToString().PadLeft(4);
                }

                Console.WriteLine(cout);
            }
        }

        private static void test26()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST26 tests LEVEL_TO_ORDER_EXPONENTIAL_SLOW.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    17 December 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int dim;
            int dim_num = 11;
            int[] level;
            int[] order;
            int r;
            int[] rule;

            Console.WriteLine("");
            Console.WriteLine("TEST26");
            Console.WriteLine("  LEVEL_TO_ORDER_EXPONENTIAL_SLOW uses a slow exponential rule to");
            Console.WriteLine("  determine the order of a rule from its level.");
            Console.WriteLine("");
            Console.WriteLine("  Since it is really only useful for fully nested rules,");
            Console.WriteLine("  we only consider rules 11, 12 and 13.");
            Console.WriteLine("");
            Console.WriteLine("RULE/LEVEL   0     1     2     3     4     5     6     7     8     9    10");
            Console.WriteLine("");

            level = new int[dim_num];
            order = new int[dim_num];
            rule = new int[dim_num];

            for (r = 11; r <= 13; r++)
            {
                for (dim = 0; dim < dim_num; dim++)
                {
                    rule[dim] = r;
                }

                for (dim = 0; dim < dim_num; dim++)
                {
                    level[dim] = dim;
                }

                LevelToOrder.level_to_order_exponential_slow(dim_num, level, rule, ref order);

                string cout = "  " + r.ToString().PadLeft(4) + "  ";

                for (dim = 0; dim < dim_num; dim++)
                {
                    cout += "  " + order[dim].ToString().PadLeft(4);
                }

                Console.WriteLine(cout);
            }
        }

        private static void test27()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST27 tests LEVEL_TO_ORDER_LINEAR.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    17 December 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int dim;
            int dim_num = 11;
            int[] level;
            int[] order;
            int r;
            int[] rule;

            Console.WriteLine("");
            Console.WriteLine("TEST27");
            Console.WriteLine("  LEVEL_TO_ORDER_LINEAR uses a linear rule to");
            Console.WriteLine("  determine the order of a rule from its level.");
            Console.WriteLine("");
            Console.WriteLine("RULE/LEVEL   0     1     2     3     4     5     6     7     8     9    10");
            Console.WriteLine("");

            level = new int[dim_num];
            order = new int[dim_num];
            rule = new int[dim_num];

            for (r = 1; r <= 10; r++)
            {
                for (dim = 0; dim < dim_num; dim++)
                {
                    rule[dim] = r;
                }

                for (dim = 0; dim < dim_num; dim++)
                {
                    level[dim] = dim;
                }

                LevelToOrder.level_to_order_linear(dim_num, level, rule, ref order);

                string cout = "  " + r.ToString().PadLeft(4) + "  ";

                for (dim = 0; dim < dim_num; dim++)
                {
                    cout += "  " + order[dim].ToString().PadLeft(4);
                }

                Console.WriteLine(cout);
            }
        }

        private static void test28()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST28 tests PATTERSON_LOOKUP against LEGENDRE_INTEGRAL.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    17 October 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double error;
            double estimate;
            double exact;
            double[] f;
            int i;
            int level;
            int level_max = 8;
            int n;
            int order;
            int p;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST28");
            Console.WriteLine("  PATTERSON_LOOKUP computes a Gauss-Patterson rule");
            Console.WriteLine("  which is appropriate for integrands of the form");
            Console.WriteLine("    Integral ( -1 <= x <= +1 ) f(x)  dx.");
            Console.WriteLine("");
            Console.WriteLine("  LEGENDRE_INTEGRAL determines the exact value of");
            Console.WriteLine("  this integal when f(x) = x^n.");
            Console.WriteLine("");
            Console.WriteLine("  A rule of order ORDER should be exact for monomials X^N");
            Console.WriteLine("  up to N = (3*ORDER+1)/2");
            Console.WriteLine("");
            Console.WriteLine("  In the following table, for each order, the LAST THREE estimates");
            Console.WriteLine("  are made on monomials that exceed the exactness limit for the rule.");
            Console.WriteLine("");
            Console.WriteLine("     Order         N       Estimate       Exact            Error");

            for (level = 0; level <= level_max; level++)
            {
                order = (int) Math.Pow(2, level + 1) - 1;

                Console.WriteLine("");

                f = new double[order];
                w = new double[order];
                x = new double[order];

                PattersonQuadrature.patterson_lookup(order, ref x, ref w);

                p = p switch
                {
                    //
                    //  Truncate at 50.
                    //
                    > 50 => 50,
                    _ => order switch
                    {
                        1 => 1,
                        _ => (3 * order + 1) / 2
                    }
                };

                for (n = 0; n <= p + 3; n += 1)
                {
                    exact = Integral.legendre_integral(n);

                    switch (n)
                    {
                        case 0:
                        {
                            for (i = 0; i < order; i++)
                            {
                                f[i] = 1.0;
                            }

                            break;
                        }
                        default:
                        {
                            for (i = 0; i < order; i++)
                            {
                                f[i] = Math.Pow(x[i], n);
                            }

                            break;
                        }
                    }

                    estimate = 0.0;
                    for (i = 0; i < order; i++)
                    {
                        estimate += w[i] * f[i];
                    }

                    error = typeMethods.r8_abs(exact - estimate);

                    Console.WriteLine("  " + order.ToString().PadLeft(8)
                                           + "  " + n.ToString().PadLeft(8)
                                           + "  " + estimate.ToString("0.######").PadLeft(14)
                                           + "  " + exact.ToString("0.######").PadLeft(14)
                                           + "  " + error.ToString("0.######").PadLeft(14) + "");
                }
            }
        }

        private static void test285()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST285 tests POINT_RADIAL_TOL_UNIQUE_INDEX_INC1, INC2 and INC3.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    09 October 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int M = 2;
            int N1 = 11;
            int N2 = 8;

            int m = M;
            int n1 = N1;
            int n2 = N2;

            double[] a1 =
            {
                0.0, 0.0,
                0.5, 0.5,
                1.0, 0.0,
                0.0, 1.0,
                1.0, 1.0,
                0.500000001, 0.5,
                0.0, 0.0,
                0.0, 0.5,
                0.5, 0.0,
                1.0, 0.5,
                0.5, 1.0
            };
            double[] a2 =
            {
                0.4999999999, 0.5,
                0.75, 0.25,
                0.500000001, 0.9999999999,
                0.500000001, 0.0000000001,
                0.25, 0.75,
                0.75, 0.25,
                0.250000001, 0.7499999999,
                0.75, 0.75
            };
            double[] a3;
            int i1;
            int i2;
            int i3;
            int[] indx1;
            int[] indx2;
            int[] indx3;
            int n3 = 0;
            double[] r1;
            double[] r2;
            double[] r3;
            int seed;
            double tol;
            int undx_value;
            int[] undx1;
            int[] undx2;
            int[] undx3;
            int unique_num1 = 0;
            int unique_num2 = 0;
            int unique_num3 = 0;
            bool[] unique1;
            bool[] unique2;
            bool[] unique3;
            int[] xdnu1;
            int[] xdnu2;
            int[] xdnu3;
            double[] z;

            seed = 123456789;

            Console.WriteLine("");
            Console.WriteLine("TEST285");
            Console.WriteLine("  POINT_RADIAL_TOL_UNIQUE_INDEX_INC1 can index unique");
            Console.WriteLine("  points in a \"permanent\" point set;");
            Console.WriteLine("  POINT_RADIAL_TOL_UNIQUE_INDEX_INC2 can incremented by");
            Console.WriteLine("  \"temporary\" points.");
            Console.WriteLine("  POINT_RADIAL_TOL_UNIQUE_INDEX_INC3 can merge permanent");
            Console.WriteLine("  and temporary points.");

            tol = Math.Sqrt(typeMethods.r8_epsilon());
            Console.WriteLine("  Using tolerance TOL = " + tol + "");
            //
            //  Step 1
            //
            indx1 = new int[n1];
            r1 = new double[n1];
            undx1 = new int[n1];
            unique1 = new bool[n1];
            xdnu1 = new int[n1];
            z = new double[m];

            Dataset.point_radial_tol_unique_index_inc1(m, n1, a1, tol, ref seed, z, r1,
                indx1, unique1, ref unique_num1, undx1, xdnu1);

            Console.WriteLine("");
            Console.WriteLine("  UNIQUE_NUM1 = " + unique_num1 + "");
            Console.WriteLine("  Expected   =  " + 9 + "");

            Console.WriteLine("");
            Console.WriteLine("  Item I1, unique index XDNU1[I1], representative location UNDX1[XDNU1[I1]]:");
            Console.WriteLine("");
            Console.WriteLine("           I1  XDNU1  UNDX1");
            Console.WriteLine("");
            for (i1 = 0; i1 < n1; i1++)
            {
                Console.WriteLine("  " + "    "
                                       + "  " + i1.ToString().PadLeft(4)
                                       + "  " + xdnu1[i1].ToString().PadLeft(4)
                                       + "  " + undx1[xdnu1[i1]].ToString().PadLeft(4) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("  Unique item I1, location UNDX1(I1), value A1(:,UNDX1(I1)):");
            Console.WriteLine("");
            Console.WriteLine("          I1 UNDX1  --A1(1,*)---  --A1(2,*)---");
            Console.WriteLine("");
            for (i1 = 0; i1 < unique_num1; i1++)
            {
                Console.WriteLine("  " + "    "
                                       + "  " + i1.ToString().PadLeft(4)
                                       + "  " + undx1[i1].ToString().PadLeft(4)
                                       + "  " + a1[0 + undx1[i1] * m].ToString().PadLeft(12)
                                       + "  " + a1[1 + undx1[i1] * m].ToString().PadLeft(12) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("  Unique item I1, location UNDX1(I1), value A1(:,UNDX1(I1)):");
            Console.WriteLine("");
            Console.WriteLine("          I1 UNIQUE1       R1     --A1(1,I1)--  --A1(2,I1)--");
            Console.WriteLine("");
            for (i1 = 0; i1 < n1; i1++)
            {
                Console.WriteLine("  " + "    "
                                       + "  " + i1.ToString().PadLeft(4)
                                       + "  " + unique1[i1].ToString().PadLeft(4)
                                       + "  " + r1[i1].ToString().PadLeft(12)
                                       + "  " + a1[0 + i1 * m].ToString().PadLeft(12)
                                       + "  " + a1[1 + i1 * m].ToString().PadLeft(12) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("          I1   INDX1   R1(INDX1)");
            Console.WriteLine("");
            for (i1 = 0; i1 < n1; i1++)
            {
                Console.WriteLine("  " + "    "
                                       + "  " + i1.ToString().PadLeft(4)
                                       + "  " + indx1[i1].ToString().PadLeft(4)
                                       + "  " + r1[indx1[i1]].ToString().PadLeft(12) + "");
            }

            //
            //  Step 2
            //
            indx2 = new int[n2];
            r2 = new double[n2];
            undx2 = new int[n2];
            unique2 = new bool[n2];
            xdnu2 = new int[n2];

            Dataset.point_radial_tol_unique_index_inc2(m, n1, a1, n2, a2, tol, z,
                r1, indx1, unique1, unique_num1, undx1, xdnu1,
                ref r2, ref indx2, ref unique2, ref unique_num2, ref undx2, ref xdnu2);

            Console.WriteLine("");
            Console.WriteLine("  UNIQUE_NUM2 = " + unique_num2 + "");
            Console.WriteLine("  Expected   =  " + 3 + "");

            Console.WriteLine("");
            Console.WriteLine("  Item I2, unique index XDNU2[I2], representative location UNDX2[XDNU2[I2]]:");
            Console.WriteLine("");
            Console.WriteLine("    I2 XDNU2 UNDX2");
            Console.WriteLine("");
            Console.WriteLine("  (Temporary data)");
            Console.WriteLine("");
            for (i2 = 0; i2 < n2; i2++)
            {
                if (xdnu2[i2] < unique_num1)
                {
                    undx_value = undx1[xdnu2[i2]];
                }
                else
                {
                    undx_value = undx2[xdnu2[i2] - unique_num1];
                }

                Console.WriteLine("  " + i2.ToString().PadLeft(4)
                                       + "  " + (i2 + n1).ToString().PadLeft(4)
                                       + "  " + xdnu2[i2].ToString().PadLeft(4)
                                       + "  " + undx_value.ToString().PadLeft(4) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("  Unique item I2, location UNDX2(I2), value A2(:,UNDX2(I2)):");
            Console.WriteLine("");
            Console.WriteLine("    I2 UNDX2  --A2(1,*)---  --A2(2,*)---");
            Console.WriteLine("");
            Console.WriteLine("  (Temporary data)");
            Console.WriteLine("");
            for (i2 = 0; i2 < unique_num2; i2++)
            {
                Console.WriteLine("  " + i2.ToString().PadLeft(4)
                                       + "  " + (i2 + unique_num1).ToString().PadLeft(4)
                                       + "  " + undx2[i2].ToString().PadLeft(4)
                                       + "  " + a2[0 + (undx2[i2] - n1) * m].ToString().PadLeft(12)
                                       + "  " + a2[1 + (undx2[i2] - n1) * m].ToString().PadLeft(12) + "");
            }

            //
            //  Step 3.
            //
            Console.WriteLine("");
            Console.WriteLine("  Merge the temporary data with the permanent data.");

            a3 = new double[m * (n1 + n2)];
            indx3 = new int[n1 + n2];
            r3 = new double[n1 + n2];
            undx3 = new int[n1 + n2];
            unique3 = new bool[n1 + n2];
            xdnu3 = new int[n1 + n2];

            Dataset.point_radial_tol_unique_index_inc3(m,
                n1, a1, r1, indx1, unique1, unique_num1, undx1, xdnu1,
                n2, a2, r2, indx2, unique2, unique_num2, undx2, xdnu2,
                ref n3, ref a3, ref r3, ref indx3, ref unique3, ref unique_num3, ref undx3, ref xdnu3);

            Console.WriteLine("");
            Console.WriteLine("  UNIQUE_NUM3 = " + unique_num3 + "");
            Console.WriteLine("  Expected   =  12");

            Console.WriteLine("");
            Console.WriteLine("  Item I3, unique index XDNU3[I3], representative location UNDX3[XDNU3[I3]]:");
            Console.WriteLine("");
            Console.WriteLine("          I3 XDNU3 UNDX3");
            Console.WriteLine("");
            for (i3 = 0; i3 < n3; i3++)
            {
                Console.WriteLine("  " + "    "
                                       + "  " + i3.ToString().PadLeft(4)
                                       + "  " + xdnu3[i3].ToString().PadLeft(4)
                                       + "  " + undx3[xdnu3[i3]].ToString().PadLeft(4) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("  Unique item I3, location UNDX3(I3), value A3(:,UNDX3(I3)):");
            Console.WriteLine("");
            Console.WriteLine("          I3 UNDX3  --A3(1,*)---  --A3(2,*)---");
            Console.WriteLine("");
            for (i3 = 0; i3 < unique_num3; i3++)
            {
                Console.WriteLine("  " + "    "
                                       + "  " + i3.ToString().PadLeft(4)
                                       + "  " + undx3[i3].ToString().PadLeft(4)
                                       + "  " + a3[0 + undx3[i3] * m].ToString().PadLeft(12)
                                       + "  " + a3[1 + undx3[i3] * m].ToString().PadLeft(12) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("  Unique item I3, location UNDX3(I3), value A3(:,UNDX3(I3)):");
            Console.WriteLine("");
            Console.WriteLine("          I3 UNIQUE3       R3     --A3(1,I3)--  --A3(2,I3)--");
            Console.WriteLine("");
            for (i3 = 0; i3 < n3; i3++)
            {
                Console.WriteLine("  " + "    "
                                       + "  " + i3.ToString().PadLeft(4)
                                       + "  " + unique3[i3].ToString().PadLeft(4)
                                       + "  " + r3[i3].ToString().PadLeft(12)
                                       + "  " + a3[0 + i3 * m].ToString().PadLeft(12)
                                       + "  " + a3[1 + i3 * m].ToString().PadLeft(12) + "");
            }

            Console.WriteLine("");
            Console.WriteLine("          I3   INDX3   R3(INDX3)");
            Console.WriteLine("");
            for (i3 = 0; i3 < n3; i3++)
            {
                Console.WriteLine("  " + "    "
                                       + "  " + i3.ToString().PadLeft(4)
                                       + "  " + indx3[i3].ToString().PadLeft(4)
                                       + "  " + r3[indx3[i3]].ToString().PadLeft(12) + "");
            }
        }

        private static void test29()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST29 tests R8COL_TOL_UNDEX.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    19 July 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int M = 3;
            int N = 22;

            double[] a =
            {
                1.9, 0.0, 10.0,
                2.0, 6.0, 10.0,
                4.0, 8.0, 12.0,
                1.0, 5.0, 9.0,
                3.0, 7.0, 11.0,
                2.0, 6.0, 0.0,
                2.0, 0.0, 10.1,
                2.0, 0.1, 10.0,
                3.0, 4.0, 18.0,
                1.9, 8.0, 10.0,
                0.0, 0.0, 0.0,
                0.0, 6.0, 10.0,
                2.1, 0.0, 10.0,
                2.0, 6.0, 10.0,
                3.0, 7.0, 11.0,
                2.0, 0.0, 10.0,
                2.0, 0.0, 10.0,
                2.0, 6.0, 10.0,
                1.0, 5.0, 9.0,
                2.0, 0.0, 10.1,
                1.0, 5.0, 9.1,
                1.0, 5.1, 9.0
            };
            double[] au;
            int i;
            int j;
            int m = M;
            int n = N;
            int n_unique;
            double tol;
            int[] undx;
            int[] xdnu;

            Console.WriteLine("");
            Console.WriteLine("TEST29");
            Console.WriteLine("  R8COL_TOL_UNDEX produces index vectors which create a sorted");
            Console.WriteLine("  list of the tolerably unique columns of an R8COL,");
            Console.WriteLine("  and a map from the original R8COL to the (implicit)");
            Console.WriteLine("  R8COL of sorted tolerably unique elements.");

            typeMethods.r8mat_transpose_print(m, n, a, "  The unsorted R8COL (transposed):");

            tol = 0.25;

            Console.WriteLine("");
            Console.WriteLine("  Using tolerance = " + tol + "");

            n_unique = typeMethods.r8col_tol_unique_count(m, n, a, tol);

            Console.WriteLine("");
            Console.WriteLine("  Number of tolerably unique columns is " + n_unique + "");

            au = new double[m * n_unique];
            undx = new int[n_unique];
            xdnu = new int[n];

            typeMethods.r8col_tol_undex(m, n, a, n_unique, tol, ref undx, ref xdnu);

            Console.WriteLine("");
            Console.WriteLine("  XDNU points to the representative for each item.");
            Console.WriteLine("  UNDX selects the representatives.");
            Console.WriteLine("");
            Console.WriteLine("     I  XDNU  UNDX");
            Console.WriteLine("");
            for (i = 0; i < n_unique; i++)
            {
                Console.WriteLine("  " + i.ToString().PadLeft(4)
                                       + "  " + xdnu[i].ToString().PadLeft(4)
                                       + "  " + undx[i].ToString().PadLeft(4) + "");
            }

            for (i = n_unique; i < n; i++)
            {
                Console.WriteLine("  " + i.ToString().PadLeft(4)
                                       + "  " + xdnu[i].ToString().PadLeft(4) + "");
            }

            for (j = 0; j < n_unique; j++)
            {
                for (i = 0; i < m; i++)
                {
                    au[i + j * m] = a[i + undx[j] * m];
                }
            }

            typeMethods.r8mat_transpose_print(m, n_unique, au,
                "  The tolerably unique R8COL (transposed):");

        }

        private static void test30()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST30 tests R8VEC_SORT_HEAP_INDEX_A_NEW.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    02 October 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double[] a;
            double[] b;
            int i;
            int[] indx;
            int n = 20;
            int seed;

            Console.WriteLine("");
            Console.WriteLine("TEST30");
            Console.WriteLine("  R8VEC_SORT_HEAP_INDEX_A_NEW creates an ascending");
            Console.WriteLine("  sort index for a R8VEC.");

            seed = 123456789;

            a = UniformRNG.r8vec_uniform_01_new(n, ref seed);

            typeMethods.r8vec_print(n, a, "  The unsorted array:");

            indx = typeMethods.r8vec_sort_heap_index_a_new(n, a);

            typeMethods.i4vec_print(n, indx, "  The index vector:");

            b = new double[n];

            for (i = 0; i < n; i++)
            {
                b[i] = a[indx[i]];
            }

            typeMethods.r8vec_print(n, b, "  The sorted array A(INDX(:)):");

        }

        private static void test31()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST31 tests HCE_COMPUTE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    08 March 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int M = 11;
            int N = (2 * M);

            int i;
            int m = M;
            int n = N;
            double[] r = new double[N];
            double[] w = new double[N];
            double[] x = new double[N];

            Console.WriteLine("");
            Console.WriteLine("TEST31:");
            Console.WriteLine("  HCE_COMPUTE returns a quadrature rule");
            Console.WriteLine("  for piecewise Hermite cubic splines which are based");
            Console.WriteLine("  on equally spaced function and derivative data.");
            Console.WriteLine("");
            Console.WriteLine("  Here we compute a rule of order N = " + n + "");

            HermiteQuadrature.hce_compute(n, ref x, ref w);

            Console.WriteLine("");
            Console.WriteLine("     I        X(I)        W(I)");
            Console.WriteLine("");

            for (i = 0; i < n; i++)
            {
                Console.WriteLine("  " + i.ToString().PadLeft(4)
                                       + "  " + x[i].ToString().PadLeft(10)
                                       + "  " + w[i].ToString().PadLeft(10) + "");
            }

        }

        private static void test32()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST32 tests HCC_COMPUTE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    17 March 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int M = 11;
            int N = (2 * M);

            int i;
            int m = M;
            int n = N;
            double[] r = new double[N];
            double[] w = new double[N];
            double[] x = new double[N];

            Console.WriteLine("");
            Console.WriteLine("TEST32:");
            Console.WriteLine("  HCC_COMPUTE returns a quadrature rule");
            Console.WriteLine("  for piecewise Hermite cubic splines which are based");
            Console.WriteLine("  on Chebyshev-spaced function and derivative data.");
            Console.WriteLine("");
            Console.WriteLine("  Here we compute a rule of order N = " + n + "");

            HermiteQuadrature.hcc_compute(n, ref x, ref w);

            Console.WriteLine("");
            Console.WriteLine("     I        X(I)        W(I)");
            Console.WriteLine("");

            for (i = 0; i < n; i++)
            {
                Console.WriteLine("  " + i.ToString().PadLeft(4)
                                       + "  " + x[i].ToString().PadLeft(10)
                                       + "  " + w[i].ToString().PadLeft(10) + "");
            }

        }

        private static void test33()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST33 tests HC_COMPUTE_WEIGHTS_FROM_POINTS.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    31 March 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int N = 11;

            double dn = 0;
            double fn = 0;
            int j;
            int n = N;
            double q;
            double q_exact;
            double[] r;
            double s = 0;
            int seed;
            double t = 0;
            int test;
            double[] w;
            double[] x = new double[N];

            Console.WriteLine("");
            Console.WriteLine("TEST33:");
            Console.WriteLine("  HC_COMPUTE_WEIGHTS_FROM_POINTS returns quadrature weights");
            Console.WriteLine("  given the points.");

            seed = 123456789;

            for (test = 1; test <= 3; test++)
            {
                r = UniformRNG.r8vec_uniform_01_new(n, ref seed);

                x[0] = r[0];
                for (j = 1; j < n; j++)
                {
                    x[j] = x[j - 1] + r[j];
                }

                Console.WriteLine("");
                Console.WriteLine("  Trial #" + test + ":");
                Console.WriteLine("  Random spacing");
                Console.WriteLine("  Number of points N = " + n + "");
                Console.WriteLine("  Interval = [" + x[0] + ", " + x[n - 1] + "]");

                w = new double [2 * n];

                HermiteQuadrature.hc_compute_weights_from_points(n, x, ref w);

                q = 0.0;

                for (j = 0; j < n; j++)
                {
                    cubic_value(x[j], ref fn, ref dn, ref s, ref t);
                    q = q + w[0 + j * 2] * fn + w[1 + j * 2] * dn;
                }

                q_exact = cubic_integrate(x[0], x[n - 1]);

                Console.WriteLine("");
                Console.WriteLine("  Q         = " + q + "");
                Console.WriteLine("  Q (exact) = " + q_exact + "");

            }
        }

        private static void test34()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST34 uses HERMITE_INTERPOLANT on the Runge function using equally spaced data.
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
            double[] ysp;
            double yt;

            Console.WriteLine("");
            Console.WriteLine("TEST34");
            Console.WriteLine("  HERMITE_INTERPOLANT computes the Hermite interpolant to data.");
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

                ys = new double[ns];
                ysp = new double[ns];

                Hermite.hermite_interpolant_value(nd, xd, yd, xdp, ydp, ns, xs, ref ys, ref ysp);

                max_dif = 0.0;
                for (i = 0; i < ns; i++)
                {
                    xt = xs[i];
                    yt = 1.0 / (1.0 + xt * xt);
                    max_dif = Math.Max(max_dif, typeMethods.r8_abs(ys[i] - yt));
                }

                Console.WriteLine("  " + n.ToString().PadLeft(4)
                                       + "  " + max_dif.ToString().PadLeft(14) + "");

            }
        }

        private static void test35()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST35 uses HERMITE_INTERPOLANT on the Runge function using Chebyshev spaced data.
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
            double[] ysp;
            double yt;

            Console.WriteLine("");
            Console.WriteLine("TEST35");
            Console.WriteLine("  HERMITE_INTERPOLANT computes the Hermite interpolant to data.");
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

                ys = new double[ns];
                ysp = new double[ns];

                Hermite.hermite_interpolant_value(nd, xd, yd, xdp, ydp, ns, xs, ref ys, ref ysp);

                max_dif = 0.0;
                for (i = 0; i < ns; i++)
                {
                    xt = xs[i];
                    yt = 1.0 / (1.0 + xt * xt);
                    max_dif = Math.Max(max_dif, typeMethods.r8_abs(ys[i] - yt));
                }

                Console.WriteLine("  " + n.ToString().PadLeft(4)
                                       + "  " + max_dif.ToString().PadLeft(14) + "");

            }
        }

        private static void test36()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST36 tests HERMITE_INTERPOLANT_RULE.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 October 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double a;
            double b;
            int e;
            double error;
            double exact;
            int i;
            int k;
            int n;
            double q;
            double[] w;
            double[] x;

            Console.WriteLine("");
            Console.WriteLine("TEST36:");
            Console.WriteLine("  HERMITE_INTERPOLANT_RULE");
            Console.WriteLine("  is given a set of N abscissas for a Hermite interpolant");
            Console.WriteLine("  and returns N pairs of quadrature weights");
            Console.WriteLine("  for function and derivative values at the abscissas.");
            //
            //  1: Behavior with increasing N.
            //
            a = 0.0;
            b = 1.0;

            Console.WriteLine("");
            Console.WriteLine("  Observe behavior of quadrature weights for increasing N");
            Console.WriteLine("  We are working in " + a + " <= X <= " + b + "");

            for (n = 3; n <= 11; n += 2)
            {
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
            }

            //
            //  2: Integral estimates with equally spaced points.
            //
            a = -5.0;
            b = +5.0;
            n = 11;

            Console.WriteLine("");
            Console.WriteLine("  Use the rule with N = " + n + " to estimate integrals.");
            Console.WriteLine("  Points are equally spaced.");
            Console.WriteLine("  We are working in " + a + " <= X <= " + b + "");

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

            q = 0.0;
            k = 0;
            for (i = 0; i < n; i++)
            {
                q = q + w[k] * 1 + w[k + 1] * 0.0;
                k += 2;
            }

            Console.WriteLine("");
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

            //
            //  3: Integral estimates with Chebyshev spaced points.
            //
            a = -5.0;
            b = +5.0;
            n = 11;

            Console.WriteLine("");
            Console.WriteLine("  Use the rule with N = " + n + " to estimate integrals.");
            Console.WriteLine("  Points are Chebyshev spaced.");
            Console.WriteLine("  We are working in " + a + " <= X <= " + b + "");

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

            q = 0.0;
            k = 0;
            for (i = 0; i < n; i++)
            {
                q = q + w[k] * 1 + w[k + 1] * 0.0;
                k += 2;
            }

            Console.WriteLine("");
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

            //
            //  4: Integral estimates with Legendre spaced points.
            //
            a = -5.0;
            b = +5.0;
            n = 11;

            Console.WriteLine("");
            Console.WriteLine("  Use the rule with N = " + n + " to estimate integrals.");
            Console.WriteLine("  Points are Legendre spaced.");
            Console.WriteLine("  We are working in " + a + " <= X <= " + b + "");

            x = typeMethods.r8vec_legendre_new(n, a, b);
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

            q = 0.0;
            k = 0;
            for (i = 0; i < n; i++)
            {
                q = q + w[k] * 1 + w[k + 1] * 0.0;
                k += 2;
            }

            Console.WriteLine("");
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

            //
            //  5: Integral estimates with Chebyshev spaced points on 1/(1+x^2), increasing N.
            //
            a = -5.0;
            b = +5.0;

            Console.WriteLine("");
            Console.WriteLine("  Approximate integral of 1/(1+x^2) with increasing N.");
            Console.WriteLine("  Points are Chebyshev spaced.");
            Console.WriteLine("  We are working in " + a + " <= X <= " + b + "");
            Console.WriteLine("");
            Console.WriteLine("     N     Estimate         Error");
            Console.WriteLine("");

            for (n = 2; n <= 11; n++)
            {
                x = typeMethods.r8vec_chebyshev_new(n, a, b);
                w = Hermite.hermite_interpolant_rule(n, a, b, x);

                q = 0.0;
                k = 0;
                for (i = 0; i < n; i++)
                {
                    q = q + w[k] / (1.0 + x[i] * x[i])
                        - w[k + 1] * 2.0 * x[i] / Math.Pow(1.0 + x[i] * x[i], 2);
                    k += 2;
                }

                exact = Math.Atan(b) - Math.Atan(a);
                error = typeMethods.r8_abs(q - exact);
                Console.WriteLine("  " + n.ToString().PadLeft(4)
                                       + "  " + q.ToString().PadLeft(14)
                                       + "  " + error.ToString().PadLeft(14) + "");

            }

            //
            //  6: Integral estimates, with Chebyshev spaced points, for monomials, using N = 11.
            //
            a = -1.0;
            b = 1.0;

            Console.WriteLine("");
            Console.WriteLine("  Compute integral estimates for monomials X^0 through X^15.");
            Console.WriteLine("  Use N = 5, 9, 13, 17, 21 point rules.");
            Console.WriteLine("  Points are Chebyshev spaced.");
            Console.WriteLine("  We are working in " + a + " <= X <= " + b + "");

            for (n = 5; n <= 21; n += 4)
            {
                x = typeMethods.r8vec_chebyshev_new(n, a, b);
                w = Hermite.hermite_interpolant_rule(n, a, b, x);

                Console.WriteLine("");
                Console.WriteLine("  Estimates are made using N = " + n + "");
                Console.WriteLine("  F(X)         Integral        Estimate           Error");
                Console.WriteLine("");
                for (e = 0; e <= 15; e++)
                {
                    q = 0.0;
                    k = 0;
                    for (i = 0; i < n; i++)
                    {
                        switch (e)
                        {
                            case 0:
                                q += w[k];
                                break;
                            default:
                                q = q + w[k] * Math.Pow(x[i], e) + w[k + 1] * e * Math.Pow(x[i], e - 1);
                                break;
                        }

                        k += 2;
                    }

                    exact = (Math.Pow(b, e + 1) - Math.Pow(a, e + 1)) / (double) (e + 1);
                    Console.WriteLine("  X^" + e.ToString().PadLeft(7)
                                             + "  " + exact.ToString().PadLeft(14)
                                             + "  " + q.ToString().PadLeft(14)
                                             + "  " + typeMethods.r8_abs(exact - q).ToString().PadLeft(14) + "");
                }

                q = 0.0;
                k = 0;
                for (i = 0; i < n; i++)
                {
                    q = q + w[k] / (1.0 + x[i] * x[i])
                        - w[k + 1] * 2.0 * x[i] / Math.Pow(1.0 + x[i] * x[i], 2);
                    k += 2;
                }

                exact = Math.Atan(b) - Math.Atan(a);
                Console.WriteLine("  1/(1+x^2)"
                                  + "  " + exact.ToString().PadLeft(14)
                                  + "  " + q.ToString().PadLeft(14)
                                  + "  " + typeMethods.r8_abs(exact - q).ToString().PadLeft(14) + "");

            }
        }

        private static void test37()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST37 checks that the HGK weights are correctly scaled.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    04 October 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int o;
            int[] order = {1, 3, 9, 19, 35, 37, 41, 43};
            double sqrtpi = 1.7724538509055159;
            int rule;
            double s;
            double[] w;

            Console.WriteLine("");
            Console.WriteLine("TEST37");
            Console.WriteLine("  HERMITE_GENZ_KEISTER_LOOKUP_WEIGHTS looks up weights");
            Console.WriteLine("  for Genz-Keister quadrature rules for the Hermite weight function.");
            Console.WriteLine("");
            Console.WriteLine("  This test simply checks that, for each rule, the quadrature");
            Console.WriteLine("  weights correctly sum to sqrt(pi).");
            Console.WriteLine("");
            Console.WriteLine(" Index     Order         Sum of Weights");
            Console.WriteLine("");

            for (rule = 0; rule < 8; rule++)
            {
                o = order[rule];

                w = new double[o];

                Burkardt.PolynomialNS.Hermite.hermite_genz_keister_lookup_weights(o, ref w);

                s = typeMethods.r8vec_sum(o, w);

                Console.WriteLine("  " + rule.ToString().PadLeft(4)
                                       + "  " + o.ToString().PadLeft(8)
                                       + "  " + s.ToString("0.######").PadLeft(14) + "");

            }

            Console.WriteLine("");
            Console.WriteLine(" Correct sum:            " + sqrtpi + "");
        }

        private static void test38()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST38 tabulates the Hermite interpolant and its derivative. 
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
            Console.WriteLine("TEST38");
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

        private static void test39()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST39 tests the LEVEL_TO_ORDER_** functions.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    31 December 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int growth;
            int level;
            int order;

            Console.WriteLine("");
            Console.WriteLine("TEST39:");
            Console.WriteLine("  Test the various LEVEL_TO_ORDER_** functions,");
            Console.WriteLine("  which, for a given type of quadrature rule, accept");
            Console.WriteLine("  LEVEL, the index of the rule in the family, and");
            Console.WriteLine("  GROWTH (0=slow, 1=moderate, 2=full), a growth rate, and");
            Console.WriteLine("  return the appropriate corresponding order for the rule.");
            Console.WriteLine("");
            Console.WriteLine("  LEVEL_TO_ORDER_EXP_CC:");
            Console.WriteLine("  Slow/moderate/full exponential growth typical of");
            Console.WriteLine("  the Clenshaw Curtis rule.");
            Console.WriteLine("");
            Console.WriteLine("  LEVEL  ORDER  ORDER  ORDER");
            Console.WriteLine("         G = 0  G = 1  G = 2");
            Console.WriteLine("");

            for (level = 0; level <= 8; level++)
            {
                string cout = "  " + level.ToString().PadLeft(5);
                for (growth = 0; growth <= 2; growth++)
                {
                    order = LevelToOrder.level_to_order_exp_cc(level, growth);
                    cout += "  " + order.ToString().PadLeft(5);
                }

                Console.WriteLine(cout);
            }

            Console.WriteLine("");
            Console.WriteLine("  LEVEL_TO_ORDER_EXP_F2:");
            Console.WriteLine("  Slow/moderate/full exponential growth typical of");
            Console.WriteLine("  the Fejer Type 2 rule.");
            Console.WriteLine("");
            Console.WriteLine("  LEVEL  ORDER  ORDER  ORDER");
            Console.WriteLine("         G = 0  G = 1  G = 2");
            Console.WriteLine("");

            for (level = 0; level <= 8; level++)
            {
                string cout = "  " + level.ToString().PadLeft(5);
                for (growth = 0; growth <= 2; growth++)
                {
                    order = LevelToOrder.level_to_order_exp_f2(level, growth);
                    cout += "  " + order.ToString().PadLeft(5);
                }

                Console.WriteLine(cout);
            }

            Console.WriteLine("");
            Console.WriteLine("  LEVEL_TO_ORDER_EXP_GAUSS:");
            Console.WriteLine("  Slow/moderate/full exponential growth typical of");
            Console.WriteLine("  a Gauss rule.");
            Console.WriteLine("");
            Console.WriteLine("  LEVEL  ORDER  ORDER  ORDER");
            Console.WriteLine("         G = 0  G = 1  G = 2");
            Console.WriteLine("");

            for (level = 0; level <= 8; level++)
            {
                string cout = "  " + level.ToString().PadLeft(5);
                for (growth = 0; growth <= 2; growth++)
                {
                    order = LevelToOrder.level_to_order_exp_gauss(level, growth);
                    cout += "  " + order.ToString().PadLeft(5);
                }

                Console.WriteLine(cout);
            }

            Console.WriteLine("");
            Console.WriteLine("  LEVEL_TO_ORDER_EXP_GP:");
            Console.WriteLine("  Slow/moderate/full exponential growth typical of");
            Console.WriteLine("  a Gauss-Patterson rule.");
            Console.WriteLine("");
            Console.WriteLine("  LEVEL  ORDER  ORDER  ORDER");
            Console.WriteLine("         G = 0  G = 1  G = 2");
            Console.WriteLine("");

            for (level = 0; level <= 8; level++)
            {
                string cout = "  " + level.ToString().PadLeft(5);
                for (growth = 0; growth <= 2; growth++)
                {
                    order = LevelToOrder.level_to_order_exp_gp(level, growth);
                    cout += "  " + order.ToString().PadLeft(5);
                }

                Console.WriteLine(cout);
            }

            Console.WriteLine("");
            Console.WriteLine("  LEVEL_TO_ORDER_EXP_HGK:");
            Console.WriteLine("  Slow/moderate/full exponential growth typical of");
            Console.WriteLine("  a Hermite Genz-Keister rule.");
            Console.WriteLine("");
            Console.WriteLine("  LEVEL  ORDER  ORDER  ORDER");
            Console.WriteLine("         G = 0  G = 1  G = 2");
            Console.WriteLine("");

            for (level = 0; level <= 8; level++)
            {
                for (growth = 0; growth <= 2; growth++)
                {
                    switch (growth)
                    {
                        case 2 when 5 < level:
                            break;
                        default:
                            order = LevelToOrder.level_to_order_exp_hgk(level, growth);
                            break;
                    }
                }

                Console.WriteLine("");
            }

            Console.WriteLine("");
            Console.WriteLine("  LEVEL_TO_ORDER_LINEAR_NN:");
            Console.WriteLine("  Slow/moderate linear growth typical of");
            Console.WriteLine("  a non-nested Gauss rule.");
            Console.WriteLine("");
            Console.WriteLine("  LEVEL  ORDER  ORDER");
            Console.WriteLine("         G = 0  G = 1");
            Console.WriteLine("");

            for (level = 0; level <= 8; level++)
            {
                string cout = "  " + level.ToString().PadLeft(5);
                for (growth = 0; growth <= 1; growth++)
                {
                    order = LevelToOrder.level_to_order_linear_nn(level, growth);
                    cout += "  " + order.ToString().PadLeft(5);
                }

                Console.WriteLine(cout);
            }

            Console.WriteLine("");
            Console.WriteLine("  LEVEL_TO_ORDER_LINEAR_WN:");
            Console.WriteLine("  Slow/moderate linear growth typical of");
            Console.WriteLine("  a weakly-nested Gauss rule.");
            Console.WriteLine("");
            Console.WriteLine("  LEVEL  ORDER  ORDER");
            Console.WriteLine("         G = 0  G = 1");
            Console.WriteLine("");

            for (level = 0; level <= 8; level++)
            {
                string cout = "  " + level.ToString().PadLeft(5);
                for (growth = 0; growth <= 1; growth++)
                {
                    order = LevelToOrder.level_to_order_linear_wn(level, growth);
                    cout += "  " + order.ToString().PadLeft(5);
                }

                Console.WriteLine(cout);
            }
        }

        private static double cubic_antiderivative(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CUBIC_ANTIDERIVATIVE evaluates the antiderivative function of a cubic.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    28 January 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double X, the argument.
            //
            //    Output, double CUBIC_ANTIDERIVATIVE, the value.
            //
        {
            double value = 0;

            value = x * x * (5.0 + x * (-7.0 / 3.0 + x * 1.0 / 4.0));

            return value;
        }

        private static double cubic_integrate(double a, double b)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CUBIC_INTEGRATE integrates the cubic from A to B.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    15 February 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double A, B, the integration interval.
            //
            //    Output, double Q, the integral from A to B.
            //
        {
            double q;

            q = cubic_antiderivative(b) - cubic_antiderivative(a);

            return q;
        }

        private static void cubic_value(double x, ref double f, ref double d, ref double s, ref double t)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CUBIC_VALUE evaluates a cubic function.
            //
            //  Discussion:
            //
            //    f(x) =   x^3 -  7 x^2 + 10 x
            //    d(x) = 3 x^2 - 14 x   + 10
            //    s(x) = 6 x   - 14
            //    t(x) = 6
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    15 February 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double X, the argument.
            //
            //    Output, double F, D, S, T, the value and first three
            //    derivatives of the cubic function.
            //
        {
            f = 0.0 + x * (10.0 + x * (-7.0 + x * 1.0));
            d = 10.0 + x * (-14.0 + x * 3.0);
            s = -14.0 + x * 6.0;
            t = 6.0;

        }
    }
}