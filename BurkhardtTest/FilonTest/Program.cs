using System;
using Burkardt.IntegralNS;

namespace FilonTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for FILON_TEST.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    20 May 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            Console.WriteLine("");
            Console.WriteLine("FILON_TEST");
            Console.WriteLine("  Test the FILON library.");

            test01();
            test02();
            test03();
            test04();
            test05();
            test06();

            Console.WriteLine("");
            Console.WriteLine("FILON_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void test01()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST01 tests FILON_TAB_COS.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    20 May 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double a;
            double b;
            double exact = 0;
            double[] ftab = new double[1];
            int i;
            int k;
            int n;
            const double r8_pi = 3.141592653589793;
            double result;
            double t = 0;
            double[] x;

            a = 0.0;
            b = 2.0 * r8_pi;

            n = 11;
            //
            //  Set the X values.
            //
            x = new double[n];
            for (i = 0; i < n; i++)
            {
                x[i] = ((double) (n - i - 1) * a
                        + (double) (i) * b)
                       / (double) (n - 1);
            }

            Console.WriteLine("");
            Console.WriteLine("TEST01");
            Console.WriteLine("  FILON_TAB_COS estimates the integral of.");
            Console.WriteLine("  F(X) * COS ( T * X )");
            Console.WriteLine("  Use integrands F(X)=1, X, X^2.");
            Console.WriteLine("");
            Console.WriteLine("  A = " + a + "");
            Console.WriteLine("  B = " + b + "");
            Console.WriteLine("  N = " + n + "");
            Console.WriteLine("");
            Console.WriteLine("       T                      Approximate             Exact");

            for (k = 1; k <= 3; k++)
            {
                if (k == 1)
                {
                    t = 1.0;
                }
                else if (k == 2)
                {
                    t = 2.0;
                }
                else if (k == 3)
                {
                    t = 10.0;
                }

                Console.WriteLine("");

                for (i = 1; i <= 3; i++)
                {
                    if (i == 1)
                    {
                        ftab = zero_integrand(n, x);
                    }
                    else if (i == 2)
                    {
                        ftab = one_integrand(n, x);
                    }
                    else if (i == 3)
                    {
                        ftab = two_integrand(n, x);
                    }

                    result = Filon.filon_tab_cos(n, ftab, a, b, t);

                    if (i == 1)
                    {
                        exact = (Math.Sin(t * b) - Math.Sin(t * a)) / t;
                    }
                    else if (i == 2)
                    {
                        exact = ((Math.Cos(t * b) + t * b * Math.Sin(t * b))
                                 - (Math.Cos(t * a) + t * a * Math.Sin(t * a))) / t / t;
                    }
                    else if (i == 3)
                    {
                        exact = ((2.0 * t * b * Math.Cos(t * b)
                                  + (t * t * b * b - 2.0) * Math.Sin(t * b))
                                 - (2.0 * t * a * Math.Cos(t * a)
                                    + (t * t * a * a - 2.0) * Math.Sin(t * a))) / t / t / t;
                    }

                    Console.WriteLine(t.ToString().PadLeft(24) + "  "
                        + result.ToString().PadLeft(24) + "  "
                        + exact.ToString().PadLeft(24) + "");
                }
            }
        }

        static void test02()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST02 tests FILON_TAB_COS.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    20 May 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double a;
            double b;
            double error;
            double exact;
            double[] ftab;
            int i;
            int j;
            int n;
            const double r8_pi = 3.141592653589793;
            double result;
            double t;
            double[] x;
            //
            //  Example suggested by James Roedder.
            //
            Console.WriteLine("");
            Console.WriteLine("TEST02");
            Console.WriteLine("  Integrate F(X) = log(1+X)*Math.Cos(T*X):");
            Console.WriteLine("  Supply integrand as a table.");
            Console.WriteLine("  T = 10, and N increases");
            Console.WriteLine("");
            Console.WriteLine("       N    Approximate             Exact                   Error");
            Console.WriteLine("");

            a = 0.0;
            b = 2.0 * r8_pi;

            for (j = 1; j <= 6; j++)
            {
                n = (int)Math.Pow(2, j) * 10 + 1;
                //
                //  Set the X values.
                //
                x = new double[n];
                for (i = 0; i < n; i++)
                {
                    x[i] = ((double) (n - i - 1) * a
                            + (double) (i) * b)
                           / (double) (n - 1);
                }

                ftab = log_integrand(n, x);

                t = 10.0;

                result = Filon.filon_tab_cos(n, ftab, a, b, t);

                exact = -0.008446594405;
                error = result - exact;

                Console.WriteLine(n.ToString().PadLeft(6) + "  "
                    + result.ToString().PadLeft(24) + "  "
                    + exact.ToString().PadLeft(24) + "  "
                    + error.ToString().PadLeft(24) + "");

            }
        }

        static void test03()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST03 tests FILON_FUN_COS.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    20 May 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double a;
            double b;
            double error;
            double exact;
            int j;
            int n;
            const double r8_pi = 3.141592653589793;
            double result;
            double t;
            //
            //  Example suggested by James Roedder.
            //
            Console.WriteLine("");
            Console.WriteLine("TEST03");
            Console.WriteLine("  Integrate F(X)=log(1+X)*Math.Cos(T*X):");
            Console.WriteLine("  Supply integrand as a function.");
            Console.WriteLine("  T = 10, and N increases");
            Console.WriteLine("");
            Console.WriteLine("       N    Approximate             Exact                   Error");
            Console.WriteLine("");

            a = 0.0;
            b = 2.0 * r8_pi;

            for (j = 1; j <= 6; j++)
            {
                n = (int)Math.Pow(2, j) * 10 + 1;

                t = 10.0;

                result = Filon.filon_fun_cos(n, log_integrand, a, b, t);

                exact = -0.008446594405;
                error = result - exact;

                Console.WriteLine(n.ToString().PadLeft(6) + "  "
                    + result.ToString().PadLeft(24) + "  "
                    + exact.ToString().PadLeft(24) + "  "
                    + error.ToString().PadLeft(24) + "");
            }
        }

        static void test04()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST04 tests FILON_TAB_SIN.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    20 May 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double a;
            double b;
            double exact = 0;
            double[] ftab = new double[1];
            int i;
            int k;
            int n;
            double r8_pi = 3.141592653589793;
            double result;
            double t = 0;
            double[] x;

            a = 0.0;
            b = 2.0 * r8_pi;
            n = 11;
            //
            //  Set the X values.
            //
            x = new double[n];
            for (i = 0; i < n; i++)
            {
                x[i] = ((double) (n - i - 1) * a
                        + (double) (i) * b)
                       / (double) (n - 1);
            }

            Console.WriteLine("");
            Console.WriteLine("TEST04");
            Console.WriteLine("  FILON_TAB_SIN estimates the integral of.");
            Console.WriteLine("  F(X) * SIN ( T * X )");
            Console.WriteLine("  Use integrands 1, X, X^2.");
            Console.WriteLine("");
            Console.WriteLine("  A = " + a + "");
            Console.WriteLine("  B = " + b + "");
            Console.WriteLine("  N = " + n + "");
            Console.WriteLine("");
            Console.WriteLine("       T                      Approximate             Exact");
            Console.WriteLine("");

            for (k = 1; k <= 3; k++)
            {
                if (k == 1)
                {
                    t = 1.0;
                }
                else if (k == 2)
                {
                    t = 2.0;
                }
                else if (k == 3)
                {
                    t = 10.0;
                }

                Console.WriteLine("");

                for (i = 1; i <= 3; i++)
                {
                    if (i == 1)
                    {
                        ftab = zero_integrand(n, x);
                    }
                    else if (i == 2)
                    {
                        ftab = one_integrand(n, x);
                    }
                    else if (i == 3)
                    {
                        ftab = two_integrand(n, x);
                    }

                    result = Filon.filon_tab_sin(n, ftab, a, b, t);

                    if (i == 1)
                    {
                        exact = (-Math.Cos(t * b) + Math.Cos(t * a)) / t;
                    }
                    else if (i == 2)
                    {
                        exact = ((Math.Sin(t * b) - t * b * Math.Cos(t * b))
                                 - (Math.Sin(t * a) - t * a * Math.Cos(t * a))) / t / t;
                    }
                    else if (i == 3)
                    {
                        exact = ((2.0 * t * b * Math.Sin(t * b)
                                  + (2.0 - t * t * b * b) * Math.Cos(t * b))
                                 - (2.0 * t * a * Math.Sin(t * a)
                                    + (2.0 - t * t * a * a) * Math.Cos(t * a))) / t / t / t;
                    }

                    Console.WriteLine(t.ToString().PadLeft(24) + "  "
                        + result.ToString().PadLeft(24) + "  "
                        + exact.ToString().PadLeft(24) + "");
                }
            }


            return;
        }
        //****************************************************************************80

        static void test05()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST05 tests FILON_TAB_COS.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    20 May 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double a;
            double b;
            double error;
            double exact;
            double[] ftab;
            int i;
            int j;
            int n;
            const double r8_pi = 3.141592653589793;
            double result;
            double t;
            double[] x;
            //
            //  Example suggested by James Roedder.
            //
            Console.WriteLine("");
            Console.WriteLine("TEST05");
            Console.WriteLine("  Integrate F(X)=log(1+X)*Math.Sin(T*X):");
            Console.WriteLine("  Supply integrand as a table.");
            Console.WriteLine("  T = 10, and N increases");
            Console.WriteLine("");
            Console.WriteLine("       N    Approximate             Exact                   Error");
            Console.WriteLine("");

            a = 0.0;
            b = 2.0 * r8_pi;

            for (j = 1; j <= 6; j++)
            {
                n = (int)Math.Pow(2, j) * 10 + 1;
                //
                //  Set the X values.
                //
                x = new double[n];
                for (i = 0; i < n; i++)
                {
                    x[i] = ((double) (n - i - 1) * a
                            + (double) (i) * b)
                           / (double) (n - 1);
                }

                ftab = log_integrand(n, x);

                t = 10.0;

                result = Filon.filon_tab_sin(n, ftab, a, b, t);

                exact = -0.19762680771872;
                error = result - exact;

                Console.WriteLine(n.ToString().PadLeft(6) + "  "
                    + result.ToString().PadLeft(24) + "  "
                    + exact.ToString().PadLeft(24) + "  "
                    + error.ToString().PadLeft(24) + "");

            }
        }

        static void test06()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST06 tests FILON_FUN_COS.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    20 May 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double a;
            double b;
            double error;
            double exact;
            int j;
            int n;
            const double r8_pi = 3.141592653589793;
            double result;
            double t;
            //
            //  Example suggested by James Roedder.
            //
            Console.WriteLine("");
            Console.WriteLine("TEST06");
            Console.WriteLine("  Integrate F(X)=log(1+X)*Math.Sin(T*X):");
            Console.WriteLine("  Supply integrand as a function.");
            Console.WriteLine("  T = 10, and N increases");
            Console.WriteLine("");
            Console.WriteLine("       N    Approximate             Exact                   Error");
            Console.WriteLine("");

            a = 0.0;
            b = 2.0 * r8_pi;

            for (j = 1; j <= 6; j++)
            {
                n = (int)Math.Pow(2, j) * 10 + 1;

                t = 10.0;

                result = Filon.filon_fun_sin(n, log_integrand, a, b, t);

                exact = -0.19762680771872;
                error = result - exact;

                Console.WriteLine(n.ToString().PadLeft(6) + "  "
                    + result.ToString().PadLeft(24) + "  "
                    + exact.ToString().PadLeft(24) + "  "
                    + error.ToString().PadLeft(24) + "");
            }
        }

        static double[] zero_integrand(int n, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ZERO_INTEGRAND evaluates the integrand x^0.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    20 May 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of points.
            //
            //    Input, double X[N], the evaluation points.
            //
            //    Output, double ZERO_INTEGRAND[N], the function values.
            //
        {
            double[] fx;
            int i;

            fx = new double[n];

            for (i = 0; i < n; i++)
            {
                fx[i] = 1.0;
            }

            return fx;
        }

        static double[] one_integrand(int n, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ONE_INTEGRAND evaluates the integrand X.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    20 May 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of points.
            //
            //    Input, double X[N], the evaluation points.
            //
            //    Output, double ONE_INTEGRAND[N], the function values.
            //
        {
            double[] fx;
            int i;

            fx = new double[n];

            for (i = 0; i < n; i++)
            {
                fx[i] = x[i];
            }

            return fx;
        }

        static double[] two_integrand(int n, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TWO_INTEGRAND evaluates the integrand X^2.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    20 May 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of points.
            //
            //    Input, double X[N], the evaluation points.
            //
            //    Output, double TWO_INTEGRAND[N], the function values.
            //
        {
            double[] fx;
            int i;

            fx = new double[n];

            for (i = 0; i < n; i++)
            {
                fx[i] = x[i] * x[i];
            }

            return fx;
        }

        static double[] log_integrand(int n, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LOG_INTEGRAND evaluates the logarithmic integrand.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    20 May 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of points.
            //
            //    Input, double X[N], the evaluation points.
            //
            //    Output, double LOG_INTEGRAND[N], the function values.
            //
        {
            double[] fx;
            int i;

            fx = new double[n];

            for (i = 0; i < n; i++)
            {
                fx[i] = Math.Log(1.0 + x[i]);
            }

            return fx;
        }
    }
}