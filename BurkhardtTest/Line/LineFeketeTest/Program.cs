using System;
using Burkardt.LineNS;
using Burkardt.Types;

namespace LineFeketeTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for LINE_FEKETE_RULE_TEST.
            //
            //  Discussion:
            //
            //    LINE_FEKETE_RULE_TEST tests the LINE_FEKETE_RULE library.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 April 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            int TEST_NUM = 3;

            int m;
            int[] m_test = {5, 11, 21};
            int test;
            int test_num = TEST_NUM;

            Console.WriteLine("");
            Console.WriteLine("LINE_FEKETE_RULE_TEST");
            Console.WriteLine("  Test the LINE_FEKETE_RULE library.");

            for (test = 0; test < test_num; test++)
            {
                m = m_test[test];
                test01(m);
            }

            for (test = 0; test < test_num; test++)
            {
                m = m_test[test];
                test02(m);
            }

            for (test = 0; test < test_num; test++)
            {
                m = m_test[test];
                test03(m);
            }

            Console.WriteLine("");
            Console.WriteLine("LINE_FEKETE_RULE_TEST");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void test01(int m)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST01 seeks Fekete points in [-1,+1].
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 April 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Alvise Sommariva, Marco Vianello,
            //    Computing approximate Fekete points by QR factorizations of Vandermonde 
            //    matrices,
            //    Computers and Mathematics with Applications,
            //    Volume 57, 2009, pages 1324-1336.
            //
            //  Parameters:
            //
            //    Input, int M, the dimension of the polynomial space.
            //
        {
            int N = 5001;

            double a;
            double b;
            int n = N;
            int nf = 0;
            double[] wf;
            double wf_sum;
            double[] x;
            double[] xf;

            a = -1.0;
            b = +1.0;
            x = typeMethods.r8vec_linspace_new(n, a, b);

            Console.WriteLine("");
            Console.WriteLine("TEST01:");
            Console.WriteLine("  Seek Fekete points in [" + a + "," + b + "]");
            Console.WriteLine("  using " + n + " equally spaced sample points");
            Console.WriteLine("  for polynomials of degree M = " + m + "");
            Console.WriteLine("  using the monomial basis and uniform weight.");

            wf = new double[m];
            xf = new double[m];
            LineFekete.line_fekete_monomial(m, a, b, n, x, ref nf, ref xf, ref wf);

            Console.WriteLine("");
            Console.WriteLine("  NF = " + nf + "");
            typeMethods.r8vec_print(nf, xf, "  Estimated Fekete points XF:");

            wf_sum = typeMethods.r8vec_sum(nf, wf);
            Console.WriteLine("");
            Console.WriteLine("  Sum(WF) = " + wf_sum + "");
        }

        static void test02(int m)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST02 seeks Fekete points in [-1,+1].
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 April 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    L Bos, N Levenberg,
            //    On the calculation of approximate Fekete points: the univariate case,
            //    Electronic Transactions on Numerical Analysis,
            //    Volume 30, pages 377-397, 2008.
            //
            //  Parameters:
            //
            //    Input, int M, the dimension of the polynomial space.
            //
        {
            int N = 5001;

            double a;
            double b;
            int n = N;
            int nf = 0;
            double[] wf;
            double wf_sum;
            double[] x;
            double[] xf;

            a = -1.0;
            b = +1.0;
            x = typeMethods.r8vec_linspace_new(n, a, b);

            Console.WriteLine("");
            Console.WriteLine("TEST02:");
            Console.WriteLine("  Seek Fekete points in [" + a + "," + b + "]");
            Console.WriteLine("  using " + n + " equally spaced sample points");
            Console.WriteLine("  for polynomials of degree M = " + m + "");
            Console.WriteLine("  with the Chebyshev basis.");

            wf = new double[m];
            xf = new double[m];
            LineFekete.line_fekete_chebyshev(m, a, b, n, x, ref nf, ref xf, ref wf);

            Console.WriteLine("");
            Console.WriteLine("  NF = " + nf + "");
            typeMethods.r8vec_print(nf, xf, "  Estimated Fekete points XF:");
            wf_sum = typeMethods.r8vec_sum(nf, wf);
            Console.WriteLine("");
            Console.WriteLine("  Sum(WF) = " + wf_sum + "");
        }

        static void test03(int m)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TEST03 seeks Fekete points in [-1,+1].
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 April 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the dimension of the polynomial space.
            //
        {
            int N = 5001;

            double a;
            double b;
            int n = N;
            int nf = 0;
            double[] wf;
            double wf_sum;
            double[] x;
            double[] xf;

            a = -1.0;
            b = +1.0;
            x = typeMethods.r8vec_linspace_new(n, a, b);

            Console.WriteLine("");
            Console.WriteLine("TEST03:");
            Console.WriteLine("  Seek Fekete points in [" + a + "," + b + "]");
            Console.WriteLine("  using " + n + " equally spaced sample points");
            Console.WriteLine("  for polynomials of degree M = " + m + "");
            Console.WriteLine("  with the Legendre basis and uniform weight.");

            wf = new double[m];
            xf = new double[m];
            LineFekete.line_fekete_legendre(m, a, b, n, x, ref nf, ref xf, ref wf);

            Console.WriteLine("");
            Console.WriteLine("  NF = " + nf + "");
            typeMethods.r8vec_print(nf, xf, "  Estimated Fekete points XF:");
            wf_sum = typeMethods.r8vec_sum(nf, wf);
            Console.WriteLine("");
            Console.WriteLine("  Sum(WF) = " + wf_sum + "");
        }
    }
}