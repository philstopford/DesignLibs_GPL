using System;
using Burkardt;
using Burkardt.Function;

namespace GlominTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    glomin_test() tests glomin().
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 June 2021
            //
            //  Author:
            //
            //    John Burkardt
            //
        {
            double a;
            double b;
            double c;
            double e;
            double m;
            double t;

            Console.WriteLine("");
            Console.WriteLine("glomin_test():");
            Console.WriteLine("  glomin() seeks a global minimizer of a function F(X)");
            Console.WriteLine("  in an interval [A,B], given some upper bound M ");
            Console.WriteLine("  for the second derivative of F.");

            e = Math.Sqrt(double.Epsilon);
            t = Math.Sqrt(double.Epsilon);

            Console.WriteLine("");
            Console.WriteLine("  Tolerances:");
            Console.WriteLine("  e = " + e + "");
            Console.WriteLine("  t = " + t + "");

            a = 7.0;
            b = 9.0;
            c = (a + b) / 2.0;
            m = 0.0;
            glomin_example(a, b, c, m, e, t, h_01,
                "h_01(x) = 2 - x");

            a = 7.0;
            b = 9.0;
            c = (a + b) / 2.0;
            m = 100.0;
            glomin_example(a, b, c, m, e, t, h_01,
                "h_01(x) = 2 - x");

            a = -1.0;
            b = +2.0;
            c = (a + b) / 2.0;
            m = 2.0;
            glomin_example(a, b, c, m, e, t, h_02,
                "h_02(x) = x * x");

            a = -1.0;
            b = +2.0;
            c = (a + b) / 2.0;
            m = 2.1;
            glomin_example(a, b, c, m, e, t, h_02,
                "h_02(x) = x * x");

            a = -0.5;
            b = +2.0;
            c = (a + b) / 2.0;
            m = 14.0;
            glomin_example(a, b, c, m, e, t, h_03,
                "h_03(x) = x^3 + x^2");

            a = -0.5;
            b = +2.0;
            c = (a + b) / 2.0;
            m = 28.0;
            glomin_example(a, b, c, m, e, t, h_03,
                "h_03(x) = x^3 + x^2");

            a = -10.0;
            b = +10.0;
            c = (a + b) / 2.0;
            m = 72.0;
            glomin_example(a, b, c, m, e, t, h_04,
                "h_04(x) = ( x + sin(x) ) * exp(-x*x)");

            a = -10.0;
            b = +10.0;
            c = (a + b) / 2.0;
            m = 72.0;
            glomin_example(a, b, c, m, e, t, h_05,
                "h_05(x) = ( x - sin(x) ) * exp(-x*x)");

            Console.WriteLine("");
            Console.WriteLine("glomin_test():");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void glomin_example(double a, double b, double c, double m,
                double e, double t, Func<double, double> f, string title)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    glomin_example() tests glomin() on one test function.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    27 May 2021
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Input:
            //
            //    double A, B, the endpoints of the interval.
            //
            //    double C, an initial guess for the global
            //    minimizer.  If no good guess is known, C = A or B is acceptable.
            //
            //    double M, the bound on the second derivative.
            //
            //    double E, a positive tolerance, a bound for the
            //    absolute error in the evaluation of F(X) for any X in [A,B].
            //
            //    double T, a positive absolute error tolerance.
            //
            //    double F ( double x ), the name of a user-supplied
            //    function whose global minimum is being sought.
            //
            //    string TITLE, a title for the problem.
            //
        {
            int calls = 0;
            double fa;
            double fb;
            double fx;
            double x = 0;

            fx = Glomin.glomin(a, b, c, m, e, t, f, ref x, ref calls);
            fa = f(a);
            fb = f(b);

            Console.WriteLine("");
            Console.WriteLine("  " + title + "");
            Console.WriteLine("  M = " + m + "");
            Console.WriteLine("");
            Console.WriteLine("           A                 X             B");
            Console.WriteLine("         F(A)              F(X)          F(B)");
            Console.WriteLine("  " + a.ToString("0.######").PadLeft(14)
                                   + "  " + x.ToString().PadLeft(14)
                                   + "  " + b.ToString().PadLeft(14) + "");
            Console.WriteLine("  " + fa.ToString().PadLeft(14)
                                   + "  " + fx.ToString().PadLeft(14)
                                   + "  " + fb.ToString().PadLeft(14) + "");
            Console.WriteLine("  Number of calls to F = " + calls + "");

            return;
        }

        static double h_01(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    h_01() evaluates 2 - x.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 April 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Input:
            //
            //    double X, the evaluation point.
            //
            //  Output:
            //
            //    double double H_01, the value of the function at X.
            //
        {
            double value;

            value = 2.0 - x;

            return value;
        }

        static double h_02(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    h_02() evaluates x^2.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 April 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Input:
            //
            //    double X, the evaluation point.
            //
            //  Output:
            //
            //    double H_02, the value of the function at X.
            //
        {
            double value;

            value = x * x;

            return value;
        }

        static double h_03(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    h_03() evaluates x^3+x^2.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 April 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Input:
            //
            //    double X, the evaluation point.
            //
            //  Output:
            //
            //    double H_03, the value of the function at X.
            //
        {
            double value;

            value = x * x * (x + 1.0);

            return value;
        }

        static double h_04(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    h_04() evaluates ( x + sin ( x ) ) * exp ( - x * x ).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 April 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Input:
            //
            //    double X, the evaluation point.
            //
            //  Output:
            //
            //    double H_04, the value of the function at X.
            //
        {
            double value;

            value = (x + Math.Sin(x)) * Math.Exp(-x * x);

            return value;
        }

        static double h_05(double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    h_05() evaluates ( x - sin ( x ) ) * exp ( - x * x ).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 April 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Input:
            //
            //    double X, the evaluation point.
            //
            //  Output:
            //
            //    double H_05, the value of the function at X.
            //
        {
            double value;

            value = (x - Math.Sin(x)) * Math.Exp(-x * x);

            return value;
        }
    }
}