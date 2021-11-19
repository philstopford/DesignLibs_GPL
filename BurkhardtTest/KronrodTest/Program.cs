using System;
using Burkardt.Kronrod;

namespace KronrodTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for KRONROD_TEST.
        //
        //  Discussion:
        //
        //    KRONROD_TEST tests KRONROD.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 August 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("KRONROD_TEST:");
        Console.WriteLine("  Test the KRONROD library.");

        test01();
        test02();
        test03();

        Console.WriteLine("");
        Console.WriteLine("KRONROD_TEST:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests the code for the odd case N = 3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 August 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double eps;
        int i;
        int i2;
        int n = 3;
        double s;
        double[] w1;
        double[] w2;
        double[] wg =
        {
            0.555555555555555555556,
            0.888888888888888888889,
            0.555555555555555555556
        };
        double[] wk =
        {
            0.104656226026467265194,
            0.268488089868333440729,
            0.401397414775962222905,
            0.450916538658474142345,
            0.401397414775962222905,
            0.268488089868333440729,
            0.104656226026467265194
        };
        double[] x;
        double[] xg =
        {
            -0.77459666924148337704,
            0.0,
            0.77459666924148337704
        };
        double[] xk =
        {
            -0.96049126870802028342,
            -0.77459666924148337704,
            -0.43424374934680255800,
            0.0,
            0.43424374934680255800,
            0.77459666924148337704,
            0.96049126870802028342
        };

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Request KRONROD to compute the Gauss rule");
        Console.WriteLine("  of order 3, and the Kronrod extension of");
        Console.WriteLine("  order 3+4=7.");
        Console.WriteLine("");
        Console.WriteLine("  Compare to exact data.");

        eps = 0.000001;
        w1 = new double[n + 1];
        w2 = new double[n + 1];
        x = new double[n + 1];

        Quadrature.kronrod(n, eps, ref x, ref w1, ref w2);

        Console.WriteLine("");
        Console.WriteLine("  KRONROD returns 3 vectors of length " + n + 1 + "");
        Console.WriteLine("");
        Console.WriteLine("     I      X               WK              WG");
        Console.WriteLine("");
        for (i = 1; i <= n + 1; i++)
        {
            Console.WriteLine("  " + i.ToString().PadLeft(4)
                                   + "  " + x[i - 1].ToString().PadLeft(14)
                                   + "  " + w1[i - 1].ToString().PadLeft(14)
                                   + "  " + w2[i - 1].ToString().PadLeft(14) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("               Gauss Abscissas");
        Console.WriteLine("            Exact           Computed");
        Console.WriteLine("");
        for (i = 1; i <= n; i++)
        {
            if (2 * i <= n + 1)
            {
                i2 = 2 * i;
                s = -1.0;
            }
            else
            {
                i2 = 2 * (n + 1) - 2 * i;
                s = +1.0;
            }

            Console.WriteLine("  " + i.ToString().PadLeft(4)
                                   + "  " + xg[i - 1].ToString().PadLeft(14)
                                   + "  " + (s * x[i2 - 1]).ToString().PadLeft(14) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("               Gauss Weights");
        Console.WriteLine("            Exact           Computed");
        Console.WriteLine("");
        for (i = 1; i <= n; i++)
        {
            if (2 * i <= n + 1)
            {
                i2 = 2 * i;
            }
            else
            {
                i2 = 2 * (n + 1) - 2 * i;
            }

            Console.WriteLine("  " + i.ToString().PadLeft(4)
                                   + "  " + wg[i - 1].ToString().PadLeft(14)
                                   + "  " + w2[i2 - 1].ToString().PadLeft(14) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("             Gauss Kronrod Abscissas");
        Console.WriteLine("            Exact           Computed");
        Console.WriteLine("");
        for (i = 1; i <= 2 * n + 1; i++)
        {
            if (i <= n + 1)
            {
                i2 = i;
                s = -1.0;
            }
            else
            {
                i2 = 2 * (n + 1) - i;
                s = +1.0;
            }

            Console.WriteLine("  " + i.ToString().PadLeft(4)
                                   + "  " + xk[i - 1].ToString().PadLeft(14)
                                   + "  " + (s * x[i2 - 1]).ToString().PadLeft(14) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("             Gauss Kronrod Weights");
        Console.WriteLine("            Exact           Computed");
        Console.WriteLine("");
        for (i = 1; i <= 2 * n + 1; i++)
        {
            if (i <= n + 1)
            {
                i2 = i;
            }
            else
            {
                i2 = 2 * (n + 1) - i;
            }

            Console.WriteLine("  " + i.ToString().PadLeft(4)
                                   + "  " + wk[i - 1].ToString().PadLeft(14)
                                   + "  " + w1[i2 - 1].ToString().PadLeft(14) + "");
        }
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests the code for the even case N = 4.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 August 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double eps;
        int i;
        int n = 4;
        double[] w1;
        double[] w2;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  Request KRONROD to compute the Gauss rule");
        Console.WriteLine("  of order 4, and the Kronrod extension of");
        Console.WriteLine("  order 4+5=9.");

        eps = 0.000001;
        w1 = new double[n + 1];
        w2 = new double[n + 1];
        x = new double[n + 1];

        Quadrature.kronrod(n, eps, ref x, ref w1, ref w2);

        Console.WriteLine("");
        Console.WriteLine("  KRONROD returns 3 vectors of length " + n + 1 + "");
        Console.WriteLine("");
        Console.WriteLine("     I      X               WK              WG");
        Console.WriteLine("");
        for (i = 1; i <= n + 1; i++)
        {
            Console.WriteLine("  " + i.ToString().PadLeft(4)
                                   + "  " + x[i - 1].ToString().PadLeft(14)
                                   + "  " + w1[i - 1].ToString().PadLeft(14)
                                   + "  " + w2[i - 1].ToString().PadLeft(14) + "");
        }
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 uses the program to estimate an integral.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    24 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double eps;
        double exact = 1.5643964440690497731;
        int i;
        double i1;
        double i2;
        int n;
        double[] w1;
        double[] w2;
        double[] x;

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  Call Kronrod to estimate the integral of a function.");
        Console.WriteLine("  Keep trying until the error is small.");
        //
        //  EPS just tells KRONROD how carefully it must compute X, W1 and W2.
        //  It is NOT a statement about the accuracy of your integral estimate!
        //
        eps = 0.000001;
        //
        //  Start the process with a 1 point rule.
        //
        n = 1;

        for (;;)
        {
            //
            //  Make space.
            //
            w1 = new double[n + 1];
            w2 = new double[n + 1];
            x = new double[n + 1];

            Quadrature.kronrod(n, eps, ref x, ref w1, ref w2);
            //
            //  Compute the estimates.
            //  There are two complications here:
            //
            //  1) Both rules use all the points.  However, the lower order rule uses
            //     a zero weight for the points it doesn't need.
            //
            //  2) The points X are all positive, and are listed in descending order.
            //     this means that 0 is always in the list, and always occurs as the
            //     last member.  Therefore, the integral estimates should use the
            //     function value at 0 once, and the function values at the other
            //     X values "twice", that is, once at X and once at -X.
            //
            i1 = w1[n] * f(x[n]);
            i2 = w2[n] * f(x[n]);

            for (i = 0; i < n; i++)
            {
                i1 += w1[i] * (f(-x[i]) + f(x[i]));
                i2 += w2[i] * (f(-x[i]) + f(x[i]));
            }

            if (Math.Abs(i1 - i2) < 0.0001)
            {
                Console.WriteLine("");
                Console.WriteLine("  Error tolerance satisfied with N = " + n + "");
                Console.WriteLine("  Coarse integral estimate = " + i1.ToString("0.########") + "");
                Console.WriteLine("  Fine   integral estimate = " + i2 + "");
                Console.WriteLine("  Error estimate = " + Math.Abs(i2 - i1) + "");
                Console.WriteLine("  Actual error = " + Math.Abs(exact - i2) + "");
                break;
            }

            if (25 < n)
            {
                Console.WriteLine("");
                Console.WriteLine("  Error tolerance failed even for n = " + n + "");
                Console.WriteLine("  Canceling iteration, and accepting bad estimates!");
                Console.WriteLine("  Coarse integral estimate = " + i1 + "");
                Console.WriteLine("  Fine   integral estimate = " + i2 + "");
                Console.WriteLine("  Error estimate = " + Math.Abs(i2 - i1) + "");
                Console.WriteLine("  Actual error = " + Math.Abs(exact - i2) + "");
                break;
            }

            n = 2 * n + 1;
        }

    }

    private static double f(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F is a function whose integral from -1 to +1 is to be estimated.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    24 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the argument.
        //
        //    Output, double F, the value.
        //
    {
        double value = 0;

        value = 1.0 / (x * x + 1.005);

        return value;
    }
}