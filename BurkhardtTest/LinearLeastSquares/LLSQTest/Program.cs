using System;
using Burkardt.SolveNS;

namespace LLSQTest;

internal class Program
{
    private static void Main(string[] args)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for LLSQ_TEST.
        //
        //  Discussion:
        //
        //    LLSQ_TEST tests the LLSQ library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 January 2019
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("LLSQ_TEST");
        Console.WriteLine("  Test the LLSQ library.");

        test01();
        test02();
        //
        //  Terminate.
        //
        Console.WriteLine("");
        Console.WriteLine("LLSQ_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 calls LLSQ to match 15 data values.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 July 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a = 0;
        double b = 0;
        double error;
        int i;
        int n = 15;
        double[] x =
        {
            1.47, 1.50, 1.52, 1.55, 1.57, 1.60, 1.63, 1.65, 1.68, 1.70,
            1.73, 1.75, 1.78, 1.80, 1.83
        };
        double[] y =
        {
            52.21, 53.12, 54.48, 55.84, 57.20, 58.57, 59.93, 61.29, 63.11, 64.47,
            66.28, 68.10, 69.92, 72.19, 74.46
        };

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  LLSQ can compute the formula for a line of the form");
        Console.WriteLine("    y = A * x + B");
        Console.WriteLine("  which minimizes the RMS error to a set of N data values.");

        LinearLeastSquares.llsq(n, x, y, ref a, ref b);

        Console.WriteLine("");
        Console.WriteLine("  Estimated relationship is y = " + a + " * x + " + b + "");
        Console.WriteLine("  Expected value is         y = 61.272 * x - 39.062");
        Console.WriteLine("");
        Console.WriteLine("     I      X       Y      B+A*X    |error|");
        Console.WriteLine("");
        error = 0.0;
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString().PadLeft(4)
                                   + "  " + x[i].ToString().PadLeft(7)
                                   + "  " + y[i].ToString().PadLeft(7)
                                   + "  " + (b + a * x[i]).ToString().PadLeft(7)
                                   + "  " + (b + a * x[i] - y[i]).ToString().PadLeft(7) + "");
            error += Math.Pow(b + a * x[i] - y[i], 2);
        }

        error = Math.Sqrt(error / n);
        Console.WriteLine("");
        Console.WriteLine("  RMS error =                      " + error + "");

    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 calls LLSQ0 to match 14 data values.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 January 2019
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        double a = 0;
        double error;
        int i;
        int n = 14;
        double[] x =
        {
            0.00, 0.10, 0.15, 0.20, 0.25,
            0.30, 0.35, 0.40, 0.45, 0.50,
            0.55, 0.60, 0.65, 0.70
        };
        double[] y =
        {
            0.0000, 0.0865, 0.1015, 0.1106, 0.1279,
            0.1892, 0.2695, 0.2888, 0.2425, 0.3465,
            0.3225, 0.3764, 0.4263, 0.4562
        };

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  LLSQ0 can compute the formula for a line of the form");
        Console.WriteLine("    y = A * x");
        Console.WriteLine("  which minimizes the RMS error to a set of N data values.");

        LinearLeastSquares.llsq0(n, x, y, ref a);

        Console.WriteLine("");
        Console.WriteLine("  Estimated relationship is y = " + a + " * x");
        Console.WriteLine("");
        Console.WriteLine("     I      X       Y        A*X    |error|");
        Console.WriteLine("");
        error = 0.0;
        for (i = 0; i < n; i++)
        {
            Console.WriteLine("  " + i.ToString().PadLeft(4)
                                   + "  " + x[i]
                                   + "  " + y[i]
                                   + "  " + a * x[i]
                                   + "  " + (a * x[i] - y[i]) + "");
            error += Math.Pow(a * x[i] - y[i], 2);
        }

        error = Math.Sqrt(error / n);
        Console.WriteLine("");
        Console.WriteLine("  RMS error =                      " + error + "");
    }
}