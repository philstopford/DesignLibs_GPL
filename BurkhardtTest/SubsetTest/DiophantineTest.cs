using System;
using Burkardt.SolveNS;

namespace SubsetTestNS;

public static class DiophantineTest
{
    public static void diophantine_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIOPHANTINE_TEST tests DIOPHANTINE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int TEST_NUM = 20;

        int[] a_test =
        {
            1027, 1027, 1027, 1027, -1027,
            -1027, -1027, -1027, 6, 0,
            0, 0, 1, 1, 1,
            1024, 0, 0, 5, 2
        };
        int[] b_test =
        {
            712, 712, -712, -712, 712,
            712, -712, -712, 8, 0,
            1, 1, 0, 0, 1,
            -15625, 0, 3, 0, 4
        };
        int[] c_test =
        {
            7, -7, 7, -7, 7,
            -7, 7, -7, 50, 0,
            0, 1, 0, 1, 0,
            11529, 1, 11, 19, 7
        };

        bool error = false;
        int test_i;
        int x = 0;
        int y = 0;

        Console.WriteLine("");
        Console.WriteLine("DIOPHANTINE_TEST");
        Console.WriteLine("  DIOPHANTINE solves a Diophantine equation:");
        Console.WriteLine("    A * X + B * Y = C");
        Console.WriteLine("");
        Console.WriteLine("        A         B         C         X     Y     Residual");
        Console.WriteLine("");

        for (test_i = 0; test_i < TEST_NUM; test_i++)
        {
            int a = a_test[test_i];
            int b = b_test[test_i];
            int c = c_test[test_i];

            Diophantine.diophantine(a, b, c, ref error, ref x, ref y);

            switch (error)
            {
                case true:
                    Console.WriteLine(a.ToString().PadLeft(10) + "  "
                                                               + b.ToString().PadLeft(10) + "  "
                                                               + c.ToString().PadLeft(10) + "  "
                                                               + "(Error occurred!)" + "");
                    break;
                default:
                    int r = a * x + b * y - c;
                    Console.WriteLine(a.ToString().PadLeft(10) + "  "
                                                               + b.ToString().PadLeft(10) + "  "
                                                               + c.ToString().PadLeft(10) + "  "
                                                               + x.ToString().PadLeft(10) + "  "
                                                               + y.ToString().PadLeft(10) + "  "
                                                               + r.ToString().PadLeft(10) + "");
                    break;
            }

        }
    }

    public static void diophantine_solution_minimize_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIOPHANTINE_SOLUTION_MINIMIZE_TEST tests DIOPHANTINE_SOLUTION_MINIMIZE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int a = 4096;
        const int b = -15625;
        const int c = 46116;

        Console.WriteLine("");
        Console.WriteLine("DIOPHANTINE_SOLUTION_MINIMIZE_TEST");
        Console.WriteLine("  DIOPHANTINE_SOLUTION_MINIMIZE computes a minimal");
        Console.WriteLine("  Euclidean norm solution of a Diophantine equation:");
        Console.WriteLine("    A * X + B * Y = C");

        int x = 665499996;
        int y = 174456828;

        int r = a * x + b * y - c;

        Console.WriteLine("");
        Console.WriteLine("  Coefficients:");
        Console.WriteLine("    A = " + a.ToString().PadLeft(12) + "");
        Console.WriteLine("    B = " + b.ToString().PadLeft(12) + "");
        Console.WriteLine("    C = " + c.ToString().PadLeft(12) + "");
        Console.WriteLine("  Solution:");
        Console.WriteLine("    X = " + x.ToString().PadLeft(12) + "");
        Console.WriteLine("    Y = " + y.ToString().PadLeft(12) + "");
        Console.WriteLine("  Residual R = A * X + B * Y - C:");
        Console.WriteLine("    R = " + r.ToString().PadLeft(12) + "");

        Diophantine.diophantine_solution_minimize(a, b, ref x, ref y);

        r = a * x + b * y - c;

        Console.WriteLine("");
        Console.WriteLine("  DIOPHANTINE_SOLUTION_MINIMIZE returns");
        Console.WriteLine("  the minimized solution:");
        Console.WriteLine("    X = " + x.ToString().PadLeft(12) + "");
        Console.WriteLine("    Y = " + y.ToString().PadLeft(12) + "");
        Console.WriteLine("  Residual R = A * X + B * Y - C:");
        Console.WriteLine("    R = " + r.ToString().PadLeft(12) + "");

        x = 15621;
        y = 4092;

        r = a * x + b * y - c;

        Console.WriteLine("");
        Console.WriteLine("  Here is the minimal positive solution:");
        Console.WriteLine("    X = " + x.ToString().PadLeft(12) + "");
        Console.WriteLine("    Y = " + y.ToString().PadLeft(12) + "");
        Console.WriteLine("  Residual R = A * X + B * Y - C:");
        Console.WriteLine("    R = " + r.ToString().PadLeft(12) + "");

    }

}