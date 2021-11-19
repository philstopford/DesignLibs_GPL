using System;
using Burkardt.AppliedStatistics;
using Burkardt.Types;

namespace ASA144Test;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for ASA144_TEST.
        //
        //  Discussion:
        //
        //    ASA144_TEST tests the ASA144 library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 January 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("ASA144_TEST");
        Console.WriteLine("  Test the ASA144 library.");

        test01();

        Console.WriteLine("");
        Console.WriteLine("ASA144_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests RCONT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 January 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int NROW = 5;
        const int NCOL = 5;

        int ifault = 0;
        int[] matrix = new int[NROW * NCOL];
        int[] ncolt =  {
            2, 2, 2, 2, 1
        };
        int[] nrowt =  {
            3, 2, 2, 1, 1
        };
        int[] nsubt = new int[NCOL];
        int test;
        const int test_num = 10;

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  RCONT constructs a random matrix with");
        Console.WriteLine("  given row and column sums.");

        typeMethods.i4vec_print(NROW, nrowt, "  The rowsum vector:");
        typeMethods.i4vec_print(NCOL, ncolt, "  The columnsum vector: ");

        bool key = false;

        Algorithms.RContData data = new();

        for (test = 1; test <= test_num; test++)
        {
            Algorithms.rcont(ref data, NROW, NCOL, nrowt, ncolt, ref nsubt, ref matrix, ref key, ref ifault);

            if (ifault != 0)
            {
                Console.WriteLine("");
                Console.WriteLine("  RCONT returned IFAULT = " + ifault + "");
                return;
            }

            typeMethods.i4mat_print(NROW, NCOL, matrix, "  The rowcolsum matrix:");
        }
    }
}