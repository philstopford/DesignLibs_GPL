using System;
using Burkardt.Types;

namespace SearchTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for SEARCH_SERIAL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int a = 1;
        int b = typeMethods.i4_huge();
        const int c = 45;

        Console.WriteLine("");
        Console.WriteLine("SEARCH_SERIAL:");
        Console.WriteLine("  Search the integers from A to B");
        Console.WriteLine("  for a value J such that F(J) = C.");
        Console.WriteLine("");
        Console.WriteLine("  A           = " + a + "");
        Console.WriteLine("  B           = " + b + "");
        Console.WriteLine("  C           = " + c + "");

        DateTime wtime = DateTime.Now;

        int j = search(a, b, c);

        DateTime wtime2 = DateTime.Now;

        switch (j)
        {
            case -1:
                Console.WriteLine("");
                Console.WriteLine("  No solution was found.");
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("  Found     J = " + j + "");
                Console.WriteLine("  Verify F(J) = " + f(j) + "");
                break;
        }

        Console.WriteLine("  Elapsed CPU time is " + (wtime2 - wtime).TotalSeconds + "");

        Console.WriteLine("");
        Console.WriteLine("SEARCH_SERIAL:");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static int f(int i)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F is the function we are analyzing.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int I, the argument.
        //
        //    Input, int F, the value.
        //
    {
        int j;

        int value = i;

        for (j = 1; j <= 5; j++)
        {
            int k = value / 127773;

            value = 16807 * (value - k * 127773) - k * 2836;

            switch (value)
            {
                case <= 0:
                    value += typeMethods.i4_huge();
                    break;
            }
        }

        return value;
    }

    private static int search(int a, int b, int c)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SEARCH searches integers in [A,B] for a J so that F(J) = C.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int A, B, the search range.
        //
        //    Input, int C, the desired function value.
        //
        //    Output, int SEARCH, the computed solution, or -1
        //    if no solution was found.
        //
    {
        int i;

        int j = -1;

        for (i = a; i <= b; i++)
        {
            int fi = f(i);

            if (fi != c)
            {
                continue;
            }

            j = i;
            break;
        }

        return j;
    }
}