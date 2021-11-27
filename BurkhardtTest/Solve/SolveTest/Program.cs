using System;
using Burkardt.SolveNS;
using Burkardt.Types;

namespace SolveTest;

internal static class Program
{
    private static void Main()
//****************************************************************************80
//
//  Discussion:
//
//    MAIN is the main program for SOLVE_TEST.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 May 2014
//
//  Author:
//
//    John Burkardt
//
    {
        Console.WriteLine("");
        Console.WriteLine("SOLVE_TEST");
        Console.WriteLine("  Test the SOLVE library.");

        test01();

        Console.WriteLine("");
        Console.WriteLine("SOLVE_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

//****************************************************************************80
//
//  Purpose:
//
//    TEST01 demonstrates how a 3X3 linear system can be set up and solved.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    27 May 2014
//
//  Author:
//
//    John Burkardt
//
    {
        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Set up a linear system, and solve it by calling a function.");
        Console.WriteLine("  The linear system is to be accessed using the A[i][j] notation.");
        Console.WriteLine("  It is usually difficult to use this notation and still be able");
        Console.WriteLine("  to pass the array A to a function.  An ordinary doubly indexed");
        Console.WriteLine("  array would require the receiving function to know, IN ADVANCE,");
        Console.WriteLine("  the exact fixed value of the second dimension of A, which defeats");
        Console.WriteLine("  the goal of writing general usable library software.");
        Console.WriteLine("");
        Console.WriteLine("  Here, I think I have made it easy, at the cost of:");
        Console.WriteLine("  * declaring the array as a 'double **a'");
        Console.WriteLine("  * creating the array with r8rmat_new() or r8rmat_zero()");
        Console.WriteLine("  * solving the system with r8rmat_fs_new()");
        Console.WriteLine("  * deleting the array with r8rmat_delete()");
//
//  Define the array size.
//
        int n = 3;
//
//  Create the array that will contain the matrix.
//
        double[][] a = typeMethods.r8rmat_new(n, n);
//
//  Set the array values.
//
        a[0][0] = 1;
        a[0][1] = 2;
        a[0][2] = 3;

        a[1][0] = 4;
        a[1][1] = 5;
        a[1][2] = 6;

        a[2][0] = 7;
        a[2][1] = 8;
        a[2][2] = 0;
//
//  Create the right hand side.
//
        double[] b = new double[n];
//
//  Set the right hand side values.
//
        b[0] = 14;
        b[1] = 32;
        b[2] = 23;
//
//  Request the solution of A*x=b.
//
        double[] x = Solve.r8rmat_fs_new(n, a, b);

        typeMethods.r8vec_print(n, x, "  Solution:");
    }
}