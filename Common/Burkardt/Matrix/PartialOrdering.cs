using System;

namespace Burkardt.MatrixNS;

public static class PartialOrdering
{
    public static bool pord_check ( int n, int[] a )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PORD_CHECK checks a matrix representing a partial ordering.
        //
        //  Discussion:
        //
        //    The array A is supposed to represent a partial ordering of
        //    the elements of a set of N objects.
        //
        //    For distinct indices I and J, the value of A(I,J) is:
        //
        //      1, if I + J
        //      0, otherwise ( I and J may be unrelated, or perhaps J + I).
        //
        //    Diagonal elements of A are ignored.
        //
        //    This routine checks that the values of A do represent
        //    a partial ordering.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 May 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of elements in the set.
        //
        //    Input, int A[N*N], the partial ordering.  A[I+J*N] is
        //    1 if I is less than J in the partial ordering,
        //    0 otherwise.
        //
        //    Output, bool PORD_CHECK, is true if an error was detected.
        //
    {
        int i;

        switch (n)
        {
            case <= 0:
                Console.WriteLine("");
                Console.WriteLine("PORD_CHECK - Fatal error!");
                Console.WriteLine("  N must be positive, but N = " + n + "");
                return true;
        }

        for ( i = 0; i < n; i++ )
        {
            int j;
            for ( j = i + 1; j < n; j++ )
            {
                switch (a[i+j*n])
                {
                    case > 0 when 0 < a[j+i*n]:
                        Console.WriteLine("");
                        Console.WriteLine("PORD_CHECK - Fatal error!");
                        Console.WriteLine("  For indices I = " + i + "");
                        Console.WriteLine("  and J = " + j + "");
                        Console.WriteLine("  A(I,J) = " + a[i+j*n] + "");
                        Console.WriteLine("  A(J,I) = " + a[j+i*n] + "");
                        return true;
                }
            }
        }

        return false;
    }
}