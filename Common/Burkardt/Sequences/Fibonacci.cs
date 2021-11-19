using System;
using Burkardt.Types;

namespace Burkardt.Sequence;

public static class Fibonacci
{
    public static int fibonacci_direct(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FIBONACCI_DIRECT computes the N-th Fibonacci number directly.
        //
        //  Discussion:
        //
        //    The formula is:
        //
        //      F(N) = ( PHIP**N - PHIM**N ) / sqrt(5)
        //
        //    where
        //
        //      PHIP = ( 1 + sqrt(5) ) / 2,
        //      PHIM = ( 1 - sqrt(5) ) / 2.
        //
        //  Example:
        //
        //     N   F
        //    --  --
        //     0   0
        //     1   1
        //     2   1
        //     3   2
        //     4   3
        //     5   5
        //     6   8
        //     7  13
        //     8  21
        //     9  34
        //    10  55
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 May 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the index of the Fibonacci number to compute.
        //    N should be nonnegative.
        //
        //    Output, int FIBONACCI_DIRECT, the value of the N-th Fibonacci number.
        //
    {
        const double sqrt5 = 2.236068;
        const double phim = (1.0 - sqrt5) / 2.0;
        const double phip = (1.0 + sqrt5) / 2.0;

        int f = n switch
        {
            < 0 => 0,
            _ => (int) typeMethods.r8_nint((Math.Pow(phip, n) - Math.Pow(phim, n)) / Math.Sqrt(5.0))
        };

        return f;
    }

    public static void fibonacci_floor(int n, ref int f, ref int i)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FIBONACCI_FLOOR returns the largest Fibonacci number less or equal to N.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 May 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the positive integer whose Fibonacci "floor" is desired.
        //
        //    Output, int *F, the largest Fibonacci number less than or equal to N.
        //
        //    Output, int *I, the index of the F.
        //
    {
        switch (n)
        {
            case <= 0:
                i = 0;
                f = 0;
                break;
            default:
            {
                i = (int)(
                    Math.Log(0.5 * (2 * n + 1) * Math.Sqrt(5.0))
                    / Math.Log(0.5 * (1.0 + Math.Sqrt(5.0))));

                f = fibonacci_direct(i);

                if (n < f)
                {
                    i -= 1;
                    f = fibonacci_direct(i);
                }

                break;
            }
        }
    }

    public static void fibonacci_recursive(int n, ref int[] f)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FIBONACCI_RECURSIVE computes the first N Fibonacci numbers.
        //
        //  Discussion:
        //
        //    The 'golden ratio' PHI = (1+sqrt(5))/2 satisfies the equation
        //
        //      X*X-X-1=0
        //
        //    which is often written as:
        //
        //       X        1
        //      --- =  ------
        //       1      X - 1
        //
        //    expressing the fact that a rectangle, whose sides are in proportion X:1,
        //    is similar to the rotated rectangle after a square of side 1 is removed.
        //
        //      <----X---->
        //
        //      +-----*---*
        //      |     |   |  1
        //      |     |   |
        //      +-----*---+
        //      <--1-><X-1>
        //
        //    The formula is:
        //
        //      PHIP = ( 1 + sqrt(5) ) / 2
        //      PHIM = ( 1 - sqrt(5) ) / 2
        //      F(N) = ( PHIP^N + PHIM^N ) / sqrt(5)
        //
        //    Moreover, F(N) can be computed by computing PHIP**N / sqrt(5) and rounding
        //    to the nearest whole number.
        //
        //  First terms:
        //
        //      1
        //      1
        //      2
        //      3
        //      5
        //      8
        //     13
        //     21
        //     34
        //     55
        //     89
        //    144
        //
        //    The 40th number is                  102,334,155.
        //    The 50th number is               12,586,269,025.
        //    The 100th number is 354,224,848,179,261,915,075.
        //
        //  Recursion:
        //
        //    F(1) = 1
        //    F(2) = 1
        //
        //    F(N) = F(N-1) + F(N-2)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 May 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the highest Fibonacci number to compute.
        //
        //    Output, int F[N], the first N Fibonacci numbers.
        //
    {
        int i;

        switch (n)
        {
            case <= 0:
                return;
        }

        f[0] = 1;

        switch (n)
        {
            case <= 1:
                return;
        }

        f[1] = 1;

        for (i = 2; i < n; i++)
        {
            f[i] = f[i - 1] + f[i - 2];
        }
    }
}