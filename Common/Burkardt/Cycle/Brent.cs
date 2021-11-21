using System;

namespace Burkardt.Cycle;

public static class Brent
{
    public static void cycle_brent( Func < int, int  > f, int x0, ref int lam, ref int mu )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CYCLE_BRENT finds a cycle in an iterated mapping using Brent's method.
        //
        //  Discussion:
        //
        //    Suppose we a repeatedly apply a function f(), starting with the argument
        //    x0, then f(x0), f(f(x0)) and so on.  Suppose that the range of f is finite.
        //    Then eventually the iteration must reach a cycle.  Once the cycle is reached,
        //    succeeding values stay within that cycle.
        //
        //    Starting at x0, there is a "nearest element" of the cycle, which is
        //    reached after MU applications of f.
        //
        //    Once the cycle is entered, the cycle has a length LAM, which is the number
        //    of steps required to first return to a given value.
        //
        //    This function uses Brent's method to determine the values of MU and LAM,
        //    given F and X0.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 June 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Richard Brent,
        //    An improved Monte Carlo factorization algorithm,
        //    BIT,
        //    Volume 20, Number 2, 1980, pages 176-184.
        //
        //  Parameters:
        //
        //    Input, int F ( int i ), the name of the function 
        //    to be analyzed.
        //
        //    Input, int X0, the starting point.
        //
        //    Output, int &LAM, the length of the cycle.
        //
        //    Output, int &MU, the index in the sequence starting
        //    at X0, of the first appearance of an element of the cycle.
        //
    {
        int i;

        int power = 1;
        lam = 1;
        int tortoise = x0;
        int hare = f(x0);

        while (tortoise != hare)
        {
            if (power == lam)
            {
                tortoise = hare;
                power *= 2;
                lam = 0;
            }

            hare = f(hare);
            lam += 1;
        }

        mu = 0;
        tortoise = x0;
        hare = x0;

        for (i = 0; i < lam; i++)
        {
            hare = f(hare);
        }

        while (tortoise != hare)
        {
            tortoise = f(tortoise);
            hare = f(hare);
            mu += 1;
        }
    }
}