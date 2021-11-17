using System;
using Burkardt.Types;

namespace Burkardt.Symbol;

public static class Legendre
{
    public static int legendre_symbol(int q, int p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEGENDRE_SYMBOL evaluates the Legendre symbol (Q/P).
        //
        //  Definition:
        //
        //    Let P be an odd prime.  Q is a QUADRATIC RESIDUE modulo P
        //    if there is an integer R such that R^2 = Q ( mod P ).
        //    The Legendre symbol ( Q / P ) is defined to be:
        //
        //      + 1 if Q ( mod P ) /= 0 and Q is a quadratic residue modulo P,
        //      - 1 if Q ( mod P ) /= 0 and Q is not a quadratic residue modulo P,
        //        0 if Q ( mod P ) == 0.
        //
        //    We can also define ( Q / P ) for P = 2 by:
        //
        //      + 1 if Q ( mod P ) /= 0
        //        0 if Q ( mod P ) == 0
        //
        //  Example:
        //
        //    (0/7) =   0
        //    (1/7) = + 1  ( 1^2 = 1 mod 7 )
        //    (2/7) = + 1  ( 3^2 = 2 mod 7 )
        //    (3/7) = - 1
        //    (4/7) = + 1  ( 2^2 = 4 mod 7 )
        //    (5/7) = - 1
        //    (6/7) = - 1
        //
        //  Discussion:
        //
        //    For any prime P, exactly half of the integers from 1 to P-1
        //    are quadratic residues.
        //
        //    ( 0 / P ) = 0.
        //
        //    ( Q / P ) = ( mod ( Q, P ) / P ).
        //
        //    ( Q / P ) = ( Q1 / P ) * ( Q2 / P ) if Q = Q1 * Q2.
        //
        //    If Q is prime, and P is prime and greater than 2, then:
        //
        //      if ( Q == 1 ) then
        //
        //        ( Q / P ) = 1
        //
        //      else if ( Q == 2 ) then
        //
        //        ( Q / P ) = + 1 if mod ( P, 8 ) = 1 or mod ( P, 8 ) = 7,
        //        ( Q / P ) = - 1 if mod ( P, 8 ) = 3 or mod ( P, 8 ) = 5.
        //
        //      else
        //
        //        ( Q / P ) = - ( P / Q ) if Q = 3 ( mod 4 ) and P = 3 ( mod 4 ),
        //                  =   ( P / Q ) otherwise.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 March 2001
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Charles Pinter,
        //    A Book of Abstract Algebra,
        //    McGraw Hill, 1982, pages 236-237.
        //
        //    Daniel Zwillinger,
        //    CRC Standard Mathematical Tables and Formulae,
        //    30th Edition,
        //    CRC Press, 1996, pages 86-87.
        //
        //  Parameters:
        //
        //    Input, int Q, an integer whose Legendre symbol with
        //    respect to P is desired.
        //
        //    Input, int P, a prime number, greater than 1, with respect
        //    to which the Legendre symbol of Q is desired.
        //
        //    Output, int LEGENDRE_SYMBOL, the Legendre symbol (Q/P).
        //    Ordinarily, this will be -1, 0 or 1.
        //    L = -2, P is less than or equal to 1.
        //    L = -3, P is not prime.
        //    L = -4, the internal stack of factors overflowed.
        //    L = -5, not enough factorization space.
        //
    {
        int FACTOR_MAX = 20;
        int STACK_MAX = 50;

        int[] factor = new int[FACTOR_MAX];
        int i;
        int l;
        int nfactor = 0;
        int nleft = 0;
        int nmore;
        int nstack;
        int[] power = new int[FACTOR_MAX];
        int[] pstack = new int[STACK_MAX];
        int[] qstack = new int[STACK_MAX];
        switch (p)
        {
            //
            //  P must be greater than 1.
            //
            case <= 1:
                Console.WriteLine("");
                Console.WriteLine("LEGENDRE_SYMBOL - Fatal error!");
                Console.WriteLine("  P must be greater than 1.");
                return 1;
        }

        //
        //  P must be prime.
        //
        if (!typeMethods.i4_is_prime(p))
        {
            Console.WriteLine("");
            Console.WriteLine("LEGENDRE_SYMBOL - Fatal error!");
            Console.WriteLine("  P is not prime.");
            return 1;
        }

        switch (q % p)
        {
            //
            //  ( k*P / P ) = 0.
            //
            case 0:
                return 0;
        }

        switch (p)
        {
            //
            //  For the special case P = 2, (Q/P) = 1 for all odd numbers.
            //
            case 2:
                return 1;
        }

        //
        //  Force Q to be nonnegative.
        //
        while (q < 0)
        {
            q += p;
        }

        nstack = 0;
        l = 1;

        for (;;)
        {
            q %= p;
            //
            //  Decompose Q into factors of prime powers.
            //
            typeMethods.i4_factor(q, FACTOR_MAX, ref nfactor, ref factor, ref power, ref nleft);

            if (nleft != 1)
            {
                Console.WriteLine("");
                Console.WriteLine("LEGENDRE_SYMBOL - Fatal error!");
                Console.WriteLine("  Not enough factorization space.");
                return 1;
            }

            //
            //  Each factor which is an odd power is added to the stack.
            //
            nmore = 0;

            for (i = 0; i < nfactor; i++)
            {
                switch (power[i] % 2)
                {
                    case 1:
                    {
                        nmore += 1;

                        if (STACK_MAX <= nstack)
                        {
                            Console.WriteLine("");
                            Console.WriteLine("LEGENDRE_SYMBOL - Fatal error!");
                            Console.WriteLine("  Stack overflow!");
                            return 1;
                        }

                        pstack[nstack] = p;
                        qstack[nstack] = factor[i];
                        nstack += 1;
                        break;
                    }
                }
            }

            if (nmore != 0)
            {
                nstack -= 1;
                q = qstack[nstack];
                switch (q)
                {
                    //
                    //  Check for a Q of 1 or 2.
                    //
                    case 1:
                    case 2 when (p % 8 == 1 || p % 8 == 7):
                        l = +1 * l;
                        break;
                    case 2 when (p % 8 == 3 || p % 8 == 5):
                        l = -1 * l;
                        break;
                    default:
                    {
                        l = (p % 4) switch
                        {
                            3 when q % 4 == 3 => -1 * l,
                            _ => l
                        };

                        typeMethods.i4_swap(ref p, ref q);

                        continue;
                    }
                }

            }

            //
            //  If the stack is empty, we're done.
            //
            if (nstack == 0)
            {
                break;
            }

            //
            //  Otherwise, get the last P and Q from the stack, and process them.
            //
            nstack -= 1;
            p = pstack[nstack];
            q = qstack[nstack];
        }

        return l;
    }

}