using System;
using Burkardt.Types;

namespace Burkardt.Symbol
{
    public class Jacobi
    {
        public static int jacobi_symbol(int q, int p)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    JACOBI_SYMBOL evaluates the Jacobi symbol (Q/P).
            //
            //  Discussion:
            //
            //    If P is prime, then
            //
            //      Jacobi Symbol (Q/P) = Legendre Symbol (Q/P)
            //
            //    Else
            //
            //      let P have the prime factorization
            //
            //        P = Product ( 1 <= I <= N ) P(I)^E(I)
            //
            //      Jacobi Symbol (Q/P) =
            //
            //        Product ( 1 <= I <= N ) Legendre Symbol (Q/P(I))^E(I)
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    30 June 2000
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Daniel Zwillinger,
            //    CRC Standard Mathematical Tables and Formulae,
            //    30th Edition,
            //    CRC Press, 1996, pages 86-87.
            //
            //  Parameters:
            //
            //    Input, int Q, an integer whose Jacobi symbol with
            //    respect to P is desired.
            //
            //    Input, int P, the number with respect to which the Jacobi
            //    symbol of Q is desired.  P should be 2 or greater.
            //
            //    Output, int JACOBI_SYMBOL, the Jacobi symbol (Q/P).
            //    Ordinarily, L will be -1, 0 or 1.
            //    If JACOBI_SYMBOL is -10, an error occurred.
            //
        {
            int FACTOR_MAX = 20;

            int[] factor = new int[FACTOR_MAX];
            int i;
            int l;
            int nfactor = 0;
            int nleft = 0;
            int[] power = new int[FACTOR_MAX];
            int value;
            //
            //  P must be greater than 1.
            //
            if (p <= 1)
            {
                Console.WriteLine("");
                Console.WriteLine("JACOBI_SYMBOL - Fatal error!");
                Console.WriteLine("  P must be greater than 1.");
                return (1);
            }

            //
            //  Decompose P into factors of prime powers.
            //
            typeMethods.i4_factor(p, FACTOR_MAX, ref nfactor, ref factor, ref power, ref nleft);

            if (nleft != 1)
            {
                Console.WriteLine("");
                Console.WriteLine("JACOBI_SYMBOL - Fatal error!");
                Console.WriteLine("  Not enough factorization space.");
                return (1);
            }

            //
            //  Force Q to be nonnegative.
            //
            while (q < 0)
            {
                q = q + p;
            }

            //
            //  For each prime factor, compute the Legendre symbol, and
            //  multiply the Jacobi symbol by the appropriate factor.
            //
            value = 1;

            for (i = 0; i < nfactor; i++)
            {
                l = Legendre.legendre_symbol(q, factor[i]);

                if (l < -1)
                {
                    Console.WriteLine("");
                    Console.WriteLine("JACOBI_SYMBOL - Fatal error!");
                    Console.WriteLine("  Error during Legendre symbol calculation.");
                    return (1);
                }

                value = value * (int)Math.Pow((double)l, power[i]);
            }

            return value;
        }
        
    }
}