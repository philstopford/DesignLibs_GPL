using System;
using Burkardt.Types;

namespace Burkardt
{
    public static class Josephus
    {
        public static int josephus(int n, int m, int k)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    JOSEPHUS returns the position X of the K-th man to be executed.
            //
            //  Discussion:
            //
            //    The classic Josephus problem concerns a circle of 41 men.
            //    Every third man is killed and removed from the circle.  Counting
            //    and executing continues until all are dead.  Where was the last
            //    survivor sitting?
            //
            //    Note that the first person killed was sitting in the third position.
            //    Moreover, when we get down to 2 people, and we need to count the
            //    "third" one, we just do the obvious thing, which is to keep counting
            //    around the circle until our count is completed.
            //
            //    The process may be regarded as generating a permutation of
            //    the integers from 1 to N.  The permutation would be the execution
            //    list, that is, the list of the executed men, by position number.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    29 May 2003
            //
            //  Author:
            //
            //     John Burkardt.
            //
            //  Reference:
            //
            //    W W Rouse Ball,
            //    Mathematical Recreations and Essays,
            //    Macmillan, 1962, pages 32-36.
            //
            //    Donald Knuth,
            //    The Art of Computer Programming,
            //    Volume 1, Fundamental Algorithms,
            //    Addison Wesley, 1968, pages 158-159.
            //
            //    Donald Knuth,
            //    The Art of Computer Programming,
            //    Volume 3, Sorting and Searching,
            //    Addison Wesley, 1968, pages 18-19.
            //
            //  Parameters:
            //
            //    Input, int N, the number of men.
            //    N must be positive.
            //
            //    Input, int M, the counting index.
            //    M must not be zero.  Ordinarily, M is positive, and no greater than N.
            //
            //    Input, int K, the index of the executed man of interest.
            //    K must be between 1 and N.
            //
            //    Output, int JOSEPHUS, the position of the K-th man.
            //    The value will be between 1 and N.
            //
        {
            int m2;
            int x;

            if (n <= 0)
            {
                Console.WriteLine("");
                Console.WriteLine("JOSEPHUS - Fatal error!");
                Console.WriteLine("  N <= 0.");
                return (1);
            }

            if (m == 0)
            {
                Console.WriteLine("");
                Console.WriteLine("JOSEPHUS - Fatal error!");
                Console.WriteLine("  M = 0.");
                return (1);
            }

            if (k <= 0 || n < k)
            {
                Console.WriteLine("");
                Console.WriteLine("JOSEPHUS - Fatal error!");
                Console.WriteLine("  J <= 0 or N < K.");
                return (1);
            }

            //
            //  In case M is bigger than N, or negative, get the
            //  equivalent positive value between 1 and N.
            //  You can skip this operation if 1 <= M <= N.
            //
            m2 = typeMethods.i4_modp(m, n);

            x = k * m2;

            while (n < x)
            {
                x = (m2 * (x - n) - 1) / (m2 - 1);
            }

            return x;
        }
    }
}