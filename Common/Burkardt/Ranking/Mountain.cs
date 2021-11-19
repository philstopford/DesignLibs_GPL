using System;
using Burkardt.Types;

namespace Burkardt.RankingNS;

public static partial class Ranking
{
    public static int mountain ( int n, int x, int y )

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    MOUNTAIN enumerates the mountains.
        // 
        //  Licensing:
        // 
        //    This code is distributed under the GNU LGPL license.
        // 
        //  Modified:
        // 
        //    24 July 2011
        // 
        //  Author:
        // 
        //    John Burkardt
        // 
        //  Reference:
        // 
        //    Donald Kreher, Douglas Simpson,
        //    Combinatorial Algorithms,
        //    CRC Press, 1998,
        //    ISBN: 0-8493-3988-X,
        //    LC: QA164.K73.
        // 
        //  Parameters:
        // 
        //    Input, int N, ...
        //    N must be positive.
        // 
        //    Input, int X, Y, ...
        //    0 <= X <= 2 * N,
        // 
        //    Output, int MOUNTAIN, the value of the "mountain function"
        //    M ( N, X, Y ), which is the number of all mountain ranges from
        //    (X,Y) to (2*N,0) which do not drop below sea level.
        // 
    {
        int value;
        switch (n)
        {
            // 
            //  Check.
            // 
            case <= 0:
                Console.WriteLine("");
                Console.WriteLine("MOUNTAIN - Fatal error!");
                Console.WriteLine("  N <= 0.");
                Console.WriteLine("  N = " + n + "");
                return 1;
        }

        switch (x)
        {
            case < 0:
                Console.WriteLine("");
                Console.WriteLine("MOUNTAIN - Fatal error!");
                Console.WriteLine("  X < 0.");
                Console.WriteLine("  X = " + x + "");
                return 1;
        }
        if ( 2 * n < x )
        {
            Console.WriteLine("");
            Console.WriteLine("MOUNTAIN - Fatal error!");
            Console.WriteLine("  2 * N < X.");
            Console.WriteLine("  X = " + x + "");
            Console.WriteLine("  N = " + n + "");
            return 1;
        }

        switch (y)
        {
            // 
            //  Special cases.
            // 
            case < 0:
                value = 0;
                break;
            default:
            {
                if ( 2 * n < x + y )
                {
                    value = 0;
                }
                else
                {
                    switch (( x + y ) % 2)
                    {
                        case 1:
                            value = 0;
                            break;
                        default:
                            int a = 2 * n - x;
                            int b = n - ( x + y ) / 2;
                            int c = n - 1 - ( x + y ) / 2;
                            value = typeMethods.i4_choose ( a, b ) - typeMethods.i4_choose ( a, c );
                            break;
                    }
                }

                break;
            }
        }
        return value;
    }
}