using System;

namespace Burkardt.Function
{
    public static class Partition
    {
        public static int plane_partition_num ( int n )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PLANE_PARTITION_NUM returns the number of plane partitions of the integer N.
            //
            //  Discussion:
            //
            //    A plane partition of a positive integer N is a partition of N in which
            //    the parts have been arranged in a 2D array that is nonincreasing across
            //    rows and columns.  There are six plane partitions of 3:
            //
            //      3   2 1   2   1 1 1   1 1   1
            //                1           1     1
            //                                  1
            //
            //  First Values:
            //
            //     N PP(N)
            //     0    1
            //     1    1
            //     2    3
            //     3    6
            //     4   13
            //     5   24
            //     6   48
            //     7   86
            //     8  160
            //     9  282
            //    10  500
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    27 April 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Frank Olver, Daniel Lozier, Ronald Boisvert, Charles Clark,
            //    NIST Handbook of Mathematical Functions,
            //    Cambridge University Press, 2010,
            //    ISBN: 978-0521140638,
            //    LC: QA331.N57.
            //  
            //  Parameters:
            //
            //    Input, int N, the number, which must be at least 0.
            //
            //    Output, int PLANE_PARTITION_NUM, the number of 
            //    plane partitions of N.
            //
        {
            int j;
            int k;
            int nn;
            int[] pp;
            int s2;
            int value;

            if ( n < 0 )
            {
                Console.WriteLine("");
                Console.WriteLine("PLANE_PARTITION_NUM - Fatal error!");
                Console.WriteLine("  0 <= N is required.");
                return ( 1 );
            }

            pp = new int[n+1];
            nn = 0;
            pp[nn] = 1;

            nn = 1;
            if ( nn <= n )
            {
                pp[nn] = 1;
            }

            for ( nn = 2; nn <= n; nn++ )
            {
                pp[nn] = 0;
                for ( j = 1; j <= nn; j++ )
                {
                    s2 = 0;
                    for ( k = 1; k <= j; k++ )
                    {
                        if ( ( j % k ) == 0 )
                        {
                            s2 = s2 + k * k;
                        }
                    }
                    pp[nn] = pp[nn] + pp[nn-j] * s2;
                }
                pp[nn] = pp[nn] / nn;
            }

            value = pp[n];

            return value;
        }

    }
}