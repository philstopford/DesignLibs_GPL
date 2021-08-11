using System;
using Burkardt.RandomNS;
using Burkardt.SubsetNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt
{
    public class CompNZData
    {
        public int h = 0;
        public int t = 0;
    }

    public static class Comp
    {
        public static int comp_enum(int n, int k)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    COMP_ENUM returns the number of compositions of the integer N into K parts.
            //
            //  Discussion:
            //
            //    A composition of the integer N into K parts is an ordered sequence
            //    of K nonnegative integers which sum to N.  The compositions (1,2,1)
            //    and (1,1,2) are considered to be distinct.
            //
            //    The 28 compositions of 6 into three parts are:
            //
            //      6 0 0,  5 1 0,  5 0 1,  4 2 0,  4 1 1,  4 0 2,
            //      3 3 0,  3 2 1,  3 1 2,  3 0 3,  2 4 0,  2 3 1,
            //      2 2 2,  2 1 3,  2 0 4,  1 5 0,  1 4 1,  1 3 2,
            //      1 2 3,  1 1 4,  1 0 5,  0 6 0,  0 5 1,  0 4 2,
            //      0 3 3,  0 2 4,  0 1 5,  0 0 6.
            //
            //    The formula for the number of compositions of N into K parts is
            //
            //      Number = ( N + K - 1 )! / ( N! * ( K - 1 )! )
            //
            //    Describe the composition using N '1's and K-1 dividing lines '|'.
            //    The number of distinct permutations of these symbols is the number
            //    of compositions.  This is equal to the number of permutations of
            //    N+K-1 things, with N identical of one kind and K-1 identical of another.
            //
            //    Thus, for the above example, we have:
            //
            //      Number = ( 6 + 3 - 1 )! / ( 6! * (3-1)! ) = 28
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 December 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Albert Nijenhuis, Herbert Wilf,
            //    Combinatorial Algorithms for Computers and Calculators,
            //    Second Edition,
            //    Academic Press, 1978,
            //    ISBN: 0-12-519260-6,
            //    LC: QA164.N54.
            //
            //  Parameters:
            //
            //    Input, int N, the integer whose compositions are desired.
            //
            //    Input, int K, the number of parts in the composition.
            //
            //    Output, int COMP_ENUM, the number of compositions of N
            //    into K parts.
            //
        {
            int number;

            number = typeMethods.i4_choose(n + k - 1, n);

            return number;
        }

        public static void comp_next_grlex(int kc, ref int[] xc)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    COMP_NEXT_GRLEX returns the next composition in grlex order.
            //
            //  Discussion:
            //
            //    Example:
            //
            //    KC = 3
            //
            //    #   XC(1) XC(2) XC(3)  Degree
            //      +------------------------
            //    1 |  0     0     0        0
            //      |
            //    2 |  0     0     1        1
            //    3 |  0     1     0        1
            //    4 |  1     0     0        1
            //      |
            //    5 |  0     0     2        2
            //    6 |  0     1     1        2
            //    7 |  0     2     0        2
            //    8 |  1     0     1        2
            //    9 |  1     1     0        2
            //   10 |  2     0     0        2
            //      |
            //   11 |  0     0     3        3
            //   12 |  0     1     2        3
            //   13 |  0     2     1        3
            //   14 |  0     3     0        3
            //   15 |  1     0     2        3
            //   16 |  1     1     1        3
            //   17 |  1     2     0        3
            //   18 |  2     0     1        3
            //   19 |  2     1     0        3
            //   20 |  3     0     0        3
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 December 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int KC, the number of parts of the composition.
            //    1 <= KC.
            //
            //    Input/output, int XC[KC], the current composition.
            //    Each entry of XC must be nonnegative.
            //    On return, XC has been replaced by the next composition in the
            //    grlex order.
            //
        {
            int i;
            int im1 = 0;
            int j;
            int t = 0;
            //
            //  Ensure that 1 <= KC.
            //
            if (kc < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("COMP_NEXT_GRLEX - Fatal error!");
                Console.WriteLine("  KC < 1");
                return;
            }

            //
            //  Ensure that 0 <= XC(I).
            //
            for (i = 0; i < kc; i++)
            {
                if (xc[i] < 0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("COMP_NEXT_GRLEX - Fatal error!");
                    Console.WriteLine("  XC[I] < 0");
                    return;
                }
            }

            //
            //  Find I, the index of the rightmost nonzero entry of X.
            //
            i = 0;
            for (j = kc; 1 <= j; j--)
            {
                if (0 < xc[j - 1])
                {
                    i = j;
                    break;
                }
            }

            //
            //  set T = X(I)
            //  set XC(I) to zero,
            //  increase XC(I-1) by 1,
            //  increment XC(KC) by T-1.
            //
            if (i == 0)
            {
                xc[kc - 1] = 1;
                return;
            }
            else if (i == 1)
            {
                t = xc[0] + 1;
                im1 = kc;
            }
            else if (1 < i)
            {
                t = xc[i - 1];
                im1 = i - 1;
            }

            xc[i - 1] = 0;
            xc[im1 - 1] = xc[im1 - 1] + 1;
            xc[kc - 1] = xc[kc - 1] + t - 1;
        }

        public static int[] comp_random_grlex(int kc, int rank1, int rank2, ref int seed, ref int rank )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COMP_RANDOM_GRLEX: random composition with degree less than or equal to NC.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 September 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int KC, the number of parts in the composition.
        //
        //    Input, int RANK1, RANK2, the minimum and maximum ranks.
        //    1 <= RANK1 <= RANK2.
        //
        //    Input/output, int &SEED, the random number seed.
        //
        //    Output, int &RANK, the rank of the composition.
        //
        //    Output, int COMP_RANDOM_GRLEX[KC], the random composition.
        //
        {
            int[] xc;

            //
            //  Ensure that 1 <= KC.
            //
            if (kc < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("COMP_RANDOM_GRLEX - Fatal error!");
                Console.WriteLine("  KC < 1");
                return null;
            }

            //
            //  Ensure that 1 <= RANK1.
            //
            if (rank1 < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("COMP_RANDOM_GRLEX - Fatal error!");
                Console.WriteLine("  RANK1 < 1");
                return null;
            }

            //
            //  Ensure that RANK1 <= RANK2.
            //
            if (rank2 < rank1)
            {
                Console.WriteLine("");
                Console.WriteLine("COMP_RANDOM_GRLEX - Fatal error!");
                Console.WriteLine("  RANK2 < RANK1");
                return null;
            }

            //
            //  Choose RANK between RANK1 and RANK2.
            //
            rank = UniformRNG.i4_uniform_ab(rank1, rank2, ref seed);
            //
            //  Recover the composition of given RANK.
            //
            xc = comp_unrank_grlex(kc, rank);

            return xc;
        }

        public static int comp_rank_grlex(int kc, int[] xc)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    COMP_RANK_GRLEX computes the graded lexicographic rank of a composition.
            //
            //  Discussion:
            //
            //    The graded lexicographic ordering is used, over all KC-compositions
            //    for NC = 0, 1, 2, ...
            //
            //    For example, if KC = 3, the ranking begins:
            //
            //    Rank  Sum    1  2  3
            //    ----  ---   -- -- --
            //       1    0    0  0  0
            //
            //       2    1    0  0  1
            //       3    1    0  1  0
            //       4    1    1  0  1
            //
            //       5    2    0  0  2
            //       6    2    0  1  1
            //       7    2    0  2  0
            //       8    2    1  0  1
            //       9    2    1  1  0
            //      10    2    2  0  0
            //
            //      11    3    0  0  3
            //      12    3    0  1  2
            //      13    3    0  2  1
            //      14    3    0  3  0
            //      15    3    1  0  2
            //      16    3    1  1  1
            //      17    3    1  2  0
            //      18    3    2  0  1
            //      19    3    2  1  0
            //      20    3    3  0  0
            //
            //      21    4    0  0  4
            //      ..   ..   .. .. ..
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 December 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int KC, the number of parts in the composition.
            //    1 <= KC.
            //
            //    Input, int XC[KC], the composition.
            //    For each 1 <= I <= KC, we have 0 <= XC(I).
            //
            //    Output, int COMP_RANK_GRLEX, the rank of the composition.
            //
        {
            int i;
            int j;
            int ks;
            int n;
            int nc;
            int ns;
            int rank;
            int tim1;
            int[] xs;
            //
            //  Ensure that 1 <= KC.
            //
            if (kc < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("COMP_RANK_GRLEX - Fatal error!");
                Console.WriteLine("  KC < 1");
                return (1);
            }

            //
            //  Ensure that 0 <= XC(I).
            //
            for (i = 0; i < kc; i++)
            {
                if (xc[i] < 0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("COMP_RANK_GRLEX - Fatal error!");
                    Console.WriteLine("  XC[I] < 0");
                    return (1);
                }
            }

            //
            //  NC = sum ( XC )
            //
            nc = typeMethods.i4vec_sum(kc, xc);
            //
            //  Convert to KSUBSET format.
            //
            ns = nc + kc - 1;
            ks = kc - 1;
            xs = new int[ks];
            xs[0] = xc[0] + 1;
            for (i = 2; i < kc; i++)
            {
                xs[i - 1] = xs[i - 2] + xc[i - 1] + 1;
            }

            //
            //  Compute the rank.
            //
            rank = 1;

            for (i = 1; i <= ks; i++)
            {
                if (i == 1)
                {
                    tim1 = 0;
                }
                else
                {
                    tim1 = xs[i - 2];
                }

                if (tim1 + 1 <= xs[i - 1] - 1)
                {
                    for (j = tim1 + 1; j <= xs[i - 1] - 1; j++)
                    {
                        rank = rank + typeMethods.i4_choose(ns - j, ks - i);
                    }
                }
            }

            for (n = 0; n < nc; n++)
            {
                rank = rank + typeMethods.i4_choose(n + kc - 1, n);
            }

            return rank;
        }
        
        public static void comp_to_ksub ( int nc, int kc, int[] ac, ref int ns, ref int ks, ref int[] as_ )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COMP_TO_KSUB converts a composition to a K-subset.
        //
        //  Discussion:
        //
        //    There is a bijection between K subsets and compositions.
        //
        //    Because we allow a composition to have entries that are 0, we need
        //    to implicitly add 1 to each entry before establishing the bijection.
        //
        //    Let AC be a composition of NC into KC parts.
        //
        //    Then let
        //      NS = NC + KC - 1
        //      KS = KC - 1
        //    and define
        //      AS(I) = sum ( AC(1:I) + 1 ), for I = 1 : KS.
        //      
        //    Then AS is a KS subset of the integers 1 through NS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 December 2013
        //
        //  Parameters:
        //
        //    Input, int NC, the composition sum.
        //
        //    Input, int KC, the number of parts of the composition.
        //
        //    Input, int AC[KC], the parts of the composition.
        //
        //    Output, int &NS, the size of the set.
        //
        //    Output, int &KS, the size of the subset.
        //
        //    Output, int AS[KS], the entries of the K-subset, 
        //    in increasing order.
        //
        {
            int i;

            ns = nc + kc - 1;
            ks = kc - 1;
                as_[0] = ac[0] + 1;
            for ( i = 1; i < kc - 1; i++ )
            {
                as_[i] = as_[i-1] + ac[i] + 1;
            }
        }

        public static int[] comp_unrank_grlex(int kc, int rank)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    COMP_UNRANK_GRLEX computes the composition of given grlex rank.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 December 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int KC, the number of parts of the composition.
            //    1 <= KC.
            //
            //    Input, int RANK, the rank of the composition.
            //    1 <= RANK.
            //
            //    Output, int COMP_UNRANK_GRLEX[KC], the composition XC of the given rank.
            //    For each I, 0 <= XC[I] <= NC, and 
            //    sum ( 1 <= I <= KC ) XC[I] = NC.
            //
        {
            int i;
            int j;
            int ks;
            int nc;
            int ns;
            int r;
            int rank1;
            int rank2;
            int[] xc;
            int[] xs;
            //
            //  Ensure that 1 <= KC.
            //
            if (kc < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("COMP_UNRANK_GRLEX - Fatal error!");
                Console.WriteLine("  KC < 1");
                return null;
            }

            //
            //  Ensure that 1 <= RANK.
            //
            if (rank < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("COMP_UNRANK_GRLEX - Fatal error!");
                Console.WriteLine("  RANK < 1");
                return null;
            }

            //
            //  Determine the appropriate value of NC.
            //  Do this by adding up the number of compositions of sum 0, 1, 2, 
            //  ..., without exceeding RANK.  Moreover, RANK - this sum essentially
            //  gives you the rank of the composition within the set of compositions
            //  of sum NC.  And that's the number you need in order to do the
            //  unranking.
            //
            rank1 = 1;
            nc = -1;
            for (;;)
            {
                nc = nc + 1;
                r = typeMethods.i4_choose(nc + kc - 1, nc);
                if (rank < rank1 + r)
                {
                    break;
                }

                rank1 = rank1 + r;
            }

            rank2 = rank - rank1;
            //
            //  Convert to KSUBSET format.
            //  Apology: an unranking algorithm was available for KSUBSETS,
            //  but not immediately for compositions.  One day we will come back
            //  and simplify all this.
            //
            ks = kc - 1;
            ns = nc + kc - 1;
            xs = new int[ks];

            j = 1;

            for (i = 1; i <= ks; i++)
            {
                r = typeMethods.i4_choose(ns - j, ks - i);

                while (r <= rank2 && 0 < r)
                {
                    rank2 = rank2 - r;
                    j = j + 1;
                    r = typeMethods.i4_choose(ns - j, ks - i);
                }

                xs[i - 1] = j;
                j = j + 1;
            }

            //
            //  Convert from KSUBSET format to COMP format.
            //
            xc = new int[kc];
            xc[0] = xs[0] - 1;
            for (i = 2; i < kc; i++)
            {
                xc[i - 1] = xs[i - 1] - xs[i - 2] - 1;
            }

            xc[kc - 1] = ns - xs[ks - 1];

            return xc;
        }

        public static void comp_next(int n, int k, ref int[] a, ref bool more, ref int h, ref int t)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    COMP_NEXT computes the compositions of the integer N into K parts.
            //
            //  Discussion:
            //
            //    A composition of the integer N into K parts is an ordered sequence
            //    of K nonnegative integers which sum to N.  The compositions (1,2,1)
            //    and (1,1,2) are considered to be distinct.
            //
            //    The routine computes one composition on each call until there are no more.
            //    For instance, one composition of 6 into 3 parts is
            //    3+2+1, another would be 6+0+0.
            //
            //    On the first call to this routine, set MORE = FALSE.  The routine
            //    will compute the first element in the sequence of compositions, and
            //    return it, as well as setting MORE = TRUE.  If more compositions
            //    are desired, call again, and again.  Each time, the routine will
            //    return with a new composition.
            //
            //    However, when the LAST composition in the sequence is computed 
            //    and returned, the routine will reset MORE to FALSE, signaling that
            //    the end of the sequence has been reached.
            //
            //    This routine originally used a SAVE statement to maintain the
            //    variables H and T.  I have decided that it is safer
            //    to pass these variables as arguments, even though the user should
            //    never alter them.  This allows this routine to safely shuffle
            //    between several ongoing calculations.
            //
            //
            //    There are 28 compositions of 6 into three parts.  This routine will
            //    produce those compositions in the following order:
            //
            //     I         A
            //     -     ---------
            //     1     6   0   0
            //     2     5   1   0
            //     3     4   2   0
            //     4     3   3   0
            //     5     2   4   0
            //     6     1   5   0
            //     7     0   6   0
            //     8     5   0   1
            //     9     4   1   1
            //    10     3   2   1
            //    11     2   3   1
            //    12     1   4   1
            //    13     0   5   1
            //    14     4   0   2
            //    15     3   1   2
            //    16     2   2   2
            //    17     1   3   2
            //    18     0   4   2
            //    19     3   0   3
            //    20     2   1   3
            //    21     1   2   3
            //    22     0   3   3
            //    23     2   0   4
            //    24     1   1   4
            //    25     0   2   4
            //    26     1   0   5
            //    27     0   1   5
            //    28     0   0   6
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 July 2008
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Albert Nijenhuis, Herbert Wilf,
            //    Combinatorial Algorithms for Computers and Calculators,
            //    Second Edition,
            //    Academic Press, 1978,
            //    ISBN: 0-12-519260-6,
            //    LC: QA164.N54.
            //
            //  Parameters:
            //
            //    Input, int N, the integer whose compositions are desired.
            //
            //    Input, int K, the number of parts in the composition.
            //
            //    Input/output, int A[K], the parts of the composition.
            //
            //    Input/output, bool *MORE.
            //    Set MORE = FALSE on first call.  It will be reset to TRUE on return
            //    with a new composition.  Each new call returns another composition until
            //    MORE is set to FALSE when the last composition has been computed
            //    and returned.
            //
            //    Input/output, int *H, *T, two internal parameters needed for the
            //    computation.  The user should allocate space for these in the calling
            //    program, include them in the calling sequence, but never alter them!
            //
        {
            if (!(more))
            {
                t = n;
                h = 0;
                a[0] = n;
                for (int i = 1; i < k; i++)
                {
                    a[i] = 0;
                }
            }
            else
            {
                if (1 < t)
                {
                    h = 0;
                }

                h = h + 1;
                t = a[h - 1];
                a[h - 1] = 0;
                a[0] = t - 1;
                a[h] = a[h] + 1;
            }

            more = (a[k - 1] != n);
        }

        public static int compnz_enum ( int n, int k )

            /******************************************************************************/
            /*
              Purpose:
            
                COMPNZ_ENUM returns the number of nonzero compositions of the N into K parts.
            
              Discussion:
            
                A composition of the integer N into K nonzero parts is an ordered sequence
                of K positive integers which sum to N.  The compositions (1,2,1)
                and (1,1,2) are considered to be distinct.
            
                The 10 compositions of 6 into three nonzero parts are:
            
                  4 1 1,  3 2 1,  3 1 2,  2 3 1,  2 2 2,  2 1 3,
                  1 4 1,  1 3 2,  1 2 3,  1 1 4.
            
                The formula for the number of compositions of N into K nonzero
                parts is
            
                  Number = ( N - 1 )! / ( ( N - K )! * ( K - 1 )! )
            
                (Describe the composition using N-K '1's and K-1 dividing lines '|'.
                The number of distinct permutations of these symbols is the number
                of compositions into nonzero parts.  This is equal to the number of
                permutations of  N-1 things, with N-K identical of one kind
                and K-1 identical of another.)
            
                Thus, for the above example, we have:
            
                  Number = ( 6 - 1 )! / ( ( 6 - 3 )! * ( 3 - 1 )! ) = 10
            
              Licensing:
            
                This code is distributed under the GNU LGPL license.
            
              Modified:
            
                05 December 2013
            
              Author:
            
                John Burkardt
            
              Reference:
            
                Albert Nijenhuis, Herbert Wilf,
                Combinatorial Algorithms for Computers and Calculators,
                Second Edition,
                Academic Press, 1978,
                ISBN: 0-12-519260-6,
                LC: QA164.N54.
            
              Parameters:
            
                Input, int N, the integer whose compositions are desired.
            
                Input, int K, the number of parts in the composition.
            
                Output, int COMPNZ_ENUM, the number of compositions of N into
                K nonzero parts.
            */
        {
            int number;

            number = typeMethods.i4_choose ( n - 1, n - k );

            return number;
        }
        public static void compnz_next(ref CompNZData data, int n, int k, ref int[] a, ref bool more)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    COMPNZ_NEXT computes the compositions of integer N into K nonzero parts.
            //
            //  Discussion:
            //
            //    A composition of the integer N into K nonzero parts is an ordered sequence
            //    of K positive integers which sum to N.  The compositions (1,2,1)
            //    and (1,1,2) are considered to be distinct.
            //
            //    The routine computes one composition on each call until there are no more.
            //    For instance, one composition of 6 into 3 parts is 3+2+1, another would
            //    be 4+1+1 but 5+1+0 is not allowed since it includes a zero part.
            //
            //    On the first call to this routine, set MORE = FALSE.  The routine
            //    will compute the first element in the sequence of compositions, and
            //    return it, as well as setting MORE = TRUE.  If more compositions
            //    are desired, call again, and again.  Each time, the routine will
            //    return with a new composition.
            //
            //    However, when the LAST composition in the sequence is computed
            //    and returned, the routine will reset MORE to FALSE, signaling that
            //    the end of the sequence has been reached.
            //
            //  Example:
            //
            //    The 10 compositions of 6 into three nonzero parts are:
            //
            //      4 1 1,  3 2 1,  3 1 2,  2 3 1,  2 2 2,  2 1 3,
            //      1 4 1,  1 3 2,  1 2 3,  1 1 4.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    01 December 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Albert Nijenhuis, Herbert Wilf,
            //    Combinatorial Algorithms,
            //    Academic Press, 1978, second edition,
            //    ISBN 0-12-519260-6.
            //
            //  Parameters:
            //
            //    Input, int N, the integer whose compositions are desired.
            //
            //    Input, int K, the number of parts in the composition.
            //    K must be no greater than N.
            //
            //    Input/output, int A[K], the parts of the composition.
            //
            //    Input/output, bool *MORE.
            //    Set MORE = FALSE on first call.  It will be reset to TRUE on return
            //    with a new composition.  Each new call returns another composition until
            //    MORE is set to FALSE when the last composition has been computed
            //    and returned.
            //
        {
            int i;
            //
            //  We use the trick of computing ordinary compositions of (N-K)
            //  into K parts, and adding 1 to each part.
            //
            if (n < k)
            {
                more = false;
                for (i = 0; i < k; i++)
                {
                    a[i] = -1;
                }

                return;
            }

            //
            //  The first computation.
            //
            if (!(more))
            {
                data.t = n - k;
                data.h = 0;
                a[0] = n - k;
                for (i = 1; i < k; i++)
                {
                    a[i] = 0;
                }
            }
            else
            {
                for (i = 0; i < k; i++)
                {
                    a[i] = a[i] - 1;
                }

                if (1 < data.t)
                {
                    data.h = 0;
                }

                data.h = data.h + 1;
                data.t = a[data.h - 1];
                a[data.h - 1] = 0;
                a[0] = data.t - 1;
                a[data.h] = a[data.h] + 1;

            }

            more = (a[k - 1] != (n - k));

            for (i = 0; i < k; i++)
            {
                a[i] = a[i] + 1;
            }
        }

        public static void compnz_random(int n, int k, ref int seed, ref int[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    COMPNZ_RANDOM selects a random composition of the integer N into K nonzero parts.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    01 December 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Albert Nijenhuis, Herbert Wilf,
            //    Combinatorial Algorithms for Computers and Calculators,
            //    Second Edition,
            //    Academic Press, 1978,
            //    ISBN: 0-12-519260-6,
            //    LC: QA164.N54.
            //
            //  Parameters:
            //
            //    Input, int N, the integer to be decomposed.
            //
            //    Input, int K, the number of parts in the composition.
            //    K must be no greater than N.
            //
            //    Input/output, int &SEED, a seed for the random number generator.
            //
            //    Output, int A[K], the parts of the composition.
            //
        {
            int i;
            int l;
            int m;

            if (n < k)
            {
                for (i = 0; i < k; i++)
                {
                    a[i] = -1;
                }

                return;
            }

            if (1 < k)
            {
                Ksub.ksub_random(n - 1, k - 1, ref seed, ref a);
            }

            a[k - 1] = n;
            l = 0;

            for (i = 0; i < k; i++)
            {
                m = a[i];
                a[i] = a[i] - l - 1;
                l = m;
            }

            for (i = 0; i < k; i++)
            {
                a[i] = a[i] + 1;
            }

        }

        public static void compnz_to_ksub(int nc, int kc, int[] ac, ref int ns, ref int ks, ref int[] as_)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    COMPNZ_TO_KSUB converts a nonzero composition to a K-subset.
            //
            //  Discussion:
            //
            //    There is a bijection between K subsets and compositions.
            //
            //    Let AC be a composition of NC into KC parts.
            //
            //    Then let
            //      NS = NC - 1
            //      KS = KC - 1
            //    and define
            //      AS(I) = sum ( AC(1:I) ), for I = 1 : KS.
            //      
            //    Then AS is a KS subset of the integers 1 through NS.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    24 May 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NC, the composition sum.
            //
            //    Input, int KC, the number of parts of the composition.
            //
            //    Input, int AC[KC], the parts of the composition.
            //
            //    Output, int &NS, the size of the set.
            //
            //    Output, int &KS, the size of the subset.
            //
            //    Output, int AS[KS], the entries of the K-subset, 
            //    in increasing order.
            //
        {
            int i;

            ns = nc - 1;
            ks = kc - 1;
            as_[0] = ac[0];
            for (i = 1; i < kc - 1; i++)
            {
                as_[i] = as_[i - 1] + ac[i];
            }
        }

    }

    public class SubCompData
    {
        public bool more2;
        public int n2;
    }

    public static class SubComp
    {
        public static void subcomp_next(ref SubCompData data, int n, int k, ref int[] a, ref bool more, ref int h,
                ref int t)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SUBCOMP_NEXT computes the next subcomposition of N into K parts.
            //
            //  Discussion:
            //
            //    A composition of the integer N into K parts is an ordered sequence
            //    of K nonnegative integers which sum to a value of N.
            //
            //    A subcomposition of the integer N into K parts is a composition
            //    of M into K parts, where 0 <= M <= N.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 July 2008
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Albert Nijenhuis, Herbert Wilf.
            //    C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Albert Nijenhuis, Herbert Wilf,
            //    Combinatorial Algorithms for Computers and Calculators,
            //    Second Edition,
            //    Academic Press, 1978,
            //    ISBN: 0-12-519260-6,
            //    LC: QA164.N54.
            //
            //  Parameters:
            //
            //    Input, int N, the integer whose subcompositions are desired.
            //
            //    Input, int K, the number of parts in the subcomposition.
            //
            //    Input/output, int A[K], the parts of the subcomposition.
            //
            //    Input/output, bool *MORE, set by the user to start the computation,
            //    and by the routine to terminate it.
            //
            //    Input/output, int *H, *T, two internal parameters needed for the
            //    computation.  The user should allocate space for these in the calling
            //    program, include them in the calling sequence, but never alter them!
            //
        {
            int i;
            //
            //  The first computation.
            //
            if (!(more))
            {
                data.n2 = 0;

                for (i = 0; i < k; i++)
                {
                    a[i] = 0;
                }

                data.more2 = false;
                h = 0;
                t = 0;

                more = true;
            }
            //
            //  Do the next element at the current value of N.
            //
            else if (data.more2)
            {
                Comp.comp_next(data.n2, k, ref a, ref data.more2, ref h, ref t);
            }
            else
            {
                data.more2 = false;
                data.n2 = data.n2 + 1;

                Comp.comp_next(data.n2, k, ref a, ref data.more2, ref h, ref t);
            }

            //
            //  Termination occurs if MORE2 = FALSE and N2 = N.
            //
            if (!data.more2 && data.n2 == n)
            {
                more = false;
            }
        }
    }
}