using System;
using Burkardt.Types;

namespace Burkardt
{
    public class CompNZData
    {
        public int h = 0;
        public int t = 0;
    }

    public static class Comp
    {
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


    }

    public class SubCompData
    {
        public bool more2;
        public int n2;
    }

    public static class SubComp
    {
        public static void subcomp_next(ref SubCompData data, int n, int k, ref int[] a, ref bool more, ref int h, ref int t )

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