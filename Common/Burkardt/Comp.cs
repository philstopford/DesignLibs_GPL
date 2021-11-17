﻿using System;
using Burkardt.SubsetNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.Composition;

public class CompNZData
{
    public int h;
    public int t;
}

public static class Comp
{
    public static int[] get_seq(int d, int norm, int seq_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GET_SEQ generates all positive integer D-vectors that sum to NORM.
        //
        //  Discussion:
        //
        //    This function computes a list, in reverse dictionary order, of
        //    all D-vectors of positive values that sum to NORM.
        //
        //    For example, call get_seq ( 3, 5, 6, fs ) returns
        //
        //      3  1  1
        //      2  2  1
        //      2  1  2
        //      1  3  1
        //      1  2  2
        //      1  1  3
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 December 2012
        //
        //  Author:
        //
        //    Original MATLAB version by Florian Heiss, Viktor Winschel.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Florian Heiss, Viktor Winschel,
        //    Likelihood approximation by numerical integration on sparse grids,
        //    Journal of Econometrics,
        //    Volume 144, 2008, pages 62-80.
        //
        //  Parameters:
        //
        //    Input, int D, the dimension.
        //    1 <= D.
        //
        //    Input, int NORM, the value that each row must sum to.
        //    D <= NORM.
        //
        //    Input, int SEQ_NUM, the number of rows of FS.
        //
        //    Output, int GET_SEQ[SEQ_NUM*D].  Each row of FS represents 
        //    one vector with all elements positive and summing to NORM.
        //
    {
        int i;

        int[] seq = new int[d];
        int[] fs = new int[seq_num * d];

        if (norm < d)
        {
            Console.WriteLine("");
            Console.WriteLine("GET_SEQ - Fatal error!");
            Console.WriteLine("  NORM = " + norm + " < D = " + "");
            return null;
        }

        for (i = 0; i < d; i++)
        {
            seq[i] = 0;
        }

        //
        //  The algorithm is written to work with vectors whose minimum value is
        //  allowed to be zero.  So we subtract D from NORM at the beginning and
        //  then increment the result vectors by 1 at the end!
        //
        int a = norm - d;
        seq[0] = a;

        int row = 0;
        for (i = 0; i < d; i++)
        {
            fs[row + i * seq_num] = seq[i] + 1;
        }

        int c = 0;

        while (seq[d - 1] < a)
        {
            if (c == d - 1)
            {
                for (i = c - 1; 0 <= i; i--)
                {
                    c = i;
                    if (seq[i] != 0)
                    {
                        break;
                    }
                }
            }

            seq[c] -= 1;
            c += 1;
            seq[c] = a;
            for (i = 0; i < c; i++)
            {
                seq[c] -= seq[i];
            }

            if (c < d - 1)
            {
                for (i = c + 1; i < d; i++)
                {
                    seq[i] = 0;
                }
            }

            row += 1;
            for (i = 0; i < d; i++)
            {
                fs[row + i * seq_num] = seq[i] + 1;
            }
        }

        return fs;
    }

    public static int num_seq ( int n, int k )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NUM_SEQ returns the number of compositions of the integer N into K parts.
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
        //    08 December 2012
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
        //    Output, int NUM_SEQ, the number of compositions of N
        //    into K parts.
        //
    {
        int value = typeMethods.i4_choose ( n + k - 1, n );

        return value;
    }
        
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
        int number = typeMethods.i4_choose(n + k - 1, n);

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
        switch (kc)
        {
            //
            //  Ensure that 1 <= KC.
            //
            case < 1:
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
            switch (xc[i])
            {
                case < 0:
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
            if (0 >= xc[j - 1])
            {
                continue;
            }

            i = j;
            break;
        }

        switch (i)
        {
            //
            //  set T = X(I)
            //  set XC(I) to zero,
            //  increase XC(I-1) by 1,
            //  increment XC(KC) by T-1.
            //
            case 0:
                xc[kc - 1] = 1;
                return;
            case 1:
                t = xc[0] + 1;
                im1 = kc;
                break;
            case > 1:
                t = xc[i - 1];
                im1 = i - 1;
                break;
        }

        xc[i - 1] = 0;
        xc[im1 - 1] += 1;
        xc[kc - 1] = xc[kc - 1] + t - 1;
    }

    public static void comp_random(int n, int k, ref int seed, ref int[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COMP_RANDOM selects a random composition of the integer N into K parts.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 May 2015
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
        //    Input, int N, the integer to be decomposed.
        //
        //    Input, int K, the number of parts in the composition.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, int A[K], the parts of the composition.
        //
    {
        int i;

        int[] b = new int[k - 1];

        Ksub.ksub_random2(n + k - 1, k - 1, ref seed, ref b);

        for (i = 0; i < k - 1; i++)
        {
            a[i] = b[i];
        }

        a[k - 1] = n + k;

        int l = 0;

        for (i = 0; i < k; i++)
        {
            int m = a[i];
            a[i] = a[i] - l - 1;
            l = m;
        }
    }

    public static int[] comp_random_grlex(int kc, int rank1, int rank2, ref int seed, ref int rank)

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
        switch (kc)
        {
            //
            //  Ensure that 1 <= KC.
            //
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("COMP_RANDOM_GRLEX - Fatal error!");
                Console.WriteLine("  KC < 1");
                return null;
        }

        switch (rank1)
        {
            //
            //  Ensure that 1 <= RANK1.
            //
            case < 1:
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
        int[] xc = comp_unrank_grlex(kc, rank);

        return xc;
    }
        
    public static int[] comp_random_new ( int n, int k, ref int seed )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    COMP_RANDOM_NEW selects a random composition of the integer N into K parts.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 April 2003
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
        //    Input, int N, the integer to be decomposed.
        //
        //    Input, int K, the number of parts in the composition.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, int COMP_RANDOM_NEW[K], the parts of the composition.
        //
    {
        int i;

        int[] a = new int[k];

        Ksub.ksub_random2 ( n+k-1, k-1, ref seed, ref a );

        a[k-1] = n + k;
        int l = 0;

        for ( i = 0; i < k; i++ )
        {
            int m = a[i];
            a[i] = a[i] - l - 1;
            l = m;
        }

        return a;
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
        int n;
        switch (kc)
        {
            //
            //  Ensure that 1 <= KC.
            //
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("COMP_RANK_GRLEX - Fatal error!");
                Console.WriteLine("  KC < 1");
                return 1;
        }

        //
        //  Ensure that 0 <= XC(I).
        //
        for (i = 0; i < kc; i++)
        {
            switch (xc[i])
            {
                case < 0:
                    Console.WriteLine("");
                    Console.WriteLine("COMP_RANK_GRLEX - Fatal error!");
                    Console.WriteLine("  XC[I] < 0");
                    return 1;
            }
        }

        //
        //  NC = sum ( XC )
        //
        int nc = typeMethods.i4vec_sum(kc, xc);
        //
        //  Convert to KSUBSET format.
        //
        int ns = nc + kc - 1;
        int ks = kc - 1;
        int[] xs = new int[ks];
        xs[0] = xc[0] + 1;
        for (i = 2; i < kc; i++)
        {
            xs[i - 1] = xs[i - 2] + xc[i - 1] + 1;
        }

        //
        //  Compute the rank.
        //
        int rank = 1;

        for (i = 1; i <= ks; i++)
        {
            int tim1 = i switch
            {
                1 => 0,
                _ => xs[i - 2]
            };

            if (tim1 + 1 > xs[i - 1] - 1)
            {
                continue;
            }

            int j;
            for (j = tim1 + 1; j <= xs[i - 1] - 1; j++)
            {
                rank += typeMethods.i4_choose(ns - j, ks - i);
            }
        }

        for (n = 0; n < nc; n++)
        {
            rank += typeMethods.i4_choose(n + kc - 1, n);
        }

        return rank;
    }

    public static void comp_to_ksub(int nc, int kc, int[] ac, ref int ns, ref int ks, ref int[] as_)

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
        for (i = 1; i < kc - 1; i++)
        {
            as_[i] = as_[i - 1] + ac[i] + 1;
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
        int r;
        switch (kc)
        {
            //
            //  Ensure that 1 <= KC.
            //
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("COMP_UNRANK_GRLEX - Fatal error!");
                Console.WriteLine("  KC < 1");
                return null;
        }

        switch (rank)
        {
            //
            //  Ensure that 1 <= RANK.
            //
            case < 1:
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
        int rank1 = 1;
        int nc = -1;
        for (;;)
        {
            nc += 1;
            r = typeMethods.i4_choose(nc + kc - 1, nc);
            if (rank < rank1 + r)
            {
                break;
            }

            rank1 += r;
        }

        int rank2 = rank - rank1;
        //
        //  Convert to KSUBSET format.
        //  Apology: an unranking algorithm was available for KSUBSETS,
        //  but not immediately for compositions.  One day we will come back
        //  and simplify all this.
        //
        int ks = kc - 1;
        int ns = nc + kc - 1;
        int[] xs = new int[ks];

        int j = 1;

        for (i = 1; i <= ks; i++)
        {
            r = typeMethods.i4_choose(ns - j, ks - i);

            while (r <= rank2 && 0 < r)
            {
                rank2 -= r;
                j += 1;
                r = typeMethods.i4_choose(ns - j, ks - i);
            }

            xs[i - 1] = j;
            j += 1;
        }

        //
        //  Convert from KSUBSET format to COMP format.
        //
        int[] xc = new int[kc];
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
        switch (more)
        {
            case false:
            {
                t = n;
                h = 0;
                a[0] = n;
                for (int i = 1; i < k; i++)
                {
                    a[i] = 0;
                }

                break;
            }
            default:
            {
                h = t switch
                {
                    > 1 => 0,
                    _ => h
                };

                h += 1;
                t = a[h - 1];
                a[h - 1] = 0;
                a[0] = t - 1;
                a[h] += 1;
                break;
            }
        }

        more = a[k - 1] != n;
    }

    public static int compnz_enum(int n, int k)

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
        int number = typeMethods.i4_choose(n - 1, n - k);

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

        switch (more)
        {
            //
            //  The first computation.
            //
            case false:
            {
                data.t = n - k;
                data.h = 0;
                a[0] = n - k;
                for (i = 1; i < k; i++)
                {
                    a[i] = 0;
                }

                break;
            }
            default:
            {
                for (i = 0; i < k; i++)
                {
                    a[i] -= 1;
                }

                data.h = data.t switch
                {
                    > 1 => 0,
                    _ => data.h
                };

                data.h += 1;
                data.t = a[data.h - 1];
                a[data.h - 1] = 0;
                a[0] = data.t - 1;
                a[data.h] += 1;
                break;
            }
        }

        more = a[k - 1] != n - k;

        for (i = 0; i < k; i++)
        {
            a[i] += 1;
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

        if (n < k)
        {
            for (i = 0; i < k; i++)
            {
                a[i] = -1;
            }

            return;
        }

        switch (k)
        {
            case > 1:
                Ksub.ksub_random(n - 1, k - 1, ref seed, ref a);
                break;
        }

        a[k - 1] = n;
        int l = 0;

        for (i = 0; i < k; i++)
        {
            int m = a[i];
            a[i] = a[i] - l - 1;
            l = m;
        }

        for (i = 0; i < k; i++)
        {
            a[i] += 1;
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
        switch (more)
        {
            //
            //  The first computation.
            //
            case false:
            {
                data.n2 = 0;

                int i;
                for (i = 0; i < k; i++)
                {
                    a[i] = 0;
                }

                data.more2 = false;
                h = 0;
                t = 0;

                more = true;
                break;
            }
            //
            default:
            {
                switch (data.more2)
                {
                    case true:
                        Comp.comp_next(data.n2, k, ref a, ref data.more2, ref h, ref t);
                        break;
                    default:
                        data.more2 = false;
                        data.n2 += 1;

                        Comp.comp_next(data.n2, k, ref a, ref data.more2, ref h, ref t);
                        break;
                }

                break;
            }
        }

        more = data.more2 switch
        {
            //
            //  Termination occurs if MORE2 = FALSE and N2 = N.
            //
            false when data.n2 == n => false,
            _ => more
        };
    }

    public static void subcompnz_next(ref CompNZData data, int n, int k, ref int[] a, ref bool more,
            ref int h, ref int t, ref int n2, ref bool more2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SUBCOMPNZ_NEXT computes the next subcomposition of N into K nonzero parts.
        //
        //  Discussion:
        //
        //    A composition of the integer N into K nonzero parts is an ordered sequence
        //    of K positive integers which sum to a value of N.
        //
        //    A subcomposition of the integer N into K nonzero parts is a composition
        //    of M into K nonzero parts, where 0 < M <= N.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 June 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the integer whose subcompositions are desired.
        //
        //    Input, int K, the number of parts in the subcomposition.
        //    K must be no greater than N.
        //
        //    Input/output, int A[K], the parts of the subcomposition.
        //
        //    Input/output, bool &MORE, set by the user to start the computation,
        //    and by the routine to terminate it.
        //
        //    Input/output, int &H, &T, &N2, internal parameters needed for the
        //    computation.  The user should allocate space for these in the calling
        //    program, include them in the calling sequence, but never alter them!
        //
        //    Input/output, bool &MORE2, an internal parameter needed for the
        //    computation.  The user should allocate space for it in the calling
        //    program, include it in the calling sequence, but never alter it!
        //
    {
        int i;

        if (n < k)
        {
            for (i = 0; i < k; i++)
            {
                a[i] = -1;
            }

            return;
        }

        switch (more)
        {
            //
            //  The first computation.
            //
            case false:
            {
                for (i = 0; i < k; i++)
                {
                    a[i] = 1;
                }

                more = true;
                h = 0;
                t = 0;
                n2 = k;
                more2 = false;
                break;
            }
            //
            default:
            {
                switch (more2)
                {
                    case true:
                        Comp.compnz_next(ref data, n2, k, ref a, ref more2);
                        break;
                    default:
                        more2 = false;
                        n2 += 1;

                        Comp.compnz_next(ref data, n2, k, ref a, ref more2);
                        break;
                }

                break;
            }
        }

        more = more2 switch
        {
            //
            //  Termination occurs if MORE2 = FALSE and N2 = N.
            //
            false when n2 == n => false,
            _ => more
        };
    }

    public static void subcompnz2_next(ref CompNZData data, int n_lo, int n_hi, int k, ref int[] a, ref bool more,
            ref int h, ref int t, ref int n2, ref bool more2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SUBCOMPNZ2_NEXT computes the next subcomposition of N into K nonzero parts.
        //
        //  Discussion:
        //
        //    A composition of the integer N into K nonzero parts is an ordered sequence
        //    of K positive integers which sum to a value of N.
        //
        //    A subcomposition of the integer N into K nonzero parts is a composition
        //    of M into K nonzero parts, where 0 < M <= N.
        //
        //    This routine computes all compositions of K into nonzero parts which sum
        //    to values between N_LO and N_HI.
        //
        //    The routine SUBCOMPNZ_NEXT can be regarded as a special case 
        //    where N_LO = K.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 June 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N_LO, N_HI, the range of values of N for which compositions
        //    are desired.
        //    N_LO must be no greater than N_HI.
        //
        //    Input, int K, the number of parts in the subcomposition.
        //    K must be no greater than N_HI.
        //
        //    Input/output, int A[K], the parts of the subcomposition.
        //
        //    Input/output, bool &MORE, set by the user to start the computation,
        //    and by the routine to terminate it.
        //
        //    Input/output, int &H, &T, &N2, internal parameters needed for the
        //    computation.  The user should allocate space for these in the calling
        //    program, include them in the calling sequence, but never alter them!
        //
        //    Input/output, bool &MORE2, an internal parameter needed for the
        //    computation.  The user should allocate space for it in the calling
        //    program, include it in the calling sequence, but never alter it!
        //
    {
        int i;

        if (n_hi < k)
        {
            for (i = 0; i < k; i++)
            {
                a[i] = -1;
            }

            return;
        }

        if (n_hi < n_lo)
        {
            for (i = 0; i < k; i++)
            {
                a[i] = -1;
            }

            return;
        }

        switch (more)
        {
            //
            //  The first computation.
            //
            case false:
                more = true;
                h = 0;
                t = 0;
                n2 = Math.Max(k, n_lo);
                more2 = false;

                Comp.compnz_next(ref data, n2, k, ref a, ref more2);
                break;
            //
            default:
            {
                switch (more2)
                {
                    case true:
                        Comp.compnz_next(ref data, n2, k, ref a, ref more2);
                        break;
                    default:
                        n2 += 1;

                        Comp.compnz_next(ref data, n2, k, ref a, ref more2);
                        break;
                }

                break;
            }
        }

        more = more2 switch
        {
            //
            //  Termination occurs if MORE2 = FALSE and N2 = N_HI.
            //
            false when n2 == n_hi => false,
            _ => more
        };
    }
}