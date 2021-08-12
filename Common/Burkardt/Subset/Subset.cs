using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.SubsetNS
{
    public static class Subset
    {
        public static bool subset_check(int n, int[] t)

            //****************************************************************************80
            // 
            //  Purpose:
            //
            //    SUBSET_CHECK checks a subset.
            // 
            //  Licensing:
            // 
            //    This code is distributed under the GNU LGPL license.
            // 
            //  Modified:
            // 
            //    28 July 2011
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
            //    Input, int N, the number of elements in the master set.
            //    N must be positive.
            // 
            //    Input, int T[N], the subset.  If T(I) = 0, item I is
            //    not in the subset; if T(I) = 1, item I is in the subset.
            // 
            //    Output, bool SUBSET_CHECK.
            //    TRUE, the data is legal.
            //    FALSE, the data is not legal.
        {
            bool check;
            int i;

            check = true;

            if (n < 1)
            {
                check = false;
                return check;
            }

            for (i = 0; i < n; i++)
            {
                if (t[i] != 0 && t[i] != 1)
                {
                    check = false;
                    return check;
                }
            }

            return check;
        }
        
                public static void subset_colex_successor(int n, ref int[] t, ref int rank )

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    SUBSET_COLEX_SUCCESSOR computes the subset colexicographic successor.
        // 
        //  Discussion:
        // 
        //    In the original code, there is a last element with no successor.
        // 
        //  Licensing:
        // 
        //    This code is distributed under the GNU LGPL license.
        // 
        //  Modified:
        // 
        //    26 July 2011
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
        //    Input, int N, the number of elements in the master set.
        //    N must be positive.
        // 
        //    Input/output, int T[N], describes a subset.  T(I) is 0 if
        //    the I-th element of the master set is not in the subset, and is
        //    1 if the I-th element is part of the subset.
        //    On input, T describes a subset.
        //    On output, T describes the next subset in the ordering.
        //    If the input T was the last in the ordering, then the output T
        //    will be the first.
        // 
        //    Input/output, int &RANK, the rank.
        //    If RANK = -1 on input, then the routine understands that this is
        //    the first call, and that the user wishes the routine to supply
        //    the first element in the ordering, which has RANK = 0.
        //    In general, the input value of RANK is increased by 1 for output,
        //    unless the very last element of the ordering was input, in which
        //    case the output value of RANK is 0.
        // 
        {
            bool check;
            int i;
            // 
            //  Return the first element.
            // 
            if (rank == -1)
            {
                for (i = 0; i < n; i++)
                {
                    t[i] = 0;
                }

                rank = 0;
                return;
            }

            // 
            //  Check.
            // 
            check = subset_check(n, t);

            if (!check)
            {
                Console.WriteLine("");
                Console.WriteLine("SUBSET_COLEX_SUCCESSOR - Fatal error!");
                Console.WriteLine("  The subset is illegal.");
                return;
            }

            for (i = 0; i < n; i++)
            {
                if (t[i] == 0)
                {
                    t[i] = 1;
                    rank = rank + 1;
                    return;
                }
                else
                {
                    t[i] = 0;
                }
            }

            rank = 0;
        }
                
                        public static int[] subset_complement(int n, int[] a)

            //****************************************************************************80
            // 
            //  Purpose:
            //
            //    SUBSET_COMPLEMENT computes the complement of a set.
            // 
            //  Licensing:
            // 
            //    This code is distributed under the GNU LGPL license.
            // 
            //  Modified:
            // 
            //    25 July 2011
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
            //    Input, int N, the order of the master set, of which A is
            //    a subset.  N must be positive.
            // 
            //    Input, int A[N], a subset of the master set.
            //    A(I) = 0 if the I-th element is in the subset A, and is
            //    1 otherwise.
            // 
            //    Output, int SUBSET_COMPLEMENT[N], the complement of A.
            // 
        {
            int[] b;
            bool check;
            int i;
            // 
            //  Check.
            // 
            check = subset_check(n, a);

            if (!check)
            {
                Console.WriteLine("");
                Console.WriteLine("SUBSET_COMPLEMENT - Fatal error!");
                Console.WriteLine("  The subset is illegal.");
                return null;
            }

            b = new int[n];

            for (i = 0; i < n; i++)
            {
                b[i] = 1 - a[i];
            }

            return b;
        }

        public static int subset_distance(int n, int[] t1, int[] t2 )

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    SUBSET_DISTANCE computes the Hamming distance between two sets.
        // 
        //  Discussion:
        // 
        //    The sets T1 and T2 are assumed to be subsets of a set of N elements.
        // 
        //  Licensing:
        // 
        //    This code is distributed under the GNU LGPL license.
        // 
        //  Modified:
        // 
        //    25 July 2011
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
        //    Input, int N, the order of the master set, of which T1 and
        //    T2 are subsets.  N must be positive.
        // 
        //    Input, int T1[N], T2[N], two subsets of the master set.
        //    T1(I) = 0 if the I-th element is in the subset T1, and is
        //    1 otherwise; T2 is defined similarly.
        // 
        //    Output, int SUBSET_DISTANCE, the Hamming distance between T1 and T2,
        //    defined as the number of elements of the master set which are
        //    in either T1 or T2 but not both.
        // 
        {
            bool check;
            int dist;
            int i;
            // 
            //  Check.
            // 
            check = subset_check(n, t1);

            if (!check)
            {
                Console.WriteLine("");
                Console.WriteLine("SUBSET_DISTANCE - Fatal error!");
                Console.WriteLine("  The subset is illegal.");
                return (1);
            }

            check = subset_check(n, t2);

            if (!check)
            {
                Console.WriteLine("");
                Console.WriteLine("SUBSET_DISTANCE - Fatal error!");
                Console.WriteLine("  The subset is illegal.");
                return (1);
            }

            dist = 0;

            for (i = 0; i < n; i++)
            {
                if ((t1[i] == 0 && t2[i] != 0) || (t1[i] != 0 && t2[i] == 0))
                {
                    dist = dist + 1;
                }
            }

            return dist;
        }

        public static int subset_enum(int n)

            //****************************************************************************80
            // 
            //  Purpose:
            //
            //    SUBSET_ENUM enumerates the subsets of a set with N elements.
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
            //  Parameters:
            // 
            //    Input, int N, the number of elements in the set.
            //    N must be at least 0.
            // 
            //    Output, int SUBSET_ENUM, the number of distinct elements.
            // 
        {
            int value = (int)Math.Pow(2, n);

            return value;
        }

        public static int[] subset_intersect(int n, int[] a, int[] b )

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    SUBSET_INTERSECT computes the intersection of two sets.
        // 
        //  Licensing:
        // 
        //    This code is distributed under the GNU LGPL license.
        // 
        //  Modified:
        // 
        //    25 July 2011
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
        //    Input, int N, the order of the master set, of which A and
        //    B are subsets.  N must be positive.
        // 
        //    Input, int A[N], B[N], two subsets of the master set.
        //    A(I) = 0 if the I-th element is in the subset A, and is
        //    1 otherwise; B is defined similarly.
        // 
        //    Output, int SUBSET_INTERSECT[N], the intersection of A and B.
        // 
        {
            int[] c;
            bool check;
            int i;
            // 
            //  Check.
            // 
            check = subset_check(n, a);

            if (!check)
            {
                Console.WriteLine("");
                Console.WriteLine("SUBSET_INTERSECTION - Fatal error!");
                Console.WriteLine("  The subset is illegal.");
                return null;
            }

            check = subset_check(n, b);

            if (!check)
            {
                Console.WriteLine("");
                Console.WriteLine("SUBSET_INTERSECTION - Fatal error!");
                Console.WriteLine("  The subset is illegal.");
                return null;
            }

            c = new int[n];

            for (i = 0; i < n; i++)
            {
                c[i] = Math.Min(a[i], b[i]);
            }

            return c;
        }

                public static void subset_lex_successor(int n, ref int[] t, ref int rank )

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    SUBSET_LEX_SUCCESSOR computes the subset lexicographic successor.
        // 
        //  Discussion:
        // 
        //    In the original code, there is a last element with no successor.
        // 
        //  Licensing:
        // 
        //    This code is distributed under the GNU LGPL license.
        // 
        //  Modified:
        // 
        //    27 July 2011
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
        //    Input, int N, the number of elements in the master set.
        //    N must be positive.
        // 
        //    Input/output, int T[N], describes a subset.  T(I) is 0 if
        //    the I-th element of the master set is not in the subset, and is
        //    1 if the I-th element is part of the subset.
        //    On input, T describes a subset.
        //    On output, T describes the next subset in the ordering.
        //    If the input T was the last in the ordering, then the output T
        //    will be the first.
        // 
        //    Input/output, int &RANK, the rank.
        //    If RANK = -1 on input, then the routine understands that this is
        //    the first call, and that the user wishes the routine to supply
        //    the first element in the ordering, which has RANK = 0.
        //    In general, the input value of RANK is increased by 1 for output,
        //    unless the very last element of the ordering was input, in which
        //    case the output value of RANK is 0.
        // 
        {
            bool check;
            int i;
            // 
            //  Return the first element.
            // 
            if (rank == -1)
            {
                for (i = 0; i < n; i++)
                {
                    t[i] = 0;
                }

                rank = 0;
                return;
            }

            // 
            //  Check.
            // 
            check = subset_check(n, t);

            if (!check)
            {
                Console.WriteLine("");
                Console.WriteLine("SUBSET_LEX_SUCCESSOR - Fatal error!");
                Console.WriteLine("  The subset is illegal.");
                return;
            }

            for (i = n - 1; 0 <= i; i--)
            {
                if (t[i] == 0)
                {
                    t[i] = 1;
                    rank = rank + 1;
                    return;
                }
                else
                {
                    t[i] = 0;
                }
            }

            rank = 0;
        }

        public static void subset_next ( int n, ref int[] t, ref int rank )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SUBSET_NEXT computes the subset lexicographic successor.
        //
        //  Discussion:
        //
        //    This is a lightly modified version of "subset_lex_successor()" from COMBO.
        //
        //  Example:
        //
        //    On initial call, N is 5 and the input value of RANK is -1.
        //    Then here are the successive outputs from the program:
        //
        //   Rank   T1   T2   T3   T4   T5
        //   ----   --   --   --   --   --
        //      0    0    0    0    0    0
        //      1    0    0    0    0    1
        //      2    0    0    0    1    0
        //      3    0    0    0    1    1
        //     ..   ..   ..   ..   ..   ..
        //     30    1    1    1    1    0
        //     31    1    1    1    1    1
        //     -1    0    0    0    0    0  <-- Reached end of cycle.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 May 2012
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
        //    Input, int N, the number of elements in the master set.
        //    N must be positive.
        //
        //    Input/output, int T[N], describes a subset.  T(I) is 0 if
        //    the I-th element of the master set is not in the subset, and is
        //    1 if the I-th element is part of the subset.
        //    On input, T describes a subset.
        //    On output, T describes the next subset in the ordering.
        //
        //    Input/output, int &RANK, the rank.
        //    If RANK = -1 on input, then the routine understands that this is
        //    the first call, and that the user wishes the routine to supply
        //    the first element in the ordering, which has RANK = 0.
        //    In general, the input value of RANK is increased by 1 for output,
        //    unless the very last element of the ordering was input, in which
        //    case the output value of RANK is -1.
        //
        {
            int i;
            //
            //  Return the first element.
            //
            if (rank == -1)
            {
                for (i = 0; i < n; i++)
                {
                    t[i] = 0;
                }

                rank = 0;
                return;
            }

            for (i = n - 1; 0 <= i; i--)
            {
                if (t[i] == 0)
                {
                    t[i] = 1;
                    rank = rank + 1;
                    return;
                }
                else
                {
                    t[i] = 0;
                }
            }

            rank = -1;
        }

        public static void subset_by_size_next(int n, ref int[] a, ref int subsize, ref bool more,
                ref bool more2, ref int m, ref int m2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SUBSET_BY_SIZE_NEXT returns all subsets of an N set, in order of size.
            //
            //  Example:
            //
            //    N = 4:
            //
            //    1 2 3 4
            //    1 2 3
            //    1 2 4
            //    1 3 4
            //    1 3
            //    1 4
            //    2 3
            //    1
            //    2
            //    3
            //    (the empty set)
            //
            //  Discussion:
            //
            //    The subsets are returned in decreasing order of size, with the
            //    empty set last.
            //
            //    For a given size K, the K subsets are returned in lexicographic order.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    09 June 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the size of the set.
            //
            //    Input/output, int A[N].  The entries A(1:SUBSIZE) contain
            //    the elements of the subset.  The elements are given in ascending
            //    order.
            //
            //    Output, int &SUBSIZE, the number of elements in the subset.
            //
            //    Input/output, bool &MORE.  Set MORE = FALSE before first call
            //    for a new sequence of subsets.  It then is set and remains
            //    TRUE as long as the subset computed on this call is not the
            //    final one.  When the final subset is computed, MORE is set to
            //    FALSE as a signal that the computation is done.
            //
            //    Input/output, bool &MORE2, a variable for bookkeeping.
            //    The user should declare this variable, but need not initialize it.
            //    The output value from one call must be the input value for the next.
            //
            //    Input/output, int &M, &M2, variables for bookkeeping.
            //    The user should declare this variable, but need not initialize it.
            //    The output value from one call must be the input value for the next.
            //
        {
            if (!more)
            {
                subsize = n;
                more = true;
                more2 = false;
                m = 0;
                m2 = 0;
            }
            else if (!more2)
            {
                subsize = subsize - 1;
            }

            //
            //  Compute the next subset of size SUBSIZE.
            //
            if (0 < subsize)
            {
                Ksub.ksub_next(n, subsize, ref a, ref more2, ref m, ref m2);
            }
            else if (subsize == 0)
            {
                more = false;
            }
        }

        public static void subset_lex_next(int n, bool jmp, int ndim, ref int k, ref int[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SUBSET_LEX_NEXT generates the subsets of a set of N elements, one at a time.
            //
            //  Discussion:
            //
            //    The subsets are generated in lexicographical order.  
            //
            //    The routine can also be forced to generate only those subsets whose 
            //    size is no greater than some user-specified maximum.
            //
            //  Example:
            //
            //    N = 5, JMP = ( K == 3 )
            //
            //    1
            //    1 2
            //    1 2 3
            //    1 2 4
            //    1 2 5
            //    1 3
            //    1 3 4
            //    1 3 5
            //    1 4
            //    1 4 5
            //    1 5
            //    2
            //    2 3
            //    2 3 4
            //    2 3 5
            //    2 4
            //    2 4 5
            //    2 5
            //    3
            //    3 4
            //    3 4 5
            //    3 5
            //    4
            //    4 5
            //    5
            //    empty set.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 November 2004
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
            //    Input, int N, the order of the main set from which subsets
            //    are chosen.
            //
            //    Input, bool JMP.  In the simplest case, set JMP = .FALSE. for
            //    a normal computation.  But to jump over supersets of the input set,
            //    set JMP = TRUE.  Setting JMP = ( K == 3 ) before every new call
            //    will, for example, force all the subsets returned
            //    to have cardinality 3 or less.
            //
            //    Input, int NDIM, the allowed storage for A.  If NDIM < N,
            //    JMP must be used to avoid creation of a subset too large to store in A.
            //
            //    Input/output, int &K.  On first call, the user must set K = 0 as
            //    a startup signal to the program.  Thereafter, the routine returns
            //    the size of the computed subset in K.  On the last return,
            //    the empty set is returned and K is 0, which is a signal to
            //    the user that the computation is complete.
            //
            //    Input/output, int A[NDIM].  A(I) is the I-th element of the
            //    subset, listed in increasing order, with 0's in entries
            //    beyond entry K.
            //
        {
            int is_;

            if (k <= 0)
            {
                if (jmp)
                {
                    return;
                }
                is_ = 0;
                k = 1;
                a[0] = 1;
            }
            else if (a[k - 1] != n)
            {
                is_ = a[k - 1];

                if (!jmp)
                {
                    k = k + 1;
                }

                a[k - 1] =  is_ +1;
            }
            else
            {
                k = k - 1;

                if (k != 0)
                {
                    a[k - 1] = a[k - 1] + 1;
                }
            }
        }

        public static void subset_gray_next(int n, ref int[] a, ref bool more, ref int ncard, ref int iadd)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SUBSET_GRAY_NEXT generates all subsets of a set of order N, one at a time.
            //
            //  Discussion:
            //
            //    It generates the subsets one at a time, by adding or subtracting
            //    exactly one element on each step.
            //
            //    The user should set MORE = .FALSE. and the value of N before
            //    the first call.  On return, the user may examine A which contains
            //    the definition of the new subset, and must check .MORE., because
            //    as soon as it is .FALSE. on return, all the subsets have been
            //    generated and the user probably should cease calling.
            //
            //    The first set returned is the empty set.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 May 2003
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
            //    Input, int N, the order of the total set from which
            //    subsets will be drawn.
            //
            //    Input/output, int A[N].  On each return, the Gray code for the newly
            //    generated subset.  A[I] = 0 if element I is in the subset, 1 otherwise.
            //
            //    Input/output, bool &MORE.  Set this variable FALSE before
            //    the first call.  Normally, MORE will be returned TRUE but once
            //    all the subsets have been generated, MORE will be
            //    reset FALSE on return and you should stop calling the program.
            //
            //    Input/output, int &NCARD, the cardinality of the set returned,
            //    which may be any value between 0 (the empty set) and N (the
            //    whole set).
            //
            //    Output, int &IADD, the element which was added or removed to the
            //    previous subset to generate the current one.  Exception:
            //    the empty set is returned on the first call, and IADD is set to -1.
        {
            int i;
            //
            //  First set returned is the empty set.
            //
            if (!more)
            {
                for (i = 0; i < n; i++)
                {
                    a[i] = 0;
                }

                iadd = 0;
                ncard = 0;
                more = true;
            }
            else
            {
                iadd = 1;

                if ((ncard % 2) != 0)
                {
                    for (;;)
                    {
                        iadd = iadd + 1;
                        if (a[iadd - 2] != 0)
                        {
                            break;
                        }
                    }
                }

                a[iadd - 1] = 1 - a[iadd - 1];
                ncard = ncard + 2 * a[iadd - 1] - 1;
                //
                //  Last set returned is the singleton A(N).
                //
                if (ncard == a[n - 1])
                {
                    more = false;
                }
            }
        }
        
                public static int[] subset_random(int n, ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SUBSET_RANDOM returns a random subset.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    24 December 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the size of the set.
            //
            //    Input/output, int &SEED, a seed for the random number
            //    generator.
            //
            //    Output, int SUBSET_RANDOM[N], defines the subset using 0 and 1 values.
            //
        {
            int[] s = UniformRNG.i4vec_uniform_ab_new(n, 0, 1, ref seed);

            return s;
        }

        public static int[] subset_union(int n, int[] a, int[] b )

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    SUBSET_UNION computes the union of two sets.
        // 
        //  Licensing:
        // 
        //    This code is distributed under the GNU LGPL license.
        // 
        //  Modified:
        // 
        //    25 July 2011
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
        //    Input, int N, the order of the master set, of which A and
        //    B are subsets.  N must be positive.
        // 
        //    Input, int A[N], B[N], two subsets of the master set.
        //    A(I) = 0 if the I-th element is in the subset A, and is
        //    1 otherwise; B is defined similarly.
        // 
        //    Output, int SUBSET_UNION[N], the union of A and B.
        // 
        {
            int[] c;
            bool check;
            int i;
            // 
            //  Check.
            // 
            check = subset_check(n, a);

            if (!check)
            {
                Console.WriteLine("");
                Console.WriteLine("SUBSET_UNION - Fatal error!");
                Console.WriteLine("  The subset is illegal.");
                return null;
            }

            check = subset_check(n, b);

            if (!check)
            {
                Console.WriteLine("");
                Console.WriteLine("SUBSET_UNION - Fatal error!");
                Console.WriteLine("  The subset is illegal.");
                return null;
            }

            c = new int[n];

            for (i = 0; i < n; i++)
            {
                c[i] = Math.Max(a[i], b[i]);
            }

            return c;
        }

        public static int subset_weight(int n, int[] t)

            //****************************************************************************80
            // 
            //  Purpose:
            //
            //    SUBSET_WEIGHT computes the Hamming weight of a set.
            // 
            //  Discussion:
            // 
            //    The Hamming weight is simply the number of elements in the set.
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
            //    Input, int N, the order of the master set, of which T
            //    is a subset.  N must be positive.
            // 
            //    Input, int T[N], defines the subset T.
            //    T(I) is 1 if I is an element of T, and 0 otherwise.
            // 
            //    Output, int SUBSET_WEIGHT, the Hamming weight of the subset T.
            // 
        {
            bool check;
            int weight;
            // 
            //  Check.
            // 
            check = subset_check(n, t);

            if (!check)
            {
                Console.WriteLine("");
                Console.WriteLine("SUBSET_WEIGHT - Fatal error!");
                Console.WriteLine("  The subset is illegal.");
                return (1);
            }

            weight = typeMethods.i4vec_sum(n, t);

            return weight;
        }

        public static int[] subset_xor(int n, int[] a, int[] b )

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    SUBSET_XOR computes the symmetric difference of two sets.
        // 
        //  Licensing:
        // 
        //    This code is distributed under the GNU LGPL license.
        // 
        //  Modified:
        // 
        //    25 July 2011
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
        //    Input, int N, the order of the master set, of which A and
        //    B are subsets.  N must be positive.
        // 
        //    Input, int A[N], B[N], two subsets of the master set.
        //    A(I) = 0 if the I-th element is in the subset A, and is
        //    1 otherwise; B is defined similarly.
        // 
        //    Output, int SUBSET_XOR[N], the symmetric difference of A and B.
        // 
        {
            int[] c;
            bool check;
            int i;
            // 
            //  Check.
            // 
            check = subset_check(n, a);

            if (!check)
            {
                Console.WriteLine("");
                Console.WriteLine("SUBSET_XOR - Fatal error!");
                Console.WriteLine("  The subset is illegal.");
                return null;
            }

            check = subset_check(n, b);

            if (!check)
            {
                Console.WriteLine("");
                Console.WriteLine("SUBSET_XOR - Fatal error!");
                Console.WriteLine("  The subset is illegal.");
                return null;
            }

            c = new int[n];

            for (i = 0; i < n; i++)
            {
                c[i] = Math.Max(a[i], b[i]) - Math.Min(a[i], b[i]);
            }

            return c;
        }

        public static int subsetsum_swap(int n, ref int[] a, int sum_desired, ref int[] index )

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    SUBSETSUM_SWAP seeks a solution of the subset sum problem by swapping.
        // 
        //  Discussion:
        // 
        //    Given a collection of N not necessarily distinct positive integers A(I),
        //    and a positive integer SUM_DESIRED, select a subset of the values so that
        //    their sum is as close as possible to SUM_DESIRED without exceeding it.
        // 
        //  Algorithm:
        // 
        //    Start with no values selected, and SUM_ACHIEVED = 0.
        // 
        //    Consider each element A(I):
        // 
        //      If A(I) is not selected and SUM_ACHIEVED + A(I) <= SUM_DESIRED,
        //        select A(I).
        // 
        //      If A(I) is still not selected, and there is a selected A(J)
        //      such that SUM_GOT < SUM_ACHIEVED + A(I) - A(J),
        //        select A(I) and deselect A(J).
        // 
        //      If no items were selected on this sweep,
        //        exit.
        //      Otherwise,
        //        repeat the search.
        // 
        //  Licensing:
        // 
        //    This code is distributed under the GNU LGPL license.
        // 
        //  Modified:
        // 
        //    25 July 2011
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
        //    Input, int N, the number of values.  N must be positive.
        // 
        //    Input/output, int A[N], a collection of positive values.
        //    On output, A has been sorted into descending order.
        // 
        //    Input, int SUM_DESIRED, the desired sum.
        // 
        //    Output, int INDEX[N]; INDEX(I) is 1 if A(I) is part of the
        //    sum, and 0 otherwise.
        // 
        //    Output, int SUBSETSUM_SWAP, the sum of the selected
        //    elements.
        // 
        {
            int i;
            int j;
            int nmove;
            int sum_achieved;
            // 
            //  Initialize.
            // 
            sum_achieved = 0;

            for (i = 0; i < n; i++)
            {
                index[i] = 0;
            }

            // 
            //  Sort into descending order.
            // 
            typeMethods.i4vec_sort_insert_d(n, ref a);

            for (;;)
            {
                nmove = 0;

                for (i = 0; i < n; i++)
                {
                    if (index[i] == 0)
                    {
                        if (sum_achieved + a[i] <= sum_desired)
                        {
                            index[i] = 1;
                            sum_achieved = sum_achieved + a[i];
                            nmove = nmove + 1;
                            continue;
                        }
                    }

                    if (index[i] == 0)
                    {
                        for (j = 0; j < n; j++)
                        {
                            if (index[j] == 1)
                            {
                                if (sum_achieved < sum_achieved + a[i] - a[j] &&
                                    sum_achieved + a[i] - a[j] <= sum_desired)
                                {
                                    index[j] = 0;
                                    index[i] = 1;
                                    nmove = nmove + 2;
                                    sum_achieved = sum_achieved + a[i] - a[j];
                                    break;
                                }
                            }
                        }
                    }
                }

                if (nmove <= 0)
                {
                    break;
                }
            }

            return sum_achieved;
        }

    }
}