using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.RankingNS
{
    public static partial class Ranking
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

        public static int subset_colex_rank(int n, int[] t)

            //****************************************************************************80
            // 
            //  Purpose:
            //
            //    SUBSET_COLEX_RANK computes the colexicographic rank of a subset.
            // 
            //  Licensing:
            // 
            //    This code is distributed under the GNU LGPL license.
            // 
            //  Modified:
            // 
            //    22 August 2011
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
            //    Input, int N, the number of items in the master set.
            //    N must be positive.
            // 
            //    Input, int T[N], the subset.  If T(I) = 0, item I is
            //    not in the subset; if T(I) = 1, item I is in the subset.
            // 
            //    Output, int SUBSET_COLEX_RANK, the rank of the subset.
            // 
        {
            bool check;
            int i;
            int rank;
            // 
            //  Check.
            // 
            check = subset_check(n, t);

            if (!check)
            {
                Console.WriteLine("");
                Console.WriteLine("SUBSET_COLEX_RANK - Fatal error!");
                Console.WriteLine("  The subset is illegal.");
                return(1);
            }

            rank = 0;

            for (i = 0; i < n; i++)
            {
                if (t[i] == 1)
                {
                    rank = rank + (int)Math.Pow(2, i);
                }
            }

            return rank;
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

        public static int[] subset_colex_unrank(int rank, int n)

            //****************************************************************************80
            // 
            //  Purpose:
            //
            //    SUBSET_COLEX_UNRANK computes the subset of given colexicographic rank.
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
            //    Input, int RANK, the rank of the subset.
            // 
            //    Input, int N, the number of items in the master set.
            //    N must be positive.
            // 
            //    Output, int SUBSET_COLEX_UNRANK[N], the subsetof the given rank.
            //    If T(I) = 0, item I is not in the subset; if T(I) = 1, item I is
            //    in the subset.
            // 
        {
            int i;
            int nsub;
            int rank_copy;
            int[] t;
            // 
            //  Check.
            // 
            if (n < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("SUBSET_COLEX_UNRANK - Fatal error!");
                Console.WriteLine("  Input N is illegal.");
                return null;
            }

            nsub = subset_enum(n);

            if (rank < 0 || nsub < rank)
            {
                Console.WriteLine("");
                Console.WriteLine("SUBSET_COLEX_UNRANK - Fatal error!");
                Console.WriteLine("  The input rank is illegal.");
                return null;
            }

            rank_copy = rank;
            t = new int[n];

            for (i = 0; i < n; i++)
            {
                if ((rank_copy % 2) == 1)
                {
                    t[i] = 1;
                }
                else
                {
                    t[i] = 0;
                }

                rank_copy = rank_copy / 2;
            }

            return t;
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

        public static int subset_lex_rank(int n, int[] t)

            //****************************************************************************80
            // 
            //  Purpose:
            //
            //    SUBSET_LEX_RANK computes the lexicographic rank of a subset.
            // 
            //  Licensing:
            // 
            //    This code is distributed under the GNU LGPL license.
            // 
            //  Modified:
            // 
            //    22 August 2011
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
            //    Input, int N, the number of items in the master set.
            //    N must be positive.
            // 
            //    Input, int T[N], the subset.  If T(I) = 0, item I is
            //    not in the subset; if T(I) = 1, item I is in the subset.
            // 
            //    Output, int SUBSET_LEX_RANK, the rank of the subset.
            // 
        {
            bool check;
            int i;
            int rank;
            // 
            //  Check.
            // 
            check = subset_check(n, t);

            if (!check)
            {
                Console.WriteLine("");
                Console.WriteLine("SUBSET_LEX_RANK - Fatal error!");
                Console.WriteLine("  The subset is illegal.");
                return (1);
            }

            rank = 0;

            for (i = 0; i < n; i++)
            {
                if (t[i] == 1)
                {
                    rank = rank + (int)Math.Pow(2, n - i - 1);
                }
            }

            return rank;
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

        public static int[] subset_lex_unrank(int rank, int n)

            //****************************************************************************80
            // 
            //  Purpose:
            //
            //    SUBSET_LEX_UNRANK computes the subset of given lexicographic rank.
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
            //    Input, int RANK, the rank of the subset.
            // 
            //    Input, int N, the number of items in the master set.
            //    N must be positive.
            // 
            //    Output, int SUBSET_LEX_UNRANK[N], the subset of the given rank.
            //    If T(I) = 0, item I is not in the subset; if T(I) = 1, item I is in
            //    the subset.
            // 
        {
            int i;
            int nsub;
            int rank_copy;
            int[] t;
            // 
            //  Check.
            // 
            if (n < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("SUBSET_LEX_UNRANK - Fatal error!");
                Console.WriteLine("  Input N is illegal.");
                return null;
            }

            nsub = subset_enum(n);

            if (rank < 0 || nsub < rank)
            {
                Console.WriteLine("");
                Console.WriteLine("SUBSET_LEX_UNRANK - Fatal error!");
                Console.WriteLine("  The input rank is illegal.");
                return null;
            }

            rank_copy = rank;
            t = new int[n];

            for (i = n - 1; 0 <= i; i--)
            {
                if ((rank_copy % 2) == 1)
                {
                    t[i] = 1;
                }
                else
                {
                    t[i] = 0;
                }

                rank_copy = rank_copy / 2;
            }

            return t;
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