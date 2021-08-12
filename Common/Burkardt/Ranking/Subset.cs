using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.RankingNS
{
    public static partial class Ranking
    {


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
            check = SubsetNS.Subset.subset_check(n, t);

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

            nsub = SubsetNS.Subset.subset_enum(n);

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

        public static int subset_gray_rank(int n, int[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SUBSET_GRAY_RANK ranks a subset of an N set, using the Gray code ordering.
            //
            //  Example:
            //
            //    N = 4
            //
            //       A       Rank
            //    -------   -----
            //
            //    0 0 0 0       1
            //    0 0 0 1       2
            //    0 0 1 1       3
            //    0 0 1 0       4
            //    0 1 1 0       5
            //    0 1 1 1       6
            //    0 1 0 1       7
            //    0 1 0 0       8
            //    1 1 0 0       9
            //    1 1 0 1      10
            //    1 1 1 1      11
            //    1 1 1 0      12
            //    1 0 1 0      13
            //    1 0 1 1      14
            //    1 0 0 1      15
            //    1 0 0 0      16
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    26 May 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of the total set from which
            //    subsets will be drawn.
            //
            //    Input, int A[N]; A(I) is 1 if element I is in the set,
            //    and 0 otherwise.
            //
            //    Output, int SUBSET_GRAY_RANK, the rank of the subset in the 
            //    Gray code ordering.
            //
        {
            int gray;
            int i;
            int rank;
            uint[] ua;
            uint ugray;

            ua = new uint[n];

            for (i = 0; i < n; i++)
            {
                ua[i] = (uint)a[i];
            }

            ugray = typeMethods.ubvec_to_ui4(n, ua);

            gray = (int)ugray;

            rank = gray_rank(gray);

            rank = rank + 1;

            return rank;
        }

        public static void subset_gray_unrank(int rank, int n, ref int[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SUBSET_GRAY_UNRANK produces a subset of an N set of the given Gray code rank.
            //
            //  Example:
            //
            //    N = 4
            //
            //     Rank     A    
            //    -----  -------
            //
            //        1  0 0 0 0
            //        2  0 0 0 1
            //        3  0 0 1 1
            //        4  0 0 1 0
            //        5  0 1 1 0
            //        6  0 1 1 1
            //        7  0 1 0 1
            //        8  0 1 0 0
            //        9  1 1 0 0
            //       10  1 1 0 1
            //       11  1 1 1 1
            //       12  1 1 1 0
            //       13  1 0 1 0
            //       14  1 0 1 1
            //       15  1 0 0 1
            //       16  1 0 0 0
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    28 June 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int RANK, the rank of the subset in the Gray code ordering.
            //
            //    Input, int N, the order of the total set from which
            //    subsets will be drawn.
            //
            //    Output, int A[N]; A(I) is 1 if element I is in the set,
            //    and 0 otherwise.
            //
        {
            int gray;
            int i;
            uint[] ua;
            uint ugray;

            gray = gray_unrank(rank - 1);

            ugray = (uint)gray;

            ua = new uint [n];

            typeMethods.ui4_to_ubvec(ugray, n, ref ua);

            for (i = 0; i < n; i++)
            {
                a[i] = (int)ua[i];
            }
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
            check = SubsetNS.Subset.subset_check(n, t);

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

            nsub = SubsetNS.Subset.subset_enum(n);

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

    }
}