using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.RankingNS
{
    public static partial class Ranking
    {
        public static bool perm_check(int n, int[] p)

            //****************************************************************************80
            // 
            //  Purpose:
            //
            //    PERM_CHECK checks a representation of a permutation.
            // 
            //  Discussion:
            // 
            //    The routine is given N and P, a vector of length N.
            //    P is a legal represention of a permutation of the integers from
            //    1 to N if and only if every integer from 1 to N occurs
            //    as a value of P(I) for some I between 1 and N.
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
            //  Parameters:
            // 
            //    Input, int N, the number of values being permuted.
            //    N must be positive.
            // 
            //    Input, int P[N], the array to check.
            //
            //    Output, bool PERM_CHECK.
            //    TRUE, the data is legal.
            //    FALSE, the data is not legal.
        {
            bool check;
            int i;
            int ifind;
            int iseek;

            check = true;

            if (n < 1)
            {
                check = false;
                return check;
            }

            for (i = 0; i < n; i++)
            {
                if (p[i] < 1 || n < p[i])
                {
                    check = false;
                    return check;
                }
            }

            for (iseek = 1; iseek <= n; iseek++)
            {
                ifind = -1;
                for (i = 0; i < n; i++)
                {
                    if (p[i] == iseek)
                    {
                        ifind = i;
                        break;
                    }
                }

                if (ifind == -1)
                {
                    check = false;
                    return check;
                }
            }

            return check;
        }

        public static int perm_enum(int n)

            //****************************************************************************80
            // 
            //  Purpose:
            //
            //    PERM_ENUM enumerates the permutations on N digits.
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
            //    Input, int N, the number of values being permuted.
            //    N must be nonnegative.
            // 
            //    Output, int PERM_ENUM, the number of distinct elements.
            // 
        {
            int value = typeMethods.i4_factorial(n);

            return value;
        }

        public static int[] perm_inv(int n, int[] p)

            //****************************************************************************80
            // 
            //  Purpose:
            //
            //    PERM_INV computes the inverse of a permutation.
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
            //    Input, int N, the number of values being permuted.
            //    N must be positive.
            // 
            //    Input, int P[N], describes the permutation.
            //    P(I) is the item which is permuted into the I-th place
            //    by the permutation.
            // 
            //    Output, int PERM_INV[N], the inverse permutation.
            // 
        {
            bool check;
            int i;
            int[] pinv;
            // 
            //  Check.
            // 
            check = perm_check(n, p);

            if (!check)
            {
                Console.WriteLine("");
                Console.WriteLine("PERM_INV - Fatal error!");
                Console.WriteLine("  Permutation is illegal.");
                return null;
            }

            pinv = new int[n];

            for (i = 0; i < n; i++)
            {
                pinv[p[i] - 1] = i + 1;
            }

            return pinv;
        }

        public static int perm_lex_rank(int n, int[] p)

            //****************************************************************************80
            // 
            //  Purpose:
            //
            //    PERM_LEX_RANK computes the lexicographic rank of a permutation.
            // 
            //  Discussion:
            // 
            //    The original code altered the input permutation.
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
            //    Input, int N, the number of values being permuted.
            //    N must be positive.
            // 
            //    Input, int P[N], describes the permutation.
            //    P[I] is the item which is permuted into the I-th place
            //    by the permutation.
            // 
            //    Output, int PERM_LEX_RANK, the rank of the permutation.
            // 
        {
            bool check;
            int i;
            int j;
            int[] pcopy;
            int rank;
            // 
            //  Check.
            // 
            check = perm_check(n, p);

            if (!check)
            {
                Console.WriteLine("");
                Console.WriteLine("PERM_LEX_RANK - Fatal error!");
                Console.WriteLine("  Permutation is illegal.");
                return (1);
            }

            rank = 0;
            pcopy = new int[n];

            for (i = 0; i < n; i++)
            {
                pcopy[i] = p[i];
            }

            for (j = 0; j < n; j++)
            {
                rank = rank + (pcopy[j] - 1) * typeMethods.i4_factorial(n - 1 - j);
                for (i = j + 1; i < n; i++)
                {
                    if (pcopy[j] < pcopy[i])
                    {
                        pcopy[i] = pcopy[i] - 1;
                    }
                }
            }

            return rank;
        }

        public static void perm_lex_successor(int n, ref int[] p, ref int rank )

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    PERM_LEX_SUCCESSOR computes the lexicographic permutation successor.
        // 
        //  Example:
        // 
        //    RANK  Permutation
        // 
        //       0  1 2 3 4
        //       1  1 2 4 3
        //       2  1 3 2 4
        //       3  1 3 4 2
        //       4  1 4 2 3
        //       5  1 4 3 2
        //       6  2 1 3 4
        //       ...
        //      23  4 3 2 1
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
        //    Input, int N, the number of values being permuted.
        //    N must be positive.
        // 
        //    Input/output, int P[N], describes the permutation.
        //    P(I) is the item which is permuted into the I-th place
        //    by the permutation.
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
            int j;
            int temp;
            // 
            //  Return the first element.
            // 
            if (rank == -1)
            {
                typeMethods.i4vec_indicator1(n, ref p);
                rank = 0;
                return;
            }

            // 
            //  Check.
            // 
            check = perm_check(n, p);

            if (!check)
            {
                Console.WriteLine("");
                Console.WriteLine("PERM_LEX_SUCCESSOR - Fatal error!");
                Console.WriteLine("  Permutation is illegal.");
                return;
            }

            // 
            //  Seek I, the highest index for which the next element is bigger.
            // 
            i = n - 1;

            for (;;)
            {
                if (i <= 0)
                {
                    break;
                }

                if (p[i - 1] <= p[i])
                {
                    break;
                }

                i = i - 1;
            }

            // 
            //  If no I could be found, then we have reach the final permutation,
            //  N, N-1, ..., 2, 1.  Time to start over again.
            // 
            if (i == 0)
            {
                typeMethods.i4vec_indicator1(n, ref p);
                rank = 0;
            }
            else
            {
                // 
                //  Otherwise, look for the the highest index after I whose element
                //  is bigger than I''s.  We know that I+1 is one such value, so the
                //  loop will never fail.
                // 
                j = n;
                while (p[j - 1] < p[i - 1])
                {
                    j = j - 1;
                }

                // 
                //  Interchange elements I and J.
                // 
                temp = p[i - 1];
                p[i - 1] = p[j - 1];
                p[j - 1] = temp;
                // 
                //  Reverse the elements from I+1 to N.
                // 
                typeMethods.i4vec_reverse(n - i, ref p, aIndex: + i);

                rank = rank + 1;
            }
        }

        public static int[] perm_lex_unrank(int rank, int n)

            //****************************************************************************80
            // 
            //  Purpose:
            //
            //    PERM_LEX_UNRANK computes the permutation of given lexicographic rank.
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
            //    Input, int RANK, the rank of the permutation.
            // 
            //    Input, int N, the number of values being permuted.
            //    N must be positive.
            // 
            //    Output, int PERM_LEX_UNRANK[N], describes the permutation.
            // 
        {
            int d;
            int i;
            int j;
            int nperm;
            int[] p;
            int rank_copy;
            // 
            //  Check.
            // 
            if (n < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("PERM_LEX_UNRANK - Fatal error!");
                Console.WriteLine("  Input N is illegal.");
                return null;
            }

            nperm = perm_enum(n);

            if (rank < 0 || nperm < rank)
            {
                Console.WriteLine("");
                Console.WriteLine("PERM_LEX_UNRANK - Fatal error!");
                Console.WriteLine("  The input rank is illegal.");
                return null;
            }

            rank_copy = rank;

            p = new int[n];

            p[n - 1] = 1;

            for (j = 1; j <= n - 1; j++)
            {
                d = (rank_copy % typeMethods.i4_factorial(j + 1)) / typeMethods.i4_factorial(j);
                rank_copy = rank_copy - d * typeMethods.i4_factorial(j);
                p[n - j - 1] = d + 1;

                for (i = n - j + 1; i <= n; i++)
                {
                    if (d < p[i - 1])
                    {
                        p[i - 1] = p[i - 1] + 1;
                    }
                }
            }

            return p;
        }

        public static int[] perm_mul(int n, int[] p, int[] q )

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    PERM_MUL computes the product of two permutations.
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
        //    Donald Kreher, Douglas Simpson,inson,
        //    Combinatorial Algorithms,
        //    CRC Press, 1998,
        //    ISBN: 0-8493-3988-X,
        //    LC: QA164.K73.
        // 
        //  Parameters:
        // 
        //    Input, int N, the number of values being permuted.
        //    N must be positive.
        // 
        //    Input, int P[N], Q[N], describes the permutation factors.
        // 
        //    Output, int PERMN_MUL[N], the product permutation P * Q.
        //    R(I) = P(Q(I)).
        // 
        {
            bool check;
            int i;
            int[] r;
            // 
            //  Check.
            // 
            check = perm_check(n, p);

            if (!check)
            {
                Console.WriteLine("");
                Console.WriteLine("PERM_MUL - Fatal error!");
                Console.WriteLine("  Permutation is illegal.");
                return null;
            }

            check = perm_check(n, q);

            if (!check)
            {
                Console.WriteLine("");
                Console.WriteLine("PERM_MUL - Fatal error!");
                Console.WriteLine("  Permutation is illegal.");
                return null;
            }

            // 
            //  Use a temporary vector for the result, to avoid problems if
            //  some arguments are actually identified.
            // 
            r = new int[n];

            for (i = 0; i < n; i++)
            {
                r[i] = p[q[i] - 1];
            }

            return r;
        }

        public static int perm_parity(int n, int[] p)

            //****************************************************************************80
            // 
            //  Purpose:
            //
            //    PERM_PARITY computes the parity of a permutation.
            // 
            //  Discussion:
            // 
            //    The routine requires the use of a temporary array.
            // 
            //    A permutation is called "even" or "odd", depending on whether
            //    it is equivalent to an even or odd number of pairwise
            //    transpositions.  This is known as the "parity" of the
            //    permutation.
            // 
            //    The "sign" of a permutation is +1 if it has even parity,
            //    and -1 if it has odd parity.
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
            //    Input, int N, the number of values being permuted.
            //    N must be positive.
            // 
            //    Input, int P[N], describes the permutation.
            //    P(I) is the item which is permuted into the I-th place
            //    by the permutation.
            // 
            //    Output, int PERM_PARITY, the parity of the permutation.
            //    0, the permutation has even parity.
            //    1, the permutation has odd parity.
            // 
        {
            int[] a;
            int c;
            bool check;
            int i;
            int j;
            int parity;
            // 
            //  Check.
            // 
            check = perm_check(n, p);

            if (!check)
            {
                Console.WriteLine("");
                Console.WriteLine("PERM_PARITY - Fatal error!");
                Console.WriteLine("  Permutation is illegal.");
                return (1);
            }

            a = new int[n];

            for (i = 0; i < n; i++)
            {
                a[i] = 0;
            }

            c = 0;

            for (j = 1; j <= n; j++)
            {
                if (a[j - 1] == 0)
                {
                    c = c + 1;
                    a[j - 1] = 1;
                    i = j;

                    while (p[i - 1] != j)
                    {
                        i = p[i - 1];
                        a[i - 1] = 1;
                    }
                }
            }

            parity = (n - c) % 2;
            
            return parity;
        }

        public static void perm_print(int n, int[] p, string title )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PERM_PRINT prints a permutation.
        //
        //  Discussion:
        //
        //    The permutation is assumed to be zero-based.
        //
        //  Example:
        //
        //    Input:
        //
        //      P = 6 1 2 0 4 2 5
        //
        //    Printed output:
        //
        //      "This is the permutation:"
        //
        //      0 1 2 3 4 5 6
        //      6 1 2 0 4 2 5
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 May 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of objects permuted.
        //
        //    Input, int P[N], the permutation, in standard index form.
        //
        //    Input, string TITLE, a title.
        //    If no title is supplied, then only the permutation is printed.
        //
        {
            int i;
            int ihi;
            int ilo;
            int inc = 20;

            if (typeMethods.s_len_trim(title) != 0)
            {
                Console.WriteLine("");
                Console.WriteLine(title + "");

                for (ilo = 0; ilo < n; ilo = ilo + inc)
                {
                    ihi = ilo + inc;
                    if (n < ihi)
                    {
                        ihi = n;
                    }

                    Console.WriteLine("");
                    string cout = "  ";
                    for (i = ilo; i < ihi; i++)
                    {
                        cout += i.ToString().PadLeft(4);
                    }

                    Console.WriteLine(cout);
                    cout = "  ";
                    for (i = ilo; i < ihi; i++)
                    {
                        cout += p[i].ToString().PadLeft(4);
                    }

                    Console.WriteLine(cout);
                }
            }
            else
            {
                for (ilo = 0; ilo < n; ilo = ilo + inc)
                {
                    ihi = ilo + inc;
                    if (n < ihi)
                    {
                        ihi = n;
                    }

                    string cout = "  ";
                    for (i = ilo; i < ihi; i++)
                    {
                        cout += p[i].ToString().PadLeft(4);
                    }

                    Console.WriteLine(cout);
                }
            }

            return;
        }

        public static int[] perm_random(int n, ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM_RANDOM selects a random permutation of 1, ..., N.
            //
            //  Discussion:
            //
            //    The algorithm is known as the Fisher-Yates or Knuth shuffle.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    11 January 2016
            //
            //  Author:
            //
            //    John Burkardt.
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
            //    Input, int N, the number of objects to be permuted.
            //
            //    Input/output, int &SEED, a seed for the random number generator.
            //
            //    Output, int PERM_RANDOM[N], a permutation of ( 1, 2, ..., N ), in standard
            //    index form.
            //
        {
            int i;
            int j;
            int[] p;
            int t;

            p = new int[n];

            for (i = 0; i < n; i++)
            {
                p[i] = i + 1;
            }

            for (i = 0; i < n - 1; i++)
            {
                j = UniformRNG.i4_uniform_ab(i, n - 1, ref seed);

                t = p[i];
                p[i] = p[j];
                p[j] = t;
            }

            return p;
        }

        public static int perm_tj_rank(int n, int[] p)

            //****************************************************************************80
            // 
            //  Purpose:
            //
            //    PERM_TJ_RANK computes the Trotter-Johnson rank of a permutation.
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
            //    Input, int N, the number of values being permuted.
            //    N must be positive.
            // 
            //    Input, int P[N], describes the permutation.
            //    P(I) is the item which is permuted into the I-th place
            //    by the permutation.
            // 
            //    Output, int PERM_TJ_RANK, the rank of the permutation.
            // 
        {
            bool check;
            int i;
            int j;
            int k;
            int rank;
            // 
            //  Check.
            // 
            check = perm_check(n, p);

            if (!check)
            {
                Console.WriteLine("");
                Console.WriteLine("PERM_TJ_RANK - Fatal error!");
                Console.WriteLine("  Permutation is illegal.");
                return(1);
            }

            rank = 0;

            for (j = 2; j <= n; j++)
            {
                k = 1;
                i = 1;

                while (p[i - 1] != j)
                {
                    if (p[i - 1] < j)
                    {
                        k = k + 1;
                    }

                    i = i + 1;
                }

                if ((rank % 2) == 0)
                {
                    rank = j * rank + j - k;
                }
                else
                {
                    rank = j * rank + k - 1;
                }
            }

            return rank;
        }

        public static void perm_tj_successor(int n, ref int[] p, ref int rank )

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    PERM_TJ_SUCCESSOR computes the Trotter-Johnson permutation successor.
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
        //    Input, int N, the number of values being permuted.
        //    N must be positive.
        // 
        //    Input/output, int P[N], describes the permutation.
        //    P(I) is the item which is permuted into the I-th place
        //    by the permutation.
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
            int d;
            bool done;
            int i;
            int m;
            int par;
            int[] q;
            int st;
            int temp;
            // 
            //  Return the first element.
            // 
            if (rank == -1)
            {
                typeMethods.i4vec_indicator1(n, ref p);
                rank = 0;
                return;
            }

            // 
            //  Check.
            // 
            check = perm_check(n, p);

            if (!check)
            {
                Console.WriteLine("");
                Console.WriteLine("PERM_TJ_SUCCESSOR - Fatal error!");
                Console.WriteLine("  Permutation is illegal.");
                return;
            }

            q = new int[n];

            st = 0;
            for (i = 0; i < n; i++)
            {
                q[i] = p[i];
            }

            done = false;
            m = n;

            while (1 < m && !done)
            {
                d = 1;
                while (q[d - 1] != m)
                {
                    d = d + 1;
                }

                for (i = d; i < m; i++)
                {
                    q[i - 1] = q[i];
                }

                par = perm_parity(m - 1, q);

                if (par == 1)
                {
                    if (d == m)
                    {
                        m = m - 1;
                    }
                    else
                    {
                        temp = p[st + d - 1];
                        p[st + d - 1] = p[st + d];
                        p[st + d] = temp;
                        done = true;
                    }
                }
                else
                {
                    if (d == 1)
                    {
                        m = m - 1;
                        st = st + 1;
                    }
                    else
                    {
                        temp = p[st + d - 1];
                        p[st + d - 1] = p[st + d - 2];
                        p[st + d - 2] = temp;
                        done = true;
                    }
                }
            }

            // 
            //  Last element was input.  Return first one.
            // 
            if (m == 1)
            {
                typeMethods.i4vec_indicator1(n, ref p);
                rank = 0;
                return;
            }

            rank = rank + 1;
        }

        public static int[] perm_tj_unrank(int rank, int n)

            //****************************************************************************80
            // 
            //  Purpose:
            //
            //    PERM_TJ_UNRANK computes the permutation of given Trotter-Johnson rank.
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
            //    Input, int RANK, the rank of the permutation.
            // 
            //    Input, int N, the number of values being permuted.
            //    N must be positive.
            // 
            //    Output, int PERM_TJ_UNRANK[N], describes the permutation.
            // 
        {
            int i;
            int j;
            int k;
            int jhi;
            int nperm;
            int[] p;
            int r1;
            int r2;
            // 
            //  Check.
            // 
            if (n < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("PERM_TJ_UNRANK - Fatal error!");
                Console.WriteLine("  Input N is illegal.");
                return null;
            }

            nperm = perm_enum(n);

            if (rank < 0 || nperm < rank)
            {
                Console.WriteLine("");
                Console.WriteLine("PERM_TJ_UNRANK - Fatal error!");
                Console.WriteLine("  The input rank is illegal.");
                return null;
            }

            p = new int[n];

            p[0] = 1;
            r2 = 0;

            for (j = 2; j <= n; j++)
            {
                // 
                //  Replace this ratio of factorials!
                // 
                r1 = (rank * typeMethods.i4_factorial(j)) / typeMethods.i4_factorial(n);
                k = r1 - j * r2;

                if ((r2 % 2) == 0)
                {
                    jhi = j - k;
                }
                else
                {
                    jhi = k + 1;
                }

                for (i = j - 1; jhi <= i; i--)
                {
                    p[i] = p[i - 1];
                }

                p[jhi - 1] = j;

                r2 = r1;
            }

            return p;
        }

        public static void perm_to_cycle(int n, int[] p, ref int ncycle, ref int[] t, ref int[] index )

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    PERM_TO_CYCLE converts a permutation from array to cycle form.
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
        //    Input, int N, the number of values being permuted.
        //    N must be positive.
        // 
        //    Input, int P[N], describes the permutation using a
        //    single array.  For each index I, I -> P(I).
        // 
        //    Output, int &NCYCLE, the number of cycles.
        //    1 <= NCYCLE <= N.
        // 
        //    Output, int T[N], INDEX[N], describes the permutation
        //    as a collection of NCYCLE cycles.  The first cycle is
        //    T(1) -> T(2) -> ... -> T(INDEX(1)) -> T(1).
        // 
        {
            bool check;
            int i;
            int j;
            int nset;
            // 
            //  Check.
            // 
            check = perm_check(n, p);

            if (!check)
            {
                Console.WriteLine("");
                Console.WriteLine("PERM_TO_CYCLE - Fatal error!");
                Console.WriteLine("  Permutation is illegal.");
                return;
            }

            // 
            //  Initialize.
            // 
            ncycle = 0;
            for (i = 0; i < n; i++)
            {
                index[i] = 0;
            }

            for (i = 0; i < n; i++)
            {
                t[i] = 0;
            }

            nset = 0;
            // 
            //  Find the next unused entry.
            //
            for (i = 1; i <= n; i++)
            {
                if (0 < p[i - 1])
                {
                    ncycle = ncycle + 1;
                    index[ncycle - 1] = 1;

                    nset = nset + 1;
                    t[nset - 1] = p[i - 1];
                    p[i - 1] = -p[i - 1];

                    for (;;)
                    {
                        j = t[nset - 1];

                        if (p[j - 1] < 0)
                        {
                            break;
                        }

                        index[ncycle - 1] = index[ncycle - 1] + 1;

                        nset = nset + 1;
                        t[nset - 1] = p[j - 1];
                        p[j - 1] = -p[j - 1];
                    }
                }
            }

            // 
            //  If no unused entries remain, we are done.
            //  Restore the sign of the permutation and return.
            //
            for (i = 0; i < n; i++)
            {
                p[i] = -p[i];
            }
        }
    }
}