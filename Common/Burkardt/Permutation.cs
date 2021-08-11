using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt
{
    public static partial class Permutation
    {
        public static void euler_row ( int n, ref int[] ieuler )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    EULER_ROW returns the N-th row of Euler's triangle.
            //
            //  Discussion:
            //
            //    E(N,K) counts the number of permutations of the N digits that have
            //    exactly K "ascents", that is, K places where the Ith digit is
            //    less than the (I+1)th digit.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    10 June 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the row of Euler's triangle desired.
            //
            //    Output, int IEULER[N+1], the N-th row of Euler's
            //    triangle, IEULER[K] contains the value of E(N,K).  Note
            //    that IEULER[0] should be 1 and IEULER[N] should be 0.
            //
        {
            int irow;
            int k;

            ieuler[0] = 1;

            if ( 0 < n )
            {
                ieuler[1] = 0;

                for ( irow = 2; irow <= n; irow++ )
                {
                    ieuler[irow] = 0;

                    for ( k = irow-1; 1 <= k; k-- )
                    {
                        ieuler[k] = ( k + 1 ) * ieuler[k] + ( irow - k ) * ieuler[k-1];
                    }
                    ieuler[0] = 1;
                }
            }
            return;
        }
        
        public static void inversion_to_perm0 ( int n, int[] ins, ref int[] p )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INVERSION_TO_PERM0: inversion sequence to permutation (0,...,N-1).
        //
        //  Discussion:
        //
        //    For a given permutation P acting on objects 0 through N-1, the
        //    inversion sequence INS is defined as:
        //
        //      INS(1) = 0
        //      INS(I) = number of values J < I for which P(I) < P(J).
        //
        //  Example:
        //
        //    Input:
        //
        //      ( 0, 0, 2, 1, 3 )
        //
        //    Output:
        //
        //      ( 2, 4, 0, 3, 1 )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 June 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Dennis Stanton, Dennis White,
        //    Constructive Combinatorics,
        //    Springer, 1986,
        //    ISBN: 0387963472,
        //    LC: QA164.S79.
        //
        //  Parameters:
        //
        //    Input, int N, the number of objects being permuted.
        //
        //    Input, int INS[N], the inversion sequence of a permutation.
        //    It must be the case that 0 <= INS(I) < I for 1 <= I <= N.
        //
        //    Output, int P[N], the permutation.
        //
        {
            int i;
            int itemp;
            int j;

            typeMethods.i4vec_indicator0 ( n, ref p );

            for ( i = n - 1; 1 <= i; i-- )
            {
                itemp = p[i-ins[i]];

                for ( j = i-ins[i]; j <= i-1; j++ )
                {
                    p[j] = p[j+1];
                }

                p[i] = itemp;
            }
        }
        
        public static void perm_ascend(int n, int[] a, ref int length, ref int[] sub)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM_ASCEND computes the longest ascending subsequence of permutation.
            //
            //  Discussion:
            //
            //    Although this routine is intended to be applied to a permutation,
            //    it will work just as well for an arbitrary vector.
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
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of the permutation.
            //
            //    Input, int A[N], the permutation to be examined.
            //
            //    Output, int &LENGTH, the length of the longest increasing subsequence.
            //
            //    Output, int SUB[N], contains in entries 1 through LENGTH
            //    a longest increasing subsequence of A.
            //
        {
            int i;
            int j;
            int k;
            int[] top;
            int[] top_prev;

            if (n <= 0)
            {
                length = 0;
                return;
            }

            top = new int[n];
            for (i = 0; i < n; i++)
            {
                top[i] = 0;
            }

            top_prev = new int[n];
            for (i = 0; i < n; i++)
            {
                top_prev[i] = 0;
            }

            for (i = 0; i < n; i++)
            {
                sub[i] = 0;
            }

            length = 0;

            for (i = 1; i <= n; i++)
            {
                k = 0;

                for (j = 1; j <= length; j++)
                {
                    if (a[i - 1] <= a[top[j - 1] - 1])
                    {
                        k = j;
                        break;
                    }
                }

                if (k == 0)
                {
                    length = length + 1;
                    k = length;
                }

                top[k - 1] = i;

                if (1 < k)
                {
                    top_prev[i - 1] = top[k - 2];
                }
                else
                {
                    top_prev[i - 1] = 0;
                }
            }

            j = top[length - 1];
            sub[length - 1] = a[j - 1];

            for (i = length - 1; 1 <= i; i--)
            {
                j = top_prev[j - 1];
                sub[i - 1] = a[j - 1];
            }
        }

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

        public static int perm_fixed_enum(int n, int m)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM_FIXED_ENUM enumerates the permutations of N objects with M fixed.
            //
            //  Discussion:
            //
            //    A permutation of N objects with M fixed is a permutation in which
            //    exactly M of the objects retain their original positions.  If
            //    M = 0, the permutation is a "derangement".  If M = N, the
            //    permutation is the identity.
            //
            //    The formula is:
            //
            //      F(N,M) = ( N! / M! ) * ( 1 - 1/1! + 1/2! - 1/3! ... 1/(N-M)! )
            //             = COMB(N,M) * D(N-M)
            //
            //    where D(N-M) is the number of derangements of N-M objects.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    02 June 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of objects to be permuted.
            //    N should be at least 1.
            //
            //    Input, int M, the number of objects that retain their
            //    position.  M should be between 0 and N.
            //
            //    Output, int PERM_FIXED_ENUM, the number of derangements of N objects
            //    in which M objects retain their positions.
            //
        {
            int fnm;

            if (n <= 0)
            {
                fnm = 1;
            }
            else if (m < 0)
            {
                fnm = 0;
            }
            else if (n < m)
            {
                fnm = 0;
            }
            else if (m == n)
            {
                fnm = 1;
            }
            else if (n == 1)
            {
                if (m == 1)
                {
                    fnm = 1;
                }
                else
                {
                    fnm = 0;
                }
            }
            else
            {
                fnm = typeMethods.i4_choose(n, m) * Derange.derange_enum(n - m);
            }

            return fnm;
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

        public static void perm_lex_successor(int n, ref int[] p, ref int rank)

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
                typeMethods.i4vec_reverse(n - i, ref p, aIndex: +i);

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

        public static int[] perm_mul(int n, int[] p, int[] q)

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

        public static void perm_print(int n, int[] p, string title)

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
                return (1);
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

        public static void perm_tj_successor(int n, ref int[] p, ref int rank)

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

        public static void perm_to_cycle(int n, int[] p, ref int ncycle, ref int[] t, ref int[] index)

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

        public static int multiperm_enum(int n, int k, int[] counts)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MULTIPERM_ENUM enumerates multipermutations.
            //
            //  Discussion:
            //
            //    A multipermutation is a permutation of objects, some of which are
            //    identical.
            //
            //    While there are 6 permutations of the distinct objects A,B,C, there
            //    are only 3 multipermutations of the objects A,B,B.
            //
            //    In general, there are N! permutations of N distinct objects, but
            //    there are N! / ( ( M1! ) ( M2! ) ... ( MK! ) ) multipermutations
            //    of N objects, in the case where the N objects consist of K
            //    types, with M1 examples of type 1, M2 examples of type 2 and so on,
            //    and for which objects of the same type are indistinguishable.
            //
            //  Example:
            //
            //    Input:
            //
            //      N = 5, K = 3, COUNTS = (/ 1, 2, 2 /)
            //
            //    Output:
            //
            //      Number = 30
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 July 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of items in the multipermutation.
            //
            //    Input, int K, the number of types of items.
            //    1 <= K.  Ordinarily, K <= N, but we allow any positive K, because
            //    we also allow entries in COUNTS to be 0.
            //
            //    Input, int COUNTS[K], the number of items of each type.
            //    0 <= COUNTS(1:K) <= N and sum ( COUNTS(1:K) ) = N.
            //
            //    Output, int MULTIPERM_ENUM, the number of multipermutations.
            //
        {
            int i;
            int j;
            int number;
            int sum;
            int top;

            if (n < 0)
            {
                number = -1;
                return number;
            }

            if (n == 0)
            {
                number = 1;
                return number;
            }

            if (k < 1)
            {
                number = -1;
                return number;
            }

            for (i = 0; i < k; i++)
            {
                if (counts[i] < 0)
                {
                    number = -1;
                    return number;
                }
            }

            sum = 0;
            for (i = 0; i < k; i++)
            {
                sum = sum + counts[i];
            }

            if (sum != n)
            {
                number = -1;
                return number;
            }

            //
            //  Ready for computation.
            //  By design, the integer division should never have a remainder.
            //
            top = 0;
            number = 1;

            for (i = 0; i < k; i++)
            {
                for (j = 1; j <= counts[i]; j++)
                {
                    top = top + 1;
                    number = (number * top) / j;
                }
            }

            return number;
        }

        public static void multiperm_next(int n, ref int[] a, ref bool more)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MULTIPERM_NEXT returns the next multipermutation.
            //
            //  Discussion:
            //
            //    To begin the computation, the user must set up the first multipermutation.
            //    To compute ALL possible multipermutations, this first permutation should
            //    list the values in ascending order.
            //
            //    The routine will compute, one by one, the next multipermutation,
            //    in lexicographical order.  On the call after computing the last 
            //    multipermutation, the routine will return MORE = FALSE (and will
            //    reset the multipermutation to the FIRST one again.)
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 March 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of items in the multipermutation.
            //
            //    Input/output, int A[N]; on input, the current multipermutation.
            //    On output, the next multipermutation.
            //
            //    Output, bool &MORE, is TRUE if the next multipermutation
            //    was computed, or FALSE if no further multipermutations could
            //    be computed.
            //
        {
            int i;
            int m;
            int temp;
            //
            //  Step 1:
            //  Find M, the last location in A for which A(M) < A(M+1).
            //
            m = 0;
            for (i = 1; i <= n - 1; i++)
            {
                if (a[i - 1] < a[i])
                {
                    m = i;
                }
            }

            //
            //  Step 2:
            //  If no M was found, we've run out of multipermutations.
            //
            if (m == 0)
            {
                more = false;
                typeMethods.i4vec_sort_heap_a(n, ref a);
                return;
            }
            else
            {
                more = true;
            }

            //
            //  Step 3:
            //  Ascending sort A(M+1:N).
            //
            if (m + 1 < n)
            {
                typeMethods.i4vec_sort_heap_a(n - m, ref a, aIndex: + m);
            }

            //
            //  Step 4:
            //  Locate the first larger value after A(M).
            //
            i = 1;
            for (;;)
            {
                if (a[m - 1] < a[m + i - 1])
                {
                    break;
                }

                i = i + 1;
            }

            //
            //  Step 5:
            //  Interchange A(M) and the next larger value.
            //
            temp = a[m - 1];
            a[m - 1] = a[m + i - 1];
            a[m + i - 1] = temp;
        }

        public static int perm0_break_count(int n, int[] p)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM0_BREAK_COUNT counts the breaks in a permutation of (0,...,N-1).
            //
            //  Discussion:
            //
            //    We begin with a permutation of order N.  We prepend an element
            //    labeled "-1" and append an element labeled "N".  There are now
            //    N+1 pairs of neighbors.  A "break" is a pair of neighbors whose
            //    value differs by more than 1.  
            //
            //    The identity permutation has a break count of 0.  The maximum
            //    break count is N+1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    15 June 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of the permutation.
            //
            //    Input, int P[N], a permutation, in standard index form.
            //
            //    Output, int PERM0_BREAK_COUNT, the number of breaks in the permutation.
            //
        {
            int i;
            int value;

            value = 0;
            //
            //  Make sure the permutation is a legal one.
            //  (This is not an efficient way to do so!)
            //
            if (!perm0_check(n, p))
            {
                Console.WriteLine("");
                Console.WriteLine("PERM0_BREAK_COUNT - Fatal error!");
                Console.WriteLine("  PERM0_CHECK rejects permutation.");
                return 1;
            }

            if (p[0] != 0)
            {
                value = value + 1;
            }

            for (i = 1; i <= n - 1; i++)
            {
                if (Math.Abs(p[i] - p[i - 1]) != 1)
                {
                    value = value + 1;
                }
            }

            if (p[n - 1] != n - 1)
            {
                value = value + 1;
            }

            return value;
        }

        public static bool perm0_check(int n, int[] p)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM0_CHECK checks a permutation of ( 0, ..., N-1 ).
            //
            //  Discussion:
            //
            //    The routine verifies that each of the integers from 0 to
            //    to N-1 occurs among the N entries of the permutation.
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
            //    Input, int N, the number of entries.
            //
            //    Input, int P[N], the array to check.
            //
            //    Output, bool PERM0_CHECK, is 
            //    TRUE if P is a legal permutation of 0,...,N-1.
            //    FALSE if P is not a legal permuation of 0,...,N-1.
            //
        {
            bool check;
            int location;
            int value;

            check = true;

            for (value = 0; value < n; value++)
            {
                check = false;

                for (location = 0; location < n; location++)
                {
                    if (p[location] == value)
                    {
                        check = true;
                        break;
                    }
                }

                if (!check)
                {
                    Console.WriteLine("");
                    Console.WriteLine("PERM0_CHECK - Warning!");
                    Console.WriteLine("  Permutation is missing value " + value + "");
                    break;
                }

            }

            return check;
        }

        public static void perm0_cycle(int n, int[] p, ref int isgn, ref int ncycle, int iopt )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PERM0_CYCLE analyzes a permutation of (0,...,N-1).
        //
        //  Discussion:
        //
        //    The routine will count cycles, find the sign of a permutation,
        //    and tag a permutation.
        //
        //  Example:
        //
        //    Input:
        //
        //      N = 9
        //      IOPT = 1
        //      P = 2, 3, 9, 6, 7, 8, 5, 4, 1
        //
        //    Output:
        //
        //      NCYCLE = 3
        //      ISGN = +1
        //      P = -2, 3, 9, -6, -7, 8, 5, 4, 1
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
        //    Input, int N, the number of objects being permuted.
        //
        //    Input/output, int P[N].  On input, P describes a
        //    permutation, in the sense that entry I is to be moved to P[I].
        //    If IOPT = 0, then P will not be changed by this routine.
        //    If IOPT = 1, then on output, P will be "tagged".  That is,
        //    one element of every cycle in P will be negated.  In this way,
        //    a user can traverse a cycle by starting at any entry I1 of P
        //    which is negative, moving to I2 = ABS(P[I1]), then to
        //    P[I2], and so on, until returning to I1.
        //
        //    Output, int &ISGN, the "sign" of the permutation, which is
        //    +1 if the permutation is even, -1 if odd.  Every permutation
        //    may be produced by a certain number of pairwise switches.
        //    If the number of switches is even, the permutation itself is
        //    called even.
        //
        //    Output, int &NCYCLE, the number of cycles in the permutation.
        //
        //    Input, int IOPT, requests tagging.
        //    0, the permutation will not be tagged.
        //    1, the permutation will be tagged.
        //
        {
            int i;
            int i1;
            int i2;
            int is_;

            if (!perm0_check(n, p))
            {
                Console.WriteLine("");
                Console.WriteLine("PERM0_CYCLE - Fatal error!");
                Console.WriteLine("  PERM0_CHECK rejects permutation.");
                return;
            }

            //
            //  Increment.
            //
            for (i = 0; i < n; i++)
            {
                p[i] = p[i] + 1;
            }

            is_ = 1;
            ncycle = n;

            for (i = 1; i <= n; i++)
            {
                i1 = p[i - 1];

                while (i < i1)
                {
                    ncycle = ncycle - 1;
                    i2 = p[i1 - 1];
                    p[i1 - 1] = -i2;
                    i1 = i2;
                }

                if (iopt != 0)
                {
                    is_ = -typeMethods.i4_sign(p[i - 1]);
                }

                p[i - 1] = Math.Abs(p[i - 1]) * typeMethods.i4_sign( is_ );
            }

            isgn = 1 - 2 * ((n - ncycle) % 2);
            //
            //  Decrement.
            //
            for (i = 0; i < n; i++)
            {
                p[i] = p[i] - 1;
            }

            return;
        }

        public static int perm0_distance(int n, int[] a, int[] b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PERM0_DISTANCE computes the distance of two permutations of (0,...,N-1).
        //
        //  Discussion:
        //
        //    The distance is known as the Ulam metric.
        //
        //    If we let N be the order of the permutations A and B, and L(P) be
        //    the length of the longest ascending subsequence of a permutation P,
        //    then the Ulam metric distance between A and B is
        //
        //      N - L ( A * inverse ( B ) ).
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
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the permutation.
        //
        //    Input, int A[N], B[N], the permutations to be examined.
        //
        //    Output, int PERM0_DISTANCE, the Ulam metric distance between A and B.
        //
        {
            int[] binv;
            int[] c;
            int length = 0;
            int[] sub;
            int value;

            c = new int[n];
            sub = new int[n];

            binv = perm0_inverse(n, b);

            perm0_mul(n, a, binv, ref c);

            perm_ascend(n, c, ref length, ref sub);

            value = n - length;

            return value;
        }

        public static void perm0_free(int npart, int[] ipart, int nfree, ref int[] ifree )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PERM0_FREE reports unused items in a partial permutation of (0,...,N-1).
        //
        //  Discussion:
        //
        //    It is assumed that the N objects being permuted are the integers
        //    from 0 to N-1, and that IPART contains a "partial" permutation, that
        //    is, the NPART entries of IPART represent the beginning of a
        //    permutation of all N items.
        //
        //    The routine returns in IFREE the items that have not been used yet.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 June 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NPART, the number of entries in IPART.  NPART may be 0.
        //
        //    Input, int IPART[NPART], the partial permutation, which should
        //    contain, at most once, some of the integers between 1 and
        //    NPART+NFREE.
        //
        //    Input, int NFREE, the number of integers that have not been
        //    used in IPART.  This is simply N - NPART.  NFREE may be zero.
        //
        //    Output, int IFREE[NFREE], the integers between 1 and NPART+NFREE
        //    that were not used in IPART.
        //
        {
            int i;
            int j;
            int k;
            int match;
            int n;

            n = npart + nfree;

            if (npart < 0)
            {
                Console.WriteLine("");
                Console.WriteLine("PERM0_FREE - Fatal error!");
                Console.WriteLine("  NPART < 0.");
                return;
            }
            else if (npart == 0)
            {
                typeMethods.i4vec_indicator0(n, ref ifree);
            }
            else if (nfree < 0)
            {
                Console.WriteLine("");
                Console.WriteLine("PERM0_FREE - Fatal error!");
                Console.WriteLine("  NFREE < 0.");
                return;
            }
            else if (nfree == 0)
            {
                return;
            }
            else
            {
                k = 0;

                for (i = 0; i < n; i++)
                {
                    match = -1;

                    for (j = 0; j < npart; j++)
                    {
                        if (ipart[j] == i)
                        {
                            match = j;
                            break;
                        }
                    }

                    if (match == -1)
                    {
                        k = k + 1;

                        if (nfree < k)
                        {
                            Console.WriteLine("");
                            Console.WriteLine("PERM0_FREE - Fatal error!");
                            Console.WriteLine("  The partial permutation is illegal.");
                            Console.WriteLine("  Technically, because NFREE < K.");
                            Console.WriteLine("  N     = " + n + "");
                            Console.WriteLine("  NPART = " + npart + "");
                            Console.WriteLine("  NFREE = " + nfree + "");
                            Console.WriteLine("  K =     " + k + "");
                            Console.WriteLine("");
                            Console.WriteLine("  The partial permutation:");
                            Console.WriteLine("");
                            string cout = "";
                            for (i = 0; i < npart; i++)
                            {
                                cout += ipart[i].ToString().PadLeft(2) + "  ";
                            }

                            Console.WriteLine(cout);
                            return;
                        }

                        ifree[k - 1] = i;
                    }
                }
            }
        }

        public static int[] perm0_inverse(int n, int[] p1)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM0_INVERSE inverts a permutation of (0,...,N-1).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    08 June 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of objects being permuted.
            //
            //    Input, int P1[N], the permutation.
            //
            //    Output, int PERM0_INVERSE[N], the inverse permutation.
            //
        {
            int i;
            int i0;
            int i1;
            int i2;
            int is_;
            int[] p2;

            if (n <= 0)
            {
                Console.WriteLine("");
                Console.WriteLine("PERM0_INVERSE - Fatal error!");
                Console.WriteLine("  Input value of N = " + n + "");
                return null;
            }

            if (!perm0_check(n, p1))
            {
                Console.WriteLine("");
                Console.WriteLine("PERM0_INVERSE - Fatal error!");
                Console.WriteLine("  PERM0_CHECK rejects permutation.");
                return null;
            }

            p2 = new int[n];
            for (i = 0; i < n; i++)
            {
                p2[i] = p1[i] + 1;
            }

            is_ = 1;

            for (i = 1; i <= n; i++)
            {
                i1 = p2[i - 1];

                while (i < i1)
                {
                    i2 = p2[i1 - 1];
                    p2[i1 - 1] = -i2;
                    i1 = i2;
                }

                is_ = -typeMethods.i4_sign(p2[i - 1]);
                p2[i - 1] = Math.Abs(p2[i - 1]) * typeMethods.i4_sign( is_ );

            }

            for (i = 1; i <= n; i++)
            {
                i1 = -p2[i - 1];

                if (0 <= i1)
                {
                    i0 = i;

                    for (;;)
                    {
                        i2 = p2[i1 - 1];
                        p2[i1 - 1] = i0;

                        if (i2 < 0)
                        {
                            break;
                        }

                        i0 = i1;
                        i1 = i2;
                    }
                }
            }

            typeMethods.i4vec_decrement(n, ref p2);

            return p2;
        }

        public static void perm0_inverse2(int n, int[] p)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM0_INVERSE2 inverts a permutation of (0,...,N-1).
            //
            //  Discussion:
            //
            //    The routine needs no extra vector storage in order to compute the
            //    inverse of a permutation.
            //
            //    This feature might be useful if the permutation is large.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    12 June 2015
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
            //    Input, int N, the number of objects in the permutation.
            //
            //    Input/output, int P[N], the permutation, in standard index form.
            //    On output, the inverse permutation.
            //
        {
            int i;
            int ii;
            int j;
            int k;
            int m;

            if (!perm0_check(n, p))
            {
                Console.WriteLine("");
                Console.WriteLine("PERM0_INVERSE2 - Fatal error!");
                Console.WriteLine("  PERM0_CHECK rejects permutation.");
                return;
            }

            for (i = 0; i < n; i++)
            {
                p[i] = p[i] + 1;
            }

            for (ii = 1; ii <= n; ii++)
            {
                m = n + 1 - ii;
                i = p[m - 1];

                if (i < 0)
                {
                    p[m - 1] = -i;
                }
                else if (i != m)
                {
                    k = m;

                    for (;;)
                    {
                        j = p[i - 1];
                        p[i - 1] = -k;

                        if (j == m)
                        {
                            p[m - 1] = i;
                            break;
                        }

                        k = i;
                        i = j;
                    }
                }
            }

            for (i = 0; i < n; i++)
            {
                p[i] = p[i] - 1;
            }
        }

        public static int[] perm0_inverse3_new(int n, int[] p)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM0_INVERSE3 produces the inverse of a permutation of (0,...,N-1).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    12 June 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of items permuted.
            //
            //    Input, int P[N], a permutation.
            //
            //    Output, int PERM0_INVERSE3_NEW[N], the inverse permutation.
            //
        {
            int i;
            int[] p_inv;

            if (!perm0_check(n, p))
            {
                Console.WriteLine("");
                Console.WriteLine("PERM0_INVERSE3 - Fatal error!");
                Console.WriteLine("  PERM0_CHECK rejects permutation.");
                return null;
            }

            p_inv = new int[n];

            for (i = 0; i < n; i++)
            {
                p_inv[p[i]] = i;
            }

            return p_inv;
        }

        public static void perm0_lex_next(int n, int[] p, ref bool more )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PERM0_LEX_NEXT generates permutations of (0,...,N-1) in lexical order.
        //
        //  Example:
        //
        //    N = 3
        //
        //    1   0 1 2
        //    2   0 2 1
        //    3   1 0 2
        //    4   1 2 0
        //    5   2 0 1
        //    6   2 1 0
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
        //    Original FORTRAN77 version by Mok-Kong Shen.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Mok-Kong Shen,
        //    Algorithm 202: Generation of Permutations in Lexicographical Order,
        //    Communications of the ACM,
        //    Volume 6, September 1963, page 517.
        //
        //  Parameters:
        //
        //    Input, int N, the number of elements being permuted.
        //
        //    Input/output, int P[N], the permutation, in standard index form.
        //    The user should not alter the elements of Pbetween successive
        //    calls.  The routine uses the input value of P to determine
        //    the output value.
        //
        //    Input/output, bool &MORE.
        //    On the first call, the user should set MORE =FALSE which signals
        //    the routine to do initialization.
        //    On return, if MORE is TRUE then another permutation has been
        //    computed and returned, while if MORE is FALSE there are no more
        //    permutations.
        //
        {
            int j;
            int k;
            int temp;
            int u;
            int w;
            //
            //  Initialization.
            //
            if (!more)
            {
                typeMethods.i4vec_indicator0(n, ref p);
                more = true;
            }
            else
            {
                if (n <= 1)
                {
                    more = false;
                    return;
                }

                w = n;

                while (p[w - 1] < p[w - 2])
                {
                    if (w == 2)
                    {
                        more = false;
                        return;
                    }

                    w = w - 1;
                }

                u = p[w - 2];

                for (j = n; w <= j; j--)
                {
                    if (u < p[j - 1])
                    {
                        p[w - 2] = p[j - 1];
                        p[j - 1] = u;

                        for (k = 0; k <= (n - w - 1) / 2; k++)
                        {
                            temp = p[n - k - 1];
                            p[n - k - 1] = p[w + k - 1];
                            p[w + k - 1] = temp;
                        }

                        return;
                    }
                }
            }
        }

        public static void perm0_mul(int n, int[] p1, int[] p2, ref int[] p3 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PERM0_MUL "multiplies" two permutations of (0,...,N-1).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 June 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the permutations.
        //
        //    Input, int P1[N], P2[N], the permutations.
        //
        //    Output, int P3[N], the product permutation.
        //
        {
            int i;

            if (!perm0_check(n, p1))
            {
                Console.WriteLine("");
                Console.WriteLine("PERM0_MUL - Fatal error!");
                Console.WriteLine("  PERM0_CHECK rejects permutation.");
                return;
            }

            if (!perm0_check(n, p2))
            {
                Console.WriteLine("");
                Console.WriteLine("PERM0_MUL - Fatal error!");
                Console.WriteLine("  PERM0_CHECK rejects permutation.");
                return;
            }

            for (i = 0; i < n; i++)
            {
                p3[i] = p2[p1[i]];
            }
        }

        public static void perm0_next(int n, int[] p, ref bool more, ref bool even )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PERM0_NEXT computes permutations of (0,...,N-1), one at a time.
        //
        //  Discussion:
        //
        //    The routine is initialized by calling with MORE = TRUE, in which case
        //    it returns the identity permutation.
        //
        //    If the routine is called with MORE = FALSE, then the successor of the
        //    input permutation is computed.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 June 2015
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
        //    Input, int N, the number of objects being permuted.
        //
        //    Input/output, int P[N], the permutation, in standard index form.
        //    On the first call, the input value is unimportant.
        //    On subsequent calls, the input value should be the same
        //    as the output value from the previous call.  In other words, the
        //    user should just leave P alone.
        //    On output, contains the "next" permutation.
        //
        //    Input/output, bool &MORE.
        //    Set MORE = FALSE before the first call.
        //    MORE will be reset to TRUE and a permutation will be returned.
        //    Each new call produces a new permutation until
        //    MORE is returned FALSE.
        //
        //    Input/output, bool &EVEN.
        //    The input value of EVEN should simply be its output value from the
        //    previous call; (the input value on the first call doesn't matter.)
        //    On output, EVEN is TRUE if the output permutation is even, that is,
        //    involves an even number of transpositions.
        //
        {
            int i = 0;
            int i1;
            int ia = 0;
            int id;
            int is_;
            int j;
            int l = 0;
            int m;

            if (!more)
            {
                typeMethods.i4vec_indicator0(n, ref p);

                more = true;
                even = true;

                if (n == 1)
                {
                    more = false;
                    return;
                }

                if (p[n - 1] != 0 || p[0] != 1 + (n % 2))
                {
                    return;
                }

                for (i = 1; i <= n - 3; i++)
                {
                    if (p[i] != p[i - 1] + 1)
                    {
                        return;
                    }
                }

                more = false;
            }
            else
            {
                if (n == 1)
                {
                    p[0] = 0;
                    more = false;
                    return;
                }

                if (even)
                {
                    ia = p[0];
                    p[0] = p[1];
                    p[1] = ia;
                    even = false;

                    if (p[n - 1] != 0 || p[0] != 1 + (n % 2))
                    {
                        return;
                    }

                    for (i = 1; i <= n - 3; i++)
                    {
                        if (p[i] != p[i - 1] + 1)
                        {
                            return;
                        }
                    }

                    more = false;
                    return;
                }
                else
                {
                    more = false;

                        is_ = 0;

                    for (i1 = 2; i1 <= n; i1++)
                    {
                        ia = p[i1 - 1];
                        i = i1 - 1;
                        id = 0;

                        for (j = 1; j <= i; j++)
                        {
                            if (ia < p[j - 1])
                            {
                                id = id + 1;
                            }
                        }

                        is_ = id + is_;

                        if (id != i * ( is_  % 2) )
                        {
                            more = true;
                            break;
                        }
                    }

                    if (!more)
                    {
                        p[0] = 0;
                        return;
                    }
                }

                m = (( is_ +1 ) % 2 ) *(n + 1);

                for (j = 1; j <= i; j++)
                {
                    if (typeMethods.i4_sign(p[j - 1] - ia) != typeMethods.i4_sign(p[j - 1] - m))
                    {
                        m = p[j - 1];
                        l = j;
                    }
                }

                p[l - 1] = ia;
                p[i1 - 1] = m;
                even = true;
            }
        }

        public class PermNext2Data
        {
            public int[] active = null;
            public int[] idir = null;
            public int[] invers = null;
            
        }

        public static void perm0_next2(ref PermNext2Data data, int n, int[] p, ref bool done )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PERM0_NEXT2 generates permutations of (0,...,N-1).
        //
        //  Discussion:
        //
        //    The routine generates the permutations one at a time.  It uses a
        //    particular ordering of permutations, generating them from the first
        //    (which is the identity permutation) to the N!-th.  The same ordering
        //    is used by the routines PERM0_RANK and PERM0_UNRANK.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 October 2006
        //
        //  Author:
        //
        //    John Burkardt.
        //
        //  Reference:
        //
        //    Dennis Stanton, Dennis White,
        //    Constructive Combinatorics,
        //    Springer, 1986,
        //    ISBN: 0387963472,
        //    LC: QA164.S79.
        //
        //  Parameters:
        //
        //    Input, int N, the number of elements in the set to be permuted.
        //
        //    Input/output, int P[N], the permutation, in standard index form.
        //
        //    Input/output, bool &DONE.  The user should set the input value of
        //    DONE only once, before the first call to compute the permutations.
        //    The user should set DONE to TRUE, which signals the routine
        //    that it is to initialize itself.
        //    Thereafter, the routine will set DONE to FALSE and will
        //    compute a new permutation on each call.
        //    However, when there are no more permutations to compute, the
        //    routine will not return a new permutation, but instead will
        //    return DONE with the value TRUE.  At this point, all the
        //    permutations have been computed.
        //
        {
            int i;
            int j;
            int nactiv;
            //
            //  An input value of FALSE for DONE is assumed to mean a new
            //  computation is beginning.
            //
            if (done)
            {
                typeMethods.i4vec_indicator1(n, ref p);

                data.active = new int[n];

                data.idir = new int[n];

                data.invers = new int[n];

                for (i = 0; i < n; i++)
                {
                    data.invers[i] = p[i];
                }

                for (i = 0; i < n; i++)
                {
                    data.idir[i] = -1;
                }

                data.active[0] = 0;
                for (i = 1; i < n; i++)
                {
                    data.active[i] = 1;
                }

                //
                //  Set the DONE flag to FALSE, signifying there are more permutations
                //  to come.  Except, of course, that we must take care of the trivial case!
                //
                if (1 < n)
                {
                    done = false;
                }
                else
                {
                    done = true;
                }
            }
            //
            //  Otherwise, assume we are in a continuing computation
            //
            else
            {
                nactiv = 0;

                for (i = 1; i <= n; i++)
                {
                    if (data.active[i - 1] != 0)
                    {
                        nactiv = i;
                    }
                }

                if (nactiv <= 0)
                {
                    done = true;
                }
                else
                {
                    j = data.invers[nactiv - 1];

                    p[j - 1] = p[j + data.idir[nactiv - 1] - 1];
                    p[j + data.idir[nactiv - 1] - 1] = nactiv;

                    data.invers[nactiv - 1] = data.invers[nactiv - 1] + data.idir[nactiv - 1];
                    data.invers[p[j - 1] - 1] = j;

                    if (j + 2 * data.idir[nactiv - 1] < 1 || n < j + 2 * data.idir[nactiv - 1])
                    {
                        data.idir[nactiv - 1] = -data.idir[nactiv - 1];
                        data.active[nactiv - 1] = 0;
                    }
                    else if (nactiv < p[j + 2 * data.idir[nactiv - 1] - 1])
                    {
                        data.idir[nactiv - 1] = -data.idir[nactiv - 1];
                        data.active[nactiv - 1] = 0;
                    }

                    for (i = nactiv; i < n; i++)
                    {
                        data.active[i] = 1;
                    }
                }
            }

            if (done)
            {
                data.active = null;
                data.idir = null;
                data.invers = null;
            }
        }

        public static void perm0_next3(int n, int[] p, ref bool more, ref int rank )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PERM0_NEXT3 computes permutations of (0,...,N-1).
        //
        //  Discussion:
        //
        //    The routine is initialized by calling with MORE = TRUE in which case
        //    it returns the identity permutation.
        //
        //    If the routine is called with MORE = FALSE then the successor of the
        //    input permutation is computed.
        //
        //    Trotter's algorithm is used.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 November 2018
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Hale Trotter,
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Hale Trotter,
        //    PERM, Algorithm 115,
        //    Communications of the Association for Computing Machinery,
        //    Volume 5, 1962, pages 434-435.
        //
        //  Parameters:
        //
        //    Input, int N, the number of objects being permuted.
        //
        //    Input/output, int P[N], the permutation, in standard index form.
        //    If MORE is TRUE, then P is assumed to contain the
        //    "previous" permutation, and on P(I) is the value
        //    of the I-th object under the next permutation.
        //    Otherwise, P will be set to the "first" permutation.
        //
        //    Input/output, bool &MORE.
        //    Set MORE = FALSE before first calling this routine.
        //    MORE will be reset to TRUE and a permutation will be returned.
        //    Each new call produces a new permutation until MORE is returned FALSE
        //
        //    Input/output, int *RANK, the rank of the current permutation.
        //
        {
            int i;
            int m2;
            int n2;
            int q;
            int s;
            int t;
            int temp;

            if (!more)
            {
                for (i = 0; i < n; i++)
                {
                    p[i] = i;
                }

                more = true;
                rank = 1;
            }
            else
            {
                n2 = n;
                m2 = rank;
                s = n;

                for (;;)
                {
                    q = m2 % n2;
                    t = m2 % (2 * n2);

                    if (q != 0)
                    {
                        break;
                    }

                    if (t == 0)
                    {
                        s = s - 1;
                    }

                    m2 = m2 / n2;
                    n2 = n2 - 1;
                    if (n2 == 0)
                    {
                        for (i = 0; i < n; i++)
                        {
                            p[i] = i;
                        }

                        more = false;
                        rank = 1;
                        break;
                    }
                }

                if (n2 != 0)
                {
                    if (q == t)
                    {
                        s = s - q;
                    }
                    else
                    {
                        s = s + q - n2;
                    }

                    temp = p[s - 1];
                    p[s - 1] = p[s];
                    p[s] = temp;

                    rank = rank + 1;
                }
            }
        }

        public static void perm0_print(int n, int[] p, string title )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PERM0_PRINT prints a permutation of (0,...,N-1).
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
        //    23 June 2015
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
            string cout = "";

            if (0 < title.Length)
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
                    cout = "  ";
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

                    cout = "  ";
                    for (i = ilo; i < ihi; i++)
                    {
                        cout += p[i].ToString().PadLeft(4);
                    }

                    Console.WriteLine(cout);
                }
            }

            return;
        }

        public static void perm0_random(int n, ref int seed, ref int[] p )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PERM0_RANDOM selects a random permutation of (0,...,N-1).
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
        //    29 May 2003
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
        //    Input, int N, the number of objects to be permuted.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, int P[N], a permutation of (0, 1, 2, ..., N-1).
        //
        {
            int i;
            int j;
            int k;

            for (i = 0; i < n; i++)
            {
                p[i] = i;
            }

            for (i = 0; i < n - 1; i++)
            {
                j = UniformRNG.i4_uniform_ab(i, n - 1, ref seed);

                k = p[i];
                p[i] = p[j];
                p[j] = k;
            }
        }

        public static void perm0_random2(int n, ref int seed, ref int[] p )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PERM0_RANDOM2 selects a random permutation of (0,...,N-1).
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
        //    Original FORTRAN77 version by James Filliben.
        //    C++ version by John Burkardt
        //
        //  Reference:
        //
        //    K L Hoffman, D R Shier,
        //    Algorithm 564,
        //    A Test Problem Generator for Discrete Linear L1 Approximation Problems,
        //    ACM Transactions on Mathematical Software,
        //    Volume 6, Number 4, December 1980, pages 615-617.
        //
        //  Parameters:
        //
        //    Input, int N, the number of elements of the array.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, int P[N], a permutation, in standard index form.
        //
        {
            int i;
            int j;
            int temp;

            if (n < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("PERM0_RANDOM2 - Fatal error!");
                Console.WriteLine("  Illegal input value of N  = " + n + "");
                Console.WriteLine("  N must be at least 1!");
                return;
            }

            if (n == 1)
            {
                p[0] = 0;
                return;
            }

            typeMethods.i4vec_indicator0(n, ref p);

            for (i = 1; i <= n; i++)
            {
                j = i + UniformRNG.i4_uniform_ab(1, n, ref seed);

                if (n < j)
                {
                    j = j - n;
                }

                if (i != j)
                {
                    temp = p[j - 1];
                    p[j - 1] = p[i - 1];
                    p[i - 1] = temp;
                }
            }

            return;
        }

        public static int perm0_rank(int n, int[] p, ref int[] invers )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PERM0_RANK ranks a permutation of (0,...,N-1).
        //
        //  Discussion:
        //
        //    This is the same as asking for the step at which PERM_NEXT2
        //    would compute the permutation.  The value of the rank will be
        //    between 1 and N!.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 June 2015
        //
        //  Author:
        //
        //    John Burkardt.
        //
        //  Reference:
        //
        //    Dennis Stanton, Dennis White,
        //    Constructive Combinatorics,
        //    Springer, 1986,
        //    ISBN: 0387963472,
        //    LC: QA164.S79.
        //
        //  Parameters:
        //
        //    Input, int N, the number of elements in the set that
        //    is permuted by P.
        //
        //    Input, int P[N], a permutation, in standard index form.
        //
        //    Output, int INVERS[N], the inverse permutation of P.
        //    It is computed as part of the algorithm, and may be of use
        //    to the user.  INVERS(P(I)) = I for each entry I.
        //
        //    Output, int PERM0_RANK, the rank of the permutation.  This
        //    gives the order of the given permutation in the set of all
        //    the permutations on N elements.
        //
        {
            int count;
            int i;
            int j;
            int rank;
            int rem;
            //
            //  Make sure the permutation is a legal one.
            //  (This is not an efficient way to do so!)
            //
            if (!perm0_check(n, p))
            {
                Console.WriteLine("");
                Console.WriteLine("PERM0_RANK - Fatal error!");
                Console.WriteLine("  PERM0_CHECK rejects permutation.");
                return 1;
            }

            //
            //  Compute the inverse permutation.
            //
            for (i = 0; i < n; i++)
            {
                invers[i] = p[i];
            }

            perm0_inverse2(n, invers);

            rank = 0;

            for (i = 1; i <= n; i++)
            {

                count = 0;

                for (j = 0; j < invers[i - 1]; j++)
                {
                    if (p[j] < i)
                    {
                        count = count + 1;
                    }
                }

                if ((rank % 2) == 1)
                {
                    rem = count;
                }
                else
                {
                    rem = i - 1 - count;
                }

                rank = i * rank + rem;
            }

            rank = rank + 1;

            return rank;
        }

        public static int perm0_sign(int n, int[] p)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM0_SIGN returns the sign of a permutation of (0,...,N-1).
            //
            //  Discussion:
            //
            //    A permutation can always be replaced by a sequence of pairwise
            //    transpositions.  A given permutation can be represented by
            //    many different such transposition sequences, but the number of
            //    such transpositions will always be odd or always be even.
            //    If the number of transpositions is even or odd, the permutation is
            //    said to be even or odd.
            //
            //  Example:
            //
            //    Input:
            //
            //      N = 9
            //      P = 2, 3, 9, 6, 7, 8, 5, 4, 1
            //
            //    Output:
            //
            //      PERM0_SIGN = +1
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 November 2012
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
            //    Input, int N, the number of objects permuted.
            //
            //    Input, int P[N], a permutation, in standard index form.
            //
            //    Output, int PERM0_SIGN, the "sign" of the permutation.
            //    +1, the permutation is even,
            //    -1, the permutation is odd.
            //
        {
            int i;
            int j;
            int p_sign;
            int[] q;
            int temp;

            if (!perm0_check(n, p))
            {
                Console.WriteLine("");
                Console.WriteLine("PERM0_SIGN - Fatal error!");
                Console.WriteLine("  PERM0_CHECK rejects permutation.");
                return 1;
            }

            //
            //  Make a temporary copy of the permutation.
            //
            q = new int[n];
            for (i = 0; i < n; i++)
            {
                q[i] = p[i];
            }

            //
            //  Start with P_SIGN indicating an even permutation.
            //  Restore each element of the permutation to its correct position,
            //  updating P_SIGN as you go.
            //
            p_sign = 1;

            for (i = 1; i <= n - 1; i++)
            {
                j = typeMethods.i4vec_index(n, q, i);

                if (j != i - 1)
                {
                    temp = q[i - 1];
                    q[i - 1] = q[j];
                    q[j] = temp;

                    p_sign = -p_sign;
                }
            }
            
            return p_sign;
        }

        public static void perm0_to_equiv(int n, int[] p, ref int npart, ref int[] jarray, ref int[] iarray )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PERM0_TO_EQUIV finds the partition induced by a permutation of (0,...,N-1).
        //
        //  Example:
        //
        //    Input:
        //
        //      N = 9
        //      P = 1, 2, 8, 5, 6, 7, 4, 3, 0
        //
        //    Output:
        //
        //      NPART = 3
        //      JARRAY = 4, 3, 2
        //      IARRAY = 1, 1, 1, 2  3  2  3  2, 1
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
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of objects being permuted.
        //
        //    Input, int P[N], a permutation, in standard index form.
        //
        //    Output, int &NPART, number of subsets in the partition.
        //
        //    Output, int JARRAY[N].  JARRAY(I) is the number of elements
        //    in the I-th subset of the partition.
        //
        //    Output, int IARRAY[N].  IARRAY(I) is the class to which
        //    element I belongs.
        //
        {
            int i;
            int j;
            int k;

            if (!perm0_check(n, p))
            {
                Console.WriteLine("");
                Console.WriteLine("PERM0_TO_EQUIV - Fatal error!");
                Console.WriteLine("  PERM0_CHECK rejects permutation.");
                return;
            }

            //
            //  Initialize.
            //
            for (i = 0; i < n; i++)
            {
                iarray[i] = -1;
            }

            for (i = 0; i < n; i++)
            {
                jarray[i] = -1;
            }

            npart = 0;
            //
            //  Search for the next item J which has not been assigned a subset/orbit.
            //
            for (j = 1; j <= n; j++)
            {
                if (iarray[j - 1] != -1)
                {
                    continue;
                }

                //
                //  Begin a new subset/orbit.
                //
                npart = npart + 1;
                k = j;
                //
                //  Add the item to the subset/orbit.
                //
                for (;;)
                {
                    jarray[npart - 1] = jarray[npart - 1] + 1;
                    iarray[k - 1] = npart;
                    //
                    //  Apply the permutation.  If the permuted object isn't already in the
                    //  subset/orbit, add it.
                    //
                    k = p[k - 1];

                    if (iarray[k - 1] != -1)
                    {
                        break;
                    }
                }
            }

            return;
        }

        public static void perm0_to_inversion(int n, int[] p, int[] ins )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PERM0_TO_INVERSION: permutation (0,...,N-1) to inversion sequence.
        //
        //  Discussion:
        //
        //    For a given permutation P acting on objects 1 through N, the inversion
        //    sequence INS is defined as:
        //
        //      INS(1) = 0
        //      INS(I) = number of values J < I for which P(I) < P(J).
        //
        //  Example:
        //
        //    Input:
        //
        //      ( 2, 4, 0, 3, 1 )
        //
        //    Output:
        //
        //      ( 0, 0, 2, 1, 3 )
        //
        //    The original permutation can be recovered from the inversion sequence.
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
        //    John Burkardt
        //
        //  Reference:
        //
        //    Dennis Stanton, Dennis White,
        //    Constructive Combinatorics,
        //    Springer, 1986,
        //    ISBN: 0387963472,
        //    LC: QA164.S79.
        //
        //  Parameters:
        //
        //    Input, int N, the number of objects being permuted.
        //
        //    Input, int P[N], the permutation, in standard index form.
        //    The I-th item has been mapped to P(I).
        //
        //    Output, int INS[N], the inversion sequence of the permutation.
        //
        {
            int i;
            int j;

            if (!perm0_check(n, p))
            {
                Console.WriteLine("");
                Console.WriteLine("PERM0_TO_INVERSION - Fatal error!");
                Console.WriteLine("  PERM0_CHECK rejects permutation.");
                return;
            }

            for (i = 0; i < n; i++)
            {
                ins[i] = 0;
            }

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < i; j++)
                {
                    if (p[i] < p[j])
                    {
                        ins[i] = ins[i] + 1;
                    }
                }
            }
        }

        public static void perm0_to_ytb(int n, int[] p, ref int[] lambda, ref int[] a )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PERM0_TO_YTB converts a permutation of (0,...,N-1) to a Young tableau.
        //
        //  Discussion:
        //
        //    The mapping is not invertible.  In most cases, several permutations
        //    correspond to the same tableau.
        //
        //  Example:
        //
        //    N = 7
        //    P = 6 1 3 0 4 2 5
        //
        //    YTB =
        //      1 2 3 6
        //      4 5
        //      7
        //
        //    LAMBDA = 4 2 1 0 0 0 0
        //
        //    A = 1 1 1 2 2 1 3
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
        //    Input, int N, the integer to be partitioned.
        //
        //    Input, int P[N], a permutation, in standard index form.
        //
        //    Output, int LAMBDA[N].  LAMBDA[I] is the length of the I-th row.
        //
        //    Output, int A[N].  A[I] is the row containing I.
        //
        {
            bool another;
            int compare;
            int i;
            int put_index;
            int put_row;
            int put_value;
            //
            //  Initialize.
            //
            for (i = 0; i < n; i++)
            {
                lambda[i] = 0;
            }

            for (i = 0; i < n; i++)
            {
                a[i] = 0;
            }

            //
            //  Now insert each item of the permutation.
            //
            for (put_index = 1; put_index <= n; put_index++)
            {
                put_value = p[put_index - 1];
                put_row = 1;

                for (;;)
                {
                    another = false;

                    for (compare = put_value + 1; compare <= n; compare++)
                    {
                        if (a[compare - 1] == put_row)
                        {
                            another = true;
                            a[put_value - 1] = put_row;
                            a[compare - 1] = 0;
                            put_value = compare;
                            put_row = put_row + 1;
                            break;
                        }
                    }

                    if (!another)
                    {
                        break;
                    }
                }

                a[put_value - 1] = put_row;
                lambda[put_row - 1] = lambda[put_row - 1] + 1;

            }

            return;
        }

        public static void perm0_unrank(int n, int rank, int[] p)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM0_UNRANK "unranks" a permutation of (0,...,N-1).
            //
            //  Discussion:
            //
            //    That is, given a rank, it computes the corresponding permutation.
            //    This is the same as asking for the permutation which PERM0_NEXT2
            //    would compute at the RANK-th step.
            //
            //    The value of the rank should be between 1 and N!.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    12 June 2015
            //
            //  Author:
            //
            //    John Burkardt.
            //
            //  Reference:
            //
            //    Dennis Stanton, Dennis White,
            //    Constructive Combinatorics,
            //    Springer, 1986,
            //    ISBN: 0387963472,
            //    LC: QA164.S79.
            //
            //  Parameters:
            //
            //    Input, int N, the number of elements in the set.
            //
            //    Input, int RANK, the desired rank of the permutation.  This
            //    gives the order of the given permutation in the set of all
            //    the permutations on N elements, using the ordering of PERM0_NEXT2.
            //
            //    Output, int P[N], the permutation, in standard index form.
            //
        {
            int i;
            int icount;
            int iprev;
            int irem;
            int j;
            int jdir;
            int jrank;

            for (i = 0; i < n; i++)
            {
                p[i] = -1;
            }

            if (rank < 1 || typeMethods.i4_factorial(n) < rank)
            {
                Console.WriteLine("");
                Console.WriteLine("PERM0_UNRANK - Fatal error!");
                Console.WriteLine("  Illegal input value for RANK.");
                Console.WriteLine("  RANK must be between 1 and N!,");
                Console.WriteLine("  but the input value is " + rank + "");
                return;
            }

            jrank = rank - 1;

            for (i = 1; i <= n; i++)
            {
                iprev = n + 1 - i;
                irem = jrank % iprev;
                jrank = jrank / iprev;

                if ((jrank % 2) == 1)
                {
                    j = 0;
                    jdir = 1;
                }
                else
                {
                    j = n + 1;
                    jdir = -1;
                }

                icount = 0;

                for (;;)
                {
                    j = j + jdir;

                    if (p[j - 1] == -1)
                    {
                        icount = icount + 1;
                    }

                    if (irem < icount)
                    {
                        break;
                    }
                }

                p[j - 1] = iprev - 1;
            }
        }

        public static void perm1_canon_to_cycle(int n, int[] p1, ref int[] p2 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PERM1_CANON_TO_CYCLE: permutation of (1,...,N) from canonical to cycle form.
        //
        //  Discussion:
        //
        //    This algorithm relies on the assumption that the entres in the permutation
        //    vector are strictly positive, and uses the sign of the entries as a temporary
        //    marker.  Thus it cannot be implemented directly with a 0-based permutation.
        //
        //  Example:
        //
        //    Input:
        //
        //      4 5 2 1 6 3
        //
        //    Output:
        //
        //      -4 5 -2 -1 6 3,
        //      indicating the cycle structure
        //      ( 4, 5 ) ( 2 ) ( 1, 6, 3 )
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
        //    John Burkardt
        //
        //  Reference:
        //
        //    Donald Knuth,
        //    The Art of Computer Programming,
        //    Volume 1, Fundamental Algorithms,
        //    Addison Wesley, 1968, page 176.
        //
        //  Parameters:
        //
        //    Input, int N, the number of objects permuted.
        //
        //    Input, int P1[N], the permutation, in canonical form.
        //
        //    Output, int P2[N], the permutation, in cycle form.
        //
        {
            int i;
            int pmin;

            for (i = 0; i < n; i++)
            {
                p2[i] = p1[i];
            }

            pmin = p2[0] + 1;

            for (i = 0; i < n; i++)
            {
                if (p2[i] < pmin)
                {
                    pmin = p2[i];
                    p2[i] = -p2[i];
                }
            }

            return;
        }

        public static bool perm1_check(int n, int[] p)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PERM1_CHECK checks a permutation of (1, ..., N ).
            //
            //  Discussion:
            //
            //    The routine verifies that each of the integers from 0 to
            //    to N-1 occurs among the N entries of the permutation.
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
            //    Input, int N, the number of entries.
            //
            //    Input, int P[N], the array to check.
            //
            //    Output, bool PERM1_CHECK, is 
            //    TRUE if P is a legal permutation of 1,...,N.
            //    FALSE if P is not a legal permuation of 1,...,N.
            //
        {
            bool check;
            int location;
            int value;

            check = true;

            for (value = 1; value <= n; value++)
            {
                check = false;

                for (location = 0; location < n; location++)
                {
                    if (p[location] == value)
                    {
                        check = true;
                        break;
                    }
                }

                if (!check)
                {
                    Console.WriteLine("");
                    Console.WriteLine("PERM1_CHECK - Warning!");
                    Console.WriteLine("  Permutation is missing value " + value + "");
                    break;
                }

            }

            return check;
        }

        public static void perm1_cycle_to_canon(int n, int[] p1, ref int[] p2 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PERM1_CYCLE_TO_CANON: permutation of (1,...,N) from cycle to canonical form.
        //
        //  Discussion:
        //
        //    This algorithm relies on the assumption that the entres in the permutation
        //    vector are strictly positive, and uses the sign of the entries as a temporary
        //    marker.  Thus it cannot be implemented directly with a 0-based permutation.
        //
        //  Example:
        //
        //    Input:
        //
        //      -6 3 1 -5, 4 -2,
        //      indicating the cycle structure
        //      ( 6, 3, 1 ) ( 5, 4 ) ( 2 )
        //
        //    Output:
        //
        //      4 5 2 1 6 3
        //
        //  Discussion:
        //
        //    The procedure is to "rotate" the elements of each cycle so that
        //    the smallest element is first:
        //
        //      ( 1, 6, 3 ) ( 4, 5 ) ( 2 )
        //
        //    and then to sort the cycles in decreasing order of their first
        //    (and lowest) element:
        //
        //      ( 4, 5 ) ( 2 ) ( 1, 6, 3 )
        //
        //    and then to drop the parentheses:
        //
        //      4 5 2 1 6 3
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    24 June 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Donald Knuth,
        //    The Art of Computer Programming,
        //    Volume 1, Fundamental Algorithms,
        //    Addison Wesley, 1968, pages 176.
        //
        //  Parameters:
        //
        //    Input, int N, the number of objects permuted.
        //
        //    Input, int P1[N], the permutation, in cycle form.
        //
        //    Output, int P2[N], the permutation, in canonical form.
        //
        {
            int[] hi;
            int i;
            int[] indx;
            int j;
            int k;
            int[] lo;
            int ncycle;
            int next;
            int nhi;
            int nlo;
            int nmin;
            int[] pmin;
            int[] ptemp;

            hi = new int[n];
            lo = new int[n];
            pmin = new int[n];
            ptemp = new int[n];

            for (i = 0; i < n; i++)
            {
                p2[i] = p1[i];
            }

            //
            //  Work on the next cycle.
            //
            nlo = 1;
            ncycle = 0;

            while (nlo <= n)
                //
                //  Identify NHI, the last index in this cycle.
                //
            {
                ncycle = ncycle + 1;

                nhi = nlo;

                while (nhi < n)
                {
                    if (p2[nhi] < 0)
                    {
                        break;
                    }

                    nhi = nhi + 1;
                }

                //
                //  Identify the smallest value in this cycle.
                //
                p2[nlo - 1] = -p2[nlo - 1];
                pmin[ncycle - 1] = p2[nlo - 1];
                nmin = nlo;

                for (i = nlo + 1; i <= nhi; i++)
                {
                    if (p2[i - 1] < pmin[ncycle - 1])
                    {
                        pmin[ncycle - 1] = p2[i - 1];
                        nmin = i;
                    }
                }

                //
                //  Rotate the cycle so A_MIN occurs first.
                //
                for (i = nlo; i <= nmin - 1; i++)
                {
                    ptemp[i + nhi - nmin] = p2[i - 1];
                }

                for (i = nmin; i <= nhi; i++)
                {
                    ptemp[i + nlo - nmin - 1] = p2[i - 1];
                }

                lo[ncycle - 1] = nlo;
                hi[ncycle - 1] = nhi;
                //
                //  Prepare to operate on the next cycle.
                //
                nlo = nhi + 1;
            }

            //
            //  Compute a sorting index for the cycle minima.
            //
            indx = typeMethods.i4vec_sort_heap_index_d(ncycle, pmin);
            //
            //  Copy the cycles out of the temporary array in sorted order.
            //
            j = 0;

            for (i = 0; i < ncycle; i++)
            {
                next = indx[i];
                nlo = lo[next];
                nhi = hi[next];

                for (k = nlo; k <= nhi; k++)
                {
                    j = j + 1;
                    p2[j - 1] = ptemp[k - 1];
                }
            }
        }

        public static void perm1_cycle_to_index(int n, int[] p1, int[] p2 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PERM1_CYCLE_TO_INDEX: permutation of (1,...,N) from cycle to index form.
        //
        //  Discussion:
        //
        //    This algorithm relies on the assumption that the entres in the permutation
        //    vector are strictly positive, and uses the sign of the entries as a temporary
        //    marker.  Thus it cannot be implemented directly with a 0-based permutation.
        //
        //  Example:
        //
        //    Input:
        //
        //      N = 9
        //      P1 = -1, 2, 3, 9, -4, 6, 8, -5, 7
        //
        //    Output:
        //
        //      P2 = 2, 3, 9, 6, 7, 8, 5, 4, 1
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
        //    Input, int N, the number of objects being permuted.
        //
        //    Input, int P1[N], the permutation, in cycle form.
        //
        //    Output, int P2[N], the permutation, in standard index form.
        //
        {
            int j;
            int k1 = 0;
            int k2 = 0;
            int k3 = 0;

            for (j = 1; j <= n; j++)
            {
                k1 = p1[j - 1];

                if (k1 < 0)
                {
                    k1 = -k1;
                    k3 = k1;
                }

                if (j + 1 <= n)
                {
                    k2 = p1[j];
                    if (k2 < 0)
                    {
                        k2 = k3;
                    }
                }
                else
                {
                    k2 = k3;
                }

                p2[k1 - 1] = k2;
            }
        }

        public static void perm1_index_to_cycle(int n, int[] p1, int[] p2 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PERM1_INDEX_TO_CYCLE: permutation of (1,...,N) from index to cycle form.
        //
        //  Discussion:
        //
        //    This algorithm relies on the assumption that the entres in the permutation
        //    vector are strictly positive, and uses the sign of the entries as a temporary
        //    marker.  Thus it cannot be implemented directly with a 0-based permutation.
        //
        //  Example:
        //
        //    Input:
        //
        //      N = 9
        //      P1 = 2, 3, 9, 6, 7, 8, 5, 4, 1
        //
        //    Output:
        //
        //      P2 = -1, 2, 3, 9, -4, 6, 8, -5, 7
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    24 June 2015
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
        //    Input, int N, the number of objects being permuted.
        //
        //    Input, int P1[N], the permutation, in standard index form.
        //
        //    Output, int P2[N], the permutation, in cycle form.
        //
        {
            int i;
            int j;
            int k;

            i = 0;
            j = 1;

            while (j <= n)
            {
                if (p1[j - 1] < 0)
                {
                    j = j + 1;
                }
                else
                {
                    k = j;

                    i = i + 1;
                    p2[i - 1] = -k;

                    while (p1[k - 1] != j)
                    {
                        i = i + 1;
                        p2[i - 1] = p1[k - 1];
                        p1[k - 1] = -p1[k - 1];
                        k = Math.Abs(p1[k - 1]);
                    }

                    p1[k - 1] = -p1[k - 1];
                }
            }

            for (i = 0; i < n; i++)
            {
                p1[i] = Math.Abs(p1[i]);
            }
        }

        public static void perm1_print(int n, int[] p, string title )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PERM1_PRINT prints a permutation of (1,...,N).
        //
        //  Example:
        //
        //    Input:
        //
        //      P = 7 2 3 1 5 3 6
        //
        //    Printed output:
        //
        //      "This is the permutation:"
        //
        //      1 2 3 4 5 6 7
        //      7 2 3 1 5 3 6
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 June 2015
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
            string cout = "";

            if (0 < title.Length)
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
                    cout = "  ";
                    for (i = ilo; i < ihi; i++)
                    {
                        cout += (i + 1).ToString().PadLeft(4);
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

                    cout = "  ";
                    for (i = ilo; i < ihi; i++)
                    {
                        cout += p[i].ToString().PadLeft(4);
                    }

                    Console.WriteLine(cout);
                }
            }
        }

        public static void random_permutation(int n, ref double[] x, ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    RANDOM_PERMUTATION applies a random permutation to an array.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    19 February 2016
            //
            //  Author:
            //
            //    Original C version by Warren Smith.
            //    This C++ version by John Burkardt.
            //
            //  Parameters:
            //
            //    Input, unsigned int N, indicates the size of X.
            //
            //    Input/output, double X[N+2].  On output, entries X[1] through
            //    X[N] have been randomly permuted.
            //
            //    Input/output, int &SEED, a seed for the random number generator.
            //
        {
            int i;
            int j;
            double t;

            for (i = 1; i < n; i++)
            {
                j = UniformRNG.i4_uniform_ab(i, n, ref seed);

                t = x[i];
                x[i] = x[j];
                x[j] = t;
            }
        }

    }
}