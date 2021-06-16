using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.RankingNS
{
    public static partial class Ranking
    {
        public static int npart_enum(int n, int npart)

            //****************************************************************************80
            // 
            //  Purpose:
            //
            //    NPART_ENUM enumerates the number of partitions of N with NPART parts.
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
            //    Input, int N, the integer to be partitioned.
            //    Normally N must be positive, but for this routine any
            //    N is allowed.
            // 
            //    Input, int NPART, the number of parts of the partition.
            //    Normally, 1 <= NPART <= N is required,
            //    but for this routine any value of NPART is allowed.
            // 
            //    Output, int NPART_ENUM is the number of partitions of N
            //    with NPART parts.
            // 
        {
            int[] p;
            int value;

            if (n <= 0)
            {
                value = 0;
            }
            else if (npart <= 0 || n < npart)
            {
                value = 0;
            }
            else
            {
                p = npart_table(n, npart);

                value = p[n + npart * (n + 1)];
            }

            return value;
        }

        public static int[] npart_rsf_lex_random(int n, int npart, ref int seed)

            //****************************************************************************80
            // 
            //  Purpose:
            //
            //    NPART_RSF_LEX_RANDOM returns a random RSF NPART partition.
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
            //    Input, int N, the integer to be partitioned.
            //    N must be positive.
            // 
            //    Input, int NPART, the number of parts of the partition.
            //    1 <= NPART <= N.
            // 
            //    Input/output, int &SEED, a seed for the random number
            //    generator.
            // 
            //    Output, int NPART_RSF_LEX_RANDOM[NPART], contains the partition.
            //    A(1) through A(NPART) contain the nonzero integers which
            //    sum to N.
            // 
        {
            int[] a;
            int npartitions;
            int rank;

            npartitions = npart_enum(n, npart);

            rank = UniformRNG.i4_uniform_ab(1, npartitions, ref seed);

            a = npart_rsf_lex_unrank(rank, n, npart);

            return a;
        }

        public static int npart_rsf_lex_rank(int n, int npart, int[] a)

            //****************************************************************************80
            // 
            //  Purpose:
            //
            //    NPART_RSF_LEX_RANK computes the lex rank of an RSF NPART partition.
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
            //    Input, int N, the integer to be partitioned.
            //    N must be positive.
            // 
            //    Input, int NPART, the number of parts of the partition.
            //    1 <= NPART <= N.
            // 
            //    Input, int A[NPART], contains the partition.
            //    A(1) through A(NPART) contain the nonzero integers which
            //    sum to N.
            // 
            //    Output, int NPART_RSF_LEX_RANK, the rank of the partition.
            // 
        {
            int[] b;
            bool check;
            int i;
            int ncopy;
            int npartcopy;
            int[] p;
            int rank;
            // 
            //  Check.
            // 
            check = part_rsf_check(n, npart, a);

            if (!check)
            {
                Console.WriteLine("");
                Console.WriteLine("NPART_RSF_LEX_RANK - Fatal error!");
                Console.WriteLine("  The input array is illegal.");
                return (1);
            }

            // 
            //  Get the table of partitions of N with NPART parts.
            // 
            p = npart_table(n, npart);
            // 
            //  Copy the partition "backwards".
            //
            b = new int[npart];

            for (i = 1; i <= npart; i++)
            {
                b[i - 1] = a[npart - i];
            }

            rank = 0;
            ncopy = n;
            npartcopy = npart;

            while (0 < ncopy && 0 < npartcopy)
            {
                if (b[npartcopy - 1] == 1)
                {
                    ncopy = ncopy - 1;
                    npartcopy = npartcopy - 1;
                }
                else
                {
                    for (i = 0; i < npartcopy; i++)
                    {
                        b[i] = b[i] - 1;
                    }

                    rank = rank + p[ncopy - 1 + (npartcopy - 1) * (n + 1)];
                    ncopy = ncopy - npartcopy;
                }
            }

            return rank;
        }

        public static void npart_rsf_lex_successor(int n, int npart, ref int[] a, ref int rank )

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    NPART_RSF_LEX_SUCCESSOR computes the RSF lex successor for NPART partitions.
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
        //    Input, int N, the integer to be partitioned.
        //    N must be at least 1.
        // 
        //    Input, int NPART, the number of parts of the partition.
        //    1 <= NPART <= N.
        // 
        //    Input/output, int A[NPART], contains the partition.
        //    A(1) through A(NPART) contain the nonzero integers which
        //    sum to N.
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
            int i;
            int j;
            // 
            //  Return the first element.
            // 
            if (rank == -1)
            {
                if (npart < 1)
                {
                    Console.WriteLine("");
                    Console.WriteLine("NPART_RSF_LEX_SUCCESSOR - Fatal error!");
                    Console.WriteLine("  NPART < 1.");
                    return;
                }

                for (i = 0; i < npart - 1; i++)
                {
                    a[i] = 1;
                }

                a[npart - 1] = n - (npart - 1);

                rank = 0;
                return;
            }

            // 
            //  Check.
            // 
            check = part_rsf_check(n, npart, a);

            if (!check)
            {
                Console.WriteLine("");
                Console.WriteLine("NPART_RSF_LEX_SUCCESSOR - Fatal error!");
                Console.WriteLine("  The input array is illegal.");
                return;
            }

            // 
            //  Find the first index I for which A(NPART+1-I) + 1 < A(NPART).
            // 
            i = 2;

            for (;;)
            {
                if (npart < i)
                {
                    break;
                }

                if (a[npart - i] + 1 < a[npart - 1])
                {
                    break;
                }

                i = i + 1;
            }

            // 
            //  If no such index, we''ve reached the end of the line.
            // 
            if (i == npart + 1)
            {
                for (i = 0; i < npart - 1; i++)
                {
                    a[i] = 1;
                }

                a[npart - 1] = n - (npart - 1);

                rank = 0;
                return;
            }
            // 
            //  Otherwise, increment A(NPART+1-I), and adjust other entries.
            // 
            else
            {
                a[npart - i] = a[npart - i] + 1;
                d = -1;

                for (j = i - 1; 2 <= j; j--)
                {
                    d = d + a[npart - j] - a[npart - i];
                    a[npart - j] = a[npart - i];
                }

                a[npart - 1] = a[npart - 1] + d;
            }

            rank = rank + 1;

            return;
        }

        public static int[] npart_rsf_lex_unrank(int rank, int n, int npart)

            //****************************************************************************80
            // 
            //  Purpose:
            //
            //    NPART_RSF_LEX_UNRANK unranks an RSF NPART partition in the lex ordering.
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
            //    Input, int RANK, the rank of the partition.
            // 
            //    Input, int N, the integer to be partitioned.
            //    N must be positive.
            // 
            //    Input, int NPART, the number of parts of the partition.
            //    1 <= NPART <= N.
            // 
            //    Output, int NPART_RSF_LEX_UNRANK[NPART], contains the partition.
            //    A(1) through A(NPART) contain the nonzero integers which
            //    sum to N.
            // 
        {
            int[] a;
            int i;
            int ncopy;
            int npartcopy;
            int npartitions;
            int[] p;
            int rank_copy;
            // 
            //  Check.
            // 
            if (n <= 0)
            {
                Console.WriteLine("");
                Console.WriteLine("NPART_RSF_LEX_UNRANK - Fatal error!");
                Console.WriteLine("  The input N is illegal.");
                return null;
            }

            if (npart < 1 || n < npart)
            {
                Console.WriteLine("");
                Console.WriteLine("NPART_RSF_LEX_UNRANK - Fatal error!");
                Console.WriteLine("  The input NPART is illegal.");
                return null;
            }

            npartitions = npart_enum(n, npart);

            if (rank < 0 || npartitions < rank)
            {
                Console.WriteLine("");
                Console.WriteLine("NPART_RSF_LEX_UNRANK - Fatal error!");
                Console.WriteLine("  The input rank is illegal.");
                return null;
            }

            // 
            //  Get the table of partitions of N with NPART parts.
            // 
            p = npart_table(n, npart);

            a = new int[npart];

            for (i = 0; i < npart; i++)
            {
                a[i] = 0;
            }

            rank_copy = rank;
            ncopy = n;
            npartcopy = npart;

            while (0 < ncopy)
            {
                if (rank_copy < p[ncopy - 1 + (npartcopy - 1) * (n + 1)])
                {
                    a[npart - npartcopy] = a[npart - npartcopy] + 1;
                    ncopy = ncopy - 1;
                    npartcopy = npartcopy - 1;
                }
                else
                {
                    for (i = 1; i <= npartcopy; i++)
                    {
                        a[npart - i] = a[npart - i] + 1;
                    }

                    rank_copy = rank_copy - p[ncopy - 1 + (npartcopy - 1) * (n + 1)];
                    ncopy = ncopy - npartcopy;
                }
            }

            return a;
        }

        public static void npart_sf_lex_successor(int n, int npart, ref int[] a, ref int rank )

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    NPART_SF_LEX_SUCCESSOR computes SF NPART partition.
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
        //    Input, int N, the integer to be partitioned.
        //    N must be positive.
        // 
        //    Input, int NPART, the number of parts of the partition.
        //    1 <= NPART <= N.
        // 
        //    Input/output, int A[NPART], contains the partition.
        //    A(1) through A(NPART) contain the nonzero integers which
        //    sum to N.  The values in A must be in DESCENDING order.
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
            int indx = 0;
            int temp;
            // 
            //  Return the first element.
            // 
            if (rank == -1)
            {
                typeMethods.i4vec_part2(n, npart, ref a);
                rank = 0;
                return;
            }

            // 
            //  Check.
            // 
            check = part_sf_check(n, npart, a);

            if (!check)
            {
                Console.WriteLine("");
                Console.WriteLine("NPART_SF_LEX_SUCCESSOR - Fatal error!");
                Console.WriteLine("  The partition is illegal.");
                return;
            }

            // 
            //  Find the last entry that is 2 or more.
            // 
            for (i = npart; 1 <= i; i--)
            {
                if (1 < a[i - 1])
                {
                    indx = i;
                    break;
                }
            }

            // 
            //  As long as the last nonunit occurs after the first position,
            //  have it donate 1 to the left.
            // 
            if (1 < indx)
            {
                a[indx - 1] = a[indx - 1] - 1;
                a[indx - 2] = a[indx - 2] + 1;
                indx = indx - 1;

                for (;;)
                {
                    if (indx <= 1)
                    {
                        break;
                    }

                    if (a[indx - 1] <= a[indx - 2])
                    {
                        break;
                    }

                    temp = a[indx - 1];
                    a[indx - 1] = a[indx - 2];
                    a[indx - 2] = temp;

                    indx = indx - 1;
                }

                // 
                //  Sum the tail.
                // 
                temp = 0;
                for (i = indx; i < npart; i++)
                {
                    temp = temp + a[i];
                }

                // 
                //  Partition the tail sum equally over the tail.
                // 
                typeMethods.i4vec_part2(temp, npart - indx, ref a, xIndex: + indx);

                rank = rank + 1;
            }
            // 
            //  If A(2) through A(NPART) are 1, then this is the last element.
            //  Return the first one.
            // 
            else
            {
                typeMethods.i4vec_part2(n, npart, ref a);
                rank = 0;
            }
        }

        public static int[] npart_table(int n, int npart)

            //****************************************************************************80
            // 
            //  Purpose:
            //
            //    NPART_TABLE tabulates the number of partitions of N having NPART parts.
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
            //    Input, int N, the integer to be partitioned.
            //    N must be positive.
            // 
            //    Input, int NPART, the number of parts of the partition.
            //    1 <= NPART <= N.
            // 
            //    Output, int NPART_TABLE[(N+1)*(NPART+1)], P(I,J) is the number of
            //    partitions of I having J parts.
            // 
        {
            int i;
            int j;
            int[] p;

            p = new int[(n + 1) * (npart + 1)];

            p[0 + 0 * (n + 1)] = 1;
            for (i = 1; i <= n; i++)
            {
                p[i + 0 * (n + 1)] = 0;
            }

            for (i = 1; i <= n; i++)
            {
                for (j = 1; j <= npart; j++)
                {
                    if (i < j)
                    {
                        p[i + j * (n + 1)] = 0;
                    }
                    else if (i < 2 * j)
                    {
                        p[i + j * (n + 1)] = p[i - 1 + (j - 1) * (n + 1)];
                    }
                    else
                    {
                        p[i + j * (n + 1)] = p[i - 1 + (j - 1) * (n + 1)] + p[i - j + j * (n + 1)];
                    }
                }
            }

            return p;
        }
    }
}