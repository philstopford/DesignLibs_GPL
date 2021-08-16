using System;
using Burkardt.Function;
using Burkardt.Uniform;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static void i4_partition_conj(int n, int[] iarray1, int[] mult1, int npart1,
                ref int[] iarray2, ref int[] mult2, ref int npart2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4_PARTITION_CONJ computes the conjugate of a partition.
            //
            //  Discussion:
            //
            //    A partition of an integer N is a set of positive integers which
            //    add up to N.  The conjugate of a partition P1 of N is another partition
            //    P2 of N obtained in the following way:
            //
            //      The first element of P2 is the number of parts of P1 greater than
            //      or equal to 1.
            //
            //      The K-th element of P2 is the number of parts of P1 greater than
            //      or equal to K.
            //
            //    Clearly, P2 will have no more than N elements; it may be surprising
            //    to find that P2 is guaranteed to be a partition of N.  However, if
            //    we symbolize the initial partition P1 by rows of X's, then we can
            //    see that P2 is simply produced by grouping by columns:
            //
            //        6 3 2 2 1
            //      5 X X X X X
            //      4 X X X X
            //      2 X X
            //      1 X
            //      1 X
            //      1 X
            //
            //  Example:
            //
            //    14 = 5 + 4 + 2 + 1 + 1 + 1
            //
            //    The conjugate partition is:
            //
            //    14 = 6 + 3 + 2 + 2 + 1
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    28 May 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters
            //
            //    Input, int N, the integer to be partitioned.
            //
            //    Input, int IARRAY1[NPART1], contains the parts of
            //    the partition.  The value of N is represented by
            //
            //      sum ( 1 <= I <= NPART1 ) MULT1(I) * IARRAY1(I).
            //
            //    Input, int MULT1[NPART1], counts the multiplicity of
            //    the parts of the partition.  MULT1(I) is the multiplicity
            //    of the part IARRAY1(I), for 1 <= I <= NPART1.
            //
            //    Input, int NPART1, the number of "parts" in the partition.
            //
            //    Output, int IARRAY2[N], contains the parts of
            //    the conjugate partition in entries 1 through NPART2.
            //
            //    Output, int MULT2[N], counts the multiplicity of
            //    the parts of the conjugate partition in entries 1 through NPART2.
            //
            //    Output, int &NPART2, the number of "parts" in the conjugate partition.
            //
        {
            int i;
            int itemp;
            int itest;

            for (i = 0; i < n; i++)
            {
                iarray2[i] = 0;
            }

            for (i = 0; i < n; i++)
            {
                mult2[i] = 0;
            }

            npart2 = 0;

            itest = 0;

            for (;;)
            {
                itest = itest + 1;

                itemp = 0;

                for (i = 0; i < npart1; i++)
                {
                    if (itest <= iarray1[i])
                    {
                        itemp = itemp + mult1[i];
                    }
                }

                if (itemp <= 0)
                {
                    break;
                }

                if (0 < npart2)
                {
                    if (itemp == iarray2[npart2 - 1])
                    {
                        mult2[npart2 - 1] = mult2[npart2 - 1] + 1;
                    }
                    else
                    {
                        npart2 = npart2 + 1;
                        iarray2[npart2 - 1] = itemp;
                        mult2[npart2 - 1] = 1;
                    }
                }
                else
                {
                    npart2 = npart2 + 1;
                    iarray2[npart2 - 1] = itemp;
                    mult2[npart2 - 1] = 1;
                }
            }
        }

        public static void i4_partition_count(int n, ref int[] p)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4_PARTITION_COUNT computes the number of partitions of an integer.
            //
            //  Discussion:
            //
            //    Partition numbers are difficult to compute.  This routine uses
            //    Euler's method, which observes that:
            //
            //      P(0) = 1
            //      P(N) =   P(N-1)  + P(N-2)
            //             - P(N-5)  - P(N-7)
            //             + P(N-12) + P(N-15)
            //             - ...
            //
            //      where the numbers 1, 2, 5, 7, ... to be subtracted from N in the
            //      indices are the successive pentagonal numbers, (with both positive 
            //      and negative indices) with the summation stopping when a negative 
            //      index is reached.
            //
            //  First values:
            //
            //    N   P
            //
            //    0   1
            //    1   1
            //    2   2
            //    3   3
            //    4   5
            //    5   7
            //    6  11
            //    7  15
            //    8  22
            //    9  30
            //   10  42
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    28 May 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    John Conway, Richard Guy,
            //    The Book of Numbers,
            //    Springer Verlag, 1996, page 95.
            //
            //  Parameters:
            //
            //    Input, int N, the index of the highest partition number desired.
            //
            //    Output, int P[N+1], the partition numbers.
            //
        {
            int i;
            int j;
            int pj;
            int sgn;

            p[0] = 1;

            for (i = 1; i <= n; i++)
            {
                p[i] = 0;

                j = 0;
                sgn = 1;

                for (;;)
                {
                    j = j + 1;
                    pj = Pentagon.pentagon_num(j);

                    if (i < pj)
                    {
                        break;
                    }

                    p[i] = p[i] + sgn * p[i - pj];
                    sgn = -sgn;
                }

                j = 0;
                sgn = 1;

                for (;;)
                {
                    j = j - 1;
                    pj = Pentagon.pentagon_num(j);

                    if (i < pj)
                    {
                        break;
                    }

                    p[i] = p[i] + sgn * p[i - pj];
                    sgn = -sgn;

                }
            }

        }

        public static int[] i4_partition_count2(int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4_PARTITION_COUNT2 computes the number of partitions of an integer.
            //
            //  First values:
            //
            //    N   P
            //
            //    0   1
            //    1   1
            //    2   2
            //    3   3
            //    4   5
            //    5   7
            //    6  11
            //    7  15
            //    8  22
            //    9  30
            //   10  42
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    01 August 2004
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
            //    Input, int N, the largest integer to be considered.
            //
            //    Output, int I4_PARTITION_COUNT2[0:N], the partition numbers.
            //
        {
            int i;
            int j;
            int[] p;
            int s;
            int t;
            int total;

            if (n < 0)
            {
                return null;
            }

            p = new int[n + 1];

            p[0] = 1;

            if (n < 1)
            {
                return p;
            }

            p[1] = 1;

            for (i = 2; i <= n; i++)
            {
                total = 0;

                for (t = 1; t <= i; t++)
                {
                    s = 0;
                    j = i;

                    for (;;)
                    {
                        j = j - t;

                        if (0 < j)
                        {
                            s = s + p[j];
                        }
                        else
                        {
                            if (j == 0)
                            {
                                s = s + 1;
                            }

                            break;
                        }
                    }

                    total = total + s * t;
                }

                p[i] = (int)(total / i);
            }

            return p;
        }

        public static void i4_partition_count_values(ref int n_data, ref int n, ref int c)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4_PARTITION_COUNT_VALUES returns some values of the int partition count.
            //
            //  Discussion:
            //
            //    A partition of an integer I is a representation of the integer
            //    as the sum of nonzero positive integers.  The order of the summands
            //    does not matter.  Thus, the number 5 has the following partitions
            //    and no more:
            //
            //    5 = 5
            //      = 4 + 1 
            //      = 3 + 2 
            //      = 3 + 1 + 1 
            //      = 2 + 2 + 1 
            //      = 2 + 1 + 1 + 1 
            //      = 1 + 1 + 1 + 1 + 1
            //
            //    so the number of partitions of 5 is 7.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    06 February 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Milton Abramowitz, Irene Stegun,
            //    Handbook of Mathematical Functions,
            //    US Department of Commerce, 1964,
            //    ISBN: 0-486-61272-4,
            //    LC: QA47.A34.
            //
            //  Parameters:
            //
            //    Input/output, int &N_DATA.
            //    On input, if N_DATA is 0, the first test data is returned, and N_DATA
            //    is set to 1.  On each subsequent call, the input value of N_DATA is
            //    incremented and that test data item is returned, if available.  When 
            //    there is no more test data, N_DATA is set to 0.
            //
            //    Output, int &N, the integer.
            //
            //    Output, int &C, the number of partitions of the integer.
            //
        {
            int N_MAX = 21;

            int[] c_vec =
            {
                1,
                1, 2, 3, 5, 7, 11, 15, 22, 30, 42,
                56, 77, 101, 135, 176, 231, 297, 385, 490, 627
            };
            int[] n_vec =
            {
                0,
                1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
                11, 12, 13, 14, 15, 16, 17, 18, 19, 20
            };

            if (n_data < 0)
            {
                n_data = 0;
            }

            if (N_MAX <= n_data)
            {
                n_data = 0;
                n = 0;
                c = 0;
            }
            else
            {
                n = n_vec[n_data];
                c = c_vec[n_data];
                n_data = n_data + 1;
            }
        }

        public static void i4_partition_next(ref bool done, int[] a, int[] mult, int n, ref int npart)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4_PARTITION_NEXT generates the partitions of an integer, one at a time.
            //
            //  Discussion:
            //
            //    The number of partitions of N is:
            //
            //      1     1
            //      2     2
            //      3     3
            //      4     5
            //      5     7
            //      6    11
            //      7    15
            //      8    22
            //      9    30
            //     10    42
            //     11    56
            //     12    77
            //     13   101
            //     14   135
            //     15   176
            //     16   231
            //     17   297
            //     18   385
            //     19   490
            //     20   627
            //     21   792
            //     22  1002
            //     23  1255
            //     24  1575
            //     25  1958
            //     26  2436
            //     27  3010
            //     28  3718
            //     29  4565
            //     30  5604
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    28 May 2003
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
            //    Input/output, bool &DONE.
            //    On first call, the user should set DONE to TRUE to signal
            //    that the program should initialize data.
            //    On each return, the programs sets DONE to FALSE if it
            //    has another partition to return.  If the program returns
            //    with DONE TRUE, then there are no more partitions.
            //
            //    Output, int A[N].  A contains the parts of
            //    the partition.  The value of N is represented by
            //      N = sum ( 1 <= I <= NPART ) MULT(I) * A(I).
            //
            //    Output, int MULT[N].  MULT counts the multiplicity of
            //    the parts of the partition.  MULT(I) is the multiplicity
            //    of the part A(I), for 1 <= I <= NPART.
            //
            //    Input, int N, the integer to be partitioned.
            //
            //    Output, int &NPART, the number of "parts" in the partition.
            //
        {
            int i;
            int is_;
            int iu;
            int iv;
            int iw;
            int k;
            int k1;

            if (n <= 0)
            {
                Console.WriteLine("");
                Console.WriteLine("I4_PARTITION_NEXT - Fatal error!");
                Console.WriteLine("  N must be positive.");
                Console.WriteLine("  The input value of N was " + n + "");
                return;
            }

            if (done)
            {
                a[0] = n;
                for (i = 1; i < n; i++)
                {
                    a[i] = 0;
                }

                mult[0] = 1;
                for (i = 1; i < n; i++)
                {
                    mult[i] = 0;
                }

                npart = 1;
                done = false;
            }
            else
            {
                if (1 < a[(npart) - 1] || 1 < npart)
                {
                    done = false;

                    if (a[(npart) - 1] == 1)
                    {
                        is_ = a[npart - 2] + mult[npart - 1];
                        k = npart - 1;
                    }
                    else
                    {
                        is_ = a[npart - 1];
                        k = npart;
                    }

                    iw = a[k - 1] - 1;
                    iu = is_ / iw;
                    iv = is_ % iw;
                    mult[k - 1] = mult[k - 1] - 1;

                    if (mult[k - 1] == 0)
                    {
                        k1 = k;
                    }
                    else
                    {
                        k1 = k + 1;
                    }

                    mult[k1 - 1] = iu;
                    a[k1 - 1] = iw;

                    if (iv == 0)
                    {
                        npart = k1;
                    }
                    else
                    {
                        mult[k1] = 1;
                        a[k1] = iv;
                        npart = k1 + 1;
                    }
                }
                else
                {
                    done = true;
                }
            }
        }

        public static void i4_partition_next2(int n, ref int[] a, ref int[] mult, ref int npart, ref bool more)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4_PARTITION_NEXT2 computes the partitions of an integer one at a time.
            //
            //  Discussion:
            //
            //    Unlike compositions, order is not important in a partition.  Thus
            //    the sequences 3+2+1 and 1+2+3 represent distinct compositions, but
            //    not distinct partitions.  Also 0 is never returned as one of the
            //    elements of the partition.
            //
            //  Example:
            //
            //    Sample partitions of 6 include:
            //
            //      6 = 4+1+1 = 3+2+1 = 2+2+2
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
            //    John Burkardt.
            //
            //  Parameters:
            //
            //    Input, int N, the integer whose partitions are desired.
            //
            //    Output, int A[N].  A(I) is the I-th distinct part
            //    of the partition, for I = 1, NPART.  Note that if a certain number
            //    shows up several times in the partition, it is listed only
            //    once in A, and its multiplicity is counted in MULT.
            //
            //    Output, int MULT[N].  MULT(I) is the multiplicity of A(I)
            //    in the partition, for I = 1, NPART; that is, the number of repeated
            //    times that A(I) is used in the partition.
            //
            //    Output, int &NPART, the number of distinct, nonzero parts in the
            //    output partition.
            //
            //    Input/output, bool &MORE.  Set MORE = FALSE on first call.  It
            //    will be reset TRUE on return with the first partition.
            //    Keep calling for more partitions until MORE
            //    is returned FALSE
            //
        {
            int iff;
            int is_;
            int isum;

            if (!more)
            {
                npart = 1;
                a[npart - 1] = n;
                mult[npart - 1] = 1;
                more = mult[npart - 1] != n;
                return;
            }

            isum = 1;

            if (a[npart - 1] <= 1)
            {
                isum = mult[npart - 1] + 1;
                npart = npart - 1;
            }

            iff = a[npart - 1] - 1;

            if (mult[npart - 1] != 1)
            {
                mult[npart - 1] = mult[npart - 1] - 1;
                npart = npart + 1;
            }

            a[npart - 1] = iff;
            mult[npart - 1] = 1 + (isum / iff);
            is_ = isum % iff;

            if (0 < is_)
            {
                npart = npart + 1;
                a[npart - 1] = is_;
                mult[npart - 1] = 1;
            }

            //
            //  There are more partitions, as long as we haven't just computed
            //  the last one, which is N copies of 1.
            //
            more = mult[npart - 1] != n;

        }

        public static void i4_partition_print(int n, int npart, int[] a, int[] mult)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4_PARTITION_PRINT prints a partition of an integer.
            //
            //  Discussion:
            //
            //    A partition of an int N is a representation of the integer as
            //    the sum of nonzero integers:
            //
            //      N = A1 + A2 + A3 + ...
            //
            //    It is standard practice to gather together all the values that 
            //    are equal, and replace them in the sum by a single term, multiplied
            //    by its "multiplicity":
            //
            //      N = M1 * A1 + M2 * A2 + ... + M(NPART) * A(NPART)
            //    
            //    In this representation, every A is a unique positive number, and 
            //    no M is zero (or negative).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    29 June 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the integer to be partitioned.
            //
            //    Input, int NPART, the number of "parts" in the partition.
            //
            //    Input, int A[NPART], the parts of the partition.  
            //
            //    Input, int MULT[NPART], the multiplicities of the parts.
            //
        {
            int i;

            string cout = "  " + n + " = ";

            for (i = 0; i < npart; i++)
            {
                if (0 < i)
                {
                    cout += " + ";
                }

                cout += mult[i] + " * " + a[i];
            }

            Console.WriteLine(cout);
        }

        public static void i4_partition_random(int n, int[] table, ref int seed, ref int[] a, ref int[] mult,
                ref int npart)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4_PARTITION_RANDOM selects a random partition of the int N.
            //
            //  Discussion:
            //
            //    Note that some elements of the partition may be 0.  The partition is
            //    returned as (MULT(I),I), with NPART nonzero entries in MULT, and
            //
            //      N = sum ( 1 <= I <= N ) MULT(I) * I.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    05 May 2003
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
            //    Input, int N, the integer to be partitioned.
            //
            //    Input, int TABLE[N], the number of partitions of each integer 
            //    from 1 to N.  This table may be computed by I4_PARTITION_COUNT2.
            //
            //    Input/output, int &SEED, a seed for the random number generator.
            //
            //    Output, int A[N], contains in A(1:NPART) the parts of the partition.
            //
            //    Output, int MULT[N], contains in MULT(1:NPART) the multiplicity
            //    of the parts.
            //
            //    Output, int &NPART, the number of parts in the partition chosen,
            //    that is, the number of integers I with nonzero multiplicity MULT(I).
            //
        {
            int i;
            int i1;
            int id;
            int j;
            int m;
            double z;

            m = n;
            npart = 0;
            for (i = 0; i < n; i++)
            {
                mult[i] = 0;
            }

            while (0 < m)
            {
                z = UniformRNG.r8_uniform_01(ref seed);
                z = m * table[m - 1] * z;
                id = 1;
                i1 = m;
                j = 0;

                for (;;)
                {
                    j = j + 1;
                    i1 = i1 - id;

                    if (i1 < 0)
                    {
                        id = id + 1;
                        i1 = m;
                        j = 0;
                        continue;
                    }

                    if (i1 == 0)
                    {
                        z = z - id;
                        if (0.0 < z)
                        {
                            id = id + 1;
                            i1 = m;
                            j = 0;
                            continue;
                        }
                        else
                        {
                            break;
                        }
                    }

                    if (0 < i1)
                    {
                        z = z - id * table[i1 - 1];
                        if (z <= 0.0)
                        {
                            break;
                        }
                    }
                }

                mult[id - 1] = mult[id - 1] + j;
                npart = npart + j;
                m = i1;
            }

            //
            //  Reformulate the partition in the standard form.
            //  NPART is the number of distinct parts.
            //
            npart = 0;

            for (i = 1; i <= n; i++)
            {
                if (mult[i - 1] != 0)
                {
                    npart = npart + 1;
                    a[npart - 1] = i;
                    mult[npart - 1] = mult[i - 1];
                }
            }

            for (i = npart + 1; i <= n; i++)
            {
                mult[i - 1] = 0;
            }

        }

        public static void i4_partitions_next(int s, ref int[] m)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    I4_PARTITIONS_NEXT: next partition into S parts.
            //
            //  Discussion:
            //
            //    This function generates, one at a time, entries from the list of
            //    nondecreasing partitions of the integers into S or fewer parts.
            //
            //    The list is ordered first by the integer that is partitioned
            //    (the sum of the entries), and second by decreasing lexical order
            //    in the partition vectors.
            //
            //    The first value returned is the only such partition of 0.
            //
            //    Next comes the only partition of 1.
            //
            //    There follow two partitions of 2, and so on.
            //
            //    Typical use of this function begins with an initialization call,
            //    and then repeated calls in which the output from the previous call
            //    is used as input to the next call:
            //
            //    m = [ 0, 0, 0 ];
            //
            //    while ( condition )
            //      m = i4_partitions_next ( s, m );
            //    end
            //
            //  Example:
            //
            //    S = 3
            //
            //    P  D    M
            //    _  _  _____
            //    1  0  0 0 0
            //    2  1  1 0 0
            //    3  2  2 0 0
            //    4  2  1 1 0
            //    5  3  3 0 0
            //    6  3  2 1 0
            //    7  3  1 1 1
            //    8  4  4 0 0
            //    9  4  3 1 0
            //   10  4  2 2 0
            //   11  4  2 1 1
            //   12  5  5 0 0
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 August 2010
            //
            //  Author:
            //
            //    Original MATLAB version by Alan Genz.
            //    C++ version by John Burkardt.
            //
            //  Parameters:
            //
            //    Input, int S, the number of items in the partition.
            //
            //    Input/output, int M[S].  On input, the current partition.  
            //    On first call, this should be a nondecreasing partition.  Thereafter, it 
            //    should be the output partition from the previous call.  On output, the
            //    next partition.
            //
        {
            int i;
            int j;
            int msum;

            msum = m[0];

            for (i = 1; i < s; i++)
            {
                msum = msum + m[i];

                if (m[0] <= m[i] + 1)
                {
                    m[i] = 0;
                }
                else
                {
                    m[0] = msum - i * (m[i] + 1);
                    for (j = 1; j <= i; j++)
                    {
                        m[j] = m[i] + 1;
                    }

                    return;
                }
            }

            //
            //  If we failed to find a suitable index I, put
            //  the entire sum into M(1), increment by 1, and
            //  prepare to partition the next integer.
            //
            m[0] = msum + 1;

        }
    }
}