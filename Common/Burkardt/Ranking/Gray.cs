using System;
using Burkardt.Types;

namespace Burkardt.RankingNS
{
    public static partial class Ranking
    {
        public static bool gray_code_check(int n, int[] t)

            //****************************************************************************80
            // 
            //  Purpose:
            //
            //    GRAY_CODE_CHECK checks a Gray code element.
            // 
            //  Licensing:
            // 
            //    This code is distributed under the GNU LGPL license.
            // 
            //  Modified:
            // 
            //    25 November 2015
            // 
            //  Author:
            // 
            //    John Burkardt
            // 
            //  Parameters:
            // 
            //    Input, int N, the number of digits in each element.
            //    N must be positive.
            // 
            //    Input, int T[N], an element of the Gray code.
            //    Each entry T(I) is either 0 or 1.
            // 
            //    Output, bool GRAY_CODE_CHECK.
            //    TRUE, the data is legal.
            //    FALSE, the data is not legal.
            // 
        {
            int i;

            bool check = true;

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

        public static int gray_code_enum(int n)

            //****************************************************************************80
            // 
            //  Purpose:
            //
            //    GRAY_CODE_ENUM enumerates the Gray codes on N digits.
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
            //    Input, int N, the number of digits in each element.
            //    N must be nonnegative.
            // 
            //    Output, int GRAY_CODE_ENUM, the number of distinct elements.
            // 
        {
            int value = (int)Math.Pow(2, n);

            return value;
        }

        public static int gray_code_rank(int n, int[] t)

            //****************************************************************************80
            // 
            //  Purpose:
            //
            //    GRAY_CODE_RANK computes the rank of a Gray code element.
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
            //    Input, int N, the number of digits in each element.
            //    N must be positive.
            // 
            //    Input, int T[N], an element of the Gray code.
            //    Each entry is either 0 or 1.
            // 
            //    Output, int GRAY_CODE_RANK, the rank of the element.
            // 
        {
            int b;
            bool check;
            int i;
            int rank;
            // 
            //  Check.
            // 
            check = gray_code_check(n, t);

            if (!check)
            {
                Console.WriteLine("");
                Console.WriteLine("GRAY_CODE_RANK - Fatal error!");
                Console.WriteLine("  The input array is illegal.");
                return (1);
            }

            rank = 0;
            b = 0;

            for (i = n - 1; 0 <= i; i--)
            {
                if (t[n - i - 1] != 0)
                {
                    b = 1 - b;
                }

                if (b == 1)
                {
                    rank = rank + (int)Math.Pow(2, i);
                }
            }

            return rank;
        }

        public static void gray_code_successor(int n, ref int[] t, ref int rank)

            //****************************************************************************80
            // 
            //  Purpose:
            //
            //    GRAY_CODE_SUCCESSOR computes the binary reflected Gray code successor.
            // 
            //  Example:
            // 
            //    000, 001, 011, 010, 110, 111, 101, 100,
            //    after which the sequence repeats.
            // 
            //  Discussion:
            // 
            //    In the original code, the successor of the element that has an
            //    initial 1 followed by N-1 zeroes is undefined.  In this version,
            //    the successor is the element with N zeroes.
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
            //    Input, int N, the number of digits in each element.
            //    N must be positive.
            // 
            //    Input/output, int T[N].
            //    On input, T contains an element of the Gray code, that is,
            //    each entry T(I) is either 0 or 1.
            //    On output, T contains the successor to the input value; this
            //    is an element of the Gray code, which differs from the input
            //    value in a single position.
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
            int weight;
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
            check = gray_code_check(n, t);

            if (!check)
            {
                Console.WriteLine("");
                Console.WriteLine("GRAY_CODE_SUCCESSOR - Fatal error!");
                Console.WriteLine("  The input array is illegal.");
                return;
            }

            weight = typeMethods.i4vec_sum(n, t);

            if ((weight % 2) == 0)
            {
                if (t[n - 1] == 0)
                {
                    t[n - 1] = 1;
                }
                else
                {
                    t[n - 1] = 0;
                }

                rank = rank + 1;
                return;
            }
            else
            {
                for (i = n - 1; 1 <= i; i--)
                {
                    if (t[i] == 1)
                    {
                        if (t[i - 1] == 0)
                        {
                            t[i - 1] = 1;
                        }
                        else
                        {
                            t[i - 1] = 0;
                        }

                        rank = rank + 1;
                        return;
                    }
                }

                // 
                //  The final element was input.
                //  Return the first element.
                // 
                for (i = 0; i < n; i++)
                {
                    t[i] = 0;
                }

                rank = 0;
            }
        }

        public static int[] gray_code_unrank(int rank, int n)

            //****************************************************************************80
            // 
            //  Purpose:
            //
            //    GRAY_CODE_UNRANK computes the Gray code element of given rank.
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
            //    Input, int RANK, the rank of the element.
            //    0 <= RANK <= 2^N.
            // 
            //    Input, int N, the number of digits in each element.
            //    N must be positive.
            // 
            //    Output, int GRAY_CODE_UNRANK[N], the element of the Gray code which has
            //    the given rank.
            // 
        {
            int b;
            int bprime;
            int i;
            int ngray;
            int rank_copy;
            int[] t;
            // 
            //  Check.
            // 
            if (n < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("GRAY_CODE_UNRANK - Fatal error!");
                Console.WriteLine("  Input N is illegal.");
                return null;
            }

            ngray = gray_code_enum(n);

            if (rank < 0 || ngray < rank)
            {
                Console.WriteLine("");
                Console.WriteLine("GRAY_CODE_UNRANK - Fatal error!");
                Console.WriteLine("  The input rank is illegal.");
                return null;
            }

            t = new int[n];

            rank_copy = rank;
            for (i = 0; i < n; i++)
            {
                t[i] = 0;
            }

            bprime = 0;

            for (i = n - 1; 0 <= i; i--)
            {
                b = rank_copy / (int)Math.Pow(2, i);

                if (b != bprime)
                {
                    t[n - i - 1] = 1;
                }

                bprime = b;
                rank_copy = rank_copy - b * (int)Math.Pow(2, i);
            }

            return t;
        }

        public static void gray_next(int n, ref int change, ref int k, ref int[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GRAY_NEXT generates the next Gray code by switching one item at a time.
            //
            //  Discussion:
            //
            //    On the first call only, the user must set CHANGE = -N.
            //    This initializes the routine to the Gray code for N zeroes.
            //
            //    Each time it is called thereafter, it returns in CHANGE the index
            //    of the item to be switched in the Gray code.  The sign of CHANGE
            //    indicates whether the item is to be added or subtracted (or
            //    whether the corresponding bit should become 1 or 0).  When
            //    CHANGE is equal to N+1 on output, all the Gray codes have been
            //    generated.
            //
            //  Example:
            //
            //    N  CHANGE         Subset in/out   Binary Number
            //                      Interpretation  Interpretation
            //                       1 2 4 8
            //   --  ---------      --------------  --------------
            //
            //    4   -4 / 0         0 0 0 0         0
            //
            //        +1             1 0 0 0         1
            //          +2           1 1 0 0         3
            //        -1             0 1 0 0         2
            //            +3         0 1 1 0         6
            //        +1             1 1 1 0         7
            //          -2           1 0 1 0         5
            //        -1             0 0 1 0         4
            //              +4       0 0 1 1        12
            //        +1             1 0 1 1        13
            //          +2           1 1 1 1        15
            //        -1             0 1 1 1        14
            //            -3         0 1 0 1        10
            //        +1             1 1 0 1        11
            //          -2           1 0 0 1         9
            //        -1             0 0 0 1         8
            //              -4       0 0 0 0         0
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
            //    Input/output, int &CHANGE.  This item is used for input only
            //    on the first call for a particular sequence of Gray codes,
            //    at which time it must be set to -N.  This corresponds to
            //    all items being excluded, or all bits being 0, in the Gray code.
            //    On output, CHANGE indicates which of the N items must be "changed", 
            //    and the sign indicates whether the item is to be added or removed
            //    (or the bit is to become 1 or 0).  Note that on return from the 
            //    first call, CHANGE is set to 0, indicating that we begin with
            //    the empty set.
            //
            //    Input/output, int &K, a bookkeeping variable.
            //    The user must declare this variable before the first call.
            //    The output value from one call should be the input value for the next.
            //    The user should not change this variable.
            //
            //    Input/output, int A[N], a bookkeeping variable.
            //    The user must declare this variable before the first call.
            //    The output value from one call should be the input value for the next.
            //    The user should not change this variable.
            //
        {
            int i;

            if (n <= 0)
            {
                Console.WriteLine("");
                Console.WriteLine("GRAY_NEXT - Fatal error!");
                Console.WriteLine("  Input value of N <= 0.");
                return;
            }

            if (change == -n)
            {
                change = 0;
                k = 1;
                for (i = 0; i < n; i++)
                {
                    a[i] = 0;
                }

                return;
            }

            //
            //  First determine WHICH item is to be changed.
            //
            if ((k % 2) == 1)
            {
                change = 1;
            }
            else
            {
                for (i = 1; i <= n; i++)
                {
                    if (a[i - 1] == 1)
                    {
                        change = i + 1;
                        break;
                    }
                }
            }

            //
            //  Take care of the terminal case CHANGE = N + 1.
            //
            if (change == n + 1)
            {
                change = n;
            }

            //
            //  Now determine HOW the item is to be changed.
            //
            if (a[change - 1] == 0)
            {
                a[change - 1] = 1;
            }
            else if (a[change - 1] == 1)
            {
                a[change - 1] = 0;
                change = -(change);
            }

            //
            //  Update the counter.
            //
            k = k + 1;
            //
            //  If the output CHANGE = -N, then we're done.
            //
            if (change == -n)
            {
                k = 0;
            }

        }

        public static int gray_rank(int gray)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GRAY_RANK ranks a Gray code.
            //
            //  Discussion:
            //
            //    Given the number GRAY, its ranking is the order in which it would be
            //    visited in the Gray code ordering.  The Gray code ordering begins
            //
            //    Rank  Gray  Gray
            //          (Dec) (Bin)
            //
            //       0     0  0000
            //       1     1  0001
            //       2     3  0011
            //       3     2  0010
            //       4     6  0110
            //       5     7  0111
            //       6     5  0101
            //       7     4  0100
            //       8    12  0110
            //       etc
            //
            //   This routine is given a Gray code, and has to return the rank.
            //
            //  Example:
            //
            //    Gray  Gray  Rank
            //    (Dec) (Bin)
            //
            //     0       0     0
            //     1       1     1
            //     2      10     3
            //     3      11     2
            //     4     100     7
            //     5     101     6
            //     6     110     4
            //     7     111     5
            //     8    1000    15
            //     9    1001    14
            //    10    1010    12
            //    11    1011    13
            //    12    1100     8
            //    13    1101     9
            //    14    1110    11
            //    15    1111    10
            //    16   10000    31
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 May 2003
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
            //    Input, int GRAY, the Gray code to be ranked.
            //
            //    Output, int GRAY_RANK, the rank of GRAY, and the integer whose Gray
            //    code is GRAY.
            //
        {
            int i;
            int nbits = 32;
            int rank;

            rank = 0;

            if (typeMethods.i4_btest(gray, nbits - 1))
            {
                rank = typeMethods.i4_bset(rank, nbits - 1);
            }

            for (i = nbits - 2; 0 <= i; i--)
            {
                if (typeMethods.i4_btest(rank, i + 1) != typeMethods.i4_btest(gray, i))
                {
                    rank = typeMethods.i4_bset(rank, i);
                }
            }

            return rank;
        }

        public static int gray_rank2(int gray)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GRAY_RANK2 ranks a Gray code.
            //
            //  Discussion:
            //
            //    In contrast to GRAY_RANK, this routine is entirely arithmetical,
            //    and does not require access to bit testing and setting routines.
            //
            //
            //    Given the number GRAY, its ranking is the order in which it would be
            //    visited in the Gray code ordering.  The Gray code ordering begins
            //
            //    Rank  Gray  Gray
            //          (Dec) (Bin)
            //
            //       0     0  0000
            //       1     1  0001
            //       2     3  0011
            //       3     2  0010
            //       4     6  0110
            //       5     7  0111
            //       6     5  0101
            //       7     4  0100
            //       8    12  0110
            //       etc
            //
            //   This routine is given a Gray code, and has to return the rank.
            //
            //  Example:
            //
            //    Gray  Gray  Rank
            //    (Dec) (Bin)
            //
            //     0       0     0
            //     1       1     1
            //     2      10     3
            //     3      11     2
            //     4     100     7
            //     5     101     6
            //     6     110     4
            //     7     111     5
            //     8    1000    15
            //     9    1001    14
            //    10    1010    12
            //    11    1011    13
            //    12    1100     8
            //    13    1101     9
            //    14    1110    11
            //    15    1111    10
            //    16   10000    31
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    17 July 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int GRAY, the Gray code to be ranked.
            //
            //    Output, int GRAY_RANK, the rank of GRAY, and the integer whose Gray
            //    code is GRAY.
            //
        {
            int k;
            bool last;
            bool next;
            int rank;
            int two_k;

            if (gray < 0)
            {
                Console.WriteLine("");
                Console.WriteLine("GRAY_RANK2 - Fatal error!");
                Console.WriteLine("  Input value of GRAY < 0.");
                return (1);
            }

            if (gray == 0)
            {
                rank = 0;
                return rank;
            }

            //
            //  Find TWO_K, the largest power of 2 less than or equal to GRAY.
            //
            k = 0;
            two_k = 1;
            while (2 * two_k <= gray)
            {
                two_k = two_k * 2;
                k = k + 1;
            }

            rank = two_k;
            last = true;
            gray = gray - two_k;

            while (0 < k)
            {
                two_k = two_k / 2;
                k = k - 1;

                next = (two_k <= gray && gray < two_k * 2);

                if (next)
                {
                    gray = gray - two_k;
                }

                if (next != last)
                {
                    rank = rank + two_k;
                    last = true;
                }
                else
                {
                    last = false;
                }
            }

            return rank;
        }

        public static int gray_unrank(int rank)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GRAY_UNRANK unranks a Gray code.
            //
            //  Discussion:
            //
            //    The binary values of the Gray codes of successive integers differ in
            //    just one bit.
            //
            //    The sequence of Gray codes for 0 to (2**N)-1 can be interpreted as a
            //    Hamiltonian cycle on a graph of the cube in N dimensions.
            //
            //  Example:
            //
            //    Rank  Gray  Gray
            //          (Dec) (Bin)
            //
            //     0     0       0
            //     1     1       1
            //     2     3      11
            //     3     2      10
            //     4     6     110
            //     5     7     111
            //     6     5     101
            //     7     4     100
            //     8    12    1100
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 May 2003
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
            //    Input, int RANK, the integer whose Gray code is desired.
            //
            //    Output, int GRAY_UNRANK, the Gray code of the given rank.
            //
        {
            int gray;
            int i;
            int nbits = 32;

            gray = 0;

            if (typeMethods.i4_btest(rank, nbits - 1))
            {
                gray = typeMethods.i4_bset(gray, nbits - 1);
            }

            for (i = nbits - 2; 0 <= i; i--)
            {
                if (typeMethods.i4_btest(rank, i + 1) != typeMethods.i4_btest(rank, i))
                {
                    gray = typeMethods.i4_bset(gray, i);
                }
            }

            return gray;
        }

        public static int gray_unrank2(int rank)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    GRAY_UNRANK2 unranks a Gray code.
            //
            //  Discussion:
            //
            //    In contrast to GRAY_UNRANK, this routine is entirely arithmetical,
            //    and does not require access to bit testing and setting routines.
            //
            //    The binary values of the Gray codes of successive integers differ in
            //    just one bit.
            //
            //    The sequence of Gray codes for 0 to (2**N)-1 can be interpreted as a
            //    Hamiltonian cycle on a graph of the cube in N dimensions.
            //
            //  Example:
            //
            //    Rank  Gray  Gray
            //          (Dec) (Bin)
            //
            //     0     0       0
            //     1     1       1
            //     2     3      11
            //     3     2      10
            //     4     6     110
            //     5     7     111
            //     6     5     101
            //     7     4     100
            //     8    12    1100
            //     9    14    1001
            //    10    12    1010
            //    11    13    1011
            //    12     8    1100
            //    13     9    1101
            //    14    11    1110
            //    15    10    1111
            //    16    31   10000
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    17 July 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int RANK, the integer whose Gray code is desired.
            //
            //    Output, int GRAY_UNRANK2, the Gray code of the given rank.
            //
        {
            int gray;
            int k;
            bool last;
            bool next;
            int two_k;

            if (rank <= 0)
            {
                gray = 0;
                return gray;
            }

            k = 0;
            two_k = 1;
            while (2 * two_k <= rank)
            {
                two_k = two_k * 2;
                k = k + 1;
            }

            gray = two_k;
            rank = rank - two_k;
            next = true;

            while (0 < k)
            {
                two_k = two_k / 2;
                k = k - 1;

                last = next;
                next = (two_k <= rank && rank <= two_k * 2);

                if (next != last)
                {
                    gray = gray + two_k;
                }

                if (next)
                {
                    rank = rank - two_k;
                }
            }

            return gray;
        }


    }
}