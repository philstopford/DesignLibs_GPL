using System;
using Burkardt.Types;

namespace Burkardt.RankingNS;

public static partial class Ranking
{


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
        int i;
        // 
        //  Check.
        // 
        bool check = typeMethods.gray_code_check(n, t);

        switch (check)
        {
            case false:
                Console.WriteLine("");
                Console.WriteLine("GRAY_CODE_RANK - Fatal error!");
                Console.WriteLine("  The input array is illegal.");
                return 1;
        }

        int rank = 0;
        int b = 0;

        for (i = n - 1; 0 <= i; i--)
        {
            if (t[n - i - 1] != 0)
            {
                b = 1 - b;
            }

            switch (b)
            {
                case 1:
                    rank += (int)Math.Pow(2, i);
                    break;
            }
        }

        return rank;
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
        int i;
        switch (n)
        {
            // 
            //  Check.
            // 
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("GRAY_CODE_UNRANK - Fatal error!");
                Console.WriteLine("  Input N is illegal.");
                return null;
        }

        int ngray = typeMethods.gray_code_enum(n);

        if (rank < 0 || ngray < rank)
        {
            Console.WriteLine("");
            Console.WriteLine("GRAY_CODE_UNRANK - Fatal error!");
            Console.WriteLine("  The input rank is illegal.");
            return null;
        }

        int[] t = new int[n];

        int rank_copy = rank;
        for (i = 0; i < n; i++)
        {
            t[i] = 0;
        }

        int bprime = 0;

        for (i = n - 1; 0 <= i; i--)
        {
            int b = rank_copy / (int)Math.Pow(2, i);

            if (b != bprime)
            {
                t[n - i - 1] = 1;
            }

            bprime = b;
            rank_copy -= b * (int)Math.Pow(2, i);
        }

        return t;
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
        const int nbits = 32;

        int rank = 0;

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
        int rank;

        switch (gray)
        {
            case < 0:
                Console.WriteLine("");
                Console.WriteLine("GRAY_RANK2 - Fatal error!");
                Console.WriteLine("  Input value of GRAY < 0.");
                return 1;
            case 0:
                rank = 0;
                return rank;
        }

        //
        //  Find TWO_K, the largest power of 2 less than or equal to GRAY.
        //
        int k = 0;
        int two_k = 1;
        while (2 * two_k <= gray)
        {
            two_k *= 2;
            k += 1;
        }

        rank = two_k;
        bool last = true;
        gray -= two_k;

        while (0 < k)
        {
            two_k /= 2;
            k -= 1;

            bool next = two_k <= gray && gray < two_k * 2;

            switch (next)
            {
                case true:
                    gray -= two_k;
                    break;
            }

            if (next != last)
            {
                rank += two_k;
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
        int i;
        const int nbits = 32;

        int gray = 0;

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

        switch (rank)
        {
            case <= 0:
                gray = 0;
                return gray;
        }

        int k = 0;
        int two_k = 1;
        while (2 * two_k <= rank)
        {
            two_k *= 2;
            k += 1;
        }

        gray = two_k;
        rank -= two_k;
        bool next = true;

        while (0 < k)
        {
            two_k /= 2;
            k -= 1;

            bool last = next;
            next = two_k <= rank && rank <= two_k * 2;

            if (next != last)
            {
                gray += two_k;
            }

            switch (next)
            {
                case true:
                    rank -= two_k;
                    break;
            }
        }

        return gray;
    }

    public static int vec_gray_rank(int n, int[] base_, ref int[] a )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VEC_GRAY_RANK computes the rank of a product space element.
        //
        //  Discussion:
        //
        //    The rank applies only to the elements as produced by the routine
        //    VEC_GRAY_NEXT.
        //
        //  Example:
        //
        //    N = 2, BASE = ( 2, 3 ), A = ( 1, 2 ),
        //
        //    RANK = 4.
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
        //    Input, int N, the number of components.
        //
        //    Input, int BASE[N], contains the number of degrees of
        //    freedom of each component.  The output values of A will
        //    satisfy 0 <= A[I] < BASE[I].
        //
        //    Input, int A[N], the product space element, with the
        //    property that 0 <= A[I] < BASE[I] for each entry I.
        //
        //    Output, int VEC_RANK, the rank, or order, of the element in
        //    the list of all elements.  The rank count begins at 1.
        //
    {
        int i;

        int rank = 0;

        for (i = 0; i < n; i++)
        {
            int c = (rank % 2) switch
            {
                1 => base_[i] - a[i] - 1,
                _ => a[i]
            };

            rank = base_[i] * rank + c;
        }

        rank += 1;

        return rank;
    }

    public static void vec_gray_unrank(int n, int[] base_, int rank, ref int[] a )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    VEC_GRAY_UNRANK computes the product space element of a given rank.
        //
        //  Discussion:
        //
        //    The rank applies only to the elements as produced by the routine
        //    VEC_GRAY_NEXT.
        //
        //  Example:
        //
        //    N = 2, BASE = ( 2, 3 ), RANK = 4.
        //
        //    A = ( 1, 2 ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 October 2006
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
        //    Input, int N, the number of components.
        //
        //    Input, int BASE[N], contains the number of degrees of
        //    freedom of each component.  The output values of A will
        //    satisfy 0 <= A[I] < BASE[I].
        //
        //    Input, int RANK, the desired rank, or order, of the element in
        //    the list of all elements.  The rank count begins at 1 and extends
        //    to MAXRANK = Product ( 0 <= I <= N ) BASE[I].
        //
        //    Output, int A[N], the product space element of the given rank.
        //
    {
        int i;

        int s = rank - 1;

        for (i = n - 1; 0 <= i; i--)
        {
            a[i] = s % base_[i];
            s /= base_[i];

            a[i] = (s % 2) switch
            {
                1 => base_[i] - a[i] - 1,
                _ => a[i]
            };
        }
    }
}