using System;
using Burkardt.Types;

namespace Burkardt.RankingNS;

public static partial class Ranking
{
    public static bool rgf_check(int m, int[] f)

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    RGF_CHECK checks a restricted growth function.
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
        //    Input, int M, the domain of the RGF is the integers
        //    from 1 to M.  M must be positive.
        // 
        //    Input, int F(M), the restricted growth function.
        // 
        //    Output, bool RGF_CHECK.
        //    TRUE, the data is legal.
        //    FALSE, the data is not legal.
    {
        bool check;
        int fmax;
        int i;

        check = true;

        switch (m)
        {
            case <= 0:
                check = false;
                return check;
        }

        fmax = 0;
        for (i = 0; i < m; i++)
        {
            if (f[i] <= 0 || fmax + 1 < f[i])
            {
                check = false;
                return check;
            }

            fmax = Math.Max(fmax, f[i]);
        }

        return check;
    }

    public static int rgf_enum(int m)

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    RGF_ENUM enumerates the restricted growth functions on M.
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
        //    Input, int M, the domain of the RGF is the integers
        //    from 1 to M.  M must be positive.  However, for the enumeration routine
        //    only, it is legal to call with any value of M.
        // 
        //    Output, int RGF_ENUM, the number of restricted growth
        //    functions.
        // 
    {
        int[] b;
        int i;
        int j;
        int value;

        switch (m)
        {
            case < 0:
                value = 0;
                break;
            case 0:
                value = 1;
                break;
            default:
            {
                b = new int[m + 1];
                b[0] = 1;
                for (j = 1; j <= m; j++)
                {
                    b[j] = 0;
                    for (i = 0; i < j; i++)
                    {
                        b[j] += typeMethods.i4_choose(j - 1, i) * b[i];
                    }
                }

                value = b[m];
                break;
            }
        }

        return value;
    }

    public static int[] rgf_g_table(int m)

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    RGF_G_TABLE tabulates the generalized restricted growth functions.
        // 
        //  Example:
        // 
        //    M = 6
        // 
        //    D =  1    1    1    1    1    1    1
        //         1    2    3    4    5    6    0
        //         2    5   10   17   26    0    0
        //         5   15   37   77    0    0    0
        //        15   52  151    0    0    0    0
        //        52  203    0    0    0    0    0
        //       203    0    0    0    0    0    0
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
        //    Input, int M, indicates how many rows and columns are to
        //    be computed.  M must be nonnegative.
        // 
        //    Output, int RGF_G_TABLE[(M+1)*(M+1)], the first M+1 rows and
        //    M+1 columns of the table of the number of generalized restricted growth
        //    functions.  D(I,J) is the number of GRGF''s of length I with restriction
        //    parameter J.
        // 
    {
        int[] d;
        int i;
        int j;

        d = new int[(m + 1) * (m + 1)];

        for (j = 0; j <= m; j++)
        {
            d[0 + j * (m + 1)] = 1;
        }

        for (i = 1; i <= m; i++)
        {
            for (j = 0; j <= m; j++)
            {
                if (j <= m - i)
                {
                    d[i + j * (m + 1)] = j * d[i - 1 + j * (m + 1)] + d[i - 1 + (j + 1) * (m + 1)];
                }
                else
                {
                    d[i + j * (m + 1)] = 0;
                }
            }
        }

        return d;
    }

    public static int rgf_rank(int m, int[] f)

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    RGF_RANK ranks a restricted growth function.
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
        //    Input, int M, the domain of the RGF is the integers
        //    from 1 to M.  M must be positive.
        // 
        //    Input, int F[M], the restricted growth function.
        // 
        //    Output, int RGF_RANK, the rank of the restricted growth
        //    function.
        // 
    {
        bool check;
        int[] d;
        int i;
        int j;
        int rank;
        // 
        //  Check.
        // 
        check = rgf_check(m, f);

        switch (check)
        {
            case false:
                Console.WriteLine("");
                Console.WriteLine("RGF_RANK - Fatal error!");
                Console.WriteLine("  The input array is illegal!");
                return 1;
        }

        // 
        //  Get the generalized restricted growth function table.
        // 
        d = rgf_g_table(m);

        rank = 0;
        j = 1;
        for (i = 2; i <= m; i++)
        {
            rank += (f[i - 1] - 1) * d[m - i + j * (m + 1)];
            j = Math.Max(j, f[i - 1]);
        }

        return rank;
    }

    public static void rgf_successor(int m, ref int[] f, ref int rank)

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    RGF_SUCCESSOR generates the next restricted growth function.
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
        //    Input, int M, the domain of the RGF is the integers
        //    from 1 to M.  M must be positive.
        // 
        //    Input/output, int F[M], the restricted growth function.
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
        int fmax;
        int i;
        int j;
        switch (rank)
        {
            // 
            //  Return the first element.
            // 
            case -1:
            {
                for (i = 0; i < m; i++)
                {
                    f[i] = 1;
                }

                rank = 0;
                return;
            }
        }

        // 
        //  Check.
        // 
        check = rgf_check(m, f);

        switch (check)
        {
            case false:
                Console.WriteLine("");
                Console.WriteLine("RGF_SUCCESSOR - Fatal error!");
                Console.WriteLine("  The input array is illegal!");
                return;
        }

        // 
        //  Find the first position from the right which can be incremented.
        // 
        for (i = m; 2 <= i; i--)
        {
            fmax = 1;
            for (j = 2; j < i; j++)
            {
                fmax = Math.Max(fmax, f[j - 1]);
            }

            // 
            //  Increment the function at this position, and set later entries to 1.
            // 
            if (f[i - 1] != fmax + 1)
            {
                f[i - 1] += 1;
                for (j = i + 1; j <= m; j++)
                {
                    f[j - 1] = 1;
                }

                rank += 1;
                return;
            }
        }

        // 
        //  The final element was input.
        //  Return the first element.
        // 
        for (i = 0; i < m; i++)
        {
            f[i] = 1;
        }

        rank = 0;
    }

    public static void rgf_to_setpart(int m, int[] f, ref int nsub, ref int[] s, ref int[] index)

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    RGF_TO_SETPART converts a restricted growth function to a set partition.
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
        //    Input, int M, the domain of the RGF is the integers
        //    from 1 to M.  M must be positive.
        // 
        //    Input, int F[M], the restricted growth function.
        // 
        //    Output, int NSUB, the number of nonempty subsets into
        //    which the set is partitioned.
        // 
        //    Output, int S[M], describes the partition of a set of
        //    M objects into NSUB nonempty subsets.  If element I of the
        //    superset belongs to subset J, then S(I) = J.
        // 
        //    Output, int INDEX[M], lists the location in S of the last
        //    element of each subset.  Thus, the elements of subset 1
        //    are S(1) through S(INDEX(1)), the elements of subset 2
        //    are S(INDEX(1)+1) through S(INDEX(2)) and so on.
        // 
    {
        bool check;
        int i;
        int j;
        int k;
        // 
        //  Check.
        // 
        check = rgf_check(m, f);

        switch (check)
        {
            case false:
                Console.WriteLine("");
                Console.WriteLine("RGF_TO_SETPART - Fatal error!");
                Console.WriteLine("  The input array is illegal!");
                return;
        }

        // 
        //  Determine the number of subsets.
        // 
        nsub = typeMethods.i4vec_max(m, f);
        // 
        //  Initialize.
        //
        for (i = 0; i < m; i++)
        {
            s[i] = 0;
        }

        for (i = 0; i < nsub; i++)
        {
            index[i] = 0;
        }

        // 
        //  For each subset I, collect the indices of F which have value I.
        //  These are the elements of the I-th subset.
        // 
        k = 0;
        for (i = 1; i <= nsub; i++)
        {
            for (j = 1; j <= m; j++)
            {
                if (f[j - 1] == i)
                {
                    k += 1;
                    s[k - 1] = j;
                }
            }

            index[i - 1] = k;
        }
    }

    public static int[] rgf_unrank(int rank, int m)

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    RGF_UNRANK returns the restricted growth function of a given rank.
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
        //    Input, int RANK, the rank of the restricted growth
        //    function.
        // 
        //    Input, int M, the domain of the RGF is the integers
        //    from 1 to M.  M must be positive.
        // 
        //    Output, int RGF_UNRANK[M], the restricted growth function.
        // 
    {
        int[] d;
        int[] f;
        int i;
        int j;
        int nrgf;
        int rank_copy;
        switch (m)
        {
            // 
            //  Check.
            // 
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("RGF_UNRANK - Fatal error!");
                Console.WriteLine("  Input M is illegal.");
                return null;
        }

        nrgf = rgf_enum(m);

        if (rank < 0 || nrgf < rank)
        {
            Console.WriteLine("");
            Console.WriteLine("RGF_UNRANK - Fatal error!");
            Console.WriteLine("  The input rank is illegal.");
            return null;
        }

        // 
        //  Get the generalized restricted growth function table.
        // 
        d = rgf_g_table(m);

        f = new int[m];

        rank_copy = rank;
        j = 1;
        f[0] = 1;

        for (i = 2; i <= m; i++)
        {
            if (j * d[m - i + j * (m + 1)] <= rank_copy)
            {
                f[i - 1] = j + 1;
                rank_copy -= j * d[m - i + j * (m + 1)];
                j += 1;
            }
            else
            {
                f[i - 1] = 1 + rank_copy / d[m - i + j * (m + 1)];
                rank_copy %= d[m - i + j * (m + 1)];
            }
        }

        return f;
    }
}