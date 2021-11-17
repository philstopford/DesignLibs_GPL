using System;

namespace Burkardt.RankingNS;

public static partial class Ranking
{
    public static int[] stirling_numbers1(int m, int n)

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    STIRLING_NUMBERS1 computes Stirling numbers of the first kind.
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
        //    Input, int M, the maximum row to compute.
        //    M must be nonnegative.
        // 
        //    Input, int N, the maximum column to compute.
        //    N must be nonnegative.
        // 
        //    Output, int S[(M+1)*(N+1)], the first M+1 rows and N+1 columns
        //    of the table of Stirling numbers of the first kind.
        // 
    {
        int i;
        int j;
        int[] s;

        s = new int[(m + 1) * (n + 1)];

        s[0 + 0 * (m + 1)] = 1;
        for (j = 1; j <= n; j++)
        {
            s[0 + j * (m + 1)] = 0;
        }

        for (i = 1; i <= m; i++)
        {
            s[i + 0 * (m + 1)] = 0;
        }

        // 
        //  This loop may be extraneous.
        // 
        for (i = 0; i <= Math.Min(m, n - 1); i++)
        {
            s[i + (i + 1) * (m + 1)] = 0;
        }

        for (i = 1; i <= m; i++)
        {
            for (j = 1; j <= n; j++)
            {
                if (j <= i)
                {
                    s[i + j * (m + 1)] = s[i - 1 + (j - 1) * (m + 1)] - (i - 1) * s[i - 1 + j * (m + 1)];
                }
                else
                {
                    s[i + j * (m + 1)] = 0;
                }
            }
        }

        return s;
    }

    public static int[] stirling_numbers2(int m, int n)

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    STIRLING_NUMBERS2 computes Stirling numbers of the second kind.
        // 
        //  Discussion:
        // 
        //    The reference has a typographical error, referring to
        //    S(I-J,J-1) instead of S(I-1,J-1).
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
        //    Input, int M, the maximum row to compute.
        //    M must be nonnegative.
        // 
        //    Input, int N, the maximum column to compute.
        //    N must be nonnegative.
        // 
        //    Output, int S[(M+1)*(N+1)], the first M+1 rows and N+1 columns
        //    of the table of Stirling numbers of the second kind.
        // 
    {
        int i;
        int j;
        int[] s;

        s = new int[(m + 1) * (n + 1)];

        s[0 + 0 * (m + 1)] = 1;
        for (j = 1; j <= n; j++)
        {
            s[0 + j * (m + 1)] = 0;
        }

        for (i = 1; i <= m; i++)
        {
            s[i + 0 * (m + 1)] = 0;
        }

        // 
        //  This loop may be extraneous.
        // 
        for (i = 0; i <= Math.Min(m, n - 1); i++)
        {
            s[i + (i + 1) * (m + 1)] = 0;
        }

        for (i = 1; i <= m; i++)
        {
            for (j = 1; j <= n; j++)
            {
                if (j <= i)
                {
                    s[i + j * (m + 1)] = j * s[i - 1 + j * (m + 1)] + s[i - 1 + (j - 1) * (m + 1)];
                }
                else
                {
                    s[i + j * (m + 1)] = 0;
                }
            }
        }

        return s;
    }
}