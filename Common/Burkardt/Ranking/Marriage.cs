namespace Burkardt.RankingNS;

public static partial class Ranking
{
    public static void marriage(int n, int[] prefer, int[] rank, ref int[] fiancee, ref int[] next )

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    MARRIAGE finds a stable set of marriages for given preferences.
        // 
        //  Discussion:
        // 
        //    Given a set of N men and N women who must be married in pairs,
        //    and information defining the relative rankings that each person
        //    assigns to the candidates of the opposite sex, this routine finds
        //    a stable set of marriages for them.
        // 
        //    A stable set of marriages is a pairing of the men and women with
        //    the stability property: if M1 marries W1 and M2 marries W2, then
        //    it is never the case that M1 and W2 would both prefer to be married
        //    to each other.
        // 
        //    An important application of stable marriage algorithms occurs in
        //    the annual matching of medical residents to hospitals.
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
        //    Robert Sedgewick,
        //    Algorithms in C,
        //    Addison-Wesley, 1990,
        //    ISBN: 0-201-51425-7,
        //    LC: QA76.73.C15S43.
        // 
        //  Parameters:
        // 
        //    Input, int N, the number of pairs of men and women.
        // 
        //    Input, int PREFER[N*N]; for man I, the value of
        //    PREFER(I,J) represents his J-th preference for a wife.
        // 
        //    Input, int RANK[N*N]; for woman I, the value of RANK(I,J)
        //    represents her ranking of man number J.  A value of 1 for RANK(I,J)
        //    means woman I ranks man J most preferable, while a value of N
        //    would mean she ranked him least preferable.
        // 
        //    Output, int FIANCEE[N]; for woman I, FIANCEE(I) is the
        //    man to whom she is now engaged.
        // 
        //    Output, int NEXT[N]; for man I, NEXT(I) is his preference
        //    ranking for the woman to whom he is now engaged.  A value of 1 represents
        //    his first choice, a value of N his last.
        // 
    {
        int i;
        int m;
        int temp;
        int w;
        // 
        //  For man I, NEXT(I) is the woman I has most recently proposed to,
        //  and hence NEXT(I)+1 is the next one to try.
        // 
        for (i = 0; i < n; i++)
        {
            next[i] = 0;
        }

        // 
        //  For woman I, FIANCEE(I) is the man she has agree to marry,
        //  or 0 if she has not agreed to any man yet.
        // 
        for (i = 0; i < n; i++)
        {
            fiancee[i] = -1;
        }

        // 
        //  Start with an unengaged man, and end with an engaged woman.
        // 
        for (i = 1; i <= n; i++)
        {
            m = i;

            for (;;)
            {
                next[m - 1] += 1;

                w = prefer[m - 1 + (next[m - 1] - 1) * n];

                if (fiancee[w - 1] == -1)
                {
                    fiancee[w - 1] = m;
                    break;
                }

                if (rank[w - 1 + (m - 1) * n] < rank[w - 1 + (fiancee[w - 1] - 1) * n])
                {
                    temp = fiancee[w - 1];
                    fiancee[w - 1] = m;
                    m = temp;
                }
            }
        }
    }
}