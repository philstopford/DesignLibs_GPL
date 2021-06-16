namespace Burkardt.RankingNS
{
    public static partial class Ranking
    {
        public static void queens(int n, int[] iarray, int k, ref int nstack, ref int[] istack,
        int maxstack )

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    QUEENS finds possible positions for the K-th nonattacking queen.
        // 
        //  Discussion:
        // 
        //    The chessboard is N by N, and is being filled one column at a time,
        //    with a tentative solution to the nonattacking queen problem.  So
        //    far, K-1 rows have been chosen, and we now need to provide a list
        //    of all possible rows that might be used in column K.
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
        //  Parameters:
        // 
        //    Input, int N, the total number of queens to place, and
        //    the length of a side of the chessboard.
        // 
        //    Input, int IARRAY[N].  The first K-1 entries of IARRAY
        //    record the rows into which queens have already been placed.
        // 
        //    Input, int K, the column for which we need possible
        //    row positions for the next queen.
        // 
        //    Input/output, int &NSTACK, the current length of stack.
        //    On output, this has been updated.
        // 
        //    Input/output, int ISTACK[MAXSTACK].  On output, we
        //    have added the candidates, and the number of candidates, to the end
        //    of the stack.
        // 
        //    Input, int MAXSTACK, maximum dimension of ISTACK.
        // 
        {
            bool diag;
            int irow;
            int jcol;
            int ncan;
            bool row;

            ncan = 0;

            for (irow = 1; irow <= n; irow++)
            {
                // 
                //  If row IROW has already been used, that is it.
                // 
                row = false;

                for (jcol = 1; jcol <= k - 1; jcol++)
                {
                    if (iarray[jcol - 1] == irow)
                    {
                        row = true;
                    }
                }

                if (!row)
                {
                    diag = false;

                    for (jcol = 1; jcol <= k - 1; jcol++)
                    {
                        if (irow == iarray[jcol - 1] + k - jcol ||
                            irow == iarray[jcol - 1] - (k - jcol))
                        {
                            diag = true;
                        }
                    }

                    if (!diag)
                    {
                        ncan = ncan + 1;
                        istack[nstack] = irow;
                        nstack = nstack + 1;
                    }
                }
            }

            istack[nstack] = ncan;
            nstack = nstack + 1;
        }
    }
}