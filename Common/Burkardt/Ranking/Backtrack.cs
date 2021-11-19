namespace Burkardt.RankingNS;

public static partial class Ranking
{
    public static void backtrack(int l, ref int[] iarray, ref int indx, ref int k, ref int nstack,
            ref int[] stack, int maxstack )

        //****************************************************************************80
        // 
        //  Purpose:
        //
        //    BACKTRACK supervises a backtrack search.
        // 
        //  Discussion:
        // 
        //    The routine builds a vector, one element at a time, which is
        //    required to satisfy some condition.
        // 
        //    At any time, the partial vector may be discovered to be
        //    unsatisfactory, but the routine records information about where the
        //    last arbitrary choice was made, so that the search can be
        //    carried out efficiently, rather than starting out all over again.
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
        //    Input, int L, the length of the completed candidate vector.
        // 
        //    Input/output, int IARRAY[L], the candidate vector.
        // 
        //    Input/output, int &INDX.
        //    On input, set INDX = 0 to start a search.
        //    On output:
        //    1, a complete output vector has been determined.
        //    2, candidates are needed.
        //    3, no more possible vectors exist.
        // 
        //    Input/output, int &K, the current length of the candidate
        //    vector.
        // 
        //    Input/output, int &NSTACK, the current length of the stack.
        // 
        //    Input/output, int STACK[MAXSTACK], a list of candidates
        //    for positions 1 through K.
        // 
        //    Input, int MAXSTACK, the maximum length of the stack.
        // 
    {
        switch (indx)
        {
            // 
            //  If this is the first call, request a candidate for position 1.
            // 
            case 0:
                k = 1;
                nstack = 0;
                indx = 2;
                return;
        }

        // 
        //  Examine the stack.
        // 
        for (;;)
        {
            nstack -= 1;
            // 
            //  If there are candidates for position K, take the first available
            //  one off the stack, and increment K.
            // 
            //  This may cause K to reach the desired value of L, in which case
            //  we need to signal the user that a complete set of candidates
            //  is being returned.
            // 
            if (stack[nstack] != 0)
            {
                iarray[k - 1] = stack[nstack - 1];
                stack[nstack - 1] = stack[nstack] - 1;

                if (k != l)
                {
                    k += 1;
                    indx = 2;
                }
                else
                {
                    indx = 1;
                }

                break;
            }
            // 
            //  If there are no candidates for position K, then decrement K.
            //  If K is still positive, repeat the examination of the stack.
            // 

            k -= 1;

            if (k > 0)
            {
                continue;
            }

            indx = 3;
            break;
        }
    }
}