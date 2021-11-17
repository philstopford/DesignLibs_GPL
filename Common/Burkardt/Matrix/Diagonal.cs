using System;

namespace Burkardt.MatrixNS;

public static class Diagonal
{
    public static int[] diag_index ( int m, int[] ia, int[] ja, int n )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIAG_INDEX determines where the diagonal matrix entries are stored.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 January 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of adjacencies.
        //
        //    Input, int IA[M], JA[M], the row and column indices of adjacencies.
        //
        //    Input, int N, the number of nodes.
        //
        //    Output, int DIAG_INDEX[N], contains for each index 0 <= I < N, the unique
        //    index J such that IA[J] = JA[J] = I.
        //
    {
        int[] diag;
        int i;
        int j;

        diag = new int[n];

        for ( i = 0; i < n; i++ )
        {
            diag[i] = -1;
        }

        for ( j = 0; j < m; j++ )
        {
            if ( ia[j] == ja[j] )
            {
                i = ia[j] - 1;
                if ( diag[i] != -1 )
                {
                    Console.WriteLine("");
                    Console.WriteLine("DIAG_INDEX - Fatal error!");
                    Console.WriteLine("  Multiple occurrences of diagonal pairs.");
                    Console.WriteLine("  IA[" + j       + "] = JA[" + j       + "] = " + ia[j] + " and");
                    Console.WriteLine("  IA[" + diag[i] + "] = JA[" + diag[i] + "] = " + ia[j] + "");
                    return null;
                }
                diag[i] = j;
            }
        }

        for ( i = 0; i < n; i++ )
        {
            switch (diag[i])
            {
                case -1:
                    Console.WriteLine("");
                    Console.WriteLine("DIAG_INDEX - Fatal error!");
                    Console.WriteLine("  DIAG[" + i + "] = -1.");
                    return null;
            }
        }

        return diag;
    }        
}