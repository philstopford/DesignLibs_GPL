namespace Burkardt.CompressedRow;

public static class DiagonalPointerCR
{
    public static void diagonal_pointer_cr ( int n, int nz_num, int[] ia, int[] ja, ref int[] ua )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIAGONAL_POINTER_CR finds diagonal entries in a sparse compressed row matrix.
        //
        //  Discussion:
        //
        //    The matrix A is assumed to be stored in compressed row format.  Only
        //    the nonzero entries of A are stored.  The vector JA stores the
        //    column index of the nonzero value.  The nonzero values are sorted
        //    by row, and the compressed row vector IA then has the property that
        //    the entries in A and JA that correspond to row I occur in indices
        //    IA[I] through IA[I+1]-1.
        //
        //    The array UA can be used to locate the diagonal elements of the matrix.
        //
        //    It is assumed that every row of the matrix includes a diagonal element,
        //    and that the elements of each row have been ascending sorted.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 July 2007
        //
        //  Author:
        //
        //    Original C version by Lili Ju.
        //    C++ version by John Burkardt.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the system.
        //
        //    Input, int NZ_NUM, the number of nonzeros.
        //
        //    Input, int IA[N+1], JA[NZ_NUM], the row and column indices
        //    of the matrix values.  The row vector has been compressed.  On output,
        //    the order of the entries of JA may have changed because of the sorting.
        //
        //    Output, int UA[N], the index of the diagonal element of each row.
        //
    {
        int i;
        int j;
        int j1;
        int j2;

        for ( i = 0; i < n; i++ )
        {
            ua[i] = -1;
            j1 = ia[i];
            j2 = ia[i+1];

            for ( j = j1; j < j2; j++ )
            {
                if ( ja[j] == i ) 
                {
                    ua[i] = j;
                }
            }

        }
    }
}