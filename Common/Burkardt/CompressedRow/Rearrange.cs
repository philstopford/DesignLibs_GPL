namespace Burkardt.CompressedRow
{
    public static class RearrangeCR
    {
        public static void rearrange_cr(int n, int nz_num, int[] ia, ref int[] ja, ref double[] a )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    REARRANGE_CR sorts a sparse compressed row matrix.
        //
        //  Discussion:
        //
        //    This routine guarantees that the entries in the CR matrix
        //    are properly sorted.
        //
        //    After the sorting, the entries of the matrix are rearranged in such
        //    a way that the entries of each column are listed in ascending order
        //    of their column values.
        //
        //    The matrix A is assumed to be stored in compressed row format.  Only
        //    the nonzero entries of A are stored.  The vector JA stores the
        //    column index of the nonzero value.  The nonzero values are sorted
        //    by row, and the compressed row vector IA then has the property that
        //    the entries in A and JA that correspond to row I occur in indices
        //    IA[I] through IA[I+1]-1.
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
        //    Input, int IA[N+1], the compressed row index.
        //
        //    Input/output, int JA[NZ_NUM], the column indices.  On output,
        //    the order of the entries of JA may have changed because of the sorting.
        //
        //    Input/output, double A[NZ_NUM], the matrix values.  On output, the
        //    order of the entries may have changed because of the sorting.
        //
        {
            double dtemp;
            int i;
            int is_;
            int itemp;
            int j;
            int j1;
            int j2;
            int k;

            for (i = 0; i < n; i++)
            {
                j1 = ia[i];
                j2 = ia[i + 1];
                    is_ = j2 - j1;

                for (k = 1; k < is_; k++)
                {
                    for (j = j1; j < j2 - k; j++)
                    {
                        if (ja[j + 1] < ja[j])
                        {
                            itemp = ja[j + 1];
                            ja[j + 1] = ja[j];
                            ja[j] = itemp;

                            dtemp = a[j + 1];
                            a[j + 1] = a[j];
                            a[j] = dtemp;
                        }
                    }
                }
            }
        }
    }
}