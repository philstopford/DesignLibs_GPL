namespace Burkardt.CompressedRow
{
    public static class LUSCR
    {
        public static void lus_cr(int n, int nz_num, int[] ia, int[] ja, double[] l, int[] ua,
        double[] r, ref double[] z, int rIndex = 0, int zIndex = 0 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LUS_CR applies the incomplete LU preconditioner.
        //
        //  Discussion:
        //
        //    The linear system M * Z = R is solved for Z.  M is the incomplete
        //    LU preconditioner matrix, and R is a vector supplied by the user.
        //    So essentially, we're solving L * U * Z = R.
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
        //    Input, int IA[N+1], JA[NZ_NUM], the row and column indices
        //    of the matrix values.  The row vector has been compressed.
        //
        //    Input, double L[NZ_NUM], the matrix values.
        //
        //    Input, int UA[N], the index of the diagonal element of each row.
        //
        //    Input, double R[N], the right hand side.
        //
        //    Output, double Z[N], the solution of the system M * Z = R.
        //
        {
            int i;
            int j;
            double[] w;

            w = new double[n];
            //
            //  Copy R in.
            //
            for (i = 0; i < n; i++)
            {
                w[i] = r[rIndex + i];
            }

            //
            //  Solve L * w = w where L is unit lower triangular.
            //
            for (i = 1; i < n; i++)
            {
                for (j = ia[i]; j < ua[i]; j++)
                {
                    w[i] = w[i] - l[j] * w[ja[j]];
                }
            }

            //
            //  Solve U * w = w, where U is upper triangular.
            //
            for (i = n - 1; 0 <= i; i--)
            {
                for (j = ua[i] + 1; j < ia[i + 1]; j++)
                {
                    w[i] = w[i] - l[j] * w[ja[j]];
                }

                w[i] = w[i] / l[ua[i]];
            }

            //
            //  Copy Z out.
            //
            for (i = 0; i < n; i++)
            {
                z[zIndex + i] = w[i];
            }
        }
    }
}