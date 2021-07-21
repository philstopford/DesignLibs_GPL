namespace Burkardt.CompressedRow
{
    public static class ATXCR
    {
        public static void atx_cr(int n, int nz_num, int[] ia, int[] ja, double[] a, double[] x,
                ref double[] w, int aIndex = 0, int xIndex = 0, int wIndex = 0)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ATX_CR computes A'*x for a matrix stored in sparse compressed row form.
            //
            //  Discussion:
            //
            //    The Sparse Compressed Row storage format is used.
            //
            //    The matrix A is assumed to be sparse.  To save on storage, only
            //    the nonzero entries of A are stored.  The vector JA stores the
            //    column index of the nonzero value.  The nonzero values are sorted
            //    by row, and the compressed row vector IA then has the property that
            //    the entries in A and JA that correspond to row I occur in indices
            //    IA[I] through IA[I+1]-1.
            //
            //    For this version of MGMRES, the row and column indices are assumed
            //    to use the C/C++ convention, in which indexing begins at 0.
            //
            //    If your index vectors IA and JA are set up so that indexing is based 
            //    at 1, then each use of those vectors should be shifted down by 1.
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
            //  Reference:
            //
            //    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
            //    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
            //    Charles Romine, Henk van der Vorst,
            //    Templates for the Solution of Linear Systems:
            //    Building Blocks for Iterative Methods,
            //    SIAM, 1994,
            //    ISBN: 0898714710,
            //    LC: QA297.8.T45.
            //
            //    Tim Kelley,
            //    Iterative Methods for Linear and Nonlinear Equations,
            //    SIAM, 2004,
            //    ISBN: 0898713528,
            //    LC: QA297.8.K45.
            //
            //    Yousef Saad,
            //    Iterative Methods for Sparse Linear Systems,
            //    Second Edition,
            //    SIAM, 2003,
            //    ISBN: 0898715342,
            //    LC: QA188.S17.
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
            //    Input, double A[NZ_NUM], the matrix values.
            //
            //    Input, double X[N], the vector to be multiplied by A'.
            //
            //    Output, double W[N], the value of A'*X.
            //
        {
            int i;
            int k;
            int k1;
            int k2;

            for (i = 0; i < n; i++)
            {
                w[wIndex + i] = 0.0;
                k1 = ia[i];
                k2 = ia[i + 1];
                for (k = k1; k < k2; k++)
                {
                    w[wIndex + ja[k]] = w[wIndex + ja[k]] + a[aIndex + k] * x[xIndex + i];
                }
            }
        }

    }
}