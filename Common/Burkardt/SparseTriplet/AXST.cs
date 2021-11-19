namespace Burkardt.SparseTripletNS;

public static class AXST
{
    public static void ax_st(int n, int nz_num, int[] ia, int[] ja, double[] a, double[] x,
            ref double[] w, int aIndex = 0, int xIndex = 0, int wIndex = 0 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    AX_ST computes A*x for a matrix stored in sparse triplet form.
        //
        //  Discussion:
        //
        //    The matrix A is assumed to be sparse.  To save on storage, only
        //    the nonzero entries of A are stored.  For instance, the K-th nonzero
        //    entry in the matrix is stored by:
        //
        //      A(K) = value of entry,
        //      IA(K) = row of entry,
        //      JA(K) = column of entry.
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
        //    Input, int IA[NZ_NUM], JA[NZ_NUM], the row and column indices
        //    of the matrix values.
        //
        //    Input, double A[NZ_NUM], the matrix values.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Output, double W[N], the value of A*X.
        //
    {
        int i;
        int k;

        for (i = 0; i < n; i++)
        {
            w[wIndex + i] = 0.0;
        }

        for (k = 0; k < nz_num; k++)
        {
            i = ia[k];
            int j = ja[k];
            w[wIndex + i] += a[aIndex + k] * x[xIndex + j];
        }
    }

}