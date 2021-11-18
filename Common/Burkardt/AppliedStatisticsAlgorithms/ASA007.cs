namespace Burkardt.AppliedStatistics;

public static partial class Algorithms
{
    public static void syminv ( double[] a, int n, ref double[] c, double[] w, ref int nullty, 
            ref int ifault )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SYMINV computes the inverse of a symmetric matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 February 2008
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Michael Healy.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Michael Healy,
        //    Algorithm AS 7:
        //    Inversion of a Positive Semi-Definite Symmetric Matrix,
        //    Applied Statistics,
        //    Volume 17, Number 2, 1968, pages 198-199.
        //
        //  Parameters:
        //
        //    Input, double A((N*(N+1))/2), a positive definite matrix stored 
        //    by rows in lower triangular form as a one dimensional array, in the sequence
        //    A(1,1),
        //    A(2,1), A(2,2),
        //    A(3,1), A(3,2), A(3,3), and so on.
        //
        //    Input, int N, the order of A.
        //
        //    Output, double C((N*(N+1))/2), the inverse of A, or generalized
        //    inverse if A is singular, stored using the same storage scheme employed
        //    for A.  The program is written in such a way that A and U can share storage.
        //
        //    Workspace, double W(N).
        //
        //    Output, int *NULLTY, the rank deficiency of A.  If NULLTY is zero,
        //    the matrix is judged to have full rank.
        //
        //    Output, int *IFAULT, error indicator.
        //    0, no error detected.
        //    1, N < 1.
        //    2, A is not positive semi-definite.
        //
    {
        ifault = 0;

        switch (n)
        {
            case <= 0:
                ifault = 1;
                return;
        }

        int nrow = n;
        //
        //  Compute the Cholesky factorization of A.
        //  The result is stored in C.
        //
        int nn = n * (n + 1) / 2;

        cholesky(a, n, nn, ref c, ref nullty, ref ifault);

        if (ifault != 0)
        {
            return;
        }

        //
        //  Invert C and form the product (Cinv)' * Cinv, where Cinv is the inverse
        //  of C, row by row starting with the last row.
        //  IROW = the row number, 
        //  NDIAG = location of last element in the row.
        //
        int irow = nrow;
        int ndiag = nn;
        //
        //  Special case, zero diagonal element.
        //
        for (;;)
        {
            int l;
            switch (c[ndiag - 1])
            {
                case 0.0:
                {
                    l = ndiag;
                    for (int j = irow; j <= nrow; j++)
                    {
                        c[l - 1] = 0.0;
                        l += j;
                    }

                    break;
                }
                default:
                {
                    l = ndiag;
                    for (int i = irow; i <= nrow; i++)
                    {
                        w[i - 1] = c[l - 1];
                        l += i;
                    }

                    int icol = nrow;
                    int jcol = nn;
                    int mdiag = nn;

                    for (;;)
                    {
                        l = jcol;

                        double x;
                        if (icol == irow)
                        {
                            x = 1.0 / w[irow - 1];
                        }
                        else
                        {
                            x = 0.0;
                        }

                        int k = nrow;

                        while (irow < k)
                        {
                            x -= w[k - 1] * c[l - 1];
                            k -= 1;
                            l -= 1;

                            if (mdiag < l)
                            {
                                l = l - k + 1;
                            }
                        }

                        c[l - 1] = x / w[irow - 1];

                        if (icol <= irow)
                        {
                            break;
                        }

                        mdiag -= icol;
                        icol -= 1;
                        jcol -= 1;
                    }

                    break;
                }
            }

            ndiag -= irow;
            irow -= 1;

            if (irow <= 0)
            {
                break;
            }
        }
    }
}