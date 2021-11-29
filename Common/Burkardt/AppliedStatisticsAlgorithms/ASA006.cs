using System;

namespace Burkardt.AppliedStatistics;

public static partial class Algorithms
{
    public static void cholesky ( double[] a, int n, int nn, ref double[] u, ref int nullty, 
            ref int ifault )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CHOLESKY computes the Cholesky factorization of a PDS matrix.
        //
        //  Discussion:
        //
        //    For a positive definite symmetric matrix A, the Cholesky factor U
        //    is an upper triangular matrix such that A = U' * U.
        //
        //    This routine was originally named "CHOL", but that conflicted with
        //    a built in MATLAB routine name.
        //
        //    The missing initialization "II = 0" has been added to the code.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 February 2008
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Michael Healy.
        //    Modifications by AJ Miller.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Michael Healy,
        //    Algorithm AS 6:
        //    Triangular decomposition of a symmetric matrix,
        //    Applied Statistics,
        //    Volume 17, Number 2, 1968, pages 195-197.
        //
        //  Parameters:
        //
        //    Input, double A((N*(N+1))/2), a positive definite matrix 
        //    stored by rows in lower triangular form as a one dimensional array, 
        //    in the sequence
        //    A(1,1),
        //    A(2,1), A(2,2),
        //    A(3,1), A(3,2), A(3,3), and so on.
        //
        //    Input, int N, the order of A.
        //
        //    Input, int NN, the dimension of the array used to store A, 
        //    which should be at least (N*(N+1))/2.
        //
        //    Output, double U((N*(N+1))/2), an upper triangular matrix,
        //    stored by columns, which is the Cholesky factor of A.  The program is
        //    written in such a way that A and U can share storage.
        //
        //    Output, int NULLTY, the rank deficiency of A.  If NULLTY is zero,
        //    the matrix is judged to have full rank.
        //
        //    Output, int IFAULT, an error indicator.
        //    0, no error was detected;
        //    1, if N < 1;
        //    2, if A is not positive semi-definite.
        //    3, if NN < (N*(N+1))/2.
        //
        //  Local Parameters:
        //
        //    Local, double ETA, should be set equal to the smallest positive
        //    value such that 1.0 + ETA is calculated as being greater than 1.0 in the
        //    accuracy being used.
        //
    {
        const double eta = 1.0E-09;
        double w = 0;

        ifault = 0;
        nullty = 0;

        switch (n)
        {
            case <= 0:
                ifault = 1;
                return;
        }

        if ( nn < n * ( n + 1 ) / 2 )
        {
            ifault = 3;
            return;
        }

        int j = 1;
        int k = 0;
        int ii = 0;
        //
        //  Factorize column by column, ICOL = column number.
        //
        for (int icol = 1; icol <= n; icol++ )
        {
            ii += icol;
            double x = eta * eta * a[ii-1];
            int l = 0;
            //
            //  IROW = row number within column ICOL.
            //
            for (int irow = 1; irow <= icol; irow++ )
            {
                k += 1;
                w = a[k-1];
                int m = j;

                int i;
                for ( i = 1; i < irow; i++ )
                {
                    l += 1;
                    w -= u[l-1] * u[m-1];
                    m += 1;
                }

                l += 1;

                if ( irow == icol ) 
                {
                    break;
                }

                if ( u[l-1] != 0.0 )
                {
                    u[k-1] = w / u[l-1];
                }
                else
                {
                    u[k-1] = 0.0;

                    if (!(Math.Abs(x * a[k - 1]) < w * w))
                    {
                        continue;
                    }

                    ifault = 2;
                    return;
                }
            }
            //
            //  End of row, estimate relative accuracy of diagonal element.
            //
            if ( Math.Abs ( w ) <= Math.Abs ( eta * a[k-1] ) )
            {
                u[k-1] = 0.0;
                nullty += 1;
            }
            else
            {
                switch (w)
                {
                    case < 0.0:
                        ifault = 2;
                        return;
                    default:
                        u[k-1] = Math.Sqrt ( w );
                        break;
                }
            }
            j += icol;
        }
    }

    public static void subchl(double[] a, int[] b, int n, ref double[] u, ref int nullty,
            ref int ifault, int ndim, ref double det)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SUBCHL computes the Cholesky factorization of a (subset of a) PDS matrix.
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
        //    Original FORTRAN77 version by Michael Healy, PR Freeman.
        //    C++ version by  John Burkardt.
        //
        //  Reference:
        //
        //    PR Freeman,
        //    Remark AS R44:
        //    A Remark on AS 6 and AS7: Triangular decomposition of a symmetric matrix
        //    and Inversion of a positive semi-definite symmetric matrix,
        //    Applied Statistics,
        //    Volume 31, Number 3, 1982, pages 336-339.
        //
        //    Michael Healy,
        //    Algorithm AS 6:
        //    Triangular decomposition of a symmetric matrix,
        //    Applied Statistics,
        //    Volume 17, Number 2, 1968, pages 195-197.
        //
        //  Parameters:
        //
        //    Input, double A((M*(M+1))/2), a positive definite matrix
        //    stored by rows in lower triangular form as a one dimensional array,
        //    in the sequence
        //    A(1,1),
        //    A(2,1), A(2,2),
        //    A(3,1), A(3,2), A(3,3), and so on.
        //    In the simplest case, M, the order of A, is equal to N.
        //
        //    Input, int B(N), indicates the order in which the
        //    rows and columns of A are to be used.  In the simplest case,
        //    B = (1,2,3...,N).
        //
        //    Input, int N, the order of the matrix, that is,
        //    the matrix formed by using B to select N rows and columns of A.
        //
        //    Output, double U((N*(N+1))/2), an upper triangular matrix,
        //    stored by columns, which is the Cholesky factor of A.  The program is
        //    written in such a way that A and U can share storage.
        //
        //    Output, int *NULLTY, the rank deficiency of A.
        //    If NULLTY is zero, the matrix is judged to have full rank.
        //
        //    Output, int *IFAULT, an error indicator.
        //    0, no error was detected;
        //    1, if N < 1;
        //    2, if A is not positive semi-definite.
        //
        //    Input, int NDIM, the dimension of A and U, which might
        //    be presumed to be (N*(N+1))/2.
        //
        //    Output, double *DET, the determinant of the matrix.
        //
    {
        const double eta = 1.0E-09;

        ifault = 0;
        nullty = 0;
        det = 1.0;

        switch (n)
        {
            case <= 0:
                ifault = 1;
                return;
        }

        ifault = 2;
        int j = 1;
        int k = 0;

        for (int icol = 1; icol <= n; icol++)
        {
            int ij = b[icol - 1] * (b[icol - 1] - 1) / 2;
            int ii = ij + b[icol - 1];
            double x = eta * eta * a[ii - 1];
            int l = 0;

            double w = 0;
            int kk = 0;
            for (int irow = 1; irow <= icol; irow++)
            {
                kk = b[irow - 1] * (b[irow - 1] + 1) / 2;
                k += 1;
                int jj = ij + b[irow - 1];
                w = a[jj - 1];
                int m = j;

                for (int i = 1; i <= irow - 1; i++)
                {
                    l += 1;
                    w -= u[l - 1] * u[m - 1];
                    m += 1;
                }

                l += 1;

                if (irow == icol)
                {
                    break;
                }

                if (u[l - 1] != 0.0)
                {
                    u[k - 1] = w / u[l - 1];
                }
                else
                {
                    if (Math.Abs(x * a[kk - 1]) < w * w)
                    {
                        ifault = 2;
                        return;
                    }

                    u[k - 1] = 0.0;
                }
            }

            if (Math.Abs(eta * a[kk - 1]) <= Math.Abs(w))
            {
                switch (w)
                {
                    case < 0.0:
                        ifault = 2;
                        return;
                    default:
                        u[k - 1] = Math.Sqrt(w);
                        break;
                }
            }
            else
            {
                u[k - 1] = 0.0;
                nullty += 1;
            }

            j += icol;
            det = det * u[k - 1] * u[k - 1];
        }
    }
}