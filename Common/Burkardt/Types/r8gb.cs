using System;
using Burkardt.Uniform;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static double r8gb_det(int n, int ml, int mu, double[] a_lu, int[] pivot)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GB_DET computes the determinant of a matrix factored by R8GB_FA or R8GB_TRF.
        //
        //  Discussion:
        //
        //    The R8GB storage format is used for an M by N banded matrix, with lower
        //    bandwidth ML and upper bandwidth MU.  Storage includes room for ML extra
        //    superdiagonals, which may be required to store nonzero entries generated 
        //    during Gaussian elimination.
        //
        //    The original M by N matrix is "collapsed" downward, so that diagonals
        //    become rows of the storage array, while columns are preserved.  The
        //    collapsed array is logically 2*ML+MU+1 by N.  
        //
        //    The two dimensional array can be further reduced to a one dimensional
        //    array, stored by columns.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be positive.
        //
        //    Input, int ML, MU, the lower and upper bandwidths.
        //    ML and MU must be nonnegative, and no greater than N-1.
        //
        //    Input, double A_LU[(2*ML+MU+1)*N], the LU factors from R8GB_FA or R8GB_TRF.
        //
        //    Input, int PIVOT[N], the pivot vector, as computed by R8GB_FA
        //    or R8GB_TRF.
        //
        //    Output, double R8GB_DET, the determinant of the matrix.
        //
    {
        int col = 2 * ml + mu + 1;
        int i;

        double det = 1.0;

        for (i = 0; i < n; i++)
        {
            det *= a_lu[ml + mu + i * col];
        }

        for (i = 0; i < n; i++)
        {
            if (pivot[i] != i + 1)
            {
                det = -det;
            }
        }

        return det;
    }

    public static double[] r8gb_dif2(int m, int n, int ml, int mu)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GB_DIF2 sets up an R8GB second difference matrix.
        //
        //  Discussion:
        //
        //    The R8GB storage format is used for an M by N banded matrix, with lower bandwidth ML
        //    and upper bandwidth MU.  Storage includes room for ML extra superdiagonals, 
        //    which may be required to store nonzero entries generated during Gaussian 
        //    elimination.
        //
        //    The original M by N matrix is "collapsed" downward, so that diagonals
        //    become rows of the storage array, while columns are preserved.  The
        //    collapsed array is logically 2*ML+MU+1 by N.  
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 June 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows of the matrix.
        //    M must be positive.
        //
        //    Input, int N, the number of columns of the matrix.
        //    N must be positive.
        //
        //    Input, int ML, MU, the lower and upper bandwidths.
        //    ML and MU must be nonnegative, and no greater than min(M,N)-1.
        //
        //    Output, double R8GB_DIF2[(2*ML+MU+1)*N], the R8GB matrix.
        //
    {
        int col = 2 * ml + mu + 1;
        int j;

        double[] a = r8vec_zeros_new((2 * ml + mu + 1) * n);

        for (j = 1; j <= n; j++)
        {
            int diag;
            for (diag = 1; diag <= 2 * ml + mu + 1; diag++)
            {
                int i = diag + j - ml - mu - 1;

                if (i == j)
                {
                    a[diag - 1 + (j - 1) * col] = 2.0;
                }
                else if (i == j - 1 || i == j + 1)
                {
                    a[diag - 1 + (j - 1) * col] = -1.0;
                }
            }
        }

        return a;
    }

    public static int r8gb_fa(int n, int ml, int mu, ref double[] a, ref int[] pivot)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GB_FA performs a LINPACK-style PLU factorization of an R8GB matrix.
        //
        //  Discussion:
        //
        //    The R8GB storage format is used for an M by N banded matrix, with lower
        //    bandwidth ML and upper bandwidth MU.  Storage includes room for ML extra
        //    superdiagonals, which may be required to store nonzero entries generated 
        //    during Gaussian elimination.
        //
        //    The original M by N matrix is "collapsed" downward, so that diagonals
        //    become rows of the storage array, while columns are preserved.  The
        //    collapsed array is logically 2*ML+MU+1 by N.  
        //
        //    The two dimensional array can be further reduced to a one dimensional
        //    array, stored by columns.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 September 2003
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
        //    LINPACK User's Guide,
        //    SIAM, 1979,
        //    ISBN13: 978-0-898711-72-1,
        //    LC: QA214.L56.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be positive.
        //
        //    Input, int ML, MU, the lower and upper bandwidths.
        //    ML and MU must be nonnegative, and no greater than N-1.
        //
        //    Input/output, double A[(2*ML+MU+1)*N], the matrix in band storage.  
        //    On output, A has been overwritten by the LU factors.
        //
        //    Output, int PIVOT[N], the pivot vector.
        //
        //    Output, int R8GB_FA, singularity flag.
        //    0, no singularity detected.
        //    nonzero, the factorization failed on the INFO-th step.
        //
    {
        int col = 2 * ml + mu + 1;
        int i;
        int jz;
        int k;

        int m = ml + mu + 1;
        //
        //  Zero out the initial fill-in columns.
        //
        int j0 = mu + 2;
        int j1 = Math.Min(n, m) - 1;

        for (jz = j0; jz <= j1; jz++)
        {
            int i0 = m + 1 - jz;
            for (i = i0; i <= ml; i++)
            {
                a[i - 1 + (jz - 1) * col] = 0.0;
            }
        }

        jz = j1;
        int ju = 0;

        for (k = 1; k <= n - 1; k++)
        {
            //
            //  Zero out the next fill-in column.
            //
            jz += 1;
            if (jz <= n)
            {
                for (i = 1; i <= ml; i++)
                {
                    a[i - 1 + (jz - 1) * col] = 0.0;
                }
            }

            //
            //  Find L = pivot index.
            //
            int lm = Math.Min(ml, n - k);
            int l = m;

            int j;
            for (j = m + 1; j <= m + lm; j++)
            {
                if (Math.Abs(a[l - 1 + (k - 1) * col]) < Math.Abs(a[j - 1 + (k - 1) * col]))
                {
                    l = j;
                }
            }

            pivot[k - 1] = l + k - m;
            switch (a[l - 1 + (k - 1) * col])
            {
                //
                //  Zero pivot implies this column already triangularized.
                //
                case 0.0:
                    Console.WriteLine("");
                    Console.WriteLine("R8GB_FA - Fatal error!");
                    Console.WriteLine("  Zero pivot on step " + k + "");
                    return 1;
            }

            //
            //  Interchange if necessary.
            //
            double t = a[l - 1 + (k - 1) * col];
            a[l - 1 + (k - 1) * col] = a[m - 1 + (k - 1) * col];
            a[m - 1 + (k - 1) * col] = t;
            //
            //  Compute multipliers.
            //
            for (i = m + 1; i <= m + lm; i++)
            {
                a[i - 1 + (k - 1) * col] = -a[i - 1 + (k - 1) * col] / a[m - 1 + (k - 1) * col];
            }

            //
            //  Row elimination with column indexing.
            //
            ju = Math.Max(ju, mu + pivot[k - 1]);
            ju = Math.Min(ju, n);
            int mm = m;

            for (j = k + 1; j <= ju; j++)
            {
                l -= 1;
                mm -= 1;

                if (l != mm)
                {
                    t = a[l - 1 + (j - 1) * col];
                    a[l - 1 + (j - 1) * col] = a[mm - 1 + (j - 1) * col];
                    a[mm - 1 + (j - 1) * col] = t;
                }

                for (i = 1; i <= lm; i++)
                {
                    a[mm + i - 1 + (j - 1) * col] += a[mm - 1 + (j - 1) * col] * a[m + i - 1 + (k - 1) * col];
                }
            }
        }

        pivot[n - 1] = n;

        switch (a[m - 1 + (n - 1) * col])
        {
            case 0.0:
                Console.WriteLine("");
                Console.WriteLine("R8GB_FA - Fatal error!");
                Console.WriteLine("  Zero pivot on step " + n + "");
                return 1;
            default:
                return 0;
        }
    }

    public static double[] r8gb_indicator(int m, int n, int ml, int mu)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GB_INDICATOR sets up an R8GB indicator matrix.
        //
        //  Discussion:
        //
        //    Note that the R8GB storage format includes extra room for
        //    fillin entries that occur during Gauss elimination.  This routine
        //    will leave those values at 0.
        //
        //    The R8GB storage format is used for an M by N banded matrix, with lower
        //    bandwidth ML and upper bandwidth MU.  Storage includes room for ML extra
        //    superdiagonals, which may be required to store nonzero entries generated 
        //    during Gaussian elimination.
        //
        //    The original M by N matrix is "collapsed" downward, so that diagonals
        //    become rows of the storage array, while columns are preserved.  The
        //    collapsed array is logically 2*ML+MU+1 by N.  
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 March 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows of the matrix.
        //    M must be positive.
        //
        //    Input, int N, the number of columns of the matrix.
        //    N must be positive.
        //
        //    Input, int ML, MU, the lower and upper bandwidths.
        //    ML and MU must be nonnegative, and no greater than min(M,N)-1.
        //
        //    Output, double R8GB_INDICATOR[(2*ML+MU+1)*N], the R8GB matrix.
        //
    {
        int col = 2 * ml + mu + 1;
        int j;

        double[] a = new double[(2 * ml + mu + 1) * n];

        int fac = (int) Math.Pow(10, (int) Math.Log10(n) + 1);
        int k = 0;

        for (j = 1; j <= n; j++)
        {
            int diag;
            for (diag = 1; diag <= 2 * ml + mu + 1; diag++)
            {
                int i = diag + j - ml - mu - 1;

                switch (i)
                {
                    case >= 1 when i <= m && i - ml <= j && j <= i + mu:
                        a[diag - 1 + (j - 1) * col] = fac * i + j;
                        break;
                    case >= 1 when i <= m && i - ml <= j && j <= i + mu + ml:
                        a[diag - 1 + (j - 1) * col] = 0.0;
                        break;
                    default:
                        k += 1;
                        a[diag - 1 + (j - 1) * col] = -(double) k;
                        break;
                }
            }
        }

        return a;
    }

    public static double[] r8gb_ml(int n, int ml, int mu, double[] a_lu, int[] pivot, double[] x,
            int job)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GB_ML computes A * x or A' * X, using R8GB_FA factors.
        //
        //  Discussion:
        //
        //    The R8GB storage format is used for an M by N banded matrix, with lower
        //    bandwidth ML and upper bandwidth MU.  Storage includes room for ML extra
        //    superdiagonals, which may be required to store nonzero entries generated 
        //    during Gaussian elimination.
        //
        //    The original M by N matrix is "collapsed" downward, so that diagonals
        //    become rows of the storage array, while columns are preserved.  The
        //    collapsed array is logically 2*ML+MU+1 by N.  
        //
        //    The two dimensional array can be further reduced to a one dimensional
        //    array, stored by columns.
        //
        //    It is assumed that R8GB_FA has overwritten the original matrix
        //    information by LU factors.  R8GB_ML is able to reconstruct the
        //    original matrix from the LU factor data.
        //
        //    R8GB_ML allows the user to check that the solution of a linear
        //    system is correct, without having to save an unfactored copy
        //    of the matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be positive.
        //
        //    Input, int ML, MU, the lower and upper bandwidths.
        //    ML and MU must be nonnegative, and no greater than N-1.
        //
        //    Input, double A_LU[(2*ML+MU+1)*N], the LU factors from R8GB_FA.
        //
        //    Input, int PIVOT[N], the pivot vector computed by R8GB_FA.
        //
        //    Input, double X[N], the vector to be multiplied.
        //
        //    Input, int JOB, specifies the operation to be done:
        //    JOB = 0, compute A * x.
        //    JOB nonzero, compute A' * X.
        //
        //    Output, double R8GB_ML[N], the result of the multiplication.
        //
    {
        int col = 2 * ml + mu + 1;
        int i;
        int j;
        int k;
        double temp;

        double[] b = new double[n];

        for (i = 0; i < n; i++)
        {
            b[i] = x[i];
        }

        switch (job)
        {
            case 0:
            {
                //
                //  Y = U * X.
                //
                for (j = 1; j <= n; j++)
                {
                    int ilo = Math.Max(1, j - ml - mu);
                    for (i = ilo; i <= j - 1; i++)
                    {
                        b[i - 1] += a_lu[i - j + ml + mu + (j - 1) * col] * b[j - 1];
                    }

                    b[j - 1] = a_lu[j - j + ml + mu + (j - 1) * col] * b[j - 1];
                }

                //
                //  B = PL * Y = PL * U * X = A * x.
                //
                for (j = n - 1; 1 <= j; j--)
                {
                    int ihi = Math.Min(n, j + ml);
                    for (i = j + 1; i <= ihi; i++)
                    {
                        b[i - 1] -= a_lu[i - j + ml + mu + (j - 1) * col] * b[j - 1];
                    }

                    k = pivot[j - 1];

                    if (k != j)
                    {
                        temp = b[k - 1];
                        b[k - 1] = b[j - 1];
                        b[j - 1] = temp;
                    }
                }

                break;
            }
            default:
            {
                //
                //  Y = ( PL )' * X.
                //
                int jhi;
                for (j = 1; j <= n - 1; j++)
                {
                    k = pivot[j - 1];

                    if (k != j)
                    {
                        temp = b[k - 1];
                        b[k - 1] = b[j - 1];
                        b[j - 1] = temp;
                    }

                    jhi = Math.Min(n, j + ml);
                    for (i = j + 1; i <= jhi; i++)
                    {
                        b[j - 1] -= b[i - 1] * a_lu[i - j + ml + mu + (j - 1) * col];
                    }

                }

                //
                //  B = U' * Y = ( PL * U )' * X = A' * X.
                //
                for (i = n; 1 <= i; i--)
                {
                    jhi = Math.Min(n, i + ml + mu);
                    for (j = i + 1; j <= jhi; j++)
                    {
                        b[j - 1] += b[i - 1] * a_lu[i - j + ml + mu + (j - 1) * col];
                    }

                    b[i - 1] *= a_lu[i - i + ml + mu + (i - 1) * col];
                }

                break;
            }
        }

        return b;
    }

    public static double[] r8gb_mu(int n, int ml, int mu, double[] a_lu, int[] pivot, double[] x,
            int job)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GB_MU computes A * x or A' * X, using R8GB_TRF factors.
        //
        //  Warning:
        //
        //    This routine must be updated to allow for rectangular matrices.
        //
        //  Discussion:
        //
        //    The R8GB storage format is used for an M by N banded matrix, with lower
        //    bandwidth ML and upper bandwidth MU.  Storage includes room for ML extra
        //    superdiagonals, which may be required to store nonzero entries generated 
        //    during Gaussian elimination.
        //
        //    The original M by N matrix is "collapsed" downward, so that diagonals
        //    become rows of the storage array, while columns are preserved.  The
        //    collapsed array is logically 2*ML+MU+1 by N.  
        //
        //    It is assumed that R8GB_TRF has overwritten the original matrix
        //    information by LU factors.  R8GB_MU is able to reconstruct the
        //    original matrix from the LU factor data.
        //
        //    R8GB_MU allows the user to check that the solution of a linear
        //    system is correct, without having to save an unfactored copy
        //    of the matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 October 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Edward Anderson, Zhaojun Bai, Christian Bischof, Susan Blackford, 
        //    James Demmel, Jack Dongarra, Jeremy Du Croz, Anne Greenbaum, 
        //    Sven Hammarling, Alan McKenney, Danny Sorensen,
        //    LAPACK User's Guide,
        //    Second Edition,
        //    SIAM, 1995.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be positive.
        //
        //    Input, int ML, MU, the lower and upper bandwidths.
        //    ML and MU must be nonnegative, and no greater than N-1.
        //
        //    Input, double A_LU[(2*ML+MU+1)*N], the LU factors from R8GB_TRF.
        //
        //    Input, int PIVOT[N], the pivot vector computed by R8GB_TRF.
        //
        //    Input, double X[N], the vector to be multiplied.
        //
        //    Input, int JOB, specifies the operation to be done:
        //    JOB = 0, compute A * x.
        //    JOB nonzero, compute A' * X.
        //
        //    Output, double R8GB_MU[N], the result of the multiplication.
        //
    {
        int i;
        int j;
        int k;
        double t;

        double[] b = new double[n];

        for (i = 0; i < n; i++)
        {
            b[i] = x[i];
        }

        switch (job)
        {
            case 0:
            {
                //
                //  Y = U * X.
                //
                for (j = 1; j <= n; j++)
                {
                    int ilo = Math.Max(1, j - ml - mu);
                    for (i = ilo; i <= j - 1; i++)
                    {
                        b[i - 1] += a_lu[i - j + ml + mu + (j - 1) * (2 * ml + mu + 1)] * b[j - 1];
                    }

                    b[j - 1] = a_lu[j - j + ml + mu + (j - 1) * (2 * ml + mu + 1)] * b[j - 1];
                }

                //
                //  B = PL * Y = PL * U * X = A * x.
                //
                for (j = n - 1; 1 <= j; j--)
                {
                    int ihi = Math.Min(n, j + ml);
                    for (i = j + 1; i <= ihi; i++)
                    {
                        b[i - 1] += a_lu[i - j + ml + mu + (j - 1) * (2 * ml + mu + 1)] * b[j - 1];
                    }

                    k = pivot[j - 1];

                    if (k != j)
                    {
                        t = b[k - 1];
                        b[k - 1] = b[j - 1];
                        b[j - 1] = t;
                    }
                }

                break;
            }
            default:
            {
                //
                //  Y = ( PL )' * X.
                //
                int jhi;
                for (j = 1; j <= n - 1; j++)
                {
                    k = pivot[j - 1];

                    if (k != j)
                    {
                        t = b[k - 1];
                        b[k - 1] = b[j - 1];
                        b[j - 1] = t;
                    }

                    jhi = Math.Min(n, j + ml);
                    for (i = j + 1; i <= jhi; i++)
                    {
                        b[j - 1] += b[i - 1] * a_lu[i - j + ml + mu + (j - 1) * (2 * ml + mu + 1)];
                    }
                }

                //
                //  B = U' * Y = ( PL * U )' * X = A' * X.
                //
                for (i = n; 1 <= i; i--)
                {
                    jhi = Math.Min(n, i + ml + mu);
                    for (j = i + 1; j <= jhi; j++)
                    {
                        b[j - 1] += b[i - 1] * a_lu[i - j + ml + mu + (j - 1) * (2 * ml + mu + 1)];
                    }

                    b[i - 1] *= a_lu[i - i + ml + mu + (i - 1) * (2 * ml + mu + 1)];
                }

                break;
            }
        }

        return b;
    }

    public static double[] r8gb_mtv(int m, int n, int ml, int mu, double[] a, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GB_MTV multilies a vector times an R8GB matrix.
        //
        //  Discussion:
        //
        //    The R8GB storage format is used for an M by N banded matrix, with lower
        //    bandwidth ML and upper bandwidth MU.  Storage includes room for ML extra
        //    superdiagonals, which may be required to store nonzero entries generated 
        //    during Gaussian elimination.
        //
        //    The original M by N matrix is "collapsed" downward, so that diagonals
        //    become rows of the storage array, while columns are preserved.  The
        //    collapsed array is logically 2*ML+MU+1 by N.  
        //
        //    The two dimensional array can be further reduced to a one dimensional
        //    array, stored by columns.
        //
        //    For our purposes, X*A and A'*X mean the same thing.
        //
        //    LINPACK and LAPACK storage of general band matrices requires
        //    an extra ML upper diagonals for possible fill in entries during
        //    Gauss elimination.  This routine does not access any entries
        //    in the fill in diagonals, because it assumes that the matrix
        //    has NOT had Gauss elimination applied to it.  If the matrix
        //    has been Gauss eliminated, then the routine R8GB_MU must be
        //    used instead.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows of the matrix.
        //    M must be positive.
        //
        //    Input, int N, the number of columns of the matrix.
        //    N must be positive.
        //
        //    Input, int ML, MU, the lower and upper bandwidths.
        //    ML and MU must be nonnegative, and no greater than min(M,N)-1.
        //
        //    Input, double A[(2*ML+MU+1)*N], the R8GB matrix.
        //
        //    Input, double X[M], the vector to be multiplied by A.
        //
        //    Output, double R8GB_MTV[N], the product X*A or A'*X.
        //
    {
        int col = 2 * ml + mu + 1;
        int j;

        double[] b = r8vec_zeros_new(n);

        for (j = 1; j <= n; j++)
        {
            int ilo = Math.Max(1, j - mu);
            int ihi = Math.Min(m, j + ml);
            int i;
            for (i = ilo; i <= ihi; i++)
            {
                b[j - 1] += x[i - 1] * a[i - j + ml + mu + (j - 1) * col];
            }
        }

        return b;
    }

    public static double[] r8gb_mv(int m, int n, int ml, int mu, double[] a, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GB_MV multiplies an R8GB matrix times a vector.
        //
        //  Discussion:
        //
        //    The R8GB storage format is used for an M by N banded matrix, with lower
        //    bandwidth ML and upper bandwidth MU.  Storage includes room for ML extra
        //    superdiagonals, which may be required to store nonzero entries generated 
        //    during Gaussian elimination.
        //
        //    The original M by N matrix is "collapsed" downward, so that diagonals
        //    become rows of the storage array, while columns are preserved.  The
        //    collapsed array is logically 2*ML+MU+1 by N.  
        //
        //    The two dimensional array can be further reduced to a one dimensional
        //    array, stored by columns.
        //
        //    LINPACK and LAPACK storage of general band matrices requires
        //    an extra ML upper diagonals for possible fill in entries during
        //    Gauss elimination.  This routine does not access any entries
        //    in the fill in diagonals, because it assumes that the matrix
        //    has NOT had Gauss elimination applied to it.  If the matrix
        //    has been Gauss eliminated, then the routine R8GB_MU must be
        //    used instead.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows of the matrix.
        //    M must be positive.
        //
        //    Input, int N, the number of columns of the matrix.
        //    N must be positive.
        //
        //    Input, int ML, MU, the lower and upper bandwidths.
        //    ML and MU must be nonnegative, and no greater than min(M,N)-1.
        //
        //    Input, double A[(2*ML+MU+1)*N], the R8GB matrix.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Output, double R8GB_MV[M], the product A * x.
        //
    {
        int col = 2 * ml + mu + 1;
        int i;

        double[] b = r8vec_zeros_new(m);

        for (i = 1; i <= m; i++)
        {
            int jlo = Math.Max(1, i - ml);
            int jhi = Math.Min(n, i + mu);
            int j;
            for (j = jlo; j <= jhi; j++)
            {
                b[i - 1] += a[i - j + ml + mu + (j - 1) * col] * x[j - 1];
            }
        }

        return b;
    }

    public static int r8gb_nz_num(int m, int n, int ml, int mu, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GB_NZ_NUM counts the nonzeroes in an R8GB matrix.
        //
        //  Discussion:
        //
        //    The R8GB storage format is used for an M by N banded matrix, with lower
        //    bandwidth ML and upper bandwidth MU.  Storage includes room for ML extra
        //    superdiagonals, which may be required to store nonzero entries generated 
        //    during Gaussian elimination.
        //
        //    The original M by N matrix is "collapsed" downward, so that diagonals
        //    become rows of the storage array, while columns are preserved.  The
        //    collapsed array is logically 2*ML+MU+1 by N.  
        //
        //    LINPACK and LAPACK band storage requires that an extra ML
        //    superdiagonals be supplied to allow for fillin during Gauss
        //    elimination.  Even though a band matrix is described as
        //    having an upper bandwidth of MU, it effectively has an
        //    upper bandwidth of MU+ML.  This routine will examine
        //    values it finds in these extra bands, so that both unfactored
        //    and factored matrices can be handled.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 October 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows of the matrix.
        //    M must be positive.
        //
        //    Input, int N, the number of columns of the matrix.
        //    N must be positive.
        //
        //    Input, int ML, MU, the lower and upper bandwidths.
        //    ML and MU must be nonnegative, and no greater than min(M,N)-1.
        //
        //    Input, double A[(2*ML+MU+1)*N], the R8GB matrix.
        //
        //    Output, int R8GB_NZ_NUM, the number of nonzero entries in A.
        //
    {
        int i;

        int nz_num = 0;

        for (i = 0; i < m; i++)
        {
            int jlo = Math.Max(0, i - ml);
            int jhi = Math.Min(n - 1, i + mu + ml);
            int j;
            for (j = jlo; j <= jhi; j++)
            {
                if (a[i - j + ml + mu + j * (2 * ml + mu + 1)] != 0.0)
                {
                    nz_num += 1;
                }
            }
        }

        return nz_num;
    }

    public static void r8gb_print(int m, int n, int ml, int mu, double[] a, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GB_PRINT prints an R8GB matrix.
        //
        //  Discussion:
        //
        //    The R8GB storage format is used for an M by N banded matrix, with lower
        //    bandwidth ML and upper bandwidth MU.  Storage includes room for ML extra
        //    superdiagonals, which may be required to store nonzero entries generated 
        //    during Gaussian elimination.
        //
        //    The original M by N matrix is "collapsed" downward, so that diagonals
        //    become rows of the storage array, while columns are preserved.  The
        //    collapsed array is logically 2*ML+MU+1 by N.  
        //
        //    The two dimensional array can be further reduced to a one dimensional
        //    array, stored by columns.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 April 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows of the matrix.
        //    M must be positive.
        //
        //    Input, int N, the number of columns of the matrix.
        //    N must be positive.
        //
        //    Input, int ML, MU, the lower and upper bandwidths.
        //    ML and MU must be nonnegative, and no greater than min(M,N)-1..
        //
        //    Input, double A[(2*ML+MU+1)*N], the R8GB matrix.
        //
        //    Input, string TITLE, a title.
        //
    {
        r8gb_print_some(m, n, ml, mu, a, 1, 1, m, n, title);
    }

    public static void r8gb_print_some(int m, int n, int ml, int mu, double[] a, int ilo,
            int jlo, int ihi, int jhi, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GB_PRINT_SOME prints some of an R8GB matrix.
        //
        //  Discussion:
        //
        //    The R8GB storage format is used for an M by N banded matrix, with lower
        //    bandwidth ML and upper bandwidth MU.  Storage includes room for ML extra
        //    superdiagonals, which may be required to store nonzero entries generated 
        //    during Gaussian elimination.
        //
        //    The original M by N matrix is "collapsed" downward, so that diagonals
        //    become rows of the storage array, while columns are preserved.  The
        //    collapsed array is logically 2*ML+MU+1 by N.  
        //
        //    The two dimensional array can be further reduced to a one dimensional
        //    array, stored by columns.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 June 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows of the matrix.
        //    M must be positive.
        //
        //    Input, int N, the number of columns of the matrix.
        //    N must be positive.
        //
        //    Input, int ML, MU, the lower and upper bandwidths.
        //    ML and MU must be nonnegative, and no greater than min(M,N)-1..
        //
        //    Input, double A[(2*ML+MU+1)*N], the R8GB matrix.
        //
        //    Input, int ILO, JLO, IHI, JHI, designate the first row and
        //    column, and the last row and column to be printed.
        //
        //    Input, string TITLE, a title.
        //
    {
        const int INCX = 5;

        int col = 2 * ml + mu + 1;
        int j2lo;

        Console.WriteLine("");
        Console.WriteLine(title + "");
        //
        //  Print the columns of the matrix, in strips of 5.
        //
        for (j2lo = jlo; j2lo <= jhi; j2lo += INCX)
        {
            int j2hi = j2lo + INCX - 1;
            j2hi = Math.Min(j2hi, n);
            j2hi = Math.Min(j2hi, jhi);

            Console.WriteLine("");
            string cout = "  Col: ";
            int j;
            for (j = j2lo; j <= j2hi; j++)
            {
                cout += j.ToString().PadLeft(7) + "       ";
            }

            Console.WriteLine(cout);
            Console.WriteLine("  Row");
            Console.WriteLine("  ---");
            //
            //  Determine the range of the rows in this strip.
            //
            int i2lo = Math.Max(ilo, 1);
            i2lo = Math.Max(i2lo, j2lo - mu);

            int i2hi = Math.Min(ihi, m);
            i2hi = Math.Min(i2hi, j2hi + ml);

            int i;
            for (i = i2lo; i <= i2hi; i++)
            {
                //
                //  Print out (up to) 5 entries in row I, that lie in the current strip.
                //
                cout = i.ToString().PadLeft(6) + "  ";
                for (j = j2lo; j <= j2hi; j++)
                {
                    if (mu < j - i || ml < i - j)
                    {
                        cout += "            ";
                    }
                    else
                    {
                        cout += a[i - j + ml + mu + (j - 1) * col].ToString().PadLeft(10) + "  ";
                    }
                }

                Console.WriteLine(cout);
            }
        }

    }

    public static double[] r8gb_random(int m, int n, int ml, int mu, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GB_RANDOM randomizes an R8GB matrix.
        //
        //  Discussion:
        //
        //    The R8GB storage format is used for an M by N banded matrix, with lower
        //    bandwidth ML and upper bandwidth MU.  Storage includes room for ML extra
        //    superdiagonals, which may be required to store nonzero entries generated 
        //    during Gaussian elimination.
        //
        //    The original M by N matrix is "collapsed" downward, so that diagonals
        //    become rows of the storage array, while columns are preserved.  The
        //    collapsed array is logically 2*ML+MU+1 by N.  
        //
        //    The two dimensional array can be further reduced to a one dimensional
        //    array, stored by columns.
        //
        //    LINPACK and LAPACK band storage requires that an extra ML
        //    superdiagonals be supplied to allow for fillin during Gauss
        //    elimination.  Even though a band matrix is described as
        //    having an upper bandwidth of MU, it effectively has an
        //    upper bandwidth of MU+ML.  This routine assumes it is setting
        //    up an unfactored matrix, so it only uses the first MU upper bands,
        //    and does not place nonzero values in the fillin bands.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows of the matrix.
        //    M must be positive.
        //
        //    Input, int N, the number of columns of the matrix.
        //    N must be positive.
        //
        //    Input, int ML, MU, the lower and upper bandwidths.
        //    ML and MU must be nonnegative, and no greater than min(M,N)-1.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double R8GB_RANDOM[(2*ML+MU+1)*N], the R8GB matrix.
        //
    {
        int col = 2 * ml + mu + 1;
        int i;
        int j;

        double[] a = new double[(2 * ml + mu + 1) * n];

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < col; i++)
            {
                a[i + j * col] = 0.0;
            }
        }

        for (j = 1; j <= n; j++)
        {
            int row;
            for (row = 1; row <= col; row++)
            {
                i = row + j - ml - mu - 1;
                if (ml < row && 1 <= i && i <= m)
                {
                    a[row - 1 + (j - 1) * col] = UniformRNG.r8_uniform_01(ref seed);
                }
                else
                {
                    a[row - 1 + (j - 1) * col] = 0.0;
                }
            }
        }

        return a;
    }

    public static double[] r8gb_sl(int n, int ml, int mu, double[] a_lu, int[] pivot,
            double[] b, int job)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GB_SL solves a system factored by R8GB_FA.
        //
        //  Discussion:
        //
        //    The R8GB storage format is used for an M by N banded matrix, with lower
        //    bandwidth ML and upper bandwidth MU.  Storage includes room for ML extra
        //    superdiagonals, which may be required to store nonzero entries generated 
        //    during Gaussian elimination.
        //
        //    The original M by N matrix is "collapsed" downward, so that diagonals
        //    become rows of the storage array, while columns are preserved.  The
        //    collapsed array is logically 2*ML+MU+1 by N.  
        //
        //    The two dimensional array can be further reduced to a one dimensional
        //    array, stored by columns.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 September 2003
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Dongarra, Bunch, Moler, Stewart.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
        //    LINPACK User's Guide,
        //    SIAM, 1979,
        //    ISBN13: 978-0-898711-72-1,
        //    LC: QA214.L56.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be positive.
        //
        //    Input, int ML, MU, the lower and upper bandwidths.
        //    ML and MU must be nonnegative, and no greater than N-1.
        //
        //    Input, double A_LU[(2*ML+MU+1)*N], the LU factors from R8GB_FA.
        //
        //    Input, int PIVOT[N], the pivot vector from R8GB_FA.
        //
        //    Input, double B[N], the right hand side vector.
        //
        //    Input, int JOB.
        //    0, solve A * x = b.
        //    nonzero, solve A' * x = b.
        //
        //    Output, double R8GB_SL[N], the solution.
        //
    {
        int col = 2 * ml + mu + 1;
        int i;
        int k;
        int l;
        int la;
        int lb;
        int lm;
        double t;

        double[] x = new double[n];

        for (i = 0; i < n; i++)
        {
            x[i] = b[i];
        }

        int m = mu + ml + 1;
        switch (job)
        {
            //
            //  Solve A * x = b.
            //
            case 0:
            {
                switch (ml)
                {
                    //
                    //  Solve L * Y = B.
                    //
                    case >= 1:
                    {
                        for (k = 1; k <= n - 1; k++)
                        {
                            lm = Math.Min(ml, n - k);
                            l = pivot[k - 1];

                            if (l != k)
                            {
                                t = x[l - 1];
                                x[l - 1] = x[k - 1];
                                x[k - 1] = t;
                            }

                            for (i = 1; i <= lm; i++)
                            {
                                x[k + i - 1] += x[k - 1] * a_lu[m + i - 1 + (k - 1) * col];
                            }
                        }

                        break;
                    }
                }

                //
                //  Solve U * X = Y.
                //
                for (k = n; 1 <= k; k--)
                {
                    x[k - 1] /= a_lu[m - 1 + (k - 1) * col];
                    lm = Math.Min(k, m) - 1;
                    la = m - lm;
                    lb = k - lm;
                    for (i = 0; i <= lm - 1; i++)
                    {
                        x[lb + i - 1] -= x[k - 1] * a_lu[la + i - 1 + (k - 1) * col];
                    }
                }

                break;
            }
            //
            default:
            {
                //
                //  Solve U' * Y = B.
                //
                for (k = 1; k <= n; k++)
                {
                    lm = Math.Min(k, m) - 1;
                    la = m - lm;
                    lb = k - lm;
                    for (i = 0; i <= lm - 1; i++)
                    {
                        x[k - 1] -= x[lb + i - 1] * a_lu[la + i - 1 + (k - 1) * col];
                    }

                    x[k - 1] /= a_lu[m - 1 + (k - 1) * col];
                }

                switch (ml)
                {
                    //
                    //  Solve L' * X = Y.
                    //
                    case >= 1:
                    {
                        for (k = n - 1; 1 <= k; k--)
                        {
                            lm = Math.Min(ml, n - k);
                            for (i = 1; i <= lm; i++)
                            {
                                x[k - 1] += x[k + i - 1] * a_lu[m + i - 1 + (k - 1) * col];
                            }

                            l = pivot[k - 1];

                            if (l != k)
                            {
                                t = x[l - 1];
                                x[l - 1] = x[k - 1];
                                x[k - 1] = t;
                            }
                        }

                        break;
                    }
                }

                break;
            }
        }

        return x;
    }

    public static double[] r8gb_to_r8vec(int m, int n, int ml, int mu, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GB_TO_R8VEC copies an R8GB matrix to a real vector.
        //
        //  Discussion:
        //
        //    In C++  and FORTRAN, this routine is not really needed.  In MATLAB,
        //    a data item carries its dimensionality implicitly, and so cannot be
        //    regarded sometimes as a vector and sometimes as an array.
        //
        //    The R8GB storage format is used for an M by N banded matrix, with lower
        //    bandwidth ML and upper bandwidth MU.  Storage includes room for ML extra
        //    superdiagonals, which may be required to store nonzero entries generated 
        //    during Gaussian elimination.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 March 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns in the array.
        //
        //    Input, int ML, MU, the lower and upper bandwidths.
        //
        //    Input, double A[(2*ML+MU+1)*N], the array.
        //
        //    Output, double R8GB_TO_R8VEC[(2*ML+MU+1)*N], the vector.
        //
    {
        int j;

        double[] x = new double[(2 * ml + mu + 1) * n];

        for (j = 1; j <= n; j++)
        {
            int ihi = Math.Min(ml + mu, ml + mu + 1 - j);
            int i;
            for (i = 1; i <= ihi; i++)
            {
                x[i - 1 + (j - 1) * (2 * ml + mu + 1)] = 0.0;
            }

            int ilo = Math.Max(ihi + 1, 1);
            ihi = Math.Min(2 * ml + mu + 1, ml + mu + m + 1 - j);
            for (i = ilo; i <= ihi; i++)
            {
                x[i - 1 + (j - 1) * (2 * ml + mu + 1)] = a[i - 1 + (j - 1) * (2 * ml + mu + 1)];
            }

            ilo = ihi + 1;
            ihi = 2 * ml + mu + 1;
            for (i = ilo; i <= ihi; i++)
            {
                x[i - 1 + (j - 1) * (2 * ml + mu + 1)] = 0.0;
            }

        }

        return x;
    }

    public static double[] r8gb_to_r8ge(int m, int n, int ml, int mu, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GB_TO_R8GE copies an R8GB matrix to an R8GE matrix.
        //
        //  Discussion:
        //
        //    The R8GB storage format is used for an M by N banded matrix, with lower
        //    bandwidth ML and upper bandwidth MU.  Storage includes room for ML extra
        //    superdiagonals, which may be required to store nonzero entries generated 
        //    during Gaussian elimination.
        //
        //    The original M by N matrix is "collapsed" downward, so that diagonals
        //    become rows of the storage array, while columns are preserved.  The
        //    collapsed array is logically 2*ML+MU+1 by N.  
        //
        //    The two dimensional array can be further reduced to a one dimensional
        //    array, stored by columns.
        //
        //    LINPACK and LAPACK band storage requires that an extra ML
        //    superdiagonals be supplied to allow for fillin during Gauss
        //    elimination.  Even though a band matrix is described as
        //    having an upper bandwidth of MU, it effectively has an
        //    upper bandwidth of MU+ML.  This routine will copy nonzero
        //    values it finds in these extra bands, so that both unfactored
        //    and factored matrices can be handled.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 May 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows of the matrices.
        //    M must be positive.
        //
        //    Input, int N, the number of columns of the matrices.
        //    N must be positive.
        //
        //    Input, int ML, MU, the lower and upper bandwidths of A1.
        //    ML and MU must be nonnegative, and no greater than min(M,N)-1.
        //
        //    Input, double A[(2*ML+MU+1)*N], the R8GB matrix.
        //
        //    Output, double R8GB_TO_R8GE[M*N], the R8GE matrix.
        //
    {
        int i;
        int j;

        double[] b = new double[m * n];

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < m; i++)
            {
                b[i + j * m] = 0.0;
            }
        }

        for (i = 1; i <= m; i++)
        {
            for (j = 1; j <= n; j++)
            {
                if (i - ml <= j && j <= i + mu)
                {
                    b[i - 1 + (j - 1) * m] = a[ml + mu + i - j + (j - 1) * (2 * ml + mu + 1)];
                }
                else
                {
                    b[i - 1 + (j - 1) * m] = 0.0;
                }
            }
        }

        return b;
    }

    public static int r8gb_trf(int m, int n, int ml, int mu, ref double[] a, ref int[] pivot)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GB_TRF performs a LAPACK-style PLU factorization of an R8GB matrix.
        //
        //  Discussion:
        //
        //    The R8GB storage format is used for an M by N banded matrix, with lower
        //    bandwidth ML and upper bandwidth MU.  Storage includes room for ML extra
        //    superdiagonals, which may be required to store nonzero entries generated 
        //    during Gaussian elimination.
        //
        //    The original M by N matrix is "collapsed" downward, so that diagonals
        //    become rows of the storage array, while columns are preserved.  The
        //    collapsed array is logically 2*ML+MU+1 by N.  
        //
        //    The two dimensional array can be further reduced to a one dimensional
        //    array, stored by columns.
        //
        //    This is a simplified, standalone version of the LAPACK
        //    routine SGBTRF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 November 2003
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Anderson, Bai, Bischof, Blackford,
        //    Demmel, Dongarra, DuCroz, Greenbaum, Hammarling, McKenney, Sorensen.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Edward Anderson, Zhaojun Bai, Christian Bischof, Susan Blackford, 
        //    James Demmel, Jack Dongarra, Jeremy Du Croz, Anne Greenbaum, 
        //    Sven Hammarling, Alan McKenney, Danny Sorensen,
        //    LAPACK User's Guide,
        //    Second Edition,
        //    SIAM, 1995.
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows of the matrix A.  0 <= M.
        //
        //    Input, int N, the number of columns of the matrix A.  0 <= N.
        //
        //    Input, int ML, the number of subdiagonals within the band of A.
        //    0 <= ML.
        //
        //    Input, int MU, the number of superdiagonals within the band of A.
        //    0 <= MU.
        //
        //    Input/output, double A[(2*ML+MU+1)*N].  On input, the matrix A in band 
        //    storage, and on output, information about the PLU factorization.
        //
        //    Output, int PIVOT(min(M,N)), the pivot indices;
        //    for 1 <= i <= min(M,N), row i of the matrix was interchanged with
        //    row IPIV(i).
        //
        //    Output, int R8GB_TRF, error flag.
        //    = 0: successful exit;
        //    < 0: an input argument was illegal;
        //    > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
        //         has been completed, but the factor U is exactly
        //         singular, and division by zero will occur if it is used
        //         to solve a system of equations.
        //
    {
        int i;
        int j;

        int info = 0;
        //
        //  KV is the number of superdiagonals in the factor U, allowing for fill-in.
        //
        int kv = mu + ml;
        //
        //  Set fill-in elements in columns MU+2 to KV to zero.
        //
        for (j = mu + 2; j <= Math.Min(kv, n); j++)
        {
            for (i = kv - j + 1; i <= ml; i++)
            {
                a[i - 1 + (j - 1) * (2 * ml + mu + 1)] = 0.0;
            }
        }

        //
        //  JU is the index of the last column affected by the current stage
        //  of the factorization.
        //
        int ju = 1;

        for (j = 1; j <= Math.Min(m, n); j++)
        {
            //
            //  Set the fill-in elements in column J+KV to zero.
            //
            if (j + kv <= n)
            {
                for (i = 1; i <= ml; i++)
                {
                    a[i - 1 + (j + kv - 1) * (2 * ml + mu + 1)] = 0.0;
                }
            }

            //
            //  Find the pivot and test for singularity.
            //  KM is the number of subdiagonal elements in the current column.
            //
            int km = Math.Min(ml, m - j);

            double piv = Math.Abs(a[kv + (j - 1) * (2 * ml + mu + 1)]);
            int jp = kv + 1;

            for (i = kv + 2; i <= kv + km + 1; i++)
            {
                if (!(piv < Math.Abs(a[i - 1 + (j - 1) * (2 * ml + mu + 1)])))
                {
                    continue;
                }

                piv = Math.Abs(a[i - 1 + (j - 1) * (2 * ml + mu + 1)]);
                jp = i;
            }

            jp -= kv;

            pivot[j - 1] = jp + j - 1;

            if (a[kv + jp - 1 + (j - 1) * (2 * ml + mu + 1)] != 0.0)
            {
                ju = Math.Max(ju, Math.Min(j + mu + jp - 1, n));
                //
                //  Apply interchange to columns J to JU.
                //
                if (jp != 1)
                {
                    for (i = 0; i <= ju - j; i++)
                    {
                        (a[kv + jp - i - 1 + (j + i - 1) * (2 * ml + mu + 1)], a[kv + 1 - i - 1 + (j + i - 1) * (2 * ml + mu + 1)]) = (a[kv + 1 - i - 1 + (j + i - 1) * (2 * ml + mu + 1)], a[kv + jp - i - 1 + (j + i - 1) * (2 * ml + mu + 1)]);
                    }
                }

                switch (km)
                {
                    //
                    //  Compute the multipliers.
                    //
                    case > 0:
                    {
                        for (i = kv + 2; i <= kv + km + 1; i++)
                        {
                            a[i - 1 + (j - 1) * (2 * ml + mu + 1)] /= a[kv + (j - 1) * (2 * ml + mu + 1)];
                        }

                        //
                        //  Update the trailing submatrix within the band.
                        //
                        if (j < ju)
                        {
                            int k;
                            for (k = 1; k <= ju - j; k++)
                            {
                                if (a[kv - k + (j + k - 1) * (2 * ml + mu + 1)] == 0.0)
                                {
                                    continue;
                                }

                                for (i = 1; i <= km; i++)
                                {
                                    a[kv + i - k + (j + k - 1) * (2 * ml + mu + 1)] -= a[kv + i + (j - 1) * (2 * ml + mu + 1)] *
                                        a[kv - k + (j + k - 1) * (2 * ml + mu + 1)];
                                }
                            }
                        }

                        break;
                    }
                }
            }
            else
                //
                //  If pivot is zero, set INFO to the index of the pivot
                //  unless a zero pivot has already been found.
                //
            {
                info = info switch
                {
                    0 => j,
                    _ => info
                };
            }
        }

        return info;
    }

    public static double[] r8gb_trs(int n, int ml, int mu, int nrhs, char trans, double[] a,
            int[] pivot, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GB_TRS solves a linear system factored by R8GB_TRF.
        //
        //  Discussion:
        //
        //    The R8GB storage format is used for an M by N banded matrix, with lower
        //    bandwidth ML and upper bandwidth MU.  Storage includes room for ML extra
        //    superdiagonals, which may be required to store nonzero entries generated 
        //    during Gaussian elimination.
        //
        //    The original M by N matrix is "collapsed" downward, so that diagonals
        //    become rows of the storage array, while columns are preserved.  The
        //    collapsed array is logically 2*ML+MU+1 by N.  
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    26 January 2004
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Anderson, Bai, Bischof, Blackford,
        //    Demmel, Dongarra, DuCroz, Greenbaum, Hammarling, McKenney, Sorensen.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Edward Anderson, Zhaojun Bai, Christian Bischof, Susan Blackford, 
        //    James Demmel, Jack Dongarra, Jeremy Du Croz, Anne Greenbaum, 
        //    Sven Hammarling, Alan McKenney, Danny Sorensen,
        //    LAPACK User's Guide,
        //    Second Edition,
        //    SIAM, 1995.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix A.
        //    N must be positive.
        //
        //    Input, int ML, the number of subdiagonals within the band of A.
        //    ML must be at least 0, and no greater than N - 1.
        //
        //    Input, int MU, the number of superdiagonals within the band of A.
        //    MU must be at least 0, and no greater than N - 1.
        //
        //    Input, int NRHS, the number of right hand sides and the number of
        //    columns of the matrix B.  NRHS must be positive.
        //
        //    Input, char TRANS, specifies the form of the system.
        //    'N':  A * x = b  (No transpose)
        //    'T':  A'* X = B  (Transpose)
        //    'C':  A'* X = B  (Conjugate transpose = Transpose)
        //
        //    Input, double A[(2*ML+MU+1)*N], the LU factorization of the band matrix
        //    A, computed by R8GB_TRF.  
        //
        //    Input, int PIVOT[N], the pivot indices; for 1 <= I <= N, row I
        //    of the matrix was interchanged with row PIVOT(I).
        //
        //    Input, double B[N*NRHS], the right hand side vectors.
        //
        //    Output, double R8GB_TRS[N*NRHS], the solution vectors.
        //
    {
        int i;
        int j;
        int k;
        int l;
        int lm;
        double temp;
        //
        //  Test the input parameters.
        //
        if (trans != 'N' && trans != 'n' &&
            trans != 'T' && trans != 't' &&
            trans != 'C' && trans != 'c')
        {
            return null;
        }

        switch (n)
        {
            case <= 0:
                return null;
        }
        switch (ml)
        {
            case < 0:
                return null;
        }
        switch (mu)
        {
            case < 0:
                return null;
        }
        switch (nrhs)
        {
            case <= 0:
                return null;
        }

        double[] x = new double[n * nrhs];

        for (k = 0; k < nrhs; k++)
        {
            for (i = 0; i < n; i++)
            {
                x[i + k * n] = b[i + k * n];
            }
        }

        int kd = mu + ml + 1;
        switch (trans)
        {
            //
            //  Solve A * x = b.
            //
            //  Solve L * x = b, overwriting b with x.
            //
            //  L is represented as a product of permutations and unit lower
            //  triangular matrices L = P(1) * L(1) * ... * P(n-1) * L(n-1),
            //  where each transformation L(i) is a rank-one modification of
            //  the identity matrix.
            //
            case 'N':
            case 'n':
            {
                switch (ml)
                {
                    case > 0:
                    {
                        for (j = 1; j <= n - 1; j++)
                        {
                            lm = Math.Min(ml, n - j);
                            l = pivot[j - 1];

                            for (k = 0; k < nrhs; k++)
                            {
                                temp = x[l - 1 + k * n];
                                x[l - 1 + k * n] = x[j - 1 + k * n];
                                x[j - 1 + k * n] = temp;
                            }

                            for (k = 0; k < nrhs; k++)
                            {
                                if (x[j - 1 + k * n] != 0.0)
                                {
                                    for (i = 1; i <= lm; i++)
                                    {
                                        x[j + i - 1 + k * n] -= a[kd + i - 1 + (j - 1) * (2 * ml + mu + 1)] *
                                                                x[j - 1 + k * n];
                                    }
                                }
                            }
                        }

                        break;
                    }
                }

                //
                //  Solve U * x = b, overwriting b with x.
                //
                for (k = 0; k < nrhs; k++)
                {
                    for (j = n; 1 <= j; j--)
                    {
                        if (x[j - 1 + k * n] == 0.0)
                        {
                            continue;
                        }

                        l = ml + mu + 1 - j;
                        x[j - 1 + k * n] /= a[ml + mu + (j - 1) * (2 * ml + mu + 1)];
                        for (i = j - 1; Math.Max(1, j - ml - mu) <= i; i--)
                        {
                            x[i - 1 + k * n] -= a[l + i - 1 + (j - 1) * (2 * ml + mu + 1)] * x[j - 1 + k * n];
                        }
                    }
                }

                break;
            }
            default:
            {
                //
                //  Solve A' * x = b.
                //
                //  Solve U' * x = b, overwriting b with x.
                //
                for (k = 0; k < nrhs; k++)
                {
                    for (j = 1; j <= n; j++)
                    {
                        temp = x[j - 1 + k * n];
                        l = ml + mu + 1 - j;
                        for (i = Math.Max(1, j - ml - mu); i <= j - 1; i++)
                        {
                            temp -= a[l + i - 1 + (j - 1) * (2 * ml + mu + 1)] * x[i - 1 + k * n];
                        }

                        x[j - 1 + k * n] = temp / a[ml + mu + (j - 1) * (2 * ml + mu + 1)];
                    }
                }

                switch (ml)
                {
                    //
                    //  Solve L' * x = b, overwriting b with x.
                    //
                    case > 0:
                    {
                        for (j = n - 1; 1 <= j; j--)
                        {
                            lm = Math.Min(ml, n - j);

                            for (k = 0; k < nrhs; k++)
                            {
                                for (i = 1; i <= lm; i++)
                                {
                                    x[j - 1 + k * n] -= x[j + i - 1 + k * n] * a[kd + i - 1 + (j - 1) * (2 * ml + mu + 1)];
                                }
                            }

                            l = pivot[j - 1];

                            for (k = 0; k < nrhs; k++)
                            {
                                temp = x[l - 1 + k * n];
                                x[l - 1 + k * n] = x[j - 1 + k * n];
                                x[j - 1 + k * n] = temp;
                            }
                        }

                        break;
                    }
                }

                break;
            }
        }

        return x;
    }

    public static double[] r8gb_zeros(int m, int n, int ml, int mu)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GB_ZEROS zeros an R8GB matrix.
        //
        //  Discussion:
        //
        //    The R8GB storage format is used for an M by N banded matrix, with lower
        //    bandwidth ML and upper bandwidth MU.  Storage includes room for ML extra
        //    superdiagonals, which may be required to store nonzero entries generated 
        //    during Gaussian elimination.
        //
        //    The original M by N matrix is "collapsed" downward, so that diagonals
        //    become rows of the storage array, while columns are preserved.  The
        //    collapsed array is logically 2*ML+MU+1 by N.  
        //
        //    The two dimensional array can be further reduced to a one dimensional
        //    array, stored by columns.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows of the matrix.
        //    M must be nonnegative.
        //
        //    Input, int N, the number of columns of the matrix.
        //    N must be nonnegative.
        //
        //    Input, int ML, MU, the lower and upper bandwidths.
        //    ML and MU must be nonnegative and no greater than min(M,N)-1.
        //
        //    Output, double R8GB_ZERO[(2*ML+MU+1)*N], the R8GB matrix.
        //
    {
        int col = 2 * ml + mu + 1;
        int j;

        double[] a = new double[col * n];

        for (j = 0; j < n; j++)
        {
            int row;
            for (row = 0; row < col; row++)
            {
                a[row + j * col] = 0.0;
            }
        }

        return a;
    }

}