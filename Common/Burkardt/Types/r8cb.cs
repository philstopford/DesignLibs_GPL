using System;
using Burkardt.Uniform;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static double r8cb_det(int n, int ml, int mu, double[] a_lu)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8CB_DET computes the determinant of an R8CB matrix factored by R8CB_NP_FA.
            //
            //  Discussion:
            //
            //    The R8CB storage format is appropriate for a compact banded matrix.
            //    It is assumed that the matrix has lower and upper bandwidths ML and MU,
            //    respectively.  The matrix is stored in a way similar to that used
            //    by LINPACK and LAPACK for a general banded matrix, except that in
            //    this mode, no extra rows are set aside for possible fillin during pivoting.
            //    Thus, this storage format is suitable if you do not intend to factor
            //    the matrix, or if you can guarantee that the matrix can be factored
            //    without pivoting.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    25 March 2004
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
            //    Input, double A_LU[(ML+MU+1)*N], the LU factors from R8CB_FA.
            //
            //    Output, double R8CB_DET, the determinant of the matrix.
            //
        {
            double det;
            int j;

            det = 1.0;
            for (j = 0; j < n; j++)
            {
                det = det * a_lu[mu + j * (ml + mu + 1)];
            }

            return det;
        }

        public static double[] r8cb_dif2(int m, int n, int ml, int mu)

            //*****************************************************************************/
            //
            //  Purpose:
            //
            //    R8CB_DIF2 sets up an R8CB second difference matrix.
            //
            //  Discussion:
            //
            //    The R8CB storage format is appropriate for a compact banded matrix.
            //    It is assumed that the matrix has lower and upper bandwidths ML and MU,
            //    respectively.  The matrix is stored in a way similar to that used
            //    by LINPACK and LAPACK for a general banded matrix, except that in
            //    this mode, no extra rows are set aside for possible fillin during pivoting.
            //    Thus, this storage format is suitable if you do not intend to factor
            //    the matrix, or if you can guarantee that the matrix can be factored
            //    without pivoting.
            //
            //    The original M by N matrix is "collapsed" downward, so that diagonals
            //    become rows of the storage array, while columns are preserved.  The
            //    collapsed array is logically ML+MU+1 by N.  
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    27 July 2016
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
            //    Output, double R8CB_DIF2[(ML+MU+1)*N], the matrix.
            //
        {
            double[] a;
            int diag;
            int i;
            int j;

            a = r8vec_zeros_new((ml + mu + 1) * n);

            for (j = 0; j < n; j++)
            {
                for (diag = 0; diag < ml + mu + 1; diag++)
                {
                    i = diag + j - mu;

                    if (i == j)
                    {
                        a[diag + j * (ml + mu + 1)] = 2.0;
                    }
                    else if (i == j + 1 || i == j - 1)
                    {
                        a[diag + j * (ml + mu + 1)] = -1.0;
                    }
                }
            }

            return a;
        }

        public static double[] r8cb_indicator(int m, int n, int ml, int mu)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8CB_INDICATOR sets up an R8CB indicator matrix.
            //
            //  Discussion:
            //
            //    The R8CB storage format is appropriate for a compact banded matrix.
            //    It is assumed that the matrix has lower and upper bandwidths ML and MU,
            //    respectively.  The matrix is stored in a way similar to that used
            //    by LINPACK and LAPACK for a general banded matrix, except that in
            //    this mode, no extra rows are set aside for possible fillin during pivoting.
            //    Thus, this storage format is suitable if you do not intend to factor
            //    the matrix, or if you can guarantee that the matrix can be factored
            //    without pivoting.
            //
            //    The original M by N matrix is "collapsed" downward, so that diagonals
            //    become rows of the storage array, while columns are preserved.  The
            //    collapsed array is logically ML+MU+1 by N.  
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    23 January 2004
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
            //    Output, double R8CB_INDICATOR[(ML+MU+1)*N], the R8CB matrix.
            //
        {
            double[] a;
            int col = ml + mu + 1;
            int diag;
            int fac;
            int i;
            int j;
            int k;

            a = r8vec_zeros_new((ml + mu + 1) * n);

            fac = (int) Math.Pow(10, (int) Math.Log10(n) + 1);
            k = 0;

            for (j = 1; j <= n; j++)
            {
                for (diag = 1; diag <= ml + mu + 1; diag++)
                {
                    i = diag + j - mu - 1;

                    if (1 <= i && i <= m && i - ml <= j && j <= i + mu)
                    {
                        a[diag - 1 + (j - 1) * col] = (double) (fac * i + j);
                    }
                    else
                    {
                        k = k + 1;
                        a[diag - 1 + (j - 1) * col] = -(double) k;
                    }
                }
            }

            return a;
        }

        public static double[] r8cb_ml(int n, int ml, int mu, double[] a_lu, double[] x, int job)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8CB_ML computes A * x or A' * X, using R8CB_NP_FA factors.
            //
            //  Discussion:
            //
            //    The R8CB storage format is appropriate for a compact banded matrix.
            //    It is assumed that the matrix has lower and upper bandwidths ML and MU,
            //    respectively.  The matrix is stored in a way similar to that used
            //    by LINPACK and LAPACK for a general banded matrix, except that in
            //    this mode, no extra rows are set aside for possible fillin during pivoting.
            //    Thus, this storage format is suitable if you do not intend to factor
            //    the matrix, or if you can guarantee that the matrix can be factored
            //    without pivoting.
            //
            //    It is assumed that R8CB_NP_FA has overwritten the original matrix
            //    information by LU factors.  R8CB_ML is able to reconstruct the
            //    original matrix from the LU factor data.
            //
            //    R8CB_ML allows the user to check that the solution of a linear
            //    system is correct, without having to save an unfactored copy
            //    of the matrix.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    19 November 2003
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
            //    Input, double A_LU[(ML+MU+1)*N], the LU factors from R8CB_NP_FA.
            //
            //    Input, double X[N], the vector to be multiplied.
            //
            //    Input, int JOB, specifies the operation to be done:
            //    JOB = 0, compute A * x.
            //    JOB nonzero, compute A' * x.
            //
            //    Output, double R8CB_ML[N], the result of the multiplication.
            //
        {
            double[] b;
            int i;
            int ihi;
            int ilo;
            int j;
            int jhi;
            int nrow = ml + mu + 1;

            b = r8vec_zeros_new(n);

            for (i = 0; i < n; i++)
            {
                b[i] = x[i];
            }

            if (job == 0)
            {
                //
                //  Y = U * X.
                //
                for (j = 0; j < n; j++)
                {
                    ilo = Math.Max(0, j - mu);
                    for (i = ilo; i < j; i++)
                    {
                        b[i] = b[i] + a_lu[i - j + mu + j * nrow] * b[j];
                    }

                    b[j] = a_lu[j - j + mu + j * nrow] * b[j];
                }

                //
                //  B = PL * Y = PL * U * X = A * x.
                //
                for (j = n - 2; 0 <= j; j--)
                {
                    ihi = Math.Min(n - 1, j + ml);
                    for (i = j + 1; i <= ihi; i++)
                    {
                        b[i] = b[i] - a_lu[i - j + mu + j * nrow] * b[j];
                    }
                }
            }
            else
            {
                //
                //  Y = ( PL )' * X.
                //
                for (j = 0; j < n - 1; j++)
                {
                    jhi = Math.Min(n - 1, j + ml);
                    for (i = j + 1; i <= jhi; i++)
                    {
                        b[j] = b[j] - b[i] * a_lu[i - j + mu + j * nrow];
                    }
                }

                //
                //  B = U' * Y = ( PL * U )' * X = A' * X.
                //
                for (i = n - 1; 0 <= i; i--)
                {
                    jhi = Math.Min(n - 1, i + mu);
                    for (j = i + 1; j <= jhi; j++)
                    {
                        b[j] = b[j] + b[i] * a_lu[i - j + mu + j * nrow];
                    }

                    b[i] = b[i] * a_lu[i - i + mu + i * nrow];
                }
            }

            return b;
        }

        public static double[] r8cb_mtv(int m, int n, int ml, int mu, double[] a, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8CB_MTV multiplies a vector by an R8CB matrix.
            //
            //  Discussion:
            //
            //    The R8CB storage format is appropriate for a compact banded matrix.
            //    It is assumed that the matrix has lower and upper bandwidths ML and MU,
            //    respectively.  The matrix is stored in a way similar to that used
            //    by LINPACK and LAPACK for a general banded matrix, except that in
            //    this mode, no extra rows are set aside for possible fillin during pivoting.
            //    Thus, this storage format is suitable if you do not intend to factor
            //    the matrix, or if you can guarantee that the matrix can be factored
            //    without pivoting.
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
            //    Input, int M, N, the order of the matrix.
            //
            //    Input, int ML, MU, the lower and upper bandwidths.
            //    ML and MU must be nonnegative, and no greater than N-1.
            //
            //    Input, double A[(ML+MU+1)*N], the R8CB matrix.
            //
            //    Input, double X[M], the vector to be multiplied by A.
            //
            //    Output, double R8CB_MTV[N], the product X*A.
            //
        {
            double[] b;
            int i;
            int j;
            int jhi;
            int jlo;

            b = r8vec_zeros_new(n);

            for (i = 0; i < m; i++)
            {
                jlo = Math.Max(0, i - ml);
                jhi = Math.Min(n - 1, i + mu);
                for (j = jlo; j <= jhi; j++)
                {
                    b[j] = b[j] + x[i] * a[i - j + mu + j * (ml + mu + 1)];
                }
            }

            return b;
        }

        public static double[] r8cb_mv(int m, int n, int ml, int mu, double[] a, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8CB_MV multiplies an R8CB matrix times a vector.
            //
            //  Discussion:
            //
            //    The R8CB storage format is appropriate for a compact banded matrix.
            //    It is assumed that the matrix has lower and upper bandwidths ML and MU,
            //    respectively.  The matrix is stored in a way similar to that used
            //    by LINPACK and LAPACK for a general banded matrix, except that in
            //    this mode, no extra rows are set aside for possible fillin during pivoting.
            //    Thus, this storage format is suitable if you do not intend to factor
            //    the matrix, or if you can guarantee that the matrix can be factored
            //    without pivoting.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 October 2003
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
            //    Input, double A[(ML+MU+1)*N], the R8CB matrix.
            //
            //    Input, double X[N], the vector to be multiplied by A.
            //
            //    Output, double R8CB_MV[M], the product A * x.
            //
        {
            double[] b;
            int i;
            int j;
            int jhi;
            int jlo;

            b = r8vec_zeros_new(m);

            for (i = 0; i < m; i++)
            {
                jlo = Math.Max(0, i - ml);
                jhi = Math.Min(n - 1, i + mu);
                for (j = jlo; j <= jhi; j++)
                {
                    b[i] = b[i] + a[i - j + mu + j * (ml + mu + 1)] * x[j];
                }
            }

            return b;
        }

        public static int r8cb_np_fa(int n, int ml, int mu, ref double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8CB_NP_FA factors an R8CB matrix by Gaussian elimination.
            //
            //  Discussion:
            //
            //    The R8CB storage format is appropriate for a compact banded matrix.
            //    It is assumed that the matrix has lower and upper bandwidths ML and MU,
            //    respectively.  The matrix is stored in a way similar to that used
            //    by LINPACK and LAPACK for a general banded matrix, except that in
            //    this mode, no extra rows are set aside for possible fillin during pivoting.
            //    Thus, this storage format is suitable if you do not intend to factor
            //    the matrix, or if you can guarantee that the matrix can be factored
            //    without pivoting.
            //
            //    R8CB_NP_FA is a version of the LINPACK routine SGBFA, modifed to use
            //    no pivoting, and to be applied to the R8CB compressed band matrix storage
            //    format.  It will fail if the matrix is singular, or if any zero
            //    pivot is encountered.
            //
            //    If R8CB_NP_FA successfully factors the matrix, R8CB_NP_SL may be called
            //    to solve linear systems involving the matrix.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    14 October 2003
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
            //    Input/output, double A[(ML+MU+1)*N], the compact band matrix.
            //    On input, the coefficient matrix of the linear system.
            //    On output, the LU factors of the matrix.
            //
            //    Output, int R8CB_NP_FA, singularity flag.
            //    0, no singularity detected.
            //    nonzero, the factorization failed on the INFO-th step.
            //
        {
            int i;
            int j;
            int ju;
            int k;
            int lm;
            int m;
            int mm;
            //
            //  The value of M is MU + 1 rather than ML + MU + 1.
            //
            m = mu + 1;
            ju = 0;

            for (k = 1; k <= n - 1; k++)
            {
                //
                //  If our pivot entry A(MU+1,K) is zero, then we must give up.
                //
                if (a[m - 1 + (k - 1) * (ml + mu + 1)] == 0.0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("R8CB_FA - Fatal error!");
                    Console.WriteLine("  Zero pivot on step " + k + "");
                    return 1;
                }

                //
                //  LM counts the number of nonzero elements that lie below the current
                //  diagonal entry, A(K,K).
                //
                //  Multiply the LM entries below the diagonal by -1/A(K,K), turning
                //  them into the appropriate "multiplier" terms in the L matrix.
                //
                lm = Math.Min(ml, n - k);
                for (i = m + 1; i <= m + lm; i++)
                {
                    a[i - 1 + (k - 1) * (ml + mu + 1)] =
                        -a[i - 1 + (k - 1) * (ml + mu + 1)] / a[m - 1 + (k - 1) * (ml + mu + 1)];
                }

                //
                //  MM points to the row in which the next entry of the K-th row is, A(K,J).
                //  We then add L(I,K)*A(K,J) to A(I,J) for rows I = K+1 to K+LM.
                //
                ju = Math.Max(ju, mu + k);
                ju = Math.Min(ju, n);
                mm = m;

                for (j = k + 1; j <= ju; j++)
                {
                    mm = mm - 1;
                    for (i = 1; i <= lm; i++)
                    {
                        a[mm + i - 1 + (j - 1) * (ml + mu + 1)] = a[mm + i - 1 + (j - 1) * (ml + mu + 1)]
                                                                  + a[mm - 1 + (j - 1) * (ml + mu + 1)] *
                                                                  a[m + i - 1 + (k - 1) * (ml + mu + 1)];
                    }
                }
            }

            if (a[m - 1 + (n - 1) * (ml + mu + 1)] == 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("R8CB_FA - Fatal error!");
                Console.WriteLine("  Zero pivot on step " + n + "");
                return (1);
            }

            return 0;
        }

        public static double[] r8cb_np_sl(int n, int ml, int mu, double[] a_lu, double[] b, int job)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8CB_NP_SL solves an R8CB system factored by R8CB_NP_FA.
            //
            //  Discussion:
            //
            //    The R8CB storage format is appropriate for a compact banded matrix.
            //    It is assumed that the matrix has lower and upper bandwidths ML and MU,
            //    respectively.  The matrix is stored in a way similar to that used
            //    by LINPACK and LAPACK for a general banded matrix, except that in
            //    this mode, no extra rows are set aside for possible fillin during pivoting.
            //    Thus, this storage format is suitable if you do not intend to factor
            //    the matrix, or if you can guarantee that the matrix can be factored
            //    without pivoting.
            //
            //    R8CB_NP_SL can also solve the related system A' * x = b.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    18 January 2004
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
            //    Input, double A_LU[(ML+MU+1)*N], the LU factors from R8CB_NP_FA.
            //
            //    Input, double B[N], the right hand side of the linear system.
            //
            //    Input, int JOB.
            //    If JOB is zero, the routine will solve A * x = b.
            //    If JOB is nonzero, the routine will solve A' * x = b.
            //
            //    Output, double R8CB_NP_SL[N], the solution of the linear system, X.
            //
        {
            int i;
            int k;
            int la;
            int lb;
            int lm;
            int m;
            double[] x;

            x = r8vec_zeros_new(n);

            for (i = 0; i < n; i++)
            {
                x[i] = b[i];
            }

            m = mu + 1;
            //
            //  Solve A * x = b.
            //
            if (job == 0)
            {
                //
                //  Solve PL * Y = B.
                //
                if (0 < ml)
                {
                    for (k = 1; k <= n - 1; k++)
                    {
                        lm = Math.Min(ml, n - k);
                        for (i = 0; i < lm; i++)
                        {
                            x[k + i] = x[k + i] + x[k - 1] * a_lu[m + i + (k - 1) * (ml + mu + 1)];
                        }
                    }
                }

                //
                //  Solve U * X = Y.
                //
                for (k = n; 1 <= k; k--)
                {
                    x[k - 1] = x[k - 1] / a_lu[m - 1 + (k - 1) * (ml + mu + 1)];
                    lm = Math.Min(k, m) - 1;
                    la = m - lm;
                    lb = k - lm;
                    for (i = 0; i <= lm - 1; i++)
                    {
                        x[lb + i - 1] = x[lb + i - 1] - x[k - 1] * a_lu[la + i - 1 + (k - 1) * (ml + mu + 1)];
                    }
                }
            }
            //
            //  Solve A' * X = B.
            //
            else
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
                        x[k - 1] = x[k - 1] - a_lu[la + i - 1 + (k - 1) * (ml + mu + 1)] * x[lb + i - 1];
                    }

                    x[k - 1] = x[k - 1] / a_lu[m - 1 + (k - 1) * (ml + mu + 1)];

                }

                //
                //  Solve ( PL )' * X = Y.
                //
                if (0 < ml)
                {
                    for (k = n - 1; 1 <= k; k--)
                    {
                        lm = Math.Min(ml, n - k);
                        for (i = 0; i < lm; i++)
                        {
                            x[k - 1] = x[k - 1] + a_lu[m + i + (k - 1) * (ml + mu + 1)] * x[k + i];
                        }
                    }
                }
            }

            return x;
        }

        public static void r8cb_print(int m, int n, int ml, int mu, double[] a, string title)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8CB_PRINT prints an R8CB matrix.
            //
            //  Discussion:
            //
            //    The R8CB storage format is appropriate for a compact banded matrix.
            //    It is assumed that the matrix has lower and upper bandwidths ML and MU,
            //    respectively.  The matrix is stored in a way similar to that used
            //    by LINPACK and LAPACK for a general banded matrix, except that in
            //    this mode, no extra rows are set aside for possible fillin during pivoting.
            //    Thus, this storage format is suitable if you do not intend to factor
            //    the matrix, or if you can guarantee that the matrix can be factored
            //    without pivoting.
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
            //    Input, int M, N, the number of rows and columns of the matrix.
            //
            //    Input, int ML, MU, the lower and upper bandwidths.
            //    ML and MU must be nonnegative, and no greater than min(M,N)-1..
            //
            //    Input, double A[(ML+MU+1)*N], the R8CB matrix.
            //
            //    Input, string TITLE, a title.
            //
        {
            r8cb_print_some(m, n, ml, mu, a, 1, 1, m, n, title);

            return;
        }

        public static void r8cb_print_some(int m, int n, int ml, int mu, double[] a, int ilo,
                int jlo, int ihi, int jhi, string title)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8CB_PRINT_SOME prints some of an R8CB matrix.
            //
            //  Discussion:
            //
            //    The R8CB storage format is appropriate for a compact banded matrix.
            //    It is assumed that the matrix has lower and upper bandwidths ML and MU,
            //    respectively.  The matrix is stored in a way similar to that used
            //    by LINPACK and LAPACK for a general banded matrix, except that in
            //    this mode, no extra rows are set aside for possible fillin during pivoting.
            //    Thus, this storage format is suitable if you do not intend to factor
            //    the matrix, or if you can guarantee that the matrix can be factored
            //    without pivoting.
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
            //    Input, int M, N, the number of rows and columns of the matrix.
            //
            //    Input, int ML, MU, the lower and upper bandwidths.
            //    ML and MU must be nonnegative, and no greater than min(M,N)-1.
            //
            //    Input, double A[(ML+MU+1)*N], the R8CB matrix.
            //
            //    Input, int ILO, JLO, IHI, JHI, designate the first row and
            //    column, and the last row and column to be printed.
            //
            //    Input, string TITLE, a title.
            //
        {
            int INCX = 5;

            int i;
            int i2hi;
            int i2lo;
            int j;
            int j2hi;
            int j2lo;
            string cout = "";

            Console.WriteLine("");
            Console.WriteLine(title + "");
            //
            //  Print the columns of the matrix, in strips of 5.
            //
            for (j2lo = jlo; j2lo <= jhi; j2lo = j2lo + INCX)
            {
                j2hi = j2lo + INCX - 1;
                j2hi = Math.Min(j2hi, n);
                j2hi = Math.Min(j2hi, jhi);

                Console.WriteLine("");
                cout = "  Col: ";

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
                i2lo = Math.Max(ilo, 1);
                i2lo = Math.Max(i2lo, j2lo - mu);
                i2hi = Math.Min(ihi, m);
                i2hi = Math.Min(i2hi, j2hi + ml);

                for (i = i2lo; i <= i2hi; i++)
                {
                    cout = i.ToString().PadLeft(4) + "  ";
                    //
                    //  Print out (up to) 5 entries in row I, that lie in the current strip.
                    //
                    for (j = j2lo; j <= j2hi; j++)
                    {
                        if (ml < i - j || mu < j - i)
                        {
                            cout += "              ";
                        }
                        else
                        {
                            cout += a[i - j + mu + (j - 1) * (ml + mu + 1)].ToString().PadLeft(12) + "  ";
                        }
                    }

                    Console.WriteLine(cout);
                }
            }
        }

        public static double[] r8cb_random(int m, int n, int ml, int mu, ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8CB_RANDOM randomizes an R8CB matrix.
            //
            //  Discussion:
            //
            //    The R8CB storage format is appropriate for a compact banded matrix.
            //    It is assumed that the matrix has lower and upper bandwidths ML and MU,
            //    respectively.  The matrix is stored in a way similar to that used
            //    by LINPACK and LAPACK for a general banded matrix, except that in
            //    this mode, no extra rows are set aside for possible fillin during pivoting.
            //    Thus, this storage format is suitable if you do not intend to factor
            //    the matrix, or if you can guarantee that the matrix can be factored
            //    without pivoting.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    25 July 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the order of the matrix.
            //
            //    Input, int ML, MU, the lower and upper bandwidths.
            //    ML and MU must be nonnegative, and no greater than N-1.
            //
            //    Input/output, int &SEED, a seed for the random number generator.
            //
            //    Output, double R8CB_RANDOM[(ML+MU+1)*N], the R8CB matrix.
            //
        {
            double[] a;
            int i;
            int ihi;
            int ilo;
            int j;

            a = r8vec_zeros_new((ml + mu + 1) * n);
            //
            //  Set the entries that correspond to matrix elements.
            //
            for (j = 0; j < n; j++)
            {
                ilo = Math.Max(0, j - mu);
                ihi = Math.Min(m - 1, j + ml);

                for (i = ilo; i <= ihi; i++)
                {
                    a[i - j + mu + j * (ml + mu + 1)] = UniformRNG.r8_uniform_01(ref seed);
                }
            }

            return a;
        }

        public static double[] r8cb_to_r8vec(int m, int n, int ml, int mu, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8CB_TO_R8VEC copies an R8CB matrix to a real vector.
            //
            //  Discussion:
            //
            //    In C++ and FORTRAN, this routine is not really needed.  In MATLAB,
            //    a data item carries its dimensionality implicitly, and so cannot be
            //    regarded sometimes as a vector and sometimes as an array.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    27 July 2016
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
            //    Input, double A[(ML+MU+1)*N], the array to be copied.
            //
            //    Output, double R8CB_TO_R8VEC[(ML+MU+1)*N], the vector.
            //
        {
            int i;
            int ihi;
            int ilo;
            int j;
            double[] x;

            x = r8vec_zeros_new((ml + mu + 1) * n);

            for (j = 0; j < n; j++)
            {
                ilo = Math.Max(mu - j, 0);
                ihi = mu + Math.Min(ml, m - j - 1);
                for (i = ilo; i <= ihi; i++)
                {
                    x[i + j * (ml + mu + 1)] = a[i + j * (ml + mu + 1)];
                }
            }

            return x;
        }

        public static double[] r8cb_to_r8ge(int m, int n, int ml, int mu, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8CB_TO_R8GE copies an R8CB matrix to an R8GE matrix.
            //
            //  Discussion:
            //
            //    The R8CB storage format is appropriate for a compact banded matrix.
            //    It is assumed that the matrix has lower and upper bandwidths ML and MU,
            //    respectively.  The matrix is stored in a way similar to that used
            //    by LINPACK and LAPACK for a general banded matrix, except that in
            //    this mode, no extra rows are set aside for possible fillin during pivoting.
            //    Thus, this storage format is suitable if you do not intend to factor
            //    the matrix, or if you can guarantee that the matrix can be factored
            //    without pivoting.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    17 November 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, N, the order of the matrices.
            //
            //    Input, int ML, MU, the lower and upper bandwidths of A.
            //    ML and MU must be nonnegative, and no greater than N-1.
            //
            //    Input, double A[(ML+MU+1)*N], the R8CB matrix.
            //
            //    Output, double R8CB_TO_R8GE[M*N], the R8GE matrix.
            //
        {
            double[] b;
            int i;
            int j;

            b = r8vec_zeros_new(m * n);

            for (i = 0; i < m; i++)
            {
                for (j = 0; j < n; j++)
                {
                    if (j - mu <= i && i <= j + ml)
                    {
                        b[i + j * m] = a[mu + i - j + j * (ml + mu + 1)];
                    }
                }
            }

            return b;
        }

        public static double[] r8cb_zeros(int n, int ml, int mu)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8CB_ZEROS zeros an R8CB matrix.
            //
            //  Discussion:
            //
            //    The R8CB storage format is appropriate for a compact banded matrix.
            //    It is assumed that the matrix has lower and upper bandwidths ML and MU,
            //    respectively.  The matrix is stored in a way similar to that used
            //    by LINPACK and LAPACK for a general banded matrix, except that in
            //    this mode, no extra rows are set aside for possible fillin during pivoting.
            //    Thus, this storage format is suitable if you do not intend to factor
            //    the matrix, or if you can guarantee that the matrix can be factored
            //    without pivoting.
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
            //    N must be nonnegative.
            //
            //    Input, int ML, MU, the lower and upper bandwidths.
            //    ML and MU must be nonnegative and no greater than N-1.
            //
            //    Output, double R8CB_ZERO[(ML+MU+1)*N), the R8CB matrix.
            //
        {
            double[] a;

            a = r8vec_zeros_new((ml + mu + 1) * n);

            return a;
        }

    }
}