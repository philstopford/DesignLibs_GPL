using System;
using Burkardt.MatrixNS;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static double[] r8mat_solve_svd ( int m, int n, double[] a, double[] b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_SOLVE_SVD solves a linear system A*x=b using the SVD.
        //
        //  Discussion:
        //
        //    When the system is determined, the solution is the solution in the
        //    ordinary sense, and A*x = b.
        //
        //    When the system is overdetermined, the solution minimizes the
        //    L2 norm of the residual ||A*x-b||.
        //
        //    When the system is underdetermined, ||A*x-b|| should be zero, and
        //    the solution is the solution of minimum L2 norm, ||x||.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 June 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns
        //    in the matrix A.
        //
        //    Input, double A[M,*N], the matrix.
        //
        //    Input, double B[M], the right hand side.
        //
        //    Output, double R8MAT_SOLVE_SVD[N], the solution.
        //
    {
        int i;
        int j;
        //
        //  Compute the SVD decomposition.
        //
        double[] a_copy = r8mat_copy_new(m, n, a);
        int lda = m;
        double[] sdiag = new double[Math.Max(m + 1, n)];
        double[] e = new double[Math.Max(m + 1, n)];
        double[] u = new double[m * m];
        int ldu = m;
        double[] v = new double[n * n];
        int ldv = n;
        double[] work = new double[m];
        const int job = 11;

        int info = DSVDC.dsvdc(ref a_copy, lda, m, n, ref sdiag, ref e, ref u, ldu, ref v, ldv, work, job);

        if (info != 0)
        {
            Console.WriteLine("");
            Console.WriteLine("R8MAT_SOLVE_SVD - Fatal error!");
            Console.WriteLine("  The SVD could not be calculated.");
            Console.WriteLine("  LINPACK routine DSVDC returned a nonzero");
            Console.WriteLine("  value of the error flag, INFO = " + info + "");
            return null;
        }

        double[] s = new double [m * n];

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < m; i++)
            {
                s[i + j * m] = 0.0;
            }
        }

        for (i = 0; i < Math.Min(m, n); i++)
        {
            s[i + i * m] = sdiag[i];
        }

        //
        //  Compute the pseudo inverse.
        //
        double[] sp = new double [n * m];

        for (j = 0; j < m; j++)
        {
            for (i = 0; i < n; i++)
            {
                sp[i + j * m] = 0.0;
            }
        }

        for (i = 0; i < Math.Min(m, n); i++)
        {
            if (s[i + i * m] != 0.0)
            {
                sp[i + i * n] = 1.0 / s[i + i * m];
            }
        }

        double[] a_pseudo = new double[n * m];

        for (j = 0; j < m; j++)
        {
            for (i = 0; i < n; i++)
            {
                a_pseudo[i + j * n] = 0.0;
                int k;
                for (k = 0; k < n; k++)
                {
                    int l;
                    for (l = 0; l < m; l++)
                    {
                        a_pseudo[i + j * n] += v[i + k * n] * sp[k + l * n] * u[j + l * m];
                    }
                }
            }
        }

        //
        //  Compute x = A_pseudo * b.
        //
        double[] x = r8mat_mv_new(n, m, a_pseudo, b);

        return x;
    }

    public static int r8mat_solve(int n, int rhs_num, ref double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_SOLVE uses Gauss-Jordan elimination to solve an N by N linear system.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //    Entry A(I,J) is stored as A[I+J*N]
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 August 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, int RHS_NUM, the number of right hand sides.  RHS_NUM
        //    must be at least 0.
        //
        //    Input/output, double A[N*(N+RHS_NUM)], contains in rows and columns 1
        //    to N the coefficient matrix, and in columns N+1 through
        //    N+RHS_NUM, the right hand sides.  On output, the coefficient matrix
        //    area has been destroyed, while the right hand sides have
        //    been overwritten with the corresponding solutions.
        //
        //    Output, int R8MAT_SOLVE, singularity flag.
        //    0, the matrix was not singular, the solutions were computed;
        //    J, factorization failed on step J, and the solutions could not
        //    be computed.
        //
    {
        int j;

        for (j = 0; j < n; j++)
        {
            //
            //  Choose a pivot row.
            //
            int ipivot = j;
            double apivot = a[j + j * n];

            int i;
            for (i = j; i < n; i++)
            {
                if (!(Math.Abs(apivot) < Math.Abs(a[i + j * n])))
                {
                    continue;
                }

                apivot = a[i + j * n];
                ipivot = i;
            }

            switch (apivot)
            {
                case 0.0:
                    return j;
            }

            //
            //  Interchange.
            //
            for (i = 0; i < n + rhs_num; i++)
            {
                (a[ipivot + i * n], a[j + i * n]) = (a[j + i * n], a[ipivot + i * n]);
            }

            //
            //  A(J,J) becomes 1.
            //
            a[j + j * n] = 1.0;
            int k;
            for (k = j; k < n + rhs_num; k++)
            {
                a[j + k * n] /= apivot;
            }

            //
            //  A(I,J) becomes 0.
            //
            for (i = 0; i < n; i++)
            {
                if (i == j)
                {
                    continue;
                }

                double factor = a[i + j * n];
                a[i + j * n] = 0.0;
                for (k = j; k < n + rhs_num; k++)
                {
                    a[i + k * n] -= factor * a[j + k * n];
                }
            }
        }

        return 0;
    }

    public static double[] r8mat_solve_2d(double[] a, double[] b, ref double det)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_SOLVE_2D solves a 2 by 2 linear system using Cramer's rule.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //    If the determinant DET is returned as zero, then the matrix A is
        //    singular, and does not have an inverse.  In that case, X is
        //    returned as the NULL vector.
        //
        //    If DET is nonzero, then its value is roughly an estimate
        //    of how nonsingular the matrix A is.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 November 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A[2*2], the matrix.
        //
        //    Input, double B[2], the right hand side.
        //
        //    Output, double *DET, the determinant of the system.
        //
        //    Output, double R8MAT_SOLVE_2D[2], the solution of the system,
        //    if DET is nonzero.  Otherwise, the NULL vector.
        //
    {
        //
        //  Compute the determinant.
        //
        det = a[0 + 0 * 2] * a[1 + 1 * 2] - a[0 + 1 * 2] * a[1 + 0 * 2];
        switch (det)
        {
            //
            //  If the determinant is zero, bail out.
            //
            case 0.0:
                return null;
        }

        //
        //  Compute the solution.
        //
        double[] x = new double[2];

        x[0] = (a[1 + 1 * 2] * b[0] - a[0 + 1 * 2] * b[1]) / det;
        x[1] = (-a[1 + 0 * 2] * b[0] + a[0 + 0 * 2] * b[1]) / det;

        return x;
    }

    public static double[] r8mat_solve_3d(double[] a, double[] b, ref double det)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_SOLVE_3D solves a 3 by 3 linear system using Cramer's rule.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //    If the determinant DET is returned as zero, then the matrix A is
        //    singular, and does not have an inverse.  In that case, X is
        //    returned as the NULL vector.
        //
        //    If DET is nonzero, then its value is roughly an estimate
        //    of how nonsingular the matrix A is.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 December 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A[3*3], the matrix.
        //
        //    Input, double B[3], the right hand side.
        //
        //    Output, double *DET, the determinant of the system.
        //
        //    Output, double R8MAT_SOLVE_3D[3], the solution of the system,
        //    if DET is nonzero.  Otherwise, the NULL vector.
        //
    {
        //
        //  Compute the determinant.
        //
        det = a[0 + 0 * 3] * (a[1 + 1 * 3] * a[2 + 2 * 3] - a[1 + 2 * 3] * a[2 + 1 * 3])
              + a[0 + 1 * 3] * (a[1 + 2 * 3] * a[2 + 0 * 3] - a[1 + 0 * 3] * a[2 + 2 * 3])
              + a[0 + 2 * 3] * (a[1 + 0 * 3] * a[2 + 1 * 3] - a[1 + 1 * 3] * a[2 + 0 * 3]);
        switch (det)
        {
            //
            //  If the determinant is zero, bail out.
            //
            case 0.0:
                return null;
        }

        //
        //  Compute the solution.
        //
        double[] x = new double[3];

        x[0] = ((a[1 + 1 * 3] * a[2 + 2 * 3] - a[1 + 2 * 3] * a[2 + 1 * 3]) * b[0]
                - (a[0 + 1 * 3] * a[2 + 2 * 3] - a[0 + 2 * 3] * a[2 + 1 * 3]) * b[1]
                + (a[0 + 1 * 3] * a[1 + 2 * 3] - a[0 + 2 * 3] * a[1 + 1 * 3]) * b[2]) / det;

        x[1] = (-(a[1 + 0 * 3] * a[2 + 2 * 3] - a[1 + 2 * 3] * a[2 + 0 * 3]) * b[0]
                + (a[0 + 0 * 3] * a[2 + 2 * 3] - a[0 + 2 * 3] * a[2 + 0 * 3]) * b[1]
                - (a[0 + 0 * 3] * a[1 + 2 * 3] - a[0 + 2 * 3] * a[1 + 0 * 3]) * b[2]) / det;

        x[2] = ((a[1 + 0 * 3] * a[2 + 1 * 3] - a[1 + 1 * 3] * a[2 + 0 * 3]) * b[0]
                - (a[0 + 0 * 3] * a[2 + 1 * 3] - a[0 + 1 * 3] * a[2 + 0 * 3]) * b[1]
                + (a[0 + 0 * 3] * a[1 + 1 * 3] - a[0 + 1 * 3] * a[1 + 0 * 3]) * b[2]) / det;

        return x;
    }

    public static double[] r8mat_l_solve(int n, double[] a, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_L_SOLVE solves a lower triangular linear system.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of rows and columns of the matrix A.
        //
        //    Input, double A[N*N], the N by N lower triangular matrix.
        //
        //    Input, double B[N], the right hand side of the linear system.
        //
        //    Output, double R8MAT_L_SOLVE[N], the solution of the linear system.
        //
    {
        int i;

        double[] x = new double[n];
        //
        //  Solve L * x = b.
        //
        for (i = 0; i < n; i++)
        {
            double temp = 0.0;
            int j;
            for (j = 0; j < i; j++)
            {
                temp += a[i + j * n] * x[j];
            }

            x[i] = (b[i] - temp) / a[i + i * n];
        }

        return x;
    }
        
    public static double[] r8mat_lt_solve(int n, double[] a, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_LT_SOLVE solves a transposed lower triangular linear system.
        //
        //  Discussion:
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //    Given the lower triangular matrix A, the linear system to be solved is:
        //
        //      A' * x = b
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 September 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of rows and columns of the matrix A.
        //
        //    Input, double A[N*N], the N by N lower triangular matrix.
        //
        //    Input, double B[N], the right hand side of the linear system.
        //
        //    Output, double R8MAT_LT_SOLVE[N], the solution of the linear system.
        //
    {
        int j;

        double[] x = new double[n];

        for (j = n - 1; 0 <= j; j--)
        {
            x[j] = b[j];
            int i;
            for (i = j + 1; i < n; i++)
            {
                x[j] -= x[i] * a[i + j * n];
            }

            x[j] /= a[j + j * n];
        }

        return x;
    }

    public static double[] r8mat_solve2(int n, ref double[] a, ref double[] b, ref int ierror)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_SOLVE2 computes the solution of an N by N linear system.
        //
        //  Discussion: 							    
        //
        //    An R8MAT is a doubly dimensioned array of R8 values, stored as a vector 
        //    in column-major order.
        //
        //    The linear system may be represented as
        //
        //      A*X = B
        //
        //    If the linear system is singular, but consistent, then the routine will
        //    still produce a solution.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of equations.
        //
        //    Input/output, double A[N*N].
        //    On input, A is the coefficient matrix to be inverted.
        //    On output, A has been overwritten.
        //
        //    Input/output, double B[N].
        //    On input, B is the right hand side of the system.
        //    On output, B has been overwritten.
        //
        //    Output, int *IERROR.
        //    0, no error detected.
        //    1, consistent singularity.
        //    2, inconsistent singularity.
        //
        //    Output, double R8MAT_SOLVE2[N], the solution of the linear system.
        //
    {
        int imax;

        ierror = 0;

        int[] piv = i4vec_zero_new(n);
        double[] x = r8vec_zero_new(n);
        //
        //  Process the matrix.
        //
        for (int k = 1; k <= n; k++)
        {
            //
            //  In column K:
            //    Seek the row IMAX with the properties that:
            //      IMAX has not already been used as a pivot;
            //      A(IMAX,K) is larger in magnitude than any other candidate.
            //
            double amax = 0.0;
            imax = 0;
            for (int i = 1; i <= n; i++)
            {
                switch (piv[i - 1])
                {
                    case 0:
                    {
                        if (amax < Math.Abs(a[i - 1 + (k - 1) * n]))
                        {
                            imax = i;
                            amax = Math.Abs(a[i - 1 + (k - 1) * n]);
                        }

                        break;
                    }
                }
            }

            //
            //  If you found a pivot row IMAX, then,
            //    eliminate the K-th entry in all rows that have not been used for pivoting.
            //
            if (imax == 0)
            {
                continue;
            }

            {
                piv[imax - 1] = k;
                for (int j = k + 1; j <= n; j++)
                {
                    a[imax - 1 + (j - 1) * n] /= a[imax - 1 + (k - 1) * n];
                }

                b[imax - 1] /= a[imax - 1 + (k - 1) * n];
                a[imax - 1 + (k - 1) * n] = 1.0;

                for (int i = 1; i <= n; i++)
                {
                    switch (piv[i - 1])
                    {
                        case 0:
                        {
                            for (int j = k + 1; j <= n; j++)
                            {
                                a[i - 1 + (j - 1) * n] -= a[i - 1 + (k - 1) * n] * a[imax - 1 + (j - 1) * n];
                            }

                            b[i - 1] -= a[i - 1 + (k - 1) * n] * b[imax - 1];
                            a[i - 1 + (k - 1) * n] = 0.0;
                            break;
                        }
                    }
                }
            }
        }

        //
        //  Now, every row with nonzero IPIV begins with a 1, and
        //  all other rows are all zero.  Begin solution.
        //
        for (int j = n; 1 <= j; j--)
        {
            imax = 0;
            for (int k = 1; k <= n; k++)
            {
                if (piv[k - 1] == j)
                {
                    imax = k;
                }
            }

            switch (imax)
            {
                case 0:
                {
                    x[j - 1] = 0.0;

                    switch (b[j - 1])
                    {
                        case 0.0:
                            ierror = 1;
                            Console.WriteLine("");
                            Console.WriteLine("R8MAT_SOLVE2 - Warning:");
                            Console.WriteLine("  Consistent singularity, equation = " + j + "");
                            break;
                        default:
                            ierror = 2;
                            Console.WriteLine("");
                            Console.WriteLine("R8MAT_SOLVE2 - Warning:");
                            Console.WriteLine("  Inconsistent singularity, equation = " + j + "");
                            break;
                    }

                    break;
                }
                default:
                {
                    x[j - 1] = b[imax - 1];

                    for (int i = 1; i <= n; i++)
                    {
                        if (i != imax)
                        {
                            b[i - 1] -= a[i - 1 + (j - 1) * n] * x[j - 1];
                        }
                    }

                    break;
                }
            }
        }

        return x;
    }
        
    public static double[] r8mat_ut_solve ( int n, double[] a, double[] b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_UT_SOLVE solves a transposed upper triangular linear system.
        //
        //  Discussion:
        //
        //    An R8MAT is an MxN array of R8's, stored by (I,J) -> [I+J*M].
        //
        //    Given the upper triangular matrix A, the linear system to be solved is:
        //
        //      A' * x = b
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 October 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of rows and columns of
        //    the matrix A.
        //
        //    Input, double A[N*N], the N by N upper triangular matrix.
        //
        //    Input, double B[N], the right hand side of the linear system.
        //
        //    Output, double R8MAT_UT_SOLVE[N], the solution of the linear system.
        //
    {
        int i;
        //
        //  Solve U' * x = b.
        //
        double[] x = new double[n];

        for ( i = 0; i < n; i++ )
        {
            x[i] = b[i];
            int j;
            for ( j = 0; j < i; j++ )
            {
                x[i] -= a[j+i*n] * x[j];
            }
            x[i] /= a[i+i*n];
        }

        return x;
    }
        
    public static double[] r8mat_upsol(int n, double[] r, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_UPSOL solves R * X = B, for an upper triangular matrix R.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 June 2013
        //
        //  Author:
        //
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
        //
        //    Input, double R[N*N], the upper triangular matrix.
        //
        //    Input, double B[N], the right hand side.
        //
        //    Output, double R8MAT_UPSOL[N], the solution.
        //
    {

        double[] x = new double[n];
        for (int j = 0; j < n; j++)
        {
            x[j] = b[j];
        }

        return x;
    }

    public static double[] r8mat_utsol(int n, double[] r, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_UTSOL solves R' * X = B for an upper triangular matrix R.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 June 2013
        //
        //  Author:
        //
        //    John Burkardt
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
        //
        //    Input, double R[N*N], the upper triangular matrix.
        //
        //    Input, double B[N], the right hand side.
        //
        //    Output, double X[N], the solution.
        //
    {
        double[] x = new double[n];

        for (int j = 0; j < n; j++)
        {
            x[j] = b[j];
        }

        x[0] /= r[0 + 0 * n];

        for (int j = 1; j < n; j++)
        {
            for (int i = 0; i < j; i++)
            {
                x[j] -= r[i + j * n] * x[i];
            }

            x[j] /= r[j + j * n];
        }

        return x;
    }
        
}