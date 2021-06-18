namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static void r8vm_sl(int n, double[] a, double[] b, int job, double[] x, ref int info )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VM_SL solves a R8VM system.
        //
        //  Discussion:
        //
        //    The R8VM storage format is used for an M by N Vandermonde matrix.
        //    An M by N Vandermonde matrix is defined by the values in its second
        //    row, which will be written here as X(1:N).  The matrix has a first 
        //    row of 1's, a second row equal to X(1:N), a third row whose entries
        //    are the squares of the X values, up to the M-th row whose entries
        //    are the (M-1)th powers of the X values.  The matrix can be stored
        //    compactly by listing just the values X(1:N).
        //
        //    Vandermonde systems are very close to singularity.  The singularity
        //    gets worse as N increases, and as any pair of values defining
        //    the matrix get close.  Even a system as small as N = 10 will
        //    involve the 9th power of the defining values.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 September 2003
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Golub, VanLoan.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Gene Golub, Charles Van Loan,
        //    Matrix Computations,
        //    Third Edition,
        //    Johns Hopkins, 1996.
        //
        //  Parameters:
        //
        //    Input, int N, the number of rows and columns of the matrix.
        //
        //    Input, double A[N], the R8VM matrix.
        //
        //    Input, double B[N], the right hand side.
        //
        //    Input, int JOB, specifies the system to solve.
        //    0, solve A * x = b.
        //    nonzero, solve A' * x = b.
        //
        //    Output, int &INFO.
        //    0, no error.
        //    nonzero, at least two of the values in A are equal.
        //
        //    Output, double X[N], the solution of the linear system.
        //
        {
            int i;
            int j;
            //
            //  Check for explicit singularity.
            //
            info = 0;

            for (j = 0; j < n; j++)
            {
                for (i = j + 1; i < n; i++)
                {
                    if (a[i] == a[j])
                    {
                        info = 1;
                        return;
                    }
                }
            }

            for (i = 0; i < n; i++)
            {
                x[i] = b[i];
            }

            if (job == 0)
            {
                for (j = 1; j <= n - 1; j++)
                {
                    for (i = n; j + 1 <= i; i--)
                    {
                        x[i - 1] = x[i - 1] - a[j - 1] * x[i - 2];
                    }
                }

                for (j = n - 1; 1 <= j; j--)
                {
                    for (i = j + 1; i <= n; i++)
                    {
                        x[i - 1] = x[i - 1] / (a[i - 1] - a[i - j - 1]);
                    }

                    for (i = j; i <= n - 1; i++)
                    {
                        x[i - 1] = x[i - 1] - x[i];
                    }
                }
            }
            else
            {
                for (j = 1; j <= n - 1; j++)
                {
                    for (i = n; j + 1 <= i; i--)
                    {
                        x[i - 1] = (x[i - 1] - x[i - 2]) / (a[i - 1] - a[i - j - 1]);
                    }
                }

                for (j = n - 1; 1 <= j; j--)
                {
                    for (i = j; i <= n - 1; i++)
                    {
                        x[i - 1] = x[i - 1] - x[i] * a[j - 1];
                    }
                }

            }
        }

        public static double[] r8vm_sl_new(int n, double[] a, double[] b, int job, ref int info )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VM_SL_NEW solves a R8VM system.
        //
        //  Discussion:
        //
        //    The R8VM storage format is used for an M by N Vandermonde matrix.
        //    An M by N Vandermonde matrix is defined by the values in its second
        //    row, which will be written here as X(1:N).  The matrix has a first 
        //    row of 1's, a second row equal to X(1:N), a third row whose entries
        //    are the squares of the X values, up to the M-th row whose entries
        //    are the (M-1)th powers of the X values.  The matrix can be stored
        //    compactly by listing just the values X(1:N).
        //
        //    Vandermonde systems are very close to singularity.  The singularity
        //    gets worse as N increases, and as any pair of values defining
        //    the matrix get close.  Even a system as small as N = 10 will
        //    involve the 9th power of the defining values.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 September 2003
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Golub, VanLoan.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Gene Golub, Charles Van Loan,
        //    Matrix Computations,
        //    Third Edition,
        //    Johns Hopkins, 1996.
        //
        //  Parameters:
        //
        //    Input, int N, the number of rows and columns of the matrix.
        //
        //    Input, double A[N], the R8VM matrix.
        //
        //    Input, double B[N], the right hand side.
        //
        //    Input, int JOB, specifies the system to solve.
        //    0, solve A * x = b.
        //    nonzero, solve A' * x = b.
        //
        //    Output, int &INFO.
        //    0, no error.
        //    nonzero, at least two of the values in A are equal.
        //
        //    Output, double R8VM_SL_NEW[N], the solution of the linear system.
        //
        {
            int i;
            int j;
            double[] x;
            //
            //  Check for explicit singularity.
            //
            info = 0;

            for (j = 0; j < n; j++)
            {
                for (i = j + 1; i < n; i++)
                {
                    if (a[i] == a[j])
                    {
                        info = 1;
                        return null;
                    }
                }
            }

            x = new double[n];

            for (i = 0; i < n; i++)
            {
                x[i] = b[i];
            }

            if (job == 0)
            {
                for (j = 1; j <= n - 1; j++)
                {
                    for (i = n; j + 1 <= i; i--)
                    {
                        x[i - 1] = x[i - 1] - a[j - 1] * x[i - 2];
                    }
                }

                for (j = n - 1; 1 <= j; j--)
                {
                    for (i = j + 1; i <= n; i++)
                    {
                        x[i - 1] = x[i - 1] / (a[i - 1] - a[i - j - 1]);
                    }

                    for (i = j; i <= n - 1; i++)
                    {
                        x[i - 1] = x[i - 1] - x[i];
                    }
                }
            }
            else
            {
                for (j = 1; j <= n - 1; j++)
                {
                    for (i = n; j + 1 <= i; i--)
                    {
                        x[i - 1] = (x[i - 1] - x[i - 2]) / (a[i - 1] - a[i - j - 1]);
                    }
                }

                for (j = n - 1; 1 <= j; j--)
                {
                    for (i = j; i <= n - 1; i++)
                    {
                        x[i - 1] = x[i - 1] - x[i] * a[j - 1];
                    }
                }

            }

            return x;
        }
    }
}