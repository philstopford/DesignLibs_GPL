﻿using System;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static double r83_np_det ( int n, double[] a_lu )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83_NP_DET: determinant of a tridiagonal system factored by R83_NP_FA.
        //
        //  Discussion:
        //
        //    The R83 storage format is used for a tridiagonal matrix.
        //    The superdiagonal is stored in entries (1,2:min(M+1,N)).
        //    The diagonal in entries (2,1:min(M,N)).
        //    The subdiagonal in (3,min(M-1,N)).
        //
        //  Example:
        //
        //    An R83 matrix of order 3x5 would be stored:
        //
        //       *  A12 A23 A34  *
        //      A11 A22 A33  *   *
        //      A21 A32  *   *   *
        //
        //    An R83 matrix of order 5x5 would be stored:
        //
        //       *  A12 A23 A34 A45
        //      A11 A22 A33 A44 A55
        //      A21 A32 A43 A54  *
        //
        //    An R83 matrix of order 5x3 would be stored:
        //
        //       *  A12 A23
        //      A11 A22 A33
        //      A21 A32 A43
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
        //    N must be at least 2.
        //
        //    Input, double A_LU[3*N], the LU factors from R83_NP_FA.
        //
        //    Output, double R83_NP_DET, the determinant of the matrix.
        //
    {
        int j;

        double det = 1.0;
        for ( j = 0; j < n; j++ )
        {
            det *= a_lu[1+j*3];
        }

        return det;
    }
    public static int r83_np_fa(int n, ref double[] a)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83_NP_FA factors a R83 system without pivoting.
        //
        //  Discussion:
        //
        //    The R83 storage format is used for a tridiagonal matrix.
        //    The superdiagonal is stored in entries (1,2:N), the diagonal in
        //    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
        //    original matrix is "collapsed" vertically into the array.
        //
        //    Because this routine does not use pivoting, it can fail even when
        //    the matrix is not singular, and it is liable to make larger
        //    errors.
        //
        //    R83_NP_FA and R83_NP_SL may be preferable to the corresponding
        //    LINPACK routine SGTSL for tridiagonal systems, which factors and solves
        //    in one step, and does not save the factorization.
        //
        //  Example:
        //
        //    Here is how a R83 matrix of order 5 would be stored:
        //
        //       *  A12 A23 A34 A45
        //      A11 A22 A33 A44 A55
        //      A21 A32 A43 A54  *
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 January 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be at least 2.
        //
        //    Input/output, double A[3*N].
        //    On input, the tridiagonal matrix.  On output, factorization information.
        //
        //    Output, int R83_NP_FA, singularity flag.
        //    0, no singularity detected.
        //    nonzero, the factorization failed on the INFO-th step.
        //
    {
        for (int i = 1; i <= n - 1; i++)
        {
            switch (a[1 + (i - 1) * 3])
            {
                case 0.0:
                    Console.WriteLine("");
                    Console.WriteLine("R83_NP_FA - Fatal error!");
                    Console.WriteLine("  Zero pivot on step " + i + "");
                    return i;
                default:
                    //
                    //  Store the multiplier in L.
                    //
                    a[2 + (i - 1) * 3] /= a[1 + (i - 1) * 3];
                    //
                    //  Modify the diagonal entry in the next column.
                    //
                    a[1 + i * 3] -= a[2 + (i - 1) * 3] * a[0 + i * 3];
                    break;
            }
        }

        switch (a[1 + (n - 1) * 3])
        {
            case 0.0:
                Console.WriteLine("");
                Console.WriteLine("R83_NP_FA - Fatal error!");
                Console.WriteLine("  Zero pivot on step " + n + "");
                return n;
            default:
                return 0;
        }
    }
        
    public static double[] r83_np_fs ( int n, double[] a, double[] b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83_NP_FS factors and solves an R83 system.
        //
        //  Discussion:
        //
        //    The R83 storage format is used for a tridiagonal matrix.
        //    The superdiagonal is stored in entries (1,2:min(M+1,N)).
        //    The diagonal in entries (2,1:min(M,N)).
        //    The subdiagonal in (3,min(M-1,N)).
        //
        //    This algorithm requires that each diagonal entry be nonzero.
        //    It does not use pivoting, and so can fail on systems that
        //    are actually nonsingular.
        //
        //  Example:
        //
        //    An R83 matrix of order 3x5 would be stored:
        //
        //       *  A12 A23 A34  *
        //      A11 A22 A33  *   *
        //      A21 A32  *   *   *
        //
        //    An R83 matrix of order 5x5 would be stored:
        //
        //       *  A12 A23 A34 A45
        //      A11 A22 A33 A44 A55
        //      A21 A32 A43 A54  *
        //
        //    An R83 matrix of order 5x3 would be stored:
        //
        //       *  A12 A23
        //      A11 A22 A33
        //      A21 A32 A43
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 November 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the linear system.
        //
        //    Input/output, double A[3*N].
        //    On input, the nonzero diagonals of the linear system.
        //    On output, the data in these vectors has been overwritten
        //    by factorization information.
        //
        //    Input, double B[N], the right hand side.
        //
        //    Output, double R83_NP_FS[N], the solution of the linear system.
        //    This is NULL if there was an error because one of the diagonal
        //    entries was zero.
        //
    {
        int i;
        //
        //  Check.
        //
        for (i = 0; i < n; i++)
        {
            switch (a[1 + i * 3])
            {
                case 0.0:
                    return null;
            }
        }

        double[] x = new double[n];

        for (i = 0; i < n; i++)
        {
            x[i] = b[i];
        }

        for (i = 1; i < n; i++)
        {
            double xmult = a[2 + (i - 1) * 3] / a[1 + (i - 1) * 3];
            a[1 + i * 3] -= xmult * a[0 + i * 3];
            x[i] -= xmult * x[i - 1];
        }

        x[n - 1] /= a[1 + (n - 1) * 3];
        for (i = n - 2; 0 <= i; i--)
        {
            x[i] = (x[i] - a[0 + (i + 1) * 3] * x[i + 1]) / a[1 + i * 3];
        }

        return x;
    }


    public static double[] r83_np_ml(int n, double[] a_lu, double[] x, int job)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83_NP_ML computes Ax or xA, where A has been factored by R83_NP_FA.
        //
        //  Discussion:
        //
        //    The R83 storage format is used for a tridiagonal matrix.
        //    The superdiagonal is stored in entries (1,2:min(M+1,N)).
        //    The diagonal in entries (2,1:min(M,N)).
        //    The subdiagonal in (3,min(M-1,N)).
        //
        //  Example:
        //
        //    An R83 matrix of order 3x5 would be stored:
        //
        //       *  A12 A23 A34  *
        //      A11 A22 A33  *   *
        //      A21 A32  *   *   *
        //
        //    An R83 matrix of order 5x5 would be stored:
        //
        //       *  A12 A23 A34 A45
        //      A11 A22 A33 A44 A55
        //      A21 A32 A43 A54  *
        //
        //    An R83 matrix of order 5x3 would be stored:
        //
        //       *  A12 A23
        //      A11 A22 A33
        //      A21 A32 A43
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be at least 2.
        //
        //    Input, double A_LU[3*N], the LU factors from R83_FA.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Output, double B[N], the product.
        //
        //    Input, int JOB, specifies the product to find.
        //    0, compute A * x.
        //    nonzero, compute A' * x.
        //
    {
        int i;

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
                //  Compute X := U * X
                //
                for (i = 1; i <= n; i++)
                {
                    b[i - 1] = a_lu[1 + (i - 1) * 3] * b[i - 1];

                    if (i < n)
                    {
                        b[i - 1] += a_lu[0 + i * 3] * b[i];
                    }
                }

                //
                //  Compute X: = L * X.
                //
                for (i = n; 2 <= i; i--)
                {
                    b[i - 1] += a_lu[2 + (i - 2) * 3] * b[i - 2];
                }

                break;
            }
            default:
            {
                //
                //  Compute X: = L' * X.
                //
                for (i = 1; i <= n - 1; i++)
                {
                    b[i - 1] += a_lu[2 + (i - 1) * 3] * b[i];
                }

                //
                //  Compute X: = U' * X.
                //
                for (i = n; 1 <= i; i--)
                {
                    b[i - 1] = a_lu[1 + (i - 1) * 3] * b[i - 1];
                    switch (i)
                    {
                        case > 1:
                            b[i - 1] += a_lu[0 + (i - 1) * 3] * b[i - 2];
                            break;
                    }
                }

                break;
            }
        }

        return b;
    }

    public static double[] r83_np_sl(int n, double[] a_lu, double[] b, int job)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83_NP_SL solves a R83 system factored by R83_NP_FA.
        //
        //  Discussion:
        //
        //    The R83 storage format is used for a tridiagonal matrix.
        //    The superdiagonal is stored in entries (1,2:N), the diagonal in
        //    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
        //    original matrix is "collapsed" vertically into the array.
        //
        //  Example:
        //
        //    Here is how a R83 matrix of order 5 would be stored:
        //
        //       *  A12 A23 A34 A45
        //      A11 A22 A33 A44 A55
        //      A21 A32 A43 A54  *
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 January 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be at least 2.
        //
        //    Input, double A_LU[3*N], the LU factors from R83_NP_FA.
        //
        //    Input, double B[N], the right hand side of the linear system.
        //    On output, B contains the solution of the linear system.
        //
        //    Input, int JOB, specifies the system to solve.
        //    0, solve A * x = b.
        //    nonzero, solve A' * x = b.
        //
        //    Output, double R83_NP_SL[N], the solution of the linear system.
        //
    {
        double[] x = new double[n];

        for (int i = 0; i < n; i++)
        {
            x[i] = b[i];
        }

        switch (job)
        {
            case 0:
            {
                //
                //  Solve L * Y = B.
                //
                for (int i = 1; i < n; i++)
                {
                    x[i] -= a_lu[2 + (i - 1) * 3] * x[i - 1];
                }

                //
                //  Solve U * X = Y.
                //
                for (int i = n; 1 <= i; i--)
                {
                    x[i - 1] /= a_lu[1 + (i - 1) * 3];
                    switch (i)
                    {
                        case > 1:
                            x[i - 2] -= a_lu[0 + (i - 1) * 3] * x[i - 1];
                            break;
                    }
                }

                break;
            }
            default:
            {
                //
                //  Solve U' * Y = B
                //
                for (int i = 1; i <= n; i++)
                {
                    x[i - 1] /= a_lu[1 + (i - 1) * 3];
                    if (i < n)
                    {
                        x[i] -= a_lu[0 + i * 3] * x[i - 1];
                    }
                }

                //
                //  Solve L' * X = Y.
                //
                for (int i = n - 1; 1 <= i; i--)
                {
                    x[i - 1] -= a_lu[2 + (i - 1) * 3] * x[i];
                }

                break;
            }
        }

        return x;
    }

    public static double[] r83np_fs(int n, ref double[] a, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83NP_FS factors and solves an R83NP system.
        //
        //  Discussion:
        //
        //    The R83NP storage format is used for a tridiagonal matrix.
        //    The subdiagonal   is in entries (0,1:N-1),
        //    the diagonal      is in entries (1,0:N-1),
        //    the superdiagonal is in entries (2,0:N-2).
        //
        //    This algorithm requires that each diagonal entry be nonzero.
        //    It does not use pivoting, and so can fail on systems that
        //    are actually nonsingular.
        //
        //    The "R83NP" format used for this routine is different from the R83 format.
        //    Here, we insist that the nonzero entries
        //    for a given row now appear in the corresponding column of the
        //    packed array.
        //
        //  Example:
        //
        //    Here is how a R83 matrix of order 5 would be stored:
        //
        //       *  A21 A32 A43 A54
        //      A11 A22 A33 A44 A55
        //      A12 A23 A34 A45  *
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 May 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the linear system.
        //
        //    Input/output, double A[3*N].
        //    On input, the nonzero diagonals of the linear system.
        //    On output, the data in these vectors has been overwritten
        //    by factorization information.
        //
        //    Input, double B[N], the right hand side.
        //
        //    Output, double R83NP_FS[N], the solution of the linear system.
        //
    {
        //
        //  Check.
        //
        for (int i = 0; i < n; i++)
        {
            switch (a[1 + i * 3])
            {
                case 0.0:
                    Console.WriteLine();
                    Console.WriteLine("R83NP_FS - Fatal error!");
                    Console.WriteLine("  A[1+" + i + "*3] = 0.");
                    return new double[1];
            }
        }

        double[] x = new double[n];

        for (int i = 0; i < n; i++)
        {
            x[i] = b[i];
        }

        for (int i = 1; i < n; i++)
        {
            a[1 + i * 3] -= a[2 + (i - 1) * 3] * a[0 + i * 3] / a[1 + (i - 1) * 3];
            x[i] -= x[i - 1] * a[0 + i * 3] / a[1 + (i - 1) * 3];
        }

        x[n - 1] /= a[1 + (n - 1) * 3];
        for (int i = n - 2; 0 <= i; i--)
        {
            x[i] = (x[i] - a[2 + i * 3] * x[i + 1]) / a[1 + i * 3];
        }

        return x;
    }


}