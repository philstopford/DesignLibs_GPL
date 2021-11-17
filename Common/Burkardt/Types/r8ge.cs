using System;
using Burkardt.Uniform;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static void r8ge_cg(int n, double[] a, double[] b, ref double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_CG uses the conjugate gradient method on an R8GE system.
        //
        //  Discussion:
        //
        //    The R8GE storage format is used for a general M by N matrix.  A storage 
        //    space is made for each entry.  The two dimensional logical
        //    array can be thought of as a vector of M*N entries, starting with
        //    the M entries in the column 1, then the M entries in column 2
        //    and so on.  Considered as a vector, the entry A(I,J) is then stored
        //    in vector location I+(J-1)*M.
        //
        //    The matrix A must be a positive definite symmetric band matrix.
        //
        //    The method is designed to reach the solution after N computational
        //    steps.  However, roundoff may introduce unacceptably large errors for
        //    some problems.  In such a case, calling the routine again, using
        //    the computed solution as the new starting estimate, should improve
        //    the results.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Frank Beckman,
        //    The Solution of Linear Equations by the Conjugate Gradient Method,
        //    in Mathematical Methods for Digital Computers,
        //    edited by John Ralston, Herbert Wilf,
        //    Wiley, 1967,
        //    ISBN: 0471706892,
        //    LC: QA76.5.R3.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be positive.
        //
        //    Input, double A[N*N], the matrix.
        //
        //    Input, double B[N], the right hand side vector.
        //
        //    Input/output, double X[N].
        //    On input, an estimate for the solution, which may be 0.
        //    On output, the approximate solution vector.
        //
    {
        int i;
        int it;
        //
        //  Initialize
        //    AP = A * x,
        //    R  = b - A * x,
        //    P  = b - A * x.
        //
        double[] ap = r8ge_mv(n, n, a, x);

        double[] r = new double[n];
        for (i = 0; i < n; i++)
        {
            r[i] = b[i] - ap[i];
        }

        double[] p = new double[n];
        for (i = 0; i < n; i++)
        {
            p[i] = b[i] - ap[i];
        }

        //
        //  Do the N steps of the conjugate gradient method.
        //
        for (it = 1; it <= n; it++)
        {
            //
            //  Compute the matrix*vector product AP=A*P.
            //
            ap = r8ge_mv(n, n, a, p);
            //
            //  Compute the dot products
            //    PAP = P*AP,
            //    PR  = P*R
            //  Set
            //    ALPHA = PR / PAP.
            //
            double pap = r8vec_dot_product(n, p, ap);
            double pr = r8vec_dot_product(n, p, r);

            if (pap == 0.0)
            {
                break;
            }

            double alpha = pr / pap;
            //
            //  Set
            //    X = X + ALPHA * P
            //    R = R - ALPHA * AP.
            //
            for (i = 0; i < n; i++)
            {
                x[i] += alpha * p[i];
            }

            for (i = 0; i < n; i++)
            {
                r[i] -= alpha * ap[i];
            }

            //
            //  Compute the vector dot product
            //    RAP = R*AP
            //  Set
            //    BETA = - RAP / PAP.
            //
            double rap = r8vec_dot_product(n, r, ap);

            double beta = -rap / pap;
            //
            //  Update the perturbation vector
            //    P = R + BETA * P.
            //
            for (i = 0; i < n; i++)
            {
                p[i] = r[i] + beta * p[i];
            }
        }
    }

    public static double r8ge_co(int n, ref double[] a, ref int[] pivot)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_CO factors an R8GE matrix and estimates its condition number.
        //
        //  Discussion:
        //
        //    The R8GE storage format is used for a "general" M by N matrix.  
        //    A physical storage space is made for each logical entry.  The two 
        //    dimensional logical array is mapped to a vector, in which storage is 
        //    by columns.
        //
        //    For the system A * X = B, relative perturbations in A and B
        //    of size EPSILON may cause relative perturbations in X of size
        //    EPSILON/RCOND.
        //
        //    If RCOND is so small that the logical expression
        //      1.0 + rcond == 1.0
        //    is true, then A may be singular to working precision.  In particular,
        //    RCOND is zero if exact singularity is detected or the estimate
        //    underflows.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 October 2003
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
        //    Input, int N, the order of the matrix A.
        //
        //    Input/output, double A[N*N].  On input, a matrix to be factored.
        //    On output, the LU factorization of the matrix.
        //
        //    Output, int PIVOT[N], the pivot indices.
        //
        //    Output, double R8GE_CO, an estimate of the reciprocal condition number of A.
        //
    {
        int i;
        int j;
        int k;
        int l;
        double rcond;
        double s;
        double t;
        //
        //  Compute the L1 norm of A.
        //
        double anorm = 0.0;
        for (j = 0; j < n; j++)
        {
            s = 0.0;
            for (i = 0; i < n; i++)
            {
                s += Math.Abs(a[i + j * n]);
            }

            anorm = Math.Max(anorm, s);
        }

        //
        //  Compute the LU factorization.
        //
        int info = r8ge_fa(n, ref a, ref pivot);

        if (info != 0)
        {
            rcond = 0.0;
            return rcond;
        }

        //
        //  RCOND = 1 / ( norm(A) * (estimate of norm(inverse(A))) )
        //
        //  estimate of norm(inverse(A)) = norm(Z) / norm(Y)
        //
        //  where
        //    A * Z = Y
        //  and
        //    A' * Y = E
        //
        //  The components of E are chosen to cause maximum local growth in the
        //  elements of W, where U'*W = E.  The vectors are frequently rescaled
        //  to avoid overflow.
        //
        //  Solve U' * W = E.
        //
        double ek = 1.0;
        double[] z = new double[n];
        for (i = 0; i < n; i++)
        {
            z[i] = 0.0;
        }

        for (k = 0; k < n; k++)
        {
            if (z[k] != 0.0)
            {
                ek = -r8_sign(z[k]) * Math.Abs(ek);
            }

            if (Math.Abs(a[k + k * n]) < Math.Abs(ek - z[k]))
            {
                s = Math.Abs(a[k + k * n]) / Math.Abs(ek - z[k]);
                for (i = 0; i < n; i++)
                {
                    z[i] = s * z[i];
                }

                ek = s * ek;
            }

            double wk = ek - z[k];
            double wkm = -ek - z[k];
            s = Math.Abs(wk);
            double sm = Math.Abs(wkm);

            if (a[k + k * n] != 0.0)
            {
                wk /= a[k + k * n];
                wkm /= a[k + k * n];
            }
            else
            {
                wk = 1.0;
                wkm = 1.0;
            }

            if (k + 2 <= n)
            {
                for (j = k + 1; j < n; j++)
                {
                    sm += Math.Abs(z[j] + wkm * a[k + j * n]);
                    z[j] += wk * a[k + j * n];
                    s += Math.Abs(z[j]);
                }

                if (s < sm)
                {
                    t = wkm - wk;
                    wk = wkm;
                    for (j = k + 1; j < n; j++)
                    {
                        z[j] += t * a[k + j * n];
                    }
                }
            }

            z[k] = wk;
        }

        s = 0.0;
        for (i = 0; i < n; i++)
        {
            s += Math.Abs(z[i]);
        }

        for (i = 0; i < n; i++)
        {
            z[i] /= s;
        }

        //
        //  Solve L' * Y = W
        //
        for (k = n - 1; 0 <= k; k--)
        {
            for (i = k + 1; i < n; i++)
            {
                z[k] += z[i] * a[i + k * n];
            }

            t = Math.Abs(z[k]);
            switch (t)
            {
                case > 1.0:
                {
                    for (i = 0; i < n; i++)
                    {
                        z[i] /= t;
                    }

                    break;
                }
            }

            l = pivot[k] - 1;

            t = z[l];
            z[l] = z[k];
            z[k] = t;
        }

        t = 0.0;
        for (i = 0; i < n; i++)
        {
            t += Math.Abs(z[i]);
        }

        for (i = 0; i < n; i++)
        {
            z[i] /= t;
        }

        double ynorm = 1.0;
        //
        //  Solve L * V = Y.
        //
        for (k = 0; k < n; k++)
        {
            l = pivot[k] - 1;

            t = z[l];
            z[l] = z[k];
            z[k] = t;

            for (i = k + 1; i < n; i++)
            {
                z[i] += t * a[i + k * n];
            }

            t = Math.Abs(z[k]);

            switch (t)
            {
                case > 1.0:
                {
                    ynorm /= t;
                    for (i = 0; i < n; i++)
                    {
                        z[i] /= t;
                    }

                    break;
                }
            }
        }

        s = 0.0;
        for (i = 0; i < n; i++)
        {
            s += Math.Abs(z[i]);
        }

        for (i = 0; i < n; i++)
        {
            z[i] /= s;
        }

        ynorm /= s;
        //
        //  Solve U * Z = V.
        //
        for (k = n - 1; 0 <= k; k--)
        {
            if (Math.Abs(a[k + k * n]) < Math.Abs(z[k]))
            {
                s = Math.Abs(a[k + k * n]) / Math.Abs(z[k]);
                for (i = 0; i < n; i++)
                {
                    z[i] = s * z[i];
                }

                ynorm = s * ynorm;
            }

            if (a[k + k * n] != 0.0)
            {
                z[k] /= a[k + k * n];
            }
            else
            {
                z[k] = 1.0;
            }

            for (i = 0; i < k; i++)
            {
                z[i] -= a[i + k * n] * z[k];
            }
        }

        //
        //  Normalize Z in the L1 norm.
        //
        s = 0.0;
        for (i = 0; i < n; i++)
        {
            s += Math.Abs(z[i]);
        }

        s = 1.0 / s;

        for (i = 0; i < n; i++)
        {
            z[i] = s * z[i];
        }

        ynorm = s * ynorm;

        if (anorm != 0.0)
        {
            rcond = ynorm / anorm;
        }
        else
        {
            rcond = 0.0;
        }

        return rcond;
    }

    public static void r8ge_copy(int m, int n, double[] a1, ref double[] a2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_COPY copies one R8GE to a "new" R8GE.
        //
        //  Discussion:
        //
        //    An R8GE is a doubly dimensioned array of R8's, which
        //    may be stored as a vector in column-major order.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 July 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double A1[M*N], the matrix to be copied.
        //
        //    Output, double A2[M*N], the copy of A1.
        //
    {
        int j;

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                a2[i + j * m] = a1[i + j * m];
            }
        }
    }

    public static double[] r8ge_copy_new(int m, int n, double[] a1)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_COPY_NEW copies one R8GE to a "new" R8GE.
        //
        //  Discussion:
        //
        //    An R8GE is a doubly dimensioned array of R8's, which
        //    may be stored as a vector in column-major order.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 July 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double A1[M*N], the matrix to be copied.
        //
        //    Output, double R8GE_COPY_NEW[M*N], the copy of A1.
        //
    {
        int j;

        double[] a2 = new double[m * n];

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                a2[i + j * m] = a1[i + j * m];
            }
        }

        return a2;
    }

    public static double r8ge_det(int n, double[] a_lu, int[] pivot)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_DET computes the determinant of a matrix factored by R8GE_FA or R8GE_TRF.
        //
        //  Discussion:
        //
        //    The R8GE storage format is used for a "general" M by N matrix.  
        //    A physical storage space is made for each logical entry.  The two 
        //    dimensional logical array is mapped to a vector, in which storage is 
        //    by columns.
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
        //    Input, double A_LU[N*N], the LU factors from R8GE_FA or R8GE_TRF.
        //
        //    Input, int PIVOT[N], as computed by R8GE_FA or R8GE_TRF.
        //
        //    Output, double R8GE_DET, the determinant of the matrix.
        //
    {
        int i;

        double det = 1.0;

        for (i = 1; i <= n; i++)
        {
            det *= a_lu[i - 1 + (i - 1) * n];
            if (pivot[i - 1] != i)
            {
                det = -det;
            }
        }

        return det;
    }

    public static double[] r8ge_dif2(int m, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_DIF2 returns the DIF2 matrix in R8GE format.
        //
        //  Example:
        //
        //    N = 5
        //
        //    2 -1  .  .  .
        //   -1  2 -1  .  .
        //    . -1  2 -1  .
        //    .  . -1  2 -1
        //    .  .  . -1  2
        //
        //  Properties:
        //
        //    A is banded, with bandwidth 3.
        //
        //    A is tridiagonal.
        //
        //    Because A is tridiagonal, it has property A (bipartite).
        //
        //    A is a special case of the TRIS or tridiagonal scalar matrix.
        //
        //    A is integral, therefore det ( A ) is integral, and 
        //    det ( A ) * inverse ( A ) is integral.
        //
        //    A is Toeplitz: constant along diagonals.
        //
        //    A is symmetric: A' = A.
        //
        //    Because A is symmetric, it is normal.
        //
        //    Because A is normal, it is diagonalizable.
        //
        //    A is persymmetric: A(I,J) = A(N+1-J,N+1-I).
        //
        //    A is positive definite.
        //
        //    A is an M matrix.
        //
        //    A is weakly diagonally dominant, but not strictly diagonally dominant.
        //
        //    A has an LU factorization A = L * U, without pivoting.
        //
        //      The matrix L is lower bidiagonal with subdiagonal elements:
        //
        //        L(I+1,I) = -I/(I+1)
        //
        //      The matrix U is upper bidiagonal, with diagonal elements
        //
        //        U(I,I) = (I+1)/I
        //
        //      and superdiagonal elements which are all -1.
        //
        //    A has a Cholesky factorization A = L * L', with L lower bidiagonal.
        //
        //      L(I,I) =    sqrt ( (I+1) / I )
        //      L(I,I-1) = -sqrt ( (I-1) / I )
        //
        //    The eigenvalues are
        //
        //      LAMBDA(I) = 2 + 2 * COS(I*PI/(N+1))
        //                = 4 SIN^2(I*PI/(2*N+2))
        //
        //    The corresponding eigenvector X(I) has entries
        //
        //       X(I)(J) = sqrt(2/(N+1)) * sin ( I*J*PI/(N+1) ).
        //
        //    Simple linear systems:
        //
        //      x = (1,1,1,...,1,1),   A*x=(1,0,0,...,0,1)
        //
        //      x = (1,2,3,...,n-1,n), A*x=(0,0,0,...,0,n+1)
        //
        //    det ( A ) = N + 1.
        //
        //    The value of the determinant can be seen by induction,
        //    and expanding the determinant across the first row:
        //
        //      det ( A(N) ) = 2 * det ( A(N-1) ) - (-1) * (-1) * det ( A(N-2) )
        //                = 2 * N - (N-1)
        //                = N + 1
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 July 2000
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Robert Gregory, David Karney,
        //    A Collection of Matrices for Testing Computational Algorithms,
        //    Wiley, 1969,
        //    ISBN: 0882756494,
        //    LC: QA263.68
        //
        //    Morris Newman, John Todd,
        //    Example A8,
        //    The evaluation of matrix inversion programs,
        //    Journal of the Society for Industrial and Applied Mathematics,
        //    Volume 6, Number 4, pages 466-476, 1958.
        //
        //    John Todd,
        //    Basic Numerical Mathematics,
        //    Volume 2: Numerical Algebra,
        //    Birkhauser, 1980,
        //    ISBN: 0817608117,
        //    LC: QA297.T58.
        //
        //    Joan Westlake,
        //    A Handbook of Numerical Matrix Inversion and Solution of 
        //    Linear Equations,
        //    John Wiley, 1968,
        //    ISBN13: 978-0471936756,
        //    LC: QA263.W47.
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of the matrix.
        //
        //    Output, double R8GE_DIF2[M*N], the matrix.
        //
    {
        int j;

        double[] a = new double[m * n];

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                if (j == i - 1)
                {
                    a[i + j * m] = -1.0;
                }
                else if (j == i)
                {
                    a[i * j * m] = 2.0;
                }
                else if (j == i + 1)
                {
                    a[i + j * m] = -1.0;
                }
                else
                {
                    a[i + j * m] = 0.0;
                }
            }
        }

        return a;
    }

    public static double[] r8ge_mtm_new(int n1, int n2, int n3, double[] a, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_MTM_NEW computes C = A' * B.
        //
        //  Discussion:
        //
        //    An R8GE is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //    For this routine, the result is returned as the function value.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 September 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N1, N2, N3, the order of the matrices.
        //
        //    Input, double A[N2*N1], double B[N2*N3], the matrices to multiply.
        //
        //    Output, double R8GE_MTM_NEW[N1*N3], the product matrix C = A' * B.
        //
    {
        int i;

        double[] c = new double[n1 * n3];

        for (i = 0; i < n1; i++)
        {
            int j;
            for (j = 0; j < n3; j++)
            {
                c[i + j * n1] = 0.0;
                int k;
                for (k = 0; k < n2; k++)
                {
                    c[i + j * n1] += a[k + i * n2] * b[k + j * n2];
                }
            }
        }

        return c;
    }

    public static double[] r8ge_mtv(int m, int n, double[] a, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_MTV multiplies a vector times an R8GE matrix.
        //
        //  Discussion:
        //
        //    The R8GE storage format is used for a "general" M by N matrix.  
        //    A physical storage space is made for each logical entry.  The two 
        //    dimensional logical array is mapped to a vector, in which storage is 
        //    by columns.
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
        //    M must be positive.
        //
        //    Input, int N, the number of columns of the matrix.
        //    N must be positive.
        //
        //    Input, double A[M*N], the R8GE matrix.
        //
        //    Input, double X[M], the vector to be multiplied by A.
        //
        //    Output, double R8GE_MTV[N], the product A' * x.
        //
    {
        int i;

        double[] b = r8vec_zeros_new(n);

        for (i = 0; i < n; i++)
        {
            int j;
            for (j = 0; j < m; j++)
            {
                b[i] += a[j + i * m] * x[j];
            }
        }

        return b;
    }

    public static double[] r8ge_mu(int m, int n, double[] a_lu, char trans, int[] pivot,
            double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_MU computes A * x or A' * x, using R8GE_TRF factors.
        //
        //  Discussion:
        //
        //    The R8GE storage format is used for a "general" M by N matrix.  
        //    A physical storage space is made for each logical entry.  The two 
        //    dimensional logical array is mapped to a vector, in which storage is 
        //    by columns.
        //
        //    It is assumed that R8GE_TRF has overwritten the original matrix
        //    information by PLU factors.  R8GE_MU is able to reconstruct the
        //    original matrix from the PLU factor data.
        //
        //    R8GE_MU allows the user to check that the solution of a linear
        //    system is correct, without having to save an unfactored copy
        //    of the matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    01 September 2006
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
        //    Input, int M, the number of rows in the matrix.
        //
        //    Input, int N, the number of columns in the matrix.
        //
        //    Input, double A_LU[M*N], the LU factors from R8GE_TRF.
        //
        //    Input, char TRANS, specifies the form of the system of equations:
        //    'N':  A * x = b  (No transpose)
        //    'T':  A'* X = B  (Transpose)
        //    'C':  A'* X = B  (Conjugate transpose = Transpose)
        //
        //    Input, int PIVOT[*], the pivot vector computed by R8GE_TRF.
        //
        //    Input, double X[*], the vector to be multiplied.
        //    For the untransposed case, X should have N entries.
        //    For the transposed case, X should have M entries.
        //
        //    Output, double R8GE_MU[*], the result of the multiplication.
        //    For the untransposed case, the result should have M entries.
        //    For the transposed case, the result should have N entries.
        //
    {
        double[] b;
        int i;
        int j;
        int k;
        double temp;

        int npiv = Math.Min(m - 1, n);
        int mn_max = Math.Max(m, n);
        double[] y = new double[mn_max];

        switch (trans)
        {
            case 'n':
            case 'N':
            {
                b = new double[m];
                //
                //  Y[MN] = U[MNxN] * X[N].
                //
                for (i = 0; i < n; i++)
                {
                    y[i] = 0.0;
                }

                for (j = 1; j <= n; j++)
                {
                    for (i = 1; i <= Math.Min(j, m); i++)
                    {
                        y[i - 1] += a_lu[i - 1 + (j - 1) * m] * x[j - 1];
                    }
                }

                //
                //  Z[M] = L[MxMN] * Y[MN] = L[MxMN] * U[MNxN] * X[N].
                //
                for (i = 0; i < m; i++)
                {
                    if (i < n)
                    {
                        b[i] = y[i];
                    }
                    else
                    {
                        b[i] = 0.0;
                    }
                }

                for (j = Math.Min(m - 1, n); 1 <= j; j--)
                {
                    for (i = j + 1; i <= m; i++)
                    {
                        b[i - 1] += a_lu[i - 1 + (j - 1) * m] * y[j - 1];
                    }
                }

                //
                //  B = P * Z = P * L * Y = P * L * U * X = A * x.
                //
                for (j = npiv; 1 <= j; j--)
                {
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
            case 't':
            case 'T':
            case 'c':
            case 'C':
            {
                b = new double[n];
                //
                //  Y = P' * X:
                //
                for (i = 1; i <= npiv; i++)
                {
                    k = pivot[i - 1];

                    if (k != i)
                    {
                        temp = x[k - 1];
                        x[k - 1] = x[i - 1];
                        x[i - 1] = temp;
                    }
                }

                for (i = 0; i < n; i++)
                {
                    if (i < m)
                    {
                        b[i] = x[i];
                    }
                    else
                    {
                        b[i] = 0.0;
                    }
                }

                //
                //  Z = L' * Y:
                //
                for (j = 1; j <= Math.Min(m - 1, n); j++)
                {
                    for (i = j + 1; i <= m; i++)
                    {
                        b[j - 1] += x[i - 1] * a_lu[i - 1 + (j - 1) * m];
                    }
                }

                //
                //  B = U' * Z.
                //
                for (i = m; 1 <= i; i--)
                {
                    for (j = i + 1; j <= n; j++)
                    {
                        b[j - 1] += b[i - 1] * a_lu[i - 1 + (j - 1) * m];
                    }

                    if (i <= n)
                    {
                        b[i - 1] *= a_lu[i - 1 + (i - 1) * m];
                    }
                }

                //
                //  Now restore X.
                //
                for (i = npiv; 1 <= i; i--)
                {
                    k = pivot[i - 1];

                    if (k != i)
                    {
                        temp = x[k - 1];
                        x[k - 1] = x[i - 1];
                        x[i - 1] = temp;
                    }
                }

                break;
            }
            //
            default:
                Console.WriteLine("");
                Console.WriteLine("R8GE_MU - Fatal error!");
                Console.WriteLine("  Illegal value of TRANS = \"" + trans + "\"");
                return null;
        }


        return b;
    }

    public static double[] r8ge_mv(int m, int n, double[] a, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_MV multiplies an R8GE matrix by an R8VEC.
        //
        //  Discussion:
        //
        //    The R8GE storage format is used for a general M by N matrix.  A storage 
        //    space is made for each entry.  The two dimensional logical
        //    array can be thought of as a vector of M*N entries, starting with
        //    the M entries in the column 1, then the M entries in column 2
        //    and so on.  Considered as a vector, the entry A(I,J) is then stored
        //    in vector location I+(J-1)*M.
        //
        //    R8GE storage is used by LINPACK and LAPACK.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 July 2014
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
        //    Input, double A[M*N], the matrix.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Output, double R8GE_MV[M], the product A * x.
        //
    {
        int i;

        double[] b = new double[m];

        for (i = 0; i < m; i++)
        {
            b[i] = 0.0;
            int j;
            for (j = 0; j < n; j++)
            {
                b[i] += a[i + j * m] * x[j];
            }
        }

        return b;
    }

    public static double[] r8ge_random(int m, int n, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_RANDOM randomizes an R8GE matrix.
        //
        //  Discussion:
        //
        //    The R8GE storage format is used for a "general" M by N matrix.  
        //    A physical storage space is made for each logical entry.  The two 
        //    dimensional logical array is mapped to a vector, in which storage is 
        //    by columns.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 January 2004
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
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double R8GE_RANDOM[M*N], the randomized M by N matrix, 
        //    with entries between 0 and 1.
        //
    {
        int j;

        double[] a = new double[m * n];

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                a[i + j * m] = UniformRNG.r8_uniform_01(ref seed);
            }
        }

        return a;
    }

    public static double[] r8ge_res(int m, int n, double[] a, double[] x, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_RES computes the residual R = B-A*X for R8GE matrices.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    05 June 2014
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
        //    Input, double A[M*N], the matrix.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Input, double B[M], the desired result A * x.
        //
        //    Output, double R8GE_RES[M], the residual R = B - A * X.
        //
    {
        int i;

        double[] r = r8ge_mv(m, n, a, x);
        for (i = 0; i < m; i++)
        {
            r[i] = b[i] - r[i];
        }

        return r;
    }

    public static double[] r8ge_ml(int n, double[] a_lu, int[] pivot, double[] x, int job)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_ML computes A * x or A' * x, using R8GE_FA factors.
        //
        //  Discussion:
        //
        //    The R8GE storage format is used for a "general" M by N matrix.  
        //    A physical storage space is made for each logical entry.  The two 
        //    dimensional logical array is mapped to a vector, in which storage is 
        //    by columns.
        //
        //    It is assumed that R8GE_FA has overwritten the original matrix
        //    information by LU factors.  R8GE_ML is able to reconstruct the
        //    original matrix from the LU factor data.
        //
        //    R8GE_ML allows the user to check that the solution of a linear
        //    system is correct, without having to save an unfactored copy
        //    of the matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 September 2003
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
        //    Input, double A_LU[N*N], the LU factors from R8GE_FA.
        //
        //    Input, int PIVOT[N], the pivot vector computed by R8GE_FA.
        //
        //    Input, double X[N], the vector to be multiplied.
        //
        //    Input, int JOB, specifies the operation to be done:
        //    JOB = 0, compute A * x.
        //    JOB nonzero, compute A' * X.
        //
        //    Output, double R8GE_ML[N], the result of the multiplication.
        //
    {
        int i;
        int j;
        int k;
        double temp;
        //
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
                    for (i = 1; i <= j - 1; i++)
                    {
                        b[i - 1] += a_lu[i - 1 + (j - 1) * n] * b[j - 1];
                    }

                    b[j - 1] = a_lu[j - 1 + (j - 1) * n] * b[j - 1];
                }

                //
                //  B = PL * Y = PL * U * X = A * x.
                //
                for (j = n - 1; 1 <= j; j--)
                {
                    for (i = j + 1; i <= n; i++)
                    {
                        b[i - 1] -= a_lu[i - 1 + (j - 1) * n] * b[j - 1];
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
                //  Y = (PL)' * X:
                //
                for (j = 1; j <= n - 1; j++)
                {
                    k = pivot[j - 1];

                    if (k != j)
                    {
                        temp = b[k - 1];
                        b[k - 1] = b[j - 1];
                        b[j - 1] = temp;
                    }

                    temp = 0.0;
                    for (i = j + 1; i <= n; i++)
                    {
                        temp += b[i - 1] * a_lu[i - 1 + (j - 1) * n];
                    }

                    b[j - 1] -= temp;

                }

                //
                //  B = U' * Y = ( PL * U )' * X = A' * X.
                //
                for (i = n; 1 <= i; i--)
                {
                    for (j = i + 1; j <= n; j++)
                    {
                        b[j - 1] += b[i - 1] * a_lu[i - 1 + (j - 1) * n];
                    }

                    b[i - 1] *= a_lu[i - 1 + (i - 1) * n];
                }

                break;
            }
        }

        return b;
    }

    public static double[] r8ge_mm_new(int n1, int n2, int n3, double[] a, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_MM_NEW multiplies two matrices.
        //
        //  Discussion:
        //
        //    An R8GE is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //    For this routine, the result is returned as the function value.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 December 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N1, N2, N3, the order of the matrices.
        //
        //    Input, double A[N1*N2], double B[N2*N3], the matrices to multiply.
        //
        //    Output, double R8GE_MM_NEW[N1*N3], the product matrix C = A * B.
        //
    {
        int i;

        double[] c = new double[n1 * n3];

        for (i = 0; i < n1; i++)
        {
            int j;
            for (j = 0; j < n3; j++)
            {
                c[i + j * n1] = 0.0;
                int k;
                for (k = 0; k < n2; k++)
                {
                    c[i + j * n1] += a[i + k * n1] * b[k + j * n2];
                }
            }
        }

        return c;
    }

    public static double[] r8ge_mtm(int n, double[] a, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_MTM computes C=A'*B for R8GE matrices.
        //
        //  Discussion:
        //
        //    The R8GE storage format is used for a "general" M by N matrix.  
        //    A physical storage space is made for each logical entry.  The two 
        //    dimensional logical array is mapped to a vector, in which storage is 
        //    by columns.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    04 August 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrices.
        //    N must be positive.
        //
        //    Input, double A[N*N], B[N*N], the factors.
        //
        //    Output, double C[N*N], the product.
        //
    {
        double[] c = new double[n * n];

        for (int j = 0; j < n; j++)
        {
            for (int i = 0; i < n; i++)
            {
                c[i + j * n] = 0.0;
                for (int k = 0; k < n; k++)
                {
                    c[i + j * n] += a[k + i * n] * b[k + j * n];
                }
            }
        }

        return c;
    }

    public static void r8ge_plu(int m, int n, double[] a, ref double[] p, ref double[] l, ref double[] u)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_PLU produces the PLU factors of an R8GE matrix.
        //
        //  Discussion:
        //
        //    The R8GE storage format is used for a "general" M by N matrix.  
        //    A physical storage space is made for each logical entry.  The two 
        //    dimensional logical array is mapped to a vector, in which storage is 
        //    by columns.
        //
        //    The PLU factors of the M by N matrix A are:
        //
        //      P, an M by M permutation matrix P,
        //      L, an M by M unit lower triangular matrix,
        //      U, an M by N upper triangular matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 December 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows in A.
        //
        //    Input, int N, the number of columns in A.
        //
        //    Input, double A[M,N], the M by N matrix to be factored.
        //
        //    Output, double P[M*M], the M by M permutation factor.
        //
        //    Output, double L[M*M], the M by M unit lower triangular factor.
        //
        //    Output, double U[M*N], the M by N upper triangular factor.
        //
    {
        int j;
        //
        //  Initialize:
        //
        //    P: = M by M Identity
        //    L: = M by M Identity
        //    U: = A
        //

        r8ge_identity(m, m, ref p);
        r8ge_identity(m, m, ref l);
        r8ge_copy(m, n, a, ref u);
        //
        //  On step J, find the pivot row and the pivot value.
        //
        for (j = 0; j < Math.Min(m - 2, n - 1); j++)
        {
            double pivot_value = 0.0;
            int pivot_row = -1;

            int i;
            for (i = j; i < m; i++)
            {
                if (pivot_value < Math.Abs(u[i + j * m]))
                {
                    pivot_value = Math.Abs(u[i + j * m]);
                    pivot_row = i;
                }
            }

            //
            //  If the pivot row is nonzero swap:
            //  * rows J and PIVOT_ROW of U;
            //  * rows J and PIVOT_ROW of L and cols J and PIVOT_ROW of L;
            //  * cols J and PIVOT_ROW of P.
            //
            if (pivot_row == -1)
            {
                continue;
            }

            double t;
            int k;
            for (k = 0; k < n; k++)
            {
                t = u[j + k * m];
                u[j + k * m] = u[pivot_row + k * m];
                u[pivot_row + k * m] = t;
            }

            for (k = 0; k < m; k++)
            {
                t = l[j + k * m];
                l[j + k * m] = l[pivot_row + k * m];
                l[pivot_row + k * m] = t;
            }

            for (k = 0; k < m; k++)
            {
                t = l[k + j * m];
                l[k + j * m] = l[k + pivot_row * m];
                l[k + pivot_row * m] = t;
            }

            for (k = 0; k < m; k++)
            {
                t = p[k + j * m];
                p[k + j * m] = p[k + pivot_row * m];
                p[k + pivot_row * m] = t;
            }

            //
            //  Zero out the entries in column J, from row J+1 to M.
            //
            for (i = j + 1; i < m; i++)
            {
                if (u[i + j * m] == 0.0)
                {
                    continue;
                }

                l[i + j * m] = u[i + j * m] / u[j + j * m];
                u[i + j * m] = 0.0;
                for (k = j + 1; k < n; k++)
                {
                    u[i + k * m] -= l[i + j * m] * u[j + k * m];
                }
            }
        }

    }

    public static double[] r8ge_poly(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_POLY computes the characteristic polynomial of an R8GE matrix.
        //
        //  Discussion:
        //
        //    The R8GE storage format is used for a "general" M by N matrix.  
        //    A physical storage space is made for each logical entry.  The two 
        //    dimensional logical array is mapped to a vector, in which storage is 
        //    by columns.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 November 2003
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
        //    Input, double A[N*N], the R8GE matrix.
        //
        //    Output, double R8GE_POLY[N+1], the coefficients of the characteristic
        //    polynomial of A.  P(I) contains the coefficient of X**I.
        //
    {
        int order;

        double[] p = new double[n + 1];
        double[] work2 = new double[n * n];
        //
        //  Initialize WORK1 to the identity matrix.
        //
        double[] work1 = r8ge_identity_new(n, n);

        p[n] = 1.0;

        for (order = n - 1; 0 <= order; order--)
        {
            //
            //  Work2 = A * WORK1.
            //
            int i;
            int j;
            for (j = 0; j < n; j++)
            {
                for (i = 0; i < n; i++)
                {
                    work2[i + j * n] = 0.0;
                    int k;
                    for (k = 0; k < n; k++)
                    {
                        work2[i + j * n] += a[i + k * n] * work1[k + j * n];
                    }
                }
            }

            //
            //  Take the trace.
            //
            double trace = 0.0;
            for (i = 0; i < n; i++)
            {
                trace += work2[i + i * n];
            }

            //
            //  P(ORDER) = - Trace ( WORK2 ) / ( N - ORDER )
            //
            p[order] = -trace / (n - order);
            //
            //  WORK1 := WORK2 + P(ORDER) * Identity.
            //
            for (j = 0; j < n; j++)
            {
                for (i = 0; i < n; i++)
                {
                    work1[i + j * n] = work2[i + j * n];
                }
            }

            for (j = 0; j < n; j++)
            {
                work1[j + j * n] += p[order];
            }
        }

        return p;
    }

    public static void r8ge_print(int m, int n, double[] a, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_PRINT prints an R8GE matrix.
        //
        //  Discussion:
        //
        //    The R8GE storage format is used for a "general" M by N matrix.  
        //    A physical storage space is made for each logical entry.  The two 
        //    dimensional logical array is mapped to a vector, in which storage is 
        //    by columns.
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
        //    Input, double A[M*N], the R8GE matrix.
        //
        //    Input, string TITLE, a title.
        //
    {
        r8ge_print_some(m, n, a, 1, 1, m, n, title);
    }

    public static void r8ge_print_some(int m, int n, double[] a, int ilo, int jlo, int ihi,
            int jhi, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_PRINT_SOME prints some of an R8GE matrix.
        //
        //  Discussion:
        //
        //    The R8GE storage format is used for a "general" M by N matrix.  
        //    A physical storage space is made for each logical entry.  The two 
        //    dimensional logical array is mapped to a vector, in which storage is 
        //    by columns.
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
        //    Input, double A[M*N], the R8GE matrix.
        //
        //    Input, int ILO, JLO, IHI, JHI, designate the first row and
        //    column, and the last row and column to be printed.
        //
        //    Input, string TITLE, a title.
        //
    {
        const int INCX = 5;

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
            //
            //  For each column J in the current range...
            //
            //  Write the header.
            //
            string cout = "  Col:    ";
            for (int j = j2lo; j <= j2hi; j++)
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
            int i2hi = Math.Min(ihi, m);

            for (int i = i2lo; i <= i2hi; i++)
            {
                //
                //  Print out (up to) 5 entries in row I, that lie in the current strip.
                //
                cout = i.ToString().PadLeft(5) + "  ";
                for (int j = j2lo; j <= j2hi; j++)
                {
                    cout += a[i - 1 + (j - 1) * m].ToString().PadLeft(12) + "  ";
                }

                Console.WriteLine(cout);
            }
        }
    }

    public static double[] r8ge_dilu(int m, int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_DILU produces the diagonal incomplete LU factor of an R8GE matrix.
        //
        //  Discussion:
        //
        //    The R8GE storage format is used for a "general" M by N matrix.  
        //    A physical storage space is made for each logical entry.  The two 
        //    dimensional logical array is mapped to a vector, in which storage is 
        //    by columns.
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
        //  Parameters:
        //
        //    Input, int M, the number of rows in A.
        //
        //    Input, int N, the number of columns in A.
        //
        //    Input, double A[M*N], the R8GE matrix.
        //
        //    Output, double R8GE_DILU[M], the D-ILU factor.
        //
    {
        int i;

        double[] d = new double[m];

        for (i = 0; i < m; i++)
        {
            if (i < n)
            {
                d[i] = a[i + i * m];
            }
            else
            {
                d[i] = 0.0;
            }
        }

        for (i = 0; i < m && i < n; i++)
        {
            d[i] = 1.0 / d[i];
            int j;
            for (j = i + 1; j < m && j < n; j++)
            {
                d[j] -= a[j + i * m] * d[i] * a[i + j * m];
            }
        }

        return d;
    }

    public static int r8ge_fa(int n, ref double[] a, ref int[] pivot, int aIndex = 0, int pivotIndex = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_FA performs a LINPACK-style PLU factorization of a R8GE matrix.
        //
        //  Discussion:
        //
        //    The R8GE storage format is used for a "general" M by N matrix.  
        //    A physical storage space is made for each logical entry.  The two 
        //    dimensional logical array is mapped to a vector, in which storage is 
        //    by columns.
        //
        //    R8GE_FA is a simplified version of the LINPACK routine SGEFA.
        //
        //    The two dimensional array is stored by columns in a one dimensional
        //    array.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 September 2003
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
        //    N must be positive.
        //
        //    Input/output, double A[N*N], the matrix to be factored.
        //    On output, A contains an upper triangular matrix and the multipliers
        //    which were used to obtain it.  The factorization can be written
        //    A = L * U, where L is a product of permutation and unit lower
        //    triangular matrices and U is upper triangular.
        //
        //    Output, int PIVOT[N], a vector of pivot indices.
        //
        //    Output, int R8GE_FA, singularity flag.
        //    0, no singularity detected.
        //    nonzero, the factorization failed on the INFO-th step.
        //
    {
        int k;
        //
        for (k = 1; k <= n - 1; k++)
        {
            //
            //  Find L, the index of the pivot row.
            //
            int l = k;

            int i;
            for (i = k + 1; i <= n; i++)
            {
                if (Math.Abs(a[(l - 1 + (k - 1) * n + aIndex) % a.Length]) < Math.Abs(a[(i - 1 + (k - 1) * n + aIndex) % a.Length]))
                {
                    l = i;
                }
            }

            pivot[(k - 1 + pivotIndex) % pivot.Length] = l;
            switch (a[(l - 1 + (k - 1) * n + aIndex) % a.Length])
            {
                //
                //  If the pivot index is zero, the algorithm has failed.
                //
                case 0.0:
                    Console.WriteLine("");
                    Console.WriteLine("R8GE_FA - Fatal error!");
                    Console.WriteLine("  Zero pivot on step " + k + "");
                    return 0;
            }

            //
            //  Interchange rows L and K if necessary.
            //
            double t;
            if (l != k)
            {
                t = a[(l - 1 + (k - 1) * n + aIndex) % a.Length];
                a[(l - 1 + (k - 1) * n + aIndex) % a.Length] = a[(k - 1 + (k - 1) * n + aIndex) % a.Length];
                a[(k - 1 + (k - 1) * n + aIndex) % a.Length] = t;
            }

            //
            //  Normalize the values that lie below the pivot entry A(K,K).
            //
            for (i = k + 1; i <= n; i++)
            {
                a[(i - 1 + (k - 1) * n + aIndex) % a.Length] = -a[(i - 1 + (k - 1) * n + aIndex) % a.Length] / a[(k - 1 + (k - 1) * n + aIndex) % a.Length];
            }

            //
            //  Row elimination with column indexing.
            //
            int j;
            for (j = k + 1; j <= n; j++)
            {
                if (l != k)
                {
                    t = a[(l - 1 + (j - 1) * n + aIndex) % a.Length];
                    a[(l - 1 + (j - 1) * n + aIndex) % a.Length] = a[(k - 1 + (j - 1) * n + aIndex) % a.Length];
                    a[(k - 1 + (j - 1) * n + aIndex) % a.Length] = t;
                }

                for (i = k + 1; i <= n; i++)
                {
                    a[(i - 1 + (j - 1) * n + aIndex) % a.Length] += a[(i - 1 + (k - 1) * n + aIndex) % a.Length] * a[(k - 1 + (j - 1) * n + aIndex) % a.Length];
                }

            }

        }

        pivot[(n - 1 + pivotIndex) % pivot.Length] = n;

        switch (a[(n - 1 + (n - 1) * n + aIndex) % a.Length])
        {
            case 0.0:
                Console.WriteLine("");
                Console.WriteLine("R8GE_FA - Fatal error!");
                Console.WriteLine("  Zero pivot on step " + n + "");
                return 0;
            default:
                return 0;
        }
    }

    public static void r8ge_sl(int n, double[] a_lu, int[] pivot, ref double[] x, int job)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_SL solves a R8GE system factored by R8GE_FA.
        //
        //  Discussion:
        //
        //    The R8GE storage format is used for a "general" M by N matrix.  
        //    A physical storage space is made for each logical entry.  The two 
        //    dimensional logical array is mapped to a vector, in which storage is 
        //    by columns.
        //
        //    R8GE_SL is a simplified version of the LINPACK routine SGESL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 April 2012
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
        //    Input, double A_LU[N*N], the LU factors from R8GE_FA.
        //
        //    Input, int PIVOT[N], the pivot vector from R8GE_FA.
        //
        //    Input/output, double X[N], on input, the right hand side vector.
        //    On output, the solution vector.
        //
        //    Input, int JOB, specifies the operation.
        //    0, solve A * x = b.
        //    nonzero, solve A' * x = b.
        //
    {
        int i;
        int k;
        int l;
        double t;
        switch (job)
        {
            //
            //  Solve A * x = b.
            //
            case 0:
            {
                //
                //  Solve PL * Y = B.
                //
                for (k = 1; k <= n - 1; k++)
                {
                    l = pivot[k - 1];

                    if (l != k)
                    {
                        t = x[l - 1];
                        x[l - 1] = x[k - 1];
                        x[k - 1] = t;
                    }

                    for (i = k + 1; i <= n; i++)
                    {
                        x[i - 1] += a_lu[i - 1 + (k - 1) * n] * x[k - 1];
                    }
                }

                //
                //  Solve U * X = Y.
                //
                for (k = n; 1 <= k; k--)
                {
                    x[k - 1] /= a_lu[k - 1 + (k - 1) * n];
                    for (i = 1; i <= k - 1; i++)
                    {
                        x[i - 1] -= a_lu[i - 1 + (k - 1) * n] * x[k - 1];
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
                    t = 0.0;
                    for (i = 1; i <= k - 1; i++)
                    {
                        t += x[i - 1] * a_lu[i - 1 + (k - 1) * n];
                    }

                    x[k - 1] = (x[k - 1] - t) / a_lu[k - 1 + (k - 1) * n];
                }

                //
                //  Solve ( PL )' * X = Y.
                //
                for (k = n - 1; 1 <= k; k--)
                {
                    t = 0.0;
                    for (i = k + 1; i <= n; i++)
                    {
                        t += x[i - 1] * a_lu[i - 1 + (k - 1) * n];
                    }

                    x[k - 1] += t;

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
    }

    public static double[] r8ge_fss_new(int n, ref double[] a, int nb, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_FSS_NEW factors and solves multiple R8GE systems.
        //
        //  Discussion:
        //
        //    The R8GE storage format is used for a "general" M by N matrix.  
        //    A physical storage space is made for each logical entry.  The two 
        //    dimensional logical array is mapped to a vector, in which storage is 
        //    by columns.
        //
        //    This routine does not save the LU factors of the matrix, and hence cannot
        //    be used to efficiently solve multiple linear systems, or even to
        //    factor A at one time, and solve a single linear system at a later time.
        //
        //    This routine uses partial pivoting, but no pivot vector is required.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 June 2009
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
        //    Input/output, double A[N*N].
        //    On input, A is the coefficient matrix of the linear system.
        //    On output, A is in unit upper triangular form, and
        //    represents the U factor of an LU factorization of the
        //    original coefficient matrix.
        //
        //    Input, int NB, the number of right hand sides.
        //
        //    Input, double B[N*NB], the right hand sides of the linear systems.
        //
        //    Output, double R8GE_FSS_NEW[N*NB], the solutions of the linear systems.
        //
    {
        int i;
        int j;
        int jcol;

        double[] x = new double[n * nb];

        for (j = 0; j < nb; j++)
        {
            for (i = 0; i < n; i++)
            {
                x[i + j * n] = b[i + j * n];
            }
        }

        for (jcol = 1; jcol <= n; jcol++)
        {
            //
            //  Find the maximum element in column I.
            //
            double piv = Math.Abs(a[jcol - 1 + (jcol - 1) * n]);
            int ipiv = jcol;
            for (i = jcol + 1; i <= n; i++)
            {
                if (piv < Math.Abs(a[i - 1 + (jcol - 1) * n]))
                {
                    piv = Math.Abs(a[i - 1 + (jcol - 1) * n]);
                    ipiv = i;
                }
            }

            switch (piv)
            {
                case 0.0:
                    Console.WriteLine("");
                    Console.WriteLine("R8GE_FSS_NEW - Fatal error!");
                    Console.WriteLine("  Zero pivot on step " + jcol + "");
                    return null;
            }

            //
            //  Switch rows JCOL and IPIV, and X.
            //
            double t;
            if (jcol != ipiv)
            {
                for (j = 1; j <= n; j++)
                {
                    t = a[jcol - 1 + (j - 1) * n];
                    a[jcol - 1 + (j - 1) * n] = a[ipiv - 1 + (j - 1) * n];
                    a[ipiv - 1 + (j - 1) * n] = t;
                }

                for (j = 0; j < nb; j++)
                {
                    t = x[jcol - 1 + j * n];
                    x[jcol - 1 + j * n] = x[ipiv - 1 + j * n];
                    x[ipiv - 1 + j * n] = t;
                }
            }

            //
            //  Scale the pivot row.
            //
            t = a[jcol - 1 + (jcol - 1) * n];
            a[jcol - 1 + (jcol - 1) * n] = 1.0;
            for (j = jcol + 1; j <= n; j++)
            {
                a[jcol - 1 + (j - 1) * n] /= t;
            }

            for (j = 0; j < nb; j++)
            {
                x[jcol - 1 + j * n] /= t;
            }

            //
            //  Use the pivot row to eliminate lower entries in that column.
            //
            for (i = jcol + 1; i <= n; i++)
            {
                if (a[i - 1 + (jcol - 1) * n] != 0.0)
                {
                    t = -a[i - 1 + (jcol - 1) * n];
                    a[i - 1 + (jcol - 1) * n] = 0.0;
                    for (j = jcol + 1; j <= n; j++)
                    {
                        a[i - 1 + (j - 1) * n] += t * a[jcol - 1 + (j - 1) * n];
                    }

                    for (j = 0; j < nb; j++)
                    {
                        x[i - 1 + j * n] += t * x[jcol - 1 + j * n];
                    }
                }
            }
        }

        //
        //  Back solve.
        //
        for (jcol = n; 2 <= jcol; jcol--)
        {
            for (i = 1; i < jcol; i++)
            {
                for (j = 0; j < nb; j++)
                {
                    x[i - 1 + j * n] -= a[i - 1 + (jcol - 1) * n] * x[jcol - 1 + j * n];
                }
            }
        }

        return x;
    }

    public static double[] r8ge_hilbert(int m, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_HILBERT returns the Hilbert matrix.
        //
        //  Formula:
        //
        //    A(I,J) = 1 / ( I + J - 1 )
        //
        //  Example:
        //
        //    N = 5
        //
        //    1/1 1/2 1/3 1/4 1/5
        //    1/2 1/3 1/4 1/5 1/6
        //    1/3 1/4 1/5 1/6 1/7
        //    1/4 1/5 1/6 1/7 1/8
        //    1/5 1/6 1/7 1/8 1/9
        //
        //  Properties:
        //
        //    A is a Hankel matrix: constant along anti-diagonals.
        //
        //    A is positive definite.
        //
        //    A is symmetric: A' = A.
        //
        //    Because A is symmetric, it is normal.
        //
        //    Because A is normal, it is diagonalizable.
        //
        //    A is totally positive.
        //
        //    A is a Cauchy matrix.
        //
        //    A is nonsingular.
        //
        //    A is very ill-conditioned.
        //
        //    The entries of the inverse of A are all integers.
        //
        //    The sum of the entries of the inverse of A is N*N.
        //
        //    The ratio of the absolute values of the maximum and minimum
        //    eigenvalues is roughly EXP(3.5*N).
        //
        //    The determinant of the Hilbert matrix of order 10 is
        //    2.16417... * 10^(-53).
        //
        //    If the (1,1) entry of the 5 by 5 Hilbert matrix is changed
        //    from 1 to 24/25, the matrix is exactly singular.  And there
        //    is a similar rule for larger Hilbert matrices.
        //
        //    The family of matrices is nested as a function of N.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 July 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    MD Choi,
        //    Tricks or treats with the Hilbert matrix,
        //    American Mathematical Monthly,
        //    Volume 90, 1983, pages 301-312.
        //
        //    Robert Gregory, David Karney,
        //    A Collection of Matrices for Testing Computational Algorithms,
        //    Wiley, 1969,
        //    ISBN: 0882756494,
        //    LC: QA263.68
        //
        //    Nicholas Higham,
        //    Accuracy and Stability of Numerical Algorithms,
        //    Society for Industrial and Applied Mathematics, Philadelphia, PA,
        //    USA, 1996; section 26.1.
        //
        //    Donald Knuth,
        //    The Art of Computer Programming,
        //    Volume 1, Fundamental Algorithms, Second Edition
        //    Addison-Wesley, Reading, Massachusetts, 1973, page 37.
        //
        //    Morris Newman, John Todd,
        //    Example A13,
        //    The evaluation of matrix inversion programs,
        //    Journal of the Society for Industrial and Applied Mathematics,
        //    Volume 6, 1958, pages 466-476.
        //
        //    Joan Westlake,
        //    A Handbook of Numerical Matrix Inversion and Solution of 
        //    Linear Equations,
        //    John Wiley, 1968,
        //    ISBN13: 978-0471936756,
        //    LC: QA263.W47.
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of the matrix.
        //
        //    Output, double R8GE_HILBERT[M*N], the matrix.
        //
    {
        int j;

        double[] a = new double[m * n];

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                a[i + j * m] = 1.0 / (i + j + 1);
            }
        }

        return a;
    }

    public static double[] r8ge_hilbert_inverse(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_HILBERT_INVERSE returns the inverse of the HILBERT matrix.
        //
        //  Formula:
        //
        //    A(I,J) =  (-1)^(I+J) * (N+I-1)! * (N+J-1)! /
        //           [ (I+J-1) * ((I-1)!*(J-1)!)^2 * (N-I)! * (N-J)! ]
        //
        //  Example:
        //
        //    N = 5
        //
        //       25    -300     1050    -1400     630
        //     -300    4800   -18900    26880  -12600
        //     1050  -18900    79380  -117600   56700
        //    -1400   26880  -117600   179200  -88200
        //      630  -12600    56700   -88200   44100
        //
        //  Properties:
        //
        //    A is symmetric: A' = A.
        //
        //    Because A is symmetric, it is normal.
        //
        //    Because A is normal, it is diagonalizable.
        //
        //    A is almost impossible to compute accurately by general routines
        //    that compute the inverse.
        //
        //    A is the exact inverse of the Hilbert matrix; however, if the
        //    Hilbert matrix is stored on a finite precision computer, and
        //    hence rounded, A is actually a poor approximation
        //    to the inverse of that rounded matrix.  Even though Gauss elimination
        //    has difficulty with the Hilbert matrix, it can compute an approximate
        //    inverse matrix whose residual is lower than that of the
        //    "exact" inverse.
        //
        //    All entries of A are integers.
        //
        //    The sum of the entries of A is N^2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 June 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Output, double R8GE_HILBERT_INVERSE[N*N], the matrix.
        //
    {
        int j;

        double[] a = new double[n * n];
        //
        //  Set the (1,1) entry.
        //
        a[0 + 0 * n] = n * n;
        //
        //  Define Row 1, Column J by recursion on Row 1 Column J-1
        //
        int i = 0;
        for (j = 1; j < n; j++)
        {
            a[i + j * n] = -a[i + (j - 1) * n]
                           * ((n + j) * (i + j) * (n - j))
                           / ((i + j + 1) * j * j);
        }

        //
        //  Define Row I by recursion on row I-1
        //
        for (i = 1; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                a[i + j * n] = -a[i - 1 + j * n]
                               * ((n + i) * (i + j) * (n - i))
                               / ((i + j + 1) * i * i);
            }
        }

        return a;
    }

    public static void r8ge_identity(int m, int n, ref double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_IDENTITY sets an R8GE matrix to the identity.
        //
        //  Discussion:
        //
        //    The R8GE storage format is used for a "general" M by N matrix.  
        //    A physical storage space is made for each logical entry.  The two 
        //    dimensional logical array is mapped to a vector, in which storage is 
        //    by columns.
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
        //    Input, int M, N, the order of A.
        //
        //    Output, double A[M*N], the identity matrix.
        //
    {
        int j;

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                if (i == j)
                {
                    a[i + j * m] = 1.0;
                }
                else
                {
                    a[i + j * m] = 0.0;
                }
            }
        }

    }

    public static double[] r8ge_identity_new(int m, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_IDENTITY_NEW sets an R8GE matrix to the identity.
        //
        //  Discussion:
        //
        //    The R8GE storage format is used for a "general" M by N matrix.  
        //    A physical storage space is made for each logical entry.  The two 
        //    dimensional logical array is mapped to a vector, in which storage is 
        //    by columns.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 february 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of A.
        //
        //    Output, double R8GE_IDENTITY_NEW[M*N], the identity matrix.
        //
    {
        int j;

        double[] a = new double[m * n];

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                if (i == j)
                {
                    a[i + j * m] = 1.0;
                }
                else
                {
                    a[i + j * m] = 0.0;
                }
            }
        }

        return a;
    }

    public static void r8ge_ilu(int m, int n, double[] a, ref double[] l, ref double[] u)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_ILU produces the incomplete LU factors of an R8GE matrix.
        //
        //  Discussion:
        //
        //    The R8GE storage format is used for a "general" M by N matrix.  
        //    A physical storage space is made for each logical entry.  The two 
        //    dimensional logical array is mapped to a vector, in which storage is 
        //    by columns.
        //
        //    The incomplete LU factors of the M by N matrix A are:
        //
        //      L, an M by M unit lower triangular matrix,
        //      U, an M by N upper triangular matrix
        //
        //    with the property that L and U are computed in the same way as
        //    the usual LU factors, except that, whenever an off diagonal element
        //    of the original matrix is zero, then the corresponding value of
        //    U is forced to be zero.
        //
        //    This condition means that it is no longer the case that A = L*U.
        //
        //    On the other hand, L and U will have a simple sparsity structure
        //    related to that of A.  The incomplete LU factorization is generally
        //    used as a preconditioner in iterative schemes applied to sparse
        //    matrices.  It is presented here merely for illustration.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    24 November 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows in A.
        //
        //    Input, int N, the number of columns in A.
        //
        //    Input, double A[M*N], the R8GE matrix.
        //
        //    Output, double L[M*M], the M by M unit lower triangular factor.
        //
        //    Output, double U[M*N], the M by N upper triangular factor.
        //
    {
        int j;
        //
        //  Initialize:
        //
        //    L := M by M Identity
        //    U := A
        //
        r8ge_identity(m, m, ref l);

        r8ge_copy(m, n, a, ref u);

        int jhi = Math.Min(m - 1, n);

        for (j = 0; j < jhi; j++)
        {
            //
            //  Zero out the entries in column J, from row J+1 to M.
            //
            int i;
            for (i = j + 1; i < m; i++)
            {
                if (u[i + j * m] == 0.0)
                {
                    continue;
                }

                l[i + j * m] = u[i + j * m] / u[j + j * m];
                u[i + j * m] = 0.0;

                int k;
                for (k = j + 1; k < n; k++)
                {
                    if (u[i + k * m] != 0.0)
                    {
                        u[i + k * m] -= l[i + j * m] * u[j + k * m];
                    }
                }
            }
        }

    }

    public static double[] r8ge_indicator(int m, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_INDICATOR sets up an R8GE indicator matrix.
        //
        //  Discussion:
        //
        //    The R8GE storage format is used for a "general" M by N matrix.  
        //    A physical storage space is made for each logical entry.  The two 
        //    dimensional logical array is mapped to a vector, in which storage is 
        //    by columns.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    25 January 2005
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
        //    Output, double R8GE_INDICATOR[M*N], the R8GE matrix.
        //
    {
        int i;

        double[] a = new double[m * n];

        int fac = (int) Math.Pow(10, (int) Math.Log10(n) + 1);

        for (i = 1; i <= m; i++)
        {
            int j;
            for (j = 1; j <= n; j++)
            {
                a[i - 1 + (j - 1) * m] = fac * i + j;
            }
        }

        return a;
    }


    public static double[] r8ge_inverse(int n, double[] a, int[] pivot)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_INVERSE computes the inverse of a R8GE matrix factored by R8GE_FA.
        //
        //  Discussion:
        //
        //    The R8GE storage format is used for a "general" M by N matrix.  
        //    A physical storage space is made for each logical entry.  The two 
        //    dimensional logical array is mapped to a vector, in which storage is 
        //    by columns.
        //
        //    R8GE_INVERSE is a simplified standalone version of the LINPACK routine
        //    SGEDI.
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
        //    Input, int N, the order of the matrix A.
        //
        //    Input, double A[N*N], the factor information computed by R8GE_FA.
        //
        //    Input, int PIVOT(N), the pivot vector from R8GE_FA.
        //
        //    Output, double R8GE_INVERSE[N*N], the inverse matrix.
        //
    {
        int i;
        int j;
        int k;

        double[] b = new double[n * n];
        //
        //  Compute Inverse(U).
        //
        for (k = 1; k <= n; k++)
        {
            for (i = 1; i <= k - 1; i++)
            {
                b[i - 1 + (k - 1) * n] = -b[i - 1 + (k - 1) * n] / a[k - 1 + (k - 1) * n];
            }

            b[k - 1 + (k - 1) * n] = 1.0 / a[k - 1 + (k - 1) * n];

            for (j = k + 1; j <= n; j++)
            {
                b[k - 1 + (j - 1) * n] = 0.0;
                for (i = 1; i <= k; i++)
                {
                    b[i - 1 + (j - 1) * n] += b[i - 1 + (k - 1) * n] * a[k - 1 + (j - 1) * n];
                }
            }
        }

        //
        //  Multiply Inverse(U) by Inverse(L).
        //
        for (k = n - 1; 1 <= k; k--)
        {
            for (i = k + 1; i <= n; i++)
            {
                b[i - 1 + (k - 1) * n] = 0.0;
            }

            for (j = k + 1; j <= n; j++)
            {
                for (i = 1; i <= n; i++)
                {
                    b[i - 1 + (k - 1) * n] += b[i - 1 + (j - 1) * n] * a[j - 1 + (k - 1) * n];
                }
            }

            if (pivot[k - 1] != k)
            {
                for (i = 1; i <= n; i++)
                {
                    (b[i - 1 + (k - 1) * n], b[i - 1 + (pivot[k - 1] - 1) * n]) = (b[i - 1 + (pivot[k - 1] - 1) * n], b[i - 1 + (k - 1) * n]);
                }

            }

        }

        return b;
    }

    public static void r8ge_fs(int n, ref double[] a, ref double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_FS factors and solves an R8GE system.
        //
        //  Discussion:
        //
        //    The R8GE storage format is used for a "general" M by N matrix.  
        //    A physical storage space is made for each logical entry.  The two 
        //    dimensional logical array is mapped to a vector, in which storage is 
        //    by columns.
        //
        //    The function does not save the LU factors of the matrix, and hence cannot
        //    be used to efficiently solve multiple linear systems, or even to
        //    factor A at one time, and solve a single linear system at a later time.
        //
        //    The function uses partial pivoting, but no pivot vector is required.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    04 December 2003
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
        //    Input/output, double A[N*N].
        //    On input, A is the coefficient matrix of the linear system.
        //    On output, A is in unit upper triangular form, and
        //    represents the U factor of an LU factorization of the
        //    original coefficient matrix.
        //
        //    Input/output, double X[N], on input, the right hand side of the linear system.
        //    On output, the solution of the linear system.
        //
    {
        int i;
        int jcol;

        for (jcol = 1; jcol <= n; jcol++)
        {
            //
            //  Find the maximum element in column I.
            //
            double piv = Math.Abs(a[jcol - 1 + (jcol - 1) * n]);
            int ipiv = jcol;
            for (i = jcol + 1; i <= n; i++)
            {
                if (piv < Math.Abs(a[i - 1 + (jcol - 1) * n]))
                {
                    piv = Math.Abs(a[i - 1 + (jcol - 1) * n]);
                    ipiv = i;
                }
            }

            switch (piv)
            {
                case 0.0:
                    Console.WriteLine("");
                    Console.WriteLine("R8GE_FS - Fatal error!");
                    Console.WriteLine("  Zero pivot on step " + jcol + "");
                    return;
            }

            //
            //  Switch rows JCOL and IPIV, and X.
            //
            double t;
            int j;
            if (jcol != ipiv)
            {
                for (j = 1; j <= n; j++)
                {
                    t = a[jcol - 1 + (j - 1) * n];
                    a[jcol - 1 + (j - 1) * n] = a[ipiv - 1 + (j - 1) * n];
                    a[ipiv - 1 + (j - 1) * n] = t;
                }

                t = x[jcol - 1];
                x[jcol - 1] = x[ipiv - 1];
                x[ipiv - 1] = t;
            }

            //
            //  Scale the pivot row.
            //
            t = a[jcol - 1 + (jcol - 1) * n];
            a[jcol - 1 + (jcol - 1) * n] = 1.0;
            for (j = jcol + 1; j <= n; j++)
            {
                a[jcol - 1 + (j - 1) * n] /= t;
            }

            x[jcol - 1] /= t;
            //
            //  Use the pivot row to eliminate lower entries in that column.
            //
            for (i = jcol + 1; i <= n; i++)
            {
                if (a[i - 1 + (jcol - 1) * n] != 0.0)
                {
                    t = -a[i - 1 + (jcol - 1) * n];
                    a[i - 1 + (jcol - 1) * n] = 0.0;
                    for (j = jcol + 1; j <= n; j++)
                    {
                        a[i - 1 + (j - 1) * n] += t * a[jcol - 1 + (j - 1) * n];
                    }

                    x[i - 1] += t * x[jcol - 1];
                }
            }
        }

        //
        //  Back solve.
        //
        for (jcol = n; 2 <= jcol; jcol--)
        {
            for (i = 1; i < jcol; i++)
            {
                x[i - 1] -= a[i - 1 + (jcol - 1) * n] * x[jcol - 1];
            }
        }
    }

    public static double[] r8ge_fs_new(int n, ref double[] a, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_FS_NEW factors and solves an R8GE system.
        //
        //  Discussion:
        //
        //    The R8GE storage format is used for a "general" M by N matrix.  
        //    A physical storage space is made for each logical entry.  The two 
        //    dimensional logical array is mapped to a vector, in which storage is 
        //    by columns.
        //
        //    The function does not save the LU factors of the matrix, and hence cannot
        //    be used to efficiently solve multiple linear systems, or even to
        //    factor A at one time, and solve a single linear system at a later time.
        //
        //    The function uses partial pivoting, but no pivot vector is required.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    04 December 2003
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
        //    Input/output, double A[N*N].
        //    On input, A is the coefficient matrix of the linear system.
        //    On output, A is in unit upper triangular form, and
        //    represents the U factor of an LU factorization of the
        //    original coefficient matrix.
        //
        //    Input, double B[N], the right hand side of the linear system.
        //
        //    Output, double R8GE_FS_NEW[N], the solution of the linear system.
        //
    {
        int i;
        int jcol;

        double[] x = new double[n];

        for (i = 0; i < n; i++)
        {
            x[i] = b[i];
        }

        for (jcol = 1; jcol <= n; jcol++)
        {
            //
            //  Find the maximum element in column I.
            //
            double piv = Math.Abs(a[jcol - 1 + (jcol - 1) * n]);
            int ipiv = jcol;
            for (i = jcol + 1; i <= n; i++)
            {
                if (piv < Math.Abs(a[i - 1 + (jcol - 1) * n]))
                {
                    piv = Math.Abs(a[i - 1 + (jcol - 1) * n]);
                    ipiv = i;
                }
            }

            switch (piv)
            {
                case 0.0:
                    Console.WriteLine("");
                    Console.WriteLine("R8GE_FS_NEW - Fatal error!");
                    Console.WriteLine("  Zero pivot on step " + jcol + "");
                    return null;
            }

            //
            //  Switch rows JCOL and IPIV, and X.
            //
            int j;
            double t;
            if (jcol != ipiv)
            {
                for (j = 1; j <= n; j++)
                {
                    t = a[jcol - 1 + (j - 1) * n];
                    a[jcol - 1 + (j - 1) * n] = a[ipiv - 1 + (j - 1) * n];
                    a[ipiv - 1 + (j - 1) * n] = t;
                }

                t = x[jcol - 1];
                x[jcol - 1] = x[ipiv - 1];
                x[ipiv - 1] = t;
            }

            //
            //  Scale the pivot row.
            //
            t = a[jcol - 1 + (jcol - 1) * n];
            a[jcol - 1 + (jcol - 1) * n] = 1.0;
            for (j = jcol + 1; j <= n; j++)
            {
                a[jcol - 1 + (j - 1) * n] /= t;
            }

            x[jcol - 1] /= t;
            //
            //  Use the pivot row to eliminate lower entries in that column.
            //
            for (i = jcol + 1; i <= n; i++)
            {
                if (a[i - 1 + (jcol - 1) * n] != 0.0)
                {
                    t = -a[i - 1 + (jcol - 1) * n];
                    a[i - 1 + (jcol - 1) * n] = 0.0;
                    for (j = jcol + 1; j <= n; j++)
                    {
                        a[i - 1 + (j - 1) * n] += t * a[jcol - 1 + (j - 1) * n];
                    }

                    x[i - 1] += t * x[jcol - 1];
                }
            }
        }

        //
        //  Back solve.
        //
        for (jcol = n; 2 <= jcol; jcol--)
        {
            for (i = 1; i < jcol; i++)
            {
                x[i - 1] -= a[i - 1 + (jcol - 1) * n] * x[jcol - 1];
            }
        }

        return x;
    }

    public static void r8ge_fss(int n, ref double[] a, int nb, ref double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_FSS factors and solves multiple R8GE systems.
        //
        //  Discussion:
        //
        //    The R8GE storage format is used for a "general" M by N matrix.  
        //    A physical storage space is made for each logical entry.  The two 
        //    dimensional logical array is mapped to a vector, in which storage is 
        //    by columns.
        //
        //    This routine does not save the LU factors of the matrix, and hence cannot
        //    be used to efficiently solve multiple linear systems, or even to
        //    factor A at one time, and solve a single linear system at a later time.
        //
        //    This routine uses partial pivoting, but no pivot vector is required.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 June 2009
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
        //    Input/output, double A[N*N].
        //    On input, A is the coefficient matrix of the linear system.
        //    On output, A is in unit upper triangular form, and
        //    represents the U factor of an LU factorization of the
        //    original coefficient matrix.
        //
        //    Input, int NB, the number of right hand sides.
        //
        //    Input/output, double B[N*NB], on input, the right hand sides.
        //    on output, the solutions of the linear systems.
        //
    {
        int i;
        int j;
        int jcol;

        for (jcol = 1; jcol <= n; jcol++)
        {
            //
            //  Find the maximum element in column I.
            //
            double piv = Math.Abs(a[jcol - 1 + (jcol - 1) * n]);
            int ipiv = jcol;
            for (i = jcol + 1; i <= n; i++)
            {
                if (piv < Math.Abs(a[i - 1 + (jcol - 1) * n]))
                {
                    piv = Math.Abs(a[i - 1 + (jcol - 1) * n]);
                    ipiv = i;
                }
            }

            switch (piv)
            {
                case 0.0:
                    Console.WriteLine("");
                    Console.WriteLine("R8GE_FSS - Fatal error!");
                    Console.WriteLine("  Zero pivot on step " + jcol + "");
                    return;
            }

            //
            //  Switch rows JCOL and IPIV, and X.
            //
            double t;
            if (jcol != ipiv)
            {
                for (j = 1; j <= n; j++)
                {
                    t = a[jcol - 1 + (j - 1) * n];
                    a[jcol - 1 + (j - 1) * n] = a[ipiv - 1 + (j - 1) * n];
                    a[ipiv - 1 + (j - 1) * n] = t;
                }

                for (j = 0; j < nb; j++)
                {
                    t = b[jcol - 1 + j * n];
                    b[jcol - 1 + j * n] = b[ipiv - 1 + j * n];
                    b[ipiv - 1 + j * n] = t;
                }
            }

            //
            //  Scale the pivot row.
            //
            t = a[jcol - 1 + (jcol - 1) * n];
            a[jcol - 1 + (jcol - 1) * n] = 1.0;
            for (j = jcol + 1; j <= n; j++)
            {
                a[jcol - 1 + (j - 1) * n] /= t;
            }

            for (j = 0; j < nb; j++)
            {
                b[jcol - 1 + j * n] /= t;
            }

            //
            //  Use the pivot row to eliminate lower entries in that column.
            //
            for (i = jcol + 1; i <= n; i++)
            {
                if (a[i - 1 + (jcol - 1) * n] != 0.0)
                {
                    t = -a[i - 1 + (jcol - 1) * n];
                    a[i - 1 + (jcol - 1) * n] = 0.0;
                    for (j = jcol + 1; j <= n; j++)
                    {
                        a[i - 1 + (j - 1) * n] += t * a[jcol - 1 + (j - 1) * n];
                    }

                    for (j = 0; j < nb; j++)
                    {
                        b[i - 1 + j * n] += t * b[jcol - 1 + j * n];
                    }
                }
            }
        }

        //
        //  Back solve.
        //
        for (jcol = n; 2 <= jcol; jcol--)
        {
            for (i = 1; i < jcol; i++)
            {
                for (j = 0; j < nb; j++)
                {
                    b[i - 1 + j * n] -= a[i - 1 + (jcol - 1) * n] * b[jcol - 1 + j * n];
                }
            }
        }
    }

    public static double[] r8ge_sl_it(int n, double[] a, double[] a_lu, int[] pivot, double[] b,
            int job, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_SL_IT applies one step of iterative refinement following R8GE_SL.
        //
        //  Discussion:
        //
        //    The R8GE storage format is used for a "general" M by N matrix.  
        //    A physical storage space is made for each logical entry.  The two 
        //    dimensional logical array is mapped to a vector, in which storage is 
        //    by columns.
        //
        //    It is assumed that:
        //
        //    * the original matrix A has been factored by R8GE_FA;
        //    * the linear system A * x = b has been solved once by R8GE_SL.
        //
        //    (Actually, it is not necessary to solve the system once using R8GE_SL.
        //    You may simply supply the initial estimated solution X = 0.)
        //
        //    Each time this routine is called, it will compute the residual in
        //    the linear system, apply one step of iterative refinement, and
        //    add the computed correction to the current solution.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 October 2003
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
        //    Input, double A[N*N], the original, UNFACTORED R8GE matrix.
        //
        //    Input, double A_LU[N*N], the LU factors from R8GE_FA.
        //
        //    Input, int PIVOT[N], the pivot vector from R8GE_FA.
        //
        //    Input, double B[N], the right hand side vector.
        //
        //    Input, int JOB, specifies the operation.
        //    0, solve A*X=B.
        //    nonzero, solve A'*X=B.
        //
        //    Input, double X[N], an estimate of the solution of A * x = b.
        //
        //    Output, double R8GE_SL_IT[N], the solution after one step of 
        //    iterative refinement.
        //
    {
        int i;
        //
        //  Compute the residual vector.
        //
        double[] r = r8ge_res(n, n, a, x, b);
        //
        //  Solve A * dx = r
        //
        double[] dx = r8ge_sl_new(n, a_lu, pivot, r, job);
        //
        //  Add dx to x.
        //
        double[] x_new = new double[n];

        for (i = 0; i < n; i++)
        {
            x_new[i] = x[i] + dx[i];
        }

        return x_new;
    }

    public static double[] r8ge_sl_new(int n, double[] a_lu, int[] pivot, double[] b, int job, int aluIndex = 0, int pivotIndex = 0, int bIndex = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_SL_NEW solves an R8GE system factored by R8GE_FA.
        //
        //  Discussion:
        //
        //    The R8GE storage format is used for a "general" M by N matrix.  
        //    A physical storage space is made for each logical entry.  The two 
        //    dimensional logical array is mapped to a vector, in which storage is 
        //    by columns.
        //
        //    R8GE_SL is a simplified version of the LINPACK routine SGESL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 September 2003
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
        //    Input, double A_LU[N*N], the LU factors from R8GE_FA.
        //
        //    Input, int PIVOT[N], the pivot vector from R8GE_FA.
        //
        //    Input, double B[N], the right hand side vector.
        //
        //    Input, int JOB, specifies the operation.
        //    0, solve A * x = b.
        //    nonzero, solve A' * x = b.
        //
        //    Output, double R8GE_SL[N], the solution vector.
        //
    {
        int i;
        int k;
        int l;
        double t;
        //
        double[] x = new double[n];

        for (i = 0; i < n; i++)
        {
            x[i] = b[(bIndex + i) % b.Length];
        }

        switch (job)
        {
            //
            //  Solve A * x = b.
            //
            case 0:
            {
                //
                //  Solve PL * Y = B.
                //
                for (k = 1; k <= n - 1; k++)
                {
                    l = pivot[(k - 1 + pivotIndex) % pivot.Length];

                    if (l != k)
                    {
                        t = x[l - 1];
                        x[l - 1] = x[k - 1];
                        x[k - 1] = t;
                    }

                    for (i = k + 1; i <= n; i++)
                    {
                        x[i - 1] += a_lu[(i - 1 + (k - 1) * n + aluIndex) % a_lu.Length] * x[k - 1];
                    }
                }

                //
                //  Solve U * X = Y.
                //
                for (k = n; 1 <= k; k--)
                {
                    x[k - 1] /= a_lu[(k - 1 + (k - 1) * n + aluIndex) % a_lu.Length];
                    for (i = 1; i <= k - 1; i++)
                    {
                        x[i - 1] -= a_lu[(i - 1 + (k - 1) * n + aluIndex) % a_lu.Length] * x[k - 1];
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
                    t = 0.0;
                    for (i = 1; i <= k - 1; i++)
                    {
                        t += x[i - 1] * a_lu[(i - 1 + (k - 1) * n + aluIndex) % a_lu.Length];
                    }

                    x[k - 1] = (x[k - 1] - t) / a_lu[(k - 1 + (k - 1) * n + aluIndex) % a_lu.Length];
                }

                //
                //  Solve ( PL )' * X = Y.
                //
                for (k = n - 1; 1 <= k; k--)
                {
                    t = 0.0;
                    for (i = k + 1; i <= n; i++)
                    {
                        t += x[i - 1] * a_lu[(i - 1 + (k - 1) * n + aluIndex) % a_lu.Length];
                    }

                    x[k - 1] += t;

                    l = pivot[(k - 1 + pivotIndex) % pivot.Length];

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

        return x;
    }

    public static double[] r8ge_to_r8gb(int m, int n, int ml, int mu, ref double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_TO_R8GB copies an R8GE matrix to an R8GB matrix.
        //
        //  Discussion:
        //
        //    The R8GE storage format is used for a "general" M by N matrix.  
        //    A physical storage space is made for each logical entry.  The two 
        //    dimensional logical array is mapped to a vector, in which storage is 
        //    by columns.
        //
        //    The R8GB storage format is for an M by N banded matrix, with lower
        //    bandwidth ML and upper bandwidth MU.  Storage includes room for ML
        //    extra superdiagonals, which may be required to store nonzero entries
        //    generated during Gaussian elimination.
        //
        //    It usually doesn't make sense to try to store a general matrix
        //    in a band matrix format.  You can always do it, but it will take
        //    more space, unless the general matrix is actually banded.
        //
        //    The purpose of this routine is to allow a user to set up a
        //    banded matrix in the easy-to-use general format, and have this
        //    routine take care of the compression of the data into general
        //    format.  All the user has to do is specify the bandwidths.
        //
        //    Note that this routine "believes" what the user says about the
        //    bandwidth.  It will assume that all entries in the general matrix
        //    outside of the bandwidth are zero.
        //
        //    The original M by N matrix is "collapsed" downward, so that diagonals
        //    become rows of the storage array, while columns are preserved.  The
        //    collapsed array is logically 2*ML+MU+1 by N.
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
        //    16 March 2005
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
        //    Jack Dongarra, Jim Bunch, Cleve Moler, Pete Stewart,
        //    LINPACK User's Guide,
        //    SIAM, 1979,
        //    ISBN13: 978-0-898711-72-1,
        //    LC: QA214.L56.
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
        //    Output, double A[M*N], the R8GE matrix.
        //
        //    Input, double R8GE_TO_R8GB[(2*ML+MU+1)*N], the R8GB matrix.
        //
    {
        int i;
        int k;

        double[] b = new double[(2 * ml + mu + 1) * n];

        for (k = 0; k < (2 * ml + mu + 1) * n; k++)
        {
            b[k] = 0.0;
        }

        for (i = 1; i <= m; i++)
        {
            int jlo = Math.Max(i - ml, 1);
            int jhi = Math.Min(i + mu, n);

            int j;
            for (j = jlo; j <= jhi; j++)
            {
                b[ml + mu + i - j + (j - 1) * (2 * ml + mu + 1)] = a[i - 1 + (j - 1) * m];
            }
        }

        return b;
    }

    public static double[] r8ge_to_r8lt(int m, int n, double[] a_ge)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_TO_R8LT copies an R8GE matrix to an R8LT matrix.
        //
        //  Discussion:
        //
        //    The R8GE storage format is used for a general M by N matrix.  A storage 
        //    space is made for each entry.  The two dimensional logical
        //    array can be thought of as a vector of M*N entries, starting with
        //    the M entries in the column 1, then the M entries in column 2
        //    and so on.  Considered as a vector, the entry A(I,J) is then stored
        //    in vector location I+(J-1)*M.
        //
        //    The R8LT storage format is used for an M by N lower triangular matrix,
        //    and allocates space even for the zero entries.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    21 August 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of the matrix.
        //
        //    Input, double A_GE[M,N], the R8GE matrix.
        //
        //    Output, double R8GE_TO_R8LT[M,N], the R8LT matrix.
        //
    {
        int j;

        double[] a_lt = new double[m * n];

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                if (j <= i)
                {
                    a_lt[i + j * m] = a_ge[i + j * m];
                }
                else
                {
                    a_lt[i + j * m] = 0.0;
                }
            }
        }

        return a_lt;
    }

    public static double[] r8ge_to_r8po(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_TO_R8PO copies an R8GE matrix to an R8PO matrix.
        //
        //  Discussion:
        //
        //    The R8PO format assumes the matrix is square and symmetric; it is also 
        //    typically assumed that the matrix is positive definite.  These are not
        //    required here.  The copied R8PO matrix simply zeros out the lower triangle
        //    of the R8GE matrix.
        //
        //    The R8GE storage format is used for a general M by N matrix.  A storage 
        //    space is made for each entry.  The two dimensional logical
        //    array can be thought of as a vector of M*N entries, starting with
        //    the M entries in the column 1, then the M entries in column 2
        //    and so on.  Considered as a vector, the entry A(I,J) is then stored
        //    in vector location I+(J-1)*M.
        //
        //    The R8PO storage format is used for a symmetric positive definite 
        //    matrix and its inverse.  (The Cholesky factor of an R8PO matrix is an
        //    upper triangular matrix, so it will be in R8GE storage format.)
        //
        //    Only the diagonal and upper triangle of the square array are used.
        //    This same storage scheme is used when the matrix is factored by
        //    R8PO_FA, or inverted by R8PO_INVERSE.  For clarity, the lower triangle
        //    is set to zero.
        //
        //    R8PO storage is used by LINPACK and LAPACK.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 August 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A[N,N], the R8GE matrix.
        //
        //    Output, double R8GE_TO_R8PO[N,N], the R8PO matrix.
        //
    {
        int i;

        double[] b = new double[n * n];

        for (i = 0; i < n; i++)
        {
            int j;
            for (j = 0; j < n; j++)
            {
                if (i <= j)
                {
                    b[i + j * n] = a[i + j * n];
                }
                else
                {
                    b[i + j * n] = 0.0;
                }
            }
        }

        return b;
    }

    public static void r8ge_to_r8ri(int n, double[] a, int nz, ref int[] ija, ref double[] sa)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_TO_R8RI converts an R8GE matrix to R8RI form.
        //
        //  Discussion:
        //
        //    A R8GE matrix is in general storage.
        //
        //    An R8RI matrix is in row indexed sparse storage form.
        //
        //    The size of the arrays IJA and SA can be determined by calling
        //    R8GE_TO_R8RI_SIZE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 January 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
        //    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
        //    Third Edition,
        //    Cambridge University Press, 2007,
        //    ISBN13: 978-0-521-88068-8,
        //    LC: QA297.N866.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A[N*N], the matrix stored in GE 
        //    or "general" format.
        //
        //    Input, int NZ, the size required for the RI
        //    or "row indexed" sparse storage.
        //
        //    Output, int IJA[NZ], the index vector.
        //
        //    Output, double SA[NZ], the value vector.
        //
    {
        int i;
        int j;
        int k;

        for (k = 0; k < n; k++)
        {
            i = k;
            j = k;
            sa[k] = a[i + j * n];
        }

        k = n;
        sa[k] = 0.0;

        for (i = 0; i <= n; i++)
        {
            ija[i] = 0;
        }

        int im = 0;

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                if (i != j)
                {
                    if (a[i + j * n] != 0.0)
                    {
                        k += 1;
                        switch (ija[i])
                        {
                            case 0:
                            {
                                int l;
                                for (l = im; l <= i; l++)
                                {
                                    ija[l] = k;
                                }

                                im = i + 1;
                                break;
                            }
                        }

                        ija[k] = j;
                        sa[k] = a[i + j * n];
                    }
                }
            }
        }

        ija[n] = k + 1;

    }

    public static int r8ge_to_r8ri_size(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_TO_R8RI_SIZE determines the size of an R8RI matrix.
        //
        //  Discussion:
        //
        //    N spaces are always used for the diagonal entries, plus a dummy.
        //    The remaining spaces store off-diagonal nonzeros.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 January 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    William Press, Brian Flannery, Saul Teukolsky, William Vetterling,
        //    Numerical Recipes in FORTRAN: The Art of Scientific Computing,
        //    Third Edition,
        //    Cambridge University Press, 2007,
        //    ISBN13: 978-0-521-88068-8,
        //    LC: QA297.N866.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A[N*N], the matrix stored in GE or "general" format.
        //
        //    Output, int R8GE_TO_R8RI_SIZE, the size required for the RI
        //    or "row indexed" sparse storage.
        //
    {
        int i;

        int nz = n + 1;

        for (i = 0; i < n; i++)
        {
            int j;
            for (j = 0; j < n; j++)
            {
                if (i != j)
                {
                    if (a[i + j * n] != 0.0)
                    {
                        nz += 1;
                    }
                }
            }
        }

        return nz;
    }

    public static double[] r8ge_to_r8ut(int m, int n, ref double[] a_ge)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_TO_R8UT copies an R8GE matrix to an R8UT matrix.
        //
        //  Discussion:
        //
        //    The R8GE storage format is used for a general M by N matrix.  A storage 
        //    space is made for each entry.  The two dimensional logical
        //    array can be thought of as a vector of M*N entries, starting with
        //    the M entries in the column 1, then the M entries in column 2
        //    and so on.  Considered as a vector, the entry A(I,J) is then stored
        //    in vector location I+(J-1)*M.
        //
        //    The R8UT storage format is used for an M by N upper triangular matrix,
        //    and allocates space even for the zero entries.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 August 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of the matrix.
        //
        //    Input, double A_GE[M,N], the R8GE matrix.
        //
        //    Output, double R8GE_TO_R8UT[M,N], the R8UT matrix.
        //
    {
        int j;

        double[] a_ut = new double[m * n];

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                if (i <= j)
                {
                    a_ut[i + j * m] = a_ge[i + j * m];
                }
                else
                {
                    a_ut[i + j * m] = 0.0;
                }
            }
        }

        return a_ut;
    }

    public static double[] r8ge_to_r8vec(int m, int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_TO_R8VEC copies an R8GE matrix to a real vector.
        //
        //  Discussion:
        //
        //    The R8GE storage format is used for a "general" M by N matrix.  
        //    A physical storage space is made for each logical entry.  The two 
        //    dimensional logical array is mapped to a vector, in which storage is 
        //    by columns.
        //
        //    In C++  and FORTRAN, this routine is not really needed.  In MATLAB,
        //    a data item carries its dimensionality implicitly, and so cannot be
        //    regarded sometimes as a vector and sometimes as an array.
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
        //    Input, double R8VEC_TO_R8GE[M*N], the array to be copied.
        //
        //    Output, double X[M*N], the vector.
        //
    {
        int j;

        double[] x = new double[m * n];

        int k = 0;
        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                x[k] = a[i + j * m];
                k += 1;
            }
        }

        return x;
    }

    public static double[] r8ge_to_r8vm(int m, int n, double[] a_ge)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_TO_R8VM converts an R8GE matrix to an R8VM matrix.
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
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 August 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns of the matrix.
        //
        //    Input, double A_GE[M*N], the R8GE matrix.
        //
        //    Output, double R8GE_TO_R8VM[N], the R8VM matrix.
        //
    {
        int j;

        double[] a_vm = new double[n];

        int i = 1;
        for (j = 0; j < n; j++)
        {
            a_vm[j] = a_ge[i + j * m];
        }

        return a_vm;
    }

    public static double[] r8ge_transpose_new(int m, int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_TRANSPOSE_NEW returns the transpose of an R8GE.
        //
        //  Discussion:
        //
        //    An R8GE is a doubly dimensioned array of R8 values, stored as a vector
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
        //    Input, int M, N, the number of rows and columns of the matrix A.
        //
        //    Input, double A[M*N], the matrix whose transpose is desired.
        //
        //    Output, double R8GE_TRANSPOSE_NEW[N*M], the transposed matrix.
        //
    {
        int j;

        double[] b = new double[n * m];

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                b[j + i * n] = a[i + j * m];
            }
        }

        return b;
    }

    public static void r8ge_transpose_print(int m, int n, double[] a, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_TRANSPOSE_PRINT prints an R8GE, transposed.
        //
        //  Discussion:
        //
        //    An R8GE is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 September 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double A[M*N], an M by N matrix to be printed.
        //
        //    Input, string TITLE, a title.
        //
    {
        r8ge_transpose_print_some(m, n, a, 1, 1, m, n, title);
    }

    public static void r8ge_transpose_print_some(int m, int n, double[] a, int ilo, int jlo,
            int ihi, int jhi, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_TRANSPOSE_PRINT_SOME prints some of an R8GE, transposed.
        //
        //  Discussion:
        //
        //    An R8GE is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 April 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Input, double A[M*N], an M by N matrix to be printed.
        //
        //    Input, int ILO, JLO, the first row and column to print.
        //
        //    Input, int IHI, JHI, the last row and column to print.
        //
        //    Input, string TITLE, a title.
        //
    {
        const int INCX = 5;

        int i2lo;

        Console.WriteLine("");
        Console.WriteLine(title + "");

        if (m <= 0 || n <= 0)
        {
            Console.WriteLine("");
            Console.WriteLine("  (None)");
            return;
        }

        int i2lo_lo = ilo switch
        {
            < 1 => 1,
            _ => ilo
        };

        int i2lo_hi = ihi < m ? m : ihi;

        for (i2lo = i2lo_lo; i2lo <= i2lo_hi; i2lo += INCX)
        {
            int i2hi = i2lo + INCX - 1;

            if (m < i2hi)
            {
                i2hi = m;
            }

            if (ihi < i2hi)
            {
                i2hi = ihi;
            }

            int inc = i2hi + 1 - i2lo;

            Console.WriteLine("");
            string cout = "  Row: ";
            int i;
            for (i = i2lo; i <= i2hi; i++)
            {
                cout += (i - 1).ToString().PadLeft(7) + "       ";
            }

            Console.WriteLine(cout);
            Console.WriteLine("  Col");
            Console.WriteLine("");

            int j2lo = jlo switch
            {
                < 1 => 1,
                _ => jlo
            };

            int j2hi = n < jhi ? n : jhi;

            int j;
            for (j = j2lo; j <= j2hi; j++)
            {
                cout = (j - 1).ToString().PadLeft(5) + ":";
                int i2;
                for (i2 = 1; i2 <= inc; i2++)
                {
                    i = i2lo - 1 + i2;
                    cout += a[i - 1 + (j - 1) * m].ToString().PadLeft(14);
                }

                Console.WriteLine(cout);
            }
        }
    }

    public static int r8ge_trf(int m, int n, ref double[] a, ref int[] pivot)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_TRF performs a LAPACK-style PLU factorization of an R8GE matrix.
        //
        //  Discussion:
        //
        //    The R8GE storage format is used for a "general" M by N matrix.  
        //    A physical storage space is made for each logical entry.  The two 
        //    dimensional logical array is mapped to a vector, in which storage is 
        //    by columns.
        //
        //    R8GE_TRF is a standalone version of the LAPACK routine SGETRF.
        //
        //    The factorization uses partial pivoting with row interchanges,
        //    and has the form
        //      A = P * L * U
        //    where P is a permutation matrix, L is lower triangular with unit
        //    diagonal elements (lower trapezoidal if N < M), and U is upper
        //    triangular (upper trapezoidal if M < N).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 November 2003
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Anderson, Bai, Bischof, Blackford,
        //    Demmel, Dongarra, DuCroz, Greenbaum, Hammarling, McKenney, Sorensen.
        //    C++ version by John Burkardt
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
        //    Input/output, double A[M*N].
        //    On entry, the M by N matrix to be factored.
        //    On exit, the factors L and U from the factorization
        //    A = P*L*U; the unit diagonal elements of L are not stored.
        //
        //    Output, int PIVOT[min(M,N)], the pivot indices.
        //
        //    Output, int R8GE_TRF.
        //    = 0: successful exit
        //    = -K, the K-th argument had an illegal value
        //    = K: U(K,K) is exactly zero. The factorization
        //         has been completed, but the factor U is exactly
        //         singular, and division by zero will occur if it is used
        //         to solve a system of equations.
        //
    {
        int j;
        //
        //  Test the input parameters.
        //
        int info = 0;

        switch (m)
        {
            case < 0:
                return -1;
        }

        switch (n)
        {
            case < 0:
                return -2;
        }

        if (m == 0 || n == 0)
        {
            return 0;
        }

        for (j = 1; j <= Math.Min(m, n); j++)
        {
            //
            //  Find the pivot.
            //
            double temp = Math.Abs(a[j - 1 + (j - 1) * m]);
            int jp = j;
            int i;
            for (i = j + 1; i <= m; i++)
            {
                if (temp < Math.Abs(a[i - 1 + (j - 1) * m]))
                {
                    temp = Math.Abs(a[i - 1 + (j - 1) * m]);
                    jp = i;
                }
            }

            pivot[j - 1] = jp;
            //
            //  Apply the interchange to columns 1:N.
            //  Compute elements J+1:M of the J-th column.
            //
            if (a[jp - 1 + (j - 1) * m] != 0.0)
            {
                if (jp != j)
                {
                    int jj;
                    for (jj = 1; jj <= n; jj++)
                    {
                        temp = a[j - 1 + (jj - 1) * m];
                        a[j - 1 + (jj - 1) * m] = a[jp - 1 + (jj - 1) * m];
                        a[jp - 1 + (jj - 1) * m] = temp;
                    }
                }

                if (j < m)
                {
                    for (i = j + 1; i <= m; i++)
                    {
                        a[i - 1 + (j - 1) * m] /= a[j - 1 + (j - 1) * m];
                    }
                }
            }
            else
            {
                info = info switch
                {
                    0 => j,
                    _ => info
                };
            }

            //
            //  Update the trailing submatrix.
            //
            if (j >= Math.Min(m, n))
            {
                continue;
            }

            int ii;
            for (ii = j + 1; ii <= m; ii++)
            {
                for (i = j + 1; i <= n; i++)
                {
                    a[ii - 1 + (i - 1) * m] -= a[ii - 1 + (j - 1) * m] * a[j - 1 + (i - 1) * m];
                }
            }
        }

        return info;
    }

    public static double[] r8ge_trs(int n, int nrhs, char trans, double[] a, int[] pivot,
            double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_TRS solves a system of linear equations factored by R8GE_TRF.
        //
        //  Discussion:
        //
        //    The R8GE storage format is used for a "general" M by N matrix.  
        //    A physical storage space is made for each logical entry.  The two 
        //    dimensional logical array is mapped to a vector, in which storage is 
        //    by columns.
        //
        //    R8GE_TRS is a standalone version of the LAPACK routine SGETRS.
        //
        //    R8GE_TRS solves a system of linear equations
        //      A * x = b  or  A' * X = B
        //    with a general N by N matrix A using the PLU factorization computed
        //    by R8GE_TRF.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 November 2003
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Anderson, Bai, Bischof, Blackford,
        //    Demmel, Dongarra, DuCroz, Greenbaum, Hammarling, McKenney, Sorensen.
        //    C++ version by John Burkardt
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
        //    Input, int N, the order of the matrix A.  0 <= N.
        //
        //    Input, int NRHS, the number of right hand sides.  0 <= NRHS.
        //
        //    Input, char TRANS, specifies the form of the system of equations:
        //    'N':  A * x = b  (No transpose)
        //    'T':  A'* X = B  (Transpose)
        //    'C':  A'* X = B  (Conjugate transpose = Transpose)
        //
        //    Input, double A[N*N], the factors L and U from the factorization
        //    A = P*L*U as computed by R8GE_TRF.
        //
        //    Input, int PIVOT[N], the pivot indices from R8GE_TRF.
        //
        //    Input, double B[N*NRHS], the right hand side matrix.
        //
        //    Output, double R8GE_TRS[N*NRHS], the solution matrix X.
        //
    {
        int i;
        int j;
        int k;
        double temp;

        if (trans != 'n' && trans != 'N' &&
            trans != 't' && trans != 'T' &&
            trans != 'c' && trans != 'C')
        {
            return null;
        }

        switch (n)
        {
            case < 0:
                return null;
        }

        switch (nrhs)
        {
            case < 0:
                return null;
        }

        switch (n)
        {
            case 0:
                return null;
        }

        switch (nrhs)
        {
            case 0:
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

        switch (trans)
        {
            case 'n':
            case 'N':
            {
                //
                //  Apply row interchanges to the right hand sides.
                //
                for (i = 1; i <= n; i++)
                {
                    if (pivot[i - 1] == i)
                    {
                        continue;
                    }

                    for (k = 0; k < nrhs; k++)
                    {
                        temp = x[i - 1 + k * n];
                        x[i - 1 + k * n] = x[pivot[i - 1] - 1 + k * n];
                        x[pivot[i - 1] - 1 + k * n] = temp;
                    }
                }

                //
                //  Solve L * x = b, overwriting b with x.
                //
                for (k = 0; k < nrhs; k++)
                {
                    for (j = 1; j <= n - 1; j++)
                    {
                        for (i = j + 1; i <= n; i++)
                        {
                            x[i - 1 + k * n] -= a[i - 1 + (j - 1) * n] * x[j - 1 + k * n];
                        }
                    }
                }

                //
                //  Solve U * x = b, overwriting b with x.
                //
                for (k = 0; k < nrhs; k++)
                {
                    for (j = n; 1 <= j; j--)
                    {
                        x[j - 1 + k * n] /= a[j - 1 + (j - 1) * n];
                        for (i = 1; i < j; i++)
                        {
                            x[i - 1 + k * n] -= a[i - 1 + (j - 1) * n] * x[j - 1 + k * n];
                        }
                    }
                }

                break;
            }
            default:
            {
                //
                //  Solve U' * x = b, overwriting b with x.
                //
                for (k = 0; k < nrhs; k++)
                {
                    for (j = 1; j <= n; j++)
                    {
                        x[j - 1 + k * n] /= a[j - 1 + (j - 1) * n];
                        for (i = j + 1; i <= n; i++)
                        {
                            x[i - 1 + k * n] -= a[j - 1 + (i - 1) * n] * x[j - 1 + k * n];
                        }
                    }
                }

                //
                //  Solve L' * x = b, overwriting b with x.
                //
                for (k = 0; k < nrhs; k++)
                {
                    for (j = n; 2 <= j; j--)
                    {
                        for (i = 1; i < j; i++)
                        {
                            x[i - 1 + k * n] -= a[j - 1 + (i - 1) * n] * x[j - 1 + k * n];
                        }
                    }
                }

                //
                //  Apply row interchanges to the solution vectors.
                //
                for (i = n; 1 <= i; i--)
                {
                    if (pivot[i - 1] == i)
                    {
                        continue;
                    }

                    for (k = 0; k < nrhs; k++)
                    {
                        temp = x[i - 1 + k * n];
                        x[i - 1 + k * n] = x[pivot[i - 1] - 1 + k * n];
                        x[pivot[i - 1] - 1 + k * n] = temp;
                    }
                }

                break;
            }
        }

        return x;
    }

    public static double[] r8ge_zeros(int m, int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_ZEROS zeros an R8GE matrix.
        //
        //  Discussion:
        //
        //    The R8GE storage format is used for a "general" M by N matrix.  
        //    A physical storage space is made for each logical entry.  The two 
        //    dimensional logical array is mapped to a vector, in which storage is 
        //    by columns.
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
        //    M must be positive.
        //
        //    Input, int N, the number of columns of the matrix.
        //    N must be positive.
        //
        //    Output, double R8GE_ZERO[M*N], the M by N matrix.
        //
    {
        int j;

        double[] a = new double[m * n];

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                a[i + j * m] = 0.0;
            }
        }

        return a;
    }
        
    public static double[] r8ge_zeros_new ( int m, int n )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_ZEROS_NEW returns a new zeroed R8GE.
        //
        //  Discussion:
        //
        //    An R8GE is a doubly dimensioned array of R8 values, stored as a vector
        //    in column-major order.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns.
        //
        //    Output, double R8GE_ZEROS_NEW[M*N], the new zeroed matrix.
        //
    {
        int j;

        double[] a = new double[m*n];

        for ( j = 0; j < n; j++ )
        {
            int i;
            for ( i = 0; i < m; i++ )
            {
                a[i+j*m] = 0.0;
            }
        }
        return a;
    }

}