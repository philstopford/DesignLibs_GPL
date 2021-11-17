using System;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static double[] r83_cr_fa(int n, double[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83_CR_FA decomposes a real tridiagonal matrix using cyclic reduction.
        //
        //  Discussion:
        //
        //    The R83 storage format is used for a tridiagonal matrix.
        //    The superdiagonal is stored in entries (1,2:N), the diagonal in
        //    entries (2,1:N), and the subdiagonal in (3,1:N-1).  Thus, the
        //    original matrix is "collapsed" vertically into the array.
        //
        //    Once R83_CR_FA has decomposed a matrix A, then R83_CR_SL may be used to solve
        //    linear systems A * x = b.
        //
        //    R83_CR_FA does not employ pivoting.  Hence, the results can be more
        //    sensitive to ill-conditioning than standard Gauss elimination.  In
        //    particular, R83_CR_FA will fail if any diagonal element of the matrix
        //    is zero.  Other matrices may also cause R83_CR_FA to fail.
        //
        //    R83_CR_FA can be guaranteed to work properly if the matrix is strictly
        //    diagonally dominant, that is, if the absolute value of the diagonal
        //    element is strictly greater than the sum of the absolute values of
        //    the offdiagonal elements, for each equation.
        //
        //    The algorithm may be illustrated by the following figures:
        //
        //    The initial matrix is given by:
        //
        //          D1 U1
        //          L1 D2 U2
        //             L2 D3 U3
        //                L3 D4 U4
        //                   L4 D U5
        //                      L5 D6
        //
        //    Rows and columns are permuted in an odd/even way to yield:
        //
        //          D1       U1
        //             D3    L2 U3
        //                D5    L4 U5
        //          L1 U2    D2
        //             L3 U4    D4
        //                L5       D6
        //
        //    A block LU decomposition is performed to yield:
        //
        //          D1      |U1
        //             D3   |L2 U3
        //                D5|   L4 U5
        //          --------+--------
        //                  |D2'F3
        //                  |F1 D4'F4
        //                  |   F2 D6'
        //
        //    For large systems, this reduction is repeated on the lower right hand
        //    tridiagonal subsystem until a completely upper triangular system
        //    is obtained.  The system has now been factored into the product of a
        //    lower triangular system and an upper triangular one, and the information
        //    defining this factorization may be used by R83_CR_SL to solve linear
        //    systems.
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
        //    23 March 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Roger Hockney,
        //    A fast direct solution of Poisson's equation using Fourier Analysis,
        //    Journal of the ACM,
        //    Volume 12, Number 1, pages 95-113, January 1965.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be positive.
        //
        //    Input, double A[3*N], the R83 matrix.
        //
        //    Output, double R83_CR_FA[3*(2*N+1)], factorization information 
        //    needed by R83_CR_SL.
        //
    {
        double[] a_cr;
        int iful;
        int ifulp;
        int ihaf;
        int il;
        int ilp;
        int inc;
        int incr;
        int ipnt;
        int ipntp;
        int j;

        switch (n)
        {
            case <= 0:
                Console.WriteLine("");
                Console.WriteLine("R83_CR_FA - Fatal error!");
                Console.WriteLine("  Nonpositive N = " + n + "");
                return null;
        }

        a_cr = new double[3 * (2 * n + 1)];

        switch (n)
        {
            case 1:
                a_cr[0 + 0 * 3] = 0.0;
                a_cr[0 + 1 * 3] = 0.0;
                a_cr[0 + 2 * 3] = 0.0;
                a_cr[1 + 0 * 3] = 0.0;
                a_cr[1 + 1 * 3] = 1.0 / a[1 + 0 * 3];
                a_cr[1 + 2 * 3] = 0.0;
                a_cr[2 + 0 * 3] = 0.0;
                a_cr[2 + 1 * 3] = 0.0;
                a_cr[2 + 2 * 3] = 0.0;

                return a_cr;
        }

        //
        //  Zero out the workspace entries.
        //
        a_cr[0 + 0 * 3] = 0.0;
        for (j = 1; j <= n - 1; j++)
        {
            a_cr[0 + j * 3] = a[0 + j * 3];
        }

        for (j = n; j <= 2 * n; j++)
        {
            a_cr[0 + j * 3] = 0.0;
        }

        a_cr[1 + 0 * 3] = 0.0;
        for (j = 1; j <= n; j++)
        {
            a_cr[1 + j * 3] = a[1 + (j - 1) * 3];
        }

        for (j = n + 1; j <= 2 * n; j++)
        {
            a_cr[1 + j * 3] = 0.0;
        }

        a_cr[2 + 0 * 3] = 0.0;
        for (j = 1; j <= n - 1; j++)
        {
            a_cr[2 + j * 3] = a[2 + (j - 1) * 3];
        }

        for (j = n; j <= 2 * n; j++)
        {
            a_cr[2 + j * 3] = 0.0;
        }

        il = n;
        ipntp = 0;

        while (1 < il)
        {
            ipnt = ipntp;
            ipntp += il;
            inc = (il % 2) switch
            {
                1 => il + 1,
                _ => il
            };

            incr = inc / 2;
            il /= 2;
            ihaf = ipntp + incr + 1;
            ifulp = ipnt + inc + 2;

            for (ilp = incr; 1 <= ilp; ilp--)
            {
                ifulp -= 2;
                iful = ifulp - 1;
                ihaf -= 1;

                a_cr[1 + iful * 3] = 1.0 / a_cr[1 + iful * 3];
                a_cr[2 + iful * 3] *= a_cr[1 + iful * 3];
                a_cr[0 + ifulp * 3] *= a_cr[1 + (ifulp + 1) * 3];
                a_cr[1 + ihaf * 3] = a_cr[1 + ifulp * 3]
                                     - a_cr[0 + iful * 3] * a_cr[2 + iful * 3]
                                     - a_cr[0 + ifulp * 3] * a_cr[2 + ifulp * 3];
                a_cr[2 + ihaf * 3] = -a_cr[2 + ifulp * 3] * a_cr[2 + (ifulp + 1) * 3];
                a_cr[0 + ihaf * 3] = -a_cr[0 + ifulp * 3] * a_cr[0 + (ifulp + 1) * 3];
            }
        }

        a_cr[1 + (ipntp + 1) * 3] = 1.0 / a_cr[1 + (ipntp + 1) * 3];

        return a_cr;
    }

    public static double[] r83_cr_sl(int n, double[] a_cr, double[] b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83_CR_SL solves a real linear system factored by R83_CR_FA.
        //
        //  Discussion:
        //
        //    The matrix A must be tridiagonal.  R83_CR_FA is called to compute the
        //    LU factors of A.  It does so using a form of cyclic reduction.  If
        //    the factors computed by R83_CR_FA are passed to R83_CR_SL, then one or 
        //    many linear systems involving the matrix A may be solved.
        //
        //    Note that R83_CR_FA does not perform pivoting, and so the solution 
        //    produced by R83_CR_SL may be less accurate than a solution produced 
        //    by a standard Gauss algorithm.  However, such problems can be 
        //    guaranteed not to occur if the matrix A is strictly diagonally 
        //    dominant, that is, if the absolute value of the diagonal coefficient 
        //    is greater than the sum of the absolute values of the two off diagonal 
        //    coefficients, for each row of the matrix.
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
        //    01 January 2004
        //
        //  Author:
        //
        //    C++ version by John Burkardt
        //
        //  Reference:
        //
        //    Roger Hockney,
        //    A fast direct solution of Poisson's equation using Fourier Analysis,
        //    Journal of the ACM,
        //    Volume 12, Number 1, pages 95-113, January 1965.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be positive.
        //
        //    Input, double A_CR[3*(2*N+1)], factorization information computed by R83_CR_FA.
        //
        //    Input, double B[N], the right hand side.
        //
        //    Output, double R83_CR_SL[N], the solution.
        //
    {
        int i;
        int iful;
        int ifulm;
        int ihaf = 0;
        int il;
        int ipnt;
        int ipntp;
        int ndiv;
        double[] rhs;
        double[] x;

        switch (n)
        {
            case <= 0:
                Console.WriteLine("");
                Console.WriteLine("R83_CR_SL - Fatal error!");
                Console.WriteLine("  Nonpositive N = " + n + "");
                return null;
            case 1:
                x = new double[1];
                x[0] = a_cr[1 + 1 * 3] * b[0];
                return x;
        }

        //
        //  Set up RHS.
        //
        rhs = new double[2 * n + 1];

        rhs[0] = 0.0;
        for (i = 1; i <= n; i++)
        {
            rhs[i] = b[i - 1];
        }

        for (i = n + 1; i <= 2 * n; i++)
        {
            rhs[i] = 0.0;
        }

        il = n;
        ndiv = 1;
        ipntp = 0;

        while (1 < il)
        {
            ipnt = ipntp;
            ipntp += il;
            il /= 2;
            ndiv *= 2;
            ihaf = ipntp;

            for (iful = ipnt + 2; iful <= ipntp; iful += 2)
            {
                ihaf += 1;
                rhs[ihaf] = rhs[iful]
                            - a_cr[2 + (iful - 1) * 3] * rhs[iful - 1]
                            - a_cr[0 + iful * 3] * rhs[iful + 1];
            }
        }

        rhs[ihaf] *= a_cr[1 + ihaf * 3];
        ipnt = ipntp;

        while (0 < ipnt)
        {
            ipntp = ipnt;
            ndiv /= 2;
            il = n / ndiv;
            ipnt -= il;
            ihaf = ipntp;

            for (ifulm = ipnt + 1; ifulm <= ipntp; ifulm += 2)
            {
                iful = ifulm + 1;
                ihaf += 1;
                rhs[iful] = rhs[ihaf];
                rhs[ifulm] = a_cr[1 + ifulm * 3] * (
                    rhs[ifulm]
                    - a_cr[2 + (ifulm - 1) * 3] * rhs[ifulm - 1]
                    - a_cr[0 + ifulm * 3] * rhs[iful]);
            }
        }

        x = new double[n];

        for (i = 0; i < n; i++)
        {
            x[i] = rhs[i + 1];
        }

        return x;
    }

    public static double[] r83_cr_sls(int n, double[] a_cr, int nb, double[] b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83_CR_SLS solves several real linear systems factored by R83_CR_FA.
        //
        //  Discussion:
        //
        //    The matrix A must be tridiagonal.  R83_CR_FA is called to compute the
        //    LU factors of A.  It does so using a form of cyclic reduction.  If
        //    the factors computed by R83_CR_FA are passed to R83_CR_SLS, then one or 
        //    many linear systems involving the matrix A may be solved.
        //
        //    Note that R83_CR_FA does not perform pivoting, and so the solutions
        //    produced by R83_CR_SLS may be less accurate than a solution produced 
        //    by a standard Gauss algorithm.  However, such problems can be 
        //    guaranteed not to occur if the matrix A is strictly diagonally 
        //    dominant, that is, if the absolute value of the diagonal coefficient 
        //    is greater than the sum of the absolute values of the two off diagonal 
        //    coefficients, for each row of the matrix.
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
        //    09 May 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Roger Hockney,
        //    A fast direct solution of Poisson's equation using Fourier Analysis,
        //    Journal of the ACM,
        //    Volume 12, Number 1, pages 95-113, January 1965.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be positive.
        //
        //    Input, double A_CR[3*(2*N+1)], factorization information computed by R83_CR_FA.
        //
        //    Input, int NB, the number of systems.
        //
        //    Input, double B[N*NB], the right hand sides.
        //
        //    Output, double R83_CR_SL[N*NB], the solutions.
        //
    {
        int i;
        int iful;
        int ifulm;
        int ihaf = 0;
        int il;
        int ipnt;
        int ipntp;
        int j;
        int ndiv;
        double[] rhs;
        double[] x;

        switch (n)
        {
            case <= 0:
                Console.WriteLine("");
                Console.WriteLine("R83_CR_SLS - Fatal error!");
                Console.WriteLine("  Nonpositive N = " + n + "");
                return null;
            case 1:
            {
                x = new double[1 * nb];
                for (j = 0; j < nb; j++)
                {
                    x[0 + j * n] = a_cr[1 + 1 * 3] * b[0 + j * n];
                }

                return x;
            }
        }

        //
        //  Set up RHS.
        //
        rhs = new double[(2 * n + 1) * nb];

        for (j = 0; j < nb; j++)
        {
            rhs[0 + j * (2 * n + 1)] = 0.0;
            for (i = 1; i <= n; i++)
            {
                rhs[i + j * (2 * n + 1)] = b[i - 1 + j * n];
            }

            for (i = n + 1; i <= 2 * n; i++)
            {
                rhs[i + j * (2 * n + 1)] = 0.0;
            }
        }

        il = n;
        ndiv = 1;
        ipntp = 0;

        while (1 < il)
        {
            ipnt = ipntp;
            ipntp += il;
            il /= 2;
            ndiv *= 2;

            for (j = 0; j < nb; j++)
            {
                ihaf = ipntp;
                for (iful = ipnt + 2; iful <= ipntp; iful += 2)
                {
                    ihaf += 1;
                    rhs[ihaf + j * (2 * n + 1)] = rhs[iful + j * (2 * n + 1)]
                                                  - a_cr[2 + (iful - 1) * 3] * rhs[iful - 1 + j * (2 * n + 1)]
                                                  - a_cr[0 + iful * 3] * rhs[iful + 1 + j * (2 * n + 1)];
                }
            }
        }

        for (j = 0; j < nb; j++)
        {
            rhs[ihaf + j * (2 * n + 1)] *= a_cr[1 + ihaf * 3];
        }

        ipnt = ipntp;

        while (0 < ipnt)
        {
            ipntp = ipnt;
            ndiv /= 2;
            il = n / ndiv;
            ipnt -= il;

            for (j = 0; j < nb; j++)
            {
                ihaf = ipntp;
                for (ifulm = ipnt + 1; ifulm <= ipntp; ifulm += 2)
                {
                    iful = ifulm + 1;
                    ihaf += 1;
                    rhs[iful + j * (2 * n + 1)] = rhs[ihaf + j * (2 * n + 1)];
                    rhs[ifulm + j * (2 * n + 1)] = a_cr[1 + ifulm * 3] * (
                        rhs[ifulm + j * (2 * n + 1)]
                        - a_cr[2 + (ifulm - 1) * 3] * rhs[ifulm - 1 + j * (2 * n + 1)]
                        - a_cr[0 + ifulm * 3] * rhs[iful + j * (2 * n + 1)]);
                }
            }
        }

        x = new double[n * nb];

        for (j = 0; j < nb; j++)
        {
            for (i = 0; i < n; i++)
            {
                x[i + j * n] = rhs[i + 1 + j * (2 * n + 1)];
            }
        }

        return x;
    }
}