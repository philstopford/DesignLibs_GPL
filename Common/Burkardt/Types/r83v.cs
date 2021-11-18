using System;
using Burkardt.Uniform;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static void r83v_cg(int n, double[] a1, double[] a2, double[] a3, double[] b,
            ref double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83V_CG uses the conjugate gradient method on an R83V system.
        //
        //  Discussion:
        //
        //    The R83V storage format is used for a tridiagonal matrix.
        //    The subdiagonal is in A(min(M-1,N)).
        //    The diagonal is in B(min(M,N)).
        //    The superdiagonal is in C(min(M,N-1)).
        //
        //  Example:
        //
        //    An R83V matrix of order 3x5 would be stored:
        //
        //      B1  C1  **  **  **
        //      A1  B2  C2  **  **
        //      **  A2  B3  C3  **
        //
        //    An R83 matrix of order 5x5 would be stored:
        //
        //      B1  C1  **  **  **
        //      A1  B2  C2  **  **
        //      **  A2  B3  C3  **
        //      **  **  A3  B4  C4
        //      **  **  **  A4  B5
        //
        //    An R83 matrix of order 5x3 would be stored:
        //
        //      B1  C1  **
        //      A1  B2  C2
        //      **  A2  B3
        //      **  **  A3
        //      **  **  **
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 February 2016
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
        //
        //    Input, double A1(N-1), A2(N), A3(N-1), the matrix.
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
        double[] ap = r83v_mv(n, n, a1, a2, a3, x);

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
            ap = r83v_mv(n, n, a1, a2, a3, p);
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

    public static void r83v_cg_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83V_CG_TEST tests R83V_CG.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 February 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("R83_CG_TEST");
        Console.WriteLine("  R83_CG applies CG to an R83 matrix.");

        int n = 10;
        //
        //  Let A be the -1 2 -1 matrix.
        //
        double[] a1 = new double[n - 1];
        double[] a2 = new double[n];
        double[] a3 = new double[n - 1];

        r83v_dif2(n, n, ref a1, ref a2, ref a3);
        //
        //  Choose a random solution.
        //
        int seed = 123456789;
        double[] x1 = UniformRNG.r8vec_uniform_01_new(n, ref seed);
        //
        //  Compute the corresponding right hand side.
        //
        double[] b = r83v_mv(n, n, a1, a2, a3, x1);
        //
        //  Call the CG routine.
        //
        double[] x2 = new double[n];
        for (i = 0; i < n; i++)
        {
            x2[i] = 1.0;
        }

        r83v_cg(n, a1, a2, a3, b, ref x2);
        //
        //  Compute the residual.
        //
        double[] r = r83v_res(n, n, a1, a2, a3, x2, b);
        double r_norm = r8vec_norm(n, r);
        //
        //  Compute the error.
        //
        double e_norm = r8vec_norm_affine(n, x1, x2);
        //
        //  Report.
        //
        Console.WriteLine("");
        Console.WriteLine("  Number of variables N = " + n + "");
        Console.WriteLine("  Norm of residual ||Ax-b|| = " + r_norm + "");
        Console.WriteLine("  Norm of error ||x1-x2|| = " + e_norm + "");
    }

    public static void r83v_copy(int m, int n, double[] a1, double[] a2, double[] a3,
            ref double[] b1, ref double[] b2, ref double[] b3)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83V_COPY copies a matrix in R83V format.
        //
        //  Discussion:
        //
        //    The R83V storage format is used for a tridiagonal matrix.
        //    The subdiagonal is in A(min(M-1,N)).
        //    The diagonal is in B(min(M,N)).
        //    The superdiagonal is in C(min(M,N-1)).
        //
        //  Example:
        //
        //    An R83V matrix of order 3x5 would be stored:
        //
        //      B1  C1  **  **  **
        //      A1  B2  C2  **  **
        //      **  A2  B3  C3  **
        //
        //    An R83 matrix of order 5x5 would be stored:
        //
        //      B1  C1  **  **  **
        //      A1  B2  C2  **  **
        //      **  A2  B3  C3  **
        //      **  **  A3  B4  C4
        //      **  **  **  A4  B5
        //
        //    An R83 matrix of order 5x3 would be stored:
        //
        //      B1  C1  **
        //      A1  B2  C2
        //      **  A2  B3
        //      **  **  A3
        //      **  **  **
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 February 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of the matrix.
        //
        //    Input, double A1(min(M-1,N)), A2(min(M,N)), A3(min(M,N-1)), the matrix.
        //
        //    Output, double B1(min(M-1,N)), B2(min(M,N)), B3(min(M,N-1)), the copy.
        //
    {
        int i;

        int ahi = Math.Min(m - 1, n);
        int bhi = Math.Min(m, n);
        int chi = Math.Min(m, n - 1);

        for (i = 0; i < ahi; i++)
        {
            b1[i] = a1[i];
        }

        for (i = 0; i < bhi; i++)
        {
            b2[i] = a2[i];
        }

        for (i = 0; i < chi; i++)
        {
            b3[i] = a3[i];
        }
    }

    public static void r83v_copy_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83V_COPY_TEST tests R83V_COPY.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 February 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int m = 0;
        int n = 0;

        Console.WriteLine("");
        Console.WriteLine("R83V_COPY_TEST");
        Console.WriteLine("  R83V_COPY copies an R83V matrix.");
        Console.WriteLine("  We check three cases, M<N, M=N, M>N.");

        for (i = 1; i <= 3; i++)
        {
            switch (i)
            {
                case 1:
                    m = 3;
                    n = 5;
                    break;
                case 2:
                    m = 5;
                    n = 5;
                    break;
                case 3:
                    m = 5;
                    n = 3;
                    break;
            }

            int ahi = Math.Min(m - 1, n);
            int bhi = Math.Min(m, n);
            int chi = Math.Min(m, n - 1);
            double[] a1 = new double[ahi];
            double[] a2 = new double[bhi];
            double[] a3 = new double[chi];

            r83v_indicator(m, n, ref a1, ref a2, ref a3);
            r83v_print(m, n, a1, a2, a3, "  R83V matrix A:");

            double[] b1 = new double[ahi];
            double[] b2 = new double[bhi];
            double[] b3 = new double[chi];
            r83v_copy(m, n, a1, a2, a3, ref b1, ref b2, ref b3);
            r83v_print(m, n, b1, b2, b3, "  B = copy of A:");

        }

    }

    public static double[] r83v_cr_fa(int n, double[] a1, double[] a2, double[] a3)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83V_CR_FA decomposes an R83V matrix using cyclic reduction.
        //
        //  Discussion:
        //
        //    The R83V storage format is used for a tridiagonal matrix.
        //    The subdiagonal is in A(min(M-1,N)).
        //    The diagonal is in B(min(M,N)).
        //    The superdiagonal is in C(min(M,N-1)).
        //
        //    Once R83V_CR_FA has decomposed a matrix A, then R83V_CR_SL may be used to
        //    solve linear systems A * x = b.
        //
        //    This function does not employ pivoting.  Hence, the results can be more
        //    sensitive to ill-conditioning than standard Gauss elimination.  In
        //    particular, this function will fail if any diagonal element of the matrix
        //    is zero.  Other matrices may also cause this function to fail.
        //
        //    This function can be guaranteed to work properly if the matrix is strictly
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
        //    defining this factorization may be used by R83V_CR_SL to solve linear
        //    systems.
        //
        //  Example:
        //
        //    An R83V matrix of order 3x5 would be stored:
        //
        //      B1  C1  **  **  **
        //      A1  B2  C2  **  **
        //      **  A2  B3  C3  **
        //
        //    An R83 matrix of order 5x5 would be stored:
        //
        //      B1  C1  **  **  **
        //      A1  B2  C2  **  **
        //      **  A2  B3  C3  **
        //      **  **  A3  B4  C4
        //      **  **  **  A4  B5
        //
        //    An R83 matrix of order 5x3 would be stored:
        //
        //      B1  C1  **
        //      A1  B2  C2
        //      **  A2  B3
        //      **  **  A3
        //      **  **  **
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 February 2016
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
        //
        //    Input, double A1(N-1), A2(N), A3(N-1), the matrix.
        //
        //    Output, double R83V_CR_FA[3*(2*N+1)], factorization information 
        //    needed by R83V_CR_SL.
        //
    {
        int j;

        switch (n)
        {
            case <= 0:
                Console.WriteLine("");
                Console.WriteLine("R83V_CR_FA - Fatal error!");
                Console.WriteLine("  Nonpositive N = " + n + "");
                return null;
        }

        double[] a_cr = new double[3 * (2 * n + 1)];

        for (j = 0; j < 2 * n + 1; j++)
        {
            int i;
            for (i = 0; i < 3; i++)
            {
                a_cr[i + j * 3] = 0.0;
            }
        }

        switch (n)
        {
            case 1:
                a_cr[1 + 0 * 3] = 1.0 / a2[0];
                return a_cr;
        }

        for (j = 1; j <= n - 1; j++)
        {
            a_cr[0 + j * 3] = a3[j - 1];
        }

        for (j = 1; j <= n; j++)
        {
            a_cr[1 + j * 3] = a2[j - 1];
        }

        for (j = 1; j <= n - 1; j++)
        {
            a_cr[2 + j * 3] = a1[j - 1];
        }

        int il = n;
        int ipntp = 0;

        while (1 < il)
        {
            int ipnt = ipntp;
            ipntp += il;
            int inc = (il % 2) switch
            {
                1 => il + 1,
                _ => il
            };

            int incr = inc / 2;
            il /= 2;
            int ihaf = ipntp + incr + 1;
            int ifulp = ipnt + inc + 2;

            int ilp;
            for (ilp = incr; 1 <= ilp; ilp--)
            {
                ifulp -= 2;
                int iful = ifulp - 1;
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

    public static void r83v_cr_fa_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83V_CR_FA_TEST tests R83V_CR_FA.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 February 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int j;
        const int n = 10;

        Console.WriteLine("");
        Console.WriteLine("R83V_CR_FA_TEST:");
        Console.WriteLine("  R83V_CR_FA factors an R83V matrix using cyclic reduction;");
        Console.WriteLine("");
        Console.WriteLine("  Matrix order N = " + n + "");
        Console.WriteLine("  The matrix is NOT symmetric.");
        //
        //  Set the matrix values.
        //
        double[] a1 = new double[n - 1];
        double[] a2 = new double[n];
        double[] a3 = new double[n - 1];

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = Math.Max(0, j - 1); i <= Math.Min(n - 1, j + 1); i++)
            {
                if (j == i - 1)
                {
                    a1[j] = j + 1;
                }
                else if (j == i)
                {
                    a2[j] = 4 * (j + 1);
                }
                else if (j == i + 1)
                {
                    a3[j - 1] = j + 1;
                }
            }
        }

        r83v_print(n, n, a1, a2, a3, "  The matrix:");
        //
        //  Set the desired solution.
        //
        double[] x = r8vec_indicator1_new(n);
        //
        //  Compute the corresponding right hand side.
        //
        double[] b = r83v_mv(n, n, a1, a2, a3, x);
        //
        //  Factor the matrix.
        //
        double[] a_cr = r83v_cr_fa(n, a1, a2, a3);
        //
        //  Solve the linear system.
        //
        x = r83v_cr_sl(n, a_cr, b);

        r8vec_print(n, x, "  Solution:");

    }

    public static double[] r83v_cr_sl(int n, double[] a_cr, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83V_CR_SL solves a real linear system factored by R83V_CR_FA.
        //
        //  Discussion:
        //
        //    The R83V storage format is used for a tridiagonal matrix.
        //    The subdiagonal is in A(min(M-1,N)).
        //    The diagonal is in B(min(M,N)).
        //    The superdiagonal is in C(min(M,N-1)).
        //
        //    The matrix A must be tridiagonal.  R83V_CR_FA is called to compute the
        //    LU factors of A.  It does so using a form of cyclic reduction.  If
        //    the factors computed by R83V_CR_FA are passed to R83V_CR_SL, then one or 
        //    many linear systems involving the matrix A may be solved.
        //
        //    Note that R83V_CR_FA does not perform pivoting, and so the solution 
        //    produced by R83V_CR_SL may be less accurate than a solution produced 
        //    by a standard Gauss algorithm.  However, such problems can be 
        //    guaranteed not to occur if the matrix A is strictly diagonally 
        //    dominant, that is, if the absolute value of the diagonal coefficient 
        //    is greater than the sum of the absolute values of the two off diagonal 
        //    coefficients, for each row of the matrix.
        //
        //  Example:
        //
        //    An R83V matrix of order 3x5 would be stored:
        //
        //      B1  C1  **  **  **
        //      A1  B2  C2  **  **
        //      **  A2  B3  C3  **
        //
        //    An R83 matrix of order 5x5 would be stored:
        //
        //      B1  C1  **  **  **
        //      A1  B2  C2  **  **
        //      **  A2  B3  C3  **
        //      **  **  A3  B4  C4
        //      **  **  **  A4  B5
        //
        //    An R83 matrix of order 5x3 would be stored:
        //
        //      B1  C1  **
        //      A1  B2  C2
        //      **  A2  B3
        //      **  **  A3
        //      **  **  **
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 February 2016
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
        //    Input, double A_CR[3*(2*N+1)], factorization information computed by
        //    R83V_CR_FA.
        //
        //    Input, double B[N], the right hand side.
        //
        //    Output, double R83V_CR_SL[N], the solution.
        //
    {
        int i;
        int iful;
        int ihaf = 0;
        int ipnt;
        double[] x;

        switch (n)
        {
            case <= 0:
                Console.WriteLine("");
                Console.WriteLine("R83V_CR_SL - Fatal error!");
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
        double[] rhs = new double[2 * n + 1];

        rhs[0] = 0.0;
        for (i = 1; i <= n; i++)
        {
            rhs[i] = b[i - 1];
        }

        for (i = n + 1; i <= 2 * n; i++)
        {
            rhs[i] = 0.0;
        }

        int il = n;
        int ndiv = 1;
        int ipntp = 0;

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

            int ifulm;
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

    public static void r83v_cr_sl_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83V_CR_SL_TEST tests R83V_CR_SL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 February 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int j;
        const int n = 10;

        Console.WriteLine("");
        Console.WriteLine("R83V_CR_SL_TEST:");
        Console.WriteLine("  R83V_CR_SL solves a linear system factored by R83V_CR_FA.");
        Console.WriteLine("");
        Console.WriteLine("  Matrix order N = " + n + "");
        Console.WriteLine("  The matrix is NOT symmetric.");
        //
        //  Set the matrix values.
        //
        double[] a1 = new double[n - 1];
        double[] a2 = new double[n];
        double[] a3 = new double[n - 1];

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = Math.Max(0, j - 1); i <= Math.Min(n - 1, j + 1); i++)
            {
                if (j == i - 1)
                {
                    a1[j] = j + 1;
                }
                else if (j == i)
                {
                    a2[j] = 4 * (j + 1);
                }
                else if (j == i + 1)
                {
                    a3[j - 1] = j + 1;
                }
            }
        }

        r83v_print(n, n, a1, a2, a3, "  The matrix:");
        //
        //  Set the desired solution.
        //
        double[] x = r8vec_indicator1_new(n);
        //
        //  Compute the corresponding right hand side.
        //
        double[] b = r83v_mv(n, n, a1, a2, a3, x);
        //
        //  Factor the matrix.
        //
        double[] a_cr = r83v_cr_fa(n, a1, a2, a3);
        //
        //  Solve the linear system.
        //
        x = r83v_cr_sl(n, a_cr, b);

        r8vec_print(n, x, "  Solution:");

    }

    public static double[] r83v_cr_sls(int n, double[] a_cr, int nb, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83V_CR_SLS solves several real linear systems factored by R83V_CR_FA.
        //
        //  Discussion:
        //
        //    The R83V storage format is used for a tridiagonal matrix.
        //    The subdiagonal is in A(min(M-1,N)).
        //    The diagonal is in B(min(M,N)).
        //    The superdiagonal is in C(min(M,N-1)).
        //
        //    The matrix A must be tridiagonal.  R83V_CR_FA is called to compute the
        //    LU factors of A.  It does so using a form of cyclic reduction.  If
        //    the factors computed by R83V_CR_FA are passed to R83V_CR_SLS, then one or 
        //    many linear systems involving the matrix A may be solved.
        //
        //    Note that R83V_CR_FA does not perform pivoting, and so the solutions
        //    produced by R83V_CR_SLS may be less accurate than a solution produced 
        //    by a standard Gauss algorithm.  However, such problems can be 
        //    guaranteed not to occur if the matrix A is strictly diagonally 
        //    dominant, that is, if the absolute value of the diagonal coefficient 
        //    is greater than the sum of the absolute values of the two off diagonal 
        //    coefficients, for each row of the matrix.
        //
        //  Example:
        //
        //    An R83V matrix of order 3x5 would be stored:
        //
        //      B1  C1  **  **  **
        //      A1  B2  C2  **  **
        //      **  A2  B3  C3  **
        //
        //    An R83 matrix of order 5x5 would be stored:
        //
        //      B1  C1  **  **  **
        //      A1  B2  C2  **  **
        //      **  A2  B3  C3  **
        //      **  **  A3  B4  C4
        //      **  **  **  A4  B5
        //
        //    An R83 matrix of order 5x3 would be stored:
        //
        //      B1  C1  **
        //      A1  B2  C2
        //      **  A2  B3
        //      **  **  A3
        //      **  **  **
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 February 2016
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
        //
        //    Input, double A_CR[3*(2*N+1)], factorization information computed by
        //    R83V_CR_FA.
        //
        //    Input, int NB, the number of systems.
        //
        //    Input, double B[N*NB], the right hand sides.
        //
        //    Output, double R83V_CR_SL[N*NB], the solutions.
        //
    {
        int i;
        int iful;
        int ihaf = 0;
        int ipnt;
        int j;
        double[] x;

        switch (n)
        {
            case <= 0:
                Console.WriteLine("");
                Console.WriteLine("R83V_CR_SLS - Fatal error!");
                Console.WriteLine("  Nonpositive N = " + n + "");
                return null;
            case 1:
            {
                x = new double[n * nb];
                for (j = 0; j < nb; j++)
                {
                    x[0 + j * n] = a_cr[1 + 0 * 3] * b[0 + j * n];
                }

                return x;
            }
        }

        //
        //  Set up RHS.
        //
        double[] rhs = new double[(2 * n + 1) * nb];

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

        int il = n;
        int ndiv = 1;
        int ipntp = 0;

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
                int ifulm;
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

    public static void r83v_cr_sls_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83V_CR_SLS_TEST tests R83V_CR_SLS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 February 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int j;
        const int n = 5;
        const int nb = 2;

        Console.WriteLine("");
        Console.WriteLine("R83V_CR_SLS_TEST");
        Console.WriteLine("  R83V_CR_SLS solves multiple linear systems A*x=b1:bn with R83V matrix");
        Console.WriteLine("  using cyclic reduction, after factorization by R83V_CR_FA.");
        Console.WriteLine("");
        Console.WriteLine("  Matrix order N = " + n + "");
        Console.WriteLine("  Number of linear systems = " + nb + "");
        Console.WriteLine("  Demonstrate multiple system solution method.");
        //
        //  Set the matrix values.
        //
        //
        //  Set the matrix values.
        //
        double[] a1 = new double[n - 1];
        double[] a2 = new double[n];
        double[] a3 = new double[n - 1];

        r83v_dif2(n, n, ref a1, ref a2, ref a3);

        r83v_print(n, n, a1, a2, a3, "  System matrix:");
        //
        //  Factor the matrix once.
        //
        double[] a_cr = r83v_cr_fa(n, a1, a2, a3);
        //
        //  Set up the linear systems.
        //
        double[] b = new double[n * nb];

        for (j = 0; j < nb; j++)
        {
            int i;
            for (i = 0; i < n; i++)
            {
                b[i + j * n] = 0.0;
            }
        }

        j = 0;
        b[n - 1 + j * n] = n + 1;

        j = 1;
        b[0 + j * n] = 1.0;
        b[n - 1 + j * n] = 1.0;

        r8ge_print(n, nb, b, "  RHS:");
        //
        //  Solve the linear systems.
        //
        double[] x = r83v_cr_sls(n, a_cr, nb, b);

        r8ge_print(n, nb, x, "  Solutions:");

    }

    public static void r83v_dif2(int m, int n, ref double[] a, ref double[] b, ref double[] c)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83V_DIF2 returns the DIF2 matrix in R83V format.
        //
        //  Discussion:
        //
        //    The R83V storage format is used for a tridiagonal matrix.
        //    The subdiagonal is in A(min(M-1,N)).
        //    The diagonal is in B(min(M,N)).
        //    The superdiagonal is in C(min(M,N-1)).
        //
        //  Example:
        //
        //    An R83V matrix of order 3x5 would be stored:
        //
        //      B1  C1  **  **  **
        //      A1  B2  C2  **  **
        //      **  A2  B3  C3  **
        //
        //    An R83 matrix of order 5x5 would be stored:
        //
        //      B1  C1  **  **  **
        //      A1  B2  C2  **  **
        //      **  A2  B3  C3  **
        //      **  **  A3  B4  C4
        //      **  **  **  A4  B5
        //
        //    An R83 matrix of order 5x3 would be stored:
        //
        //      B1  C1  **
        //      A1  B2  C2
        //      **  A2  B3
        //      **  **  A3
        //      **  **  **
        //
        //  Properties:
        //
        //    A is banded, with bandwidth 3.
        //    A is tridiagonal.
        //    Because A is tridiagonal, it has property A (bipartite).
        //    A is a special case of the TRIS or tridiagonal scalar matrix.
        //    A is integral, therefore det ( A ) is integral, and 
        //    det ( A ) * inverse ( A ) is integral.
        //    A is Toeplitz: constant along diagonals.
        //    A is symmetric: A' = A.
        //    Because A is symmetric, it is normal.
        //    Because A is normal, it is diagonalizable.
        //    A is persymmetric: A(I,J) = A(N+1-J,N+1-I).
        //    A is positive definite.
        //    A is an M matrix.
        //    A is weakly diagonally dominant, but not strictly diagonally dominant.
        //    A has an LU factorization A = L * U, without pivoting.
        //      The matrix L is lower bidiagonal with subdiagonal elements:
        //        L(I+1,I) = -I/(I+1)
        //      The matrix U is upper bidiagonal, with diagonal elements
        //        U(I,I) = (I+1)/I
        //      and superdiagonal elements which are all -1.
        //    A has a Cholesky factorization A = L * L', with L lower bidiagonal.
        //      L(I,I) =    sqrt ( (I+1) / I )
        //      L(I,I-1) = -sqrt ( (I-1) / I )
        //    The eigenvalues are
        //      LAMBDA(I) = 2 + 2 * COS(I*PI/(N+1))
        //                = 4 SIN^2(I*PI/(2*N+2))
        //    The corresponding eigenvector X(I) has entries
        //       X(I)(J) = sqrt(2/(N+1)) * sin ( I*J*PI/(N+1) ).
        //    Simple linear systems:
        //      x = (1,1,1,...,1,1),   A*x=(1,0,0,...,0,1)
        //      x = (1,2,3,...,n-1,n), A*x=(0,0,0,...,0,n+1)
        //    det ( A ) = N + 1.
        //    The value of the determinant can be seen by induction,
        //    and expanding the determinant across the first row:
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
        //    17 February 2016
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
        //    Output, double A(min(M-1,N)), B(min(M,N)), C(min(M,N-1)), the matrix.
        //
    {
        int i;

        int ahi = Math.Min(m - 1, n);
        int bhi = Math.Min(m, n);
        int chi = Math.Min(m, n - 1);

        for (i = 0; i < ahi; i++)
        {
            a[i] = -1.0;
        }

        for (i = 0; i < bhi; i++)
        {
            b[i] = 2.0;
        }

        for (i = 0; i < chi; i++)
        {
            c[i] = -1.0;
        }
    }

    public static void r83v_dif2_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83V_DIF2_TEST tests R83V_DIF2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 February 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int m = 0;
        int n = 0;

        Console.WriteLine("");
        Console.WriteLine("R83V_DIF2_TEST");
        Console.WriteLine("  R83V_DIF2 sets up an R83V second difference matrix.");
        Console.WriteLine("  We check three cases, M<N, M=N, M>N.");

        for (i = 1; i <= 3; i++)
        {
            switch (i)
            {
                case 1:
                    m = 3;
                    n = 5;
                    break;
                case 2:
                    m = 5;
                    n = 5;
                    break;
                case 3:
                    m = 5;
                    n = 3;
                    break;
            }

            int ahi = Math.Min(m - 1, n);
            int bhi = Math.Min(m, n);
            int chi = Math.Min(m, n - 1);
            double[] a = new double[ahi];
            double[] b = new double[bhi];
            double[] c = new double[chi];

            r83v_dif2(m, n, ref a, ref b, ref c);

            r83v_print(m, n, a, b, c, "  The R83V DIF2 matrix:");

        }

    }

    public static double[] r83v_fs(int n, double[] a1, double[] a2, double[] a3, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83V_FS solves a linear system with R83V matrix.
        //
        //  Discussion:
        //
        //    This function is based on the LINPACK SGTSL routine.
        // 
        //    The R83V storage format is used for a tridiagonal matrix.
        //    The subdiagonal is in A(min(M-1,N)).
        //    The diagonal is in B(min(M,N)).
        //    The superdiagonal is in C(min(M,N-1)).
        //
        //  Example:
        //
        //    An R83V matrix of order 3x5 would be stored:
        //
        //      B1  C1  **  **  **
        //      A1  B2  C2  **  **
        //      **  A2  B3  C3  **
        //
        //    An R83 matrix of order 5x5 would be stored:
        //
        //      B1  C1  **  **  **
        //      A1  B2  C2  **  **
        //      **  A2  B3  C3  **
        //      **  **  A3  B4  C4
        //      **  **  **  A4  B5
        //
        //    An R83 matrix of order 5x3 would be stored:
        //
        //      B1  C1  **
        //      A1  B2  C2
        //      **  A2  B3
        //      **  **  A3
        //      **  **  **
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 February 2016
        //
        //  Author:
        //
        //    John Burkardt, based on the LINPACK SGTSL function.
        //
        //  Reference:
        //
        //    Jack Dongarra, Cleve Moler, Jim Bunch and Pete Stewart,
        //    LINPACK User's Guide,
        //    SIAM, (Society for Industrial and Applied Mathematics),
        //    3600 University City Science Center,
        //    Philadelphia, PA, 19104-2688.
        //    ISBN 0-89871-172-X
        //
        //  Parameters:
        //
        //    Input, int N, the order of the tridiagonal matrix.
        //
        //    Input, double A1[N-1], A2[N], A3[N-1], the R83V matrix.
        //
        //    Input, double B[N], the right hand side.
        //
        //    Output, double R83V_FS[N], the solution.
        //
    {
        int i;
        int k;
        //
        //  Copy the input data.
        //
        double[] c = new double[n];
        double[] d = new double[n];
        double[] e = new double[n];
        double[] x = new double[n];

        c[0] = 0.0;
        for (i = 1; i < n; i++)
        {
            c[i] = a1[i - 1];
        }

        for (i = 0; i < n; i++)
        {
            d[i] = a2[i];
        }

        for (i = 0; i < n - 1; i++)
        {
            e[i] = a3[i];
        }

        e[n - 1] = 0.0;
        for (i = 0; i < n; i++)
        {
            x[i] = b[i];
        }

        //
        //  Factor.
        //
        c[0] = a2[0];

        switch (n)
        {
            case >= 2:
            {
                d[0] = e[0];
                e[0] = 0.0;
                e[n - 1] = 0.0;

                for (k = 1; k <= n - 1; k++)
                {
                    //
                    //  Find the larger of the two rows.
                    //
                    double t;
                    if (Math.Abs(c[k - 1]) <= Math.Abs(c[k]))
                    {
                        //
                        //  Interchange rows.
                        //
                        t = c[k];
                        c[k] = c[k - 1];
                        c[k - 1] = t;

                        t = d[k];
                        d[k] = d[k - 1];
                        d[k - 1] = t;

                        t = e[k];
                        e[k] = e[k - 1];
                        e[k - 1] = t;

                        t = x[k];
                        x[k] = x[k - 1];
                        x[k - 1] = t;
                    }

                    switch (c[k - 1])
                    {
                        //
                        //  Zero elements.
                        //
                        case 0.0:
                            return null;
                    }

                    t = -c[k] / c[k - 1];
                    c[k] = d[k] + t * d[k - 1];
                    d[k] = e[k] + t * e[k - 1];
                    e[k] = 0.0;
                    x[k] += t * x[k - 1];
                }

                break;
            }
        }

        switch (c[n - 1])
        {
            case 0.0:
                return null;
        }

        //
        //  Back solve.
        //
        x[n - 1] /= c[n - 1];

        switch (n)
        {
            case > 1:
            {
                x[n - 2] = (x[n - 2] - d[n - 2] * x[n - 1]) / c[n - 2];

                for (k = n - 2; 1 <= k; k--)
                {
                    x[k - 1] = (x[k - 1] - d[k - 1] * x[k] - e[k - 1] * x[k + 1]) / c[k - 1];
                }

                break;
            }
        }

        return x;
    }

    public static void r83v_fs_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83V_FS_TEST tests R83V_FS_SL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 February 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const int n = 10;

        Console.WriteLine("");
        Console.WriteLine("R83V_FS_TEST");
        Console.WriteLine("  R83V_FS factors and solves a linear system");
        Console.WriteLine("  for an R83V matrix.");
        Console.WriteLine("");
        Console.WriteLine("  Matrix order N = " + n + "");
        //
        //  Set the matrix values.
        //
        double[] a1 = new double[n - 1];
        double[] a2 = new double[n];
        double[] a3 = new double[n - 1];

        r83v_dif2(n, n, ref a1, ref a2, ref a3);
        //
        //  Set the desired solution.
        //
        double[] x1 = r8vec_indicator1_new(n);
        //
        //  Compute the corresponding right hand side.
        //
        double[] b = r83v_mv(n, n, a1, a2, a3, x1);

        r8vec_print(n, b, "  The right hand side:");
        //
        //  Solve the linear system.
        //
        double[] x2 = r83v_fs(n, a1, a2, a3, b);

        if (x2 == null || x2.Length == 0)
        {
            Console.WriteLine("");
            Console.WriteLine("  R83V_FS failed.");
        }
        else
        {
            r8vec_print(n, x2, "  Solution:");
        }

    }

    public static void r83v_gs_sl(int n, double[] a1, double[] a2, double[] a3, double[] b,
            ref double[] x, int it_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83V_GS_SL solves an R83V system using Gauss Seidel iteration.
        //
        //  Discussion:
        //
        //    The R83V storage format is used for a tridiagonal matrix.
        //    The subdiagonal is in A(min(M-1,N)).
        //    The diagonal is in B(min(M,N)).
        //    The superdiagonal is in C(min(M,N-1)).
        //
        //  Example:
        //
        //    An R83V matrix of order 3x5 would be stored:
        //
        //      B1  C1  **  **  **
        //      A1  B2  C2  **  **
        //      **  A2  B3  C3  **
        //
        //    An R83 matrix of order 5x5 would be stored:
        //
        //      B1  C1  **  **  **
        //      A1  B2  C2  **  **
        //      **  A2  B3  C3  **
        //      **  **  A3  B4  C4
        //      **  **  **  A4  B5
        //
        //    An R83 matrix of order 5x3 would be stored:
        //
        //      B1  C1  **
        //      A1  B2  C2
        //      **  A2  B3
        //      **  **  A3
        //      **  **  **
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 February 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A1(N-1), A2(N), A3(N-1), the R83V matrix.
        //
        //    Input, double B[N], the right hand side of the linear system.
        //
        //    Input/output, double X[N], an approximate solution to the system.
        //
        //    Input, int IT_MAX, the maximum number of iterations to take.
        //
    {
        int i;
        int it_num;
        //
        //  No diagonal matrix entry can be zero.
        //
        for (i = 0; i < n; i++)
        {
            switch (a2[i])
            {
                case 0.0:
                    Console.WriteLine("");
                    Console.WriteLine("R83V_GS_SL - Fatal error!");
                    Console.WriteLine("  Zero diagonal entry, index = " + i + "");
                    return;
            }
        }

        double[] x_new = new double[n];

        for (it_num = 1; it_num <= it_max; it_num++)
        {
            x_new[0] = (b[0] - a3[0] * x[1]) / a2[0];
            for (i = 1; i < n - 1; i++)
            {
                x_new[i] = (b[i] - a1[i - 1] * x_new[i - 1] - a3[i] * x[i + 1]) / a2[i];
            }

            x_new[n - 1] = (b[n - 1] - a1[n - 2] * x_new[n - 2]) / a2[n - 1];

            for (i = 0; i < n; i++)
            {
                x[i] = x_new[i];
            }
        }

    }

    public static void r83v_gs_sl_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83V_GS_SL_TEST tests R83V_GS_SL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 February 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        const int maxit = 25;
        const int n = 10;

        Console.WriteLine("");
        Console.WriteLine("R83V_GS_SL_TEST");
        Console.WriteLine("  R83V_GS_SL solves a linear system using Gauss-Seidel");
        Console.WriteLine("  iteration for an R83V matrix.");
        Console.WriteLine("");
        Console.WriteLine("  Matrix order N = " + n + "");
        Console.WriteLine("  Iterations per call = " + maxit + "");
        //
        //  Set the matrix values.
        //
        double[] a1 = new double[n - 1];
        double[] a2 = new double[n];
        double[] a3 = new double[n - 1];

        r83v_dif2(n, n, ref a1, ref a2, ref a3);
        //
        //  Set the desired solution.
        //
        double[] x = r8vec_indicator1_new(n);
        //
        //  Compute the corresponding right hand side.
        //
        double[] b = r83v_mv(n, n, a1, a2, a3, x);

        r8vec_print(n, b, "  The right hand side:");
        //
        //  Set the starting solution.
        //
        for (i = 0; i < n; i++)
        {
            x[i] = 0.0;
        }

        //
        //  Solve the linear system.
        //
        for (i = 1; i <= 3; i++)
        {
            r83v_gs_sl(n, a1, a2, a3, b, ref x, maxit);

            r8vec_print(n, x, "  Current estimated solution:");
        }

    }

    public static void r83v_indicator(int m, int n, ref double[] a, ref double[] b, ref double[] c)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83V_INDICATOR sets up an R83V indicator matrix.
        //
        //  Discussion:
        //
        //    The R83V storage format is used for a tridiagonal matrix.
        //    The subdiagonal is in A(min(M-1,N)).
        //    The diagonal is in B(min(M,N)).
        //    The superdiagonal is in C(min(M,N-1)).
        //
        //  Example:
        //
        //    An R83V matrix of order 3x5 would be stored:
        //
        //      B1  C1  **  **  **
        //      A1  B2  C2  **  **
        //      **  A2  B3  C3  **
        //
        //    An R83 matrix of order 5x5 would be stored:
        //
        //      B1  C1  **  **  **
        //      A1  B2  C2  **  **
        //      **  A2  B3  C3  **
        //      **  **  A3  B4  C4
        //      **  **  **  A4  B5
        //
        //    An R83 matrix of order 5x3 would be stored:
        //
        //      B1  C1  **
        //      A1  B2  C2
        //      **  A2  B3
        //      **  **  A3
        //      **  **  **
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 February 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of the matrix.
        //
        //    Output, double A(min(M-1,N)), B(min(M,N)), C(min(M,N-1)), the matrix.
        //
    {
        int i;

        int fac = (int) Math.Pow(10, (int) Math.Log10(n) + 1);

        int ahi = Math.Min(m - 1, n);
        int bhi = Math.Min(m, n);
        int chi = Math.Min(m, n - 1);

        for (i = 0; i < ahi; i++)
        {
            a[i] = fac * (i + 2) + i + 1;
        }

        for (i = 0; i < bhi; i++)
        {
            b[i] = fac * (i + 1) + i + 1;
        }

        for (i = 0; i < chi; i++)
        {
            c[i] = fac * (i + 1) + i + 2;
        }
    }

    public static void r83v_indicator_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83V_INDICATOR_TEST tests R83V_INDICATOR.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 February 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int m = 0;
        int n = 0;

        Console.WriteLine("");
        Console.WriteLine("R83V_INDICATOR_TEST");
        Console.WriteLine("  R83V_INDICATOR sets up an R83V indicator matrix.");
        Console.WriteLine("  We check three cases, M<N, M=N, M>N.");

        for (i = 1; i <= 3; i++)
        {
            switch (i)
            {
                case 1:
                    m = 3;
                    n = 5;
                    break;
                case 2:
                    m = 5;
                    n = 5;
                    break;
                case 3:
                    m = 5;
                    n = 3;
                    break;
            }

            int ahi = Math.Min(m - 1, n);
            int bhi = Math.Min(m, n);
            int chi = Math.Min(m, n - 1);
            double[] a = new double[ahi];
            double[] b = new double[bhi];
            double[] c = new double[chi];

            r83v_indicator(m, n, ref a, ref b, ref c);

            r83v_print(m, n, a, b, c, "  The R83V indicator matrix:");

        }

    }

    public static void r83v_jac_sl(int n, double[] a1, double[] a2, double[] a3, double[] b,
            ref double[] x, int it_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83V_JAC_SL solves an R83V system using Jacobi iteration.
        //
        //  Discussion:
        //
        //    The R83V storage format is used for a tridiagonal matrix.
        //    The subdiagonal is in A(min(M-1,N)).
        //    The diagonal is in B(min(M,N)).
        //    The superdiagonal is in C(min(M,N-1)).
        //
        //  Example:
        //
        //    An R83V matrix of order 3x5 would be stored:
        //
        //      B1  C1  **  **  **
        //      A1  B2  C2  **  **
        //      **  A2  B3  C3  **
        //
        //    An R83 matrix of order 5x5 would be stored:
        //
        //      B1  C1  **  **  **
        //      A1  B2  C2  **  **
        //      **  A2  B3  C3  **
        //      **  **  A3  B4  C4
        //      **  **  **  A4  B5
        //
        //    An R83 matrix of order 5x3 would be stored:
        //
        //      B1  C1  **
        //      A1  B2  C2
        //      **  A2  B3
        //      **  **  A3
        //      **  **  **
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 February 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A1(N-1), A2(N), A3(N-1), the R83V matrix.
        //
        //    Input, double B[N], the right hand side of the linear system.
        //
        //    Input/output, double X[N], an approximate solution to the system.
        //
        //    Input, int IT_MAX, the maximum number of iterations to take.
        //
    {
        int i;
        int it_num;

        double[] x_new = new double[n];
        //
        //  No diagonal matrix entry can be zero.
        //
        for (i = 0; i < n; i++)
        {
            switch (a2[i])
            {
                case 0.0:
                    Console.WriteLine("");
                    Console.WriteLine("R83V_JAC_SL - Fatal error!");
                    Console.WriteLine("  Zero diagonal entry, index = " + i + "");
                    return;
            }
        }

        for (it_num = 1; it_num <= it_max; it_num++)
        {
            //
            //  Solve A*x=b:
            //
            for (i = 0; i < n; i++)
            {
                x_new[i] = b[i];
            }

            for (i = 0; i < n - 1; i++)
            {
                x_new[i] -= a3[i] * x[i + 1];
            }

            for (i = 0; i < n - 1; i++)
            {
                x_new[i + 1] -= a1[i] * x[i];
            }

            //
            //  Divide by the diagonal term, and overwrite X.
            //
            for (i = 0; i < n; i++)
            {
                x[i] = x_new[i] / a2[i];
            }
        }
    }

    public static void r83v_jac_sl_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83V_JAC_SL_TEST tests R83V_JAC_SL.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 February 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        const int maxit = 25;
        const int n = 10;

        Console.WriteLine("");
        Console.WriteLine("R83V_JAC_SL_TEST");
        Console.WriteLine("  R83V_JAC_SL solves a linear system using Jacobi iteration,");
        Console.WriteLine("  for an R83V matrix.");
        Console.WriteLine("");
        Console.WriteLine("  Matrix order N = " + n + "");
        Console.WriteLine("  Iterations per call = " + maxit + "");
        //
        //  Set the matrix values.
        //
        double[] a1 = new double[n - 1];
        double[] a2 = new double[n];
        double[] a3 = new double[n - 1];

        r83v_dif2(n, n, ref a1, ref a2, ref a3);
        //
        //  Set the desired solution.
        //
        double[] x = r8vec_indicator1_new(n);
        //
        //  Compute the corresponding right hand side.
        //
        double[] b = r83v_mv(n, n, a1, a2, a3, x);

        r8vec_print(n, b, "  The right hand side:");
        //
        //  Set the starting solution.
        //
        for (i = 0; i < n; i++)
        {
            x[i] = 0.0;
        }

        //
        //  Solve the linear system.
        //
        for (i = 1; i <= 3; i++)
        {
            r83v_jac_sl(n, a1, a2, a3, b, ref x, maxit);

            r8vec_print(n, x, "  Current estimated solution:");
        }
    }

    public static double[] r83v_mtv(int m, int n, double[] a1, double[] a2, double[] a3,
            double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83V_MTV multiplies a vector times an R83V matrix.
        //
        //  Discussion:
        //
        //    The R83V storage format is used for a tridiagonal matrix.
        //    The subdiagonal is in A(min(M-1,N)).
        //    The diagonal is in B(min(M,N)).
        //    The superdiagonal is in C(min(M,N-1)).
        //
        //  Example:
        //
        //    An R83V matrix of order 3x5 would be stored:
        //
        //      B1  C1  **  **  **
        //      A1  B2  C2  **  **
        //      **  A2  B3  C3  **
        //
        //    An R83 matrix of order 5x5 would be stored:
        //
        //      B1  C1  **  **  **
        //      A1  B2  C2  **  **
        //      **  A2  B3  C3  **
        //      **  **  A3  B4  C4
        //      **  **  **  A4  B5
        //
        //    An R83 matrix of order 5x3 would be stored:
        //
        //      B1  C1  **
        //      A1  B2  C2
        //      **  A2  B3
        //      **  **  A3
        //      **  **  **
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 February 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of the linear system.
        //
        //    Input, double A1(min(M-1,N)), A2(min(M,N)), A3(min(M,N-1)), the matrix.
        //
        //    Input, double X[M], the vector to be multiplied.
        //
        //    Output, double R83V_MTV[N], the product A'*x.
        //
    {
        int j;

        double[] b = new double[n];

        for (j = 0; j < n; j++)
        {
            b[j] = 0.0;
        }

        int ahi = Math.Min(m - 1, n);
        int bhi = Math.Min(m, n);
        int chi = Math.Min(m, n - 1);

        for (j = 0; j < ahi; j++)
        {
            b[j] += a1[j] * x[j + 1];
        }

        for (j = 0; j < bhi; j++)
        {
            b[j] += a2[j] * x[j];
        }

        for (j = 0; j < chi; j++)
        {
            b[j + 1] += a3[j] * x[j];
        }

        return b;
    }

    public static void r83v_mtv_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83V_MTV_TEST tests R83V_MTV.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 February 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int m = 0;
        int n = 0;

        Console.WriteLine("");
        Console.WriteLine("R83V_MTV_TEST");
        Console.WriteLine("  R83V_MTV computes b=A'*x, where A is an R83V matrix.");
        Console.WriteLine("  We check three cases, M<N, M=N, M>N.");

        for (i = 1; i <= 3; i++)
        {
            switch (i)
            {
                case 1:
                    m = 3;
                    n = 5;
                    break;
                case 2:
                    m = 5;
                    n = 5;
                    break;
                case 3:
                    m = 5;
                    n = 3;
                    break;
            }

            int ahi = Math.Min(m - 1, n);
            int bhi = Math.Min(m, n);
            int chi = Math.Min(m, n - 1);

            double[] a1 = new double[ahi];
            double[] a2 = new double[bhi];
            double[] a3 = new double[chi];

            int seed = 123456789;
            r83v_random(m, n, ref seed, ref a1, ref a2, ref a3);
            double[] x = r8vec_indicator1_new(m);
            double[] ax = r83v_mtv(m, n, a1, a2, a3, x);

            double[] a_ge = r83v_to_r8ge(m, n, a1, a2, a3);
            double[] ax_ge = r8ge_mtv(m, n, a_ge, x);
            r8vec2_print(n, ax, ax_ge, "  Product comparison:");

        }
    }

    public static double[] r83v_mv(int m, int n, double[] a, double[] b, double[] c, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83V_MV multiplies an R83V matrix times a vector.
        //
        //  Discussion:
        //
        //    The R83V storage format is used for a tridiagonal matrix.
        //    The subdiagonal is in A(min(M-1,N)).
        //    The diagonal is in B(min(M,N)).
        //    The superdiagonal is in C(min(M,N-1)).
        //
        //  Example:
        //
        //    An R83V matrix of order 3x5 would be stored:
        //
        //      B1  C1  **  **  **
        //      A1  B2  C2  **  **
        //      **  A2  B3  C3  **
        //
        //    An R83 matrix of order 5x5 would be stored:
        //
        //      B1  C1  **  **  **
        //      A1  B2  C2  **  **
        //      **  A2  B3  C3  **
        //      **  **  A3  B4  C4
        //      **  **  **  A4  B5
        //
        //    An R83 matrix of order 5x3 would be stored:
        //
        //      B1  C1  **
        //      A1  B2  C2
        //      **  A2  B3
        //      **  **  A3
        //      **  **  **
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 February 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of the linear system.
        //
        //    Input, double A(min(M-1,N)), B(min(M,N)), C(min(M,N-1)), the R83V matrix.
        //
        //    Input, double X[N], the vector to be multiplied by A.
        //
        //    Output, double R83V_MV[M], the product A * x.
        //
    {
        int i;

        double[] ax = new double[m];

        for (i = 0; i < m; i++)
        {
            ax[i] = 0.0;
        }

        int ahi = Math.Min(m - 1, n);
        int bhi = Math.Min(m, n);
        int chi = Math.Min(m, n - 1);

        for (i = 0; i < ahi; i++)
        {
            ax[i + 1] += a[i] * x[i];
        }

        for (i = 0; i < bhi; i++)
        {
            ax[i] += b[i] * x[i];
        }

        for (i = 0; i < chi; i++)
        {
            ax[i] += c[i] * x[i + 1];
        }

        return ax;
    }

    public static void r83v_mv_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83V_MV_TEST tests R83_MV.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 February 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int m = 0;
        int n = 0;

        Console.WriteLine("");
        Console.WriteLine("R83V_MV_TEST");
        Console.WriteLine("  R83V_MV computes b=A*x, where A is an R83V matrix.");
        Console.WriteLine("  We check three cases, M<N, M=N, M>N.");

        for (i = 1; i <= 3; i++)
        {
            switch (i)
            {
                case 1:
                    m = 3;
                    n = 5;
                    break;
                case 2:
                    m = 5;
                    n = 5;
                    break;
                case 3:
                    m = 5;
                    n = 3;
                    break;
            }

            int ahi = Math.Min(m - 1, n);
            int bhi = Math.Min(m, n);
            int chi = Math.Min(m, n - 1);

            double[] a = new double[ahi];
            double[] b = new double[bhi];
            double[] c = new double[chi];

            int seed = 123456789;
            r83v_random(m, n, ref seed, ref a, ref b, ref c);
            double[] x = r8vec_indicator1_new(n);
            double[] ax = r83v_mv(m, n, a, b, c, x);

            double[] a_ge = r83v_to_r8ge(m, n, a, b, c);
            double[] ax_ge = r8ge_mv(m, n, a_ge, x);
            r8vec2_print(m, ax, ax_ge, "  Product comparison:");

        }
    }

    public static void r83v_print(int m, int n, double[] a, double[] b, double[] c,
            string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83V_PRINT prints an R83V matrix.
        //
        //  Discussion:
        //
        //    The R83V storage format is used for a tridiagonal matrix.
        //    The subdiagonal is in A(min(M-1,N)).
        //    The diagonal is in B(min(M,N)).
        //    The superdiagonal is in C(min(M,N-1)).
        //
        //  Example:
        //
        //    An R83V matrix of order 3x5 would be stored:
        //
        //      B1  C1  **  **  **
        //      A1  B2  C2  **  **
        //      **  A2  B3  C3  **
        //
        //    An R83 matrix of order 5x5 would be stored:
        //
        //      B1  C1  **  **  **
        //      A1  B2  C2  **  **
        //      **  A2  B3  C3  **
        //      **  **  A3  B4  C4
        //      **  **  **  A4  B5
        //
        //    An R83 matrix of order 5x3 would be stored:
        //
        //      B1  C1  **
        //      A1  B2  C2
        //      **  A2  B3
        //      **  **  A3
        //      **  **  **
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 February 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of the matrix.
        //
        //    Input, double A(min(M-1,N)), B(min(M,N)), C(min(M,N-1)), the R83V matrix.
        //
        //    Input, string TITLE, a title.
        //
    {
        r83v_print_some(m, n, a, b, c, 1, 1, m, n, title);
    }

    public static void r83v_print_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83V_PRINT_TEST tests R83V_PRINT.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 February 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("R83V_PRINT_TEST");
        Console.WriteLine("  R83V_PRINT prints an R83V matrix.");

        const int m = 5;
        const int n = 5;

        int ahi = Math.Min(m - 1, n);
        int bhi = Math.Min(m, n);
        int chi = Math.Min(m, n - 1);
        double[] a = new double[ahi];
        double[] b = new double[bhi];
        double[] c = new double[chi];

        r83v_indicator(m, n, ref a, ref b, ref c);

        r83v_print(m, n, a, b, c, "  The R83V  matrix:");

    }

    public static void r83v_print_some(int m, int n, double[] a, double[] b, double[] c,
            int ilo, int jlo, int ihi, int jhi, string title)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83V_PRINT_SOME prints some of an R83V matrix.
        //
        //  Discussion:
        //
        //    The R83V storage format is used for a tridiagonal matrix.
        //    The subdiagonal is in A(min(M-1,N)).
        //    The diagonal is in B(min(M,N)).
        //    The superdiagonal is in C(min(M,N-1)).
        //
        //  Example:
        //
        //    An R83V matrix of order 3x5 would be stored:
        //
        //      B1  C1  **  **  **
        //      A1  B2  C2  **  **
        //      **  A2  B3  C3  **
        //
        //    An R83 matrix of order 5x5 would be stored:
        //
        //      B1  C1  **  **  **
        //      A1  B2  C2  **  **
        //      **  A2  B3  C3  **
        //      **  **  A3  B4  C4
        //      **  **  **  A4  B5
        //
        //    An R83 matrix of order 5x3 would be stored:
        //
        //      B1  C1  **
        //      A1  B2  C2
        //      **  A2  B3
        //      **  **  A3
        //      **  **  **
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 February 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of the matrix.
        //
        //    Input, double A(min(M-1,N)), B(min(M,N)), C(min(M+1,N)), the R83V matrix.
        //
        //    Input, int ILO, JLO, IHI, JHI, designate the first row and
        //    column, and the last row and column, to be printed.
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

            int inc = j2hi + 1 - j2lo;

            Console.WriteLine("");
            string cout = "  Col: ";
            int j2;
            int j;
            for (j = j2lo; j <= j2hi; j++)
            {
                j2 = j + 1 - j2lo;
                cout += j.ToString().PadLeft(7) + "       ";
            }

            Console.WriteLine(cout);
            Console.WriteLine("  Row");
            Console.WriteLine("  ---");
            //
            //  Determine the range of the rows in this strip.
            //
            int i2lo = Math.Max(ilo, 1);
            i2lo = Math.Max(i2lo, j2lo - 1);

            int i2hi = Math.Min(ihi, m);
            i2hi = Math.Min(i2hi, j2hi + 1);

            int i;
            for (i = i2lo; i <= i2hi; i++)
            {
                //
                //  Print out (up to) 5 entries in row I, that lie in the current strip.
                //
                cout = i.ToString().PadLeft(6) + "  ";

                for (j2 = 1; j2 <= inc; j2++)
                {
                    j = j2lo - 1 + j2;

                    switch (i - j + 1)
                    {
                        case < 0:
                        case > 2:
                            cout += "              ";
                            break;
                        default:
                        {
                            if (j == i - 1)
                            {
                                cout += "  " + a[i - 2].ToString().PadLeft(12);
                            }
                            else if (j == i)
                            {
                                cout += "  " + b[i - 1].ToString().PadLeft(12);
                            }
                            else if (j == i + 1)
                            {
                                cout += "  " + c[i - 1].ToString().PadLeft(12);
                            }

                            break;
                        }
                    }
                }

                Console.WriteLine(cout);
            }
        }

    }

    public static void r83v_print_some_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83V_PRINT_SOME_TEST tests R83V_PRINT_SOME.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 February 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("R83V_PRINT_SOME_TEST");
        Console.WriteLine("  R83V_PRINT_SOME prints some of an R83V matrix.");

        const int m = 5;
        const int n = 5;

        int ahi = Math.Min(m - 1, n);
        int bhi = Math.Min(m, n);
        int chi = Math.Min(m, n - 1);
        double[] a = new double[ahi];
        double[] b = new double[bhi];
        double[] c = new double[chi];

        r83v_indicator(m, n, ref a, ref b, ref c);

        r83v_print_some(m, n, a, b, c, 2, 2, 5, 4, "  Rows 2-5, Cols 2-4:");

    }

    public static void r83v_random(int m, int n, ref int seed, ref double[] a, ref double[] b, ref double[] c)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83V_RANDOM returns a random matrix in R83V format.
        //
        //  Discussion:
        //
        //    The R83V storage format is used for a tridiagonal matrix.
        //    The subdiagonal is in A(min(M-1,N)).
        //    The diagonal is in B(min(M,N)).
        //    The superdiagonal is in C(min(M,N-1)).
        //
        //  Example:
        //
        //    An R83V matrix of order 3x5 would be stored:
        //
        //      B1  C1  **  **  **
        //      A1  B2  C2  **  **
        //      **  A2  B3  C3  **
        //
        //    An R83 matrix of order 5x5 would be stored:
        //
        //      B1  C1  **  **  **
        //      A1  B2  C2  **  **
        //      **  A2  B3  C3  **
        //      **  **  A3  B4  C4
        //      **  **  **  A4  B5
        //
        //    An R83 matrix of order 5x3 would be stored:
        //
        //      B1  C1  **
        //      A1  B2  C2
        //      **  A2  B3
        //      **  **  A3
        //      **  **  **
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 February 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of the matrix.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double A(min(M-1,N)), B(min(M,N)), C(min(M,N-1)), the matrix.
        //
    {
        int ahi = Math.Min(m - 1, n);
        int bhi = Math.Min(m, n);
        int chi = Math.Min(m, n - 1);

        UniformRNG.r8vec_uniform_01(ahi, ref seed, ref a);
        UniformRNG.r8vec_uniform_01(bhi, ref seed, ref b);
        UniformRNG.r8vec_uniform_01(chi, ref seed, ref c);

    }

    public static void r83v_random_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83V_RANDOM_TEST tests R83V_RANDOM.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 February 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int m = 0;
        int n = 0;

        Console.WriteLine("");
        Console.WriteLine("R83V_RANDOM_TEST");
        Console.WriteLine("  R83V_RANDOM sets up an R83V random matrix.");
        Console.WriteLine("  We check three cases, M<N, M=N, M>N.");

        for (i = 1; i <= 3; i++)
        {
            int seed = 123456789;

            switch (i)
            {
                case 1:
                    m = 3;
                    n = 5;
                    break;
                case 2:
                    m = 5;
                    n = 5;
                    break;
                case 3:
                    m = 5;
                    n = 3;
                    break;
            }

            int ahi = Math.Min(m - 1, n);
            int bhi = Math.Min(m, n);
            int chi = Math.Min(m, n - 1);
            double[] a = new double[ahi];
            double[] b = new double[bhi];
            double[] c = new double[chi];

            r83v_random(m, n, ref seed, ref a, ref b, ref c);

            r83v_print(m, n, a, b, c, "  The R83V random matrix:");

        }

    }

    public static double[] r83v_res(int m, int n, double[] a, double[] b, double[] c,
            double[] x, double[] ax)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83V_RES computes the residual R = B-A*X for R83V matrices.
        //
        //  Discussion:
        //
        //    The R83V storage format is used for a tridiagonal matrix.
        //    The subdiagonal is in A(min(M-1,N)).
        //    The diagonal is in B(min(M,N)).
        //    The superdiagonal is in C(min(M,N-1)).
        //
        //  Example:
        //
        //    An R83V matrix of order 3x5 would be stored:
        //
        //      B1  C1  **  **  **
        //      A1  B2  C2  **  **
        //      **  A2  B3  C3  **
        //
        //    An R83 matrix of order 5x5 would be stored:
        //
        //      B1  C1  **  **  **
        //      A1  B2  C2  **  **
        //      **  A2  B3  C3  **
        //      **  **  A3  B4  C4
        //      **  **  **  A4  B5
        //
        //    An R83 matrix of order 5x3 would be stored:
        //
        //      B1  C1  **
        //      A1  B2  C2
        //      **  A2  B3
        //      **  **  A3
        //      **  **  **
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 February 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the number of rows of the matrix.
        //
        //    Input, int N, the number of columns of the matrix.
        //
        //    Input, double A(min(M-1,N)), B(min(M,N)), C(min(M,N-1)), the matrix.
        //
        //    Input, double X[N], the vector to be multiplied.
        //
        //    Input, double AX[M], the desired result A * x.
        //
        //    Output, double R83V_RES[M], the residual R = AX - A * X.
        //
    {
        int i;

        double[] r = r83v_mv(m, n, a, b, c, x);

        for (i = 0; i < m; i++)
        {
            r[i] = ax[i] - r[i];
        }

        return r;
    }

    public static void r83v_res_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83V_RES_TEST tests R83V_RES.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 February 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int m = 0;
        int n = 0;

        Console.WriteLine("");
        Console.WriteLine("R83V_RES_TEST");
        Console.WriteLine("  R83V_RES computes b-A*x, where A is an R83V matrix.");
        Console.WriteLine("  We check three cases, M<N, M=N, M>N.");

        for (i = 1; i <= 3; i++)
        {
            switch (i)
            {
                case 1:
                    m = 3;
                    n = 5;
                    break;
                case 2:
                    m = 5;
                    n = 5;
                    break;
                case 3:
                    m = 5;
                    n = 3;
                    break;
            }

            int ahi = Math.Min(m - 1, n);
            int bhi = Math.Min(m, n);
            int chi = Math.Min(m, n - 1);
            double[] a = new double[ahi];
            double[] b = new double[bhi];
            double[] c = new double[chi];

            int seed = 123456789;
            r83v_random(m, n, ref seed, ref a, ref b, ref c);
            double[] x = r8vec_indicator1_new(n);
            double[] ax = r83v_mv(m, n, a, b, c, x);
            double[] r = r83v_res(m, n, a, b, c, x, ax);
            r8vec_print(m, r, "  Residual A*x-b:");

        }

    }

    public static double[] r83v_to_r8ge(int m, int n, double[] a1, double[] a2, double[] a3)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83V_TO_R8GE copies an R83V matrix to an R8GE matrix.
        //
        //  Discussion:
        //
        //    The R83V storage format is used for a tridiagonal matrix.
        //    The subdiagonal is in A(min(M-1,N)).
        //    The diagonal is in B(min(M,N)).
        //    The superdiagonal is in C(min(M,N-1)).
        //
        //  Example:
        //
        //    An R83V matrix of order 3x5 would be stored:
        //
        //      B1  C1  **  **  **
        //      A1  B2  C2  **  **
        //      **  A2  B3  C3  **
        //
        //    An R83 matrix of order 5x5 would be stored:
        //
        //      B1  C1  **  **  **
        //      A1  B2  C2  **  **
        //      **  A2  B3  C3  **
        //      **  **  A3  B4  C4
        //      **  **  **  A4  B5
        //
        //    An R83 matrix of order 5x3 would be stored:
        //
        //      B1  C1  **
        //      A1  B2  C2
        //      **  A2  B3
        //      **  **  A3
        //      **  **  **
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 February 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of the matrix.
        //
        //    Input, double A1(min(M-1,N)), A2(min(M,N)), A3(min(M,N-1)), the matrix.
        //
        //    Output, double R83V_TO_R8GE(M,N), the R8GE matrix.
        //
    {
        int j;
        int k;

        double[] a = new double[m * n];

        for (j = 0; j < n; j++)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                a[i + j * m] = 0.0;
            }
        }

        int ahi = Math.Min(m - 1, n);
        int bhi = Math.Min(m, n);
        int chi = Math.Min(m, n - 1);

        for (k = 0; k < ahi; k++)
        {
            a[k + 1 + k * m] = a1[k];
        }

        for (k = 0; k < bhi; k++)
        {
            a[k + k * m] = a2[k];
        }

        for (k = 0; k < chi; k++)
        {
            a[k + (k + 1) * m] = a3[k];
        }

        return a;
    }

    public static void r83v_to_r8ge_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83V_TO_R8GE_TEST tests R83V_TO_R8GE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 February 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int m = 0;
        int n = 0;

        Console.WriteLine("");
        Console.WriteLine("R83V_TO_R8GE_TEST");
        Console.WriteLine("  R83V_TO_R8GE copies an R83V matrix to an R8GE matrix.");
        Console.WriteLine("  We check three cases, M<N, M=N, M>N.");

        for (i = 1; i <= 3; i++)
        {
            switch (i)
            {
                case 1:
                    m = 3;
                    n = 5;
                    break;
                case 2:
                    m = 5;
                    n = 5;
                    break;
                case 3:
                    m = 5;
                    n = 3;
                    break;
            }

            int ahi = Math.Min(m - 1, n);
            int bhi = Math.Min(m, n);
            int chi = Math.Min(m, n - 1);
            double[] a1 = new double[ahi];
            double[] a2 = new double[bhi];
            double[] a3 = new double[chi];

            r83v_indicator(m, n, ref a1, ref a2, ref a3);
            r83v_print(m, n, a1, a2, a3, "  R83V matrix A:");

            double[] a = r83v_to_r8ge(m, n, a1, a2, a3);
            r8ge_print(m, n, a, "  R8GE version of A:");

        }

    }

    public static double[] r83v_to_r8vec(int m, int n, double[] a1, double[] a2, double[] a3)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83V_TO_R8VEC copies an R83V matrix to an R8VEC.
        //
        //  Discussion:
        //
        //    The R83V storage format is used for a tridiagonal matrix.
        //    The subdiagonal is in A(min(M-1,N)).
        //    The diagonal is in B(min(M,N)).
        //    The superdiagonal is in C(min(M,N-1)).
        //
        //  Example:
        //
        //    An R83V matrix of order 3x5 would be stored:
        //
        //      B1  C1  **  **  **
        //      A1  B2  C2  **  **
        //      **  A2  B3  C3  **
        //
        //    An R83 matrix of order 5x5 would be stored:
        //
        //      B1  C1  **  **  **
        //      A1  B2  C2  **  **
        //      **  A2  B3  C3  **
        //      **  **  A3  B4  C4
        //      **  **  **  A4  B5
        //
        //    An R83 matrix of order 5x3 would be stored:
        //
        //      B1  C1  **
        //      A1  B2  C2
        //      **  A2  B3
        //      **  **  A3
        //      **  **  **
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 February 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of the matrix.
        //
        //    Input, double A1(min(M-1,N)), A2(min(M,N)), A3(min(M,N-1)), the matrix.
        //
        //    Output, double R83V_TO_R8VEC(min(N-1,M)+min(N,M)+min(N,M-1)), the vector.
        //
    {
        int j;

        int ahi = Math.Min(m - 1, n);
        int bhi = Math.Min(m, n);
        int chi = Math.Min(m, n - 1);

        double[] a = new double[ahi + bhi + chi];

        int k = 0;
        for (j = 0; j < n; j++)
        {
            if (j < m + 1 && 1 <= j)
            {
                a[k] = a3[j - 1];
                k += 1;
            }

            if (j < m)
            {
                a[k] = a2[j];
                k += 1;
            }

            if (j < m - 1)
            {
                a[k] = a1[j];
                k += 1;
            }
        }

        return a;
    }

    public static void r83v_to_r8vec_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83V_TO_R8VEC_TEST tests R83V_TO_R8VEC.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 February 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int m = 0;
        int n = 0;

        Console.WriteLine("");
        Console.WriteLine("R83V_TO_R8VEC_TEST");
        Console.WriteLine("  R83V_TO_R8VEC copies an R83V matrix to an R8VEC.");
        Console.WriteLine("  We check three cases, M<N, M=N, M>N.");

        for (i = 1; i <= 3; i++)
        {
            switch (i)
            {
                case 1:
                    m = 3;
                    n = 5;
                    break;
                case 2:
                    m = 5;
                    n = 5;
                    break;
                case 3:
                    m = 5;
                    n = 3;
                    break;
            }

            int ahi = Math.Min(m - 1, n);
            int bhi = Math.Min(m, n);
            int chi = Math.Min(m, n - 1);
            double[] a1 = new double[ahi];
            double[] a2 = new double[bhi];
            double[] a3 = new double[chi];

            r83v_indicator(m, n, ref a1, ref a2, ref a3);
            r83v_print(m, n, a1, a2, a3, "  R83V matrix A:");

            double[] a = r83v_to_r8vec(m, n, a1, a2, a3);
            r8vec_print(ahi + bhi + chi, a, "  Vector version of A:");

        }

    }

    public static void r83v_transpose(int m, int n, double[] a1, double[] a2, double[] a3,
            ref double[] b1, ref double[] b2, ref double[] b3)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83V_TRANSPOSE makes a transposed copy of an R83V matrix.
        //
        //  Discussion:
        //
        //    The R83V storage format is used for a tridiagonal matrix.
        //    The subdiagonal is in A(min(M-1,N)).
        //    The diagonal is in B(min(M,N)).
        //    The superdiagonal is in C(min(M,N-1)).
        //
        //  Example:
        //
        //    An R83V matrix of order 3x5 would be stored:
        //
        //      B1  C1  **  **  **
        //      A1  B2  C2  **  **
        //      **  A2  B3  C3  **
        //
        //    An R83 matrix of order 5x5 would be stored:
        //
        //      B1  C1  **  **  **
        //      A1  B2  C2  **  **
        //      **  A2  B3  C3  **
        //      **  **  A3  B4  C4
        //      **  **  **  A4  B5
        //
        //    An R83 matrix of order 5x3 would be stored:
        //
        //      B1  C1  **
        //      A1  B2  C2
        //      **  A2  B3
        //      **  **  A3
        //      **  **  **
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 February 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of the matrix.
        //
        //    Input, double A1(min(M-1,N)), A2(min(M,N)), A3(min(M,N-1)), the matrix.
        //
        //    Output, double B1(min(N-1,M)), B2(min(N,M)), B3(min(N,M-1)), the copy.
        //
    {
        int i;

        int ahi = Math.Min(m - 1, n);
        int bhi = Math.Min(m, n);
        int chi = Math.Min(m, n - 1);

        for (i = 0; i < ahi; i++)
        {
            b3[i] = a1[i];
        }

        for (i = 0; i < bhi; i++)
        {
            b2[i] = a2[i];
        }

        for (i = 0; i < chi; i++)
        {
            b1[i] = a3[i];
        }
    }

    public static void r83v_transpose_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83V_TRANSPOSE_TEST tests R83V_TRANSPOSE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 February 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int m = 0;
        int n = 0;

        Console.WriteLine("");
        Console.WriteLine("R83V_TRANSPOSE_TEST");
        Console.WriteLine("  R83V_TRANSPOSE makes a transposed copy of an R83V matrix.");
        Console.WriteLine("  We check three cases, M<N, M=N, M>N.");

        for (i = 1; i <= 3; i++)
        {
            switch (i)
            {
                case 1:
                    m = 3;
                    n = 5;
                    break;
                case 2:
                    m = 5;
                    n = 5;
                    break;
                case 3:
                    m = 5;
                    n = 3;
                    break;
            }

            int ahi = Math.Min(m - 1, n);
            int bhi = Math.Min(m, n);
            int chi = Math.Min(m, n - 1);
            double[] a1 = new double[ahi];
            double[] a2 = new double[bhi];
            double[] a3 = new double[chi];

            r83v_indicator(m, n, ref a1, ref a2, ref a3);
            r83v_print(m, n, a1, a2, a3, "  R83V matrix A:");

            double[] b1 = new double[chi];
            double[] b2 = new double[bhi];
            double[] b3 = new double[ahi];
            r83v_transpose(m, n, a1, a2, a3, ref b1, ref b2, ref b3);
            r83v_print(n, m, b1, b2, b3, "  B = copy of A:");

        }
    }

    public static void r83v_zeros(int m, int n, ref double[] a, ref double[] b, ref double[] c)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83V_ZEROS returns the zero matrix in R83V format.
        //
        //  Discussion:
        //
        //    The R83V storage format is used for a tridiagonal matrix.
        //    The subdiagonal is in A(min(M-1,N)).
        //    The diagonal is in B(min(M,N)).
        //    The superdiagonal is in C(min(M,N-1)).
        //
        //  Example:
        //
        //    An R83V matrix of order 3x5 would be stored:
        //
        //      B1  C1  **  **  **
        //      A1  B2  C2  **  **
        //      **  A2  B3  C3  **
        //
        //    An R83 matrix of order 5x5 would be stored:
        //
        //      B1  C1  **  **  **
        //      A1  B2  C2  **  **
        //      **  A2  B3  C3  **
        //      **  **  A3  B4  C4
        //      **  **  **  A4  B5
        //
        //    An R83 matrix of order 5x3 would be stored:
        //
        //      B1  C1  **
        //      A1  B2  C2
        //      **  A2  B3
        //      **  **  A3
        //      **  **  **
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 February 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the order of the matrix.
        //
        //    Output, double A(min(M-1,N)), B(min(M,N)), C(min(M,N-1)), the matrix.
        //
    {
        int i;

        int ahi = Math.Min(m - 1, n);
        int bhi = Math.Min(m, n);
        int chi = Math.Min(m, n - 1);

        for (i = 0; i < ahi; i++)
        {
            a[i] = 0.0;
        }

        for (i = 0; i < bhi; i++)
        {
            b[i] = 0.0;
        }

        for (i = 0; i < chi; i++)
        {
            c[i] = 0.0;
        }
    }

    public static void r83v_zeros_test()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R83V_ZEROS_TEST tests R83V_ZEROS.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    17 February 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        int m = 0;
        int n = 0;

        Console.WriteLine("");
        Console.WriteLine("R83V_ZEROS_TEST");
        Console.WriteLine("  R83V_ZEROS sets up an R83V zero matrix.");
        Console.WriteLine("  We check three cases, M<N, M=N, M>N.");

        for (i = 1; i <= 3; i++)
        {
            switch (i)
            {
                case 1:
                    m = 3;
                    n = 5;
                    break;
                case 2:
                    m = 5;
                    n = 5;
                    break;
                case 3:
                    m = 5;
                    n = 3;
                    break;
            }

            int ahi = Math.Min(m - 1, n);
            int bhi = Math.Min(m, n);
            int chi = Math.Min(m, n - 1);
            double[] a = new double[ahi];
            double[] b = new double[bhi];
            double[] c = new double[chi];

            r83v_zeros(m, n, ref a, ref b, ref c);

            r83v_print(m, n, a, b, c, "  The R83V zero matrix:");

        }

    }

}