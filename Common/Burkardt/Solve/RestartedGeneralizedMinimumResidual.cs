using System;
using Burkardt.CompressedRow;
using Burkardt.Sparse;
using Burkardt.SparseTripletNS;
using Burkardt.Types;

namespace Burkardt.SolveNS;

public static class RestartedGeneralizedMinimumResidual
{
    public static void mgmres(double[] a, int[] ia, int[] ja, ref double[] x, double[] rhs,
            int n, int nz_num, int itr_max, int mr, double tol_abs, double tol_rel)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MGMRES applies the restarted GMRES iteration to a linear system.
        //
        //  Discussion:
        //
        //    The linear system A*X=B is solved iteratively.
        //
        //    The matrix A is assumed to be sparse.  To save on storage, only
        //    the nonzero entries of A are stored.  For instance, the K-th nonzero
        //    entry in the matrix is stored by:
        //
        //      A(K) = value of entry,
        //      IA(K) = row of entry,
        //      JA(K) = column of entry.
        //
        //    The "matrices" H and V are treated as one-dimensional vectors
        //    which store the matrix data in row major form.
        //
        //    This requires that references to H[I][J] be replaced by references
        //    to H[I+J*(MR+1)] and references to V[I][J] by V[I+J*N].
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 July 2007
        //
        //  Author:
        //
        //    Original C version by Lili Ju.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
        //    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
        //    Charles Romine, Henk van der Vorst,
        //    Templates for the Solution of Linear Systems:
        //    Building Blocks for Iterative Methods,
        //    SIAM, 1994,
        //    ISBN: 0898714710,
        //    LC: QA297.8.T45.
        //
        //    Tim Kelley,
        //    Iterative Methods for Linear and Nonlinear Equations,
        //    SIAM, 2004,
        //    ISBN: 0898713528,
        //    LC: QA297.8.K45.
        //
        //    Yousef Saad,
        //    Iterative Methods for Sparse Linear Systems,
        //    Second Edition,
        //    SIAM, 2003,
        //    ISBN: 0898715342,
        //    LC: QA188.S17.
        //
        //  Parameters:
        //
        //    Input, double A[NZ_NUM], the matrix values.
        //
        //    Input, int IA[NZ_NUM], JA[NZ_NUM], the row and column indices
        //    of the matrix values.
        //
        //    Input/output, double X[N]; on input, an approximation to
        //    the solution.  On output, an improved approximation.
        //
        //    Input, double RHS[N], the right hand side of the linear system.
        //
        //    Input, int N, the order of the linear system.
        //
        //    Input, int NZ_NUM, the number of nonzero matrix values.
        //
        //    Input, int ITR_MAX, the maximum number of (outer) iterations to take.
        //
        //    Input, int MR, the maximum number of (inner) iterations to take.
        //    MR must be less than N.
        //
        //    Input, double TOL_ABS, an absolue tolerance applied to the
        //    current residual.
        //
        //    Input, double TOL_REL, a relative tolerance comparing the
        //    current residual to the initial residual.
        //
    {
        double av;
        double[] c;
        double delta = 1.0e-03;
        double[] g;
        double[] h;
        double htmp;
        int i;
        int itr;
        int itr_used;
        int j;
        int k;
        int k_copy = 0;
        double mu;
        double[] r;
        double rho = 0;
        double rho_tol = 0;
        double[] s;
        double[] v;
        bool verbose = true;
        double[] y;

        c = new double[mr];
        g = new double[mr + 1];
        h = new double[(mr + 1) * mr];
        r = new double[n];
        s = new double[mr];
        v = new double[n * (mr + 1)];
        y = new double[mr + 1];

        itr_used = 0;

        if (n < mr)
        {
            Console.WriteLine("");
            Console.WriteLine("MGMRES - Fatal error!");
            Console.WriteLine("  N < MR.");
            Console.WriteLine("  N = " + n + "");
            Console.WriteLine("  MR = " + mr + "");
            return;
        }

        for (itr = 1; itr <= itr_max; itr++)
        {
            AX.ax(a, ia, ja, x, ref r, n, nz_num);

            for (i = 0; i < n; i++)
            {
                r[i] = rhs[i] - r[i];
            }

            rho = Math.Sqrt(typeMethods.r8vec_dot_product(n, r, r));

            switch (verbose)
            {
                case true:
                    Console.WriteLine("  ITR = " + itr + "  Residual = " + rho + "");
                    break;
            }

            rho_tol = itr switch
            {
                1 => rho * tol_rel,
                _ => rho_tol
            };

            for (i = 0; i < n; i++)
            {
                v[i + 0 * n] = r[i] / rho;
            }

            g[0] = rho;
            for (i = 1; i <= mr; i++)
            {
                g[i] = 0.0;
            }

            for (i = 0; i < mr + 1; i++)
            {
                for (j = 0; j < mr; j++)
                {
                    h[i + j * (mr + 1)] = 0.0;
                }
            }

            for (k = 1; k <= mr; k++)
            {
                k_copy = k;

                AX.ax(a, ia, ja, v, ref v, n, nz_num, xIndex:  + (k - 1) * n, wIndex:  + k * n);

                av = Math.Sqrt(typeMethods.r8vec_dot_product(n, v, v, a1Index: +k * n, a2Index: +k * n));

                for (j = 1; j <= k; j++)
                {
                    h[j - 1 + (k - 1) * (mr + 1)] =
                        typeMethods.r8vec_dot_product(n, v, v, a1Index: +k * n, a2Index: +(j - 1) * n);
                    for (i = 0; i < n; i++)
                    {
                        v[i + k * n] -= h[j - 1 + (k - 1) * (mr + 1)] * v[i + (j - 1) * n];
                    }
                }

                h[k + (k - 1) * (mr + 1)] =
                    Math.Sqrt(typeMethods.r8vec_dot_product(n, v, v, a1Index: +k * n, a2Index: +k * n));

                if (av + delta * h[k + (k - 1) * (mr + 1)] == av)
                {
                    for (j = 1; j <= k; j++)
                    {
                        htmp = typeMethods.r8vec_dot_product(n, v, v, a1Index: +k * n, a2Index: +(j - 1) * n);
                        h[j - 1 + (k - 1) * (mr + 1)] += htmp;
                        for (i = 0; i < n; i++)
                        {
                            v[i + k * n] -= htmp * v[i + (j - 1) * n];
                        }
                    }

                    h[k + (k - 1) * (mr + 1)] =
                        Math.Sqrt(typeMethods.r8vec_dot_product(n, v, v, a1Index: +k * n, a2Index: +k * n));
                }

                if (h[k + (k - 1) * (mr + 1)] != 0.0)
                {
                    for (i = 0; i < n; i++)
                    {
                        v[i + k * n] /= h[k + (k - 1) * (mr + 1)];
                    }
                }

                switch (k)
                {
                    case > 1:
                    {
                        for (i = 1; i <= k + 1; i++)
                        {
                            y[i - 1] = h[i - 1 + (k - 1) * (mr + 1)];
                        }

                        for (j = 1; j <= k - 1; j++)
                        {
                            Helpers.mult_givens(c[j - 1], s[j - 1], j - 1, ref y);
                        }

                        for (i = 1; i <= k + 1; i++)
                        {
                            h[i - 1 + (k - 1) * (mr + 1)] = y[i - 1];
                        }

                        break;
                    }
                }

                mu = Math.Sqrt(Math.Pow(h[k - 1 + (k - 1) * (mr + 1)], 2)
                               + Math.Pow(h[k + (k - 1) * (mr + 1)], 2));
                c[k - 1] = h[k - 1 + (k - 1) * (mr + 1)] / mu;
                s[k - 1] = -h[k + (k - 1) * (mr + 1)] / mu;
                h[k - 1 + (k - 1) * (mr + 1)] = c[k - 1] * h[k - 1 + (k - 1) * (mr + 1)]
                                                - s[k - 1] * h[k + (k - 1) * (mr + 1)];
                h[k + (k - 1) * (mr + 1)] = 0;
                Helpers.mult_givens(c[k - 1], s[k - 1], k - 1, ref g);

                rho = Math.Abs(g[k]);

                itr_used += 1;

                switch (verbose)
                {
                    case true:
                        Console.WriteLine("  K =   " + k + "  Residual = " + rho + "");
                        break;
                }

                if (rho <= rho_tol && rho <= tol_abs)
                {
                    break;
                }
            }

            k = k_copy - 1;
            y[k] = g[k] / h[k + k * (mr + 1)];

            for (i = k; 1 <= i; i--)
            {
                y[i - 1] = g[i - 1];
                for (j = i + 1; j <= k + 1; j++)
                {
                    y[i - 1] -= h[i - 1 + (j - 1) * (mr + 1)] * y[j - 1];
                }

                y[i - 1] /= h[i - 1 + (i - 1) * (mr + 1)];
            }

            for (i = 1; i <= n; i++)
            {
                for (j = 1; j <= k + 1; j++)
                {
                    x[i - 1] += v[i - 1 + (j - 1) * n] * y[j - 1];
                }
            }

            if (rho <= rho_tol && rho <= tol_abs)
            {
                break;
            }
        }

        switch (verbose)
        {
            case true:
                Console.WriteLine("");
                Console.WriteLine("MGMRES");
                Console.WriteLine("  Number of iterations = " + itr_used + "");
                Console.WriteLine("  Final residual = " + rho + "");
                break;
        }
    }

    public static void mgmres_st(int n, int nz_num, int[] ia, int[] ja, double[] a, ref double[] x,
            double[] rhs, int itr_max, int mr, double tol_abs, double tol_rel)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MGMRES_ST applies restarted GMRES to a matrix in sparse triplet form.
        //
        //  Discussion:
        //
        //    The linear system A*X=B is solved iteratively.
        //
        //    The matrix A is assumed to be stored in sparse triplet form.  Only
        //    the nonzero entries of A are stored.  For instance, the K-th nonzero
        //    entry in the matrix is stored by:
        //
        //      A(K) = value of entry,
        //      IA(K) = row of entry,
        //      JA(K) = column of entry.
        //
        //    The "matrices" H and V are treated as one-dimensional vectors
        //    which store the matrix data in row major form.
        //
        //    This requires that references to H[I][J] be replaced by references
        //    to H[I+J*(MR+1)] and references to V[I][J] by V[I+J*N].
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    18 July 2007
        //
        //  Author:
        //
        //    Original C version by Lili Ju.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
        //    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
        //    Charles Romine, Henk van der Vorst,
        //    Templates for the Solution of Linear Systems:
        //    Building Blocks for Iterative Methods,
        //    SIAM, 1994,
        //    ISBN: 0898714710,
        //    LC: QA297.8.T45.
        //
        //    Tim Kelley,
        //    Iterative Methods for Linear and Nonlinear Equations,
        //    SIAM, 2004,
        //    ISBN: 0898713528,
        //    LC: QA297.8.K45.
        //
        //    Yousef Saad,
        //    Iterative Methods for Sparse Linear Systems,
        //    Second Edition,
        //    SIAM, 2003,
        //    ISBN: 0898715342,
        //    LC: QA188.S17.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the linear system.
        //
        //    Input, int NZ_NUM, the number of nonzero matrix values.
        //
        //    Input, int IA[NZ_NUM], JA[NZ_NUM], the row and column indices
        //    of the matrix values.
        //
        //    Input, double A[NZ_NUM], the matrix values.
        //
        //    Input/output, double X[N]; on input, an approximation to
        //    the solution.  On output, an improved approximation.
        //
        //    Input, double RHS[N], the right hand side of the linear system.
        //
        //    Input, int ITR_MAX, the maximum number of (outer) iterations to take.
        //
        //    Input, int MR, the maximum number of (inner) iterations to take.
        //    MR must be less than N.
        //
        //    Input, double TOL_ABS, an absolute tolerance applied to the
        //    current residual.
        //
        //    Input, double TOL_REL, a relative tolerance comparing the
        //    current residual to the initial residual.
        //
    {
        double av;
        double[] c;
        double delta = 1.0e-03;
        double[] g;
        double[] h;
        double htmp;
        int i;
        int itr;
        int itr_used;
        int j;
        int k;
        int k_copy = 0;
        double mu;
        double[] r;
        double rho = 0;
        double rho_tol = 0;
        double[] s;
        double[] v;
        bool verbose = true;
        double[] y;

        c = new double[mr];
        g = new double[mr + 1];
        h = new double[(mr + 1) * mr];
        r = new double[n];
        s = new double[mr];
        v = new double[n * (mr + 1)];
        y = new double[mr + 1];

        itr_used = 0;

        if (n < mr)
        {
            Console.WriteLine("");
            Console.WriteLine("MGMRES_ST - Fatal error!");
            Console.WriteLine("  N < MR.");
            Console.WriteLine("  N = " + n + "");
            Console.WriteLine("  MR = " + mr + "");
            return;
        }

        for (itr = 1; itr <= itr_max; itr++)
        {
            AXST.ax_st(n, nz_num, ia, ja, a, x, ref r);

            for (i = 0; i < n; i++)
            {
                r[i] = rhs[i] - r[i];
            }

            rho = Math.Sqrt(typeMethods.r8vec_dot(n, r, r));

            switch (verbose)
            {
                case true:
                    Console.WriteLine("  ITR = " + itr + "  Residual = " + rho + "");
                    break;
            }

            rho_tol = itr switch
            {
                1 => rho * tol_rel,
                _ => rho_tol
            };

            for (i = 0; i < n; i++)
            {
                v[i + 0 * n] = r[i] / rho;
            }

            g[0] = rho;
            for (i = 1; i <= mr; i++)
            {
                g[i] = 0.0;
            }

            for (i = 0; i < mr + 1; i++)
            {
                for (j = 0; j < mr; j++)
                {
                    h[i + j * (mr + 1)] = 0.0;
                }
            }

            for (k = 1; k <= mr; k++)
            {
                k_copy = k;

                AXST.ax_st(n, nz_num, ia, ja, a, v, ref v, xIndex: + (k - 1) * n, wIndex: + k * n);

                av = Math.Sqrt(typeMethods.r8vec_dot(n, v, v, a1Index: + k * n, a2Index: + k * n));

                for (j = 1; j <= k; j++)
                {
                    h[j - 1 + (k - 1) * (mr + 1)] = typeMethods.r8vec_dot(n, v, v, a1Index: + k * n, a2Index: + (j - 1) * n);
                    for (i = 0; i < n; i++)
                    {
                        v[i + k * n] -= h[j - 1 + (k - 1) * (mr + 1)] * v[i + (j - 1) * n];
                    }
                }

                h[k + (k - 1) * (mr + 1)] = Math.Sqrt(typeMethods.r8vec_dot(n, v, v, a1Index: + k * n, a2Index: + k * n));

                if (av + delta * h[k + (k - 1) * (mr + 1)] == av)
                {
                    for (j = 1; j <= k; j++)
                    {
                        htmp = typeMethods.r8vec_dot(n, v, v, a1Index: + k * n, a2Index: + (j - 1) * n);
                        h[j - 1 + (k - 1) * (mr + 1)] += htmp;
                        for (i = 0; i < n; i++)
                        {
                            v[i + k * n] -= htmp * v[i + (j - 1) * n];
                        }
                    }

                    h[k + (k - 1) * (mr + 1)] = Math.Sqrt(typeMethods.r8vec_dot(n, v, v, a1Index: + k * n, a2Index: + k * n));
                }

                if (h[k + (k - 1) * (mr + 1)] != 0.0)
                {
                    for (i = 0; i < n; i++)
                    {
                        v[i + k * n] /= h[k + (k - 1) * (mr + 1)];
                    }
                }

                switch (k)
                {
                    case > 1:
                    {
                        for (i = 1; i <= k + 1; i++)
                        {
                            y[i - 1] = h[i - 1 + (k - 1) * (mr + 1)];
                        }

                        for (j = 1; j <= k - 1; j++)
                        {
                            Helpers.mult_givens(c[j - 1], s[j - 1], j - 1, ref y);
                        }

                        for (i = 1; i <= k + 1; i++)
                        {
                            h[i - 1 + (k - 1) * (mr + 1)] = y[i - 1];
                        }

                        break;
                    }
                }

                mu = Math.Sqrt(Math.Pow(h[k - 1 + (k - 1) * (mr + 1)], 2)
                               + Math.Pow(h[k + (k - 1) * (mr + 1)], 2));
                c[k - 1] = h[k - 1 + (k - 1) * (mr + 1)] / mu;
                s[k - 1] = -h[k + (k - 1) * (mr + 1)] / mu;
                h[k - 1 + (k - 1) * (mr + 1)] = c[k - 1] * h[k - 1 + (k - 1) * (mr + 1)]
                                                - s[k - 1] * h[k + (k - 1) * (mr + 1)];
                h[k + (k - 1) * (mr + 1)] = 0;
                Helpers.mult_givens(c[k - 1], s[k - 1], k - 1, ref g);

                rho = Math.Abs(g[k]);

                itr_used += 1;

                switch (verbose)
                {
                    case true:
                        Console.WriteLine("  K =   " + k + "  Residual = " + rho + "");
                        break;
                }

                if (rho <= rho_tol && rho <= tol_abs)
                {
                    break;
                }
            }

            k = k_copy - 1;
            y[k] = g[k] / h[k + k * (mr + 1)];

            for (i = k; 1 <= i; i--)
            {
                y[i - 1] = g[i - 1];
                for (j = i + 1; j <= k + 1; j++)
                {
                    y[i - 1] -= h[i - 1 + (j - 1) * (mr + 1)] * y[j - 1];
                }

                y[i - 1] /= h[i - 1 + (i - 1) * (mr + 1)];
            }

            for (i = 1; i <= n; i++)
            {
                for (j = 1; j <= k + 1; j++)
                {
                    x[i - 1] += v[i - 1 + (j - 1) * n] * y[j - 1];
                }
            }

            if (rho <= rho_tol && rho <= tol_abs)
            {
                break;
            }
        }

        switch (verbose)
        {
            case true:
                Console.WriteLine("");
                Console.WriteLine("MGMRES_ST");
                Console.WriteLine("  Number of iterations = " + itr_used + "");
                Console.WriteLine("  Final residual = " + rho + "");
                break;
        }
    }

    public static void pmgmres_ilu_cr(int n, int nz_num, int[] ia, int[] ja, double[] a,
            ref double[] x, double[] rhs, int itr_max, int mr, double tol_abs,
            double tol_rel)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PMGMRES_ILU_CR applies the preconditioned restarted GMRES algorithm.
        //
        //  Discussion:
        //
        //    The matrix A is assumed to be stored in compressed row format.  Only
        //    the nonzero entries of A are stored.  The vector JA stores the
        //    column index of the nonzero value.  The nonzero values are sorted
        //    by row, and the compressed row vector IA then has the property that
        //    the entries in A and JA that correspond to row I occur in indices
        //    IA[I] through IA[I+1]-1.
        //
        //    This routine uses the incomplete LU decomposition for the
        //    preconditioning.  This preconditioner requires that the sparse
        //    matrix data structure supplies a storage position for each diagonal
        //    element of the matrix A, and that each diagonal element of the
        //    matrix A is not zero.
        //
        //    Thanks to Jesus Pueblas Sanchez-Guerra for supplying two
        //    corrections to the code on 31 May 2007.
        //
        //    This implementation of the code stores the doubly-dimensioned arrays
        //    H and V as vectors.  However, it follows the C convention of storing
        //    them by rows, rather than my own preference for storing them by
        //    columns.   I may come back and change this some time.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    26 July 2007
        //
        //  Author:
        //
        //    Original C version by Lili Ju.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
        //    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
        //    Charles Romine, Henk van der Vorst,
        //    Templates for the Solution of Linear Systems:
        //    Building Blocks for Iterative Methods,
        //    SIAM, 1994.
        //    ISBN: 0898714710,
        //    LC: QA297.8.T45.
        //
        //    Tim Kelley,
        //    Iterative Methods for Linear and Nonlinear Equations,
        //    SIAM, 2004,
        //    ISBN: 0898713528,
        //    LC: QA297.8.K45.
        //
        //    Yousef Saad,
        //    Iterative Methods for Sparse Linear Systems,
        //    Second Edition,
        //    SIAM, 2003,
        //    ISBN: 0898715342,
        //    LC: QA188.S17.
        //
        //  Parameters:
        //
        //    Input, int N, the order of the linear system.
        //
        //    Input, int NZ_NUM, the number of nonzero matrix values.
        //
        //    Input, int IA[N+1], JA[NZ_NUM], the row and column indices
        //    of the matrix values.  The row vector has been compressed.
        //
        //    Input, double A[NZ_NUM], the matrix values.
        //
        //    Input/output, double X[N]; on input, an approximation to
        //    the solution.  On output, an improved approximation.
        //
        //    Input, double RHS[N], the right hand side of the linear system.
        //
        //    Input, int ITR_MAX, the maximum number of (outer) iterations to take.
        //
        //    Input, int MR, the maximum number of (inner) iterations to take.
        //    MR must be less than N.
        //
        //    Input, double TOL_ABS, an absolute tolerance applied to the
        //    current residual.
        //
        //    Input, double TOL_REL, a relative tolerance comparing the
        //    current residual to the initial residual.
        //
    {
        double av;
        double[] c;
        double delta = 1.0e-03;
        double[] g;
        double[] h;
        double htmp;
        int i;
        int itr;
        int itr_used;
        int j;
        int k;
        int k_copy = 0;
        double[] l;
        double mu;
        double[] r;
        double rho = 0;
        double rho_tol = 0;
        double[] s;
        int[] ua;
        double[] v;
        bool verbose = true;
        double[] y;

        itr_used = 0;

        c = new double[mr + 1];
        g = new double[mr + 1];
        h = new double[(mr + 1) * mr];
        l = new double[ia[n] + 1];
        r = new double[n];
        s = new double[mr + 1];
        ua = new int[n];
        v = new double[(mr + 1) * n];
        y = new double[mr + 1];

        RearrangeCR.rearrange_cr(n, nz_num, ia, ref ja, ref a);

        DiagonalPointerCR.diagonal_pointer_cr(n, nz_num, ia, ja, ref ua);

        ILUCR.ilu_cr(n, nz_num, ia, ja, a, ua, ref l);

        switch (verbose)
        {
            case true:
                Console.WriteLine("");
                Console.WriteLine("PMGMRES_ILU_CR");
                Console.WriteLine("  Number of unknowns = " + n + "");
                break;
        }

        for (itr = 0; itr < itr_max; itr++)
        {
            AXCR.ax_cr(n, nz_num, ia, ja, a, x, ref r);

            for (i = 0; i < n; i++)
            {
                r[i] = rhs[i] - r[i];
            }

            LUSCR.lus_cr(n, nz_num, ia, ja, l, ua, r, ref r);

            rho = Math.Sqrt(typeMethods.r8vec_dot(n, r, r));

            switch (verbose)
            {
                case true:
                    Console.WriteLine("  ITR = " + itr + "  Residual = " + rho + "");
                    break;
            }

            rho_tol = itr switch
            {
                0 => rho * tol_rel,
                _ => rho_tol
            };

            for (i = 0; i < n; i++)
            {
                v[0 * n + i] = r[i] / rho;
            }

            g[0] = rho;
            for (i = 1; i < mr + 1; i++)
            {
                g[i] = 0.0;
            }

            for (i = 0; i < mr + 1; i++)
            {
                for (j = 0; j < mr; j++)
                {
                    h[i * mr + j] = 0.0;
                }
            }

            for (k = 0; k < mr; k++)
            {
                k_copy = k;

                AXCR.ax_cr(n, nz_num, ia, ja, a, v, ref v, xIndex: + k * n, wIndex: + (k + 1) * n);

                LUSCR.lus_cr(n, nz_num, ia, ja, l, ua, v, ref v, rIndex: + (k + 1) * n, zIndex: + (k + 1) * n);

                av = Math.Sqrt(typeMethods.r8vec_dot(n, v, v, a1Index: + (k + 1) * n, a2Index: + (k + 1) * n));

                for (j = 0; j <= k; j++)
                {
                    h[j * mr + k] = typeMethods.r8vec_dot(n, v, v, a1Index: + (k + 1) * n, a2Index: + j * n);
                    for (i = 0; i < n; i++)
                    {
                        v[(k + 1) * n + i] -= h[j * mr + k] * v[j * n + i];
                    }
                }

                h[(k + 1) * mr + k] = Math.Sqrt(typeMethods.r8vec_dot(n, v, v, a1Index: + (k + 1) * n, a2Index: + (k + 1) * n));

                if (av + delta * h[(k + 1) * mr + k] == av)
                {
                    for (j = 0; j < k + 1; j++)
                    {
                        htmp = typeMethods.r8vec_dot(n, v, v, a1Index: + (k + 1) * n, a2Index: + j * n);
                        h[j * mr + k] += htmp;
                        for (i = 0; i < n; i++)
                        {
                            v[(k + 1) * n + i] -= htmp * v[j * n + i];
                        }
                    }

                    h[(k + 1) * mr + k] = Math.Sqrt(typeMethods.r8vec_dot(n, v, v, a1Index: + (k + 1) * n, a2Index: + (k + 1) * n));
                }

                if (h[(k + 1) * mr + k] != 0.0)
                {
                    for (i = 0; i < n; i++)
                    {
                        v[(k + 1) * n + i] /= h[(k + 1) * mr + k];
                    }
                }

                switch (k)
                {
                    case > 0:
                    {
                        for (i = 0; i < k + 2; i++)
                        {
                            y[i] = h[i * mr + k];
                        }

                        for (j = 0; j < k; j++)
                        {
                            Helpers.mult_givens(c[j], s[j], j, ref y);
                        }

                        for (i = 0; i < k + 2; i++)
                        {
                            h[i * mr + k] = y[i];
                        }

                        break;
                    }
                }

                mu = Math.Sqrt(h[k * mr + k] * h[k * mr + k] + h[(k + 1) * mr + k] * h[(k + 1) * mr + k]);
                c[k] = h[k * mr + k] / mu;
                s[k] = -h[(k + 1) * mr + k] / mu;
                h[k * mr + k] = c[k] * h[k * mr + k] - s[k] * h[(k + 1) * mr + k];
                h[(k + 1) * mr + k] = 0.0;
                Helpers.mult_givens(c[k], s[k], k, ref g);

                rho = Math.Abs(g[k + 1]);

                itr_used += 1;

                switch (verbose)
                {
                    case true:
                        Console.WriteLine("  K   = " + k + "  Residual = " + rho + "");
                        break;
                }

                if (rho <= rho_tol && rho <= tol_abs)
                {
                    break;
                }
            }

            k = k_copy;

            y[k] = g[k] / h[k * mr + k];
            for (i = k - 1; 0 <= i; i--)
            {
                y[i] = g[i];
                for (j = i + 1; j < k + 1; j++)
                {
                    y[i] -= h[i * mr + j] * y[j];
                }

                y[i] /= h[i * mr + i];
            }

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < k + 1; j++)
                {
                    x[i] += v[j * n + i] * y[j];
                }
            }

            if (rho <= rho_tol && rho <= tol_abs)
            {
                break;
            }
        }

        switch (verbose)
        {
            case true:
                Console.WriteLine("");
                ;
                Console.WriteLine("PMGMRES_ILU_CR:");
                Console.WriteLine("  Iterations = " + itr_used + "");
                Console.WriteLine("  Final residual = " + rho + "");
                break;
        }
    }
}