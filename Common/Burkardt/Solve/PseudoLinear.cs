using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.SolveNS;

public static class PseudoLinear
{

    public static double[] pseudo_inverse(int m, int n, double[] u, double[] s,
            double[] v)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PSEUDO_INVERSE computes the pseudoinverse.
        //
        //  Discussion:
        //
        //    Given the singular value decomposition of a real MxN matrix A:
        //
        //      A = U * S * V'
        //
        //    where 
        //
        //      U is MxM orthogonal,
        //      S is MxN, and entirely zero except for the diagonal;
        //      V is NxN orthogonal.
        //
        //    the pseudo inverse is the NxM matrix A+ with the form
        //
        //      A+ = V * S+ * U'
        //
        //    where 
        //
        //      S+ is the NxM matrix whose nonzero diagonal elements are
        //      the inverses of the corresponding diagonal elements of S.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns in the matrix A.
        //
        //    Input, double A[M*N], the matrix whose singular value
        //    decomposition we are investigating.
        //
        //    Input, double U[M*M], S[M*N], V[N*N], the factors
        //    that form the singular value decomposition of A.
        //
        //    Output, double PSEUDO_INVERSE[N*M], the pseudo_inverse of A.
        //
    {
        double[] a_pseudo;
        int i;
        int j;
        int k;
        double[] sp;
        double[] sput;

        sp = new double[n * m];
        for (j = 0; j < m; j++)
        {
            for (i = 0; i < n; i++)
            {
                if (i == j && s[i + i * m] != 0.0)
                {
                    sp[i + j * n] = 1.0 / s[i + i * m];
                }
                else
                {
                    sp[i + j * n] = 0.0;
                }
            }
        }

        sput = new double[n * m];
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < m; j++)
            {
                sput[i + j * n] = 0.0;
                for (k = 0; k < m; k++)
                {
                    sput[i + j * n] += sp[i + k * n] * u[j + k * m];
                }
            }
        }

        a_pseudo = new double[n * m];
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < m; j++)
            {
                a_pseudo[i + j * n] = 0.0;
                for (k = 0; k < n; k++)
                {
                    a_pseudo[i + j * n] += v[i + k * n] * sput[k + j * n];
                }
            }
        }


        return a_pseudo;
    }

    public static void pseudo_linear_solve_test(int m, int n, double[] a,
            double[] a_pseudo, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PSEUDO_LINEAR_SOLVE_TEST uses the pseudoinverse for linear systems.
        //
        //  Discussion:
        //
        //    Given an MxN matrix A, and its pseudoinverse A+:
        //
        //      "Solve" A  * x = b by x = A+  * b.
        //
        //      "Solve" A' * x = b by x = A+' * b.
        //
        //    When the system is overdetermined, the solution minimizes the
        //    L2 norm of the residual.  
        //
        //    When the system is underdetermined, the solution
        //    is the solution of minimum L2 norm.     
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns in the matrix A.
        //
        //    Input, double A[M*N], the matrix whose singular value
        //    decomposition we are investigating.
        //
        //    Input, double A_PSEUDO[N*M], the pseudo_inverse of A.
        //
        //    Input/output, int *SEED, a seed for the random number generator.
        //
    {
        double[] bm;
        double[] bn;
        int i;
        int j;
        double[] rm;
        double[] rn;
        double[] xm1;
        double[] xm2;
        double[] xn1;
        double[] xn2;

        Console.WriteLine("");
        Console.WriteLine("PSEUDO_LINEAR_SOLVE_TEST");
        //
        //  A * x = b, b in range of A.
        //
        xn1 = UniformRNG.r8vec_uniform_01_new(n, ref seed);
        for (i = 0; i < n; i++)
        {
            xn1[i] = typeMethods.r8_nint(10.0 * xn1[i]);
        }

        bm = new double[m];
        for (i = 0; i < m; i++)
        {
            bm[i] = 0.0;
            for (j = 0; j < n; j++)
            {
                bm[i] += a[i + j * m] * xn1[j];
            }
        }

        xn2 = new double[n];
        for (i = 0; i < n; i++)
        {
            xn2[i] = 0.0;
            for (j = 0; j < m; j++)
            {
                xn2[i] += a_pseudo[i + j * n] * bm[j];
            }
        }

        rm = new double[m];
        for (i = 0; i < m; i++)
        {
            rm[i] = bm[i];
            for (j = 0; j < n; j++)
            {
                rm[i] -= a[i + j * m] * xn2[j];
            }
        }

        Console.WriteLine("");
        Console.WriteLine("  Given:");
        Console.WriteLine("    b = A * x1");
        Console.WriteLine("  so that b is in the range of A, solve");
        Console.WriteLine("    A * x = b");
        Console.WriteLine("  using the pseudoinverse:");
        Console.WriteLine("    x2 = A+ * b.");
        Console.WriteLine("");
        Console.WriteLine("  Norm of x1 = " + typeMethods.r8vec_norm_l2(n, xn1) + "");
        Console.WriteLine("  Norm of x2 = " + typeMethods.r8vec_norm_l2(n, xn2) + "");
        Console.WriteLine("  Norm of residual = " + typeMethods.r8vec_norm_l2(m, rm) + "");

        //
        //  A * x = b, b not in range of A.
        //
        if (n < m)
        {
            Console.WriteLine("");
            Console.WriteLine("  For N < M, most systems A*x=b will not be");
            Console.WriteLine("  exactly and uniquely solvable, except in the");
            Console.WriteLine("  least squares sense.");
            Console.WriteLine("");
            Console.WriteLine("  Here is an example:");

            bm = UniformRNG.r8vec_uniform_01_new(m, ref seed);

            xn2 = new double[n];
            for (i = 0; i < n; i++)
            {
                xn2[i] = 0.0;
                for (j = 0; j < m; j++)
                {
                    xn2[i] += a_pseudo[i + j * n] * bm[j];
                }
            }

            rm = new double[m];
            for (i = 0; i < m; i++)
            {
                rm[i] = bm[i];
                for (j = 0; j < n; j++)
                {
                    rm[i] -= a[i + j * m] * xn2[j];
                }
            }

            Console.WriteLine("");
            Console.WriteLine("  Given b is NOT in the range of A, solve");
            Console.WriteLine("    A * x = b");
            Console.WriteLine("  using the pseudoinverse:");
            Console.WriteLine("    x2 = A+ * b.");
            Console.WriteLine("");
            Console.WriteLine("  Norm of x2 = " + typeMethods.r8vec_norm_l2(n, xn2) + "");
            Console.WriteLine("  Norm of residual = " + typeMethods.r8vec_norm_l2(m, rm) + "");

        }

        //
        //  A' * x = b, b is in the range of A'.
        //
        xm1 = UniformRNG.r8vec_uniform_01_new(m, ref seed);
        for (i = 0; i < m; i++)
        {
            xm1[i] = typeMethods.r8_nint(10.0 * xm1[i]);
        }

        bn = new double[n];
        for (i = 0; i < n; i++)
        {
            bn[i] = 0.0;
            for (j = 0; j < m; j++)
            {
                bn[i] += a[j + i * m] * xm1[j];
            }
        }

        xm2 = new double[m];
        for (i = 0; i < m; i++)
        {
            xm2[i] = 0.0;
            for (j = 0; j < n; j++)
            {
                xm2[i] += a_pseudo[j + i * n] * bn[j];
            }
        }

        rn = new double[n];
        for (i = 0; i < n; i++)
        {
            rn[i] = bn[i];
            for (j = 0; j < m; j++)
            {
                rn[i] -= a[j + i * m] * xm2[j];
            }
        }

        Console.WriteLine("");
        Console.WriteLine("  Given:");
        Console.WriteLine("    b = A' * x1");
        Console.WriteLine("  so that b is in the range of A', solve");
        Console.WriteLine("    A' * x = b");
        Console.WriteLine("  using the pseudoinverse:");
        Console.WriteLine("    x2 = A+' * b.");
        Console.WriteLine("");
        Console.WriteLine("  Norm of x1 = " + typeMethods.r8vec_norm_l2(m, xm1) + "");
        Console.WriteLine("  Norm of x2 = " + typeMethods.r8vec_norm_l2(m, xm2) + "");
        Console.WriteLine("  Norm of residual = " + typeMethods.r8vec_norm_l2(n, rn) + "");

        //
        //  A' * x = b, b is not in the range of A'.

        if (m < n)
        {
            Console.WriteLine("");
            Console.WriteLine("  For M < N, most systems A'*x=b will not be");
            Console.WriteLine("  exactly and uniquely solvable, except in the");
            Console.WriteLine("  least squares sense.");
            Console.WriteLine("");
            Console.WriteLine("  Here is an example:");

            bn = UniformRNG.r8vec_uniform_01_new(n, ref seed);

            xm2 = new double[m];
            for (i = 0; i < m; i++)
            {
                xm2[i] = 0.0;
                for (j = 0; j < n; j++)
                {
                    xm2[i] += a_pseudo[j + i * n] * bn[j];
                }
            }

            rn = new double[n];
            for (i = 0; i < n; i++)
            {
                rn[i] = bn[i];
                for (j = 0; j < m; j++)
                {
                    rn[i] -= a[j + i * m] * xm2[j];
                }
            }

            Console.WriteLine("");
            Console.WriteLine("  Given b is NOT in the range of A', solve");
            Console.WriteLine("    A' * x = b");
            Console.WriteLine("  using the pseudoinverse:");
            Console.WriteLine("    x2 = A+ * b.");
            Console.WriteLine("");
            Console.WriteLine("  Norm of x2 = " + typeMethods.r8vec_norm_l2(m, xm2) + "");
            Console.WriteLine("  Norm of residual = " + typeMethods.r8vec_norm_l2(n, rn) + "");

        }

    }

    public static void pseudo_product_test(int m, int n, double[] a, double[] a_pseudo)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PSEUDO_PRODUCT_TEST examines pseudoinverse products.
        //
        //  Discussion:
        //
        //    Given an MxN matrix A, and its pseudoinverse A+, we must have
        //
        //      A+ * A * A+ = A+
        //      A * A+ * A = A
        //      ( A * A+ )' = A * A+ (MxM symmetry)
        //      ( A+ * A )' = A+ * A (NxN symmetry)
        //
        //    If M <= N, A * A+ may be "interesting" (equal to or "like" the identity),
        //    if N <= M, A+ * A may be "interesting" (equal to or "like" the identity).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 September 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, N, the number of rows and columns in the matrix A.
        //
        //    Input, double A[M*N], the matrix whose singular value
        //    decomposition we are investigating.
        //
        //    Input, double A_PSEUDO[N*M], the pseudo_inverse of A.
        //
    {
        double[] bmm;
        double[] bmn;
        double[] bnm;
        double[] bnn;
        double dif1;
        double dif2;
        double dif3;
        double dif4;
        int i;
        int j;
        int k;

        Console.WriteLine("");
        Console.WriteLine("PSEUDO_PRODUCT_TEST");
        Console.WriteLine("");
        Console.WriteLine("  The following relations MUST hold:");
        Console.WriteLine("");
        Console.WriteLine("   A  * A+ * A  = A");
        Console.WriteLine("   A+ * A  * A+ = A+");
        Console.WriteLine(" ( A  * A+ ) is MxM symmetric;");
        Console.WriteLine(" ( A+ * A  ) is NxN symmetric");
        //
        //  Compute A * A+ * A.
        //
        bnn = new double[n * n];
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                bnn[i + j * n] = 0.0;
                for (k = 0; k < m; k++)
                {
                    bnn[i + j * n] += a_pseudo[i + k * n] * a[k + j * m];
                }
            }
        }

        bmn = new double[m * n];

        for (i = 0; i < m; i++)
        {
            for (j = 0; j < n; j++)
            {
                bmn[i + j * m] = 0.0;
                for (k = 0; k < n; k++)
                {
                    bmn[i + j * m] += a[i + k * m] * bnn[k + j * n];
                }
            }
        }

        dif1 = typeMethods.r8mat_dif_fro(m, n, a, bmn);

        //
        //  Compute A+ * A * A+.
        //
        bmm = new double[m * m];
        for (i = 0; i < m; i++)
        {
            for (j = 0; j < m; j++)
            {
                bmm[i + j * m] = 0.0;
                for (k = 0; k < n; k++)
                {
                    bmm[i + j * m] += a[i + k * m] * a_pseudo[k + j * n];
                }
            }
        }

        bnm = new double[n * m];

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < m; j++)
            {
                bnm[i + j * n] = 0.0;
                for (k = 0; k < m; k++)
                {
                    bnm[i + j * n] += a_pseudo[i + k * n] * bmm[k + j * m];
                }
            }
        }

        dif2 = typeMethods.r8mat_dif_fro(n, m, a_pseudo, bnm);

        //
        //  Compute norm of A * A+ - (A * A+)'.
        //
        bmm = new double[m * m];
        for (i = 0; i < m; i++)
        {
            for (j = 0; j < m; j++)
            {
                bmm[i + j * m] = 0.0;
                for (k = 0; k < n; k++)
                {
                    bmm[i + j * m] += a[i + k * m] * a_pseudo[k + j * n];
                }
            }
        }

        dif3 = 0.0;
        for (j = 0; j < m; j++)
        {
            for (i = 0; i < m; i++)
            {
                dif3 += Math.Pow(bmm[i + j * m] - bmm[j + i * m], 2);
            }
        }

        dif3 = Math.Sqrt(dif3);

        //
        //  Compute norm of A+ * A - (A+ * A)'
        //
        bnn = new double[n * n];
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                bnn[i + j * n] = 0.0;
                for (k = 0; k < m; k++)
                {
                    bnn[i + j * n] += a_pseudo[i + k * n] * a[k + j * m];
                }
            }
        }

        dif4 = 0.0;
        for (j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
            {
                dif4 += Math.Pow(bnn[i + j * n] - bnn[j + i * n], 2);
            }
        }

        dif4 = Math.Sqrt(dif4);

        //
        //  Report.
        //
        Console.WriteLine("");
        Console.WriteLine("  Here are the Frobenius norms of the errors");
        Console.WriteLine("  in these relationships:");
        Console.WriteLine("");
        Console.WriteLine("   A  * A+ * A  = A            " + dif1 + "");
        Console.WriteLine("   A+ * A  * A+ = A+           " + dif2 + "");
        Console.WriteLine(" ( A  * A+ ) is MxM symmetric; " + dif3 + "");
        Console.WriteLine(" ( A+ * A  ) is NxN symmetric; " + dif4 + "");

        Console.WriteLine("");
        Console.WriteLine("  In some cases, the matrix A * A+");
        Console.WriteLine("  may be interesting (if M <= N, then");
        Console.WriteLine("  it MIGHT look like the identity.)");
        Console.WriteLine("");
        bmm = new double[m * m];
        for (i = 0; i < m; i++)
        {
            for (j = 0; j < m; j++)
            {
                bmm[i + j * m] = 0.0;
                for (k = 0; k < n; k++)
                {
                    bmm[i + j * m] += a[i + k * m] * a_pseudo[k + j * n];
                }
            }
        }

        typeMethods.r8mat_print(m, m, bmm, "  A * A+:");


        Console.WriteLine("");
        Console.WriteLine("  In some cases, the matrix A+ * A");
        Console.WriteLine("  may be interesting (if N <= M, then");
        Console.WriteLine("  it MIGHT look like the identity.)");
        Console.WriteLine("");

        bnn = new double[n * n];
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                bnn[i + j * n] = 0.0;
                for (k = 0; k < m; k++)
                {
                    bnn[i + j * n] += a_pseudo[i + k * n] * a[k + j * m];
                }
            }
        }

        typeMethods.r8mat_print(n, n, bnn, "  A+ * A");

    }
}