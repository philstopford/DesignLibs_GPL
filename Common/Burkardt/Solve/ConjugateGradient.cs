using System;

namespace Burkardt.SolveNS;

public static class ConjugateGradient
{
    public static double[] solve_cg(int n, int[] diag, int nz_num, int[] ia, int[] ja,
            double[] a, double[] b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SOLVE_CG solves a linear system using the conjugate gradient method.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    25 January 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of nodes.
        //
        //    Input, int DIAG[N], contains for each index 0 <= I < N, the unique
        //    index J such that IA[J] = JA[J] = I.
        //
        //    Input, int NZ_NUM, the number of nonzero entries.
        //
        //    Input, int IA[NZ_NUM], JA[NZ_NUM], the row and column
        //    indices of the nonzero entries.
        //
        //    Input, double A[NZ_NUM], the nonzero entries of the matrix.
        //
        //    Input, double B[N], the right hand side.
        //
        //    Output, double SOLVE_CG[N], the solution of the linear system.
        //
    {
        double aii;
        int i;
        double rnrm2;
        ConjugateGradientData data = new();

        int it = 0;
        const int it_max = 100;
        const double tol = 1.0E-08;
        double bnrm2 = 0.0;
        for (i = 0; i < n; i++)
        {
            bnrm2 += b[i] * b[i];
        }

        bnrm2 = Math.Sqrt(bnrm2);

        double[] p = new double[n];
        double[] q = new double[n];
        double[] r = new double[n];
        double[] x = new double[n];
        double[] z = new double[n];

        for (i = 0; i < n; i++)
        {
            aii = a[diag[i]];
            x[i] = b[i] / aii;
        }

        Console.WriteLine("");
        Console.WriteLine("  Step        Residual");
        Console.WriteLine("");

        int job = 1;

        for (;;)
        {
            job = ConjugateGradientRC.cg_rc(ref data, n, b, ref x, ref r, ref z, ref p, ref q, ref job);
            //
            //  Compute q = A * p.
            //
            int k;
            int j;
            if (job == 1)
            {
                for (i = 0; i < n; i++)
                {
                    q[i] = 0.0;
                }

                for (k = 0; k < nz_num; k++)
                {
                    i = ia[k] - 1;
                    j = ja[k] - 1;
                    q[i] += a[k] * p[j];
                }
            }
            //
            //  Solve M * z = r.
            //
            else if (job == 2)
            {
                for (i = 0; i < n; i++)
                {
                    aii = a[diag[i]];
                    z[i] = r[i] / aii;
                }
            }
            //
            //  Compute r = r - A * x.
            //
            else if (job == 3)
            {
                for (k = 0; k < nz_num; k++)
                {
                    i = ia[k] - 1;
                    j = ja[k] - 1;
                    r[i] -= a[k] * x[j];
                }
            }
            //
            //  Stopping test.
            //
            else if (job == 4)
            {
                rnrm2 = 0.0;
                for (i = 0; i < n; i++)
                {
                    rnrm2 += r[i] * r[i];
                }

                rnrm2 = Math.Sqrt(rnrm2);

                if (bnrm2 == 0.0)
                {
                    if (rnrm2 <= tol)
                    {
                        break;
                    }
                }
                else
                {
                    if (rnrm2 <= tol * bnrm2)
                    {
                        break;
                    }
                }

                it += 1;
                Console.WriteLine("  " + it + "  " + rnrm2);

                if (it_max <= it)
                {
                    Console.WriteLine("");
                    Console.WriteLine("  Iteration limit exceeded.");
                    Console.WriteLine("  Terminating early.");
                    break;
                }
            }

            job = 2;
        }

        Console.WriteLine("");
        Console.WriteLine("  Number of iterations was " + it + "");
        Console.WriteLine("  Estimated error is " + rnrm2 + "");

        return x;
    }
}