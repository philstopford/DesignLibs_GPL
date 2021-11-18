using System;
using Burkardt.Types;

namespace Burkardt.SolveNS;

public static class RootsRC
{
    public static void roots_rc(int n, double[] x, double[] fx, ref double ferr, ref double[] xnew,
            ref double[] q)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ROOTS_RC solves a system of nonlinear equations using reverse communication.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 January 2013
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Gaston Gonnet.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //     
        //    Gaston Gonnet,
        //    On the Structure of Zero Finders,
        //    BIT Numerical Mathematics,
        //    Volume 17, Number 2, June 1977, pages 170-183.
        //
        //  Parameters:
        //
        //    Input, int N, the number of equations.
        //
        //    Input, double X[N].  Before the first call, the user should
        //    set X to an initial guess or estimate for the root.  Thereafter, the input
        //    value of X should be the output value of XNEW from the previous call.
        //
        //    Input, double FX[N], the value of the function at XNEW.
        //
        //    Output, double &FERR, the function error, that is, the sum of
        //    the absolute values of the most recently computed function vector.
        //
        //    Output, double XNEW[N], a new point at which a function 
        //    value is requested.
        //
        //    Workspace, double Q[(2*N+2)*(N+2)].  Before the first call 
        //    for a given problem, the user must set Q(2*N+1,1) to 0.0.
        //
    {
        double damp;
        int i;
        int j;
        int jsma;
        int jsus;
        int lda;
        double sump;
        double t;

        lda = 2 * n + 2;

        ferr = 0.0;
        for (i = 0; i < n; i++)
        {
            ferr += Math.Abs(fx[i]);
        }

        switch (q[2 * n + 1 + 0 * lda])
        {
            //
            //  Initialization if Q(2*N+1,1) = 0.0.
            //
            case 0.0:
            {
                for (i = 1; i <= n; i++)
                {
                    for (j = 1; j <= n + 1; j++)
                    {
                        q[i - 1 + (j - 1) * lda] = 0.0;
                        q[i + (j - 1) * lda] = 0.0;
                    }

                    q[i - 1 + (i - 1) * lda] = 100.0;
                    q[i + n - 1 + (i - 1) * lda] = 1.0;
                }

                for (j = 1; j <= n; j++)
                {
                    q[2 * n + (j - 1) * lda] = typeMethods.r8_huge();
                }

                for (j = 1; j <= n; j++)
                {
                    q[2 * n + 1 + (j - 1) * lda] = n;
                }

                for (i = 1; i <= n; i++)
                {
                    q[i + n - 1 + n * lda] = x[i - 1];
                }

                for (i = 1; i <= n; i++)
                {
                    q[i - 1 + n * lda] = fx[i - 1];
                }

                q[2 * n + n * lda] = ferr;
                q[2 * n + 1 + n * lda] = 0.0;
                damp = 0.99;
                break;
            }
            default:
            {
                jsus = 1;
                for (i = 2; i <= n + 1; i++)
                {
                    if (2 * n <= q[2 * n + 1 + (i - 1) * lda])
                    {
                        q[2 * n + (i - 1) * lda] = typeMethods.r8_huge();
                    }

                    if (q[2 * n + 1 + (jsus - 1) * lda] < (double)(n + 3) / 2)
                    {
                        jsus = i;
                    }

                    if ((double)(n + 3) / 2 <= q[2 * n + 1 + (i - 1) * lda] &&
                        q[2 * n + (jsus - 1) * lda] < q[2 * n + (i - 1) * lda])
                    {
                        jsus = i;
                    }
                }

                for (i = 1; i <= n; i++)
                {
                    q[i + n - 1 + (jsus - 1) * lda] = x[i - 1];
                    q[i - 1 + (jsus - 1) * lda] = fx[i - 1];
                }

                q[2 * n + (jsus - 1) * lda] = ferr;
                q[2 * n + 1 + (jsus - 1) * lda] = 0;
                jsma = 1;
                damp = 0.0;

                for (j = 1; j <= n + 1; j++)
                {
                    if (typeMethods.r8_huge() / 10.0 < q[2 * n + (j - 1) * lda])
                    {
                        damp = 0.99;
                    }

                    if (q[2 * n + (j - 1) * lda] < q[2 * n + (jsma - 1) * lda])
                    {
                        jsma = j;
                    }
                }

                if (jsma != n + 1)
                {
                    for (i = 1; i <= 2 * n + 2; i++)
                    {
                        t = q[i - 1 + (jsma - 1) * lda];
                        q[i - 1 + (jsma - 1) * lda] = q[i - 1 + n * lda];
                        q[i - 1 + n * lda] = t;
                    }
                }

                break;
            }
        }

        for (i = 1; i <= n; i++)
        {
            q[i - 1 + (n + 1) * lda] = q[i - 1 + n * lda];
        }

        //
        //  Call the linear equation solver, which should not destroy the matrix 
        //  in Q(1:N,1:N), and should overwrite the solution into Q(1:N,N+2).
        //
        typeMethods.r8mat_fs(lda, n, q, ref q, xIndex: +(n + 1) * lda);

        sump = 0.0;
        for (i = 1; i <= n; i++)
        {
            sump += q[i - 1 + (n + 1) * lda];
        }

        switch (Math.Abs(1.0 - sump))
        {
            case <= 1.0E-10:
                Console.WriteLine("");
                Console.WriteLine("ROOTS_RC - Fatal error!");
                Console.WriteLine("  SUMP almost exactly 1.");
                Console.WriteLine("  SUMP = " + sump + "");
                return;
        }

        for (i = 1; i <= n; i++)
        {
            xnew[i - 1] = q[i + n - 1 + n * lda];
            for (j = 1; j <= n; j++)
            {
                xnew[i - 1] -= q[i + n - 1 + (j - 1) * lda] * q[j - 1 + (n + 1) * lda];
            }

            //
            //  If system not complete, damp the solution.
            //
            xnew[i - 1] = xnew[i - 1] / (1.0 - sump) * (1.0 - damp) + q[i + n - 1 + n * lda] * damp;
        }

        for (j = 1; j <= n + 1; j++)
        {
            q[2 * n + 1 + (j - 1) * lda] += 1.0;
        }

    }
}