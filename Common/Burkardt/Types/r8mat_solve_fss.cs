using System;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static void r8mat_fss(int n, ref double[] a, int nb, ref double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_FSS factors and solves a system with multiple right hand sides.
        //
        //  Discussion:
        //
        //    This routine uses partial pivoting, but no pivot vector is required.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 November 2011
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
        //    Input/output, double X[N*NB], on input, the right hand sides of the
        //    linear systems.  On output, the solutions of the linear systems.
        //
    {
        int i;
        int ipiv;
        int j;
        int jcol;
        double piv;
        double t;

        for (jcol = 1; jcol <= n; jcol++)
        {
            //
            //  Find the maximum element in column I.
            //
            piv = Math.Abs(a[jcol - 1 + (jcol - 1) * n]);
            ipiv = jcol;
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
                    Console.WriteLine("R8MAT_FSS - Fatal error!");
                    Console.WriteLine("  Zero pivot on step " + jcol + "");
                    return;
            }

            //
            //  Switch rows JCOL and IPIV, and X.
            //
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
    }

    public static double[] r8mat_fss_new(int n, ref double[] a, int nb, double[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8MAT_FSS_NEW factors and solves a system with multiple right hand sides.
        //
        //  Discussion:
        //
        //    This routine uses partial pivoting, but no pivot vector is required.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 November 2011
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
        //    Output, double R8MAT_FSS_NEW[N*NB], the solutions of the linear systems.
        //
    {
        int i;
        int ipiv;
        int j;
        int jcol;
        double piv;
        double t;
        double[] x;

        x = new double[n * nb];

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
            piv = Math.Abs(a[jcol - 1 + (jcol - 1) * n]);
            ipiv = jcol;
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
                    Console.WriteLine("R8MAT_FSS_NEW - Fatal error!");
                    Console.WriteLine("  Zero pivot on step " + jcol + "");
                    return null;
            }

            //
            //  Switch rows JCOL and IPIV, and X.
            //
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
        
}