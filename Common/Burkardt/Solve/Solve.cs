using System;

namespace Burkardt
{
    public static class Solve
    {
        public static double[] r8rmat_fs_new(int n, double[][] a, double[] b)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8RMAT_FS_NEW factors and solves an R8RMAT system with one right hand side.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 July 2014
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
            //    Input, double **A, the coefficient matrix of the linear system.
            //
            //    Input, double B[N], the right hand side of the linear system.
            //
            //    Output, double R8RMAT_FS_NEW[N], the solution of the linear system.
            //
        {
            double[][] a2;
            int i;
            int j;
            int k;
            int p;
            double t;
            double[] x;
            //
            //  Create a copy of the matrix.
            //
            a2 = new double[n][];

            for (i = 0; i < n; i++)
            {
                a2[i] = new double[n];
            }

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    a2[i][j] = a[i][j];
                }
            }

            //
            //  Create X and set it to B.
            //
            x = new double[n];
            for (i = 0; i < n; i++)
            {
                x[i] = b[i];
            }

            for (k = 0; k < n; k++)
            {
                //
                //  Find the maximum element in column I.
                //
                p = k;

                for (i = k + 1; i < n; i++)
                {
                    if (Math.Abs(a2[p][k]) < Math.Abs(a2[i][k]))
                    {
                        p = i;
                    }
                }

                if (a2[p][k] == 0.0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("R8RMAT_FS_NEW - Fatal error!");
                    Console.WriteLine("  Zero pivot on step " + k + "");
                    return new double[1];
                }

                //
                //  Switch rows K and P.
                //
                if (k != p)
                {
                    for (j = 0; j < n; j++)
                    {
                        t = a2[k][j];
                        a2[k][j] = a2[p][j];
                        a2[p][j] = t;
                    }

                    t = x[k];
                    x[k] = x[p];
                    x[p] = t;
                }

                //
                //  Scale the pivot row.
                //
                t = a2[k][k];
                a2[k][k] = 1.0;
                for (j = k + 1; j < n; j++)
                {
                    a2[k][j] = a2[k][j] / t;
                }

                x[k] = x[k] / t;
                //
                //  Use the pivot row to eliminate lower entries in that column.
                //
                for (i = k + 1; i < n; i++)
                {
                    if (a2[i][k] != 0.0)
                    {
                        t = -a2[i][k];
                        a2[i][k] = 0.0;
                        for (j = k + 1; j < n; j++)
                        {
                            a2[i][j] = a2[i][j] + t * a2[k][j];
                        }

                        x[i] = x[i] + t * x[k];
                    }
                }
            }

            //
            //  Back solve.
            //
            for (j = n - 1; 1 <= j; j--)
            {
                for (i = 0; i < j; i++)
                {
                    x[i] = x[i] - a2[i][j] * x[j];
                }
            }

            return x;
        }

        public static double[] r8ge_fs ( int n, ref double[] a, double[] b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8GE_FS factors and solves a R8GE system.
        //
        //  Discussion:
        //
        //    The R8GE storage format is used for a "general" M by N matrix.  
        //    A physical storage space is made for each logical entry.  The two 
        //    dimensional logical array is mapped to a vector, in which storage is 
        //    by columns.
        //
        //    R8GE_FS does not save the LU factors of the matrix, and hence cannot
        //    be used to efficiently solve multiple linear systems, or even to
        //    factor A at one time, and solve a single linear system at a later time.
        //
        //    R8GE_FS uses partial pivoting, but no pivot vector is required.
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
        //    Output, double R8GE_FS[N], the solution of the linear system.
        //
        {
            int i;
            int ipiv;
            int j;
            int jcol;
            double piv;
            double t;
            double[] x;

            x = new double[n];

            for (i = 0; i < n; i++)
            {
                x[i] = b[i];
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

                if (piv == 0.0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("R8GE_FS - Fatal error!");
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
                    a[jcol - 1 + (j - 1) * n] = a[jcol - 1 + (j - 1) * n] / t;
                }

                x[jcol - 1] = x[jcol - 1] / t;
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
                            a[i - 1 + (j - 1) * n] = a[i - 1 + (j - 1) * n] + t * a[jcol - 1 + (j - 1) * n];
                        }

                        x[i - 1] = x[i - 1] + t * x[jcol - 1];
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
                    x[i - 1] = x[i - 1] - a[i - 1 + (jcol - 1) * n] * x[jcol - 1];
                }
            }

            return x;
        }

        public static double[] r8ge_fss_new(int n, double[] a, int nb, double[] b )

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

                if (piv == 0.0)
                {
                    Console.WriteLine();
                    Console.WriteLine("R8GE_FSS_NEW - Fatal error!");
                    Console.WriteLine("  Zero pivot on step " + jcol);
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
                    a[jcol - 1 + (j - 1) * n] = a[jcol - 1 + (j - 1) * n] / t;
                }

                for (j = 0; j < nb; j++)
                {
                    x[jcol - 1 + j * n] = x[jcol - 1 + j * n] / t;
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
                            a[i - 1 + (j - 1) * n] = a[i - 1 + (j - 1) * n] + t * a[jcol - 1 + (j - 1) * n];
                        }

                        for (j = 0; j < nb; j++)
                        {
                            x[i - 1 + j * n] = x[i - 1 + j * n] + t * x[jcol - 1 + j * n];
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
                        x[i - 1 + j * n] = x[i - 1 + j * n] - a[i - 1 + (jcol - 1) * n] * x[jcol - 1 + j * n];
                    }
                }
            }

            return x;
        }
    }
}