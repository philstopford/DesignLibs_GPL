using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static void r8mat_fs(int n, double[] a, ref double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_FS factors and solves a system with one right hand side.
            //
            //  Discussion:
            //
            //    This routine differs from R8MAT_FSS in two ways:
            //    * only one right hand side is allowed;
            //    * the input matrix A is not modified.
            //
            //    This routine uses partial pivoting, but no pivot vector is required.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    21 January 2013
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
            //    Input, double A[N*N], the coefficient matrix of the linear system.
            //
            //    Input/output, double X[N], on input, the right hand side of the
            //    linear system.  On output, the solution of the linear system.
            //
        {
            double[] a2;
            int i;
            int ipiv;
            int j;
            int jcol;
            double piv;
            double t;

            a2 = new double[n * n];

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < n; i++)
                {
                    a2[i + j * n] = a[i + j * n];
                }
            }

            for (jcol = 1; jcol <= n; jcol++)
            {
                //
                //  Find the maximum element in column I.
                //
                piv = Math.Abs(a2[jcol - 1 + (jcol - 1) * n]);
                ipiv = jcol;
                for (i = jcol + 1; i <= n; i++)
                {
                    if (piv < Math.Abs(a2[i - 1 + (jcol - 1) * n]))
                    {
                        piv = Math.Abs(a2[i - 1 + (jcol - 1) * n]);
                        ipiv = i;
                    }
                }

                if (piv == 0.0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("R8MAT_FS - Fatal error!");
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
                        t = a2[jcol - 1 + (j - 1) * n];
                        a2[jcol - 1 + (j - 1) * n] = a2[ipiv - 1 + (j - 1) * n];
                        a2[ipiv - 1 + (j - 1) * n] = t;
                    }

                    t = x[jcol - 1];
                    x[jcol - 1] = x[ipiv - 1];
                    x[ipiv - 1] = t;
                }

                //
                //  Scale the pivot row.
                //
                t = a2[jcol - 1 + (jcol - 1) * n];
                a2[jcol - 1 + (jcol - 1) * n] = 1.0;
                for (j = jcol + 1; j <= n; j++)
                {
                    a2[jcol - 1 + (j - 1) * n] = a2[jcol - 1 + (j - 1) * n] / t;
                }

                x[jcol - 1] = x[jcol - 1] / t;
                //
                //  Use the pivot row to eliminate lower entries in that column.
                //
                for (i = jcol + 1; i <= n; i++)
                {
                    if (a2[i - 1 + (jcol - 1) * n] != 0.0)
                    {
                        t = -a2[i - 1 + (jcol - 1) * n];
                        a2[i - 1 + (jcol - 1) * n] = 0.0;
                        for (j = jcol + 1; j <= n; j++)
                        {
                            a2[i - 1 + (j - 1) * n] = a2[i - 1 + (j - 1) * n] + t * a2[jcol - 1 + (j - 1) * n];
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
                    x[i - 1] = x[i - 1] - a2[i - 1 + (jcol - 1) * n] * x[jcol - 1];
                }
            }
        }

        public static double[] r8mat_fs_new(int n, double[] a, double[] b)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8MAT_FS_NEW factors and solves a system with one right hand side.
            //
            //  Discussion:
            //
            //    This routine differs from R8MAT_FSS_NEW in two ways:
            //    * only one right hand side is allowed;
            //    * the input matrix A is not modified.
            //
            //    This routine uses partial pivoting, but no pivot vector is required.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    21 January 2013
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
            //    Input, double A[N*N], the coefficient matrix of the linear system.
            //    On output, A is in unit upper triangular form, and
            //    represents the U factor of an LU factorization of the
            //    original coefficient matrix.
            //
            //    Input, double B[N], the right hand side of the linear system.
            //
            //    Output, double X[N], the solution of the linear system.
            //
        {
            double[] a2;
            int i;
            int ipiv;
            int j;
            int jcol;
            double piv;
            double t;
            double[] x;

            a2 = new double[n * n];

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < n; i++)
                {
                    a2[i + j * n] = a[i + j * n];
                }
            }

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
                piv = Math.Abs(a2[jcol - 1 + (jcol - 1) * n]);
                ipiv = jcol;
                for (i = jcol + 1; i <= n; i++)
                {
                    if (piv < Math.Abs(a2[i - 1 + (jcol - 1) * n]))
                    {
                        piv = Math.Abs(a2[i - 1 + (jcol - 1) * n]);
                        ipiv = i;
                    }
                }

                if (piv == 0.0)
                {
                    Console.WriteLine("");
                    Console.WriteLine("R8MAT_FS_NEW - Fatal error!");
                    Console.WriteLine("  Zero pivot on step " + jcol);
                    return new double[1];
                }

                //
                //  Switch rows JCOL and IPIV, and X.
                //
                if (jcol != ipiv)
                {
                    for (j = 1; j <= n; j++)
                    {
                        t = a2[jcol - 1 + (j - 1) * n];
                        a2[jcol - 1 + (j - 1) * n] = a2[ipiv - 1 + (j - 1) * n];
                        a2[ipiv - 1 + (j - 1) * n] = t;
                    }

                    t = x[jcol - 1];
                    x[jcol - 1] = x[ipiv - 1];
                    x[ipiv - 1] = t;
                }

                //
                //  Scale the pivot row.
                //
                t = a2[jcol - 1 + (jcol - 1) * n];
                a2[jcol - 1 + (jcol - 1) * n] = 1.0;
                for (j = jcol + 1; j <= n; j++)
                {
                    a2[jcol - 1 + (j - 1) * n] = a2[jcol - 1 + (j - 1) * n] / t;
                }

                x[jcol - 1] = x[jcol - 1] / t;
                //
                //  Use the pivot row to eliminate lower entries in that column.
                //
                for (i = jcol + 1; i <= n; i++)
                {
                    if (a2[i - 1 + (jcol - 1) * n] != 0.0)
                    {
                        t = -a2[i - 1 + (jcol - 1) * n];
                        a2[i - 1 + (jcol - 1) * n] = 0.0;
                        for (j = jcol + 1; j <= n; j++)
                        {
                            a2[i - 1 + (j - 1) * n] = a2[i - 1 + (j - 1) * n] + t * a2[jcol - 1 + (j - 1) * n];
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
                    x[i - 1] = x[i - 1] - a2[i - 1 + (jcol - 1) * n] * x[jcol - 1];
                }
            }

            return x;
        }
        
    }
}