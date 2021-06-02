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
    }
}