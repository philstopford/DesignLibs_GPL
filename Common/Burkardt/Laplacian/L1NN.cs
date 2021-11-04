using System;

namespace Burkardt.Laplacian
{
    public static class L1NN
    {
        public static double[] l1nn_apply(int n, double h, double[] u)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    L1NN_APPLY applies the 1D NN Laplacian to a vector.
            //
            //  Discussion:
            //
            //    The N grid points are assumed to be evenly spaced by H.
            //
            //    For N = 5, the discrete Laplacian with left Neumann and right Neumann
            //    boundary conditions on [0,6] has the matrix form L:
            //
            //       1 -1  0  0  0
            //      -1  2 -1  0  0
            //       0 -1  2 -1  0
            //       0  0 -1  2 -1
            //       0  0  0 -1  1
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    30 October 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of points.
            //    N must be at least 3.
            //
            //    Input, double H, the spacing between points.
            //
            //    Input, double U[N], the value at each point.
            //
            //    Output, double L1NN_APPLY[N], the Laplacian evaluated at each point.
            //
        {
            int i;
            double[] lu;

            if (n < 3)
            {
                Console.WriteLine("");
                Console.WriteLine("L1NN_APPLY - Fatal error!");
                Console.WriteLine("  N < 3.");
                return null;
            }

            lu = new double[n];

            i = 0;
            lu[i] = (u[i] - u[i + 1]) / h / h;
            for (i = 1; i < n - 1; i++)
            {
                lu[i] = (-u[i - 1] + 2.0 * u[i] - u[i + 1]) / h / h;
            }

            i = n - 1;
            lu[i] = (-u[i - 1] + u[i]) / h / h;

            return lu;
        }

        public static double[] l1nn_cholesky(int n, double h)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    L1NN_CHOLESKY computes the Cholesky factor of the 1D NN Laplacian.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    30 October 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of points.
            //    N must be at least 3.
            //
            //    Input, double H, the spacing between points.
            //
            //    Output, double L1NN_CHOLESKY[N*N], the Cholesky factor.
            //
        {
            double[] c;
            int i;
            int j;

            if (n < 3)
            {
                Console.WriteLine("");
                Console.WriteLine("L1NN_CHOLESKY - Fatal error!");
                Console.WriteLine("  N < 3.");
                return null;
            }

            c = new double[n * n];

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < n; i++)
                {
                    c[i + j * n] = 0.0;
                }
            }

            for (i = 0; i < n - 1; i++)
            {
                c[i + i * n] = +1.0;
                c[i + (i + 1) * n] = -1.0;
            }

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < n; i++)
                {
                    c[i + j * n] = c[i + j * n] / h;
                }
            }

            return c;
        }

        public static void l1nn_eigen(int n, double h, ref double[] v, ref double[] lambda )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    L1NN_EIGEN returns eigeninformation for the 1D NN Laplacian.
        //
        //  Discussion:
        //
        //    The grid points are assumed to be evenly spaced by H.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 October 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of points.
        //
        //    Input, double H, the spacing between points.
        //
        //    Output, double V[N*N], the eigenvectors.
        //
        //    Output, double LAMBDA[N], the eigenvalues.
        //
        {
            int i;
            double i_r8;
            int j;
            double j_r8;
            double n_r8;
            
            double theta;

            n_r8 = (double) (n);

            for (j = 0; j < n; j++)
            {
                j_r8 = (double) (j + 1);
                theta = Math.PI * (j_r8 - 1.0) / (2.0 * n_r8);
                lambda[j] = Math.Pow(2.0 * Math.Sin(theta) / h, 2);
                if (j == 0)
                {
                    for (i = 0; i < n; i++)
                    {
                        v[i + j * n] = Math.Sqrt(n_r8);
                    }
                }
                else
                {
                    for (i = 0; i < n; i++)
                    {
                        i_r8 = (double) (i + 1);
                        theta = Math.PI * (i_r8 - 0.5) * (j_r8 - 1.0) / n_r8;
                        v[i + j * n] = Math.Sqrt(2.0 / n_r8) * Math.Cos(theta);
                    }
                }
            }
        }

        public static double[] l1nn(int n, double h)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    L1NN stores the 1D NN Laplacian as a full matrix.
            //
            //  Discussion:
            //
            //    The N grid points are assumed to be evenly spaced by H.
            //
            //    For N = 5, the discrete Laplacian with Neumann boundary conditions
            //    at both ends of [0,6] has the matrix form L:
            //
            //       1 -1  0  0  0
            //      -1  2 -1  0  0
            //       0 -1  2 -1  0
            //       0  0 -1  2 -1
            //       0  0  0 -1  1
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    30 October 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of points.
            //    N must be at least 3.
            //
            //    Input, double H, the spacing between points.
            //
            //    Output, double L1NN[N*N], the Laplacian matrix.
            //
        {
            int i;
            int j;
            double[] l;

            if (n < 3)
            {
                Console.WriteLine("");
                Console.WriteLine("L1NN - Fatal error!");
                Console.WriteLine("  N < 3.");
                return null;
            }

            l = new double[n * n];

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < n; i++)
                {
                    l[i + j * n] = 0.0;
                }
            }

            i = 0;
            l[i + i * n] = 1.0 / h / h;
            l[i + (i + 1) * n] = -1.0 / h / h;

            for (i = 1; i < n - 1; i++)
            {
                l[i + (i - 1) * n] = -1.0 / h / h;
                l[i + i * n] = 2.0 / h / h;
                l[i + (i + 1) * n] = -1.0 / h / h;
            }

            i = n - 1;
            l[i + (i - 1) * n] = -1.0 / h / h;
            l[i + i * n] = 1.0 / h / h;

            return l;
        }

        public static void l1nn_lu(int n, double h, ref double[] l, ref double[] u )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    L1NN_LU computes the LU factors of the 1D NN Laplacian.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 November 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of points.
        //    N must be at least 3.
        //
        //    Input, double H, the spacing between points.
        //
        //    Output, double L[N*N], U[N*N], the LU factors.
        //
        {
            int i;
            int j;

            if (n < 3)
            {
                Console.WriteLine("");
                Console.WriteLine("L1NN_LU - Fatal error!");
                Console.WriteLine("  N < 3.");
                return;
            }

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < n; i++)
                {
                    l[i + j * n] = 0.0;
                }
            }

            for (i = 0; i < n; i++)
            {
                l[i + i * n] = 1.0;
            }

            for (i = 1; i < n; i++)
            {
                l[i + (i - 1) * n] = -1.0;
            }

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < n; i++)
                {
                    l[i + j * n] = l[i + j * n] / h;
                }
            }

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < n; i++)
                {
                    u[i + j * n] = 0.0;
                }
            }

            for (i = 0; i < n - 1; i++)
            {
                u[i + i * n] = 1.0;
            }

            i = n - 1;
            u[i + i * n] = 0.0;

            for (i = 0; i < n - 1; i++)
            {
                u[i + (i + 1) * n] = -1.0;
            }

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < n; i++)
                {
                    u[i + j * n] = u[i + j * n] / h;
                }
            }
        }
    }
}