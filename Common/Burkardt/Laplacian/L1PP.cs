using System;

namespace Burkardt.Laplacian;

public static class L1PP
{
    public static double[] l1pp_apply(int n, double h, double[] u)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    L1PP_APPLY applies the 1D PP Laplacian to a vector.
        //
        //  Discussion:
        //
        //    The N grid points are assumed to be evenly spaced by H.
        //
        //    For N = 5, the discrete Laplacian with periodic boundary conditions
        //    on [0,6] has the matrix form L:
        //
        //       2 -1  0  0 -1
        //      -1  2 -1  0  0
        //       0 -1  2 -1  0
        //       0  0 -1  2 -1
        //      -1  0  0 -1  2
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
        //    Output, double L1PP_APPLY[N], the Laplacian evaluated at each point.
        //
    {
        int i;
        double[] lu;

        switch (n)
        {
            case < 3:
                Console.WriteLine("");
                Console.WriteLine("L1PP_APPLY - Fatal error!");
                Console.WriteLine("  N < 3.");
                return null;
        }

        lu = new double[n];

        i = 0;
        lu[i] = (-u[n - 1] + 2.0 * u[i] - u[i + 1]) / h / h;
        for (i = 1; i < n - 1; i++)
        {
            lu[i] = (-u[i - 1] + 2.0 * u[i] - u[i + 1]) / h / h;
        }

        i = n - 1;
        lu[i] = (-u[i - 1] + 2.0 * u[i] - u[0]) / h / h;

        return lu;
    }

    public static double[] l1pp_cholesky(int n, double h)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    L1PP_CHOLESKY computes the Cholesky factor of the 1D PP Laplacian.
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
        //    Output, double L1PP_CHOLESKY[N*N], the Cholesky factor.
        //
    {
        double[] c;
        int i;
        double i_r8;
        int j;

        switch (n)
        {
            case < 3:
                Console.WriteLine("");
                Console.WriteLine("L1PP_CHOLESKY - Fatal error!");
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
            i_r8 = i + 1;
            c[i + i * n] = Math.Sqrt(i_r8 + 1.0) / Math.Sqrt(i_r8);
        }

        for (i = 0; i < n - 2; i++)
        {
            i_r8 = i + 1;
            c[i + (i + 1) * n] = -i_r8 / (i_r8 + 1.0) * Math.Sqrt(i_r8 + 1.0)
                                 / Math.Sqrt(i_r8);
        }

        for (i = 0; i < n - 2; i++)
        {
            i_r8 = i + 1;
            c[i + (n - 1) * n] = -1.0 / (i_r8 + 1.0) * Math.Sqrt(i_r8 + 1.0)
                                 / Math.Sqrt(i_r8);
        }

        i = n - 2;
        i_r8 = i + 1;
        c[i + (n - 1) * n] = -(double) n / (i_r8 + 1.0)
            * Math.Sqrt(i_r8 + 1.0) / Math.Sqrt(i_r8);

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
            {
                c[i + j * n] /= h;
            }
        }

        return c;
    }

    public static void l1pp_eigen(int n, double h, ref double[] v, ref double[] lambda )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    L1PP_EIGEN returns eigeninformation for the 1D PP Laplacian.
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
            
        double s;
        double theta;

        n_r8 = n;

        for (j = 0; j < n; j++)
        {
            j_r8 = j + 1;
            theta = (j % 2) switch
            {
                0 => Math.PI * (j_r8 - 1.0) / (2.0 * n_r8),
                _ => Math.PI * j_r8 / (2.0 * n_r8)
            };

            lambda[j] = Math.Pow(2.0 * Math.Sin(theta) / h, 2);

            switch (j % 2)
            {
                case 0 when j == 0:
                {
                    for (i = 0; i < n; i++)
                    {
                        v[i + j * n] = 1.0 / Math.Sqrt(n_r8);
                    }

                    break;
                }
                case 0:
                {
                    for (i = 0; i < n; i++)
                    {
                        i_r8 = i + 1;
                        theta = Math.PI * (i_r8 - 0.5) * (j_r8 - 1.0) / n_r8;
                        v[i + j * n] = Math.Sqrt(2.0 / n_r8) * Math.Cos(theta);
                    }

                    break;
                }
                default:
                {
                    if (j == n - 1)
                    {
                        s = -1.0 / Math.Sqrt(n_r8);
                        for (i = 0; i < n; i++)
                        {
                            v[i + j * n] = s;
                            s = -s;
                        }
                    }
                    else
                    {
                        for (i = 0; i < n; i++)
                        {
                            i_r8 = i + 1;
                            theta = Math.PI * (i_r8 - 0.5) * j_r8 / n_r8;
                            v[i + j * n] = Math.Sqrt(2.0 / n_r8) * Math.Sin(theta);
                        }
                    }

                    break;
                }
            }

        }
    }

    public static double[] l1pp(int n, double h)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    L1PP stores the 1D PP Laplacian as a full matrix.
        //
        //  Discussion:
        //
        //    The N grid points are assumed to be evenly spaced by H.
        //
        //    For N = 5, the discrete Laplacian with periodic boundary conditions
        //    has the matrix form L:
        //
        //       2 -1  0  0 -1
        //      -1  2 -1  0  0
        //       0 -1  2 -1  0
        //       0  0 -1  2 -1
        //      -1  0  0 -1  2
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
        //    Output, double L1PP[N*N], the Laplacian matrix.
        //
    {
        int i;
        int j;
        double[] l;

        switch (n)
        {
            case < 3:
                Console.WriteLine("");
                Console.WriteLine("L1PP - Fatal error!");
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
        l[i + i * n] = 2.0 / h / h;
        l[i + (i + 1) * n] = -1.0 / h / h;
        l[i + (n - 1) * n] = -1.0 / h / h;

        for (i = 1; i < n - 1; i++)
        {
            l[i + (i - 1) * n] = -1.0 / h / h;
            l[i + i * n] = 2.0 / h / h;
            l[i + (i + 1) * n] = -1.0 / h / h;
        }

        i = n - 1;
        l[i + 0 * n] = -1.0 / h / h;
        l[i + (i - 1) * n] = -1.0 / h / h;
        l[i + i * n] = 2.0 / h / h;

        return l;
    }

    public static void l1pp_lu(int n, double h, ref double[] l, ref double[] u )

        //****************************************************************************80
        //
        //  Discussion:
        //
        //    L1PP_LU computes the LU factors of the 1D PP Laplacian.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 November 2013
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
        double i_r8;
        int j;

        switch (n)
        {
            case < 3:
                Console.WriteLine("");
                Console.WriteLine("L1PP_LU - Fatal error!");
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

        for (i = 1; i < n - 1; i++)
        {
            i_r8 = i + 1;
            l[i + (i - 1) * n] = -(i_r8 - 1.0) / i_r8;
            l[n - 1 + (i - 1) * n] = -1.0 / i_r8;
        }

        i = n - 1;
        l[i + (i - 1) * n] = -1.0;

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
            {
                l[i + j * n] /= h;
            }
        }

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
            {
                u[i + j * n] = 0.0;
            }
        }

        for (i = 0; i < n - 2; i++)
        {
            i_r8 = i + 1;
            u[i + i * n] = (i_r8 + 1.0) / i_r8;
            u[i + (i + 1) * n] = -1.0;
            u[i + (n - 1) * n] = -1.0 / i_r8;
        }

        i = n - 2;
        i_r8 = i + 1;
        u[i + i * n] = (i_r8 + 1.0) / i_r8;
        u[i + (i + 1) * n] = -(i_r8 + 1.0) / i_r8;

        i = n - 1;
        u[i + i * n] = 0.0;

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
            {
                u[i + j * n] /= h;
            }
        }
    }
}