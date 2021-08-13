using System;
using Burkardt.Types;

namespace Burkardt.Sequence
{
    public static class Jacobi
    {

        public static double[] jacobi1(int n, double[] a, double[] b, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    JACOBI1 carries out one step of the Jacobi iteration.
            //
            //  Discussion:
            //
            //    The linear system A*x=b is to be solved.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 January 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of the matrix.
            //
            //    Input, double A[N,N], the matrix.
            //
            //    Input, double B[N], the right hand side.
            //
            //    Input, double X[N], the current solution estimate.
            //
            //    Output, double JACOBI1[N], the solution estimate updated by
            //    one step of the Jacobi iteration.
            //
        {
            int i;
            int j;
            double[] x_new;

            x_new = new double[n];

            for (i = 0; i < n; i++)
            {
                x_new[i] = b[i];
                for (j = 0; j < n; j++)
                {
                    if (j != i)
                    {
                        x_new[i] = x_new[i] - a[i + j * n] * x[j];
                    }
                }

                x_new[i] = x_new[i] / a[i + i * n];
            }

            return x_new;
        }

        public static void jacobi_eigenvalue(int n, double[] a, int it_max, ref double[] v,
                ref double[] d, ref int it_num, ref int rot_num)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    JACOBI_EIGENVALUE carries out the Jacobi eigenvalue iteration.
            //
            //  Discussion:
            //
            //    This function computes the eigenvalues and eigenvectors of a
            //    real symmetric matrix, using Rutishauser's modfications of the classical
            //    Jacobi rotation method with threshold pivoting. 
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    15 July 2013
            //
            //  Author:
            //
            //    C++ version by John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of the matrix.
            //
            //    Input, double A[N*N], the matrix, which must be square, real,
            //    and symmetric.
            //
            //    Input, int IT_MAX, the maximum number of iterations.
            //
            //    Output, double V[N*N], the matrix of eigenvectors.
            //
            //    Output, double D[N], the eigenvalues, in descending order.
            //
            //    Output, int &IT_NUM, the total number of iterations.
            //
            //    Output, int &ROT_NUM, the total number of rotations.
            //
        {
            double[] bw;
            double c;
            double g;
            double gapq;
            double h;
            int i;
            int j;
            int k;
            int l;
            int m;
            int p;
            int q;
            double s;
            double t;
            double tau;
            double term;
            double termp;
            double termq;
            double theta;
            double thresh;
            double w;
            double[] zw;

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < n; i++)
                {
                    if (i == j)
                    {
                        v[i + j * n] = 1.0;
                    }
                    else
                    {
                        v[i + j * n] = 0.0;
                    }
                }
            }

            for (i = 0; i < n; i++)
            {
                d[i] = a[i + i * n];
            }

            bw = new double[n];
            zw = new double[n];

            for (i = 0; i < n; i++)
            {
                bw[i] = d[i];
                zw[i] = 0.0;
            }

            it_num = 0;
            rot_num = 0;

            while (it_num < it_max)
            {
                it_num = it_num + 1;
                //
                //  The convergence threshold is based on the size of the elements in
                //  the strict upper triangle of the matrix.
                //
                thresh = 0.0;
                for (j = 0; j < n; j++)
                {
                    for (i = 0; i < j; i++)
                    {
                        thresh = thresh + a[i + j * n] * a[i + j * n];
                    }
                }

                thresh = Math.Sqrt(thresh) / (double) (4 * n);

                if (thresh == 0.0)
                {
                    break;
                }

                for (p = 0; p < n; p++)
                {
                    for (q = p + 1; q < n; q++)
                    {
                        gapq = 10.0 * Math.Abs(a[p + q * n]);
                        termp = gapq + Math.Abs(d[p]);
                        termq = gapq + Math.Abs(d[q]);
                        //
                        //  Annihilate tiny offdiagonal elements.
                        //
                        if (4 < it_num &&
                            termp == Math.Abs(d[p]) &&
                            termq == Math.Abs(d[q]))
                        {
                            a[p + q * n] = 0.0;
                        }
                        //
                        //  Otherwise, apply a rotation.
                        //
                        else if (thresh <= Math.Abs(a[p + q * n]))
                        {
                            h = d[q] - d[p];
                            term = Math.Abs(h) + gapq;

                            if (term == Math.Abs(h))
                            {
                                t = a[p + q * n] / h;
                            }
                            else
                            {
                                theta = 0.5 * h / a[p + q * n];
                                t = 1.0 / (Math.Abs(theta) + Math.Sqrt(1.0 + theta * theta));
                                if (theta < 0.0)
                                {
                                    t = -t;
                                }
                            }

                            c = 1.0 / Math.Sqrt(1.0 + t * t);
                            s = t * c;
                            tau = s / (1.0 + c);
                            h = t * a[p + q * n];
                            //
                            //  Accumulate corrections to diagonal elements.
                            //
                            zw[p] = zw[p] - h;
                            zw[q] = zw[q] + h;
                            d[p] = d[p] - h;
                            d[q] = d[q] + h;

                            a[p + q * n] = 0.0;
                            //
                            //  Rotate, using information from the upper triangle of A only.
                            //
                            for (j = 0; j < p; j++)
                            {
                                g = a[j + p * n];
                                h = a[j + q * n];
                                a[j + p * n] = g - s * (h + g * tau);
                                a[j + q * n] = h + s * (g - h * tau);
                            }

                            for (j = p + 1; j < q; j++)
                            {
                                g = a[p + j * n];
                                h = a[j + q * n];
                                a[p + j * n] = g - s * (h + g * tau);
                                a[j + q * n] = h + s * (g - h * tau);
                            }

                            for (j = q + 1; j < n; j++)
                            {
                                g = a[p + j * n];
                                h = a[q + j * n];
                                a[p + j * n] = g - s * (h + g * tau);
                                a[q + j * n] = h + s * (g - h * tau);
                            }

                            //
                            //  Accumulate information in the eigenvector matrix.
                            //
                            for (j = 0; j < n; j++)
                            {
                                g = v[j + p * n];
                                h = v[j + q * n];
                                v[j + p * n] = g - s * (h + g * tau);
                                v[j + q * n] = h + s * (g - h * tau);
                            }

                            rot_num = rot_num + 1;
                        }
                    }
                }

                for (i = 0; i < n; i++)
                {
                    bw[i] = bw[i] + zw[i];
                    d[i] = bw[i];
                    zw[i] = 0.0;
                }
            }

            //
            //  Restore upper triangle of input matrix.
            //
            for (j = 0; j < n; j++)
            {
                for (i = 0; i < j; i++)
                {
                    a[i + j * n] = a[j + i * n];
                }
            }

            //
            //  Ascending sort the eigenvalues and eigenvectors.
            //
            for (k = 0; k < n - 1; k++)
            {
                m = k;
                for (l = k + 1; l < n; l++)
                {
                    if (d[l] < d[m])
                    {
                        m = l;
                    }
                }

                if (m != k)
                {
                    t = d[m];
                    d[m] = d[k];
                    d[k] = t;
                    for (i = 0; i < n; i++)
                    {
                        w = v[i + m * n];
                        v[i + m * n] = v[i + k * n];
                        v[i + k * n] = w;
                    }
                }
            }
        }









    }
}