﻿using Burkardt.Types;

namespace Burkardt.HyperGeometry.Hypercube
{
    public static class Grid
    {
        public static double[] hypercube_grid(int m, int n, int[] ns, double[] a, double[] b,
        int[] c )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HYPERCUBE_GRID: grid points over the interior of a hypercube in M dimensions.
        //
        //  Discussion:
        //
        //    In M dimensional space, a logically rectangular grid is to be created.
        //    In the I-th dimension, the grid will use S(I) points.
        //    The total number of grid points is 
        //      N = product ( 1 <= I <= M ) S(I)
        //
        //    Over the interval [A(i),B(i)], we have 5 choices for grid centering:
        //      1: 0,   1/3, 2/3, 1
        //      2: 1/5, 2/5, 3/5, 4/5
        //      3: 0,   1/4, 2/4, 3/4
        //      4: 1/4, 2/4, 3/4, 1
        //      5: 1/8, 3/8, 5/8, 7/8
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    30 August 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the spatial dimension.
        //
        //    Input, int N, the number of points.
        //    N = product ( 1 <= I <= M ) NS(I).
        //
        //    Input, int NS[M], the number of points along 
        //    each dimension.
        //
        //    Input, double A[M], B[M], the endpoints for each dimension.
        //
        //    Input, int C[M], the grid centering for each dimension.
        //    1 <= C(*) <= 5.
        //
        //    Output, double HYPERCUBE_GRID[M*N] = X(M*S(1),S(2),...,S(M)), the points.
        //
        {
            int i;
            int j;
            int s;
            double[] x;
            double[] xs;

            x = new double[m * n];
            //
            //  Create the 1D grids in each dimension.
            //
            for (i = 0; i < m; i++)
            {
                s = ns[i];

                xs = new double[s];

                for (j = 0; j < s; j++)
                {
                    if (c[i] == 1)
                    {
                        if (s == 1)
                        {
                            xs[j] = 0.5 * (a[i] + b[i]);
                        }
                        else
                        {
                            xs[j] = ((double) (s - j - 1) * a[i]
                                     + (double) (j) * b[i])
                                    / (double) (s - 1);
                        }
                    }
                    else if (c[i] == 2)
                    {
                        xs[j] = ((double) (s - j) * a[i]
                                 + (double) (j + 1) * b[i])
                                / (double) (s + 1);
                    }
                    else if (c[i] == 3)
                    {
                        xs[j] = ((double) (s - j) * a[i]
                                 + (double) (j - 2) * b[i])
                                / (double) (s);
                    }
                    else if (c[i] == 4)
                    {
                        xs[j] = ((double) (s - j - 1) * a[i]
                                 + (double) (j + 1) * b[i])
                                / (double) (s);
                    }
                    else if (c[i] == 5)
                    {
                        xs[j] = ((double) (2 * s - 2 * j - 1) * a[i]
                                 + (double) (2 * j + 1) * b[i])
                                / (double) (2 * s);
                    }
                }

                typeMethods.r8vec_direct_product(i, s, xs, m, n, ref x);

            }

            return x;
        }
    }
}