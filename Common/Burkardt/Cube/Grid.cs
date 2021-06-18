using Burkardt.Types;

namespace Burkardt.Cube
{
    public static class Grid
    {
        public static double[] cube_grid(int n, int[] ns, double[] a, double[] b,
                int[] c)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CUBE_GRID: grid points over the interior of a cube in 3D.
            //
            //  Discussion:
            //
            //    In 3D, a logically rectangular grid is to be created.
            //    In the I-th dimension, the grid will use S(I) points.
            //    The total number of grid points is 
            //      N = product ( 1 <= I <= 3 ) S(I)
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
            //    31 August 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of points.
            //    N = product ( 1 <= I <= 3 ) NS(I).
            //
            //    Input, int NS[3], the number of points along 
            //    each dimension.
            //
            //    Input, double A[3], B[3], the endpoints for each dimension.
            //
            //    Input, int C[3], the grid centering for each dimension.
            //    1 <= C(*) <= 5.
            //
            //    Output, double CUBE_GRID[3*N] = X(3*S(0)*S(1)*S(2)), the points.
            //
        {
            int i;
            int j;
            int m = 3;
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