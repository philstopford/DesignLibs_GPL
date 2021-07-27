namespace Burkardt.LineNS
{
    public static class Grid
    {
        public static double[] line_grid(int n, double a, double b, int c)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    LINE_GRID: grid points over the interior of a line segment in 1D.
            //
            //  Discussion:
            //
            //    In 1D, a grid is created using N points.
            //
            //    Over the interval [A,B], we have 5 choices for grid centering:
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
            //
            //    Input, double A, B, the endpoints.
            //
            //    Input, int C, the grid centering.
            //    1 <= C <= 5.
            //
            //    Output, double LINE_GRID[N], the points.
            //
        {
            int j;
            double[] x;

            x = new double[n];
            //
            //  Create the 1D grids in each dimension.
            //
            for (j = 0; j < n; j++)
            {
                if (c == 1)
                {
                    if (n == 1)
                    {
                        x[j] = 0.5 * (a + b);
                    }
                    else
                    {
                        x[j] = ((double) (n - j - 1) * a
                                + (double) (j) * b)
                               / (double) (n - 1);
                    }
                }
                else if (c == 2)
                {
                    x[j] = ((double) (n - j) * a
                            + (double) (j + 1) * b)
                           / (double) (n + 1);
                }
                else if (c == 3)
                {
                    x[j] = ((double) (n - j) * a
                            + (double) (j - 2) * b)
                           / (double) (n);
                }
                else if (c == 4)
                {
                    x[j] = ((double) (n - j - 1) * a
                            + (double) (j + 1) * b)
                           / (double) (n);
                }
                else if (c == 5)
                {
                    x[j] = ((double) (2 * n - 2 * j - 1) * a
                            + (double) (2 * j + 1) * b)
                           / (double) (2 * n);
                }
            }

            return x;
        }
    }
}