using System;

namespace Burkardt.Ball;

public static class Grid
{
    public static double[] ball_grid(int n, double r, double[] c, int ng)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BALL_GRID computes grid points inside a ball.
        //
        //  Discussion:
        //
        //    The grid is defined by specifying the radius and center of the ball,
        //    and the number of subintervals N into which the horizontal radius
        //    should be divided.  Thus, a value of N = 2 will result in 5 points
        //    along that horizontal line.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 November 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of subintervals.
        //
        //    Input, double R, the radius of the ball.
        //
        //    Input, double C[3], the coordinates of the center of the ball.
        //
        //    Input, int NG, the number of grid points, as determined by
        //    BALL_GRID_COUNT.
        //
        //    Output, double BALL_GRID[3*NG], the grid points inside the ball.
        //
    {
        double[] bg;
        int i;
        int j;
        int k;
        int p;
        double x;
        double y;
        double z;

        bg = new double[3 * ng];

        p = 0;

        for (i = 0; i <= n; i++)
        {
            x = c[0] + r * (2 * i) / (2 * n + 1);
            for (j = 0; j <= n; j++)
            {
                y = c[1] + r * (2 * j) / (2 * n + 1);
                for (k = 0; k <= n; k++)
                {
                    z = c[2] + r * (2 * k) / (2 * n + 1);

                    if (r * r < Math.Pow(x - c[0], 2)
                        + Math.Pow(y - c[1], 2)
                        + Math.Pow(z - c[2], 2))
                    {
                        break;
                    }


                    bg[0 + p * 3] = x;
                    bg[1 + p * 3] = y;
                    bg[2 + p * 3] = z;
                    p += 1;

                    switch (i)
                    {
                        case > 0:
                            bg[0 + p * 3] = 2.0 * c[0] - x;
                            bg[1 + p * 3] = y;
                            bg[2 + p * 3] = z;
                            p += 1;
                            break;
                    }

                    switch (j)
                    {
                        case > 0:
                            bg[0 + p * 3] = x;
                            bg[1 + p * 3] = 2.0 * c[1] - y;
                            bg[2 + p * 3] = z;
                            p += 1;
                            break;
                    }

                    switch (k)
                    {
                        case > 0:
                            bg[0 + p * 3] = x;
                            bg[1 + p * 3] = y;
                            bg[2 + p * 3] = 2.0 * c[2] - z;
                            p += 1;
                            break;
                    }

                    switch (i)
                    {
                        case > 0 when 0 < j:
                            bg[0 + p * 3] = 2.0 * c[0] - x;
                            bg[1 + p * 3] = 2.0 * c[1] - y;
                            bg[2 + p * 3] = z;
                            p += 1;
                            break;
                    }

                    switch (i)
                    {
                        case > 0 when 0 < k:
                            bg[0 + p * 3] = 2.0 * c[0] - x;
                            bg[1 + p * 3] = y;
                            bg[2 + p * 3] = 2.0 * c[2] - z;
                            p += 1;
                            break;
                    }

                    switch (j)
                    {
                        case > 0 when 0 < k:
                            bg[0 + p * 3] = x;
                            bg[1 + p * 3] = 2.0 * c[1] - y;
                            bg[2 + p * 3] = 2.0 * c[2] - z;
                            p += 1;
                            break;
                    }

                    switch (i)
                    {
                        case > 0 when 0 < j && 0 < k:
                            bg[0 + p * 3] = 2.0 * c[0] - x;
                            bg[1 + p * 3] = 2.0 * c[1] - y;
                            bg[2 + p * 3] = 2.0 * c[2] - z;
                            p += 1;
                            break;
                    }
                }
            }
        }

        return bg;
    }

    public static int ball_grid_count(int n, double r, double[] c)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BALL_GRID computes grid points inside a ball.
        //
        //  Discussion:
        //
        //    The grid is defined by specifying the radius and center of the ball,
        //    and the number of subintervals N into which the horizontal radius
        //    should be divided.  Thus, a value of N = 2 will result in 5 points
        //    along that horizontal line.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 November 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of subintervals.
        //
        //    Input, double R, the radius of the ball.
        //
        //    Input, double C[3], the coordinates of the center of the ball.
        //
        //    Output, int BALL_GRID_COUNT, the number of grid points inside the ball.
        //
    {
        int i;
        int j;
        int k;
        int ng;
        double x;
        double y;
        double z;

        ng = 0;

        for (i = 0; i <= n; i++)
        {
            x = c[0] + r * (2 * i) / (2 * n + 1);
            for (j = 0; j <= n; j++)
            {
                y = c[1] + r * (2 * j) / (2 * n + 1);
                for (k = 0; k <= n; k++)
                {
                    z = c[2] + r * (2 * k) / (2 * n + 1);

                    if (r * r < Math.Pow(x - c[0], 2)
                        + Math.Pow(y - c[1], 2)
                        + Math.Pow(z - c[2], 2))
                    {
                        break;
                    }

                    ng += 1;

                    switch (i)
                    {
                        case > 0:
                            ng += 1;
                            break;
                    }

                    switch (j)
                    {
                        case > 0:
                            ng += 1;
                            break;
                    }

                    switch (k)
                    {
                        case > 0:
                            ng += 1;
                            break;
                    }

                    switch (i)
                    {
                        case > 0 when 0 < j:
                            ng += 1;
                            break;
                    }

                    switch (i)
                    {
                        case > 0 when 0 < k:
                            ng += 1;
                            break;
                    }

                    switch (j)
                    {
                        case > 0 when 0 < k:
                            ng += 1;
                            break;
                    }

                    switch (i)
                    {
                        case > 0 when 0 < j && 0 < k:
                            ng += 1;
                            break;
                    }

                }
            }
        }

        return ng;
    }
}