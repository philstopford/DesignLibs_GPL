using System;

namespace Burkardt.Ellipse;

public static class Grid
{
    public static double[] ellipse_grid(int n, double[] r, double[] c, int ng )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ELLIPSE_GRID generates grid points inside an ellipse.
        //
        //  Discussion:
        //
        //    The ellipse is specified as
        //
        //      ( ( X - C1 ) / R1 )^2 + ( ( Y - C2 ) / R2 )^2 = 1
        //
        //    The user supplies a number N.  There will be N+1 grid points along
        //    the shorter axis.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 November 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of subintervals.
        //
        //    Input, double R[2], the half axis lengths.
        //
        //    Input, double C[2], the center of the ellipse.
        //
        //    Input, int NG, the number of grid points inside the ellipse.
        //
        //    Output, double ELLIPSE_GRID[2*NG], the grid points.
        //
    {
        double h;
        int i;
        int j;
        int nj;
        int p;
        double x;
        double[] xy;
        double y;

        xy = new double[2 * ng];

        if (r[0] < r[1])
        {
            h = 2.0 * r[0] / (2 * n + 1);
            nj = (int)(Math.Ceiling(r[1] / r[0]) * n);
        }
        else
        {
            h = 2.0 * r[1] / (2 * n + 1);
            nj = n;
        }

        p = 0;

        for (j = 0; j <= nj; j++)
        {
            i = 0;
            x = c[0];
            y = c[1] + j * h;

            xy[0 + p * 2] = x;
            xy[1 + p * 2] = y;
            p += 1;

            switch (j)
            {
                case > 0:
                    xy[0 + p * 2] = x;
                    xy[1 + p * 2] = 2.0 * c[1] - y;
                    p += 1;
                    break;
            }

            for (;;)
            {
                i += 1;
                x = c[0] + i * h;

                if (1.0 < Math.Pow((x - c[0]) / r[0], 2)
                    + Math.Pow((y - c[1]) / r[1], 2))
                {
                    break;
                }

                xy[0 + p * 2] = x;
                xy[1 + p * 2] = y;
                p += 1;
                xy[0 + p * 2] = 2.0 * c[0] - x;
                xy[1 + p * 2] = y;
                p += 1;

                switch (j)
                {
                    case > 0:
                        xy[0 + p * 2] = x;
                        xy[1 + p * 2] = 2.0 * c[1] - y;
                        p += 1;
                        xy[0 + p * 2] = 2.0 * c[0] - x;
                        xy[1 + p * 2] = 2.0 * c[1] - y;
                        p += 1;
                        break;
                }
            }
        }

        return xy;
    }

    public static int ellipse_grid_count(int n, double[] r, double[] c )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ELLIPSE_GRID_COUNT counts the grid points inside an ellipse.
        //
        //  Discussion:
        //
        //    The ellipse is specified as
        //
        //      ( ( X - C1 ) / R1 )^2 + ( ( Y - C2 ) / R2 )^2 = 1
        //
        //    The user supplies a number N.  There will be N+1 grid points along
        //    the shorter axis.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 November 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of subintervals.
        //
        //    Input, double R[2], the half axis lengths.
        //
        //    Input, double C[2], the center of the ellipse.
        //
        //    Output, int ELLIPSE_GRID)_COUNT, the number of grid points inside 
        //    the ellipse.
        //
    {
        double h;
        int i;
        int j;
        int nj;
        int p;
        double x;
        double y;

        if (r[0] < r[1])
        {
            h = 2.0 * r[0] / (2 * n + 1);
            nj = (int)(Math.Ceiling(r[1] / r[0]) * n);
        }
        else
        {
            h = 2.0 * r[1] / (2 * n + 1);
            nj = n;
        }

        p = 0;

        for (j = 0; j <= nj; j++)
        {
            i = 0;
            x = c[0];
            y = c[1] + j * h;

            p += 1;

            switch (j)
            {
                case > 0:
                    p += 1;
                    break;
            }

            for (;;)
            {
                i += 1;
                x = c[0] + i * h;

                if (1.0 < Math.Pow((x - c[0]) / r[0], 2)
                    + Math.Pow((y - c[1]) / r[1], 2))
                {
                    break;
                }

                p += 1;
                p += 1;

                switch (j)
                {
                    case > 0:
                        p += 1;
                        p += 1;
                        break;
                }
            }
        }

        return p;
    }
}