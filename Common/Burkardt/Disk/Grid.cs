﻿using System;

namespace Burkardt.Disk;

public static class Grid
{
    public static double[] disk_grid(int n, double r, double[] c, int ng )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DISK_GRID computes grid points inside a disk.
        //
        //  Discussion:
        //
        //    The grid is defined by specifying the radius and center of the circle,
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
        //    09 November 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of subintervals.
        //
        //    Input, double R, the radius of the circle.
        //
        //    Input, double C[2], the coordinates of the center of the circle.
        //
        //    Input, int NG, the number of grid points, as determined by
        //    DISK_GRID_COUNT.
        //
        //    Output, double DISK_GRID[2*NG], the grid points inside the circle.
        //
    {
        int j;

        double[] cg = new double[2 * ng];

        int p = 0;

        for (j = 0; j <= n; j++)
        {
            int i = 0;
            double x = c[0];
            double y = c[1] + r * (2 * j) / (2 * n + 1);

            cg[0 + 2 * p] = x;
            cg[1 + 2 * p] = y;
            p += 1;

            switch (j)
            {
                case > 0:
                    cg[0 + 2 * p] = x;
                    cg[1 + 2 * p] = 2.0 * c[1] - y;
                    p += 1;
                    break;
            }

            for (;;)
            {
                i += 1;
                x = c[0] + r * (2 * i) / (2 * n + 1);

                if (r * r < Math.Pow(x - c[0], 2) + Math.Pow(y - c[1], 2))
                {
                    break;
                }

                cg[0 + 2 * p] = x;
                cg[1 + 2 * p] = y;
                p += 1;
                cg[0 + 2 * p] = 2.0 * c[0] - x;
                cg[1 + 2 * p] = y;
                p += 1;

                switch (j)
                {
                    case > 0:
                        cg[0 + 2 * p] = x;
                        cg[1 + 2 * p] = 2.0 * c[1] - y;
                        p += 1;
                        cg[0 + 2 * p] = 2.0 * c[0] - x;
                        cg[1 + 2 * p] = 2.0 * c[1] - y;
                        p += 1;
                        break;
                }
            }
        }

        return cg;
    }

    public static int disk_grid_count(int n, double r, double[] c )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DISK_GRID_COUNT counts the grid points inside a disk.
        //
        //  Discussion:
        //
        //    The grid is defined by specifying the radius and center of the circle,
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
        //    09 November 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of subintervals.
        //
        //    Input, double R, the radius of the circle.
        //
        //    Input, double C[2], the coordinates of the center of the circle.
        //
        //    Output, int DISK_GRID_COUNT, the number of grid points inside 
        //    the circle.
        //
    {
        int j;

        int ng = 0;

        for (j = 0; j <= n; j++)
        {
            int i = 0;
            double x = c[0];
            double y = c[1] + r * (2 * j) / (2 * n + 1);
            ng += 1;

            switch (j)
            {
                case > 0:
                    ng += 1;
                    break;
            }

            for (;;)
            {
                i += 1;
                x = c[0] + r * (2 * i) / (2 * n + 1);

                if (r * r < Math.Pow(x - c[0], 2) + Math.Pow(y - c[1], 2))
                {
                    break;
                }

                ng += 1;
                ng += 1;
                switch (j)
                {
                    case > 0:
                        ng += 1;
                        ng += 1;
                        break;
                }
            }
        }

        return ng;
    }

    public static double[] disk_grid_fibonacci(int n, double r, double[] c)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DISK_GRID_FIBONACCI computes Fibonacci grid points inside a disk.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 October 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Richard Swinbank, James Purser,
        //    Fibonacci grids: A novel approach to global modelling,
        //    Quarterly Journal of the Royal Meteorological Society,
        //    Volume 132, Number 619, July 2006 Part B, pages 1769-1793.
        //
        //  Parameters:
        //
        //    Input, int N, the number of points desired.
        //
        //    Input, double R, the radius of the circle.
        //
        //    Input, double C[2], the coordinates of the center of the circle.
        //
        //    Output, double DISK_GRID_FIBONACCI[2*N], the grid points.
        //
    {
        int i;

        double r0 = r / Math.Sqrt(n - 0.5);
        double phi = (1.0 + Math.Sqrt(5.0)) / 2.0;

        double[] g = new double[2 * n];

        for (i = 0; i < n; i++)
        {
            double gr = r0 * Math.Sqrt(i + 1 - 0.5);
            double gt = 2.0 * Math.PI * (i + 1) / phi;
            g[0 + i * 2] = c[0] + gr * Math.Cos(gt);
            g[1 + i * 2] = c[1] + gr * Math.Sin(gt);
        }

        return g;
    }
}