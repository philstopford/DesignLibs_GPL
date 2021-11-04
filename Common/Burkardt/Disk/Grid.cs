using System;

namespace Burkardt.Disk
{
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
            double[] cg;
            int i;
            int j;
            int p;
            double x;
            double y;

            cg = new double[2 * ng];

            p = 0;

            for (j = 0; j <= n; j++)
            {
                i = 0;
                x = c[0];
                y = c[1] + r * (double) (2 * j) / (double) (2 * n + 1);

                cg[0 + 2 * p] = x;
                cg[1 + 2 * p] = y;
                p = p + 1;

                if (0 < j)
                {
                    cg[0 + 2 * p] = x;
                    cg[1 + 2 * p] = 2.0 * c[1] - y;
                    p = p + 1;
                }

                for (;;)
                {
                    i = i + 1;
                    x = c[0] + r * (double) (2 * i) / (double) (2 * n + 1);

                    if (r * r < Math.Pow(x - c[0], 2) + Math.Pow(y - c[1], 2))
                    {
                        break;
                    }

                    cg[0 + 2 * p] = x;
                    cg[1 + 2 * p] = y;
                    p = p + 1;
                    cg[0 + 2 * p] = 2.0 * c[0] - x;
                    cg[1 + 2 * p] = y;
                    p = p + 1;

                    if (0 < j)
                    {
                        cg[0 + 2 * p] = x;
                        cg[1 + 2 * p] = 2.0 * c[1] - y;
                        p = p + 1;
                        cg[0 + 2 * p] = 2.0 * c[0] - x;
                        cg[1 + 2 * p] = 2.0 * c[1] - y;
                        p = p + 1;
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
            int i;
            int j;
            int ng;
            double x;
            double y;

            ng = 0;

            for (j = 0; j <= n; j++)
            {
                i = 0;
                x = c[0];
                y = c[1] + r * (double) (2 * j) / (double) (2 * n + 1);
                ng = ng + 1;

                if (0 < j)
                {
                    ng = ng + 1;
                }

                for (;;)
                {
                    i = i + 1;
                    x = c[0] + r * (double) (2 * i) / (double) (2 * n + 1);

                    if (r * r < Math.Pow(x - c[0], 2) + Math.Pow(y - c[1], 2))
                    {
                        break;
                    }

                    ng = ng + 1;
                    ng = ng + 1;
                    if (0 < j)
                    {
                        ng = ng + 1;
                        ng = ng + 1;
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
            double[] g;
            double gr;
            double gt;
            int i;
            double phi;
            
            double r0;

            r0 = r / Math.Sqrt((double) (n) - 0.5);
            phi = (1.0 + Math.Sqrt(5.0)) / 2.0;

            g = new double[2 * n];

            for (i = 0; i < n; i++)
            {
                gr = r0 * Math.Sqrt((double) (i + 1) - 0.5);
                gt = 2.0 * Math.PI * (double) (i + 1) / phi;
                g[0 + i * 2] = c[0] + gr * Math.Cos(gt);
                g[1 + i * 2] = c[1] + gr * Math.Sin(gt);
            }

            return g;
        }
    }
}