using System;
using Burkardt.Types;

namespace Burkardt.Ellipsoid;

public static class Grid
{
    public static double[] ellipsoid_grid(int n, double[] r, double[] c, int ng )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ELLIPSOID_GRID generates the grid points inside an ellipsoid.
        //
        //  Discussion:
        //
        //    The ellipsoid is specified as
        //
        //      ( ( X - C1 ) / R1 )^2 
        //    + ( ( Y - C2 ) / R2 )^2 
        //    + ( ( Z - C3 ) / R3 )^2 = 1
        //
        //    The user supplies a number N.  There will be N+1 grid points along
        //    the shortest axis.
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
        //    Input, double R[3], the half axis lengths.
        //
        //    Input, double C[3], the center of the ellipsoid.
        //
        //    Input, int NG, the number of grid points.
        //
        //    Output, double XYZ[3*NG], the grid point coordinates.
        //
    {
        double h;
        int k;
        int ni;
        int nj;
        int nk;
        double[] p = new double[3 * 8];

        int ng2 = 0;

        double[] xyz = new double[3 * ng];

        double rmin = typeMethods.r8vec_min(3, r);

        if (Math.Abs(r[0] - rmin) <= double.Epsilon)
        {
            h = 2.0 * r[0] / (2 * n + 1);
            ni = n;
            nj = (int)(Math.Ceiling(r[1] / r[0]) * n);
            nk = (int)(Math.Ceiling(r[2] / r[0]) * n);
        }
        else if (Math.Abs(r[1] - rmin) <= double.Epsilon)
        {
            h = 2.0 * r[1] / (2 * n + 1);
            nj = n;
            ni = (int)(Math.Ceiling(r[0] / r[1]) * n);
            nk = (int)(Math.Ceiling(r[2] / r[1]) * n);
        }
        else
        {
            h = 2.0 * r[2] / (2 * n + 1);
            nk = n;
            ni = (int)(Math.Ceiling(r[0] / r[2]) * n);
            nj = (int)(Math.Ceiling(r[1] / r[2]) * n);
        }

        for (k = 0; k <= nk; k++)
        {
            double z = c[2] + k * h;
            int j;
            for (j = 0; j <= nj; j++)
            {
                double y = c[1] + j * h;
                int i;
                for (i = 0; i <= ni; i++)
                {
                    double x = c[0] + i * h;
                    //
                    //  If we have left the ellipsoid, the I loop is completed.
                    //
                    if (1.0 < Math.Pow((x - c[0]) / r[0], 2)
                        + Math.Pow((y - c[1]) / r[1], 2)
                        + Math.Pow((z - c[2]) / r[2], 2))
                    {
                        break;
                    }

                    //
                    //  At least one point is generated, but more possible by symmetry.
                    //
                    int np = 0;
                    p[0 + np * 3] = x;
                    p[1 + np * 3] = y;
                    p[2 + np * 3] = z;

                    np = 1;

                    int m;
                    switch (i)
                    {
                        case > 0:
                        {
                            for (m = 0; m < np; m++)
                            {
                                p[0 + (np + m) * 3] = 2.0 * c[0] - p[0 + m * 3];
                                p[1 + (np + m) * 3] = p[1 + m * 3];
                                p[2 + (np + m) * 3] = p[2 + m * 3];
                            }

                            np = 2 * np;
                            break;
                        }
                    }

                    switch (j)
                    {
                        case > 0:
                        {
                            for (m = 0; m < np; m++)
                            {
                                p[0 + (np + m) * 3] = p[0 + m * 3];
                                p[1 + (np + m) * 3] = 2.0 * c[1] - p[1 + m * 3];
                                p[2 + (np + m) * 3] = p[2 + m * 3];
                            }

                            np = 2 * np;
                            break;
                        }
                    }

                    switch (k)
                    {
                        case > 0:
                        {
                            for (m = 0; m < np; m++)
                            {
                                p[0 + (np + m) * 3] = p[0 + m * 3];
                                p[1 + (np + m) * 3] = p[1 + m * 3];
                                p[2 + (np + m) * 3] = 2.0 * c[2] - p[2 + m * 3];
                            }

                            np = 2 * np;
                            break;
                        }
                    }

                    for (m = 0; m < np; m++)
                    {
                        int ii;
                        for (ii = 0; ii < 3; ii++)
                        {
                            xyz[ii + (ng2 + m) * 3] = p[ii + m * 3];
                        }
                    }

                    ng2 += np;
                }
            }
        }

        return xyz;
    }

    public static int ellipsoid_grid_count(int n, double[] r, double[] c )

        //****************************************************************************80
        //
        // ELLIPSOID_GRID_COUNT counts the grid points inside an ellipsoid.
        //
        //  Discussion:
        //
        //    The ellipsoid is specified as
        //
        //      ( ( X - C1 ) / R1 )^2 
        //    + ( ( Y - C2 ) / R2 )^2 
        //    + ( ( Z - C3 ) / R3 )^2 = 1
        //
        //    The user supplies a number N.  There will be N+1 grid points along
        //    the shortest axis.
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
        //    Input, double R[3], the half axis lengths.
        //
        //    Input, double C[3], the center of the ellipsoid.
        //
        //    Output, int ELLIPSOID_GRID_COUNT, the number of grid points.
        //
    {
        double h;
        int k;
        int ni;
        int nj;
        int nk;
        int np = 0;

        int ng = 0;

        double rmin = typeMethods.r8vec_min(3, r);

        if (Math.Abs(r[0] - rmin) <= double.Epsilon)
        {
            h = 2.0 * r[0] / (2 * n + 1);
            ni = n;
            nj = (int)(Math.Ceiling(r[1] / r[0]) * n);
            nk = (int)(Math.Ceiling(r[2] / r[0]) * n);
        }
        else if (Math.Abs(r[1] - rmin) <= double.Epsilon)
        {
            h = 2.0 * r[1] / (2 * n + 1);
            nj = n;
            ni = (int)(Math.Ceiling(r[0] / r[1]) * n);
            nk = (int)(Math.Ceiling(r[2] / r[1]) * n);
        }
        else
        {
            h = 2.0 * r[2] / (2 * n + 1);
            nk = n;
            ni = (int)(Math.Ceiling(r[0] / r[2]) * n);
            nj = (int)(Math.Ceiling(r[1] / r[2]) * n);
        }

        for (k = 0; k <= nk; k++)
        {
            double z = c[2] + k * h;
            int j;
            for (j = 0; j <= nj; j++)
            {
                double y = c[1] + j * h;
                int i;
                for (i = 0; i <= ni; i++)
                {
                    double x = c[0] + i * h;
                    //
                    //  If we have left the ellipsoid, the I loop is completed.
                    //
                    if (1.0 < Math.Pow((x - c[0]) / r[0], 2)
                        + Math.Pow((y - c[1]) / r[1], 2)
                        + Math.Pow((z - c[2]) / r[2], 2))
                    {
                        break;
                    }

                    np = k switch
                    {
                        > 0 => 2 * np,
                        _ => j switch
                        {
                            > 0 => 2 * np,
                            _ => i switch
                            {
                                > 0 => 2 * np,
                                //
                                //  At least one point is generated, but more possible by symmetry.
                                //
                                _ => 1
                            }
                        }
                    };

                    ng += np;
                }
            }
        }

        return ng;
    }
}