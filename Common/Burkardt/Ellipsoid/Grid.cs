using System;
using Burkardt.Types;

namespace Burkardt.Ellipsoid
{
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
            int ii;
            int i;
            int j;
            int k;
            int m;
            int ng2;
            int ni;
            int nj;
            int nk;
            int np;
            double[] p = new double[3 * 8];
            double rmin;
            double x;
            double[] xyz;
            double y;
            double z;

            ng2 = 0;

            xyz = new double[3 * ng];

            rmin = typeMethods.r8vec_min(3, r);

            if (r[0] == rmin)
            {
                h = 2.0 * r[0] / (double) (2 * n + 1);
                ni = n;
                nj = (int)(Math.Ceiling(r[1] / r[0]) * (double) (n));
                nk = (int)(Math.Ceiling(r[2] / r[0]) * (double) (n));
            }
            else if (r[1] == rmin)
            {
                h = 2.0 * r[1] / (double) (2 * n + 1);
                nj = n;
                ni = (int)(Math.Ceiling(r[0] / r[1]) * (double) (n));
                nk = (int)(Math.Ceiling(r[2] / r[1]) * (double) (n));
            }
            else
            {
                h = 2.0 * r[2] / (double) (2 * n + 1);
                nk = n;
                ni = (int)(Math.Ceiling(r[0] / r[2]) * (double) (n));
                nj = (int)(Math.Ceiling(r[1] / r[2]) * (double) (n));
            }

            for (k = 0; k <= nk; k++)
            {
                z = c[2] + (double) (k) * h;
                for (j = 0; j <= nj; j++)
                {
                    y = c[1] + (double) (j) * h;
                    for (i = 0; i <= ni; i++)
                    {
                        x = c[0] + (double) (i) * h;
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
                        np = 0;
                        p[0 + np * 3] = x;
                        p[1 + np * 3] = y;
                        p[2 + np * 3] = z;

                        np = 1;

                        if (0 < i)
                        {
                            for (m = 0; m < np; m++)
                            {
                                p[0 + (np + m) * 3] = 2.0 * c[0] - p[0 + m * 3];
                                p[1 + (np + m) * 3] = p[1 + m * 3];
                                p[2 + (np + m) * 3] = p[2 + m * 3];
                            }

                            np = 2 * np;
                        }

                        if (0 < j)
                        {
                            for (m = 0; m < np; m++)
                            {
                                p[0 + (np + m) * 3] = p[0 + m * 3];
                                p[1 + (np + m) * 3] = 2.0 * c[1] - p[1 + m * 3];
                                p[2 + (np + m) * 3] = p[2 + m * 3];
                            }

                            np = 2 * np;
                        }

                        if (0 < k)
                        {
                            for (m = 0; m < np; m++)
                            {
                                p[0 + (np + m) * 3] = p[0 + m * 3];
                                p[1 + (np + m) * 3] = p[1 + m * 3];
                                p[2 + (np + m) * 3] = 2.0 * c[2] - p[2 + m * 3];
                            }

                            np = 2 * np;
                        }

                        for (m = 0; m < np; m++)
                        {
                            for (ii = 0; ii < 3; ii++)
                            {
                                xyz[ii + (ng2 + m) * 3] = p[ii + m * 3];
                            }
                        }

                        ng2 = ng2 + np;
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
            int i;
            int j;
            int k;
            int ng;
            int ni;
            int nj;
            int nk;
            int np;
            double rmin;
            double x;
            double y;
            double z;

            ng = 0;

            rmin = typeMethods.r8vec_min(3, r);

            if (r[0] == rmin)
            {
                h = 2.0 * r[0] / (double) (2 * n + 1);
                ni = n;
                nj = (int)(Math.Ceiling(r[1] / r[0]) * (double) (n));
                nk = (int)(Math.Ceiling(r[2] / r[0]) * (double) (n));
            }
            else if (r[1] == rmin)
            {
                h = 2.0 * r[1] / (double) (2 * n + 1);
                nj = n;
                ni = (int)(Math.Ceiling(r[0] / r[1]) * (double) (n));
                nk = (int)(Math.Ceiling(r[2] / r[1]) * (double) (n));
            }
            else
            {
                h = 2.0 * r[2] / (double) (2 * n + 1);
                nk = n;
                ni = (int)(Math.Ceiling(r[0] / r[2]) * (double) (n));
                nj = (int)(Math.Ceiling(r[1] / r[2]) * (double) (n));
            }

            for (k = 0; k <= nk; k++)
            {
                z = c[2] + (double) (k) * h;
                for (j = 0; j <= nj; j++)
                {
                    y = c[1] + (double) (j) * h;
                    for (i = 0; i <= ni; i++)
                    {
                        x = c[0] + (double) (i) * h;
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
                        np = 1;
                        if (0 < i)
                        {
                            np = 2 * np;
                        }

                        if (0 < j)
                        {
                            np = 2 * np;
                        }

                        if (0 < k)
                        {
                            np = 2 * np;
                        }

                        ng = ng + np;
                    }
                }
            }

            return ng;
        }
    }
}