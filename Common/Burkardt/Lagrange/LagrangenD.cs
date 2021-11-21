using Burkardt.ClenshawCurtisNS;
using Burkardt.OrderNS;
using Burkardt.Types;

namespace Burkardt.Lagrange;

public static class LagrangenD
{
    public static double[] lagrange_interp_nd_grid(int m, int[] n_1d, double[] a, double[] b,
            int nd )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGRANGE_INTERP_ND_GRID sets an M-dimensional Lagrange interpolant grid.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 September 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the spatial dimension.
        //
        //    Input, int N_1D[M], the order of the 1D rule to be used
        //    in each dimension.
        //
        //    Input, double A[M], B[M], the lower and upper limits.
        //
        //    Input, int ND, the number of points in the product grid.
        //
        //    Output, double LAGRANGE_INTERP_ND_GRID[M*ND], the points at which data 
        //    is to be sampled.
        //
    {
        int i;
        int j;
        typeMethods.r8vecDPData data = new();
        //
        //  Compute the data points.
        //
        double[] xd = new double[m * nd];

        for (j = 0; j < nd; j++)
        {
            for (i = 0; i < m; i++)
            {
                xd[i + j * m] = 0.0;
            }
        }

        for (i = 0; i < m; i++)
        {
            int n = n_1d[i];
            double[] x_1d = ClenshawCurtis.cc_compute_points(n);
            for (j = 0; j < n; j++)
            {
                x_1d[j] = 0.5 * ((1.0 - x_1d[j]) * a[i]
                                 + (1.0 + x_1d[j]) * b[i]);
            }

            typeMethods.r8vec_direct_product(ref data, i, n, x_1d, m, nd, ref xd);
        }

        return xd;
    }

    public static double[] lagrange_interp_nd_grid2(int m, int[] ind, double[] a, double[] b,
            int nd )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGRANGE_INTERP_ND_GRID2 sets an M-dimensional Lagrange interpolant grid.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 September 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the spatial dimension.
        //
        //    Input, int IND[M], the index or level of the 1D rule 
        //    to be used in each dimension.
        //
        //    Input, double A[M], B[M], the lower and upper limits.
        //
        //    Input, int ND, the number of points in the product grid.
        //
        //    Output, double LAGRANGE_INTERP_ND_GRID2[M*ND], the points at which data 
        //    was sampled.
        //
    {
        int i;
        int j;
        typeMethods.r8vecDPData data = new();
        //
        //  Compute the data points.
        //
        double[] xd = new double[m * nd];

        for (j = 0; j < nd; j++)
        {
            for (i = 0; i < m; i++)
            {
                xd[i + j * m] = 0.0;
            }
        }

        for (i = 0; i < m; i++)
        {
            int n = Order.order_from_level_135(ind[i]);
            double[] x_1d = ClenshawCurtis.cc_compute_points(n);
            for (j = 0; j < n; j++)
            {
                x_1d[j] = 0.5 * ((1.0 - x_1d[j]) * a[i]
                                 + (1.0 + x_1d[j]) * b[i]);
            }

            typeMethods.r8vec_direct_product(ref data, i, n, x_1d, m, nd, ref xd);
        }

        return xd;
    }

    public static int lagrange_interp_nd_size(int m, int[] n_1d)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGRANGE_INTERP_ND_SIZE sizes an M-dimensional Lagrange interpolant.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 September 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the spatial dimension.
        //
        //    Input, int N_1D[M], the order of the 1D rule to be used
        //    in each dimension.
        //
        //    Output, int LAGRANGE_INTERP_ND_SIZE, the number of points in the product grid.
        //
    {
        //
        //  Determine the number of data points.
        //
        int nd = typeMethods.i4vec_product(m, n_1d);

        return nd;
    }

    public static int lagrange_interp_nd_size2(int m, int[] ind)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGRANGE_INTERP_ND_SIZE2 sizes an M-dimensional Lagrange interpolant.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 September 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the spatial dimension.
        //
        //    Input, int IND[M], the index or level of the 1D rule 
        //    to be used in each dimension.
        //
        //    Output, int LAGRANGE_INTERP_ND_SIZE2, the number of points in the product grid.
        //
    {
        int i;
        //
        //  Determine the number of data points.
        //
        int nd = 1;
        for (i = 0; i < m; i++)
        {
            int n = Order.order_from_level_135(ind[i]);
            nd *= n;
        }

        return nd;
    }

    public static double[] lagrange_interp_nd_value(int m, int[] n_1d, double[] a, double[] b,
            int nd, double[] zd, int ni, double[] xi )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGRANGE_INTERP_ND_VALUE evaluates an ND Lagrange interpolant.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 September 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the spatial dimension.
        //
        //    Input, int N_1D[M], the order of the 1D rule to be used
        //    in each dimension.
        //
        //    Input, double A[M], B[M], the lower and upper limits.
        //
        //    Input, int ND, the number of points in the product grid.
        //
        //    Input, double ZD[ND], the function evaluated at the points XD.
        //
        //    Input, int NI, the number of points at which the 
        //    interpolant is to be evaluated.
        //
        //    Input, double XI[M*NI], the points at which the interpolant is 
        //    to be evaluated.
        //
        //    Output, double LAGRANGE_INTERP_ND_VALUE[NI], the interpolant evaluated 
        //    at the points XI.
        //
    {
        int j;
        typeMethods.r8vecDPData data = new();

        double[] w = new double[nd];
        double[] zi = new double[ni];

        for (j = 0; j < ni; j++)
        {
            int i;
            for (i = 0; i < nd; i++)
            {
                w[i] = 1.0;
            }

            for (i = 0; i < m; i++)
            {
                int n = n_1d[i];
                double[] x_1d = ClenshawCurtis.cc_compute_points(n);
                int k;
                for (k = 0; k < n; k++)
                {
                    x_1d[k] = 0.5 * ((1.0 - x_1d[k]) * a[i]
                                     + (1.0 + x_1d[k]) * b[i]);
                }

                double[] value = Lagrange1D.lagrange_base_1d(n, x_1d, 1, xi, xiIndex: + i + j * m);
                typeMethods.r8vec_direct_product2(ref data, i, n, value, m, nd, ref w);
            }

            zi[j] = typeMethods.r8vec_dot_product(nd, w, zd);
        }
            
        return zi;
    }

    public static double[] lagrange_interp_nd_value2(int m, int[] ind, double[] a, double[] b,
            int nd, double[] zd, int ni, double[] xi )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGRANGE_INTERP_ND_VALUE2 evaluates an ND Lagrange interpolant.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 September 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the spatial dimension.
        //
        //    Input, int IND[M], the index or level of the 1D rule 
        //    to be used in each dimension.
        //
        //    Input, double A[M], B[M], the lower and upper limits.
        //
        //    Input, int ND, the number of points in the product grid.
        //
        //    Input, double ZD[ND], the function evaluated at the points XD.
        //
        //    Input, int NI, the number of points at which the 
        //    interpolant is to be evaluated.
        //
        //    Input, double XI[M*NI], the points at which the interpolant 
        //    is to be evaluated.
        //
        //    Output, double ZI[NI], the interpolant evaluated at the 
        //    points XI.
        //
    {
        int j;
        typeMethods.r8vecDPData data = new();

        double[] w = new double[nd];
        double[] zi = new double[ni];

        for (j = 0; j < ni; j++)
        {
            int i;
            for (i = 0; i < nd; i++)
            {
                w[i] = 1.0;
            }

            for (i = 0; i < m; i++)
            {
                int n = Order.order_from_level_135(ind[i]);
                double[] x_1d = ClenshawCurtis.cc_compute_points(n);
                int k;
                for (k = 0; k < n; k++)
                {
                    x_1d[k] = 0.5 * ((1.0 - x_1d[k]) * a[i]
                                     + (1.0 + x_1d[k]) * b[i]);
                }

                double[] value = Lagrange1D.lagrange_base_1d(n, x_1d, 1, xi, xiIndex: + i + j * m);
                typeMethods.r8vec_direct_product2(ref data, i, n, value, m, nd, ref w);
            }

            zi[j] = typeMethods.r8vec_dot_product(nd, w, zd);
        }

        return zi;
    }

}