namespace Burkardt.Lagrange;

public static class Lagrange2D
{
    public static double[] lagrange_interp_2d ( int mx, int my, double[] xd_1d, double[] yd_1d, 
            double[] zd, int ni, double[] xi, double[] yi )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGRANGE_INTERP_2D evaluates the Lagrange interpolant for a product grid.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 September 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int MX, MY, the polynomial degree in X and Y.
        //
        //    Input, double XD_1D[MX+1], YD_1D[MY+1], the 1D data locations.
        //
        //    Input, double ZD[(MX+1)*(MY+1)], the 2D array of data values.
        //
        //    Input, int NI, the number of 2D interpolation points.
        //
        //    Input, double XI[NI], YI[NI], the 2D interpolation points.
        //
        //    Output, double LAGRANGE_INTERP_2D[NI], the interpolated values.
        //
    {
        int k;

        double[] zi = new double[ni];

        for ( k = 0; k < ni; k++ )
        {
            int l = 0;
            zi[k] = 0.0;
            int j;
            for ( j = 0; j < my + 1; j++ )
            {
                int i;
                for ( i = 0; i < mx + 1; i++ )
                {
                    double lx = Lagrange1D.lagrange_basis_function_1d ( mx, xd_1d, i, xi[k] );
                    double ly = Lagrange1D.lagrange_basis_function_1d ( my, yd_1d, j, yi[k] );
                    zi[k] += zd[l] * lx * ly;
                    l += 1;
                }
            }
        }
        return zi;
    }
}