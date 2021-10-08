using Burkardt.TriangulationNS;
using Burkardt.Types;

namespace Burkardt.PiecewiseLinear
{
    public static class Interp2D
    {
        public static double[] pwl_interp_2d(int nxd, int nyd, double[] xd, double[] yd, double[] zd,
                int ni, double[] xi, double[] yi)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PWL_INTERP_2D: piecewise linear interpolant to data defined on a 2D grid.
            //
            //  Discussion:
            //
            //    Thanks to Adam Hirst for pointing out an error in the formula that
            //    chooses the interpolation triangle, 04 February 2018.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    04 February 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NXD, NYD, the number of X and Y data values.
            //
            //    Input, double XD[NXD], YD[NYD], the sorted X and Y data.
            //
            //    Input, double ZD[NXD*NYD}, the Z data.
            //
            //    Input, int NI, the number of interpolation points.
            //
            //    Input, double XI[NI], YI[NI], the coordinates of the
            //    interpolation points.
            //
            //    Output, double PWL_INTERP_2D[NI], the value of the interpolant.
            //
        {
            double alpha;
            double beta;
            double det;
            double dxa;
            double dxb;
            double dxi;
            double dya;
            double dyb;
            double dyi;
            double gamma;
            int i;
            int j;
            int k;
            double[] zi;

            zi = new double[ni];

            for (k = 0; k < ni; k++)
            {
                //
                //  For interpolation point (xi(k),yi(k)), find data intervals I and J so that:
                //
                //    xd(i) <= xi(k) <= xd(i+1),
                //    yd(j) <= yi(k) <= yd(j+1).
                //
                //  But if the interpolation point is not within a data interval, 
                //  assign the dummy interpolant value zi(k) = infinity.
                //
                i = typeMethods.r8vec_bracket5(nxd, xd, xi[k]);
                if (i == -1)
                {
                    zi[k] = typeMethods.r8_huge();
                    continue;
                }

                j = typeMethods.r8vec_bracket5(nyd, yd, yi[k]);
                if (j == -1)
                {
                    zi[k] = typeMethods.r8_huge();
                    continue;
                }

                //
                //  The rectangular cell is arbitrarily split into two triangles.
                //  The linear interpolation formula depends on which triangle 
                //  contains the data point.
                //
                //    (I,J+1)--(I+1,J+1)
                //      |\       |
                //      | \      |
                //      |  \     |
                //      |   \    |
                //      |    \   |
                //      |     \  |
                //    (I,J)---(I+1,J)
                //
                if (yi[k] < yd[j + 1] + (yd[j] - yd[j + 1]) * (xi[k] - xd[i]) / (xd[i + 1] - xd[i]))
                {
                    dxa = xd[i + 1] - xd[i];
                    dya = yd[j] - yd[j];

                    dxb = xd[i] - xd[i];
                    dyb = yd[j + 1] - yd[j];

                    dxi = xi[k] - xd[i];
                    dyi = yi[k] - yd[j];

                    det = dxa * dyb - dya * dxb;

                    alpha = (dxi * dyb - dyi * dxb) / det;
                    beta = (dxa * dyi - dya * dxi) / det;
                    gamma = 1.0 - alpha - beta;

                    zi[k] = alpha * zd[i + 1 + j * nxd] + beta * zd[i + (j + 1) * nxd] + gamma * zd[i + j * nxd];
                }
                else
                {
                    dxa = xd[i] - xd[i + 1];
                    dya = yd[j + 1] - yd[j + 1];

                    dxb = xd[i + 1] - xd[i + 1];
                    dyb = yd[j] - yd[j + 1];

                    dxi = xi[k] - xd[i + 1];
                    dyi = yi[k] - yd[j + 1];

                    det = dxa * dyb - dya * dxb;

                    alpha = (dxi * dyb - dyi * dxb) / det;
                    beta = (dxa * dyi - dya * dxi) / det;
                    gamma = 1.0 - alpha - beta;

                    zi[k] = alpha * zd[i + (j + 1) * nxd] + beta * zd[i + 1 + j * nxd] +
                            gamma * zd[i + 1 + (j + 1) * nxd];
                }
            }

            return zi;
        }
        
        public static double[] pwl_interp_2d_scattered_value ( int nd, double[] xyd, double[] zd, 
        int t_num, int[] t, int[] t_neighbor, int ni, double[] xyi )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PWL_INTERP_2D_SCATTERED_VALUE evaluates a 2d interpolant of scattered data
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    25 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ND, the number of data points.
        //
        //    Input, double XYD[2*ND], the data point coordinates.
        //
        //    Input, double ZD[ND], the data values.
        //
        //    Input, int T_NUM, the number of triangles.
        //
        //    Input, int T[3*T_NUM], the triangle information.
        //
        //    Input, int T_NEIGHBOR[3*T_NUM], the triangle neighbors.
        //
        //    Input, int NI, the number of interpolation points.
        //
        //    Input, double XYI[2*NI], the interpolation point coordinates.
        //
        //    Output, double PWL_INTERP_2D_SCATTERED_VALUE[NI], the interpolated values.
        //
        {
            double alpha = 0;
            double beta = 0;
            double gamma = 0;
            int edge = 0;
            int i;
            int j = 0;
            int step_num = 0;
            double[] zi;

            zi = new double[ni];

            DelaunaySearchData data = new DelaunaySearchData();
            
            for ( i = 0; i < ni; i++ )
            {
                Search.triangulation_search_delaunay_a (ref data, nd, xyd, 3, t_num, t, t_neighbor, 
                    xyi, ref j, ref alpha, ref beta, ref gamma, ref edge, ref step_num, pIndex: +2*i);
                
                if ( j == -1 )
                {
                    zi[i] = -1.0;
                }
                
                zi[i % zi.Length] = alpha * zd[t[(t.Length + (0+j*3)) % t.Length] % zd.Length] 
                                    + beta  * zd[t[(t.Length + (1+j*3)) % t.Length] % zd.Length] 
                                    + gamma * zd[t[(t.Length + (2+j*3)) % t.Length] % zd.Length];
            }
            return zi;
        }
    }
}