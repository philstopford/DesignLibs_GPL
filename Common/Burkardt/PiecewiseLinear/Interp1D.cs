namespace Burkardt.PiecewiseLinear;

public static class Interp1D
{
    public static double[] pwl_basis_1d(int nd, double[] xd, int ni, double[] xi)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PWL_BASIS_1D evaluates a 1D piecewise linear basis function.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 July 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ND, the number of data points.
        //
        //    Input, double XD[ND], the data points.
        //
        //    Input, int NI, the number of interpolation points.
        //
        //    Input, double XI[NI], the interpolation points.
        //
        //    Output, double PW_BASIS_1D[NI*ND], the basis function at the 
        //    interpolation points.
        //
    {
        int i;
        int j;

        double[] bk = new double[ni * nd];

        for (j = 0; j < nd; j++)
        {
            for (i = 0; i < ni; i++)
            {
                bk[i + j * ni] = 0.0;
            }
        }

        switch (nd)
        {
            case 1:
            {
                for (j = 0; j < nd; j++)
                {
                    for (i = 0; i < ni; i++)
                    {
                        bk[i + j * ni] = 1.0;
                    }
                }

                return bk;
            }
        }

        for (i = 0; i < ni; i++)
        {
            for (j = 0; j < nd; j++)
            {
                double t;
                switch (j)
                {
                    case 0 when xi[i] <= xd[j]:
                        t = (xi[i] - xd[j]) / (xd[(j + 1) % xd.Length] - xd[j]);
                        bk[(i + j * ni) % bk.Length] = 1.0 - t;
                        break;
                    default:
                    {
                        if (j == nd - 1 && xd[j] <= xi[i])
                        {
                            t = (xi[i] - xd[(xd.Length + j - 1) % xd.Length]) / (xd[j] - xd[(j - 1) % xd.Length]);
                            bk[(i + j * ni) % bk.Length] = t;
                        }
                        else if (xd[(xd.Length + j - 1) % xd.Length] < xi[i] && xi[i] <= xd[j])
                        {
                            t = (xi[i] - xd[(xd.Length + j - 1) % xd.Length]) / (xd[j] - xd[(xd.Length + j - 1) % xd.Length]);
                            bk[(i + j * ni) % bk.Length] = t;
                        }
                        else if (xd[j] <= xi[i] && xi[i] < xd[(j + 1) % xd.Length])
                        {
                            t = (xi[i] - xd[j]) / (xd[(j + 1) % xd.Length] - xd[j]);
                            bk[(i + j * ni) % bk.Length] = 1.0 - t;
                        }

                        break;
                    }
                }
            }
        }

        return bk;
    }

    public static double[] pwl_value_1d(int nd, double[] xd, double[] yd, int ni, double[] xi)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PWL_VALUE_1D evaluates the piecewise linear interpolant.
        //
        //  Discussion:
        //
        //    The piecewise linear interpolant L(ND,XD,YD)(X) is the piecewise
        //    linear function which interpolates the data (XD(I),YD(I)) for I = 1
        //    to ND.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 September 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int ND, the number of data points.
        //    ND must be at least 1.
        //
        //    Input, double XD[ND], the data points.
        //
        //    Input, double YD[ND], the data values.
        //
        //    Input, int NI, the number of interpolation points.
        //
        //    Input, double XI[NI], the interpolation points.
        //
        //    Output, double PWL_VALUE_1D[NI], the interpolated values.
        //
    {
        int i;

        double[] yi = new double[ni];

        for (i = 0; i < ni; i++)
        {
            yi[i] = 0.0;
        }

        switch (nd)
        {
            case 1:
            {
                for (i = 0; i < ni; i++)
                {
                    yi[i] = yd[0];
                }

                return yi;
            }
        }

        for (i = 0; i < ni; i++)
        {
            double t;
            if (xi[i] <= xd[0])
            {
                t = (xi[i] - xd[0]) / (xd[1] - xd[0]);
                yi[i] = (1.0 - t) * yd[0] + t * yd[1];
            }
            else if (xd[nd - 1] <= xi[i])
            {
                t = (xi[i] - xd[nd - 2]) / (xd[nd - 1] - xd[nd - 2]);
                yi[i] = (1.0 - t) * yd[nd - 2] + t * yd[nd - 1];
            }
            else
            {
                int k;
                for (k = 1; k < nd; k++)
                {
                    if (xd[k - 1] <= xi[i] && xi[i] <= xd[k])
                    {
                        t = (xi[i] - xd[k - 1]) / (xd[k] - xd[k - 1]);
                        yi[i] = (1.0 - t) * yd[k - 1] + t * yd[k];
                        break;
                    }
                }
            }
        }

        return yi;
    }
}