using System;
using Burkardt.Types;

namespace Burkardt.Interpolation
{
    public static class Shepard
    {
        public static double[] shepard_basis_1d(int nd, double[] xd, double p, int ni, double[] xi)

//****************************************************************************80
//
//  Purpose:
//
//    SHEPARD_BASIS_1D evaluates a 1D Shepard basis function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 July 2015
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Donald Shepard,
//    A two-dimensional interpolation function for irregularly spaced data,
//    ACM '68: Proceedings of the 1968 23rd ACM National Conference,
//    ACM, pages 517-524, 1969.
//
//  Parameters:
//
//    Input, int ND, the number of data points.
//
//    Input, double XD[ND], the data points.
//
//    Input, double P, the power.
//
//    Input, int NI, the number of interpolation points.
//
//    Input, double XI[NI], the interpolation points.
//
//    Output, double SHEPARD_BASIS_1D[NI*ND], the basis functions at the interpolation 
//    points.
// 
        {
            double[] w = new double[nd];
            double[] bk = new double[ni * nd];

            for (int i = 0; i < ni; i++)
            {
                if (p == 0.0)
                {
                    for (int j = 0; j < nd; j++)
                    {
                        w[j] = 1.0 / (double) (nd);
                    }
                }
                else
                {
                    int z = -1;
                    for (int j = 0; j < nd; j++)
                    {
                        w[j] = Math.Abs(xi[i] - xd[j]);
                        if (w[j] == 0.0)
                        {
                            z = j;
                            break;
                        }
                    }

                    if (z != -1)
                    {
                        for (int j = 0; j < nd; j++)
                        {
                            w[j] = 0.0;
                        }

                        w[z] = 1.0;
                    }
                    else
                    {
                        for (int j = 0; j < nd; j++)
                        {
                            w[j] = 1.0 / Math.Pow(w[j], p);
                        }

                        double s = typeMethods.r8vec_sum(nd, w);
                        for (int j = 0; j < nd; j++)
                        {
                            w[j] = w[j] / s;
                        }
                    }
                }

                for (int j = 0; j < nd; j++)
                {
                    bk[i + j * ni] = w[j];
                }
            }

            return bk;
        }

        public static double[] shepard_value_1d(int nd, double[] xd, double[] yd, double p, int ni,
            double[] xi)
//****************************************************************************80
//
//  Purpose:
//
//    SHEPARD_VALUE_1D evaluates a 1D Shepard interpolant.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    01 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Donald Shepard,
//    A two-dimensional interpolation function for irregularly spaced data,
//    ACM '68: Proceedings of the 1968 23rd ACM National Conference,
//    ACM, pages 517-524, 1969.
//
//  Parameters:
//
//    Input, int ND, the number of data points.
//
//    Input, double XD[ND], the data points.
//
//    Input, double YD[ND], the data values.
//
//    Input, double P, the power.
//
//    Input, int NI, the number of interpolation points.
//
//    Input, double XI[NI], the interpolation points.
//
//    Output, double SHEPARD_VALUE_1D[NI], the interpolated values.
//
        {
            double[] w = new double[nd];
            double[] yi = new double[ni];

            for (int i = 0; i < ni; i++)
            {
                if (p == 0.0)
                {
                    for (int j = 0; j < nd; j++)
                    {
                        w[j] = 1.0 / (double) (nd);
                    }
                }
                else
                {
                    int z = -1;
                    for (int j = 0; j < nd; j++)
                    {
                        w[j] = Math.Abs(xi[i] - xd[j]);
                        if (w[j] == 0.0)
                        {
                            z = j;
                            break;
                        }
                    }

                    if (z != -1)
                    {
                        for (int j = 0; j < nd; j++)
                        {
                            w[j] = 0.0;
                        }

                        w[z] = 1.0;
                    }
                    else
                    {
                        for (int j = 0; j < nd; j++)
                        {
                            w[j] = 1.0 / Math.Pow(w[j], p);
                        }

                        double s = typeMethods.r8vec_sum(nd, w);
                        for (int j = 0; j < nd; j++)
                        {
                            w[j] = w[j] / s;
                        }
                    }
                }

                yi[i] = typeMethods.r8vec_dot_product(nd, w, yd);
            }

            return yi;
        }

        public static double[] shepard_interp_2d(int nd, double[] xd, double[] yd, double[] zd,
        double p, int ni, double[] xi, double[] yi )
//****************************************************************************80
//
//  Purpose:
//
//    SHEPARD_INTERP_2D evaluates a 2D Shepard interpolant.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Donald Shepard,
//    A two-dimensional interpolation function for irregularly spaced data,
//    ACM '68: Proceedings of the 1968 23rd ACM National Conference,
//    ACM, pages 517-524, 1969.
//
//  Parameters:
//
//    Input, int ND, the number of data points.
//
//    Input, double XD[ND], YD[ND], the data points.
//
//    Input, double ZD[ND], the data values.
//
//    Input, double P, the power.
//
//    Input, int NI, the number of interpolation points.
//
//    Input, double XI[NI], YI[NI], the interpolation points.
//
//    Output, double SHEPARD_INTERP_2D[NI], the interpolated values.
//
        {
            double[] w = new double[nd];
            double[] zi = new double[ni];

            for (int i = 0; i < ni; i++)
            {
                if (p == 0.0)
                {
                    for (int j = 0; j < nd; j++)
                    {
                        w[j] = 1.0 / (double) (nd);
                    }
                }
                else
                {
                    int z = -1;
                    for (int j = 0; j < nd; j++)
                    {
                        w[j] = Math.Sqrt(Math.Pow(xi[i] - xd[j], 2)
                                    + Math.Pow(yi[i] - yd[j], 2));
                        if (w[j] == 0.0)
                        {
                            z = j;
                            break;
                        }
                    }

                    if (z != -1)
                    {
                        for (int j = 0; j < nd; j++)
                        {
                            w[j] = 0.0;
                        }

                        w[z] = 1.0;
                    }
                    else
                    {
                        for (int j = 0; j < nd; j++)
                        {
                            w[j] = 1.0 / Math.Pow(w[j], p);
                        }

                        double s = 0.0;
                        for (int j = 0; j < nd; j++)
                        {
                            s = s + w[j];
                        }

                        for (int j = 0; j < nd; j++)
                        {
                            w[j] = w[j] / s;
                        }
                    }
                }

                zi[i] = typeMethods.r8vec_dot_product(nd, w, zd);
            }

            return zi;
        }

        public static double[] shepard_interp_nd(int m, int nd, double[] xd, double[] zd, double p,
        int ni, double[] xi )
//****************************************************************************80
//
//  Purpose:
//
//    SHEPARD_INTERP_ND evaluates a multidimensional Shepard interpolant.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    02 October 2012
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Donald Shepard,
//    A two-dimensional interpolation function for irregularly spaced data,
//    ACM '68: Proceedings of the 1968 23rd ACM National Conference,
//    ACM, pages 517-524, 1969.
//
//  Parameters:
//
//    Input, int M, the spatial dimension.
//
//    Input, int ND, the number of data points.
//
//    Input, double XD[M*ND], the data points.
//
//    Input, double ZD[ND], the data values.
//
//    Input, double P, the power.
//
//    Input, int NI, the number of interpolation points.
//
//    Input, double XI[M*NI], the interpolation points.
//
//    Output, double SHEPARD_INTERP_ND[NI], the interpolated values.
//
        {
            double[] w = new double[nd];
            double[] zi = new double[ni];

            for (int i = 0; i < ni; i++)
            {
                if (p == 0.0)
                {
                    for (int j = 0; j < nd; j++)
                    {
                        w[j] = 1.0 / (double) (nd);
                    }
                }
                else
                {
                    int z = -1;
                    for (int j = 0; j < nd; j++)
                    {
                        double t = 0.0;
                        int i2;
                        for (i2 = 0; i2 < m; i2++)
                        {
                            t = t + Math.Pow(xi[i2 + i * m] - xd[i2 + j * m], 2);
                        }

                        w[j] = Math.Sqrt(t);
                        if (w[j] == 0.0)
                        {
                            z = j;
                            break;
                        }
                    }

                    if (z != -1)
                    {
                        for (int j = 0; j < nd; j++)
                        {
                            w[j] = 0.0;
                        }

                        w[z] = 1.0;
                    }
                    else
                    {
                        for (int j = 0; j < nd; j++)
                        {
                            w[j] = 1.0 / Math.Pow(w[j], p);
                        }

                        double s = 0.0;
                        for (int j = 0; j < nd; j++)
                        {
                            s = s + w[j];
                        }

                        for (int j = 0; j < nd; j++)
                        {
                            w[j] = w[j] / s;
                        }
                    }
                }

                zi[i] = typeMethods.r8vec_dot_product(nd, w, zd);
            }

            return zi;
        }

    }
}