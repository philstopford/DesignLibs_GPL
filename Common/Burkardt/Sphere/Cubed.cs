using System;

namespace Burkardt.SphereNS;

public static class Cubed
{
    public static double[] sphere_cubed_ijk_to_xyz_old(int n, int i, int j, int k)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_CUBED_IJK_TO_XYZ_OLD: cubed sphere IJK to XYZ coordinates.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of sections into which each 
        //    face of the cube is to be divided.
        //
        //    Input, int I, J, K, indices between 0 and N.  Normally,
        //    at least one of the indices should have the value 0 or N.
        //
        //    Output, double SPHERE_CUBED_IJK_TO_XYZ_OLD[3], coordinates of the point.
        //
    {
            
        double xc;
        double[] xyz;
        double xyzn;
        double yc;
        double zc;

        xyz = new double[3];

        switch (i)
        {
            case 0:
                xc = -1.0;
                break;
            default:
            {
                if (i == n)
                {
                    xc = +1.0;
                }
                else
                {
                    xc = Math.Tan((2 * i - n) * 0.25 * Math.PI / n);
                }

                break;
            }
        }

        switch (j)
        {
            case 0:
                yc = -1.0;
                break;
            default:
            {
                if (j == n)
                {
                    yc = +1.0;
                }
                else
                {
                    yc = Math.Tan((2 * j - n) * 0.25 * Math.PI / n);
                }

                break;
            }
        }

        switch (k)
        {
            case 0:
                zc = -1.0;
                break;
            default:
            {
                if (k == n)
                {
                    zc = +1.0;
                }
                else
                {
                    zc = Math.Tan((2 * k - n) * 0.25 * Math.PI / n);
                }

                break;
            }
        }

        xyzn = Math.Sqrt(xc * xc + yc * yc + zc * zc);

        xyz[0] = xc / xyzn;
        xyz[1] = yc / xyzn;
        xyz[2] = zc / xyzn;

        return xyz;
    }

    public static void sphere_cubed_ijk_to_xyz(int n, int i, int j, int k, ref double[] xyz, int xyzIndex = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_CUBED_IJK_TO_XYZ: cubed sphere IJK to XYZ coordinates.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of sections into which each 
        //    face of the cube is to be divided.
        //
        //    Input, int I, J, K, indices between 0 and N.  Normally,
        //    at least one of the indices should have the value 0 or N.
        //
        //    Output, double XYZ[3], coordinates of the point.
        //
    {
            
        double xc;
        double xyzn;
        double yc;
        double zc;

        switch (i)
        {
            case 0:
                xc = -1.0;
                break;
            default:
            {
                if (i == n)
                {
                    xc = +1.0;
                }
                else
                {
                    xc = Math.Tan((2 * i - n) * 0.25 * Math.PI / n);
                }

                break;
            }
        }

        switch (j)
        {
            case 0:
                yc = -1.0;
                break;
            default:
            {
                if (j == n)
                {
                    yc = +1.0;
                }
                else
                {
                    yc = Math.Tan((2 * j - n) * 0.25 * Math.PI / n);
                }

                break;
            }
        }

        switch (k)
        {
            case 0:
                zc = -1.0;
                break;
            default:
            {
                if (k == n)
                {
                    zc = +1.0;
                }
                else
                {
                    zc = Math.Tan((2 * k - n) * 0.25 * Math.PI / n);
                }

                break;
            }
        }

        xyzn = Math.Sqrt(xc * xc + yc * yc + zc * zc);

        xyz[xyzIndex + 0] = xc / xyzn;
        xyz[xyzIndex + 1] = yc / xyzn;
        xyz[xyzIndex + 2] = zc / xyzn;
    }

    public static int sphere_cubed_line_num(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_CUBED_LINE_NUM counts lines on a cubed sphere grid.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of sections into which each 
        //    face of the cube is to be divided.
        //
        //    Output, int LINE_NUM, the number of lines.
        //
    {
        int line_num;

        line_num = 0;
        switch (n)
        {
            //
            //  If N = 1, the corners form 12 lines.
            //
            case 1:
                line_num = 12;
                return line_num;
        }
        //
        //  If 1 < N, each of 8 corners connects to three neighboring edges.
        //

        line_num += 8 * 3;

        switch (n)
        {
            //
            //  If 2 < N, then each of the 12 edges includes lines.
            //
            case > 2:
                line_num += 12 * (n - 2);
                break;
        }

        switch (n)
        {
            //
            //  Lines that belong to one of the six faces.
            //
            case > 1:
                line_num += 6 * 2 * n * (n - 1);
                break;
        }

        return line_num;
    }

    public static double[] sphere_cubed_lines(int n, int line_num)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_CUBED_LINES computes the lines on a cubed sphere grid.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    15 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of sections into which 
        //    each face of the cube is to be divided.
        //
        //    Input, int LINE_NUM, the number of lines.
        //
        //    Output, double SPHERE_CUBED_LINES[3*2*LINE_NUM], distinct points 
        //    on the unit sphere generated by a cubed sphere grid.
        //
    {
        int i;
        int j;
        int l;
        double[] line_data;

        line_data = new double[3 * 2 * line_num];
        l = 0;
        switch (n)
        {
            //
            //  If N = 1, the corners form 12 lines.
            //
            case 1:
                sphere_cubed_ijk_to_xyz(n, 0, 0, 0, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
                sphere_cubed_ijk_to_xyz(n, n, 0, 0, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
                l += 1;
                sphere_cubed_ijk_to_xyz(n, n, 0, 0, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
                sphere_cubed_ijk_to_xyz(n, n, n, 0, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
                l += 1;
                sphere_cubed_ijk_to_xyz(n, n, n, 0, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
                sphere_cubed_ijk_to_xyz(n, 0, n, 0, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
                l += 1;
                sphere_cubed_ijk_to_xyz(n, 0, n, 0, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
                sphere_cubed_ijk_to_xyz(n, 0, 0, 0, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);

                l += 1;
                sphere_cubed_ijk_to_xyz(n, 0, 0, n, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
                sphere_cubed_ijk_to_xyz(n, n, 0, n, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
                l += 1;
                sphere_cubed_ijk_to_xyz(n, n, 0, n, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
                sphere_cubed_ijk_to_xyz(n, n, n, n, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
                l += 1;
                sphere_cubed_ijk_to_xyz(n, n, n, n, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
                sphere_cubed_ijk_to_xyz(n, 0, n, n, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
                l += 1;
                sphere_cubed_ijk_to_xyz(n, 0, n, n, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
                sphere_cubed_ijk_to_xyz(n, 0, 0, n, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);

                l += 1;
                sphere_cubed_ijk_to_xyz(n, 0, 0, 0, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
                sphere_cubed_ijk_to_xyz(n, 0, 0, n, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
                l += 1;
                sphere_cubed_ijk_to_xyz(n, n, 0, 0, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
                sphere_cubed_ijk_to_xyz(n, n, 0, n, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
                l += 1;
                sphere_cubed_ijk_to_xyz(n, n, n, 0, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
                sphere_cubed_ijk_to_xyz(n, n, n, n, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
                l += 1;
                sphere_cubed_ijk_to_xyz(n, 0, n, 0, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
                sphere_cubed_ijk_to_xyz(n, 0, n, n, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
                l += 1;
                return line_data;
        }
        //
        //  If 1 < N, each of 8 corners connects to three neighboring edges.
        //

        sphere_cubed_ijk_to_xyz(n, 0, 0, 0, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
        sphere_cubed_ijk_to_xyz(n, 1, 0, 0, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
        l += 1;
        sphere_cubed_ijk_to_xyz(n, 0, 0, 0, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
        sphere_cubed_ijk_to_xyz(n, 0, 1, 0, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
        l += 1;
        sphere_cubed_ijk_to_xyz(n, 0, 0, 0, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
        sphere_cubed_ijk_to_xyz(n, 0, 0, 1, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);

        l += 1;
        sphere_cubed_ijk_to_xyz(n, n, 0, 0, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
        sphere_cubed_ijk_to_xyz(n, n - 1, 0, 0, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
        l += 1;
        sphere_cubed_ijk_to_xyz(n, n, 0, 0, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
        sphere_cubed_ijk_to_xyz(n, 0, 1, 0, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
        l += 1;
        sphere_cubed_ijk_to_xyz(n, n, 0, 0, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
        sphere_cubed_ijk_to_xyz(n, n, 0, 1, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);

        l += 1;
        sphere_cubed_ijk_to_xyz(n, n, n, 0, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
        sphere_cubed_ijk_to_xyz(n, n - 1, n, 0, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
        l += 1;
        sphere_cubed_ijk_to_xyz(n, n, n, 0, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
        sphere_cubed_ijk_to_xyz(n, n, n - 1, 0, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
        l += 1;
        sphere_cubed_ijk_to_xyz(n, n, n, 0, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
        sphere_cubed_ijk_to_xyz(n, n, n, 1, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);

        l += 1;
        sphere_cubed_ijk_to_xyz(n, 0, n, 0, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
        sphere_cubed_ijk_to_xyz(n, 1, n, 0, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
        l += 1;
        sphere_cubed_ijk_to_xyz(n, 0, n, 0, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
        sphere_cubed_ijk_to_xyz(n, 0, n - 1, 0, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
        l += 1;
        sphere_cubed_ijk_to_xyz(n, 0, n, 0, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
        sphere_cubed_ijk_to_xyz(n, 0, n, 1, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);

        l += 1;
        sphere_cubed_ijk_to_xyz(n, 0, 0, n, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
        sphere_cubed_ijk_to_xyz(n, 1, 0, n, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
        l += 1;
        sphere_cubed_ijk_to_xyz(n, 0, 0, n, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
        sphere_cubed_ijk_to_xyz(n, 0, 1, n, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
        l += 1;
        sphere_cubed_ijk_to_xyz(n, 0, 0, n, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
        sphere_cubed_ijk_to_xyz(n, 0, 0, n - 1, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);

        l += 1;
        sphere_cubed_ijk_to_xyz(n, n, 0, n, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
        sphere_cubed_ijk_to_xyz(n, n - 1, 0, n, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
        l += 1;
        sphere_cubed_ijk_to_xyz(n, n, 0, n, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
        sphere_cubed_ijk_to_xyz(n, n, 1, n, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
        l += 1;
        sphere_cubed_ijk_to_xyz(n, n, 0, n, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
        sphere_cubed_ijk_to_xyz(n, n, 0, n - 1, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);

        l += 1;
        sphere_cubed_ijk_to_xyz(n, n, n, n, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
        sphere_cubed_ijk_to_xyz(n, n - 1, n, n, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
        l += 1;
        sphere_cubed_ijk_to_xyz(n, n, n, n, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
        sphere_cubed_ijk_to_xyz(n, n, n - 1, n, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
        l += 1;
        sphere_cubed_ijk_to_xyz(n, n, n, n, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
        sphere_cubed_ijk_to_xyz(n, n, n, n - 1, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);

        l += 1;
        sphere_cubed_ijk_to_xyz(n, 0, n, n, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
        sphere_cubed_ijk_to_xyz(n, 1, n, n, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
        l += 1;
        sphere_cubed_ijk_to_xyz(n, 0, n, n, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
        sphere_cubed_ijk_to_xyz(n, 0, n - 1, n, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
        l += 1;
        sphere_cubed_ijk_to_xyz(n, 0, n, n, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
        sphere_cubed_ijk_to_xyz(n, 0, n, n - 1, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
        l += 1;

        switch (n)
        {
            //
            //  If 2 < N, then each of the 12 edges includes lines.
            //
            case > 2:
            {
                for (i = 1; i <= n - 2; i++)
                {
                    sphere_cubed_ijk_to_xyz(n, i, 0, 0, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
                    sphere_cubed_ijk_to_xyz(n, i + 1, 0, 0, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
                    l += 1;
                }

                for (i = 1; i <= n - 2; i++)
                {
                    sphere_cubed_ijk_to_xyz(n, n, i, 0, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
                    sphere_cubed_ijk_to_xyz(n, n, i + 1, 0, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
                    l += 1;
                }

                for (i = 1; i <= n - 2; i++)
                {
                    sphere_cubed_ijk_to_xyz(n, n - i, n, 0, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
                    sphere_cubed_ijk_to_xyz(n, n - i - 1, n, 0, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
                    l += 1;
                }

                for (i = 1; i <= n - 2; i++)
                {
                    sphere_cubed_ijk_to_xyz(n, 0, n - i, 0, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
                    sphere_cubed_ijk_to_xyz(n, 0, n - i - 1, 0, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
                    l += 1;
                }

                for (i = 1; i <= n - 2; i++)
                {
                    sphere_cubed_ijk_to_xyz(n, i, 0, n, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
                    sphere_cubed_ijk_to_xyz(n, i + 1, 0, n, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
                    l += 1;
                }

                for (i = 1; i <= n - 2; i++)
                {
                    sphere_cubed_ijk_to_xyz(n, n, i, n, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
                    sphere_cubed_ijk_to_xyz(n, n, i + 1, n, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
                    l += 1;
                }

                for (i = 1; i <= n - 2; i++)
                {
                    sphere_cubed_ijk_to_xyz(n, n - i, n, n, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
                    sphere_cubed_ijk_to_xyz(n, n - i - 1, n, n, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
                    l += 1;
                }

                for (i = 1; i <= n - 2; i++)
                {
                    sphere_cubed_ijk_to_xyz(n, 0, n - i, n, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
                    sphere_cubed_ijk_to_xyz(n, 0, n - i - 1, n, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
                    l += 1;
                }

                for (i = 1; i <= n - 2; i++)
                {
                    sphere_cubed_ijk_to_xyz(n, 0, 0, i, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
                    sphere_cubed_ijk_to_xyz(n, 0, 0, i + 1, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
                    l += 1;
                }

                for (i = 1; i <= n - 2; i++)
                {
                    sphere_cubed_ijk_to_xyz(n, n, 0, i, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
                    sphere_cubed_ijk_to_xyz(n, n, 0, i + 1, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
                    l += 1;
                }

                for (i = 1; i <= n - 2; i++)
                {
                    sphere_cubed_ijk_to_xyz(n, n, n, i, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
                    sphere_cubed_ijk_to_xyz(n, n, n, i + 1, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
                    l += 1;
                }

                for (i = 1; i <= n - 2; i++)
                {
                    sphere_cubed_ijk_to_xyz(n, 0, n, i, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
                    sphere_cubed_ijk_to_xyz(n, 0, n, i + 1, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
                    l += 1;
                }

                break;
            }
        }

        switch (n)
        {
            //
            //  Lines that belong to one of the six faces.
            //
            case > 1:
            {
                //
                //  000, nn0
                //
                for (i = 1; i <= n - 1; i++)
                {
                    for (j = 0; j <= n - 1; j++)
                    {
                        sphere_cubed_ijk_to_xyz(n, i, j, 0, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
                        sphere_cubed_ijk_to_xyz(n, i, j + 1, 0, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
                        l += 1;
                    }
                }

                for (j = 1; j <= n - 1; j++)
                {
                    for (i = 0; i <= n - 1; i++)
                    {
                        sphere_cubed_ijk_to_xyz(n, i, j, 0, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
                        sphere_cubed_ijk_to_xyz(n, i + 1, j, 0, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
                        l += 1;
                    }
                }

                //
                //  00n, nnn
                //
                for (i = 1; i <= n - 1; i++)
                {
                    for (j = 0; j <= n - 1; j++)
                    {
                        sphere_cubed_ijk_to_xyz(n, i, j, n, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
                        sphere_cubed_ijk_to_xyz(n, i, j + 1, n, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
                        l += 1;
                    }
                }

                for (j = 1; j <= n - 1; j++)
                {
                    for (i = 0; i <= n - 1; i++)
                    {
                        sphere_cubed_ijk_to_xyz(n, i, j, n, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
                        sphere_cubed_ijk_to_xyz(n, i + 1, j, n, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
                        l += 1;
                    }
                }

                //
                //  000:n0n
                //
                for (i = 1; i <= n - 1; i++)
                {
                    for (j = 0; j <= n - 1; j++)
                    {
                        sphere_cubed_ijk_to_xyz(n, i, 0, j, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
                        sphere_cubed_ijk_to_xyz(n, i, 0, j + 1, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
                        l += 1;
                    }
                }

                for (j = 1; j <= n - 1; j++)
                {
                    for (i = 0; i <= n - 1; i++)
                    {
                        sphere_cubed_ijk_to_xyz(n, i, 0, j, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
                        sphere_cubed_ijk_to_xyz(n, i + 1, 0, j, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
                        l += 1;
                    }
                }

                //
                //  0n0:nnn
                //
                for (i = 1; i <= n - 1; i++)
                {
                    for (j = 0; j <= n - 1; j++)
                    {
                        sphere_cubed_ijk_to_xyz(n, i, n, j, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
                        sphere_cubed_ijk_to_xyz(n, i, n, j + 1, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
                        l += 1;
                    }
                }

                for (j = 1; j <= n - 1; j++)
                {
                    for (i = 0; i <= n - 1; i++)
                    {
                        sphere_cubed_ijk_to_xyz(n, i, n, j, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
                        sphere_cubed_ijk_to_xyz(n, i + 1, n, j, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
                        l += 1;
                    }
                }

                //
                //  000:0nn
                //
                for (i = 1; i <= n - 1; i++)
                {
                    for (j = 0; j <= n - 1; j++)
                    {
                        sphere_cubed_ijk_to_xyz(n, 0, i, j, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
                        sphere_cubed_ijk_to_xyz(n, 0, i, j + 1, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
                        l += 1;
                    }
                }

                for (j = 1; j <= n - 1; j++)
                {
                    for (i = 0; i <= n - 1; i++)
                    {
                        sphere_cubed_ijk_to_xyz(n, 0, i, j, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
                        sphere_cubed_ijk_to_xyz(n, 0, i + 1, j, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
                        l += 1;
                    }
                }

                //
                //  n00:nnn
                //
                for (i = 1; i <= n - 1; i++)
                {
                    for (j = 0; j <= n - 1; j++)
                    {
                        sphere_cubed_ijk_to_xyz(n, n, i, j, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
                        sphere_cubed_ijk_to_xyz(n, n, i, j + 1, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
                        l += 1;
                    }
                }

                for (j = 1; j <= n - 1; j++)
                {
                    for (i = 0; i <= n - 1; i++)
                    {
                        sphere_cubed_ijk_to_xyz(n, n, i, j, ref line_data, xyzIndex: +0 + 0 * 3 + l * 6);
                        sphere_cubed_ijk_to_xyz(n, n, i + 1, j, ref line_data, xyzIndex: +0 + 1 * 3 + l * 6);
                        l += 1;
                    }
                }

                break;
            }
        }

        if (l != line_num)
        {
            Console.WriteLine("");
            Console.WriteLine("SPHERE_CUBED_LINES - Fatal error!");
            Console.WriteLine("  LINE_NUM = " + line_num + "");
            Console.WriteLine("  L = " + l + "n");
            return null;
        }

        return line_data;
    }

    public static double[] sphere_cubed_points(int n, int ns)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_CUBED_POINTS computes the points on a cubed sphere grid.
        //
        //  Discussion:
        //
        //    For a value of N = 3, for instance, each of the 6 cube faces will
        //    be divided into 3 sections, so that a single cube face will have
        //    (3+1)x(3+1) points:
        //
        //      X---X---X---X
        //      | 1 | 4 | 7 |
        //      X---X---X---X
        //      | 2 | 5 | 8 |
        //      X---X---X---X
        //      | 3 | 6 | 9 |
        //      X---X---X---X
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of sections into which each 
        //    face of the cube is to be divided.
        //
        //    Input, int NS, the number of points.
        //
        //    Output, double SPHERE_CUBED_POINTS[3*NS], distinct points on the 
        //    unit sphere generated by a cubed sphere grid.
        //
    {
        int ns2;
        double[] xyz;

        xyz = new double[3 * ns];

        ns2 = 0;
        //
        //  Bottom full.
        //
        sphere_cubed_points_face(n, 0, 0, 0, n, n, 0, ref ns2, ref xyz);
        //
        //  To avoid repetition, draw the middles as grids of n-2 x n-1 points.
        //
        sphere_cubed_points_face(n, 0, 0, 1, 0, n - 1, n - 1, ref ns2, ref xyz);
        sphere_cubed_points_face(n, 0, n, 1, n - 1, n, n - 1, ref ns2, ref xyz);
        sphere_cubed_points_face(n, n, 1, 1, n, n, n - 1, ref ns2, ref xyz);
        sphere_cubed_points_face(n, 1, 0, 1, n, 0, n - 1, ref ns2, ref xyz);
        //
        //  Top full.
        //
        sphere_cubed_points_face(n, 0, 0, n, n, n, n, ref ns2, ref xyz);

        if (ns2 != ns)
        {
            Console.WriteLine("");
            Console.WriteLine("SPHERE_CUBED_POINTS - Fatal error");
            Console.WriteLine("  Expected to generated NS = " + ns + " points.");
            Console.WriteLine("  Generated " + ns2 + " points.");
            return null;
        }

        return xyz;
    }

    public static void sphere_cubed_points_face(int n, int i1, int j1, int k1, int i2, int j2,
            int k2, ref int ns, ref double[] xyz)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_CUBED_POINTS_FACE: points on one face of a cubed sphere grid.
        //
        //  Discussion:
        //
        //    This routine starts with NS = 0, and is called repeatedly to
        //    add points for another face.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of sections into which each face 
        //    of the cube is to be divided.
        //
        //    Input, int I1, J1, K1, I2, J2, K2, the logical indices, 
        //    between 0 and N, of two corners of the face grid.  It is guaranteed that 
        //    I1 <= I2, J1 <= J2, and K1 <= K2.  
        //
        //    Input/output, int &NS, the number of points.
        //
        //    Input/output, double XYZ[3*NS], distinct points on the unit sphere
        //    generated by a cubed sphere grid.
        //
    {
        int i;
        int j;
        int k;
            
        double xyzn;
        double xc;
        double yc;
        double zc;

        for (i = i1; i <= i2; i++)
        {
            if (i1 < i2)
            {
                xc = Math.Tan((2 * i - n) * 0.25 * Math.PI / n);
            }
            else
            {
                switch (i1)
                {
                    case 0:
                        xc = -1.0;
                        break;
                    default:
                    {
                        if (i1 == n)
                        {
                            xc = +1.0;
                        }
                        else
                        {
                            xc = 0.0;
                        }

                        break;
                    }
                }
            }

            for (j = j1; j <= j2; j++)
            {
                if (j1 < j2)
                {
                    yc = Math.Tan((2 * j - n) * 0.25 * Math.PI / n);
                }
                else
                {
                    switch (j1)
                    {
                        case 0:
                            yc = -1.0;
                            break;
                        default:
                        {
                            if (j1 == n)
                            {
                                yc = +1.0;
                            }
                            else
                            {
                                yc = 0.0;
                            }

                            break;
                        }
                    }
                }

                for (k = k1; k <= k2; k++)
                {
                    if (k1 < k2)
                    {
                        zc = Math.Tan((2 * k - n) * 0.25 * Math.PI / n);
                    }
                    else
                    {
                        switch (k1)
                        {
                            case 0:
                                zc = -1.0;
                                break;
                            default:
                            {
                                if (k1 == n)
                                {
                                    zc = +1.0;
                                }
                                else
                                {
                                    zc = 0.0;
                                }

                                break;
                            }
                        }
                    }

                    xyzn = Math.Sqrt(xc * xc + yc * yc + zc * zc);

                    xyz[0 + ns * 3] = xc / xyzn;
                    xyz[1 + ns * 3] = yc / xyzn;
                    xyz[2 + ns * 3] = zc / xyzn;
                    ns += 1;
                }
            }
        }
    }

    public static int sphere_cubed_point_num(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SPHERE_CUBED_POINT_NUM counts the points on a cubed sphere grid.
        //
        //  Discussion:
        //
        //    For a value of N = 3, for instance, each of the 6 cube faces will
        //    be divided into 3 sections, so that a single cube face will have
        //    (3+1)x(3+1) points:
        //
        //      X---X---X---X
        //      | 1 | 4 | 7 |
        //      X---X---X---X
        //      | 2 | 5 | 8 |
        //      X---X---X---X
        //      | 3 | 6 | 9 |
        //      X---X---X---X
        //
        //    The number of points is simply (N+1)^3 - (N-1)^3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 October 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of sections into which 
        //    each face of the cube is to be divided.
        //
        //    Output, int SPHERE_CUBED_POINT_NUM, the number of points.
        //
    {
        int ns;

        ns = (int)(Math.Pow(n + 1, 3) - Math.Pow(n - 1, 3));

        return ns;
    }
        
}