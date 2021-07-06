using System;

namespace Burkardt.Edge
{
    public class ThreeD
    {
        public static double fxyz1(double x, double y, double z)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FXYZ1 is the first 3D example, scalar version.
            //
            //  Discussion:
            //
            //    This function allows the user a more convenient interface when
            //    only a single input argument is supplied.  See FXYZ1_VEC for details.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 February 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
        {
            double value;
            double[] value_vec;
            double[] x_vec = new double[1];
            double[] y_vec = new double[1];
            double[] z_vec = new double[1];

            x_vec[0] = x;
            y_vec[0] = y;
            z_vec[0] = z;
            value_vec = fxyz1_vec(1, x_vec, y_vec, z_vec);
            value = value_vec[0];

            return value;
        }

        public static double[] fxyz1_vec(int n, double[] x, double[] y, double[] z )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FXYZ1_VEC is the first 3D example, vector version.
        //
        //  Discussion:
        //
        //    This example is known as the 3D Shepp-Logan phantom.
        //
        //    It should be plotted on [-1,+1] x [-1,+1.5] x [-1.5,+1.5].
        //
        //    Seventeen objects are modeled by ellipses of various gray levels,
        //    including:
        //
        //     1: Outer skull
        //     2: Inner skull
        //     3: Left eye
        //     4: Right eye
        //     5: Nose
        //     6: Mouth
        //     7: Left ear
        //     8: Right ear
        //     9: Left small tumor
        //    10: Center small tumor
        //    11: Right small tumor
        //    12: Old f
        //    13: Old g
        //    14: Old e
        //    15: Right ventricle
        //    16: Left ventricle
        //    17: Blood clot
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    14 February 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Larry Shepp,
        //    Computerized tomography and nuclear magnetic resonance,
        //    Journal of Computer Assisted Tomography,
        //    Volume 4, Number 1, February 1980, pages 94-107.
        //
        //  Parameters:
        //
        //    Input, int N, the number of points.
        //
        //    Input, double X[N], Y[N], Z[N], the arguments.
        //
        //    Output, double FXYZ1_VEC[N], the function values.
        //
        {
            double[] a1 =  {
                0.7233, 0.7008, 0.1270, 0.1270, 0.1270,
                0.4575, 0.0635, 0.0635, 0.0460, 0.0230,
                0.0230, 0.0460, 0.2100, 0.1100, 0.1600,
                0.1600, 0.0300
            }
            ;
            double[] a2 =  {
                0.9644, 0.9246, 0.1270, 0.1270, 0.3400,
                0.6099, 0.3175, 0.3175, 0.0230, 0.0230,
                0.0460, 0.0460, 0.2581, 0.2500, 0.3100,
                0.4100, 0.2000
            }
            ;
            double[] a3 =  {
                1.2700, 1.2241, 0.1270, 0.1270, 0.1700,
                0.5080, 0.3175, 0.3175, 0.0230, 0.0460,
                0.0230, 0.0460, 0.2581, 0.2300, 0.2540,
                0.3810, 0.2000
            }
            ;
            double c;
            double[] f;
            double[] g =  {
                2.0000, -0.9800, -1.0000, -1.0000, 1.5000,
                -1.0000, 1.0000, 1.0000, 0.0100, 0.0100,
                0.0100, 0.0100, 0.0100, 0.0100, -0.0200,
                -0.0200, 0.0300
            }
            ;
            int e;
            int i;
            double[] v11 =  {
                1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
                1.0000, 0.9903, -0.9903, 1.0000, 1.0000,
                1.0000, 1.0000, 1.0000, 1.0000, 0.9511,
                -0.9511, 0.9192
            }
            ;
            double[] v12 =  {
                0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                0.0000, -0.1085, -0.1085, 0.0000, 0.0000,
                0.0000, 0.0000, 0.0000, 0.0000, -0.3090,
                -0.3090, -0.3381
            }
            ;
            double[] v13 =  {
                0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                0.0000, -0.0865, -0.0865, 0.0000, 0.0000,
                0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                0.0000, 0.2020
            }
            ;
            double[] v21 =  {
                0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                0.0000, 0.1089, -0.1089, 0.0000, 0.0000,
                0.0000, 0.0000, 0.0000, 0.0000, 0.3090,
                -0.3090, 0.3452
            }
            ;
            double[] v22 =  {
                1.0000, 1.0000, 1.0000, 1.0000, 0.5446,
                1.0000, 0.9941, 0.9941, 1.0000, 1.0000,
                1.0000, 1.0000, 1.0000, 1.0000, 0.9511,
                0.9511, 0.9385
            }
            ;
            double[] v23 =  {
                0.0000, 0.0000, 0.0000, 0.0000, -0.8387,
                0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                0.0000, 0.0000
            }
            ;
            double[] v31 =  {
                0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                0.0000, 0.0860, -0.0860, 0.0000, 0.0000,
                0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                0.0000, 0.1896
            }
            ;
            double[] v32 =  {
                0.0000, 0.0000, 0.0000, 0.0000, 0.8387,
                0.0000, -0.0094, -0.0094, 0.0000, 0.0000,
                0.0000, 0.0000, 0.0000, 0.0000, 0.0000,
                0.0000, -0.0697
            }
            ;
            double[] v33 =  {
                1.0000, 1.0000, 1.0000, 1.0000, 0.5446,
                1.0000, 0.9963, 0.9963, 1.0000, 1.0000,
                1.0000, 1.0000, 1.0000, 1.0000, 1.0000,
                1.0000, -0.9794
            }
            ;
            double[] x0 =  {
                0.0000, 0.0000, 0.2583, -0.2583, 0.0000,
                0.0000, 0.7076, -0.7076, -0.0800, 0.0000,
                0.0600, 0.0000, 0.0000, 0.0000, 0.2200,
                -0.2200, 0.5600
            }
            ;
            double[] y0 =  {
                0.0000, -0.0184, 0.7534, 0.7534, 1.1398,
                0.0000, -0.1378, -0.1378, -0.6050, -0.6050,
                -0.6050, 0.1000, -0.1000, 0.3500, 0.0000,
                0.0000, -0.4000
            }
            ;
            double[] z0 =  {
                0.0000, -0.0185, 0.0000, 0.0000, -0.1957,
                -0.7620, -0.1905, -0.1905, 0.3810, 0.3810,
                0.3810, 0.3810, 0.1270, 0.3810, 0.3810,
                0.3810, 0.3810
            }
            ;

            f = new double[n];

            for (i = 0; i < n; i++)
            {
                f[i] = 0.0;

                for (e = 0; e < 17; e++)
                {
                    c = Math.Pow(((x[i] - x0[e]) * v11[e]
                             + (y[i] - y0[e]) * v12[e]
                             + (z[i] - z0[e]) * v13[e]) / a1[e], 2)
                        + Math.Pow(((x[i] - x0[e]) * v21[e]
                               + (y[i] - y0[e]) * v22[e]
                               + (z[i] - z0[e]) * v23[e]) / a2[e], 2)
                        + Math.Pow(((x[i] - x0[e]) * v31[e]
                               + (y[i] - y0[e]) * v32[e]
                               + (z[i] - z0[e]) * v33[e]) / a3[e], 2);

                    if (c <= 1.0)
                    {
                        f[i] = f[i] + g[e];
                    }

                }

            }

            return f;
        }
    }
}