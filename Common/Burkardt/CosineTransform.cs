using System;

namespace Burkardt.Transform
{
    public static class Cosine
    {
        public static double[] cosine_transform_data(int n, double[] d)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    COSINE_TRANSFORM_DATA does a cosine transform on a vector of data.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    27 August 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, integer N, the number of data points.
            //
            //    Input, double D[N], the vector of data.
            //
            //    Output, double COSINE_TRANSFORM_DATA[N], the transform coefficients.
            //
        {
            double angle;
            double[] c;
            int i;
            int j;

            c = new double[n];

            for (i = 0; i < n; i++)
            {
                c[i] = 0.0;
                for (j = 0; j < n; j++)
                {
                    angle = Math.PI * (double) (i * (2 * j + 1)) / (double) (2 * n);
                    c[i] = c[i] + Math.Cos(angle) * d[j];
                }

                c[i] = c[i] * Math.Sqrt(2.0 / (double) (n));
            }

            return c;
        }

        public static double[] cosine_transform_inverse(int n, double[] c)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    COSINE_TRANSFORM_INVERSE does an inverse cosine transform.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    27 August 2015
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, integer N, the number of data points.
            //
            //    Input, double C[N], the vector of transform coefficients
            //
            //    Output, double COSINE_TRANSFORM_INVERSE[N], the original data.
            //
        {
            double angle;
            double[] d;
            int i;
            int j;
            double r8_pi = 3.141592653589793;

            d = new double[n];

            for (i = 0; i < n; i++)
            {
                d[i] = c[0] / 2.0;
                for (j = 1; j < n; j++)
                {
                    angle = r8_pi * (double) ((2 * i + 1) * j) / (double) (2 * n);
                    d[i] = d[i] + Math.Cos(angle) * c[j];
                }

                d[i] = d[i] * Math.Sqrt(2.0 / (double) (n));
            }

            return d;
        }
    }
}