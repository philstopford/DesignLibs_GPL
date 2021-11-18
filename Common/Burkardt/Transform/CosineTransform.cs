using System;

namespace Burkardt.Transform;

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
        int i;

        double[] c = new double[n];

        for (i = 0; i < n; i++)
        {
            c[i] = 0.0;
            int j;
            for (j = 0; j < n; j++)
            {
                double angle = Math.PI * (i * (2 * j + 1)) / (2 * n);
                c[i] += Math.Cos(angle) * d[j];
            }

            c[i] *= Math.Sqrt(2.0 / n);
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
        int i;

        double[] d = new double[n];

        for (i = 0; i < n; i++)
        {
            d[i] = c[0] / 2.0;
            int j;
            for (j = 1; j < n; j++)
            {
                double angle = Math.PI * ((2 * i + 1) * j) / (2 * n);
                d[i] += Math.Cos(angle) * c[j];
            }

            d[i] *= Math.Sqrt(2.0 / n);
        }

        return d;
    }
}