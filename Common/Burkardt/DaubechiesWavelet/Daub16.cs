using Burkardt.Types;

namespace Burkardt.DaubechiesWavelet;

public static class Daub16
{
    public static double[] daub16_transform(int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DAUB16_TRANSFORM computes the DAUB16 transform of a vector.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 May 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the dimension of the vector.
        //    N must be a power of 2 and at least 4.
        //
        //    Input, double X[N], the vector to be transformed. 
        //
        //    Output, double DAUB16_TRANSFORM[N], the transformed vector.
        //
    {
        double[] c =  {
                5.441584224310400E-02,
                3.128715909142999E-01,
                6.756307362972898E-01,
                5.853546836542067E-01,
                -1.582910525634930E-02,
                -2.840155429615469E-01,
                4.724845739132827E-04,
                1.287474266204784E-01,
                -1.736930100180754E-02,
                -4.408825393079475E-02,
                1.398102791739828E-02,
                8.746094047405776E-03,
                -4.870352993451574E-03,
                -3.917403733769470E-04,
                6.754494064505693E-04,
                -1.174767841247695E-04
            }
            ;
        int i;
        int j;
        int j0;
        int j1;
        int k;
        int m;
        int p = 15;
        double[] y;
        double[] z;

        y = typeMethods.r8vec_copy_new(n, x);
        z = new double[n];

        m = n;

        while (4 <= m)
        {
            for (i = 0; i < m; i++)
            {
                z[i] = 0.0;
            }

            i = 0;

            for (j = 0; j < m - 1; j += 2)
            {
                for (k = 0; k < p; k += 2)
                {
                    j0 = typeMethods.i4_wrap(j + k, 0, m - 1);
                    j1 = typeMethods.i4_wrap(j + k + 1, 0, m - 1);
                    z[i] = z[i] + c[k] * y[j0] + c[k + 1] * y[j1];
                    z[i + m / 2] = z[i + m / 2] + c[p - k] * y[j0] - c[p - k - 1] * y[j1];
                }

                i += 1;
            }

            for (i = 0; i < m; i++)
            {
                y[i] = z[i];
            }

            m /= 2;
        }
            
        return y;
    }

    public static double[] daub16_transform_inverse(int n, double[] y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DAUB16_TRANSFORM_INVERSE inverts the DAUB16 transform of a vector.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 May 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the dimension of the vector.
        //    N must be a power of 2 and at least 4.
        //
        //    Input, double Y[N], the transformed vector. 
        //
        //    Output, double DAUB16_TRANSFORM_INVERSE[N], the original vector.
        //
    {
        double[] c =  {
                5.441584224310400E-02,
                3.128715909142999E-01,
                6.756307362972898E-01,
                5.853546836542067E-01,
                -1.582910525634930E-02,
                -2.840155429615469E-01,
                4.724845739132827E-04,
                1.287474266204784E-01,
                -1.736930100180754E-02,
                -4.408825393079475E-02,
                1.398102791739828E-02,
                8.746094047405776E-03,
                -4.870352993451574E-03,
                -3.917403733769470E-04,
                6.754494064505693E-04,
                -1.174767841247695E-04
            }
            ;
        int i;
        int i0;
        int i1;
        int j;
        int k;
        int m;
        int p = 15;
        int q;
        double[] x;
        double[] z;

        x = typeMethods.r8vec_copy_new(n, y);
        z = new double[n];

        m = 4;
        q = (p - 1) / 2;

        while (m <= n)
        {
            for (i = 0; i < n; i++)
            {
                z[i] = 0.0;
            }

            j = 0;

            for (i = -q; i < m / 2 - q; i++)
            {
                for (k = 0; k < p; k += 2)
                {
                    i0 = typeMethods.i4_wrap(i + k / 2, 0, m / 2 - 1);
                    i1 = typeMethods.i4_wrap(i + m / 2 + k / 2, m / 2, m - 1);
                    z[j] = z[j] + c[p - k - 1] * x[i0] + c[k + 1] * x[i1];
                    z[j + 1] = z[j + 1] + c[p - k] * x[i0] - c[k] * x[i1];
                }

                j += 2;
            }

            for (i = 0; i < m; i++)
            {
                x[i] = z[i];
            }

            m *= 2;
        }

        return x;
    }
}