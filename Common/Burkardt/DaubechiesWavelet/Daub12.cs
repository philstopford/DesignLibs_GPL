using System;
using Burkardt.Types;

namespace Burkardt.DaubechiesWavelet;

public static class Daub12
{
    public static double[] daub12_matrix(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DAUB12_MATRIX returns the DAUB12 matrix.
        //
        //  Discussion:
        //
        //    The DAUB12 matrix is the Daubechies wavelet transformation matrix 
        //    with 12 coefficients.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 May 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //    N must be at least 12 and a multiple of 2.
        //
        //    Output, double DAUB12_MATRIX[N*N], the matrix.
        //
    {
        int i;

        if (n < 12 || n % 2 != 0)
        {
            Console.WriteLine("");
            Console.WriteLine("DAUB12_MATRIX - Fatal error!");
            Console.WriteLine("  Order N must be at least 12 and a multiple of 2.");
            return null;
        }

        double[] a = typeMethods.r8mat_zero_new(n, n);

        double[] c = Coefficients.daub_coefficients(12);

        for (i = 0; i < n - 1; i += 2)
        {
            int j = i;
            a[i + j * n] = c[0];
            j = i + 1;
            a[i + j * n] = c[1];
            j = typeMethods.i4_wrap(i + 2, 0, n - 1);
            a[i + j * n] = c[2];
            j = typeMethods.i4_wrap(i + 3, 0, n - 1);
            a[i + j * n] = c[3];
            j = typeMethods.i4_wrap(i + 4, 0, n - 1);
            a[i + j * n] = c[4];
            j = typeMethods.i4_wrap(i + 5, 0, n - 1);
            a[i + j * n] = c[5];
            j = typeMethods.i4_wrap(i + 6, 0, n - 1);
            a[i + j * n] = c[6];
            j = typeMethods.i4_wrap(i + 7, 0, n - 1);
            a[i + j * n] = c[7];
            j = typeMethods.i4_wrap(i + 8, 0, n - 1);
            a[i + j * n] = c[8];
            j = typeMethods.i4_wrap(i + 9, 0, n - 1);
            a[i + j * n] = c[9];
            j = typeMethods.i4_wrap(i + 10, 0, n - 1);
            a[i + j * n] = c[10];
            j = typeMethods.i4_wrap(i + 11, 0, n - 1);
            a[i + j * n] = c[11];

            j = i;
            a[i + 1 + j * n] = c[11];
            j = i + 1;
            a[i + 1 + j * n] = -c[10];
            j = typeMethods.i4_wrap(i + 2, 0, n - 1);
            a[i + 1 + j * n] = c[9];
            j = typeMethods.i4_wrap(i + 3, 0, n - 1);
            a[i + 1 + j * n] = -c[8];
            j = typeMethods.i4_wrap(i + 4, 0, n - 1);
            a[i + 1 + j * n] = c[7];
            j = typeMethods.i4_wrap(i + 5, 0, n - 1);
            a[i + 1 + j * n] = -c[6];
            j = typeMethods.i4_wrap(i + 6, 0, n - 1);
            a[i + 1 + j * n] = c[5];
            j = typeMethods.i4_wrap(i + 7, 0, n - 1);
            a[i + 1 + j * n] = -c[4];
            j = typeMethods.i4_wrap(i + 8, 0, n - 1);
            a[i + 1 + j * n] = c[3];
            j = typeMethods.i4_wrap(i + 9, 0, n - 1);
            a[i + 1 + j * n] = -c[2];
            j = typeMethods.i4_wrap(i + 10, 0, n - 1);
            a[i + 1 + j * n] = c[1];
            j = typeMethods.i4_wrap(i + 11, 0, n - 1);
            a[i + 1 + j * n] = -c[0];
        }
            
        return a;
    }

    public static double[] daub12_transform(int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DAUB12_TRANSFORM computes the DAUB12 transform of a vector.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 May 2012
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
        //    Output, double DAUB12_TRANSFORM[N], the transformed vector.
        //
    {
        double[] c =  {
                0.1115407433501095E+00,
                0.4946238903984533E+00,
                0.7511339080210959E+00,
                0.3152503517091982E+00,
                -0.2262646939654400E+00,
                -0.1297668675672625E+00,
                0.0975016055873225E+00,
                0.0275228655303053E+00,
                -0.0315820393174862E+00,
                0.0005538422011614E+00,
                0.0047772575109455E+00,
                -0.0010773010853085E+00
            }
            ;
        const int p = 11;

        double[] y = typeMethods.r8vec_copy_new(n, x);
        double[] z = new double[n];

        int m = n;

        while (4 <= m)
        {
            int i;
            for (i = 0; i < m; i++)
            {
                z[i] = 0.0;
            }

            i = 0;

            int j;
            for (j = 0; j < m - 1; j += 2)
            {
                int k;
                for (k = 0; k < p; k += 2)
                {
                    int j0 = typeMethods.i4_wrap(j + k, 0, m - 1);
                    int j1 = typeMethods.i4_wrap(j + k + 1, 0, m - 1);
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

    public static double[] daub12_transform_inverse(int n, double[] y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DAUB12_TRANSFORM_INVERSE inverts the DAUB12 transform of a vector.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 May 2012
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
        //    Output, double DAUB12_TRANSFORM_INVERSE[N], the original vector.
        //
    {
        double[] c =  {
                0.1115407433501095E+00,
                0.4946238903984533E+00,
                0.7511339080210959E+00,
                0.3152503517091982E+00,
                -0.2262646939654400E+00,
                -0.1297668675672625E+00,
                0.0975016055873225E+00,
                0.0275228655303053E+00,
                -0.0315820393174862E+00,
                0.0005538422011614E+00,
                0.0047772575109455E+00,
                -0.0010773010853085E+00
            }
            ;
        const int p = 11;

        double[] x = typeMethods.r8vec_copy_new(n, y);
        double[] z = new double[n];

        int m = 4;
        const int q = (p - 1) / 2;

        while (m <= n)
        {
            int i;
            for (i = 0; i < n; i++)
            {
                z[i] = 0.0;
            }

            int j = 0;

            for (i = -q; i < m / 2 - q; i++)
            {
                int k;
                for (k = 0; k < p; k += 2)
                {
                    int i0 = typeMethods.i4_wrap(i + k / 2, 0, m / 2 - 1);
                    int i1 = typeMethods.i4_wrap(i + m / 2 + k / 2, m / 2, m - 1);
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