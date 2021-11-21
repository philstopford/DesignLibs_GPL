using System;
using Burkardt.Types;

namespace Burkardt.DaubechiesWavelet;

public static class Daub6
{
    public static double[] daub6_matrix(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DAUB6_MATRIX returns the DAUB6 matrix.
        //
        //  Discussion:
        //
        //    The DAUB6 matrix is the Daubechies wavelet transformation matrix 
        //    with 6 coefficients.
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
        //    N must be at least 6 and a multiple of 2.
        //
        //    Output, double DAUB6_MATRIX[N*N], the matrix.
        //
    {
        int i;

        if (n < 6 || n % 2 != 0)
        {
            Console.WriteLine("");
            Console.WriteLine("DAUB6_MATRIX - Fatal error!");
            Console.WriteLine("  Order N must be at least 6 and a multiple of 2.");
            return null;
        }

        double[] a = typeMethods.r8mat_zero_new(n, n);

        double[] c = Coefficients.daub_coefficients(6);

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

            j = i;
            a[i + 1 + j * n] = c[5];
            j = i + 1;
            a[i + 1 + j * n] = -c[4];
            j = typeMethods.i4_wrap(i + 2, 0, n - 1);
            a[i + 1 + j * n] = c[3];
            j = typeMethods.i4_wrap(i + 3, 0, n - 1);
            a[i + 1 + j * n] = -c[2];
            j = typeMethods.i4_wrap(i + 4, 0, n - 1);
            a[i + 1 + j * n] = c[1];
            j = typeMethods.i4_wrap(i + 5, 0, n - 1);
            a[i + 1 + j * n] = -c[0];
        }

        return a;
    }

    public static double daub6_scale(int n, double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DAUB6_SCALE recursively evaluates the DAUB6 scaling function.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 May 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the recursion level.
        //
        //    Input, double X, the point at which the function is to 
        //    be evaluated.
        //
        //    Output, double DAUB6_SCALE, the estimated value of the function.
        //
    {
        double[] c =  {
                0.3326705529500826E+00,
                0.8068915093110925E+00,
                0.4598775021184915E+00,
                -0.1350110200102545E+00,
                -0.08544127388202666E+00,
                0.03522629188570953E+00
            }
            ;
        double y;

        switch (n)
        {
            case > 0:
                y = Math.Sqrt(2.0) *
                    (c[0] * daub6_scale(n - 1, 2.0 * x)
                     + c[1] * daub6_scale(n - 1, 2.0 * x - 1.0)
                     + c[2] * daub6_scale(n - 1, 2.0 * x - 2.0)
                     + c[3] * daub6_scale(n - 1, 2.0 * x - 3.0)
                     + c[4] * daub6_scale(n - 1, 2.0 * x - 4.0)
                     + c[5] * daub6_scale(n - 1, 2.0 * x - 5.0));
                break;
            default:
            {
                y = x switch
                {
                    >= 0.0 and < 1.0 => 1.0,
                    _ => 0.0
                };

                break;
            }
        }

        return y;
    }

    public static double[] daub6_transform(int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DAUB6_TRANSFORM computes the DAUB6 transform of a vector.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 April 2012
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
        //    Output, double DAUB6_TRANSFORM[N], the transformed vector.
        //
    {
        double[] c =  {
                0.3326705529500826,
                0.8068915093110925,
                0.4598775021184915,
                -0.1350110200102545,
                -0.08544127388202666,
                0.03522629188570953
            }
            ;
        const int p = 5;

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

    public static double[] daub6_transform_inverse(int n, double[] y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DAUB6_TRANSFORM_INVERSE inverts the DAUB6 transform of a vector.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 April 2012
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
        //    Output, double DAUB6_TRANSFORM_INVERSE[N], the original vector.
        //
    {
        double[] c =  {
                0.3326705529500826,
                0.8068915093110925,
                0.4598775021184915,
                -0.1350110200102545,
                -0.08544127388202666,
                0.03522629188570953
            }
            ;
        const int p = 5;

        double[] x = typeMethods.r8vec_copy_new(n, y);
        double[] z = new double[n];

        int m = 4;
        int q = (p - 1) / 2;

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