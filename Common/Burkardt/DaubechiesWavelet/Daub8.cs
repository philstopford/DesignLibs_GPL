using System;
using Burkardt.Types;

namespace Burkardt.DaubechiesWavelet;

public static class Daub8
{
    public static double[] daub8_matrix(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DAUB8_MATRIX returns the DAUB8 matrix.
        //
        //  Discussion:
        //
        //    The DAUB8 matrix is the Daubechies wavelet transformation matrix 
        //    with 8 coefficients.
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
        //    N must be at least 8 and a multiple of 2.
        //
        //    Output, double DAUB8_MATRIX[N*N], the matrix.
        //
    {
        int i;

        if (n < 8 || n % 2 != 0)
        {
            Console.WriteLine("");
            Console.WriteLine("DAUB8_MATRIX - Fatal error!");
            Console.WriteLine("  Order N must be at least 8 and a multiple of 2.");
            return null;
        }

        double[] a = typeMethods.r8mat_zero_new(n, n);

        double[] c = Coefficients.daub_coefficients(8);

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

            j = i;
            a[i + 1 + j * n] = c[7];
            j = i + 1;
            a[i + 1 + j * n] = -c[6];
            j = typeMethods.i4_wrap(i + 2, 0, n - 1);
            a[i + 1 + j * n] = c[5];
            j = typeMethods.i4_wrap(i + 3, 0, n - 1);
            a[i + 1 + j * n] = -c[4];
            j = typeMethods.i4_wrap(i + 4, 0, n - 1);
            a[i + 1 + j * n] = c[3];
            j = typeMethods.i4_wrap(i + 5, 0, n - 1);
            a[i + 1 + j * n] = -c[2];
            j = typeMethods.i4_wrap(i + 6, 0, n - 1);
            a[i + 1 + j * n] = c[1];
            j = typeMethods.i4_wrap(i + 7, 0, n - 1);
            a[i + 1 + j * n] = -c[0];
        }

        return a;
    }

    public static double daub8_scale(int n, double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DAUB8_SCALE recursively evaluates the DAUB8 scaling function.
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
        //    Output, double DAUB8_SCALE, the estimated value of the function.
        //
    {
        double[] c =  {
                0.2303778133088964E+00,
                0.7148465705529154E+00,
                0.6308807679298587E+00,
                -0.0279837694168599E+00,
                -0.1870348117190931E+00,
                0.0308413818355607E+00,
                0.0328830116668852E+00,
                -0.0105974017850690E+00
            }
            ;
        double y;

        switch (n)
        {
            case > 0:
                y = Math.Sqrt(2.0) *
                    (c[0] * daub8_scale(n - 1, 2.0 * x)
                     + c[1] * daub8_scale(n - 1, 2.0 * x - 1.0)
                     + c[2] * daub8_scale(n - 1, 2.0 * x - 2.0)
                     + c[3] * daub8_scale(n - 1, 2.0 * x - 3.0)
                     + c[4] * daub8_scale(n - 1, 2.0 * x - 4.0)
                     + c[5] * daub8_scale(n - 1, 2.0 * x - 5.0)
                     + c[6] * daub8_scale(n - 1, 2.0 * x - 6.0)
                     + c[7] * daub8_scale(n - 1, 2.0 * x - 7.0));
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

    public static double[] daub8_transform(int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DAUB8_TRANSFORM computes the DAUB8 transform of a vector.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 April 2012
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
        //    Output, double DAUB8_TRANSFORM[N], the transformed vector.
        //
    {
        double[] c =  {
                0.2303778133088964,
                0.7148465705529154,
                0.6308807679298587,
                -0.02798376941685985,
                -0.1870348117190931,
                0.03084138183556076,
                0.03288301166688519,
                -0.01059740178506903
            }
            ;
        const int p = 7;

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

    public static double[] daub8_transform_inverse(int n, double[] y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DAUB8_TRANSFORM_INVERSE inverts the DAUB8 transform of a vector.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 April 2012
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
        //    Output, double DAUB8_TRANSFORM_INVERSE[N], the original vector.
        //
    {
        double[] c =  {
                0.2303778133088964,
                0.7148465705529154,
                0.6308807679298587,
                -0.02798376941685985,
                -0.1870348117190931,
                0.03084138183556076,
                0.03288301166688519,
                -0.01059740178506903
            }
            ;
        const int p = 7;

        double[] x = typeMethods.r8vec_copy_new(n, y);
        double[] z = new double[n];

        int m = 4;
        int q = (p - 1) / 2;

        while (m <= n)
        {
            int 
                i;
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