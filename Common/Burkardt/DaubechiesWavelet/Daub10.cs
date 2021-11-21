using System;
using Burkardt.Types;

namespace Burkardt.DaubechiesWavelet;

public static class Daub10
{
    public static double[] daub10_matrix(int n)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DAUB10_MATRIX returns the DAUB10 matrix.
        //
        //  Discussion:
        //
        //    The DAUB10 matrix is the Daubechies wavelet transformation matrix 
        //    with 10 coefficients.
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
        //    N must be at least 10 and a multiple of 2.
        //
        //    Output, double DAUB10_MATRIX[N*N], the matrix.
        //
    {
        int i;

        if (n < 10 || n % 2 != 0)
        {
            Console.WriteLine("");
            Console.WriteLine("DAUB10_MATRIX - Fatal error!");
            Console.WriteLine("  Order N must be at least 10 and a multiple of 2.");
            return null;
        }

        double[] a = typeMethods.r8mat_zero_new(n, n);

        double[] c = Coefficients.daub_coefficients(10);

        for (i = 0; i < n - 1; i += 2)
        {
            int j = i;
            a[i + j * n] = c[0];
            j = i + 1;
            a[i + j * n] = c[1];
            j =  typeMethods.i4_wrap(i + 2, 0, n - 1);
            a[i + j * n] = c[2];
            j =  typeMethods.i4_wrap(i + 3, 0, n - 1);
            a[i + j * n] = c[3];
            j =  typeMethods.i4_wrap(i + 4, 0, n - 1);
            a[i + j * n] = c[4];
            j =  typeMethods.i4_wrap(i + 5, 0, n - 1);
            a[i + j * n] = c[5];
            j =  typeMethods.i4_wrap(i + 6, 0, n - 1);
            a[i + j * n] = c[6];
            j =  typeMethods.i4_wrap(i + 7, 0, n - 1);
            a[i + j * n] = c[7];
            j =  typeMethods.i4_wrap(i + 8, 0, n - 1);
            a[i + j * n] = c[8];
            j =  typeMethods.i4_wrap(i + 9, 0, n - 1);
            a[i + j * n] = c[9];

            j = i;
            a[i + 1 + j * n] = c[9];
            j = i + 1;
            a[i + 1 + j * n] = -c[8];
            j =  typeMethods.i4_wrap(i + 2, 0, n - 1);
            a[i + 1 + j * n] = c[7];
            j =  typeMethods.i4_wrap(i + 3, 0, n - 1);
            a[i + 1 + j * n] = -c[6];
            j =  typeMethods.i4_wrap(i + 4, 0, n - 1);
            a[i + 1 + j * n] = c[5];
            j =  typeMethods.i4_wrap(i + 5, 0, n - 1);
            a[i + 1 + j * n] = -c[4];
            j =  typeMethods.i4_wrap(i + 6, 0, n - 1);
            a[i + 1 + j * n] = c[3];
            j =  typeMethods.i4_wrap(i + 7, 0, n - 1);
            a[i + 1 + j * n] = -c[2];
            j =  typeMethods.i4_wrap(i + 8, 0, n - 1);
            a[i + 1 + j * n] = c[1];
            j =  typeMethods.i4_wrap(i + 9, 0, n - 1);
            a[i + 1 + j * n] = -c[0];
        }
            
        return a;
    }

    public static double daub10_scale(int n, double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DAUB10_SCALE recursively evaluates the DAUB10 scaling function.
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
        //    Output, double DAUB10_SCALE, the estimated value of the function.
        //
    {
        double[] c =  {
                0.1601023979741929E+00,
                0.6038292697971895E+00,
                0.7243085284377726E+00,
                0.1384281459013203E+00,
                -0.2422948870663823E+00,
                -0.0322448695846381E+00,
                0.0775714938400459E+00,
                -0.0062414902127983E+00,
                -0.0125807519990820E+00,
                0.0033357252854738E+00
            }
            ;
        double y;

        switch (n)
        {
            case > 0:
                y = Math.Sqrt(2.0) *
                    (c[0] * daub10_scale(n - 1, 2.0 * x)
                     + c[1] * daub10_scale(n - 1, 2.0 * x - 1.0)
                     + c[2] * daub10_scale(n - 1, 2.0 * x - 2.0)
                     + c[3] * daub10_scale(n - 1, 2.0 * x - 3.0)
                     + c[4] * daub10_scale(n - 1, 2.0 * x - 4.0)
                     + c[5] * daub10_scale(n - 1, 2.0 * x - 5.0)
                     + c[6] * daub10_scale(n - 1, 2.0 * x - 6.0)
                     + c[7] * daub10_scale(n - 1, 2.0 * x - 7.0)
                     + c[8] * daub10_scale(n - 1, 2.0 * x - 8.0)
                     + c[9] * daub10_scale(n - 1, 2.0 * x - 9.0));
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

    public static double[] daub10_transform(int n, double[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DAUB10_TRANSFORM computes the DAUB10 transform of a vector.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 May 2012
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
        //    Output, double DAUB10_TRANSFORM[N], the transformed vector.
        //
    {
        double[] c =  {
                1.601023979741929E-01,
                6.038292697971896E-01,
                7.243085284377729E-01,
                1.384281459013207E-01,
                -2.422948870663820E-01,
                -3.224486958463837E-02,
                7.757149384004571E-02,
                -6.241490212798274E-03,
                -1.258075199908199E-02,
                3.335725285473771E-03
            }
            ;
        const int p = 9;

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

    public static double[] daub10_transform_inverse(int n, double[] y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DAUB10_TRANSFORM_INVERSE inverts the DAUB10 transform of a vector.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 May 2012
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
        //    Output, double DAUB10_TRANSFORM_INVERSE[N], the original vector.
        //
    {
        double[] c =  {
                1.601023979741929E-01,
                6.038292697971896E-01,
                7.243085284377729E-01,
                1.384281459013207E-01,
                -2.422948870663820E-01,
                -3.224486958463837E-02,
                7.757149384004571E-02,
                -6.241490212798274E-03,
                -1.258075199908199E-02,
                3.335725285473771E-03
            }
            ;
        const int p = 9;

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