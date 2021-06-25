using System;
using Burkardt.Types;

namespace Burkardt.DaubechiesWavelet
{
    public static class Daub4
    {
        public static double[] daub4_matrix(int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DAUB4_MATRIX returns the DAUB4 matrix.
            //
            //  Discussion:
            //
            //    The DAUB4 matrix is the Daubechies wavelet transformation matrix 
            //    with 4 coefficients.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    10 May 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the order of the matrix.
            //    N must be at least 4 and a multiple of 2.
            //
            //    Output, double DAUB4_MATRIX[N*N], the matrix.
            //
        {
            double[] a;
            double[] c;
            int i;
            int j;

            if (n < 4 || (n % 2) != 0)
            {
                Console.WriteLine("");
                Console.WriteLine("DAUB4_MATRIX - Fatal error!");
                Console.WriteLine("  Order N must be at least 4 and a multiple of 2.");
                return null;
            }

            a = typeMethods.r8mat_zero_new(n, n);
    
            c = Coefficients.daub_coefficients(4);

            for (i = 0; i < n - 1; i = i + 2)
            {
                j = i;
                a[i + j * n] = c[0];
                j = i + 1;
                a[i + j * n] = c[1];
                j = typeMethods.i4_wrap(i + 2, 0, n - 1);
                a[i + j * n] = c[2];
                j = typeMethods.i4_wrap(i + 3, 0, n - 1);
                a[i + j * n] = c[3];

                j = i;
                a[i + 1 + j * n] = c[3];
                j = i + 1;
                a[i + 1 + j * n] = -c[2];
                j = typeMethods.i4_wrap(i + 2, 0, n - 1);
                a[i + 1 + j * n] = c[1];
                j = typeMethods.i4_wrap(i + 3, 0, n - 1);
                a[i + 1 + j * n] = -c[0];
            }

            return a;
        }

        public static double daub4_scale(int n, double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DAUB4_SCALE recursively evaluates the DAUB4 scaling function.
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
            //    Output, double DAUB4_SCALE, the estimated value of the function.
            //
        {
            double[] c =  {
                0.4829629131445341E+00,
                0.8365163037378079E+00,
                0.2241438680420133E+00,
                -0.1294095225512603E+00
            }
            ;
            double y;

            if (0 < n)
            {
                y = Math.Sqrt(2.0) *
                    (c[0] * daub4_scale(n - 1, 2.0 * x)
                     + c[1] * daub4_scale(n - 1, 2.0 * x - 1.0)
                     + c[2] * daub4_scale(n - 1, 2.0 * x - 2.0)
                     + c[3] * daub4_scale(n - 1, 2.0 * x - 3.0));
            }
            else if (0.0 <= x && x < 1.0)
            {
                y = 1.0;
            }
            else
            {
                y = 0.0;
            }

            return y;
        }

        public static double[] daub4_transform(int n, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DAUB4_TRANSFORM computes the DAUB4 transform of a vector.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    28 April 2012
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
            //    Output, double DAUB4_TRANSFORM[N], the transformed vector.
            //
        {
            double[] c =  {
                0.4829629131445341,
                0.8365163037378079,
                0.2241438680420133,
                -0.1294095225512603
            }
            ;
            int i;
            int j;
            int j0;
            int j1;
            int j2;
            int j3;
            int m;
            double[] y;
            double[] z;

            y = typeMethods.r8vec_copy_new(n, x);
            z = new double[n];
            for (i = 0; i < n; i++)
            {
                z[i] = 0.0;
            }

            m = n;

            while (4 <= m)
            {
                i = 0;

                for (j = 0; j < m - 1; j = j + 2)
                {
                    j0 = typeMethods.i4_wrap(j, 0, m - 1);
                    j1 = typeMethods.i4_wrap(j + 1, 0, m - 1);
                    j2 = typeMethods.i4_wrap(j + 2, 0, m - 1);
                    j3 = typeMethods.i4_wrap(j + 3, 0, m - 1);

                    z[i] = c[0] * y[j0] + c[1] * y[j1]
                                        + c[2] * y[j2] + c[3] * y[j3];

                    z[i + m / 2] = c[3] * y[j0] - c[2] * y[j1]
                        + c[1] * y[j2] - c[0] * y[j3];

                    i = i + 1;
                }

                for (i = 0; i < m; i++)
                {
                    y[i] = z[i];
                }

                m = m / 2;
            }
            
            return y;
        }

        public static double[] daub4_transform_inverse(int n, double[] y)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DAUB4_TRANSFORM_INVERSE inverts the DAUB4 transform of a vector.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    28 April 2012
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
            //    Output, double DAUB4_TRANSFORM_INVERSE[N], the original vector.
            //
        {
            double[] c =  {
                0.4829629131445341,
                0.8365163037378079,
                0.2241438680420133,
                -0.1294095225512603
            }
            ;
            int i;
            int i0;
            int i1;
            int i2;
            int i3;
            int j;
            int m;
            double[] x;
            double[] z;

            x = typeMethods.r8vec_copy_new(n, y);
            z = new double[n];
            for (i = 0; i < n; i++)
            {
                z[i] = 0.0;
            }

            m = 4;

            while (m <= n)
            {
                j = 0;

                for (i = 0; i < m / 2; i++)
                {
                    i0 = typeMethods.i4_wrap(i - 1, 0, m / 2 - 1);
                    i2 = typeMethods.i4_wrap(i, 0, m / 2 - 1);

                    i1 = typeMethods.i4_wrap(i + m / 2 - 1, m / 2, m - 1);
                    i3 = typeMethods.i4_wrap(i + m / 2, m / 2, m - 1);

                    z[j] = c[2] * x[i0] + c[1] * x[i1]
                                        + c[0] * x[i2] + c[3] * x[i3];

                    z[j + 1] = c[3] * x[i0] - c[0] * x[i1]
                        + c[1] * x[i2] - c[2] * x[i3];

                    j = j + 2;
                }

                for (i = 0; i < m; i++)
                {
                    x[i] = z[i];
                }

                m = m * 2;
            }

            return x;
        }
    }
}