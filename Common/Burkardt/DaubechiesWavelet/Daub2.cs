using System;
using Burkardt.Types;

namespace Burkardt.DaubechiesWavelet
{
    public static class Daub2
    {
        public static double[] daub2_matrix(int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DAUB2_MATRIX returns the DAUB2 matrix.
            //
            //  Discussion:
            //
            //    The DAUB2 matrix is the Daubechies wavelet transformation matrix 
            //    with 2 coefficients.
            //
            //    The DAUB2 matrix is also known as the Haar matrix.
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
            //    N must be at least 2 and a multiple of 2.
            //
            //    Output, double DAUB2_MATRIX[N*N], the matrix.
            //
        {
            double[] a;
            double[] c;
            int i;

            if (n < 2 || (n % 2) != 0)
            {
                Console.WriteLine("");
                Console.WriteLine("DAUB2_MATRIX - Fatal error!");
                Console.WriteLine("  Order N must be at least 2 and a multiple of 2.");
                return null;
            }

            a = typeMethods.r8mat_zero_new(n, n);

            c = Coefficients.daub_coefficients(2);

            for (i = 0; i < n - 1; i = i + 2)
            {
                a[i + i * n] = c[0];
                a[i + (i + 1) * n] = c[1];

                a[i + 1 + i * n] = c[1];
                a[i + 1 + (i + 1) * n] = -c[0];
            }

            return a;
        }

        public static double daub2_scale(int n, double x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DAUB2_SCALE recursively evaluates the DAUB2 scaling function.
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
            //    Output, double DAUB2_SCALE, the estimated value of the function.
            //
        {
            double[] c =  {
                7.071067811865475E-01,
                7.071067811865475E-01
            }
            ;
            double y;

            if (0 < n)
            {
                y = Math.Sqrt(2.0) *
                    (c[0] * daub2_scale(n - 1, 2.0 * x)
                     + c[1] * daub2_scale(n - 1, 2.0 * x - 1.0));
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

        public static double[] daub2_transform(int n, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DAUB2_TRANSFORM computes the DAUB2 transform of a vector.
            //
            //  Discussion:
            //
            //    DAUB2 is better known as the Haar transform.
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
            //    N must be a power of 2.
            //
            //    Input, double X[N], the vector to be transformed. 
            //
            //    Output, double DAUB2_TRANSFORM[N], the transformed vector.
            //
        {
            double[] c =  {
                7.071067811865475E-01,
                7.071067811865475E-01
            }
            ;
            int i;
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

            while (2 <= m)
            {
                m = m / 2;

                for (i = 0; i < m; i++)
                {
                    z[i] = c[0] * (y[2 * i] + y[2 * i + 1]);
                    z[i + m] = c[1] * (y[2 * i] - y[2 * i + 1]);
                }

                for (i = 0; i < 2 * m; i++)
                {
                    y[i] = z[i];
                }
            }
            
            return y;
        }

        public static double[] daub2_transform_inverse(int n, double[] y)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DAUB2_TRANSFORM_INVERSE inverts the DAUB2 transform of a vector.
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
            //    N must be a power of 2.
            //
            //    Input, double Y[N], the transformed vector.
            //
            //    Output, double DAUB2_TRANSFORM_INVERSE[N], the original vector.
            //
        {
            double[] c =  {
                7.071067811865475E-01,
                7.071067811865475E-01
            }
            ;
            int i;
            int m;
            double[] x;
            double[] z;

            x = typeMethods.r8vec_copy_new(n, y);
            z = new double[n];
            for (i = 0; i < n; i++)
            {
                z[i] = 0.0;
            }

            m = 1;

            while (m * 2 <= n)
            {
                for (i = 0; i < m; i++)
                {
                    z[2 * i] = c[0] * (x[i] + x[i + m]);
                    z[2 * i + 1] = c[1] * (x[i] - x[i + m]);
                }

                for (i = 0; i < 2 * m; i++)
                {
                    x[i] = z[i];
                }

                m = m * 2;
            }
            
            return x;
        }
    }
}