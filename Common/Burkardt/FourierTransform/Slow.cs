using System;

namespace Burkardt.FourierTransform
{
    public static class Slow
    {
        public static double[] r8vec_sftb(int n, double azero, double[] a, double[] b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_SFTB computes a "slow" backward Fourier transform of real data.
        //
        //  Discussion:
        //
        //    SFTB and SFTF are inverses of each other.  If we begin with data
        //    AZERO, A, and B, and apply SFTB to it, and then apply SFTF to the
        //    resulting R vector, we should get back the original AZERO, A and B.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 February 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of data values.
        //
        //    Input, double AZERO, the constant Fourier coefficient.
        //
        //    Input, double A[N/2], B[N/2], the Fourier coefficients.
        //
        //    Output, double SFTB[N], the reconstructed data sequence.
        //
        {
            int i;
            int k;
            double[] r;
            double theta;

            r = new double[n];

            for (i = 0; i < n; i++)
            {
                r[i] = azero;
                for (k = 0; k < (n / 2); k++)
                {
                    theta = (double) ((k + 1) * i * 2) * Math.PI / (double) (n);
                    r[i] = r[i] + a[k] * Math.Cos(theta) + b[k] * Math.Sin(theta);
                }
            }

            return r;
        }

        public static void r8vec_sftf(int n, double[] r, ref double azero, ref double[] a, ref double[] b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_SFTF computes a "slow" forward Fourier transform of real data.
        //
        //  Discussion:
        //
        //    SFTF and SFTB are inverses of each other.  If we begin with data
        //    R and apply SFTB to it, and then apply SFTB to the resulting AZERO, 
        //    A, and B, we should get back the original R.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 July 2017
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of data values.
        //
        //    Input, double R[N], the data to be transformed.
        //
        //    Output, double *AZERO, = sum ( 1 <= I <= N ) R(I) / N.
        //
        //    Output, double A[N/2], B[N/2], the Fourier coefficients.
        //
        {
            int i;
            int j;
            double pi = 3.141592653589793;
            double theta;

            azero = 0.0;
            for (i = 0; i < n; i++)
            {
                azero = azero + r[i];
            }

            azero = azero / (double) (n);

            for (i = 0; i < (n / 2); i++)
            {
                a[i] = 0.0;
                b[i] = 0.0;

                for (j = 0; j < n; j++)
                {
                    theta = (double) (2 * (i + 1) * j) * pi / (double) (n);
                    a[i] = a[i] + r[j] * Math.Cos(theta);
                    b[i] = b[i] + r[j] * Math.Sin(theta);
                }

                a[i] = a[i] / (double) (n);
                b[i] = b[i] / (double) (n);

                if ((n % 2) == 1 || i != (n / 2 - 1))
                {
                    a[i] = 2.0 * a[i];
                    b[i] = 2.0 * b[i];
                }
            }

            return;
        }
    }
}