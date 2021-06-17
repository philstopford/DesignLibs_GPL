using System;
using System.Numerics;
using Burkardt.Types;

namespace Burkardt.FourierTransform
{
    public static partial class Slow
    {
        public static double[] r8vec_sftb(int n, double azero, double[] a, double[] b)

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

        public static void r8vec_sftf(int n, double[] r, ref double azero, ref double[] a, ref double[] b)

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

        public static double[] r8vec_sct(int n, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_SCT computes a "slow" cosine transform of an R8VEC.
            //
            //  Discussion:
            //
            //    This routine is provided for illustration and testing.  It is inefficient
            //    relative to optimized routines that use fast Fourier techniques.
            //
            //      Y(1) = Sum ( 1 <= J <= N ) X(J)
            //
            //      For 2 <= I <= N-1:
            //
            //        Y(I) = 2 * Sum ( 1 <= J <= N ) X(J)
            //          * cos ( PI * ( I - 1 ) * ( J - 1 ) / ( N - 1 ) )
            //
            //      Y(N) = Sum ( X(1:N:2) ) - Sum ( X(2:N:2) )
            //
            //    Applying the routine twice in succession should yield the original data,
            //    multiplied by 2 * ( N + 1 ).  This is a good check for correctness
            //    and accuracy.
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
            //    Input, double X[N], the data sequence.
            //
            //    Output, double SCT[N], the transformed data.
            //
        {
            double angle;
            int i;
            int j;
            double[] y;

            y = new double[n];

            for (i = 0; i < n; i++)
            {
                y[i] = x[0] / 2.0;

                for (j = 1; j < n - 1; j++)
                {
                    angle = Math.PI * (double) ((i * j) % (2 * (n - 1)))
                            / (double) (n - 1);
                    y[i] = y[i] + x[j] * Math.Cos(angle);
                }

                j = n - 1;

                angle = Math.PI * (double) ((i * j) % (2 * (n - 1)))
                        / (double) (n - 1);

                y[i] = y[i] + x[j] * Math.Cos(angle) / 2.0;
            }

            for (i = 0; i < n; i++)
            {
                y[i] = 2.0 * y[i] * Math.Sqrt((double) (n) / (double) (n - 1));
            }

            return y;
        }

        public static double[] r8vec_sht(int n, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_SHT computes a "slow" Hartley transform of an R8VEC.
            //
            //  Discussion:
            //
            //    The discrete Hartley transform B of a set of data A is
            //
            //      B(I) = 1/sqrt(N) * Sum ( 0 <= J <= N-1 ) A(J) * CAS(2*PI*I*J/N)
            //
            //    Here, the data and coefficients are indexed from 0 to N-1.
            //
            //    With the above normalization factor of 1/sqrt(N), the Hartley
            //    transform is its own inverse.
            //
            //    This routine is provided for illustration and testing.  It is inefficient
            //    relative to optimized routines.
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
            //  Reference:
            //
            //    Ralph Hartley,
            //    A More Symmetrical Fourier Analysis Applied to Transmission Problems,
            //    Proceedings of the Institute of Radio Engineers,
            //    Volume 30, pages 144-150, 1942.
            //
            //  Parameters:
            //
            //    Input, int N, the number of data values.
            //
            //    Input, double A[N], the data to be transformed.
            //
            //    Output, double SHT[N], the transformed data.
            //
        {
            double[] b;
            int i;
            int j;
            double theta;

            b = new double[n];

            for (i = 0; i < n; i++)
            {
                b[i] = 0.0;
            }

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    theta = 2.0 * Math.PI * (double) ((i * j) % n) / (double) (n);
                    b[i] = b[i] + a[j] * (Math.Cos(theta) + Math.Sin(theta));
                }
            }

            for (i = 0; i < n; i++)
            {
                b[i] = b[i] / Math.Sqrt((double) (n));
            }

            return b;
        }

        public static double[] r8vec_sqctb(int n, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_SQCTB computes a "slow" quarter cosine transform backward of an R8VEC.
            //
            //  Discussion:
            //
            //    This routine is provided for illustration and testing.  It is inefficient
            //    relative to optimized routines that use fast Fourier techniques.
            //
            //    For 0 <= I <= N-1,
            //
            //      Y(I) = X(0) + 2 Sum ( 1 <= J <= N-1 ) X(J) * cos ( PI * J * (I+1/2) / N )
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
            //  Reference:
            //
            //    William Briggs, Van Emden Henson,
            //    The Discrete Fourier Transform,
            //    SIAM, 1995,
            //    LC: QA403.5 B75
            //
            //  Parameters:
            //
            //    Input, int N, the number of data values.
            //
            //    Input, double X[N], the data sequence.
            //
            //    Output, double SQCTB[N], the transformed data.
            //
        {
            int i;
            int j;
            double theta;
            double[] y;

            y = new double[n];

            for (i = 0; i < n; i++)
            {
                y[i] = x[0];
            }

            for (i = 0; i < n; i++)
            {
                for (j = 1; j < n; j++)
                {
                    theta = 0.5 * Math.PI * (double) (j * (2 * i + 1)) / (double) (n);
                    y[i] = y[i] + 2.0 * x[j] * Math.Cos(theta);
                }
            }

            return y;
        }

        public static double[] r8vec_sqctf(int n, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_SQCTF computes a "slow" quarter cosine transform forward of an R8VEC.
            //
            //  Discussion:
            //
            //    This routine is provided for illustration and testing.  It is inefficient
            //    relative to optimized routines that use fast Fourier techniques.
            //
            //    For 0 <= I <= N-1,
            //
            //      Y(I) = (1/N) Sum ( 0 <= J <= N-1 ) X(J) * cos ( PI * I * (J+1/2) / N )
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
            //  Reference:
            //
            //    William Briggs, Van Emden Henson,
            //    The Discrete Fourier Transform,
            //    SIAM, 1995,
            //    QA403.5 B75
            //
            //  Parameters:
            //
            //    Input, int N, the number of data values.
            //
            //    Input, double X[N], the data sequence.
            //
            //    Output, double SQCTF[N], the transformed data.
            //
        {
            int i;
            int j;
            double theta;
            double[] y;

            y = new double[n];

            for (i = 0; i < n; i++)
            {
                y[i] = 0.0;
            }

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    theta = 0.5 * Math.PI * (double) (i * (2 * j + 1)) / (double) (n);
                    y[i] = y[i] + x[j] * Math.Cos(theta);
                }
            }

            for (i = 0; i < n; i++)
            {
                y[i] = y[i] / (double) (n);
            }

            return y;
        }

        public static double[] r8vec_sqstb(int n, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_SQSTB computes a "slow" quarter sine transform backward of an R8VEC.
            //
            //  Discussion:
            //
            //    This routine is provided for illustration and testing.  It is inefficient
            //    relative to optimized routines that use fast Fourier techniques.
            //
            //    For 0 <= I <= N-1,
            //
            //      Y(I) = -2 Sum ( 1 <= J <= N-1 ) X(J) * sin ( PI * J * (I+1/2) / N )
            //             - X(N) * cos ( pi * I )
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
            //  Reference:
            //
            //    William Briggs, Van Emden Henson,
            //    The Discrete Fourier Transform,
            //    SIAM, 1995,
            //    QA403.5 B75
            //
            //  Parameters:
            //
            //    Input, int N, the number of data values.
            //
            //    Input, double X[N], the data sequence.
            //
            //    Output, double SQSTB[N], the transformed data.
            //
        {
            int i;
            int j;
            double theta;
            double[] y;

            y = new double[n];

            for (i = 0; i < n; i++)
            {
                y[i] = 0.0;
            }

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n - 1; j++)
                {
                    theta = 0.5 * Math.PI * (double) ((j + 1) * (2 * i + 1))
                            / (double) (n);
                    y[i] = y[i] - 2.0 * x[j] * Math.Sin(theta);
                }

                theta = Math.PI * (double) (i);
                y[i] = y[i] - x[n - 1] * Math.Cos(theta);

            }

            return y;
        }

        public static double[] r8vec_sqstf(int n, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_SQSTF computes a "slow" quarter sine transform forward of an R8VEC.
            //
            //  Discussion:
            //
            //    This routine is provided for illustration and testing.  It is inefficient
            //    relative to optimized routines that use fast Fourier techniques.
            //
            //    For 1 <= I <= N,
            //
            //      Y(I) = -(1/N) Sum ( 0 <= J <= N-1 ) X(J) * sin ( PI * I * (J+1/2) / N )
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
            //  Reference:
            //
            //    William Briggs, Van Emden Henson,
            //    The Discrete Fourier Transform,
            //    SIAM, 1995,
            //    QA403.5 B75
            //
            //  Parameters:
            //
            //    Input, int N, the number of data values.
            //
            //    Input, double X[N], the data sequence.
            //
            //    Output, double SQSTF{N], the transformed data.
            //
        {
            int i;
            int j;
            double theta;
            double[] y;

            y = new double[n];

            for (i = 0; i < n; i++)
            {
                y[i] = 0.0;
            }

            for (i = 0; i < n; i++)
            {
                for (j = 0; j < n; j++)
                {
                    theta = 0.5 * Math.PI * (double) ((i + 1) * (2 * j + 1))
                            / (double) (n);
                    y[i] = y[i] + x[j] * Math.Sin(theta);
                }
            }

            for (i = 0; i < n; i++)
            {
                y[i] = -y[i] / (double) (n);
            }

            return y;
        }

        public static double[] r8vec_sst(int n, double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_SST computes a "slow" sine transform of an R8VEC.
            //
            //  Discussion:
            //
            //    This routine is provided for illustration and testing.  It is inefficient
            //    relative to optimized routines that use fast Fourier techniques.
            //
            //    For 1 <= I <= N,
            //
            //      Y(I) = Sum ( 1 <= J <= N ) X(J) * sin ( PI * I * J / ( N + 1 ) )
            //
            //    Applying the routine twice in succession should yield the original data,
            //    multiplied by N / 2.  This is a good check for correctness and accuracy.
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
            //    Input, double X[N], the data sequence.
            //
            //    Output, double SST[N], the transformed data.
            //
        {
            int i;
            int j;
            double[] theta;
            double[] y;

            theta = new double[n];

            for (i = 0; i < n; i++)
            {
                theta[i] = Math.PI * (double) (i + 1) / (double) (n + 1);
            }

            y = new double[n];

            for (i = 0; i < n; i++)
            {
                y[i] = 0.0;
            }

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < n; i++)
                {
                    y[i] = (y[i] + 2.0 * x[j] * Math.Sin((double) (j + 1) * theta[i]));
                }
            }

            return y;
        }

        public static double[] r8vec_swtb(int n, ref double[] s, ref double[] d )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_SWTB computes a "slow" backward wavelet transform of an R8VEC.
        //
        //  Discussion:
        //
        //    This function inverts the D4 Daubechies wavelet.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 June 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of data values.
        //
        //    Input/output, double S[(N+1)/2], D[(N+1)/2], the transformed data.
        //    On output, S and D have been overwritten.
        //
        //    Output, double R8VEC_SWTB[N], the original data sequence.
        //
        {
            int i;
            int im1;
            int ip1;
            int n2;
            int np1h;
            double[] x;
            double[] y;

            if ((n % 2) == 1)
            {
                n2 = n + 1;
            }
            else
            {
                n2 = n;
            }

            np1h = (n + 1) / 2;

            for (i = 0; i < np1h; i++)
            {
                d[i] = d[i] / ((Math.Sqrt(3.0) + 1.0) / Math.Sqrt(2.0));
                s[i] = s[i] / ((Math.Sqrt(3.0) - 1.0) / Math.Sqrt(2.0));
            }

            for (i = 0; i < np1h; i++)
            {
                ip1 = typeMethods.i4_wrap(i + 1, 0, np1h - 1);
                s[i] = s[i] + d[ip1];
            }

            y = new double[n2];

            for (i = 0; i < np1h; i++)
            {
                im1 = typeMethods.i4_wrap(i - 1, 0, np1h - 1);
                y[2 * i + 1] = d[i] + Math.Sqrt(3.0) / 4.0 * s[i]
                                    + (Math.Sqrt(3.0) - 2.0) / 4.0 * s[im1];
            }

            for (i = 0; i < np1h; i++)
            {
                y[2 * i] = s[i] - Math.Sqrt(3.0) * y[2 * i + 1];
            }

            x = new double[n];

            for (i = 0; i < n; i++)
            {
                x[i] = y[i];
            }
            return x;
        }

        public static void r8vec_swtf(int n, double[] x, ref double[] s, ref double[] d )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_SWTF computes a "slow" forward wavelet transform of an R8VEC.
        //
        //  Discussion:
        //
        //    This function applies the D4 Daubechies wavelet.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 June 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of data values.
        //
        //    Input, double X[N], the data sequence.
        //
        //    Output, double S[(N+1)/2], D[(N+1)/2], the transformed data.
        //
        {
            int i;
            int im1;
            int ip1;
            int n2;
            int np1h;
            double[] y;

            if ((n % 2) == 1)
            {
                n2 = n + 1;
            }
            else
            {
                n2 = n;
            }

            y = new double[n2];

            for (i = 0; i < n; i++)
            {
                y[i] = x[i];
            }

            if (n < n2)
            {
                y[n] = 0.0;
            }

            np1h = (n + 1) / 2;

            for (i = 0; i < np1h; i++)
            {
                s[i] = y[2 * i] + Math.Sqrt(3.0) * y[2 * i + 1];
            }

            for (i = 0; i < np1h; i++)
            {
                im1 = typeMethods.i4_wrap(i - 1, 0, np1h - 1);
                d[i] = y[2 * i + 1] - Math.Sqrt(3.0) / 4.0 * s[i]
                                    - (Math.Sqrt(3.0) - 2.0) / 4.0 * s[im1];
            }

            for (i = 0; i < np1h; i++)
            {
                ip1 = typeMethods.i4_wrap(i + 1, 0, np1h - 1);
                s[i] = s[i] - d[ip1];
            }

            for (i = 0; i < np1h; i++)
            {
                s[i] = (Math.Sqrt(3.0) - 1.0) / Math.Sqrt(2.0) * s[i];
                d[i] = (Math.Sqrt(3.0) + 1.0) / Math.Sqrt(2.0) * d[i];
            }
       }

    }
}