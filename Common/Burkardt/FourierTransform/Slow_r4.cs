using System;

namespace Burkardt.FourierTransform;

public static partial class Slow
{
    public static float[] r4vec_sct(int n, float[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R4VEC_SCT computes a "slow" cosine transform of an R4VEC.
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
        //    10 June 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of data values.
        //
        //    Input, float X[N], the data sequence.
        //
        //    Output, float SCT[N], the transformed data.
        //
    {
        float angle;
        int i;
        int j;
        float[] y;

        y = new float[n];

        for (i = 0; i < n; i++)
        {
            y[i] = x[0] / 2.0f;

            for (j = 1; j < n - 1; j++)
            {
                angle = (float) Math.PI * (i * j % (2 * (n - 1)))
                        / (n - 1);
                y[i] = (float) (y[i] + x[j] * Math.Cos(angle));
            }

            j = n - 1;

            angle = (float) Math.PI * (i * j % (2 * (n - 1)))
                    / (n - 1);

            y[i] = (float) (y[i] + x[j] * Math.Cos(angle) / 2.0);
        }

        for (i = 0; i < n; i++)
        {
            y[i] = (float) (2.0 * y[i] * Math.Sqrt(n / (float) (n - 1)));
        }

        return y;
    }

    public static float[] r4vec_sftb(int n, float azero, float[] a, float[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R4VEC_SFTB computes a "slow" backward Fourier transform of an R4VEC.
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
        //    10 June 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of data values.
        //
        //    Input, float AZERO, the constant Fourier coefficient.
        //
        //    Input, float A[N/2], B[N/2], the Fourier coefficients.
        //
        //    Output, float SFTB[N], the reconstructed data sequence.
        //
    {
        int i;
        int k;
        float[] r;
        float theta;

        r = new float[n];

        for (i = 0; i < n; i++)
        {
            r[i] = azero;
            for (k = 0; k < n / 2; k++)
            {
                theta = (float) ((k + 1) * i * 2 * Math.PI / (float) n);
                r[i] = (float) (r[i] + a[k] * Math.Cos(theta) + b[k] * Math.Sin(theta));
            }
        }

        return r;
    }

    public static void r4vec_sftf(int n, float[] r, ref float azero, ref float[] a, ref float[] b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R4VEC_SFTF computes a "slow" forward Fourier transform of an R4VEC.
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
        //    Input, float R[N], the data to be transformed.
        //
        //    Output, float *AZERO, = sum ( 1 <= I <= N ) R(I) / N.
        //
        //    Output, float A[N/2], B[N/2], the Fourier coefficients.
        //
    {
        int i;
        int j;
        float theta;

        azero = 0.0f;
        for (i = 0; i < n; i++)
        {
            azero += r[i];
        }

        azero /= n;

        for (i = 0; i < n / 2; i++)
        {
            a[i] = 0.0f;
            b[i] = 0.0f;

            for (j = 0; j < n; j++)
            {
                theta = (float) (2 * (i + 1) * j * Math.PI / (float) n);
                a[i] += r[j] * (float) Math.Cos(theta);
                b[i] += r[j] * (float) Math.Sin(theta);
            }

            a[i] /= n;
            b[i] /= n;

            if (n % 2 == 1 || i != n / 2 - 1)
            {
                a[i] = 2.0f * a[i];
                b[i] = 2.0f * b[i];
            }
        }
    }

    public static float[] r4vec_sht(int n, float[] a)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R4VEC_SHT computes a "slow" Hartley transform of an R4VEC.
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
        //    10 June 2010
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
        //    Input, float A[N], the data to be transformed.
        //
        //    Output, float SHT[N], the transformed data.
        //
    {
        float[] b;
        int i;
        int j;
        float theta;

        b = new float[n];

        for (i = 0; i < n; i++)
        {
            b[i] = 0.0f;
        }

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                theta = (float) (2.0f * Math.PI * (float) (i * j % n) / (float) n);
                b[i] = (float) (b[i] + a[j] * (Math.Cos(theta) + Math.Sin(theta)));
            }
        }

        for (i = 0; i < n; i++)
        {
            b[i] = (float) (b[i] / Math.Sqrt((float) n));
        }

        return b;
    }

    public static float[] r4vec_sqctb(int n, float[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R4VEC_SQCTB computes a "slow" quarter cosine transform backward of an R4VEC.
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
        //    10 June 2010
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
        //    Input, float X[N], the data sequence.
        //
        //    Output, float SQCTB[N], the transformed data.
        //
    {
        int i;
        int j;
        float theta;
        float[] y;

        y = new float[n];

        for (i = 0; i < n; i++)
        {
            y[i] = x[0];
        }

        for (i = 0; i < n; i++)
        {
            for (j = 1; j < n; j++)
            {
                theta = (float) (0.5 * Math.PI * (float) (j * (2 * i + 1)) / (float) n);
                y[i] = (float) (y[i] + 2.0 * x[j] * Math.Cos(theta));
            }
        }

        return y;
    }

    public static float[] r4vec_sqctf(int n, float[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R4VEC_SQCTF computes a "slow" quarter cosine transform forward of an R4VEC.
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
        //    10 June 2010
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
        //    Input, float X[N], the data sequence.
        //
        //    Output, float SQCTF[N], the transformed data.
        //
    {
        int i;
        int j;
        float theta;
        float[] y;

        y = new float[n];

        for (i = 0; i < n; i++)
        {
            y[i] = 0.0f;
        }

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                theta = (float) (0.5 * Math.PI * (float) (i * (2 * j + 1)) / (float) n);
                y[i] = (float) (y[i] + x[j] * Math.Cos(theta));
            }
        }

        for (i = 0; i < n; i++)
        {
            y[i] /= n;
        }

        return y;
    }

    public static float[] r4vec_sqstb(int n, float[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R4VEC_SQSTB computes a "slow" quarter sine transform backward of an R4VEC.
        //
        //  Discussion:
        //
        //    This routine is provided for illustration and testing.  It is inefficient
        //    relative to optimized routines that use fast Fourier techniques.
        //
        //    For 0 <= I <= N-1,
        //
        //      Y(I) = -2 Sum ( 1 <= J <= N-1 ) X(J) * sin ( PI * J * (I+1/2) / N )
        //             - X(N) * cos ( Math.PI * I )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 June 2010
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
        //    Input, float X[N], the data sequence.
        //
        //    Output, float SQSTB[N], the transformed data.
        //
    {
        int i;
        int j;
        float theta;
        float[] y;

        y = new float[n];

        for (i = 0; i < n; i++)
        {
            y[i] = 0.0f;
        }

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n - 1; j++)
            {
                theta = (float) (0.5 * Math.PI * (float) ((j + 1) * (2 * i + 1))
                                 / (float) n);
                y[i] = (float) (y[i] - 2.0 * x[j] * Math.Sin(theta));
            }

            theta = (float) (Math.PI * (float) i);
            y[i] = (float) (y[i] - x[n - 1] * Math.Cos(theta));

        }

        return y;
    }

    public static float[] r4vec_sqstf(int n, float[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R4VEC_SQSTF computes a "slow" quarter sine transform forward of an R4VEC.
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
        //    10 June 2010
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
        //    Input, float X[N], the data sequence.
        //
        //    Output, float SQSTF{N], the transformed data.
        //
    {
        int i;
        int j;
        float theta;
        float[] y;

        y = new float[n];

        for (i = 0; i < n; i++)
        {
            y[i] = 0.0f;
        }

        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
            {
                theta = (float) (0.5 * Math.PI * (float) ((i + 1) * (2 * j + 1))
                                 / (float) n);
                y[i] = (float) (y[i] + x[j] * Math.Sin(theta));
            }
        }

        for (i = 0; i < n; i++)
        {
            y[i] = -y[i] / n;
        }

        return y;
    }

    public static float[] r4vec_sst(int n, float[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R4VEC_SST computes a "slow" sine transform of an R4VEC.
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
        //    10 June 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of data values.
        //
        //    Input, float X[N], the data sequence.
        //
        //    Output, float SST[N], the transformed data.
        //
    {
        int i;
        int j;
        float[] theta;
        float[] y;

        theta = new float[n];

        for (i = 0; i < n; i++)
        {
            theta[i] = (float) (Math.PI * (float) (i + 1) / (float) (n + 1));
        }

        y = new float[n];

        for (i = 0; i < n; i++)
        {
            y[i] = 0.0f;
        }

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < n; i++)
            {
                y[i] = (float) (y[i] + 2.0 * x[j] * Math.Sin((j + 1) * theta[i]));
            }
        }

        return y;
    }

}