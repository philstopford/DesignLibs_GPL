using System;
using System.Numerics;

namespace Burkardt.FourierTransform;

public static partial class Slow
{
    public static Complex[] c4mat_sftb(int n1, int n2, Complex[] y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C4MAT_SFTB computes a "slow" backward Fourier transform of a C4MAT.
        //
        //  Discussion:
        //
        //    SFTF and SFTB are inverses of each other.  If we begin with data
        //    X and apply SFTF to get Y, and then apply SFTB to Y,
        //    we should get back the original X.
        //
        //    This routine is provided for illustration and testing.  It is inefficient
        //    relative to optimized routines that use fast Fourier techniques.
        //
        //    For 0 <= I1 <= N1 - 1, 
        //        0 <= I2 <= N2 - 1,
        //
        //      X(I1,I2) = Sum ( 0 <= K2 <= N2 - 1 ) Sum ( 0 <= K1 <= N1 - 1 ) 
        //        Y(K1,K2) * exp ( 2 Math.PI i I1 K1 / N1 ) * exp ( 2 Math.PI i I2 K2 / N2 )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 June 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N1, N2, the number of rows and columns of data.
        //
        //    Input, complex <float> Y[N1*N2], the Fourier coefficients.
        //
        //    Output, complex <float> C4MAT_SFTB[N1*N2], the data.
        //
    {
        int i1;
        int i2;

        Complex[] x = new Complex[n1 * n2];

        for (i2 = 0; i2 < n2; i2++)
        {
            for (i1 = 0; i1 < n1; i1++)
            {
                x[i1 + i2 * n1] = new Complex(0.0, 0.0);
            }
        }

        for (i2 = 0; i2 < n2; i2++)
        {
            int j2;
            for (j2 = 0; j2 < n2; j2++)
            {
                float theta2 = (float) (2.0 * Math.PI * (float) (i2 * j2) / (float) n2);
                Complex cs2 = new(Math.Cos(theta2), -Math.Sin(theta2));
                for (i1 = 0; i1 < n1; i1++)
                {
                    int j1;
                    for (j1 = 0; j1 < n1; j1++)
                    {
                        float theta1 = (float) (2.0 * Math.PI * (float) (i1 * j1) / (float) n1);
                        Complex cs1 = new(Math.Cos(theta1), -Math.Sin(theta1));
                        x[i1 + i2 * n1] += y[j1 + j2 * n1] * cs1 * cs2;
                    }
                }
            }
        }

        for (i2 = 0; i2 < n2; i2++)
        {
            for (i1 = 0; i1 < n1; i1++)
            {
                x[i1 + i2 * n1] /= (float) (n1 * n2);
            }
        }

        return x;
    }

    public static Complex[] c4mat_sftf(int n1, int n2, Complex[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C4MAT_SFTF computes a "slow" forward Fourier transform of a C4MAT.
        //
        //  Discussion:
        //
        //    SFTF and SFTB are inverses of each other.  If we begin with data
        //    X and apply SFTF to get Y, and then apply SFTB to Y, 
        //    we should get back the original X.
        //
        //    This routine is provided for illustration and testing.  It is inefficient
        //    relative to optimized routines that use fast Fourier techniques.
        //
        //    For 0 <= I1 <= N1 - 1, 
        //        0 <= I2 <= N2 - 1,
        //
        //      Y(I1,I2) = Sum ( 0 <= K2 <= N2 - 1 ) Sum ( 0 <= K1 <= N1 - 1 ) 
        //        X(K1,K2) * exp ( - 2 Math.PI i I1 K1 / N1 ) * exp ( - 2 Math.PI i I2 K2 / N2 )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 June 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N1, N2, the number of rows and columns of data.
        //
        //    Input, complex <float> X[N1*N2], the data to be transformed.
        //
        //    Output, complex <float> C4MAT_SFTF[N1*N2], the Fourier coefficients.
        //
    {
        int i1;
        int i2;

        Complex[] y = new Complex[n1 * n2];

        for (i2 = 0; i2 < n2; i2++)
        {
            for (i1 = 0; i1 < n1; i1++)
            {
                y[i1 + i2 * n1] = new Complex(0.0, 0.0);
            }
        }

        for (i2 = 0; i2 < n2; i2++)
        {
            int j2;
            for (j2 = 0; j2 < n2; j2++)
            {
                float theta2 = (float) (-2.0 * Math.PI * (float) (i2 * j2) / (float) n2);
                Complex cs2 = new(Math.Cos(theta2), -Math.Sin(theta2));
                for (i1 = 0; i1 < n1; i1++)
                {
                    int j1;
                    for (j1 = 0; j1 < n1; j1++)
                    {
                        float theta1 = (float) (-2.0 * Math.PI * (float) (i1 * j1) / (float) n1);
                        Complex cs1 = new(Math.Cos(theta1), -Math.Sin(theta1));
                        y[i1 + i2 * n1] += x[j1 + j2 * n1] * cs1 * cs2;
                    }
                }
            }
        }

        return y;
    }

    public static Complex[] c4vec_sftb(int n, Complex[] y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C4VEC_SFTB computes a "slow" backward Fourier transform of a C4VEC.
        //
        //  Discussion:
        //
        //    SFTF and SFTB are inverses of each other.  If we begin with data
        //    X and apply SFTF to get Y, and then apply SFTB to Y,
        //    we should get back the original X.
        //
        //    This routine is provided for illustration and testing.  It is inefficient
        //    relative to optimized routines that use fast Fourier techniques.
        //
        //    For 0 <= I <= N - 1
        //
        //      X(I) = 1/N * Sum ( 0 <= J <= N - 1 ) Y(J) * exp ( 2 Math.PI i I J / N )
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
        //    Input, complex <float> Y[N], the Fourier coefficients.
        //
        //    Output, complex <float> C4VEC_SFTB, the data.
        //
    {
        int i;

        Complex[] x = new Complex[n];

        for (i = 0; i < n; i++)
        {
            x[i] = new Complex(0.0, 0.0);
            int j;
            for (j = 0; j < n; j++)
            {
                float theta = (float) (-2.0 * Math.PI * (float) (i * j) / (float) n);
                x[i] += y[j] * new Complex(Math.Cos(theta), Math.Sin(theta));
            }

            x[i] /= (float) n;
        }

        return x;
    }

    public static Complex[] c4vec_sftf(int n, Complex[] x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    C4VEC_SFTF computes a "slow" forward Fourier transform of a C4VEC.
        //
        //  Discussion:
        //
        //    SFTF and SFTB are inverses of each other.  If we begin with data
        //    X and apply SFTF to get Y, and then apply SFTB to Y, 
        //    we should get back the original X.
        //
        //    This routine is provided for illustration and testing.  It is inefficient
        //    relative to optimized routines that use fast Fourier techniques.
        //
        //    For 0 <= I <= N - 1
        //
        //      Y(I) = Sum ( 0 <= J <= N - 1 ) X(J) * exp ( - 2 Math.PI i I J / N )
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
        //    Input, complex <float> X[N], the data to be transformed.
        //
        //    Output, complex <float> C4VEC_SFTF[N], the Fourier coefficients.
        //
    {
        int i;

        Complex[] y = new Complex[n];

        for (i = 0; i < n; i++)
        {
            y[i] = new Complex(0.0, 0.0);
            int j;
            for (j = 0; j < n; j++)
            {
                float theta = (float) (-2.0 * Math.PI * (float) (i * j) / (float) n);
                y[i] += x[j] * new Complex(Math.Cos(theta), -Math.Sin(theta));
            }
        }

        return y;
    }

}