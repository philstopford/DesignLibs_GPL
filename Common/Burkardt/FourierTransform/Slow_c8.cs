using System;
using System.Numerics;

namespace Burkardt.FourierTransform
{
    public static partial class Slow
    {
        public static Complex[] c8mat_sftb(int n1, int n2, Complex[] y)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8MAT_SFTB computes a "slow" backward Fourier transform of a C8MAT.
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
            //        Y(K1,K2) * exp ( 2 pi i I1 K1 / N1 ) * exp ( 2 pi i I2 K2 / N2 )
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
            //    Input, complex <double> Y[N1*N2], the Fourier coefficients.
            //
            //    Output, complex <double> C8MAT_SFTB[N1*N2], the data.
            //
        {
            Complex cs1;
            Complex cs2;
            int i1;
            int i2;
            int j1;
            int j2;
            double theta1;
            double theta2;
            Complex[] x;

            x = new Complex[n1 * n2];

            for (i2 = 0; i2 < n2; i2++)
            {
                for (i1 = 0; i1 < n1; i1++)
                {
                    x[i1 + i2 * n1] = new Complex(0.0, 0.0);
                }
            }

            for (i2 = 0; i2 < n2; i2++)
            {
                for (j2 = 0; j2 < n2; j2++)
                {
                    theta2 = 2.0 * Math.PI * (double) (i2 * j2) / (double) (n2);
                    cs2 = new Complex(Math.Cos(theta2), -Math.Sin(theta2));
                    for (i1 = 0; i1 < n1; i1++)
                    {
                        for (j1 = 0; j1 < n1; j1++)
                        {
                            theta1 = 2.0 * Math.PI * (double) (i1 * j1) / (double) (n1);
                            cs1 = new Complex(Math.Cos(theta1), -Math.Sin(theta1));
                            x[i1 + i2 * n1] = x[i1 + i2 * n1] + y[j1 + j2 * n1] * cs1 * cs2;
                        }
                    }
                }
            }

            for (i2 = 0; i2 < n2; i2++)
            {
                for (i1 = 0; i1 < n1; i1++)
                {
                    x[i1 + i2 * n1] = x[i1 + i2 * n1] / (double) (n1 * n2);
                }
            }

            return x;
        }

        public static Complex[] c8mat_sftf(int n1, int n2, Complex[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8MAT_SFTF computes a "slow" forward Fourier transform of a C8MAT.
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
            //        X(K1,K2) * exp ( - 2 pi i I1 K1 / N1 ) * exp ( - 2 pi i I2 K2 / N2 )
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
            //    Input, complex <double> X[N1*N2], the data to be transformed.
            //
            //    Output, complex <double> C8MAT_SFTF[N1*N2], the Fourier coefficients.
            //
        {
            Complex cs1;
            Complex cs2;
            int i1;
            int i2;
            int j1;
            int j2;
            double theta1;
            double theta2;
            Complex[] y;

            y = new Complex[n1 * n2];

            for (i2 = 0; i2 < n2; i2++)
            {
                for (i1 = 0; i1 < n1; i1++)
                {
                    y[i1 + i2 * n1] = new Complex(0.0, 0.0);
                }
            }

            for (i2 = 0; i2 < n2; i2++)
            {
                for (j2 = 0; j2 < n2; j2++)
                {
                    theta2 = -2.0 * Math.PI * (double) (i2 * j2) / (double) (n2);
                    cs2 = new Complex(Math.Cos(theta2), -Math.Sin(theta2));
                    for (i1 = 0; i1 < n1; i1++)
                    {
                        for (j1 = 0; j1 < n1; j1++)
                        {
                            theta1 = -2.0 * Math.PI * (double) (i1 * j1) / (double) (n1);
                            cs1 = new Complex(Math.Cos(theta1), -Math.Sin(theta1));
                            y[i1 + i2 * n1] = y[i1 + i2 * n1] + x[j1 + j2 * n1] * cs1 * cs2;
                        }
                    }
                }
            }

            return y;
        }

        public static Complex[] c8vec_sftb(int n, Complex[] y)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8VEC_SFTB computes a "slow" backward Fourier transform of a C8VEC.
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
            //      X(I) = 1/N * Sum ( 0 <= J <= N - 1 ) Y(J) * exp ( 2 pi i I J / N )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 November 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of data values.
            //
            //    Input, complex <double> Y[N], the Fourier coefficients.
            //
            //    Output, complex <double> C8VEC_SFTB, the data.
            //
        {
            int i;
            int j;
            double theta;
            Complex[] x;

            x = new Complex[n];

            for (i = 0; i < n; i++)
            {
                x[i] = new Complex(0.0, 0.0);
                for (j = 0; j < n; j++)
                {
                    theta = -2.0 * Math.PI * (double) (i * j) / (double) (n);
                    x[i] = x[i] + y[j] * new Complex(Math.Cos(theta), Math.Sin(theta));
                }

                x[i] = x[i] / (double) (n);
            }

            return x;
        }

        public static Complex[] c8vec_sftf(int n, Complex[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    C8VEC_SFTF computes a "slow" forward Fourier transform of a C8VEC.
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
            //      Y(I) = Sum ( 0 <= J <= N - 1 ) X(J) * exp ( - 2 pi i I J / N )
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 November 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of data values.
            //
            //    Input, complex <double> X[N], the data to be transformed.
            //
            //    Output, complex <double> C8VEC_SFTF[N], the Fourier coefficients.
            //
        {
            int i;
            int j;
            double theta;
            Complex[] y;

            y = new Complex[n];

            for (i = 0; i < n; i++)
            {
                y[i] = new Complex(0.0, 0.0);
                for (j = 0; j < n; j++)
                {
                    theta = -2.0 * Math.PI * (double) (i * j) / (double) (n);
                    y[i] = y[i] + x[j] * new Complex(Math.Cos(theta), -Math.Sin(theta));
                }
            }

            return y;
        }

    }
}