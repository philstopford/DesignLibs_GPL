using System;

namespace Burkardt.Transform;

public static class Sine
{
    public static double[] sine_transform_data(int n, double[] d)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SINE_TRANSFORM_DATA does a sine transform on a vector of data.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 February 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, integer N, the number of data points.
        //
        //    Input, double D[N], the vector of data.
        //
        //    Output, double SINE_TRANSFORM_DATA[N], the sine transform coefficients.
        //
    {
        int i;

        double[] s = new double[n];

        for (i = 0; i < n; i++)
        {
            s[i] = 0.0;
            int j;
            for (j = 0; j < n; j++)
            {
                double angle = Math.PI * ((i + 1) * (j + 1)) / (n + 1);
                s[i] += Math.Sin(angle) * d[j];
            }

            s[i] *= Math.Sqrt(2.0 / (n + 1));
        }

        return s;
    }

    public static double[] sine_transform_function(int n, double a, double b,
            Func < double, double > f)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SINE_TRANSFORM_FUNCTION does a sine transform on functional data.
        //
        //  Discussion:
        //
        //    The interval [A,B] is divided into N+1 intervals using N+2 points,
        //    which are indexed by 0 through N+1.
        //
        //    The original function F(X) is regarded as the sum of a linear function 
        //    F1 that passes through (A,F(A)) and (B,F(B)), and a function F2
        //    which is 0 at A and B.
        //
        //    The sine transform coefficients for F2 are then computed.
        //
        //    To recover the interpolant of F(X), it is necessary to combine the
        //    linear part F1 with the sine transform interpolant:
        //
        //      Interp(F)(X) = F1(X) + F2(X)
        //
        //    This can be done by calling SINE_TRANSFORM_INTERPOLANT().
        //    
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 February 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of data points.
        //
        //    Input, double A, B, the interval endpoints.
        //
        //    Input, double F ( double X ), a pointer to the function.
        //
        //    Output, SINE_TRANSFORM_FUNCTION[N], the sine transform coefficients.
        //
    {
        int i;

        double fa = f(a);
        double fb = f(b);

        double[] f2 = new double[n];

        for (i = 0; i < n; i++)
        {
            double x = ((n - i) * a
                        + (i + 1) * b)
                       / (n + 1);

            f2[i] = f(x)
                    - ((b - x) * fa
                       + (x - a) * fb)
                    / (b - a);
        }

        double[] s = new double[n];

        for (i = 0; i < n; i++)
        {
            s[i] = 0.0;

            int j;
            for (j = 0; j < n; j++)
            {
                double angle = Math.PI * ((i + 1) * (j + 1)) / (n + 1);
                s[i] += Math.Sin(angle) * f2[j];
            }

            s[i] *= Math.Sqrt(2.0 / (n + 1));
        }
            
        return s;
    }

    public static double[] sine_transform_interpolant(int n, double a, double b, double fa,
            double fb, double[] s, int nx, double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SINE_TRANSFORM_INTERPOLANT evaluates the sine transform interpolant.
        //
        //  Discussion:
        //
        //    The interval [A,B] is divided into N+1 intervals using N+2 points,
        //    which are indexed by 0 through N+1.
        //
        //    The original function F(X) is regarded as the sum of a linear function 
        //    F1 that passes through (A,F(A)) and (B,F(B)), and a function F2
        //    which is 0 at A and B.
        //
        //    The function F2 has been approximated using the sine transform,
        //    and the interpolant is then evaluated as:
        //
        //      Interp(F)(X) = F1(X) + F2(X)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 February 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of terms in the approximation.
        //
        //    Input, double A, B, the interval over which the approximant 
        //    was defined.
        //
        //    Input, double FA, FB, the function values at A and B.
        //
        //    Input, double S[N], the approximant coefficients.
        //
        //    Input, int NX, the number of evaluation points.
        //
        //    Input, double X[NX], the evaluation points.
        //
        //    Output, double SINE_TRANSFORM_INTERPOLANT[NX], the value of the interpolant.
        //
    {
        int i;

        double[] value = new double[nx];

        for (i = 0; i < nx; i++)
        {
            double f1 = ((b - x[i]) * fa
                         + (x[i] - a) * fb)
                        / (b - a);
            double f2 = 0.0;
            int j;
            for (j = 0; j < n; j++)
            {
                double angle = (j + 1) * (x[i] - a) * Math.PI / (b - a);
                f2 += s[j] * Math.Sin(angle);
            }

            f2 *= Math.Sqrt(2.0 / (n + 1));
            value[i] = f1 + f2;
        }

        return value;
    }
}