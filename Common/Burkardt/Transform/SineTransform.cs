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
        double angle;
        int i;
        int j;
            
        double[] s;

        s = new double[n];

        for (i = 0; i < n; i++)
        {
            s[i] = 0.0;
            for (j = 0; j < n; j++)
            {
                angle = Math.PI * ((i + 1) * (j + 1)) / (n + 1);
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
        double angle;
        double[] f2;
        double fa;
        double fb;
        ;
        int i;
        int j;
            
        double[] s;
        double x;

        fa = f(a);
        fb = f(b);

        f2 = new double[n];

        for (i = 0; i < n; i++)
        {
            x = ((n - i) * a
                 + (i + 1) * b)
                / (n + 1);

            f2[i] = f(x)
                    - ((b - x) * fa
                       + (x - a) * fb)
                    / (b - a);
        }

        s = new double[n];

        for (i = 0; i < n; i++)
        {
            s[i] = 0.0;

            for (j = 0; j < n; j++)
            {
                angle = Math.PI * ((i + 1) * (j + 1)) / (n + 1);
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
        double angle;
        double f1;
        double f2;
        int i;
        int j;
            
        double[] value;

        value = new double[nx];

        for (i = 0; i < nx; i++)
        {
            f1 = ((b - x[i]) * fa
                  + (x[i] - a) * fb)
                 / (b - a);
            f2 = 0.0;
            for (j = 0; j < n; j++)
            {
                angle = (j + 1) * (x[i] - a) * Math.PI / (b - a);
                f2 += s[j] * Math.Sin(angle);
            }

            f2 *= Math.Sqrt(2.0 / (n + 1));
            value[i] = f1 + f2;
        }

        return value;
    }
}