using System;
using Burkardt.FourierTransform;
using Burkardt.Types;

namespace Burkardt.Noise;

public static class Colored
{
    public static double[] f_alpha(int n, double q_d, double alpha, ref typeMethods.r8NormalData data, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    F_ALPHA generates a 1/F^ALPHA noise sequence.
        //
        //  Discussion:
        //
        //    Thanks to Miro Stoyanov for pointing out that the second half of
        //    the data returned by the inverse Fourier transform should be
        //    discarded, 24 August 2010.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    24 August 2010
        //
        //  Author:
        //
        //    Original C version by Todd Walter.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Jeremy Kasdin,
        //    Discrete Simulation of Colored Noise and Stochastic Processes
        //    and 1/f^a Power Law Noise Generation,
        //    Proceedings of the IEEE,
        //    Volume 83, Number 5, 1995, pages 802-827.
        //
        //  Parameters:
        //
        //    Input, int N, the number of samples of the noise to generate.
        //
        //    Input, double Q_D, the variance of the noise.
        //
        //    Input, double ALPHA, the exponent for the noise.
        //
        //    Input/output, int *SEED, a seed for the random number generator.
        //
        //    Output, double F_ALPHA[N], a sequence sampled with the given
        //    power law.
        //
    {
        double[] h_a;
        double h_azero = 0;
        double[] h_b;
        double[] hfa;
        int i;
        double[] w_a;
        double w_azero = 0;
        double[] w_b;
        double[] wfa;
        double wi;
        double wr;
        double[] x;
        double[] x2;
            
        //
        //  Set the deviation of the noise.
        //
        q_d = Math.Sqrt(q_d);
        //
        //  Generate the coefficients Hk.
        //
        hfa = new double[2 * n];
        hfa[0] = 1.0;
        for (i = 1; i < n; i++)
        {
            hfa[i] = hfa[i - 1]
                * (0.5 * alpha + (i - 1)) / i;
        }

        for (i = n; i < 2 * n; i++)
        {
            hfa[i] = 0.0;
        }

        //
        //  Fill Wk with white noise.
        //
        wfa = new double[2 * n];

        for (i = 0; i < n; i++)
        {
            wfa[i] = q_d * typeMethods.r8_normal_01(ref data, ref seed);
        }

        for (i = n; i < 2 * n; i++)
        {
            wfa[i] = 0.0;
        }

        //
        //  Perform the discrete Fourier transforms of Hk and Wk.
        //
        h_a = new double[n];
        h_b = new double[n];

        Slow.r8vec_sftf(2 * n, hfa, ref h_azero, ref h_a, ref h_b);

        w_a = new double[n];
        w_b = new double[n];

        Slow.r8vec_sftf(2 * n, wfa, ref w_azero, ref w_a, ref w_b);
        //
        //  Multiply the two complex vectors.
        //
        w_azero *= h_azero;

        for (i = 0; i < n; i++)
        {
            wr = w_a[i];
            wi = w_b[i];
            w_a[i] = wr * h_a[i] - wi * h_b[i];
            w_b[i] = wi * h_a[i] + wr * h_b[i];
        }

        //
        //  This scaling is introduced only to match the behavior
        //  of the Numerical Recipes code...
        //
        w_azero = w_azero * 2 * n;
        for (i = 0; i < n - 1; i++)
        {
            w_a[i] *= n;
            w_b[i] *= n;
        }

        i = n - 1;
        w_a[i] = w_a[i] * 2 * n;
        w_b[i] = w_b[i] * 2 * n;
        //
        //  Take the inverse Fourier transform of the result.
        //
        x2 = Slow.r8vec_sftb(2 * n, w_azero, w_a, w_b);
        //
        //  Only return the first N inverse Fourier transform values.
        //
        x = new double[n];
        for (i = 0; i < n; i++)
        {
            x[i] = x2[i];
        }

        return x;
    }
}