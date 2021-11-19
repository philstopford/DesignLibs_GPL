using System;
using Burkardt.Types;

namespace Burkardt.SolveNS;

public static class BlackScholes
{
    public static double[] asset_path(double s0, double mu, double sigma, double t1, int n,
            ref typeMethods.r8vecNormalData data, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ASSET_PATH simulates the behavior of an asset price over time.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 June 2016
        //
        //  Author:
        //
        //    Original MATLAB version by Desmond Higham.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Desmond Higham,
        //    Black-Scholes for Scientific Computing Students,
        //    Computing in Science and Engineering,
        //    November/December 2004, Volume 6, Number 6, pages 72-79.
        //
        //  Parameters:
        //
        //    Input, double S0, the asset price at time 0.
        //
        //    Input, double MU, the expected growth rate.
        //
        //    Input, double SIGMA, the volatility of the asset.
        //
        //    Input, double T1, the expiry date.
        //
        //    Input, integer N, the number of steps to take between 0 and T1.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double ASSET_PATH[N+1], the option values from time 0 to T1 
        //    in equal steps.
        //
    {
        int i;

        double dt = t1 / n;

        double[] r = typeMethods.r8vec_normal_01_new(n, ref data, ref seed);

        double[] s = new double[n + 1];

        s[0] = s0;
        double p = s0;
        for (i = 1; i <= n; i++)
        {
            p *= Math.Exp((mu - sigma * sigma) * dt + sigma * Math.Sqrt(dt) * r[i - 1]);
            s[i] = p;
        }

        return s;
    }

    public static double binomial(double s0, double e, double r, double sigma, double t1,
            int m)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BINOMIAL uses the binomial method for a European call.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 February 2012
        //
        //  Author:
        //
        //    Original MATLAB version by Desmond Higham.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Desmond Higham,
        //    Black-Scholes for Scientific Computing Students,
        //    Computing in Science and Engineering,
        //    November/December 2004, Volume 6, Number 6, pages 72-79.
        //
        //  Parameters:
        //
        //    Input, double S0, the asset price at time 0.
        //
        //    Input, double E, the exercise price.
        //
        //    Input, double R, the interest rate.
        //
        //    Input, double SIGMA, the volatility of the asset.
        //
        //    Input, double T1, the expiry date.
        //
        //    Input, int M, the number of steps to take 
        //    between 0 and T1.
        //
        //    Output, double BIONOMIAL, the option value at time 0.
        //
    {
        int i;
        int n;
        //
        //  Time stepsize.
        //
        double dt = t1 / m;

        double a = 0.5 * (Math.Exp(-r * dt) + Math.Exp((r + sigma * sigma) * dt));

        double d = a - Math.Sqrt(a * a - 1.0);
        double u = a + Math.Sqrt(a * a - 1.0);

        double p = (Math.Exp(r * dt) - d) / (u - d);

        double[] w = new double[m + 1];

        for (i = 0; i <= m; i++)
        {
            w[i] = Math.Max(s0 * Math.Pow(d, m - i) * Math.Pow(u, i) - e, 0.0);
        }

        //
        //  Trace backwards to get the option value at time 0.
        //
        for (n = m - 1; 0 <= n; n--)
        {
            for (i = 0; i <= n; i++)
            {
                w[i] = (1.0 - p) * w[i] + p * w[i + 1];
            }
        }

        for (i = 0; i < m + 1; i++)
        {
            w[i] = Math.Exp(-r * t1) * w[i];
        }

        double c = w[0];

        return c;
    }

    public static double bsf(double s0, double t0, double e, double r, double sigma, double t1)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BSF evaluates the Black-Scholes formula for a European call.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 February 2012
        //
        //  Author:
        //
        //    Original MATLAB version by Desmond Higham.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Desmond Higham,
        //    Black-Scholes for Scientific Computing Students,
        //    Computing in Science and Engineering,
        //    November/December 2004, Volume 6, Number 6, pages 72-79.
        //
        //  Parameters:
        //
        //    Input, double S0, the asset price at time T0.
        //
        //    Input, double T0, the time at which the asset price is known.
        //
        //    Input, double E, the exercise price.
        //
        //    Input, double R, the interest rate.
        //
        //    Input, double SIGMA, the volatility of the asset.
        //
        //    Input, double T1, the expiry date.
        //
        //    Output, double BSF, the value of the call option.
        //
    {
        double c;

        double tau = t1 - t0;

        switch (tau)
        {
            case > 0.0:
                double d1 = (Math.Log(s0 / e) + (r + 0.5 * sigma * sigma) * tau)
                            / (sigma * Math.Sqrt(tau));

                double d2 = d1 - sigma * Math.Sqrt(tau);

                double n1 = 0.5 * (1.0 + Helpers.Erf(d1 / Math.Sqrt(2.0)));
                double n2 = 0.5 * (1.0 + Helpers.Erf(d2 / Math.Sqrt(2.0)));

                c = s0 * n1 - e * Math.Exp(-r * tau) * n2;
                break;
            default:
                c = Math.Max(s0 - e, 0.0);
                break;
        }

        return c;
    }

    public static double[] forward(double e, double r, double sigma, double t1, int nx,
            int nt, double smax)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FORWARD uses the forward difference method to value a European call option.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 February 2012
        //
        //  Author:
        //
        //    Original MATLAB version by Desmond Higham.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Desmond Higham,
        //    Black-Scholes for Scientific Computing Students,
        //    Computing in Science and Engineering,
        //    November/December 2004, Volume 6, Number 6, pages 72-79.
        //
        //  Parameters:
        //
        //    Input, double E, the exercise price.
        //
        //    Input, double R, the interest rate.
        //
        //    Input, double SIGMA, the volatility of the asset.
        //
        //    Input, double T1, the expiry date.
        //
        //    Input, int NX, the number of "space" steps used to 
        //    divide the interval [0,L].
        //
        //    Input, int NT, the number of time steps.
        //
        //    Input, double SMAX, the maximum value of S to consider.
        //
        //    Output, double U[(NX-1)*(NT+1)], the value of the European 
        //    call option.
        //
    {
        int i;
        int j;

        double dt = t1 / nt;
        double dx = smax / nx;

        double[] a = new double[nx - 1];
        double[] b = new double[nx - 1];
        double[] c = new double[nx - 1];

        for (i = 0; i < nx - 1; i++)
        {
            b[i] = 1.0 - r * dt - dt * Math.Pow(sigma * (i + 1), 2);
        }

        for (i = 0; i < nx - 2; i++)
        {
            c[i] = 0.5 * dt * Math.Pow(sigma * (i + 1), 2) + 0.5 * dt * r * (i + 1);
        }

        for (i = 1; i < nx - 1; i++)
        {
            a[i] = 0.5 * dt * Math.Pow(sigma * (i + 1), 2) - 0.5 * dt * r * (i + 1);
        }

        double[] u = new double[(nx - 1) * (nt + 1)];

        double u0 = 0.0;
        for (i = 0; i < nx - 1; i++)
        {
            u0 += dx;
            u[i + 0 * (nx - 1)] = Math.Max(u0 - e, 0.0);
        }

        for (j = 0; j < nt; j++)
        {
            double t = j * t1 / nt;

            double p = 0.5 * dt * (nx - 1) * (sigma * sigma * (nx - 1) + r)
                       * (smax - e * Math.Exp(-r * t));

            for (i = 0; i < nx - 1; i++)
            {
                u[i + (j + 1) * (nx - 1)] = b[i] * u[i + j * (nx - 1)];
            }

            for (i = 0; i < nx - 2; i++)
            {
                u[i + (j + 1) * (nx - 1)] += c[i] * u[i + 1 + j * (nx - 1)];
            }

            for (i = 1; i < nx - 1; i++)
            {
                u[i + (j + 1) * (nx - 1)] += a[i] * u[i - 1 + j * (nx - 1)];
            }

            u[nx - 2 + (j + 1) * (nx - 1)] += p;
        }

        return u;
    }

    public static double[] mc(double s0, double e, double r, double sigma, double t1, int m,
            ref typeMethods.r8vecNormalData data,
            ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MC uses Monte Carlo valuation on a European call.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 June 2016
        //
        //  Author:
        //
        //    Original MATLAB version by Desmond Higham.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Desmond Higham,
        //    Black-Scholes for Scientific Computing Students,
        //    Computing in Science and Engineering,
        //    November/December 2004, Volume 6, Number 6, pages 72-79.
        //
        //  Parameters:
        //
        //    Input, double S0, the asset price at time 0.
        //
        //    Input, double E, the exercise price.
        //
        //    Input, double R, the interest rate.
        //
        //    Input, double SIGMA, the volatility of the asset.
        //
        //    Input, double T1, the expiry date.
        //
        //    Input, int M, the number of simulations.
        //
        //    Input/output, int &SEED, a seed for the random
        //    number generator.
        //
        //    Output, double MC[2], the estimated range of the valuation.
        //
    {
        int i;

        double[] u = typeMethods.r8vec_normal_01_new(m, ref data, ref seed);

        double[] svals = new double[m];

        for (i = 0; i < m; i++)
        {
            svals[i] = s0 * Math.Exp((r - 0.5 * sigma * sigma) * t1
                                     + sigma * Math.Sqrt(t1) * u[i]);
        }

        double[] pvals = new double[m];

        for (i = 0; i < m; i++)
        {
            pvals[i] = Math.Exp(-r * t1) * Math.Max(svals[i] - e, 0.0);
        }

        double pmean = 0.0;
        for (i = 0; i < m; i++)
        {
            pmean += pvals[i];
        }

        pmean /= m;

        double std = 0.0;
        for (i = 0; i < m; i++)
        {
            std += Math.Pow(pvals[i] - pmean, 2);
        }

        std = Math.Sqrt(std / (m - 1));

        double width = 1.96 * std / Math.Sqrt(m);

        double[] conf = new double[2];

        conf[0] = pmean - width;
        conf[1] = pmean + width;

        return conf;
    }
}