using Burkardt.IntegralNS;

namespace Burkardt.ODENS
{
    public static class HermitePolyChaosExpansion
    {
        public static void pce_ode_hermite(double ti, double tf, int nt, double ui, int np,
            double alpha_mu, double alpha_sigma, ref double[] t, ref double[] u )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PCE_ODE_HERMITE applies the polynomial chaos expansion to a scalar ODE.
        //
        //  Discussion:
        //
        //    The deterministic equation is
        //
        //      du/dt = - alpha * u,
        //      u(0) = u0
        //
        //    In the stochastic version, it is assumed that the decay coefficient
        //    ALPHA is a Gaussian random variable with mean value ALPHA_MU and variance
        //    ALPHA_SIGMA^2.
        //
        //    The exact expected value of the stochastic equation will be
        //
        //      u(t) = u0 * exp ( t^2/2)
        //
        //    This should be matched by the first component of the polynomial chaos
        //    expansion.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 March 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double TI, TF, the initial and final times.
        //
        //    Input, int  NT, the number of output points.
        //
        //    Input, double UI, the initial condition.
        //
        //    Input, int NP, the degree of the expansion.  Polynomials 
        //    of degree 0 through NP will be used.
        //
        //    Input, double ALPHA_MU, ALPHA_SIGMA, the mean and standard 
        //    deviation of the decay coefficient.
        //
        //    Output, double T[NT+1], U[(NT+1)*(NP+1)], the times and the PCE 
        //    coefficients at the successive time steps.
        //
        {
            double dp;
            double dt;
            int i;
            int it;
            int j;
            int k;
            double t1;
            double t2;
            double term;
            double tp;
            double[] u1;
            double[] u2;

            u1 = new double[np + 1];
            u2 = new double[np + 1];

            dt = (tf - ti) / (double)(nt);
            //
            //  Set the PCE coefficients for the initial time.
            //
            t1 = ti;

            u1[0] = ui;
            for (j = 1; j <= np; j++)
            {
                u1[j] = 0.0;
            }

            //
            //  Copy into the output arrays.
            //
            t[0] = t1;
            for (j = 0; j <= np; j++)
            {
                u[0 + j * (nt + 1)] = u1[j];
            }

            //
            //  Time integration.
            //
            for (it = 1; it <= nt; it++)
            {
                t2 = ((double)(nt - it) * ti
                      + (double)(it) * tf)
                     / (double)(nt);

                for (k = 0; k <= np; k++)
                {
                    dp = Integral.he_double_product_integral(k, k);

                    term = -alpha_mu * u1[k];

                    i = 1;
                    for (j = 0; j <= np; j++)
                    {
                        tp = Integral.he_triple_product_integral(i, j, k);
                        term = term - alpha_sigma * u1[j] * tp / dp;
                    }

                    u2[k] = u1[k] + dt * term;
                }

                //
                //  Prepare for next step.
                //
                t1 = t2;
                for (j = 0; j <= np; j++)
                {
                    u1[j] = u2[j];
                }

                //
                //  Copy into the output arrays.
                //
                t[it] = t1;
                for (j = 0; j <= np; j++)
                {
                    u[it + j * (nt + 1)] = u1[j];
                }
            }
        }

    }
}