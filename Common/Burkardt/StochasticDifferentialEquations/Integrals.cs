﻿using System;
using Burkardt.Types;

namespace Burkardt.StochasticDifferentialEquations;

public static class Integrals
{
    public static void stochastic_integral_ito(int n, ref typeMethods.r8vecNormalData data, ref int seed, ref double estimate,
            ref double exact, ref double error)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    STOCHASTIC_INTEGRAL_ITO approximates the Ito integral of W(t) dW.
        //
        //  Discussion:
        //
        //    This function estimates the Ito integral of W(t) dW over 
        //    the interval [0,1].
        //
        //    The estimates is made by taking N steps.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 September 2012
        //
        //  Author:
        //
        //    Original Matlab version by Desmond Higham.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Desmond Higham,
        //    An Algorithmic Introduction to Numerical Simulation of
        //    Stochastic Differential Equations,
        //    SIAM Review,
        //    Volume 43, Number 3, September 2001, pages 525-546.
        //
        //  Parameters:
        //
        //    Input, int N, the number of steps to take.
        //
        //    Input, int &SEED, a seed for the random number generator.
        //
        //    Output, double &ESTIMATE, the estimate of the integral.
        //
        //    Output, double &EXACT, the exact value of the integral.
        //
        //    Output, double &ERROR, the error in the integral estimate.
        //
    {
        int j;
        //
        //  Set step parameters.
        //
        const double tmax = 1.0;
        double dt = tmax / n;
        //
        //  Define the increments dW.
        //
        double[] dw = typeMethods.r8vec_normal_01_new(n, ref data, ref seed);
        for (j = 0; j < n; j++)
        {
            dw[j] = Math.Sqrt(dt) * dw[j];
        }

        //
        //  Sum the increments to get the Brownian path.
        //
        double[] w = new double[n + 1];
        w[0] = 0.0;
        for (j = 1; j <= n; j++)
        {
            w[j] = w[j - 1] + dw[j - 1];
        }

        //
        //  Approximate the Ito integral.
        //
        estimate = typeMethods.r8vec_dot_product(n, w, dw);
        //
        //  Compare with the exact solution.
        //
        exact = 0.5 * (w[n] * w[n] - tmax);
        error = Math.Abs(estimate - exact);
    }

    public static void stochastic_integral_strat(int n, ref typeMethods.r8vecNormalData data, ref int seed, ref double estimate,
            ref double exact, ref double error)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    STOCHASTIC_INTEGRAL_STRAT approximates the Stratonovich integral of W(t) dW.
        //
        //  Discussion:
        //
        //    This function estimates the Stratonovich integral of W(t) dW over 
        //    the interval [0,1].
        //
        //    The estimates is made by taking N steps.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 September 2012
        //
        //  Author:
        //
        //    Original Matlab version by Desmond Higham.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Desmond Higham,
        //    An Algorithmic Introduction to Numerical Simulation of
        //    Stochastic Differential Equations,
        //    SIAM Review,
        //    Volume 43, Number 3, September 2001, pages 525-546.
        //
        //  Parameters:
        //
        //    Input, int N, the number of steps to take.
        //
        //    Input, int &SEED, a seed for the random number generator.
        //
        //    Output, double &ESTIMATE, the estimate of the integral.
        //
        //    Output, double &EXACT, the exact value of the integral.
        //
        //    Output, double &ERROR, the error in the integral estimate.
        //
    {
        int j;
        //
        //  Set step parameters.
        //
        const double tmax = 1.0;
        double dt = tmax / n;
        //
        //  Define the increments dW.
        //
        double[] dw = typeMethods.r8vec_normal_01_new(n, ref data, ref seed);

        for (j = 0; j < n; j++)
        {
            dw[j] = Math.Sqrt(dt) * dw[j];
        }

        //
        //  Sum the increments to get the Brownian path.
        //
        double[] w = new double[n + 1];
        w[0] = 0.0;
        for (j = 1; j <= n; j++)
        {
            w[j] = w[j - 1] + dw[j - 1];
        }

        //
        //  Approximate the Stratonovich integral.
        //
        double[] u = typeMethods.r8vec_normal_01_new(n, ref data, ref seed);

        double[] v = new double[n];
        for (j = 0; j < n; j++)
        {
            v[j] = 0.5 * (w[j] + w[j + 1]) + 0.5 * Math.Sqrt(dt) * u[j];
        }

        estimate = typeMethods.r8vec_dot_product(n, v, dw);
        //
        //  Compare with the exact solution.
        //
        exact = 0.5 * w[n - 1] * w[n - 1];
        error = Math.Abs(estimate - exact);
    }
}