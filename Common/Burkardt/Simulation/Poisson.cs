﻿using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.Simulation
{
    public static class Poisson
    {
        public static void poisson_fixed_events(double lambda, int event_num, ref int seed,
                ref double[] t, ref double[] w)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POISSON_FIXED_EVENTS waits for a given number of Poisson events.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    27 September 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double LAMBDA, the average number of events per 
            //    unit time.
            //
            //    Input, int EVENT_NUM, the number of events to wait for.
            //
            //    Input/output, int &SEED, a seed for the random
            //    number generator.
            //
            //    Output, double T[EVENT_NUM+1], the time at which a total 
            //    of 0, 1, 2, ... and EVENT_NUM events were observed.
            //
            //    Output, double W[EVENT_NUM+1], the waiting time until the
            //    I-th event occurred.
            //
        {
            int i;
            double[] u;
            //
            //  Poisson waiting times follow an exponential distribution.
            //
            w[0] = 0.0;
            u = UniformRNG.r8vec_uniform_01_new(event_num, ref seed);
            for (i = 1; i <= event_num; i++)
            {
                w[i] = -Math.Log(u[i - 1]) / lambda;
            }

            //
            //  The time til event I is the sum of the waiting times 0 through I.
            //
            typeMethods.r8vec_cum(event_num + 1, w, ref t);
        }

        public static int poisson_fixed_time(double lambda, double time, ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POISSON_FIXED_TIME counts the Poisson events in a fied time.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    27 September 2012
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double LAMBDA, the average number of events 
            //    per unit time.
            //
            //    Input, double TIME, the amount of time to observe.
            //
            //    Input/output, int &SEED, a seed for the random
            //    number generator.
            //
            //    Output, int POISSON_FIXED_TIME, the number of Poisson events observed.
            //
        {
            double dt;
            int n;
            double t;
            double u;

            n = 0;
            t = 0.0;

            while (t < time)
            {
                u = UniformRNG.r8_uniform_01(ref seed);
                dt = -Math.Log(u) / lambda;
                n = n + 1;
                t = t + dt;
            }

            return n;
        }
    }
}