using System;
using Burkardt.Uniform;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static double r8_normal_01(ref int seed)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_NORMAL_01 samples the standard normal probability distribution.
            //
            //  Discussion:
            //
            //    The standard normal probability distribution function (PDF) has
            //    mean 0 and standard deviation 1.
            //
            //    The Box-Muller method is used, which is efficient, but
            //    generates two values at a time.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 June 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input/output, int &SEED, a seed for the random number generator.
            //
            //    Output, double R8_NORMAL_01, a normally distributed random value.
            //
        {
            double pi = 3.141592653589793;
            int seed2 = 0;
            int used = 0;
            double x;
            double y = 0.0;
            //
            //  If we've used an even number of values so far, generate two more, 
            //  return one, and save one.
            //
            if ((used % 2) == 0)
            {
                double r1 = UniformRNG.r8_uniform_01(ref seed);
                seed2 = seed;
                double r2 = UniformRNG.r8_uniform_01(ref seed2);

                x = Math.Sqrt(-2.0 * Math.Log(r1)) * Math.Cos(2.0 * pi * r2);
                y = Math.Sqrt(-2.0 * Math.Log(r1)) * Math.Sin(2.0 * pi * r2);
            }
            else
            {
                seed = seed2;
                x = y;
            }

            used = used + 1;

            return x;
        }        
    }
}