using System;
using Burkardt.Uniform;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public class r8NormalData
        {
            public int seed2 = 0;
            public int used = 0;
            public double y = 0.0;

        }
        public static double r8_normal_01(ref r8NormalData data, ref int seed)
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
            double x;
            //
            //  If we've used an even number of values so far, generate two more, 
            //  return one, and save one.
            //
            if ((data.used % 2) == 0)
            {
                double r1 = UniformRNG.r8_uniform_01(ref seed);
                data.seed2 = seed;
                double r2 = UniformRNG.r8_uniform_01(ref data.seed2);

                x = Math.Sqrt(-2.0 * Math.Log(r1)) * Math.Cos(2.0 * Math.PI * r2);
                data.y = Math.Sqrt(-2.0 * Math.Log(r1)) * Math.Sin(2.0 * Math.PI * r2);
            }
            else
            {
                seed = data.seed2;
                x = data.y;
            }

            data.used = data.used + 1;

            return x;
        }        
    }
}