using System;
using Burkardt.Uniform;

namespace Burkardt.Probability
{
    public static class Disk
    {
        public static double[] disk_mean(double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DISK_MEAN returns the mean of disk points.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 March 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, C, the parameters of the disk.
        //    The disk is centered at (A,B) and has radius C.
        //    0.0 < C.
        //
        //    Output, double DISK_MEAN[2], the mean of the points.
        //
        {
            double[] x = new double[2];

            x[0] = a;
            x[1] = b;

            return x;
        }

        public static double[] disk_sample(double a, double b, double c, ref int seed)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DISK_SAMPLE samples points in a disk.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 March 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, C, the parameters of the disk.
        //    The disk is centered at (A,B) and has radius C.
        //    0.0 < C.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double DISK_SAMPLE[2], a sampled point of the disk.
        //
        {
            const double r8_pi = 3.14159265358979323;

            double[] x = new double[2];

            double radius_frac = Math.Sqrt(UniformRNG.r8_uniform_01(ref seed));

            double angle = 2.0 * r8_pi * UniformRNG.r8_uniform_01(ref seed);

            x[0] = a + c * radius_frac * Math.Cos(angle);
            x[1] = b + c * radius_frac * Math.Sin(angle);

            return x;
        }

        public static double disk_variance(double a, double b, double c)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DISK_VARIANCE returns the variance of points in a disk.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 March 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, B, C, the parameters of the disk.
        //    The disk is centered at (A,B) and has radius C.
        //    0.0 < C.
        //
        //    Output, double DISK_VARIANCE, the variance of points in the disk.
        //
        {
            double value = 0.5 * c * c;

            return value;
        }
    }
}