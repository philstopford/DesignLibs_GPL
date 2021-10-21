using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.Ball
{
    public static class Geometry
    {
        public static double[] ball01_sample_2d(ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BALL01_SAMPLE_2D picks a random point in the unit ball in 2D.
            //
            //  Discussion:
            //
            //    The unit ball is the set of points such that
            //
            //      X * X + Y * Y <= 1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    25 August 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input/output, int &SEED, a seed for the random number generator.
            //
            //    Output, double BALL01_SAMPLE_2D[2], a random point in the unit ball.
            //
        {
            const double r8_pi = 3.141592653589793;
            double r;
            double theta;
            double[] x;

            r = UniformRNG.r8_uniform_01(ref seed);
            r = Math.Sqrt(r);

            theta = UniformRNG.r8_uniform_01(ref seed);
            theta = 2.0 * r8_pi * theta;

            x = new double[2];

            x[0] = r * Math.Cos(theta);
            x[1] = r * Math.Sin(theta);

            return x;
        }

        public static double[] ball01_sample_3d(ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BALL01_SAMPLE_3D picks a random point in the unit ball in 3D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    25 August 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input/output, int &SEED, a seed for the random number generator.
            //
            //    Output, double BALL01_SAMPLE_3D[3], the sample point.
            //
        {
            int DIM_NUM = 3;

            double phi;
            const double r8_pi = 3.141592653589793;
            double r;
            double theta;
            double vdot;
            double[] x;
            //
            //  Pick a uniformly random VDOT, which must be between -1 and 1.
            //  This represents the dot product of the random vector with the Z unit vector.
            //
            //  This works because the surface area of the sphere between
            //  Z and Z + dZ is independent of Z.  So choosing Z uniformly chooses
            //  a patch of area uniformly.
            //
            vdot = UniformRNG.r8_uniform_01(ref seed);
            vdot = 2.0 * vdot - 1.0;

            phi = typeMethods.r8_acos(vdot);
            //
            //  Pick a uniformly random rotation between 0 and 2 Pi around the
            //  axis of the Z vector.
            //
            theta = UniformRNG.r8_uniform_01(ref seed);
            theta = 2.0 * r8_pi * theta;
            //
            //  Pick a random radius R.
            //
            r = UniformRNG.r8_uniform_01(ref seed);
            r = Math.Pow((double) r, (double) (1.0 / 3.0));

            x = new double[DIM_NUM];

            x[0] = r * Math.Cos(theta) * Math.Sin(phi);
            x[1] = r * Math.Sin(theta) * Math.Sin(phi);
            x[2] = r * Math.Cos(phi);

            return x;
        }

        public static double[] ball01_sample_nd(int dim_num, ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BALL01_SAMPLE_ND picks a random point in the unit ball in ND.
            //
            //  Discussion:
            //
            //    DIM_NUM-1 random Givens rotations are applied to the point ( 1, 0, 0, ..., 0 ).
            //
            //    The I-th Givens rotation is in the plane of coordinate axes I and I+1,
            //    and has the form:
            //
            //     [ cos ( theta )  - sin ( theta ) ] * x(i)      = x'(i)
            //     [ sin ( theta )    cos ( theta ) ]   x(i+1)      x'(i+1)
            //
            //    Finally, a scaling is applied to set the point at a distance R
            //    from the center, in a way that results in a uniform distribution.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    25 August 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the dimension of the space.
            //
            //    Input/output, int &SEED, a seed for the random number generator.
            //
            //    Output, double BALL01_SAMPLE_ND[DIM_NUM], the random point.
            //
        {
            int i;
            double r;
            double random_cosine;
            double random_sign;
            double random_sine;
            double[] x;
            double xi;

            x = new double[dim_num];

            x[0] = 1.0;
            for (i = 1; i < dim_num; i++)
            {
                x[i] = 0.0;
            }

            for (i = 0; i < dim_num - 1; i++)
            {
                random_cosine = UniformRNG.r8_uniform_01(ref seed);
                random_cosine = 2.0 * random_cosine - 1.0;
                random_sign = UniformRNG.r8_uniform_01(ref seed);
                random_sign = (double)
                    (2 * (int) (2.0 * random_sign - 1.0));
                random_sine = random_sign
                              * Math.Sqrt(1.0 - random_cosine * random_cosine);
                xi = x[i];
                x[i] = random_cosine * xi;
                x[i + 1] = random_sine * xi;
            }

            r = UniformRNG.r8_uniform_01(ref seed);

            r = Math.Pow((double) r, 1.0 / (double) dim_num);

            for (i = 0; i < dim_num; i++)
            {
                x[i] = r * x[i];
            }

            return x;
        }

        public static double ball01_volume()

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    BALL01_VOLUME returns the volume of the unit ball in 3D.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    17 January 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Output, double BALL01_VOLUME, the volume of the unit ball.
            //
        {
            double r;
            const double r8_pi = 3.141592653589793;
            double volume;

            r = 1.0;
            volume = 4.0 * r8_pi * r * r * r / 3.0;

            return volume;
        }
    }
}