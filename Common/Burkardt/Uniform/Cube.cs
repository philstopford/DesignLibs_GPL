using Burkardt.Types;

namespace Burkardt.Uniform
{
    public static class Cube
    {
        public static double[] uniform_in_cube01(int dim_num, int n, ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    UNIFORM_IN_CUBE01 creates uniform points in the unit hypercube.
            //
            //  Discussion:
            //
            //    The unit hypercube is defined as points all of whose components are between
            //    0 and 1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    17 August 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the dimension of the space.
            //
            //    Input, int N, the number of points.
            //
            //    Input/output, int &SEED, a seed for the random number generator.
            //
            //    Output, double UNIFORM_IN_CUBE01[DIM_NUM*N], the points.
            //
        {
            double[] x;

            x = new double[dim_num * n];

            UniformRNG.r8vec_uniform_01(dim_num * n, ref seed, ref x);

            return x;
        }

        public static double[] uniform_on_cube(int m, int n, double[] c, double r, ref int seed )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UNIFORM_ON_CUBE returns random points on the surface of a cube.
        //
        //  Discussion:
        //
        //    The cube is assumed to be aligned with the coordinate axes.
        //
        //    The cube has center C and radius R.  Any point on the surface of
        //    the cube is described by
        //
        //      X = C + R * PM
        //
        //    where PM is an M-dimensional vector whose entries are between
        //    -1 and +1, and for which at least one value has norm 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 April 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the spatial dimension.
        //    1 <= M.
        //
        //    Input, int N, the number of points.
        //
        //    Input, double C[M], the coordinates of the center.
        //
        //    Input, double R, the radius.
        //
        //    Input/output, int &SEED, a seed for the random
        //    number generator.
        //
        //    Output, double UNIFORM_ON_CUBE[M*N], the coordinates of N points, chosen
        //    uniformly at random from the surface of the M-cube of center C and 
        //    radius R.
        //
        {
            int i;
            int j;
            int k;
            double[] x;
            //
            //  Choose random points within the cube of radius 1.
            //
            x = UniformRNG.r8mat_uniform_01_new(m, n, ref seed);

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    x[i + j * m] = 2.0 * x[i + j * m] - 1.0;
                }
            }

            //
            //  For each point, select a coordinate at random, and set it to +1 or -1.
            //
            for (j = 0; j < n; j++)
            {
                i = UniformRNG.i4_uniform_ab(0, m - 1, ref seed);
                k = UniformRNG.i4_uniform_ab(0, 1, ref seed);
                if (k == 0)
                {
                    x[i + j * m] = 1.0;
                }
                else
                {
                    x[i + j * m] = -1.0;
                }
            }

            //
            //  Shift by C and scale by R.
            //
            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    x[i + j * m] = c[i] + r * x[i + j * m];
                }
            }

            return x;
        }

        public static double[] uniform_on_cube01(int m, int n, ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    UNIFORM_ON_CUBE01 returns random points on the surface of the unit cube.
            //
            //  Discussion:
            //
            //    The cube is assumed to be aligned with the coordinate axes.
            //
            //    The cube has center at the origin and radius 1. Any point on the surface
            //    of the cube is described by
            //
            //      X = PM
            //
            //    where PM is an M-dimensional vector whose entries are between
            //    -1 and +1, and for which at least one value has norm 1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    20 April 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the spatial dimension.
            //    1 <= M.
            //
            //    Input, int N, the number of points.
            //
            //    Input/output, int &SEED, a seed for the random
            //    number generator.
            //
            //    Output, double UNIFORM_ON_CUBE[M*N], the coordinates of N points, chosen
            //    uniformly at random from the surface of the unit M-cube.
            //
        {
            int i;
            int j;
            int k;
            double[] x;
            //
            //  Choose random points within the cube of radius 1.
            //
            x = UniformRNG.r8mat_uniform_01_new(m, n, ref seed);

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < m; i++)
                {
                    x[i + j * m] = 2.0 * x[i + j * m] - 1.0;
                }
            }

            //
            //  For each point, select a coordinate at random, and set it to +1 or -1.
            //
            for (j = 0; j < n; j++)
            {
                i = UniformRNG.i4_uniform_ab(0, m - 1, ref seed);
                k = UniformRNG.i4_uniform_ab(0, 1, ref seed);
                if (k == 0)
                {
                    x[i + j * m] = 1.0;
                }
                else
                {
                    x[i + j * m] = -1.0;
                }
            }

            return x;
        }

    }
}