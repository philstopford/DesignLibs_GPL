using Burkardt.Sampling;
using Burkardt.Uniform;

namespace Burkardt.HyperGeometry.Hypercube;

public static class Sample
{
    public static GeometrySampleResult sample_hypercube_uniform ( int dim_num, int n, int seed )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SAMPLE_HYPERCUBE_UNIFORM returns sample points in the unit hypercube.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int N, the number of points to compute.
        //
        //    Input/output, int *SEED, a seed for the random number generator.
        //
        //    Output, double SAMPLE_HYPERCUBE_UNIFORM[DIM_NUM*N], the sample points.
        //
    {
        double[] x = UniformRNG.r8mat_uniform_01 ( dim_num, n, ref seed );

        GeometrySampleResult result = new()
        {
            seed = seed,
            result = x
        };

        return result;
    }
}