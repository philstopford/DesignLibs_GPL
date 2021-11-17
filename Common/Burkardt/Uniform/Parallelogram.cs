namespace Burkardt.Uniform;

public static class Parallelogram
{
    public static double[] uniform_in_parallelogram_map ( double[] v1, double[] v2,
            double[] v3, int n, ref int seed )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UNIFORM_IN_PARALLELOGRAM_MAP maps uniform points into a parallelogram.
        //
        //  Discussion:
        //
        //    The parallelogram is defined by three vertices, V1, V2 and V3.
        //    The missing vertex V4 is equal to V2+V3-V1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 August 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Greg Turk,
        //    Generating Random Points in a Triangle,
        //    in Graphics Gems,
        //    edited by Andrew Glassner,
        //    AP Professional, 1990, pages 24-28.
        //
        //  Parameters:
        //
        //    Input, double V1[2], V2[2], V3[2], the vertices.
        //
        //    Input, int N, the number of points.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double UNIFORM_IN_PARALLELOGRAM_MAP[2*N], the points.
        //
    {
        const int DIM_NUM = 2;

        int j;

        double[] x = new double[DIM_NUM*n];

        for ( j = 0; j < n; j++ )
        {
            double r = UniformRNG.r8_uniform_01 ( ref seed );
            double s = UniformRNG.r8_uniform_01 ( ref seed );

            int i;
            for ( i = 0; i < DIM_NUM; i++ )
            {
                x[i+j*DIM_NUM] = ( 1.0 - r - s ) * v1[i]
                                 + r       * v2[i]
                                 + s   * v3[i];
            }
        }

        return x;
    }
}