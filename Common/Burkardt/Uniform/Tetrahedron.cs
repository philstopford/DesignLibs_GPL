namespace Burkardt.Uniform
{
    public static class Tetrahedron
    {
        public static double[] uniform_in_tetrahedron ( double[] v, int n, ref int seed )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UNIFORM_IN_TETRAHEDRON returns uniform points in a tetrahedron.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 August 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Claudio Rocchini, Paolo Cignoni,
        //    Generating Random Points in a Tetrahedron,
        //    Journal of Graphics Tools,
        //    Volume 5, Number 5, 2000, pages 9-12.
        //
        //  Parameters:
        //
        //    Input, double V[3*4], the vertices of the tetrahedron.
        //
        //    Input, int N, the number of points.
        //
        //    Input/output, int &SEED, a seed for the random
        //    number generator.
        //
        //    Output, double UNIFORM_IN_TETRAHEDRON[3*N], the points.
        //
        {
            double[] c = new double[4];
            int i;
            int j;
            int k;
            double t;
            double[] x;

            x = new double[3*n];

            for ( j = 0; j < n; j++ )
            {
                UniformRNG.r8vec_uniform_01 ( 3, ref seed, ref c );

                if ( 1.0 < c[0] + c[1] )
                {
                    c[0] = 1.0 - c[0];
                    c[1] = 1.0 - c[1];
                }

                if ( 1.0 < c[1] + c[2] )
                {
                    t = c[2];
                    c[2] = 1.0 - c[0] - c[1];
                    c[1] = 1.0 - t;
                }
                else if ( 1.0 < c[0] + c[1] + c[2] )
                {
                    t = c[2];
                    c[2] = c[0] + c[1] + c[2] - 1.0;
                    c[0] = 1.0 - c[1] - t;
                }
                c[3] = 1.0 - c[0] - c[1] - c[2];

                for ( i = 0; i < 3; i++ )
                {
                    x[i+j*3] = 0.0;
                    for ( k = 0; k < 4; k++ )
                    {
                        x[i+j*3] = x[i+j*3] + v[i+k*3] * c[k];
                    }
                }
            }

            return x;
        }
    }
}