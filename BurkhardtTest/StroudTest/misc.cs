namespace StroudTest
{
    public static class misc
    {
        public static double[] setsim ( int n )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SETSIM defines a unit simplex.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 April 2008
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the spatial dimension.
            //
            //    Output, double V[N*(N+1)], the coordinates of the N+1 vertices.
            //
        {
            int i;
            int j;
            double[] v;

            v = new double[n*(n+1)];

            for ( i = 0; i < n; i++ )
            {
                for ( j = 0; j < n + 1; j++ )
                {
                    v[i+j*n] = 0.0;
                }
            }

            for ( i = 0; i < n; i++ )
            {
                j = i + 1;
                v[i+j*n] = 1.0;
            }

            return v;
        }

    }
}