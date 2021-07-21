namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static void r8vec_stutter ( int n, double[] a, int m, ref double[] am )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_STUTTER makes a "stuttering" copy of an R8VEC.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    Applying a stuttering factor M of 3, the vector A = ( 1, 5, 8 ) becomes
        //    AM = ( 1, 1, 1, 5, 5, 5, 8, 8, 8 ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 March 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the size of the input vector.
        //
        //    Input, double A[N], the vector.
        //
        //    Input, int M, the "stuttering factor".
        //
        //    Output, double AM[M*N], the stuttering vector.
        //
        {
            int i;
            int j;
            int k;

            k = 0;

            for ( i = 0; i < n; i++ )
            {
                for ( j = 0; j < m; j++ )
                {
                    am[k] = a[i];
                    k = k + 1;
                }
            }
            return;
        }

        public static double[] r8vec_stutter_new ( int n, double[] a, int m )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_STUTTER_NEW makes a "stuttering" copy of an R8VEC.
        //
        //  Discussion:
        //
        //    An R8VEC is a vector of R8's.
        //
        //    Applying a stuttering factor M of 3, the vector A = ( 1, 5, 8 ) becomes
        //    AM = ( 1, 1, 1, 5, 5, 5, 8, 8, 8 ).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 March 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the size of the input vector.
        //
        //    Input, double A[N], the vector.
        //
        //    Input, int M, the "stuttering factor".
        //
        //    Output, double R8VEC_STUTTER_NEW[M*N], the stuttering vector.
        //
        {
            double[] am;
            int i;
            int j;
            int k;

            am = new double[m*n];

            k = 0;

            for ( i = 0; i < n; i++ )
            {
                for ( j = 0; j < m; j++ )
                {
                    am[k] = a[i];
                    k = k + 1;
                }
            }
            return am;
        }

    }
}