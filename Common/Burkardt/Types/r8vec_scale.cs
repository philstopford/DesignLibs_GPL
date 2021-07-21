namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static void r8vec_scale ( double s, int n, ref double[] a )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_SCALE multiplies an R8VEC by a scale factor.
            //
            //  Discussion:
            //
            //    An R8VEC is a vector of R8's.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    22 September 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double S, the scale factor.
            //
            //    Input, int N, the number of entries in the vectors.
            //
            //    Input/output, double A[N], the vector to be scaled.
            //    On output, A[] = S * A[].
            //
        {
            int i;

            for ( i = 0; i < n; i++ )
            {
                a[i] = s * a[i];
            }
            return;
        }
        
        public static double[] r8vec_scale_01 ( int n, double[] x )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_SCALE_01 scales an R8VEC to [0,1].
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 October 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the vector.
            //
            //    Input, double X[N], the vector to be scaled.
            //
            //    Output, double R8VEC_SCALE_01[N], the scaled vector.
            //
        {
            int i;
            double xmax;
            double xmin;
            double[] xs;

            xs = r8vec_copy_new ( n, x );

            xmin = r8vec_min ( n, xs );
            xmax = r8vec_max ( n, xs );
            if ( 0.0 < xmax - xmin )
            {
                for ( i = 0; i < n; i++ )
                {
                    xs[i] = ( xs[i] - xmin ) / ( xmax - xmin );
                }
            }
            else
            {
                for ( i = 0; i < n; i++ )
                {
                    xs[i] = 0.5;
                }
            }

            return xs;
        }

        public static double[] r8vec_scale_ab ( int n, double[] x, double a, double b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8VEC_SCALE_AB scales an R8VEC to [A,B].
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 October 2018
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the number of entries in the vector.
        //
        //    Input, double X[N], the vector to be scaled.
        //
        //    Input, double A, B, define the new range;
        //
        //    Output, double R8VEC_SCALE_01[N], the scaled vector.
        //
        {
            int i;
            double xmax;
            double xmin;
            double[] xs;

            xs = r8vec_copy_new ( n, x );

            xmin = r8vec_min ( n, xs );
            xmax = r8vec_max ( n, xs );
            if ( 0.0 < xmax - xmin )
            {
                for ( i = 0; i < n; i++ )
                {
                    xs[i] = a + ( b - a ) * ( xs[i] - xmin ) / ( xmax - xmin );
                }
            }
            else
            {
                for ( i = 0; i < n; i++ )
                {
                    xs[i] = ( a + b ) / 2.0;
                }
            }

            return xs;
        }

    }
}