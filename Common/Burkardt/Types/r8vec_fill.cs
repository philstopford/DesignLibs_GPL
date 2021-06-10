namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static void r8vec_fill ( int n, double value, double[] x )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_FILL sets all entries of an R8VEC to a given value.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 December 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of elements of A.
            //
            //    Input, double VALUE, the value.
            //
            //    Output, double X[N], the array.
            //
        {
            int i;

            for ( i = 0; i < n; i++ )
            {
                x[i] = value;
            }
            return;
        }

        public static double[] r8vec_fill_new ( int n, double value )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_FILL_NEW creates an R8VEC, setting all entries to a given value.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 December 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of elements of A.
            //
            //    Input, double VALUE, the value.
            //
            //    Output, double R8VEC_FILL_NEW[N], the array.
            //
        {
            int i;
            double[] x;

            x = new double[n];

            for ( i = 0; i < n; i++ )
            {
                x[i] = value;
            }
            return x;
        }
    }
}