namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static double r8vec_dot(int n, double[] a1, double[] a2)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_DOT computes the dot product of a pair of R8VEC's.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    03 July 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the vectors.
            //
            //    Input, double A1[N], A2[N], the two vectors to be considered.
            //
            //    Output, double R8VEC_DOT, the dot product of the vectors.
            //
        {
            double value = 0.0;
            for (int i = 0; i < n; i++)
            {
                value = value + a1[i] * a2[i];
            }

            return value;
        }

        public static double r8vec_dot_product(int n, double[] a1, double[] a2)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_DOT_PRODUCT computes the dot product of a pair of R8VEC's.
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
            //    03 July 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the vectors.
            //
            //    Input, double A1[N], A2[N], the two vectors to be considered.
            //
            //    Output, double R8VEC_DOT_PRODUCT, the dot product of the vectors.
            //
        {
            double value = 0.0;
            for (int i = 0; i < n; i++)
            {
                value = value + a1[i] * a2[i];
            }

            return value;
        }
        
        public static double r8vec_i4vec_dot_product(int n, double[] r8vec, int[] i4vec)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_I4VEC_DOT_PRODUCT computes the dot product of an R8VEC and an I4VEC.
            //
            //  Discussion:
            //
            //    An R8VEC is a vector of R8's.
            //
            //    An I4VEC is a vector of I4's.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    30 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the vectors.
            //
            //    Input, double R8VEC[N], the first vector.
            //
            //    Input, int I4VEC[N], the second vector.
            //
            //    Output, double R8VEC_I4VEC_DOT_PRODUCT, the dot product of the vectors.
            //
        {
            int i;

            double value = 0.0;
            for (i = 0; i < n; i++)
            {
                value = value + r8vec[i] * (double) (i4vec[i]);
            }

            return value;
        }

    }
}