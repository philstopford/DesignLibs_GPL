using System;
using System.Linq;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static double r8vec_max(int n, double[] dvec)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_MAX returns the value of the maximum element in an R8VEC.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    15 October 1998
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the array.
            //
            //    Input, double *RVEC, a pointer to the first entry of the array.
            //
            //    Output, double R8VEC_MAX, the value of the maximum element.  This
            //    is set to 0.0 if N <= 0.
            //
        {
            if (dvec.Length <= 0)
            {
                return 0;
            }

            // Limit to the number of items in the array as a maximum
            n = Math.Min(n, dvec.Length);

            if (n == dvec.Length)
            {
                return dvec.Max();
            }

            return dvec.Take(n).Max();
        }


        public static double r8vec_min(int n, double[] dvec)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_MIN returns the value of the minimum element in an R8VEC.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    15 October 1998
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the array.
            //
            //    Input, double *RVEC, a pointer to the first entry of the array.
            //
            //    Output, double R8VEC_MIN, the value of the minimum element.  This
            //    is set to 0.0 if N <= 0.
            //
        {
            if (dvec.Length <= 0)
            {
                return 0;
            }

            // Limit to the number of items in the array as a maximum
            n = Math.Min(n, dvec.Length);

            if (n == dvec.Length)
            {
                return dvec.Min();
            }

            return dvec.Take(n).Min();
        }

        public static double r8vec_amax(int n, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_AMAX returns the maximum absolute value in an R8VEC.
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
            //    18 September 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the array.
            //
            //    Input, double A[N], the array.
            //
            //    Output, double AMAX, the value of the entry
            //    of largest magnitude.
            //
        {
            double amax;
            int i;

            amax = 0.0;
            for (i = 0; i < n; i++)
            {
                if (amax < Math.Abs(a[i]))
                {
                    amax = Math.Abs(a[i]);
                }
            }

            return amax;
        }

        public static int r8vec_amax_index(int n, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_AMAX_INDEX returns the index of the maximum absolute value in an R8VEC.
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
            //    20 September 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the array.
            //
            //    Input, double A[N], the array.
            //
            //    Output, int R8VEC_AMAX_INDEX, the index of the entry of largest magnitude.
            //
        {
            double amax;
            int amax_index;
            int i;

            if (n <= 0)
            {
                amax_index = -1;
            }
            else
            {
                amax_index = 1;
                amax = Math.Abs(a[0]);

                for (i = 2; i <= n; i++)
                {
                    if (amax < Math.Abs(a[i - 1]))
                    {
                        amax_index = i;
                        amax = Math.Abs(a[i - 1]);
                    }
                }
            }

            return amax_index;
        }

        public static double r8vec_amin(int n, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_AMIN returns the minimum absolute value in an R8VEC.
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
            //    18 September 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the array.
            //
            //    Input, double A[N], the array.
            //
            //    Output, double R8VEC_AMIN, the value of the entry
            //    of smallest magnitude.
            //
        {
            int i;
            const double r8_huge = 1.79769313486231571E+308;
            double value;

            value = r8_huge;
            for (i = 0; i < n; i++)
            {
                if (Math.Abs(a[i]) < value)
                {
                    value = Math.Abs(a[i]);
                }
            }

            return value;
        }

        public static int r8vec_amin_index(int n, double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8VEC_AMIN_INDEX returns the index of the minimum absolute value in an R8VEC.
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
            //    20 September 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of entries in the array.
            //
            //    Input, double A[N], the array.
            //
            //    Output, int R8VEC_AMIN_INDEX, the index of the entry of smallest magnitude.
            //
        {
            double amin;
            int amin_index;
            int i;

            if (n <= 0)
            {
                amin_index = -1;
            }
            else
            {
                amin_index = 1;
                amin = Math.Abs(a[0]);

                for (i = 2; i <= n; i++)
                {
                    if (Math.Abs(a[i - 1]) < amin)
                    {
                        amin_index = i;
                        amin = Math.Abs(a[i - 1]);
                    }
                }
            }

            return amin_index;
        }

    }
}