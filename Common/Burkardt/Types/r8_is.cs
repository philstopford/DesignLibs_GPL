using System;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static bool r8_is_in_01(double a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_IN_01 is TRUE if an R8 is in the range [0,1].
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    13 June 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double A, the value.
            //
            //    Output, bool R8_IS_IN_01, is TRUE if A is between 0 and 1.
            //
        {
            bool value;

            value = (0.0 <= a && a <= 1.0);

            return value;
        }

        public static bool r8_is_inf(double r)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_IS_INF determines if an R8 represents an infinite value.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 May 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double R, the number to be checked.
            //
            //    Output, bool R8_IS_INF, is TRUE if R is an infinite value.
            //
        {
            const double r8_huge = 1.79769313486231571E+308;
            bool value;

            if (r < 0.0)
            {
                value = (r < -r8_huge);
            }
            else
            {
                value = (r8_huge < r);
            }

            return value;
        }

        public static bool r8_is_insignificant(double r, double s)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_IS_INSIGNIFICANT determines if an R8 is insignificant.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    26 November 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double R, the number to be compared against.
            //
            //    Input, double S, the number to be compared.
            //
            //    Output, bool R8_IS_INSIGNIFICANT, is TRUE if S is insignificant
            //    compared to R.
            //
        {
            double t;
            double tol;
            bool value;

            value = true;

            t = r + s;
            tol = typeMethods.r8_epsilon() * Math.Abs(r);

            if (tol < Math.Abs(r - t))
            {
                value = false;
            }

            return value;
        }

        public static bool r8_is_integer(double r)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_IS_INTEGER determines if an R8 represents an integer value.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    15 April 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double R, the number to be checked.
            //
            //    Output, bool R8_IS_INTEGER, is TRUE if R is an integer value.
            //
        {
            const int i4_huge = 2147483647;
            bool value;

            if ((double) (i4_huge) < r)
            {
                value = false;
            }
            else if (r < -(double) (i4_huge))
            {
                value = false;
            }
            else if (r == (double) ((int) (r)))
            {
                value = true;
            }
            else
            {
                value = false;
            }

            return value;
        }

        public static bool r8_is_nan(double r)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    R8_IS_NAN determines if an R8 represents a NaN value.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    14 May 2016
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double R, the number to be checked.
            //
            //    Output, bool R8_IS_NAN, is TRUE if R is a NaN
            //
        {
            bool value;

            value = (r != r);

            return value;
        }
        
    }
}