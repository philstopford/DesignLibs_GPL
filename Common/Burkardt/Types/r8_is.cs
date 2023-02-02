using System;

namespace Burkardt.Types;

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
        bool value = a is >= 0.0 and <= 1.0;

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

        bool value = r switch
        {
            < 0.0 => r < -r8_huge,
            _ => r8_huge < r
        };

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
        bool value = true;

        double t = r + s;
        double tol = r8_epsilon() * Math.Abs(r);

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

        switch (r)
        {
            case > i4_huge:
            case < -(double) i4_huge:
                value = false;
                break;
            default:
            {
                value = Math.Abs(r - (int) r) switch
                {
                    <= typeMethods._r8_epsilon => true,
                    _ => false
                };

                break;
            }
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
        return double.IsNaN(r);
    }
        
}