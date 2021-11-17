using System;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static double r8_sign(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_SIGN returns the sign of an R8.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the number whose sign is desired.
        //
        //    Output, double R8_SIGN, the sign of X.
        //
    {
        return x switch
        {
            < 0.0 => -1.0,
            _ => 1.0
        };
    }

    public static double r8_sign2(double x, double y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_SIGN2 returns the first argument with the sign of the second.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 January 2002
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, Y, the input arguments.
        //
        //    Output, double R8_SIGN2, is equal to the absolute value of X, and
        //    has the sign of Y.
        //
    {
        double value = y switch
        {
            >= 0.0 => Math.Abs(x),
            _ => -Math.Abs(x)
        };

        return value;
    }

    public static double r8_sign3(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_SIGN3 returns the three-way sign of an R8.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 September 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the number whose sign is desired.
        //
        //    Output, double R8_SIGN3, the sign of X.
        //
    {
        double value = x switch
        {
            < 0.0 => -1.0,
            0.0 => 0.0,
            _ => 1.0
        };

        return value;
    }

    public static char r8_sign_char(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_SIGN_CHAR returns a character indicating the sign of an R8.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, the number whose sign is desired.
        //
        //    Output, char R8_SIGN_CHAR, the sign of X, '-', '0' or '+'.
        //
    {
        char value = x switch
        {
            < 0.0 => '-',
            0.0 => '0',
            _ => '+'
        };

        return value;
    }

    public static bool r8_sign_match(double r1, double r2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_SIGN_MATCH is TRUE if two R8's are of the same sign.
        //
        //  Discussion:
        //
        //    This test could be coded numerically as
        //
        //      if ( 0 <= r1 * r2 ) then ...
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R1, R2, the values to check.
        //
        //    Output, bool R8_SIGN_MATCH, is TRUE if ( R1 <= 0 and R2 <= 0 )
        //    or ( 0 <= R1 and 0 <= R2 ).
        //
    {
        bool value = r1 <= 0.0 && r2 <= 0.0 || 0.0 <= r1 && 0.0 <= r2;

        return value;
    }

    public static bool r8_sign_match_strict(double r1, double r2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_SIGN_MATCH_STRICT is TRUE if two R8's are of the same strict sign.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 April 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R1, R2, the values to check.
        //
        //    Output, bool R8_SIGN_MATCH_STRICT, is TRUE if the signs match.
        //
    {
        bool value = r1 < 0.0 && r2 < 0.0 ||
                     r1 == 0.0 && r2 == 0.0 ||
                     0.0 < r1 && 0.0 < r2;

        return value;
    }

    public static bool r8_sign_opposite(double r1, double r2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_SIGN_OPPOSITE is TRUE if two R8's are not of the same sign.
        //
        //  Discussion:
        //
        //    This test could be coded numerically as
        //
        //      if ( r1 * r2 <= 0.0 ) ...
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 June 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R1, R2, the values to check.
        //
        //    Output, bool R8_SIGN_OPPOSITE, is TRUE if ( R1 <= 0 and 0 <= R2 )
        //    or ( R2 <= 0 and 0 <= R1 ).
        //
    {
        bool value = r1 <= 0.0 && 0.0 <= r2 || r2 <= 0.0 && 0.0 <= r1;

        return value;
    }

    public static bool r8_sign_opposite_strict(double r1, double r2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8_SIGN_OPPOSITE_STRICT is TRUE if two R8's are strictly of opposite sign.
        //
        //  Discussion:
        //
        //    This test could be coded numerically as
        //
        //      if ( r1 * r2 < 0.0 ) ...
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 June 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R1, R2, the values to check.
        //
        //    Output, bool R8_SIGN_OPPOSITE_STRICT, is TRUE if ( R1 < 0 and 0 < R2 )
        //    or ( R2 < 0 and 0 < R1 ).
        //
    {
        bool value = r1 < 0.0 && 0.0 < r2 || r2 < 0.0 && 0.0 < r1;

        return value;
    }

}