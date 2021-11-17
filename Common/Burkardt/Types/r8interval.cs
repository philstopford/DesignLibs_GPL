namespace Burkardt.Types;

public static partial class typeMethods
{
    public static double r8int_to_r8int(double rmin, double rmax, double r, double r2min,
            double r2max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8INT_TO_R8INT maps one R8 interval to another.
        //
        //  Discussion:
        //
        //    The formula used is
        //
        //      R2 := R2MIN + ( R2MAX - R2MIN ) * ( R - RMIN ) / ( RMAX - RMIN )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 January 2001
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double RMIN, RMAX, the first range.
        //
        //    Input, double R, the number to be converted.
        //
        //    Input, double R2MAX, R2MIN, the second range.
        //
        //    Output, double R8INT_TO_R8INT, the corresponding value in
        //    the range [R2MIN,R2MAX].
        //
    {
        double r2;

        if (rmax == rmin)
        {
            r2 = (r2max + r2min) / 2.0;
        }
        else
        {
            r2 = ((rmax - r) * r2min
                  + (r - rmin) * r2max)
                 / (rmax - rmin);
        }

        return r2;
    }

    public static int r8int_to_i4int(double rmin, double rmax, double r, int imin, int imax)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8INT_TO_I4INT maps an R8 interval to an integer interval.
        //
        //  Discussion:
        //
        //    The formula used is
        //
        //      I := IMIN + ( IMAX - IMIN ) * ( R - RMIN ) / ( RMAX - RMIN )
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 January 2001
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double RMIN, RMAX, the range.
        //
        //    Input, double R, the number to be converted.
        //
        //    Input, int IMAX, IMIN, the integer range.
        //
        //    Output, int R8INT_TO_I4INT, the corresponding value in the range [IMIN,IMAX].
        //
    {
        int i;

        if (rmax == rmin)
        {
            i = (imax + imin) / 2;
        }
        else
        {
            i = (int) r8_nint(
                ((rmax - r) * imin
                 + (r - rmin) * imax)
                / (rmax - rmin));
        }

        return i;
    }
}