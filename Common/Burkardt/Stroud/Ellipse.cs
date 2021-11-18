using System;
using Burkardt.Types;

namespace Burkardt.Stroud;

public static class Ellipse
{
    public static double ellipse_area_2d(double r1, double r2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ELLIPSE_AREA_2D returns the area of an ellipse in 2D.
        //
        //  Integration region:
        //
        //    ( ( X - CENTER(1) ) / R1 )^2 + ( ( Y - CENTER(2) ) / R2 )^2 <= 1
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 March 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R1, R2, the major and minor semi-axes.
        //
        //    Output, double ELLIPSE_AREA_2D, the area of the ellipse.
        //
    {
        double value = Math.PI * r1 * r2;

        return value;
    }

    public static double ellipse_circumference_2d(double r1, double r2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ELLIPSE_CIRCUMFERENCE_2D returns the circumference of an ellipse in 2D.
        //
        //  Discussion:
        //
        //    There is no closed formula for the circumference of an ellipse.
        //
        //    Defining the eccentricity by
        //
        //      E = Math.Sqrt ( 1 - ( r2 / r1 )^2 )
        //
        //    where R1 and R2 are the major and minor axes, then
        //
        //      circumference
        //        = 4 * R1 * E(K,2*PI)
        //        = R1 * Integral ( 0 <= T <= 2*PI ) Math.Sqrt ( 1 - E^2 * sin^2 ( T ) ) dT
        //
        //    This integral can be approximated by the Gauss-Kummer formula.
        //
        //  Integration region:
        //
        //    ( ( X - CENTER(1) ) / R1 )^2 + ( ( Y - CENTER(2) ) / R2 )^2 <= 1
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 March 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    John Harris, Horst Stocker,
        //    Handbook of Mathematics and Computational Science,
        //    Springer, 1998,
        //    ISBN: 0-387-94746-9,
        //    LC: QA40.S76.
        //
        //  Parameters:
        //
        //    Input, double R1, R2, the major and minor semi-axes.
        //
        //    Output, double ELLIPSE_CIRCUMFERENCE_2D, the
        //    circumference of the ellipse.
        //
    {
        double value;

        if (Math.Abs(r1 - r2) <= double.Epsilon)
        {
            value = 2.0 * Math.PI * r1;
            return value;
        }

        //
        //  Compute the eccentricity of the ellipse.
        //
        double e = Math.Sqrt(1.0 - Math.Pow(Math.Min(r1, r2) / Math.Max(r1, r2), 2));

        value = 1.0;
        double term = value;
        int i = 0;

        for (;;)
        {
            i += 1;
            term = term * (2 * i - 3) * (2 * i - 1) * e * e
                   / (2 * 2 * i * i);

            if (Math.Abs(term) <= typeMethods.r8_epsilon() * (Math.Abs(value) + 1.0))
            {
                break;
            }

            value += term;
        }

        value = 2.0 * Math.PI * Math.Max(r1, r2) * value;

        return value;
    }

    public static double ellipse_eccentricity_2d(double r1, double r2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ELLIPSE_ECCENTRICITY_2D returns the eccentricity of an ellipse in 2D.
        //
        //  Integration region:
        //
        //    ( ( X - CENTER(1) ) / R1 )^2 + ( ( Y - CENTER(2) ) / R2 )^2 <= 1
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 March 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R1, R2, the major and minor semi-axes.
        //
        //    Output, double ELLIPSE_ECCENTRICITY_2D, the eccentricity 
        //    of the ellipse.
        //
    {
        double value;

        double minor = Math.Min(Math.Abs(r1), Math.Abs(r2));
        double major = Math.Max(Math.Abs(r1), Math.Abs(r2));

        switch (major)
        {
            case 0.0:
                value = -typeMethods.r8_huge();
                return value;
            default:
                value = Math.Sqrt(1.0 - Math.Pow(minor / major, 2));

                return value;
        }
    }

    public static double ellipsoid_volume_3d(double r1, double r2, double r3)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ELLIPSOID_VOLUME_3D returns the volume of an ellipsoid in 3d.
        //
        //  Discussion:
        //
        //    This is not a general ellipsoid, but one for which each of the 
        //    axes lies along a coordinate axis.
        //
        //  Integration region:
        //
        //      ( ( X - CENTER(1) ) / R1 )^2 
        //    + ( ( Y - CENTER(2) ) / R2 )^2
        //    + ( ( Z - CENTER(3) ) / R3 )^2 <= 1
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 March 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R1, R2, R3, the semi-axes of the ellipsoid.
        //
        //    Output, double ELLIPSOID_VOLUME_3D, the volume of the ellipsoid.
        //
    {
        double value = 4.0 / 3.0 * Math.PI * r1 * r2 * r3;

        return value;
    }

}