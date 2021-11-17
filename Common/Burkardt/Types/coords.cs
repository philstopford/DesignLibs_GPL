using System;

namespace Burkardt.Types;

public static partial class typeMethods
{
    public static double[] tp_to_xyz ( double theta, double phi )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TP_TO_XYZ converts unit spherical TP coordinates to XYZ coordinates.
        //
        //  Discussion:
        //
        //    The point is assume to lie on the unit sphere centered at the origin.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    22 September 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double THETA, PHI, the angular coordinates of a point
        //    on the unit sphere.
        //
        //    Output, double TP_TO_XYZ[3], the XYZ coordinates.
        //
    {
        double[] v;

        v = new double[3];

        v[0] = Math.Cos ( theta ) * Math.Sin ( phi );
        v[1] = Math.Sin ( theta ) * Math.Sin ( phi );
        v[2] =                 Math.Cos ( phi );

        return v;
    }
}