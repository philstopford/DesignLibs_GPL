using System;

namespace Burkardt.Torus
{
    public static class Geometry
    {
        public static double torus_area_3d ( double r1, double r2 )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TORUS_AREA_3D returns the area of a torus in 3D.
            //
            //  Discussion:
            //
            //    A torus with radii R1 and R2 is the set of points satisfying:
            //
            //      ( sqrt ( X^2 + Y^2 ) - R1 )^2 + Z^2 <= R2^2
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 September 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double R1, R2, the two radii that define the torus.
            //
            //    Output, double TORUS_AREA_3D, the area of the torus.
            //
        {
            double area;

            area = 4.0 * Math.PI * Math.PI * r1 * r2;

            return area;
        }

        public static double torus_volume_3d ( double r1, double r2 )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TORUS_VOLUME_3D computes the volume of a torus in 3D.
            //
            //  Discussion:
            //
            //    A torus with radii R1 and R2 is the set of points satisfying:
            //
            //      ( sqrt ( X^2 + Y^2 ) - R1 )^2 + Z^2 <= R2^2
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    22 April 1999
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double R1, R2, the "inner" and "outer" radii of the torus.
            //
            //    Output, double TORUS_VOLUME_3D, the volume of the torus.
            //
        {
            double volume;

            volume = 2.0 * Math.PI * Math.PI * r1 * r2 * r2;

            return volume;
        }

    }
}