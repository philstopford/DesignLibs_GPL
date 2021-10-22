using System;

namespace Burkardt.Geometry
{
    public static class XY
    {
        public static void rtp_to_xyz ( double r, double theta, double phi, ref double[] xyz )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RTP_TO_XYZ converts (R,Theta,Phi) to (X,Y,Z) coordinates.
        //
        //  Discussion:
        //
        //    R measures the distance of the point to the origin.
        //
        //    Theta measures the "longitude" of the point, between 0 and 2 PI.
        //
        //    PHI measures the angle from the "north pole", between 0 and PI.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 July 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, THETA, PHI, the radius, longitude, and
        //    declination of a point.
        //
        //    Output, double XYZ[3], the corresponding Cartesian coordinates.
        //
        {
            xyz[0] = r * Math.Cos ( theta ) * Math.Sin ( phi );
            xyz[1] = r * Math.Sin ( theta ) * Math.Sin ( phi );
            xyz[2] = r *                 Math.Cos ( phi );

        }
        
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
}