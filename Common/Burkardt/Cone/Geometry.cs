using System;

namespace Burkardt.Cone;

public static class Geometry
{
    public static double cone_area_3d(double h, double r)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CONE_AREA_3D computes the surface area of a right circular cone in 3D.
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
        //    Input, double H, R, the height of the cone, and the radius of the
        //    circle that forms the base of the cone.
        //
        //    Output, double CONE_AREA_3D, the surface area of the cone.
        //
    {
        double area;

        area = Math.PI * r * Math.Sqrt(h * h + r * r);

        return area;
    }

    public static double[] cone_centroid_3d(double r, double[] pc, double[] pt)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CONE_CENTROID_3D returns the centroid of a cone in 3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Adrian Bowyer, John Woodwark,
        //    A Programmer's Geometry,
        //    Butterworths, 1983.
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the circle at the base of the cone.
        //
        //    Input, double PC[3], the coordinates of the center of the circle.
        //
        //    Input, double PT[3], the coordinates of the tip of the cone.
        //
        //    Output, double CONE_CENTROID_3D[3], the coordinates of the centroid of the cone.
        //
    {
        double[] centroid;
        int dim_num = 3;
        int i;

        centroid = new double[3];

        for (i = 0; i < dim_num; i++)
        {
            centroid[i] = 0.75 * pc[i] + 0.25 * pt[i];
        }

        return centroid;
    }

    public static double cone_volume_3d(double h, double r)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CONE_VOLUME_3D computes the volume of a right circular cone in 3D.
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
        //    Input, double H, R, the height of the cone, and the radius of the
        //    circle that forms the base of the cone.
        //
        //    Output, double CONE_VOLUME_3D, the volume of the cone.
        //
    {
        double volume;

        volume = Math.PI * r * r * h / 3.0;

        return volume;
    }
}