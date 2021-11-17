using System;
using Burkardt.Types;

namespace Burkardt.RandomNS;

public static class Stri
{
    public static double stri_angles_to_area(double r, double a, double b, double c)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    STRI_ANGLES_TO_AREA computes the area of a spherical triangle.
        //
        //  Discussion:
        //
        //    A sphere centered at 0 in 3D satisfies the equation:
        //
        //      X**2 + Y**2 + Z**2 = R**2
        //
        //    A spherical triangle is specified by three points on the surface
        //    of the sphere.
        //
        //    The area formula is known as Girard's formula.
        //
        //    The area of a spherical triangle is:
        //
        //      AREA = ( A + B + C - PI ) * R**2
        //
        //    where A, B and C are the (surface) angles of the triangle.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the sphere.
        //
        //    Input, double A, B, C, the angles of the triangle.
        //
        //    Output, double STRI_ANGLES_TO_AREA_3D, the area of the spherical triangle.
        //
    {
        double area = r * r * (a + b + c - Math.PI);

        return area;
    }

    public static void stri_sides_to_angles(double r, double as_, double bs, double cs,
            ref double a, ref double b, ref double c)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    STRI_SIDES_TO_ANGLES computes spherical triangle angles.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the sphere.
        //
        //    Input, double AS, BS, CS, the (geodesic) length of the sides of the
        //    triangle.
        //
        //    Output, double *A, *B, *C, the spherical angles of the triangle.
        //    Angle A is opposite the side of length AS, and so on.
        //
    {
        double asu;
        double bsu;
        double csu;
        double ssu;
        double tan_a2;
        double tan_b2;
        double tan_c2;

        asu =  as_ / r;
        bsu = bs / r;
        csu = cs / r;
        ssu = (asu + bsu + csu) / 2.0;

        tan_a2 = Math.Sqrt(Math.Sin(ssu - bsu) * Math.Sin(ssu - csu) /
                           (Math.Sin(ssu) * Math.Sin(ssu - asu)));

        a = 2.0 * Math.Atan(tan_a2);

        tan_b2 = Math.Sqrt(Math.Sin(ssu - asu) * Math.Sin(ssu - csu) /
                           (Math.Sin(ssu) * Math.Sin(ssu - bsu)));

        b = 2.0 * Math.Atan(tan_b2);

        tan_c2 = Math.Sqrt(Math.Sin(ssu - asu) * Math.Sin(ssu - bsu) /
                           (Math.Sin(ssu) * Math.Sin(ssu - csu)));

        c = 2.0 * Math.Atan(tan_c2);
    }

    public static void stri_vertices_to_sides(double r, double[] v1, double[] v2,
            double[] v3, ref double as_, ref double bs, ref double cs )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    STRI_VERTICES_TO_SIDES_3D computes spherical triangle sides.
        //
        //  Discussion:
        //
        //    We can use the ACOS system call here, but the R8_ACOS routine
        //    will automatically take care of cases where the input argument is
        //    (usually slightly) out of bounds.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 January 2016
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, the radius of the sphere.
        //
        //    Input, double V1[3], V2[3], V3[3], the vertices of the spherical
        //    triangle.
        //
        //    Output, double *AS, *BS, *CS, the (geodesic) length of the sides of the
        //    triangle.
        //
    {
        as_ = r* Math.Acos(typeMethods.r8vec_dot_product
            (3, v2, v3) / (r * r) );
        bs = r * Math.Acos(typeMethods.r8vec_dot_product(3, v3, v1) / (r * r));
        cs = r * Math.Acos(typeMethods.r8vec_dot_product(3, v1, v2) / (r * r));
    }
}