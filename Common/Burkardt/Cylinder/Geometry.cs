using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.Cylinder;

public static class Geometry
{
    public static double cylinder_point_dist_3d(double[] p1, double[] p2, double r,
            double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CYLINDER_POINT_DIST_3D determines the distance from a cylinder to a point in 3D.
        //
        //  Discussion:
        //
        //    We are computing the distance to the SURFACE of the cylinder.
        //
        //    The surface of a (right) (finite) cylinder in 3D is defined by an axis,
        //    which is the line segment from point P1 to P2, and a radius R.  The points
        //    on the surface of the cylinder are:
        //    * points at a distance R from the line through P1 and P2, and whose nearest
        //      point on the line through P1 and P2 is strictly between P1 and P2,
        //    PLUS
        //    * points at a distance less than or equal to R from the line through P1
        //      and P2, whose nearest point on the line through P1 and P2 is either
        //      P1 or P2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1[3], P2[3], the first and last points
        //    on the axis line of the cylinder.
        //
        //    Input, double R, the radius of the cylinder.
        //
        //    Input, double P[3], the point.
        //
        //    Output, double CYLINDER_POINT_DIST_3D, the distance from the point
        //    to the cylinder.
        //
    {
        int DIM_NUM = 3;

        double[] axis = new double[DIM_NUM];
        double axis_length = 0.0;
        double distance = 0.0;
        int i;
        double off_axis_component = 0.0;
        double p_dot_axis = 0.0;
        double p_length = 0.0;
        double[] v1 = new double[DIM_NUM];

        for (i = 0; i < DIM_NUM; i++)
        {
            axis[i] = p2[i] - p1[i];
        }

        axis_length = typeMethods.r8vec_norm(DIM_NUM, axis);

        switch (axis_length)
        {
            case 0.0:
                distance = -typeMethods.r8_huge();
                return distance;
        }

        for (i = 0; i < DIM_NUM; i++)
        {
            axis[i] /= axis_length;
        }

        p_dot_axis = 0.0;
        for (i = 0; i < DIM_NUM; i++)
        {
            p_dot_axis += (p[i] - p1[i]) * axis[i];
        }

        switch (p_dot_axis)
        {
            //
            //  Case 1: Below bottom cap.
            //
            case <= 0.0:
                distance = Disk.Geometry.disk_point_dist_3d(p1, r, axis, p);
                break;
            //
            default:
            {
                if (p_dot_axis <= axis_length)
                {
                    for (i = 0; i < DIM_NUM; i++)
                    {
                        v1[i] = p[i] - p1[i];
                    }

                    p_length = typeMethods.r8vec_norm(DIM_NUM, v1);
                    off_axis_component = Math.Sqrt(Math.Pow(p_length, 2) - Math.Pow(p_dot_axis, 2));

                    distance = Math.Abs(off_axis_component - r);

                    if (off_axis_component < r)
                    {
                        distance = Math.Min(distance, axis_length - p_dot_axis);
                        distance = Math.Min(distance, p_dot_axis);
                    }
                }
                //
                //  Case 3: Above the top cap.
                //
                else if (axis_length < p_dot_axis)
                {
                    distance = Disk.Geometry.disk_point_dist_3d(p2, r, axis, p);
                }

                break;
            }
        }

        return distance;
    }

    public static double cylinder_point_dist_signed_3d(double[] p1, double[] p2, double r,
            double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CYLINDER_POINT_DIST_SIGNED_3D: signed distance from cylinder to point in 3D.
        //
        //  Discussion:
        //
        //    We are computing the signed distance to the SURFACE of the cylinder.
        //
        //    The surface of a (right) (finite) cylinder in 3D is defined by an axis,
        //    which is the line segment from point P1 to P2, and a radius R.  The points
        //    on the surface of the cylinder are:
        //    * points at a distance R from the line through P1 and P2, and whose nearest
        //      point on the line through P1 and P2 is strictly between P1 and P2,
        //    PLUS
        //    * points at a distance less than or equal to R from the line through P1
        //      and P2, whose nearest point on the line through P1 and P2 is either
        //      P1 or P2.
        //
        //    Points inside the surface have a negative distance.
        //    Points on the surface have a zero distance.
        //    Points outside the surface have a positive distance.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    26 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1[3], P2[3], the first and last points
        //    on the axis line of the cylinder.
        //
        //    Input, double R, the radius of the cylinder.
        //
        //    Input, double P[3], the point.
        //
        //    Output, double CYLINDER_POINT_DIST_SIGNED_3D, the signed distance from the point
        //    to the cylinder.
        //
    {
        int DIM_NUM = 3;

        double[] axis = new double[DIM_NUM];
        double axis_length = 0.0;
        double distance = 0.0;
        int i;
        double off_axis_component = 0.0;
        double p_dot_axis = 0.0;
        double p_length = 0.0;
        double[] v1 = new double[DIM_NUM];

        for (i = 0; i < DIM_NUM; i++)
        {
            axis[i] = p2[i] - p1[i];
        }

        axis_length = typeMethods.r8vec_norm(DIM_NUM, axis);

        switch (axis_length)
        {
            case 0.0:
                distance = -typeMethods.r8_huge();
                return distance;
        }

        for (i = 0; i < DIM_NUM; i++)
        {
            axis[i] /= axis_length;
        }

        p_dot_axis = 0.0;
        for (i = 0; i < DIM_NUM; i++)
        {
            p_dot_axis += (p[i] - p1[i]) * axis[i];
        }

        switch (p_dot_axis)
        {
            //
            //  Case 1: Below bottom cap.
            //
            case <= 0.0:
                distance = Disk.Geometry.disk_point_dist_3d(p1, r, axis, p);
                break;
            //
            default:
            {
                if (p_dot_axis <= axis_length)
                {
                    for (i = 0; i < DIM_NUM; i++)
                    {
                        v1[i] = p[i] - p1[i];
                    }

                    p_length = typeMethods.r8vec_norm(DIM_NUM, v1);
                    off_axis_component = Math.Sqrt(Math.Pow(p_length, 2) - Math.Pow(p_dot_axis, 2));

                    distance = off_axis_component - r;

                    switch (distance)
                    {
                        case < 0.0:
                            distance = Math.Max(distance, p_dot_axis - axis_length);
                            distance = Math.Max(distance, -p_dot_axis);
                            break;
                    }
                }
                //
                //  Case 3: Above the top cap.
                //
                else if (axis_length < p_dot_axis)
                {
                    distance = Disk.Geometry.disk_point_dist_3d(p2, r, axis, p);
                }

                break;
            }
        }

        return distance;
    }

    public static bool cylinder_point_inside_3d(double[] p1, double[] p2, double r,
            double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CYLINDER_POINT_INSIDE_3D determines if a cylinder contains a point in 3D.
        //
        //  Discussion:
        //
        //    The surface and interior of a (right) (finite) cylinder in 3D is defined
        //    by an axis, which is the line segment from point P1 to P2, and a
        //    radius R.  The points contained in the volume include:
        //    * points at a distance less than or equal to R from the line through P1
        //      and P2, whose nearest point on the line through P1 and P2 is, in fact,
        //      P1, P2, or any point between them.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1[3], P2[3], the first and last points
        //    on the axis line of the cylinder.
        //
        //    Input, double R, the radius of the cylinder.
        //
        //    Input, double P[3], the point.
        //
        //    Output, bool CYLINDER_POINT_INSIDE_3D, is TRUE if the point is inside
        //    the cylinder.
        //
    {
        int DIM_NUM = 3;

        double[] axis = new double[DIM_NUM];
        double axis_length;
        int i;
        bool inside;
        double off_axis_component;
        double p_dot_axis;
        double p_length;
        double[] v1 = new double[DIM_NUM];

        for (i = 0; i < DIM_NUM; i++)
        {
            axis[i] = p2[i] - p1[i];
        }

        axis_length = typeMethods.r8vec_norm(DIM_NUM, axis);

        switch (axis_length)
        {
            case 0.0:
                inside = false;
                return inside;
        }

        for (i = 0; i < DIM_NUM; i++)
        {
            axis[i] /= axis_length;
        }

        p_dot_axis = 0.0;
        for (i = 0; i < DIM_NUM; i++)
        {
            p_dot_axis += (p[i] - p1[i]) * axis[i];
        }

        //
        //  If the point lies below or above the "caps" of the cylinder, we're done.
        //
        if (p_dot_axis < 0.0 || axis_length < p_dot_axis)
        {
            inside = false;
        }
        //
        //  Otherwise, determine the distance from P to the axis.
        //
        else
        {
            for (i = 0; i < DIM_NUM; i++)
            {
                v1[i] = p[i] - p1[i];
            }

            p_length = typeMethods.r8vec_norm(DIM_NUM, v1);

            off_axis_component = Math.Sqrt(Math.Pow(p_length, 2) - Math.Pow(p_dot_axis, 2));

            if (off_axis_component <= r)
            {
                inside = true;
            }
            else
            {
                inside = false;
            }
        }

        return inside;
    }

    public static double[] cylinder_point_near_3d(double[] p1, double[] p2, double r,
            double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CYLINDER_POINT_NEAR_3D: nearest point on a cylinder to a point in 3D.
        //
        //  Discussion:
        //
        //    We are computing the nearest point on the SURFACE of the cylinder.
        //
        //    The surface of a (right) (finite) cylinder in 3D is defined by an axis,
        //    which is the line segment from point P1 to P2, and a radius R.  The points
        //    on the surface of the cylinder are:
        //    * points at a distance R from the line through P1 and P2, and whose nearest
        //      point on the line through P1 and P2 is strictly between P1 and P2,
        //    PLUS
        //    * points at a distance less than or equal to R from the line through P1
        //      and P2, whose nearest point on the line through P1 and P2 is either
        //      P1 or P2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    19 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1[3], P2[3], the first and last points
        //    on the axis line of the cylinder.
        //
        //    Input, double R, the radius of the cylinder.
        //
        //    Input, double P[3], the point.
        //
        //    Output, double CYLINDER_POINT_NEAR_3D[3], the nearest point on the cylinder.
        //
    {
        int DIM_NUM = 3;

        double axial_component;
        double[] axis = new double[DIM_NUM];
        double axis_length;
        double distance;
        int i;
        double[] normal;
        double[] off_axis = new double[DIM_NUM];
        double off_axis_component;
        double[] pn;

        pn = new double[DIM_NUM];

        for (i = 0; i < DIM_NUM; i++)
        {
            axis[i] = p2[i] - p1[i];
        }

        axis_length = typeMethods.r8vec_norm(DIM_NUM, axis);
        for (i = 0; i < DIM_NUM; i++)
        {
            axis[i] /= axis_length;
        }

        axial_component = 0.0;
        for (i = 0; i < DIM_NUM; i++)
        {
            axial_component += (p[i] - p1[i]) * axis[i];
        }

        for (i = 0; i < DIM_NUM; i++)
        {
            off_axis[i] = p[i] - p1[i] - axial_component * axis[i];
        }

        off_axis_component = typeMethods.r8vec_norm(DIM_NUM, off_axis);
        switch (axial_component)
        {
            //
            //  Case 1: Below bottom cap.
            //
            case <= 0.0 when off_axis_component <= r:
            {
                for (i = 0; i < DIM_NUM; i++)
                {
                    pn[i] = p1[i] + off_axis[i];
                }

                break;
            }
            case <= 0.0:
            {
                for (i = 0; i < DIM_NUM; i++)
                {
                    pn[i] = p1[i] + r / off_axis_component * off_axis[i];
                }

                break;
            }
            //
            default:
            {
                if (axial_component <= axis_length)
                {
                    switch (off_axis_component)
                    {
                        case 0.0:
                        {
                            normal = typeMethods.r8vec_any_normal(DIM_NUM, axis);
                            for (i = 0; i < DIM_NUM; i++)
                            {
                                pn[i] = p[i] + r * normal[i];
                            }

                            break;
                        }
                        default:
                        {
                            distance = Math.Abs(off_axis_component - r);

                            for (i = 0; i < DIM_NUM; i++)
                            {
                                pn[i] = p1[i] + axial_component * axis[i]
                                              + r / off_axis_component * off_axis[i];
                            }

                            if (off_axis_component < r)
                            {
                                if (axis_length - axial_component < distance)
                                {
                                    distance = axis_length - axial_component;
                                    for (i = 0; i < DIM_NUM; i++)
                                    {
                                        pn[i] = p2[i] + off_axis[i];
                                    }
                                }

                                if (axial_component < distance)
                                {
                                    distance = axial_component;
                                    for (i = 0; i < DIM_NUM; i++)
                                    {
                                        pn[i] = p1[i] + off_axis[i];
                                    }
                                }
                            }

                            break;
                        }
                    }
                }
                //
                //  Case 3: Above the top cap.
                //
                else if (axis_length < axial_component)
                {
                    if (off_axis_component <= r)
                    {
                        for (i = 0; i < DIM_NUM; i++)
                        {
                            pn[i] = p2[i] + off_axis[i];
                        }
                    }
                    else
                    {
                        for (i = 0; i < DIM_NUM; i++)
                        {
                            pn[i] = p2[i] + r / off_axis_component * off_axis[i];
                        }
                    }
                }

                break;
            }
        }

        return pn;
    }

    public static double[] cylinder_sample_3d(double[] p1, double[] p2, double r, int n,
            ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CYLINDER_SAMPLE_3D samples a cylinder in 3D.
        //
        //  Discussion:
        //
        //    We are sampling the interior of a right finite cylinder in 3D.
        //
        //    The interior of a (right) (finite) cylinder in 3D is defined by an axis,
        //    which is the line segment from point P1 to P2, and a radius R.  The points
        //    on or inside the cylinder are:
        //    * points whose distance from the line through P1 and P2 is less than
        //      or equal to R, and whose nearest point on the line through P1 and P2
        //      lies (nonstrictly) between P1 and P2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1[3], P2[3], the first and last points
        //    on the axis line of the cylinder.
        //
        //    Input, double R, the radius of the cylinder.
        //
        //    Input, int N, the number of sample points to compute.
        //
        //    Input/output, int &SEED, the random number seed.
        //
        //    Input, double CYLINDER_SAMPLE_3D[3*N], the sample points.
        //
    {
        int DIM_NUM = 3;

        double[] axis = new double[DIM_NUM];
        double axis_length;
        int i;
        int j;
        double[] p;
        double radius;
        double theta;
        double[] v2 = new double[DIM_NUM];
        double[] v3 = new double[DIM_NUM];
        double z;
        //
        //  Compute the axis vector.
        //
        for (i = 0; i < DIM_NUM; i++)
        {
            axis[i] = p2[i] - p1[i];
        }

        axis_length = typeMethods.r8vec_norm(DIM_NUM, axis);
        for (i = 0; i < DIM_NUM; i++)
        {
            axis[i] /= axis_length;
        }

        //
        //  Compute vectors V2 and V3 that form an orthogonal triple with AXIS.
        //
        Plane.Geometry.plane_normal_basis_3d(p1, axis, ref v2, ref v3);
        //
        //  Assemble the randomized information.
        //
        p = new double[DIM_NUM * n];

        for (j = 0; j < n; j++)
        {
            for (i = 0; i < DIM_NUM; i++)
            {
                radius = r * Math.Sqrt(UniformRNG.r8_uniform_01(ref seed));
                theta = 2.0 * Math.PI * UniformRNG.r8_uniform_01(ref seed);
                z = axis_length * UniformRNG.r8_uniform_01(ref seed);

                p[i + j * DIM_NUM] = p1[i]
                                     + z * axis[i]
                                     + radius * Math.Cos(theta) * v2[i]
                                     + radius * Math.Sin(theta) * v3[i];
            }
        }

        return p;
    }

    public static double cylinder_volume_3d(double[] p1, double[] p2, double r)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CYLINDER_VOLUME_3D determines the volume of a cylinder in 3D.
        //
        //  Discussion:
        //
        //    A (right) (finite) cylinder in 3D is the set of points
        //    contained on or inside a circle of radius R, whose center
        //    lies along the line segment from point P1 to P2, and whose
        //    plane is perpendicular to that line segment.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1[3], P2[3], the first and last points
        //    on the axis line of the cylinder.
        //
        //    Input, double R, the radius of the cylinder.
        //
        //    Output, double CYLINDER_VOLUME_3D, the volume of the cylinder.
        //
    {
        double h;
        double volume;

        h = typeMethods.r8vec_distance(3, p1, p2);

        volume = Math.PI * r * r * h;

        return volume;
    }

}