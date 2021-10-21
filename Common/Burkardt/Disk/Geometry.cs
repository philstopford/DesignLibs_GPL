using System;
using Burkardt.Types;

namespace Burkardt.Disk
{
    public static class Geometry
    {
        public static double disk_point_dist_3d(double[] pc, double r, double[] axis,
                double[] p)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DISK_POINT_DIST_3D determines the distance from a disk to a point in 3D.
            //
            //  Discussion:
            //
            //    A disk in 3D satisfies the equations:
            //
            //      ( P(1) - PC(1) )^2 + ( P(2) - PC(2) )^2 + ( P(3) - PC(3) <= R^2
            //
            //    and
            //
            //      P(1) * AXIS(1) + P(2) * AXIS(2) + P(3) * AXIS(3) = 0
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
            //    Input, double PC(3), the center of the disk.
            //
            //    Input, double R, the radius of the disk.
            //
            //    Input, double AXIS(3), the axis vector.
            //
            //    Input, double P(3), the point to be checked.
            //
            //    Output, double DISK_POINT_DIST_3D, the distance of the
            //    point to the disk.
            //
        {
            int DIM_NUM = 3;

            double axial_component;
            double axis_length;
            double dist;
            int i;
            double off_axis_component;
            double[] off_axis = new double[DIM_NUM];
            double[] v = new double[DIM_NUM];
            //
            //  Special case: the point is the center.
            //
            if (typeMethods.r8vec_eq(DIM_NUM, p, pc))
            {
                dist = 0.0;
                return dist;
            }

            axis_length = typeMethods.r8vec_norm(DIM_NUM, axis);

            if (axis_length <= 0.0)
            {
                dist = -typeMethods.r8_huge();
                return dist;
            }

            for (i = 0; i < DIM_NUM; i++)
            {
                v[i] = p[i] - pc[i];
            }

            axial_component = typeMethods.r8vec_dot_product(DIM_NUM, v, axis) / axis_length;
            //
            //  Special case: the point satisfies the disk equation exactly.
            //
            if (typeMethods.r8vec_norm(DIM_NUM, v) <= r && axial_component == 0.0)
            {
                dist = 0.0;
                return dist;
            }

            //
            //  Decompose P-PC into axis component and off-axis component.
            //
            for (i = 0; i < DIM_NUM; i++)
            {
                off_axis[i] = p[i] - pc[i] - axial_component * axis[i] / axis_length;
            }

            off_axis_component = typeMethods.r8vec_norm(DIM_NUM, off_axis);
            //
            //  If the off-axis component has norm less than R, the nearest point is
            //  the projection to the disk along the axial direction, and the distance
            //  is just the dot product of P-PC with unit AXIS.
            //
            if (off_axis_component <= r)
            {
                dist = Math.Abs(axial_component);
                return dist;
            }

            //
            //  Otherwise, the nearest point is along the perimeter of the disk.
            //
            dist = Math.Sqrt(Math.Pow(axial_component, 2)
                             + Math.Pow(off_axis_component - r, 2));

            return dist;
        }
    }
}