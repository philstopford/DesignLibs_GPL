using System;
using Burkardt.Types;

namespace Burkardt
{
    public static class Rotation
    {
        public static double[] rotation_axis_vector(double[] axis, double angle, double[] v)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ROTATION_AXIS_VECTOR rotates a vector around an axis vector in 3D.
            //
            //  Discussion:
            //
            //    Thanks to Cody Farnell for correcting some mistakes in an earlier
            //    version of this routine.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    30 May 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double AXIS[3], the axis vector for the rotation.
            //
            //    Input, double ANGLE, the angle, in radians, of the rotation.
            //
            //    Input, double V[3], the vector to be rotated.
            //
            //    Output, double ROTATION_AXIS_VECTOR[3], the rotated vector.
            //
        {
            double axis_norm;
            double dot;
            double norm;
            double[] normal = new double[3];
            double normal_component;
            double[] normal2;
            double[] parallel = new double[3];
            double[] rot = new double[3];
            double[] u;
            double[] w;
            //
            //  Compute the length of the rotation axis.
            //
            u = typeMethods.r8vec_copy_new(3, axis);

            axis_norm = typeMethods.r8vec_norm(3, u);

            if (axis_norm == 0.0)
            {
                w = typeMethods.r8vec_zeros_new(3);
                return w;
            }

            u[0] = u[0] / axis_norm;
            u[1] = u[1] / axis_norm;
            u[2] = u[2] / axis_norm;
            //
            //  Compute the dot product of the vector and the unit rotation axis.
            //
            dot = typeMethods.r8vec_dot_product(3, u, v);
            //
            //  Compute the parallel component of the vector.
            //
            parallel[0] = dot * u[0];
            parallel[1] = dot * u[1];
            parallel[2] = dot * u[2];
            //
            //  Compute the normal component of the vector.
            //
            normal[0] = v[0] - parallel[0];
            normal[1] = v[1] - parallel[1];
            normal[2] = v[2] - parallel[2];

            normal_component = typeMethods.r8vec_norm(3, normal);

            if (normal_component == 0.0)
            {
                w = typeMethods.r8vec_copy_new(3, parallel);
                return w;
            }

            normal[0] = normal[0] / normal_component;
            normal[1] = normal[1] / normal_component;
            normal[2] = normal[2] / normal_component;
            //
            //  Compute a second vector, lying in the plane, perpendicular
            //  to V, and forming a right-handed system.
            //
            normal2 = typeMethods.r8vec_cross_product(u, normal);

            norm = typeMethods.r8vec_norm(3, normal2);

            normal2[0] = normal2[0] / norm;
            normal2[1] = normal2[1] / norm;
            normal2[2] = normal2[2] / norm;
            //
            //  Rotate the normal component by the angle.
            //
            rot[0] = normal_component * (Math.Cos(angle) * normal[0]
                                         + Math.Sin(angle) * normal2[0]);

            rot[1] = normal_component * (Math.Cos(angle) * normal[1]
                                         + Math.Sin(angle) * normal2[1]);

            rot[2] = normal_component * (Math.Cos(angle) * normal[2]
                                         + Math.Sin(angle) * normal2[2]);

            //
            //  The rotated vector is the parallel component plus the rotated component.
            //
            w = new double[3];

            w[0] = parallel[0] + rot[0];
            w[1] = parallel[1] + rot[1];
            w[2] = parallel[2] + rot[2];


            return w;
        }

        public static double[] rotation_axis2mat(double[] axis, double angle)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ROTATION_AXIS2MAT converts a rotation from axis to matrix format in 3D.
            //
            //  Discussion:
            //
            //    The two dimensional array A is stored as a one dimensional vector, by columns.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    30 May 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    James Foley, Andries vanDam, Steven Feiner, John Hughes,
            //    Computer Graphics, Principles and Practice,
            //    Second Edition,
            //    Addison Wesley, 1990.
            //
            //  Parameters:
            //
            //    Input, double AXIS[3], the axis vector which remains unchanged by
            //    the rotation.
            //
            //    Input, double ANGLE, the angular measurement of the rotation about
            //    the axis, in radians.
            //
            //    Output, double ROTATION_AXIS2MAT[3*3], the rotation matrix.
            //
        {
            double[] a;
            double ca;
            double norm;
            double sa;
            double v1;
            double v2;
            double v3;

            v1 = axis[0];
            v2 = axis[1];
            v3 = axis[2];

            norm = Math.Sqrt(v1 * v1 + v2 * v2 + v3 * v3);

            if (norm == 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("ROTATION_AXIS2MAT - Fatal error!");
                Console.WriteLine("  Axis vector has zero norm.");
                return null;
            }

            v1 = v1 / norm;
            v2 = v2 / norm;
            v3 = v3 / norm;

            ca = Math.Cos(angle);
            sa = Math.Sin(angle);

            a = new double[3 * 3];

            a[0 + 0 * 3] = v1 * v1 + ca * (1.0 - v1 * v1);
            a[1 + 0 * 3] = (1.0 - ca) * v2 * v1 + sa * v3;
            a[2 + 0 * 3] = (1.0 - ca) * v3 * v1 - sa * v2;

            a[0 + 1 * 3] = (1.0 - ca) * v1 * v2 - sa * v3;
            a[1 + 1 * 3] = v2 * v2 + ca * (1.0 - v2 * v2);
            a[2 + 1 * 3] = (1.0 - ca) * v3 * v2 + sa * v1;

            a[0 + 2 * 3] = (1.0 - ca) * v1 * v3 + sa * v2;
            a[1 + 2 * 3] = (1.0 - ca) * v2 * v3 - sa * v1;
            a[2 + 2 * 3] = v3 * v3 + ca * (1.0 - v3 * v3);

            return a;
        }

        public static double[] rotation_axis2quat(double[] axis, double angle)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ROTATION_AXIS2QUAT converts a rotation from axis to quaternion format in 3D.
            //
            //  Discussion:
            //
            //    A rotation quaternion Q has the form:
            //
            //      Q = A + Bi + Cj + Dk
            //
            //    where A, B, C and D are double numbers, and i, j, and k are to be regarded
            //    as symbolic constant basis vectors, similar to the role of the "i"
            //    in the representation of imaginary numbers.
            //
            //    A is the cosine of half of the angle of rotation.  (B,C,D) is a
            //    unit vector pointing in the direction of the axis of rotation.
            //    Rotation multiplication and inversion can be carried out using
            //    this format and the usual rules for quaternion multiplication
            //    and inversion.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    30 May 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double AXIS[3], the axis vector which remains unchanged by
            //    the rotation.
            //
            //    Input, double ANGLE, the angular measurement of the rotation about
            //    the axis, in radians.
            //
            //    Output, double ROTATION_AXIS2QUAT[4], the quaternion representing the rotation.
            //
        {
            double norm;
            double[] q;

            norm = typeMethods.r8vec_norm(3, axis);

            if (norm == 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("ROTATION_AXIS2QUAT - Fatal error!");
                Console.WriteLine("  The axis vector is null.");
                return null;
            }

            q = new double [4];

            q[0] = Math.Cos(0.5 * angle);

            q[1] = axis[0] * Math.Sin(0.5 * angle) / norm;
            q[2] = axis[1] * Math.Sin(0.5 * angle) / norm;
            q[3] = axis[2] * Math.Sin(0.5 * angle) / norm;

            return q;
        }

        public static double[] rotation_mat_vector(double[] a, double[] v)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ROTATION_MAT_VECTOR applies a marix rotation to a vector in 3d.
            //
            //  Discussion:
            //
            //    The two dimensional array A is stored as a one dimensional vector, by columns.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    30 May 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double A[3*3], the matrix defining the rotation.
            //
            //    Input, double V[3], the vector to be rotated.
            //
            //    Output, double ROTATION_MAT_VECTOR[3], the rotated vector.
            //
        {
            double[] w;

            w = typeMethods.r8mat_mv_new(3, 3, a, v);

            return w;
        }

        public static void rotation_mat2axis(double[] a, double[] axis, ref double angle)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ROTATION_MAT2AXIS converts a rotation from matrix to axis format in 3D.
            //
            //  Discussion:
            //
            //    The two dimensional array A is stored as a one dimensional vector, by columns.
            //
            //    The computation is based on the fact that a rotation matrix must
            //    have an eigenvector corresponding to the eigenvalue of 1, hence:
            //
            //      ( A - I ) * v = 0.
            //
            //    The eigenvector V is the axis of rotation.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    30 May 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Jack Kuipers,
            //    Quaternions and Rotation Sequences,
            //    Princeton, 1998.
            //
            //  Parameters:
            //
            //    Input, double A[3*3], the rotation matrix.
            //
            //    Output, double AXIS[3], the axis vector which remains unchanged by
            //    the rotation.
            //
            //    Output, double &ANGLE, the angular measurement of the rotation about
            //    the axis, in radians.
            //
        {
            double norm;

            norm = Math.Sqrt((a[2 + 1 * 3] - a[1 + 2 * 3]) * (a[2 + 1 * 3] - a[1 + 2 * 3])
                             + (a[0 + 2 * 3] - a[2 + 0 * 3]) * (a[0 + 2 * 3] - a[2 + 0 * 3])
                             + (a[1 + 0 * 3] - a[0 + 1 * 3]) * (a[1 + 0 * 3] - a[0 + 1 * 3]));

            if (norm == 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("ROTATION_MAT2AXIS - Fatal error!");
                Console.WriteLine("  A is not a rotation matrix,");
                Console.WriteLine("  or there are multiple axes of rotation.");
                return;
            }

            axis[0] = (a[2 + 1 * 3] - a[1 + 2 * 3]) / norm;
            axis[1] = (a[0 + 2 * 3] - a[2 + 0 * 3]) / norm;
            axis[2] = (a[1 + 0 * 3] - a[0 + 1 * 3]) / norm;
            //
            //  Find the angle.
            //
            angle = typeMethods.r8_acos(0.5 *
                                        (a[0 + 0 * 3] + a[1 + 1 * 3] + a[2 + 2 * 3] - 1.0));

        }

        public static double[] rotation_mat2quat(double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ROTATION_MAT2QUAT converts a rotation from matrix to quaternion format in 3D.
            //
            //  Discussion:
            //
            //    The two dimensional array A is stored as a one dimensional vector, by columns.
            //
            //    The computation is based on the fact that a rotation matrix must
            //    have an eigenvector corresponding to the eigenvalue of 1, hence:
            //
            //      ( A - I ) * v = 0.
            //
            //    The eigenvector V is the axis of rotation.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    30 May 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Jack Kuipers,
            //    Quaternions and Rotation Sequences,
            //    Princeton, 1998.
            //
            //  Parameters:
            //
            //    Input, double A[3*3], the rotation matrix.
            //
            //    Output, double ROTATION_MAT2QUAT[4], the quaternion representing the rotation.
            //
        {
            double angle;
            double cos_phi;
            double norm;
            double[] q;
            double sin_phi;

            norm = Math.Sqrt((a[2 + 1 * 3] - a[1 + 2 * 3]) * (a[2 + 1 * 3] - a[1 + 2 * 3])
                             + (a[0 + 2 * 3] - a[2 + 0 * 3]) * (a[0 + 2 * 3] - a[2 + 0 * 3])
                             + (a[1 + 0 * 3] - a[0 + 1 * 3]) * (a[1 + 0 * 3] - a[0 + 1 * 3]));

            if (norm == 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("ROTATION_MAT2QUAT - Fatal error!");
                Console.WriteLine("  A is not a rotation matrix,");
                Console.WriteLine("  or there are multiple axes of rotation.");
                return null;
            }

            angle = typeMethods.r8_acos(0.5 *
                                        (a[0 + 0 * 3] + a[1 + 1 * 3] + a[2 + 2 * 3] - 1.0));

            cos_phi = Math.Cos(0.5 * angle);

            sin_phi = Math.Sqrt(1.0 - cos_phi * cos_phi);

            q = new double[4];

            q[0] = cos_phi;
            q[1] = sin_phi * (a[2 + 1 * 3] - a[1 + 2 * 3]) / norm;
            q[2] = sin_phi * (a[0 + 2 * 3] - a[2 + 0 * 3]) / norm;
            q[3] = sin_phi * (a[1 + 0 * 3] - a[0 + 1 * 3]) / norm;

            return q;
        }

        public static double[] rotation_quat_vector(double[] q, double[] v)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ROTATION_QUAT_VECTOR applies a quaternion rotation to a vector in 3d.
            //
            //  Discussion:
            //
            //    If Q is a unit quaternion that encodes a rotation of ANGLE
            //    radians about the vector AXIS, then for an arbitrary real
            //    vector V, the result W of the rotation on V can be written as:
            //
            //      W = Q * V * Conj(Q)
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    30 May 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double Q[4], the quaternion defining the rotation.
            //
            //    Input, double V[3], the vector to be rotated.
            //
            //    Output, double ROTATION_QUAT_VECTOR[3], the rotated vector.
            //
        {
            double[] w;

            w = new double[3];

            w[0] =
                (2.0 * (q[0] * q[0] + q[1] * q[1]) - 1.0) * v[0]
                + 2.0 * (q[1] * q[2] - q[0] * q[3]) * v[1]
                + 2.0 * (q[1] * q[3] + q[0] * q[2]) * v[2];

            w[1] =
                2.0 * (q[1] * q[2] + q[0] * q[3]) * v[0]
                + (2.0 * (q[0] * q[0] + q[2] * q[2]) - 1.0) * v[1]
                + 2.0 * (q[2] * q[3] - q[0] * q[1]) * v[2];

            w[2] =
                2.0 * (q[1] * q[3] - q[0] * q[2]) * v[0]
                + 2.0 * (q[2] * q[3] + q[0] * q[1]) * v[1]
                + (2.0 * (q[0] * q[0] + q[3] * q[3]) - 1.0) * v[2];

            return w;
        }

        public static void rotation_quat2axis(double[] q, double[] axis, ref double angle)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ROTATION_QUAT2AXIS converts a rotation from quaternion to axis format in 3D.
            //
            //  Discussion:
            //
            //    A rotation quaternion Q has the form:
            //
            //      Q = A + Bi + Cj + Dk
            //
            //    where A, B, C and D are double numbers, and i, j, and k are to be regarded
            //    as symbolic constant basis vectors, similar to the role of the "i"
            //    in the representation of imaginary numbers.
            //
            //    A is the cosine of half of the angle of rotation.  (B,C,D) is a
            //    vector pointing in the direction of the axis of rotation.
            //    Rotation multiplication and inversion can be carried out using
            //    this format and the usual rules for quaternion multiplication
            //    and inversion.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    05 August 2018
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double Q[4], the quaternion representing the rotation.
            //
            //    Output, double AXIS[3], the axis vector which remains unchanged by
            //    the rotation.
            //
            //    Output, double &ANGLE, the angular measurement of the rotation about
            //    the axis, in radians.
            //
        {
            double cos_phi;
            double sin_phi;

            sin_phi = Math.Sqrt(q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);

            cos_phi = q[0];

            angle = 2.0 * Math.Atan2(sin_phi, cos_phi);

            if (sin_phi == 0.0)
            {
                axis[0] = 1.0;
                axis[1] = 0.0;
                axis[2] = 0.0;
            }
            else
            {
                axis[0] = q[1] / sin_phi;
                axis[1] = q[2] / sin_phi;
                axis[2] = q[3] / sin_phi;
            }
        }

        public static double[] rotation_quat2mat(double[] q)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ROTATION_QUAT2MAT converts a rotation from quaternion to matrix format in 3D.
            //
            //  Discussion:
            //
            //    The two dimensional array A is stored as a one dimensional vector, by columns.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    30 May 2010
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    James Foley, Andries vanDam, Steven Feiner, John Hughes,
            //    Computer Graphics, Principles and Practice,
            //    Second Edition,
            //    Addison Wesley, 1990.
            //
            //  Parameters:
            //
            //    Input, double Q[4], the quaternion representing the rotation.
            //
            //    Output, double ROTATION_QUAT2MAT[3*3], the rotation matrix.
            //
        {
            double[] a;
            double angle;
            double ca;
            double cos_phi;
            double sa;
            double sin_phi;
            double v1;
            double v2;
            double v3;

            sin_phi = Math.Sqrt(q[1] * q[1] + q[2] * q[2] + q[3] * q[3]);

            cos_phi = q[0];

            angle = 2.0 * Math.Atan2(sin_phi, cos_phi);

            if (sin_phi == 0.0)
            {
                v1 = 1.0;
                v2 = 0.0;
                v3 = 0.0;
            }
            else
            {
                v1 = q[1] / sin_phi;
                v2 = q[2] / sin_phi;
                v3 = q[3] / sin_phi;
            }

            ca = Math.Cos(angle);
            sa = Math.Sin(angle);

            a = new double[3 * 3];

            a[0 + 0 * 3] = v1 * v1 + ca * (1.0 - v1 * v1);
            a[1 + 0 * 3] = (1.0 - ca) * v2 * v1 + sa * v3;
            a[2 + 0 * 3] = (1.0 - ca) * v3 * v1 - sa * v2;

            a[0 + 1 * 3] = (1.0 - ca) * v1 * v2 - sa * v3;
            a[1 + 1 * 3] = v2 * v2 + ca * (1.0 - v2 * v2);
            a[2 + 1 * 3] = (1.0 - ca) * v3 * v2 + sa * v1;

            a[0 + 2 * 3] = (1.0 - ca) * v1 * v3 + sa * v2;
            a[1 + 2 * 3] = (1.0 - ca) * v2 * v3 - sa * v1;
            a[2 + 2 * 3] = v3 * v3 + ca * (1.0 - v3 * v3);

            return a;
        }
    }
}