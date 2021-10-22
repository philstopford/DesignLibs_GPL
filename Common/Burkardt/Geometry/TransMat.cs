using System;

namespace Burkardt.Geometry
{
    public static class TransMat
    {
        public static void tmat_init(double[] a)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TMAT_INIT initializes the geometric transformation matrix.
            //
            //  Discussion:
            //
            //    The geometric transformation matrix can be thought of as a 4 by 4
            //    matrix "A" having components:
            //
            //      r11 r12 r13 t1
            //      r21 r22 r23 t2
            //      r31 r32 r33 t3
            //        0   0   0  1
            //
            //    This matrix encodes the rotations, scalings and translations that
            //    are applied to graphical objects.
            //
            //    A point P = (x,y,z) is rewritten in "homogeneous coordinates" as
            //    PH = (x,y,z,1).  Then to apply the transformations encoded in A to
            //    the point P, we simply compute A * PH.
            //
            //    Individual transformations, such as a scaling, can be represented
            //    by simple versions of the transformation matrix.  If the matrix
            //    A represents the current set of transformations, and we wish to
            //    apply a new transformation B, { the original points are
            //    transformed twice:  B * ( A * PH ).  The new transformation B can
            //    be combined with the original one A, to give a single matrix C that
            //    encodes both transformations: C = B * A.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    04 July 2005
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
            //    Input, double A[4*4], the geometric transformation matrix.
            //
        {
            int i;
            int j;

            for (i = 0; i < 4; i++)
            {
                for (j = 0; j < 4; j++)
                {
                    if (i == j)
                    {
                        a[i + j * 4] = 1.0;
                    }
                    else
                    {
                        a[i + j * 4] = 0.0;
                    }
                }
            }
        }

        public static void tmat_mxm(double[] a, double[] b, ref double[] c)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TMAT_MXM multiplies two geometric transformation matrices.
            //
            //  Discussion:
            //
            //    The product is accumulated in a temporary array, and { assigned
            //    to the result.  Therefore, it is legal for any two, or all three,
            //    of the arguments to share memory.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    19 October 1998
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
            //    Input, double A[4*4], the first geometric transformation matrix.
            //
            //    Input, double B[4*4], the second geometric transformation matrix.
            //
            //    Output, double C[4*4], the product A * B.
            //
        {
            double[] d = new double[4 * 4];
            int i;
            int j;
            int k;

            for (i = 0; i < 4; i++)
            {
                for (k = 0; k < 4; k++)
                {
                    d[i + k * 4] = 0.0;
                    for (j = 0; j < 4; j++)
                    {
                        d[i + k * 4] = d[i + k * 4] + a[i + j * 4] * b[j + k * 4];
                    }
                }
            }

            for (i = 0; i < 4; i++)
            {
                for (j = 0; j < 4; j++)
                {
                    c[i + j * 4] = d[i + j * 4];
                }
            }
        }

        public static void tmat_mxp(double[] a, double[] x, ref double[] y)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TMAT_MXP multiplies a geometric transformation matrix times a point.
            //
            //  Discussion:
            //
            //    The matrix will normally have the form
            //
            //      xx xy xz tx
            //      yx yy yz ty
            //      zx zy zz tz
            //       0  0  0  1
            //
            //    where the 3x3 initial block controls rotations and scalings,
            //    and the values [ tx, ty, tz ] implement a translation.
            //
            //    The matrix is stored as a vector, by COLUMNS.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 July 2005
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
            //    Input, double A[4*4], the geometric transformation matrix.
            //
            //    Input, double X[3], the point to be multiplied.  There is a
            //    "theoretical" fourth component of X, which can be assumed to
            //    equal 1.
            //
            //    Output, double Y[3], the result of A*X.  The product is accumulated in
            //    a temporary vector, and assigned to the result.  Therefore, it
            //    is legal for X and Y to share memory.
            //
        {
            int i;
            int j;
            double[] z = new double[3];

            for (i = 0; i < 3; i++)
            {
                z[i] = a[i + 3 * 4];
                for (j = 0; j < 3; j++)
                {
                    z[i] = z[i] + a[i + j * 4] * x[j];
                }
            }

            for (i = 0; i < 3; i++)
            {
                y[i] = z[i];
            }
        }

        public static void tmat_mxp2(double[] a, double[] p1, ref double[] p2, int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TMAT_MXP2 multiplies a geometric transformation matrix times N points.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    06 July 2005
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
            //    Input, double A[4*4], the geometric transformation matrix.
            //
            //    Input, double P1[3*N], the points to be multiplied.
            //
            //    Output, double P2[3*N], the transformed points.  Each product is
            //    accumulated in a temporary vector, and assigned to the
            //    result.  Therefore, it is legal for X and Y to share memory.
            //
        {
            int i;
            int j;
            int k;

            for (k = 0; k < n; k++)
            {
                for (i = 0; i < 3; i++)
                {
                    p2[i + k * 3] = a[i + 3 * 4];
                    for (j = 0; j < 3; j++)
                    {
                        p2[i + k * 3] = p2[i + k * 3] + a[i + j * 4] * p1[j + k * 3];
                    }
                }
            }
        }

        public static void tmat_mxv(double[] a, double[] x, ref double[] y)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TMAT_MXV multiplies a geometric transformation matrix times a vector.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 July 2005
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
            //    Input, double A[4*4], the geometric transformation matrix.
            //
            //    Input, double X[4], the vector to be multiplied.  The fourth component
            //    of X is implicitly assigned the value of 1.
            //
            //    Output, double Y[4], the result of A*X.  The product is accumulated in
            //    a temporary vector, and assigned to the result.  Therefore, it
            //    is legal for X and Y to share memory.
            //
        {
            int i;
            int j;
            double[] z = new double[4];

            for (i = 0; i < 3; i++)
            {
                z[i] = 0.0;
                for (j = 0; j < 3; j++)
                {
                    z[i] = z[i] + a[i + j * 4] * x[j];
                }

                z[i] = z[i] + a[i + 3 * 4];
            }

            for (i = 0; i < 3; i++)
            {
                y[i] = z[i];
            }
        }

        public static void tmat_rot_axis(double[] a, ref double[] b, double angle,
                char axis)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TMAT_ROT_AXIS applies an axis rotation to the geometric transformation matrix.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 July 2005
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
            //    Input, double A[4*4], the current geometric transformation matrix.
            //
            //    Output, double B[4*4], the modified geometric transformation matrix.
            //    A and B may share the same memory.
            //
            //    Input, double ANGLE, the angle, in degrees, of the rotation.
            //
            //    Input, character AXIS, is 'X', 'Y' or 'Z', specifying the coordinate
            //    axis about which the rotation occurs.
            //
        {
            double[] c = new double[4 * 4];
            double[] d = new double[4 * 4];
            int i;
            int j;
            double theta;

            theta = Helpers.degrees_to_radians(angle);

            tmat_init(c);

            if (axis == 'X' || axis == 'x')
            {
                c[1 + 1 * 4] = Math.Cos(theta);
                c[1 + 2 * 4] = -Math.Sin(theta);
                c[2 + 1 * 4] = Math.Sin(theta);
                c[2 + 2 * 4] = Math.Cos(theta);
            }
            else if (axis == 'Y' || axis == 'y')
            {
                c[0 + 0 * 4] = Math.Cos(theta);
                c[0 + 2 * 4] = Math.Sin(theta);
                c[2 + 0 * 4] = -Math.Sin(theta);
                c[2 + 2 * 4] = Math.Cos(theta);
            }
            else if (axis == 'Z' || axis == 'z')
            {
                c[0 + 0 * 4] = Math.Cos(theta);
                c[0 + 1 * 4] = -Math.Sin(theta);
                c[1 + 0 * 4] = Math.Sin(theta);
                c[1 + 1 * 4] = Math.Cos(theta);
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("TMAT_ROT_AXIS - Fatal error!");
                Console.WriteLine("  Illegal rotation axis: " + axis + "");
                Console.WriteLine("  Legal choices are 'X', 'Y', or 'Z'.");
                return;
            }

            tmat_mxm(c, a, ref d);

            for (i = 0; i < 4; i++)
            {
                for (j = 0; j < 4; j++)
                {
                    b[i + j * 4] = d[i + j * 4];
                }
            }
        }

        public static void tmat_rot_vector(double[] a, ref double[] b, double angle,
                double[] v)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TMAT_ROT_VECTOR applies a rotation about a vector to the geometric transformation matrix.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 July 2005
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
            //    Input, double A[4*4], the current geometric transformation matrix.
            //
            //    Output, double B[4*4], the modified geometric transformation matrix.
            //    A and B may share the same memory.
            //
            //    Input, double ANGLE, the angle, in degrees, of the rotation.
            //
            //    Input, double V[3], the coordinates of a (nonzero)
            //    point defining a vector from the origin.  The rotation will occur
            //    about this axis.
            //
        {
            double[] c = new double[4 * 4];
            double ca;
            double[] d = new double[4 * 4];
            int i;
            int j;
            double sa;
            double theta;

            if (Math.Pow(v[0], 2) + Math.Pow(v[1], 2) + Math.Pow(v[2], 2) == 0.0)
            {
                return;
            }

            theta = Helpers.degrees_to_radians(angle);

            tmat_init(c);

            ca = Math.Cos(theta);
            sa = Math.Sin(theta);

            c[0 + 0 * 4] = v[0] * v[0] + ca * (1.0 - v[0] * v[0]);
            c[0 + 1 * 4] = (1.0 - ca) * v[0] * v[1] - sa * v[2];
            c[0 + 2 * 4] = (1.0 - ca) * v[0] * v[2] + sa * v[1];

            c[1 + 0 * 4] = (1.0 - ca) * v[1] * v[0] + sa * v[2];
            c[1 + 1 * 4] = v[1] * v[1] + ca * (1.0 - v[1] * v[1]);
            c[1 + 2 * 4] = (1.0 - ca) * v[1] * v[2] - sa * v[0];

            c[2 + 0 * 4] = (1.0 - ca) * v[2] * v[0] - sa * v[1];
            c[2 + 1 * 4] = (1.0 - ca) * v[2] * v[1] + sa * v[0];
            c[2 + 2 * 4] = v[2] * v[2] + ca * (1.0 - v[2] * v[2]);

            tmat_mxm(c, a, ref d);

            for (i = 0; i < 4; i++)
            {
                for (j = 0; j < 4; j++)
                {
                    b[i + j * 4] = d[i + j * 4];
                }
            }
        }

        public static void tmat_scale(double[] a, ref double[] b, double[] s)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TMAT_SCALE applies a scaling to the geometric transformation matrix.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 July 2005
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
            //    Input, double A[4*4], the current geometric transformation matrix.
            //
            //    Output, double B[4*4], the modified geometric transformation matrix.
            //    A and B may share the same memory.
            //
            //    Input, double S[3], the scalings to be applied to the coordinates.
            //
        {
            double[] c = new double[4 * 4];
            double[] d = new double[4 * 4];
            int i;
            int j;

            tmat_init(c);

            c[0 + 0 * 4] = s[0];
            c[1 + 1 * 4] = s[1];
            c[2 + 2 * 4] = s[2];

            tmat_mxm(c, a, ref d);

            for (i = 0; i < 4; i++)
            {
                for (j = 0; j < 4; j++)
                {
                    b[i + j * 4] = d[i + j * 4];
                }
            }
        }

        public static void tmat_shear(double[] a, ref double[] b, string axis, double s)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TMAT_SHEAR applies a shear to the geometric transformation matrix.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 July 2005
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
            //    Input, double A[4*4], the current geometric transformation matrix.
            //
            //    Output, double B[4*4], the modified geometric transformation matrix.
            //    A and B may share the same memory.
            //
            //    Input, string AXIS, is 'XY', 'XZ', 'YX', 'YZ', 'ZX' or 'ZY',
            //    specifying the shear equation:
            //
            //      XY:  x' = x + s * y;
            //      XZ:  x' = x + s * z;
            //      YX:  y' = y + s * x;
            //      YZ:  y' = y + s * z;
            //      ZX:  z' = z + s * x;
            //      ZY:  z' = z + s * y.
            //
            //    Input, double S, the shear coefficient.
            //
        {
            double[] c = new double[4 * 4];
            double[] d = new double[4 * 4];
            int i;
            int j;

            tmat_init(c);

            if (axis == "XY" || axis == "xy")
            {
                c[0 + 1 * 4] = s;
            }
            else if (axis == "XZ" || axis == "xz")
            {
                c[0 + 2 * 4] = s;
            }
            else if (axis == "YX" || axis == "yx")
            {
                c[1 + 0 * 4] = s;
            }
            else if (axis == "YZ" || axis == "yz")
            {
                c[1 + 2 * 4] = s;
            }
            else if (axis == "ZX" || axis == "zx")
            {
                c[2 + 0 * 4] = s;
            }
            else if (axis == "ZY" || axis == "zy")
            {
                c[2 + 1 * 4] = s;
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("TMAT_SHEAR - Fatal error!");
                Console.WriteLine("  Illegal shear axis: " + axis + "");
                Console.WriteLine("  Legal choices are XY, XZ, YX, YZ, ZX, or ZY.");
                return;
            }

            tmat_mxm(c, a, ref d);

            for (i = 0; i < 4; i++)
            {
                for (j = 0; j < 4; j++)
                {
                    b[i + j * 4] = d[i + j * 4];
                }
            }
        }

        public static void tmat_trans(double[] a, ref double[] b, double[] v)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TMAT_TRANS applies a translation to the geometric transformation matrix.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 July 2005
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
            //    Input, double A[4*4], the current geometric transformation matrix.
            //
            //    Output, double B[4*4], the modified transformation matrix.
            //    A and B may share the same memory.
            //
            //    Input, double V[3], the translation.  This may be thought of as the
            //    point that the origin moves to under the translation.
            //
        {
            int i;
            int j;

            for (i = 0; i < 4; i++)
            {
                for (j = 0; j < 4; j++)
                {
                    b[i + j * 4] = a[i + j * 4];
                }
            }

            b[0 + 3 * 4] = b[0 + 3 * 4] + v[0];
            b[1 + 3 * 4] = b[1 + 3 * 4] + v[1];
            b[2 + 3 * 4] = b[2 + 3 * 4] + v[2];

        }
    }
}