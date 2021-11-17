using System;
using Burkardt.Geometry;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.TriangleNS;

public static class Geometry
{
    public static void triangle_angles_2d(double[] t, ref double[] angle)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_ANGLES_2D computes the angles of a triangle in 2D.
        //
        //  Discussion:
        //
        //    The law of cosines is used:
        //
        //      C * C = A * A + B * B - 2 * A * B * COS ( GAMMA )
        //
        //    where GAMMA is the angle opposite side C.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double T[2*3], the triangle vertices.
        //
        //    Output, double ANGLE[3], the angles opposite
        //    sides P1-P2, P2-P3 and P3-P1, in radians.
        //
    {
        double a;
        double b;
        double c;

        a = Math.Sqrt(Math.Pow(t[0 + 1 * 2] - t[0 + 0 * 2], 2)
                      + Math.Pow(t[1 + 1 * 2] - t[1 + 0 * 2], 2));

        b = Math.Sqrt(Math.Pow(t[0 + 2 * 2] - t[0 + 1 * 2], 2)
                      + Math.Pow(t[1 + 2 * 2] - t[1 + 1 * 2], 2));

        c = Math.Sqrt(Math.Pow(t[0 + 0 * 2] - t[0 + 2 * 2], 2)
                      + Math.Pow(t[1 + 0 * 2] - t[1 + 2 * 2], 2));
        switch (a)
        {
            //
            //  Take care of a ridiculous special case.
            //
            case 0.0 when b == 0.0 && c == 0.0:
                angle[0] = 2.0 * Math.PI / 3.0;
                angle[1] = 2.0 * Math.PI / 3.0;
                angle[2] = 2.0 * Math.PI / 3.0;
                return;
        }

        if (c == 0.0 || a == 0.0)
        {
            angle[0] = Math.PI;
        }
        else
        {
            angle[0] = typeMethods.r8_acos((c * c + a * a - b * b) / (2.0 * c * a));
        }

        if (a == 0.0 || b == 0.0)
        {
            angle[1] = Math.PI;
        }
        else
        {
            angle[1] = typeMethods.r8_acos((a * a + b * b - c * c) / (2.0 * a * b));
        }

        if (b == 0.0 || c == 0.0)
        {
            angle[2] = Math.PI;
        }
        else
        {
            angle[2] = typeMethods.r8_acos((b * b + c * c - a * a) / (2.0 * b * c));
        }

    }

    public static double[] triangle_angles_2d_new(double[] t)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_ANGLES_2D_NEW computes the angles of a triangle in 2D.
        //
        //  Discussion:
        //
        //    The law of cosines is used:
        //
        //      C * C = A * A + B * B - 2 * A * B * COS ( GAMMA )
        //
        //    where GAMMA is the angle opposite side C.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 September 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double T[2*3], the triangle vertices.
        //
        //    Output, double TRIANGLE_ANGLES_2D_NEW[3], the angles opposite
        //    sides P1-P2, P2-P3 and P3-P1, in radians.
        //
    {
        double a;
        double[] angle;
        double b;
        double c;

        angle = new double[3];

        a = Math.Sqrt(Math.Pow(t[0 + 1 * 2] - t[0 + 0 * 2], 2)
                      + Math.Pow(t[1 + 1 * 2] - t[1 + 0 * 2], 2));

        b = Math.Sqrt(Math.Pow(t[0 + 2 * 2] - t[0 + 1 * 2], 2)
                      + Math.Pow(t[1 + 2 * 2] - t[1 + 1 * 2], 2));

        c = Math.Sqrt(Math.Pow(t[0 + 0 * 2] - t[0 + 2 * 2], 2)
                      + Math.Pow(t[1 + 0 * 2] - t[1 + 2 * 2], 2));
        switch (a)
        {
            //
            //  Take care of a ridiculous special case.
            //
            case 0.0 when b == 0.0 && c == 0.0:
                angle[0] = 2.0 * Math.PI / 3.0;
                angle[1] = 2.0 * Math.PI / 3.0;
                angle[2] = 2.0 * Math.PI / 3.0;
                return angle;
        }

        if (c == 0.0 || a == 0.0)
        {
            angle[0] = Math.PI;
        }
        else
        {
            angle[0] = typeMethods.r8_acos((c * c + a * a - b * b) / (2.0 * c * a));
        }

        if (a == 0.0 || b == 0.0)
        {
            angle[1] = Math.PI;
        }
        else
        {
            angle[1] = typeMethods.r8_acos((a * a + b * b - c * c) / (2.0 * a * b));
        }

        if (b == 0.0 || c == 0.0)
        {
            angle[2] = Math.PI;
        }
        else
        {
            angle[2] = typeMethods.r8_acos((b * b + c * c - a * a) / (2.0 * b * c));
        }

        return angle;
    }

    public static void triangle_angles_3d(double[] t, ref double[] angle, int angleIndex = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_ANGLES_3D computes the angles of a triangle in 3D.
        //
        //  Discussion:
        //
        //    The law of cosines is used:
        //
        //      C * C = A * A + B * B - 2 * A * B * COS ( GAMMA )
        //
        //    where GAMMA is the angle opposite side C.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double T[3*3], the triangle vertices.
        //
        //    Output, double ANGLE[3], the angles opposite
        //    sides P1-P2, P2-P3 and P3-P1, in radians.
        //
    {
        int DIM_NUM = 3;

        double a;
        double b;
        double c;

        a = Math.Sqrt(Math.Pow(t[0 + 1 * DIM_NUM] - t[0 + 0 * DIM_NUM], 2)
                      + Math.Pow(t[1 + 1 * DIM_NUM] - t[1 + 0 * DIM_NUM], 2)
                      + Math.Pow(t[2 + 1 * DIM_NUM] - t[2 + 0 * DIM_NUM], 2));

        b = Math.Sqrt(Math.Pow(t[0 + 2 * DIM_NUM] - t[0 + 1 * DIM_NUM], 2)
                      + Math.Pow(t[1 + 2 * DIM_NUM] - t[1 + 1 * DIM_NUM], 2)
                      + Math.Pow(t[2 + 2 * DIM_NUM] - t[2 + 1 * DIM_NUM], 2));

        c = Math.Sqrt(Math.Pow(t[0 + 0 * DIM_NUM] - t[0 + 2 * DIM_NUM], 2)
                      + Math.Pow(t[1 + 0 * DIM_NUM] - t[1 + 2 * DIM_NUM], 2)
                      + Math.Pow(t[2 + 0 * DIM_NUM] - t[2 + 2 * DIM_NUM], 2));
        switch (a)
        {
            //
            //  Take care of a ridiculous special case.
            //
            case 0.0 when b == 0.0 && c == 0.0:
                angle[(0 + angleIndex) % angle.Length] = 2.0 * Math.PI / 3.0;
                angle[(1 + angleIndex) % angle.Length] = 2.0 * Math.PI / 3.0;
                angle[(2 + angleIndex) % angle.Length] = 2.0 * Math.PI / 3.0;
                return;
        }

        if (c == 0.0 || a == 0.0)
        {
            angle[(0 + angleIndex) % angle.Length] = Math.PI;
        }
        else
        {
            angle[(0 + angleIndex) % angle.Length] = typeMethods.r8_acos((c * c + a * a - b * b) / (2.0 * c * a));
        }

        if (a == 0.0 || b == 0.0)
        {
            angle[(1 + angleIndex) % angle.Length] = Math.PI;
        }
        else
        {
            angle[(1 + angleIndex) % angle.Length] = typeMethods.r8_acos((a * a + b * b - c * c) / (2.0 * a * b));
        }

        if (b == 0.0 || c == 0.0)
        {
            angle[(2 + angleIndex) % angle.Length] = Math.PI;
        }
        else
        {
            angle[(2 + angleIndex) % angle.Length] = typeMethods.r8_acos((b * b + c * c - a * a) / (2.0 * b * c));
        }
            
    }

    public static double[] triangle_angles_3d_new(double[] t)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_ANGLES_3D_NEW computes the angles of a triangle in 3D.
        //
        //  Discussion:
        //
        //    The law of cosines is used:
        //
        //      C * C = A * A + B * B - 2 * A * B * COS ( GAMMA )
        //
        //    where GAMMA is the angle opposite side C.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 September 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double T[3*3], the triangle vertices.
        //
        //    Output, double TRIANGLE_ANGLES_3D_NEW[3], the angles opposite
        //    sides P1-P2, P2-P3 and P3-P1, in radians.
        //
    {
        int DIM_NUM = 3;

        double a;
        double[] angle;
        double b;
        double c;

        angle = new double[3];

        a = Math.Sqrt(Math.Pow(t[0 + 1 * DIM_NUM] - t[0 + 0 * DIM_NUM], 2)
                      + Math.Pow(t[1 + 1 * DIM_NUM] - t[1 + 0 * DIM_NUM], 2)
                      + Math.Pow(t[2 + 1 * DIM_NUM] - t[2 + 0 * DIM_NUM], 2));

        b = Math.Sqrt(Math.Pow(t[0 + 2 * DIM_NUM] - t[0 + 1 * DIM_NUM], 2)
                      + Math.Pow(t[1 + 2 * DIM_NUM] - t[1 + 1 * DIM_NUM], 2)
                      + Math.Pow(t[2 + 2 * DIM_NUM] - t[2 + 1 * DIM_NUM], 2));

        c = Math.Sqrt(Math.Pow(t[0 + 0 * DIM_NUM] - t[0 + 2 * DIM_NUM], 2)
                      + Math.Pow(t[1 + 0 * DIM_NUM] - t[1 + 2 * DIM_NUM], 2)
                      + Math.Pow(t[2 + 0 * DIM_NUM] - t[2 + 2 * DIM_NUM], 2));
        switch (a)
        {
            //
            //  Take care of a ridiculous special case.
            //
            case 0.0 when b == 0.0 && c == 0.0:
                angle[0] = 2.0 * Math.PI / 3.0;
                angle[1] = 2.0 * Math.PI / 3.0;
                angle[2] = 2.0 * Math.PI / 3.0;
                return angle;
        }

        if (c == 0.0 || a == 0.0)
        {
            angle[0] = Math.PI;
        }
        else
        {
            angle[0] = typeMethods.r8_acos((c * c + a * a - b * b) / (2.0 * c * a));
        }

        if (a == 0.0 || b == 0.0)
        {
            angle[1] = Math.PI;
        }
        else
        {
            angle[1] = typeMethods.r8_acos((a * a + b * b - c * c) / (2.0 * a * b));
        }

        if (b == 0.0 || c == 0.0)
        {
            angle[2] = Math.PI;
        }
        else
        {
            angle[2] = typeMethods.r8_acos((b * b + c * c - a * a) / (2.0 * b * c));
        }

        return angle;

    }

    public static double triangle_area_2d(double[] t)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_AREA_2D computes the area of a triangle in 2D.
        //
        //  Discussion:
        //
        //    If the triangle's vertices are given in counter clockwise order,
        //    the area will be positive.  If the triangle's vertices are given
        //    in clockwise order, the area will be negative!
        //
        //    An earlier version of this routine always returned the absolute
        //    value of the computed area.  I am convinced now that that is
        //    a less useful result!  For instance, by returning the signed
        //    area of a triangle, it is possible to easily compute the area
        //    of a nonconvex polygon as the sum of the (possibly negative)
        //    areas of triangles formed by node 1 and successive pairs of vertices.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double T[2*3], the vertices of the triangle.
        //
        //    Output, double TRIANGLE_AREA_2D, the area of the triangle.
        //
    {
        double area;

        area = 0.5 * (
            t[0 + 0 * 2] * (t[1 + 1 * 2] - t[1 + 2 * 2]) +
            t[0 + 1 * 2] * (t[1 + 2 * 2] - t[1 + 0 * 2]) +
            t[0 + 2 * 2] * (t[1 + 0 * 2] - t[1 + 1 * 2]));

        return area;
    }

    public static double triangle_area_3d(double[] t)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_AREA_3D computes the area of a triangle in 3D.
        //
        //  Discussion:
        //
        //    This routine uses the fact that the norm of the cross product vector
        //    is the area of the parallelogram they form.  The triangle they
        //    form has half that area.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 October 2005
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
        //    Input, double T[3*3], the vertices of the triangle.
        //
        //    Output, double TRIANGLE_AREA_3D, the area of the triangle.
        //
    {
        int DIM_NUM = 3;

        double area;
        double[] cross;
        int i;
        //
        //  Compute the cross product vector.
        //
        cross = new double[DIM_NUM];

        cross[0] = (t[1 + 1 * DIM_NUM] - t[1 + 0 * DIM_NUM])
                   * (t[2 + 2 * DIM_NUM] - t[2 + 0 * DIM_NUM])
                   - (t[2 + 1 * DIM_NUM] - t[2 + 0 * DIM_NUM])
                   * (t[1 + 2 * DIM_NUM] - t[1 + 0 * DIM_NUM]);

        cross[1] = (t[2 + 1 * DIM_NUM] - t[2 + 0 * DIM_NUM])
                   * (t[0 + 2 * DIM_NUM] - t[0 + 0 * DIM_NUM])
                   - (t[0 + 1 * DIM_NUM] - t[0 + 0 * DIM_NUM])
                   * (t[2 + 2 * DIM_NUM] - t[2 + 0 * DIM_NUM]);

        cross[2] = (t[0 + 1 * DIM_NUM] - t[0 + 0 * DIM_NUM])
                   * (t[1 + 2 * DIM_NUM] - t[1 + 0 * DIM_NUM])
                   - (t[1 + 1 * DIM_NUM] - t[1 + 0 * DIM_NUM])
                   * (t[0 + 2 * DIM_NUM] - t[0 + 0 * DIM_NUM]);

        area = 0.0;
        for (i = 0; i < DIM_NUM; i++)
        {
            area += Math.Pow(cross[i], 2);
        }

        area = 0.5 * Math.Sqrt(area);

        return area;

    }

    public static double triangle_area_3d_2(double[] t)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_AREA_3D_2 computes the area of a triangle in 3D.
        //
        //  Discussion:
        //
        //    This routine computes the area "the hard way".
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double T[3*3], the vertices of the triangle.
        //
        //    Output, double TRIANGLE_AREA_3D_2, the area of the triangle.
        //
    {
        int DIM_NUM = 3;

        double alpha;
        double area;
        double base_;
        double dot;
        double height;
        double[] ph = new double[DIM_NUM];
        //
        //  Find the projection of (P3-P1) onto (P2-P1).
        //
        dot = (t[0 + 1 * DIM_NUM] - t[0 + 0 * DIM_NUM])
              * (t[0 + 2 * DIM_NUM] - t[0 + 0 * DIM_NUM])
              + (t[1 + 1 * DIM_NUM] - t[1 + 0 * DIM_NUM])
              * (t[1 + 2 * DIM_NUM] - t[1 + 0 * DIM_NUM])
              + (t[2 + 1 * DIM_NUM] - t[2 + 0 * DIM_NUM])
              * (t[2 + 2 * DIM_NUM] - t[2 + 0 * DIM_NUM]);

        base_ = Math.Sqrt(Math.Pow(t[0 + 1 * DIM_NUM] - t[0 + 0 * DIM_NUM], 2)
                          + Math.Pow(t[1 + 1 * DIM_NUM] - t[1 + 0 * DIM_NUM], 2)
                          + Math.Pow(t[2 + 1 * DIM_NUM] - t[2 + 0 * DIM_NUM], 2));
        switch (base_)
        {
            //
            //  The height of the triangle is the length of (P3-P1) after its
            //  projection onto (P2-P1) has been subtracted.
            //
            case 0.0:
                height = 0.0;
                break;
            default:
                alpha = dot / (base_ * base_);

                ph[0] = t[0 + 0 * DIM_NUM] + alpha * (t[0 + 1 * DIM_NUM] - t[0 + 0 * DIM_NUM]);
                ph[1] = t[1 + 0 * DIM_NUM] + alpha * (t[1 + 1 * DIM_NUM] - t[1 + 0 * DIM_NUM]);
                ph[2] = t[2 + 0 * DIM_NUM] + alpha * (t[2 + 1 * DIM_NUM] - t[2 + 0 * DIM_NUM]);

                height = Math.Sqrt(Math.Pow(ph[0] - t[0 + 2 * DIM_NUM], 2)
                                   + Math.Pow(ph[1] - t[1 + 2 * DIM_NUM], 2)
                                   + Math.Pow(ph[2] - t[2 + 2 * DIM_NUM], 2));
                break;
        }

        area = 0.5 * base_ * height;

        return area;

    }

    public static double triangle_area_3d_3(double[] t)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_AREA_3D_3 computes the area of a triangle in 3D.
        //
        //  Discussion:
        //
        //    This routine uses Heron's formula
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 July 2005
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
        //    Input, double T[3*3], the triangle vertices.
        //
        //    Output, double AREA, the area of the triangle.
        //
    {
        int DIM_NUM = 3;

        double area;
        int i;
        int j;
        int jp1;
        double[] s = new double[DIM_NUM];

        for (j = 0; j < DIM_NUM; j++)
        {
            jp1 = (j + 1) % DIM_NUM;
            s[j] = 0.0;
            for (i = 0; i < DIM_NUM; i++)
            {
                s[j] += Math.Pow(t[i + j * DIM_NUM] - t[i + jp1 * DIM_NUM], 2);
            }

            s[j] = Math.Sqrt(s[j]);
        }

        area = (s[0] + s[1] + s[2])
               * (-s[0] + s[1] + s[2])
               * (s[0] - s[1] + s[2])
               * (s[0] + s[1] - s[2]);

        switch (area)
        {
            case < 0.0:
                area = -1.0;
                return area;
            default:
                area = 0.25 * Math.Sqrt(area);

                return area;
        }
    }

    public static double triangle_area_heron(double[] s)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_AREA_HERON computes the area of a triangle using Heron's formula.
        //
        //  Discussion:
        //
        //    The formula is valid for any spatial dimension, depending only
        //    on the lengths of the sides, and not the coordinates of the vertices.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double S[3], the lengths of the three sides.
        //
        //    Output, double TRIANGLE_AREA_HERON, the area of the triangle, or -1.0 if the
        //    sides cannot constitute a triangle.
        //
    {
        double area;

        area = (s[0] + s[1] + s[2])
               * (-s[0] + s[1] + s[2])
               * (s[0] - s[1] + s[2])
               * (s[0] + s[1] - s[2]);

        switch (area)
        {
            case < 0.0:
                area = -1.0;
                return area;
            default:
                area = 0.25 * Math.Sqrt(area);

                return area;
        }
    }

    public static double[] triangle_area_vector_3d(double[] t)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_AREA_VECTOR_3D computes the area vector of a triangle in 3D.
        //
        //  Discussion:
        //
        //    The "area vector" of a triangle is simply a cross product of,
        //    for instance, the vectors (V2-V1) and (V3-V1), where V1, V2
        //    and V3 are the vertices of the triangle.
        //
        //    The norm of the cross product vector of two vectors is the area
        //    of the parallelogram they form.
        //
        //    Therefore, the area of the triangle is half of the norm of the
        //    area vector:
        //
        //      area = 0.5 * Math.Sqrt ( sum ( area_vector(1:3)^2 ) )
        //
        //    The reason for looking at the area vector rather than the area
        //    is that this makes it possible to compute the area of a flat
        //    polygon in 3D by summing the areas of the triangles that form
        //    a decomposition of the polygon, while allowing for both positive
        //    and negative areas.  (Sum the vectors, THEN take the norm and
        //    multiply by 1/2).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 October 2005
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
        //    Input, double T[3*3], the vertices of the triangle.
        //
        //    Output, double TRIANGLE_AREA_VECTOR_3D[3], the area vector of the triangle.
        //
    {
        int DIM_NUM = 3;

        double[] cross;
        //
        //  Compute the cross product vector.
        //
        cross = new double[DIM_NUM];

        cross[0] = (t[1 + 1 * DIM_NUM] - t[1 + 0 * DIM_NUM])
                   * (t[2 + 2 * DIM_NUM] - t[2 + 0 * DIM_NUM])
                   - (t[2 + 1 * DIM_NUM] - t[2 + 0 * DIM_NUM])
                   * (t[1 + 2 * DIM_NUM] - t[1 + 0 * DIM_NUM]);

        cross[1] = (t[2 + 1 * DIM_NUM] - t[2 + 0 * DIM_NUM])
                   * (t[0 + 2 * DIM_NUM] - t[0 + 0 * DIM_NUM])
                   - (t[0 + 1 * DIM_NUM] - t[0 + 0 * DIM_NUM])
                   * (t[2 + 2 * DIM_NUM] - t[2 + 0 * DIM_NUM]);

        cross[2] = (t[0 + 1 * DIM_NUM] - t[0 + 0 * DIM_NUM])
                   * (t[1 + 2 * DIM_NUM] - t[1 + 0 * DIM_NUM])
                   - (t[1 + 1 * DIM_NUM] - t[1 + 0 * DIM_NUM])
                   * (t[0 + 2 * DIM_NUM] - t[0 + 0 * DIM_NUM]);

        return cross;

    }

    public static double[] triangle_barycentric_2d(double[] t, double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_BARYCENTRIC_2D finds the barycentric coordinates of a point in 2D.
        //
        //  Discussion:
        //
        //    The barycentric coordinate of point X related to vertex A can be
        //    interpreted as the ratio of the area of the triangle with
        //    vertex A replaced by vertex X to the area of the original
        //    triangle.
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
        //  Parameters:
        //
        //    Input, double T[2*3], the vertices of the triangle.
        //
        //    Input, double P[2], the point to be checked.
        //
        //    Output, double C[3], the barycentric coordinates of the point with respect
        //    to the triangle.
        //
    {
        int N = 2;
        int RHS_NUM = 1;

        double[] a = new double[N * (N + RHS_NUM)];
        double[] c;
        int info;
        //
        //  Set up the linear system
        //
        //    ( X2-X1  X3-X1 ) C1  = X-X1
        //    ( Y2-Y1  Y3-Y1 ) C2    Y-Y1
        //
        //  which is satisfied by the barycentric coordinates.
        //
        a[0 + 0 * N] = t[0 + 1 * 2] - t[0 + 0 * 2];
        a[1 + 0 * N] = t[1 + 1 * 2] - t[1 + 0 * 2];

        a[0 + 1 * N] = t[0 + 2 * 2] - t[0 + 0 * 2];
        a[1 + 1 * N] = t[1 + 2 * 2] - t[1 + 0 * 2];

        a[0 + 2 * N] = p[0] - t[0 + 0 * 2];
        a[1 + 2 * N] = p[1] - t[1 + 0 * 2];
        //
        //  Solve the linear system.
        //
        info = typeMethods.r8mat_solve(N, RHS_NUM, ref a);

        if (info != 0)
        {
            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_BARYCENTRIC_2D - Fatal error!");
            Console.WriteLine("  The linear system is singular.");
            Console.WriteLine("  The input data does not form a proper triangle.");
            return null;
        }

        c = new double[3];

        c[0] = a[0 + 2 * N];
        c[1] = a[1 + 2 * N];
        c[2] = 1.0 - c[0] - c[1];

        return c;
    }

    public static double[] triangle_centroid_2d(double[] t)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_CENTROID_2D computes the centroid of a triangle in 2D.
        //
        //  Discussion:
        //
        //    The centroid of a triangle can also be considered the center
        //    of gravity, assuming that the triangle is made of a thin uniform
        //    sheet of massy material.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 July 2005
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
        //    Input, double T[2*3], the vertices of the triangle.
        //
        //    Output, double TRIANGLE_CENTROID_2D[2], the coordinates of the centroid of the triangle.
        //
    {
        int DIM_NUM = 2;

        double[] centroid;

        centroid = new double[DIM_NUM];

        centroid[0] = (t[0 + 0 * DIM_NUM] + t[0 + 1 * DIM_NUM] + t[0 + 2 * DIM_NUM]) / 3.0;
        centroid[1] = (t[1 + 0 * DIM_NUM] + t[1 + 1 * DIM_NUM] + t[1 + 2 * DIM_NUM]) / 3.0;

        return centroid;

    }

    public static double[] triangle_centroid_3d(double[] t)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_CENTROID_3D computes the centroid of a triangle in 3D.
        //
        //  Discussion:
        //
        //    The centroid of a triangle can also be considered the center
        //    of gravity, assuming that the triangle is made of a thin uniform
        //    sheet of massy material.
        //
        //    Thanks to Gordon Griesel for pointing out a typographical
        //    error in an earlier version of this program, and for pointing
        //    out a second oversight, as well.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double T[3*3], the vertices of the triangle.
        //
        //    Output, double TRIANGLE_CENTROID_3D[3], the coordinates of the centroid.
        //
    {
        int DIM_NUM = 3;

        double[] centroid;

        centroid = new double[DIM_NUM];

        centroid[0] = (t[0 + 0 * DIM_NUM] + t[0 + 1 * DIM_NUM] + t[0 + 2 * DIM_NUM]) / 3.0;
        centroid[1] = (t[1 + 0 * DIM_NUM] + t[1 + 1 * DIM_NUM] + t[1 + 2 * DIM_NUM]) / 3.0;
        centroid[2] = (t[2 + 0 * DIM_NUM] + t[2 + 1 * DIM_NUM] + t[2 + 2 * DIM_NUM]) / 3.0;

        return centroid;

    }

    public static double[] triangle_circumcenter_2d(double[] t)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_CIRCUMCENTER_2D computes the circumcenter of a triangle in 2D.
        //
        //  Discussion:
        //
        //    The circumcenter of a triangle is the center of the circumcircle, the
        //    circle that passes through the three vertices of the triangle.
        //
        //    The circumcircle contains the triangle, but it is not necessarily the
        //    smallest triangle to do so.
        //
        //    If all angles of the triangle are no greater than 90 degrees, then
        //    the center of the circumscribed circle will lie inside the triangle.
        //    Otherwise, the center will lie outside the triangle.
        //
        //    The circumcenter is the intersection of the perpendicular bisectors
        //    of the sides of the triangle.
        //
        //    In geometry, the circumcenter of a triangle is often symbolized by "O".
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 February 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double T[2*3], the triangle vertices.
        //
        //    Output, double *TRIANGLE_CIRCUMCENTER_2D[2], the circumcenter.
        //
    {
        int DIM_NUM = 2;

        double asq;
        double bot;
        double[] pc;
        double csq;
        double top1;
        double top2;

        pc = new double[DIM_NUM];

        asq = (t[0 + 1 * 2] - t[0 + 0 * 2]) * (t[0 + 1 * 2] - t[0 + 0 * 2])
              + (t[1 + 1 * 2] - t[1 + 0 * 2]) * (t[1 + 1 * 2] - t[1 + 0 * 2]);

        csq = (t[0 + 2 * 2] - t[0 + 0 * 2]) * (t[0 + 2 * 2] - t[0 + 0 * 2])
              + (t[1 + 2 * 2] - t[1 + 0 * 2]) * (t[1 + 2 * 2] - t[1 + 0 * 2]);

        top1 = (t[1 + 1 * 2] - t[1 + 0 * 2]) * csq - (t[1 + 2 * 2] - t[1 + 0 * 2]) * asq;
        top2 = -(t[0 + 1 * 2] - t[0 + 0 * 2]) * csq + (t[0 + 2 * 2] - t[0 + 0 * 2]) * asq;

        bot = (t[1 + 1 * 2] - t[1 + 0 * 2]) * (t[0 + 2 * 2] - t[0 + 0 * 2])
              - (t[1 + 2 * 2] - t[1 + 0 * 2]) * (t[0 + 1 * 2] - t[0 + 0 * 2]);

        pc[0] = t[0 + 0 * 2] + 0.5 * top1 / bot;
        pc[1] = t[1 + 0 * 2] + 0.5 * top2 / bot;

        return pc;


    }

    public static double[] triangle_circumcenter_2d_2(double[] t)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_CIRCUMCENTER_2D_2 computes the circumcenter of a triangle in 2D.
        //
        //  Discussion:
        //
        //    The circumcenter of a triangle is the center of the circumcircle, the
        //    circle that passes through the three vertices of the triangle.
        //
        //    The circumcircle contains the triangle, but it is not necessarily the
        //    smallest triangle to do so.
        //
        //    If all angles of the triangle are no greater than 90 degrees, then
        //    the center of the circumscribed circle will lie inside the triangle.
        //    Otherwise, the center will lie outside the triangle.
        //
        //    The circumcenter is the intersection of the perpendicular bisectors
        //    of the sides of the triangle.
        //
        //    Surprisingly, the diameter of the circle can be found by solving
        //    a 2 by 2 linear system.  If we label the vertices of the triangle
        //    P1, P2 and P3, then the vectors P2 - P1 and P3 - P1 are secants of
        //    the circle, and each forms a right triangle with the diameter.
        //    Hence, the dot product of P2 - P1 with the diameter vector is equal
        //    to the square of the length of P2 - P1, and similarly for P3 - P1.
        //    This determines the diameter vector originating at P1.
        //
        //    In geometry, the circumcenter of a triangle is often symbolized by "O".
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 August 2003
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
        //    Input, double T[2*3], the triangle vertices.
        //
        //    Output, double *TRIANGLE_CIRCUMCENTER_2D[2], the circumcenter.
        //
    {
        int N = 2;
        int RHS_NUM = 1;

        double[] a = new double[N * (N + RHS_NUM)];
        double[] pc;
        int info;
        //
        //  Set up the linear system.
        //
        a[0 + 0 * N] = t[0 + 1 * 2] - t[0 + 0 * 2];
        a[0 + 1 * N] = t[1 + 1 * 2] - t[1 + 0 * 2];
        a[0 + 2 * N] = Math.Pow(t[0 + 1 * 2] - t[0 + 0 * 2], 2)
                       + Math.Pow(t[1 + 1 * 2] - t[1 + 0 * 2], 2);

        a[1 + 0 * N] = t[0 + 2 * 2] - t[0 + 0 * 2];
        a[1 + 1 * N] = t[1 + 2 * 2] - t[1 + 0 * 2];
        a[1 + 2 * N] = Math.Pow(t[0 + 2 * 2] - t[0 + 0 * 2], 2)
                       + Math.Pow(t[1 + 2 * 2] - t[1 + 0 * 2], 2);
        //
        //  Solve the linear system.
        //
        info = typeMethods.r8mat_solve(N, RHS_NUM, ref a);
        //
        //  Compute the center.
        //
        pc = new double[2];

        if (info != 0)
        {
            pc[0] = 0.0;
            pc[1] = 0.0;
        }
        else
        {
            pc[0] = t[0 + 0 * 2] + 0.5 * a[0 + N * N];
            pc[1] = t[1 + 0 * 2] + 0.5 * a[1 + N * N];
        }

        return pc;
    }

    public static double[] triangle_circumcenter(int n, double[] t)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_CIRCUMCENTER computes the circumcenter of a triangle in ND.
        //
        //  Discussion:
        //
        //    Three ND points A, B and C lie on a circle.
        //
        //    The circumcenter P has the formula
        //
        //      P = ( Area ( PBC ) * A + Area ( APC) * B + Area ( ABP ) * C )
        //        / ( Area ( PBC )     + Area ( APC )    + Area ( ABP ) )
        //
        //    The details of the formula rely on information supplied
        //    by Oscar Lanzi III.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 October 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the spatial dimension.
        //
        //    Input, double T[N*3], the triangle vertices.
        //
        //    Output, double TRIANGLE_CIRCUMCENTER[N], the circumcenter.
        //
    {
        double a;
        double abp;
        double apc;
        double b;
        double c;
        int i;
        double[] p;
        double pbc;

        a = typeMethods.r8vec_normsq_affine(n, t, t, +1 * n, +2 * n);
        b = typeMethods.r8vec_normsq_affine(n, t, t, +2 * n, +0 * n);
        c = typeMethods.r8vec_normsq_affine(n, t, t, +0 * n, +1 * n);

        pbc = a * (-a + b + c);
        apc = b * (a - b + c);
        abp = c * (a + b - c);

        p = new double[n];

        for (i = 0; i < n; i++)
        {
            p[i] = (pbc * t[i + 0 * n] + apc * t[i + 1 * n] + abp * t[i + 2 * n])
                   / (pbc + apc + abp);
        }

        return p;
    }

    public static void triangle_circumcircle_2d(double[] t, ref double r, ref double[] pc)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_CIRCUMCIRCLE_2D computes the circumcircle of a triangle in 2D.
        //
        //  Discussion:
        //
        //    The circumcenter of a triangle is the center of the circumcircle, the
        //    circle that passes through the three vertices of the triangle.
        //
        //    The circumcircle contains the triangle, but it is not necessarily the
        //    smallest triangle to do so.
        //
        //    If all angles of the triangle are no greater than 90 degrees, then
        //    the center of the circumscribed circle will lie inside the triangle.
        //    Otherwise, the center will lie outside the triangle.
        //
        //    The circumcenter is the intersection of the perpendicular bisectors
        //    of the sides of the triangle.
        //
        //    In geometry, the circumcenter of a triangle is often symbolized by "O".
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
        //  Parameters:
        //
        //    Input, double T[2*3], the triangle vertices.
        //
        //    Output, double *R, PC[2], the circumradius, and the coordinates of the
        //    circumcenter of the triangle.
        //
    {
        double a;
        double b;
        double bot;
        double c;
        double top1;
        double top2;
        //
        //  Circumradius.
        //
        a = Math.Sqrt(Math.Pow(t[0 + 1 * 2] - t[0 + 0 * 2], 2) + Math.Pow(t[1 + 1 * 2] - t[1 + 0 * 2], 2));
        b = Math.Sqrt(Math.Pow(t[0 + 2 * 2] - t[0 + 1 * 2], 2) + Math.Pow(t[1 + 2 * 2] - t[1 + 1 * 2], 2));
        c = Math.Sqrt(Math.Pow(t[0 + 0 * 2] - t[0 + 2 * 2], 2) + Math.Pow(t[1 + 0 * 2] - t[1 + 2 * 2], 2));

        bot = (a + b + c) * (-a + b + c) * (a - b + c) * (a + b - c);

        switch (bot)
        {
            case <= 0.0:
                r = -1.0;
                pc[0] = 0.0;
                pc[1] = 0.0;
                return;
        }

        r = a * b * c / Math.Sqrt(bot);
        //
        //  Circumcenter.
        //
        top1 = (t[1 + 1 * 2] - t[1 + 0 * 2]) * c * c - (t[1 + 2 * 2] - t[1 + 0 * 2]) * a * a;
        top2 = (t[0 + 1 * 2] - t[0 + 0 * 2]) * c * c - (t[0 + 2 * 2] - t[0 + 0 * 2]) * a * a;
        bot = (t[1 + 1 * 2] - t[1 + 0 * 2]) * (t[0 + 2 * 2] - t[0 + 0 * 2])
              - (t[1 + 2 * 2] - t[1 + 0 * 2]) * (t[0 + 1 * 2] - t[0 + 0 * 2]);

        pc[0] = t[0 + 0 * 2] + 0.5 * top1 / bot;
        pc[1] = t[1 + 0 * 2] - 0.5 * top2 / bot;

    }

    public static void triangle_circumcircle_2d_2(double[] t, ref double r, ref double[] pc)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_CIRCUMCIRCLE_2D_2 computes the circumcircle of a triangle in 2D.
        //
        //  Discussion:
        //
        //    The circumscribed circle of a triangle is the circle that passes through
        //    the three vertices of the triangle.  The circumscribed circle contains
        //    the triangle, but it is not necessarily the smallest triangle to do so.
        //
        //    Surprisingly, the diameter of the circle can be found by solving
        //    a 2 by 2 linear system.  This is because the vectors P2 - P1
        //    and P3 - P1 are secants of the circle, and each forms a right
        //    triangle with the diameter.  Hence, the dot product of
        //    P2 - P1 with the diameter is equal to the square of the length
        //    of P2 - P1, and similarly for P3 - P1.  This determines the
        //    diameter vector originating at P1.
        //
        //    If all angles of the triangle are no greater than 90 degrees, then
        //    the center of the circumscribed circle will lie inside the triangle.
        //    Otherwise, the center will lie outside the triangle.
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
        //    Adrian Bowyer, John Woodwark,
        //    A Programmer's Geometry,
        //    Butterworths, 1983.
        //
        //  Parameters:
        //
        //    Input, double T[2*3], the triangle vertices.
        //
        //    Output, double *R, PC[2], the radius and coordinates of the center of the
        //    circumscribed circle.  If the linear system is
        //    singular, then R = -1, PC = 0.
        //
    {
        int N = 2;
        int RHS_NUM = 1;

        double[] a = new double[N * (N + RHS_NUM)];
        int info;
        //
        //  Set up the linear system.
        //
        a[0 + 0 * N] = t[0 + 1 * 2] - t[0 + 0 * 2];
        a[1 + 0 * N] = t[0 + 2 * 2] - t[0 + 0 * 2];

        a[0 + 1 * N] = t[1 + 1 * 2] - t[1 + 0 * 2];
        a[1 + 1 * N] = t[1 + 2 * 2] - t[1 + 0 * 2];

        a[0 + 2 * N] = Math.Pow(t[0 + 1 * 2] - t[0 + 0 * 2], 2) + Math.Pow(t[1 + 1 * 2] - t[1 + 0 * 2], 2);
        a[1 + 2 * N] = Math.Pow(t[0 + 2 * 2] - t[0 + 0 * 2], 2) + Math.Pow(t[1 + 2 * 2] - t[1 + 0 * 2], 2);
        //
        //  Solve the linear system.
        //
        info = typeMethods.r8mat_solve(N, RHS_NUM, ref a);
        //
        //  If the system was singular, return a consolation prize.
        //
        if (info != 0)
        {
            r = -1.0;
            pc[0] = 0.0;
            pc[1] = 0.0;
            return;
        }

        //
        //  Compute the radius and center.
        //
        r = 0.5 * Math.Sqrt(a[0 + N * N] * a[0 + N * N] + a[1 + N * N] * a[1 + N * N]);
        pc[0] = t[0 + 0 * 2] + 0.5 * a[0 + N * N];
        pc[1] = t[1 + 0 * 2] + 0.5 * a[1 + N * N];

    }

    public static double triangle_circumradius_2d(double[] t)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_CIRCUMRADIUS_2D computes the circumradius of a triangle in 2D.
        //
        //  Discussion:
        //
        //    The circumscribed circle of a triangle is the circle that passes through
        //    the three vertices of the triangle.  The circumscribed circle contains
        //    the triangle, but it is not necessarily the smallest triangle to do so.
        //
        //    The circumradius of a triangle is the radius of the circumscribed
        //    circle.
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
        //  Parameters:
        //
        //    Input, double T[2*3], the triangle vertices.
        //
        //    Output, double TRIANGLE_CIRCUMRADIUS_2D, the circumradius of the
        //    circumscribed circle.
        //
    {
        double a;
        double b;
        double bot;
        double c;
        double r;

        a = Math.Sqrt(Math.Pow(t[0 + 1 * 2] - t[0 + 0 * 2], 2) + Math.Pow(t[1 + 1 * 2] - t[1 + 0 * 2], 2));
        b = Math.Sqrt(Math.Pow(t[0 + 2 * 2] - t[0 + 1 * 2], 2) + Math.Pow(t[1 + 2 * 2] - t[1 + 1 * 2], 2));
        c = Math.Sqrt(Math.Pow(t[0 + 0 * 2] - t[0 + 2 * 2], 2) + Math.Pow(t[1 + 0 * 2] - t[1 + 2 * 2], 2));

        bot = (a + b + c) * (-a + b + c) * (a - b + c) * (a + b - c);

        switch (bot)
        {
            case <= 0.0:
                r = -1.0;
                return r;
            default:
                r = a * b * c / Math.Sqrt(bot);

                return r;
        }
    }

    public static void triangle_contains_line_exp_3d(double[] t, double[] p1,
            double[] p2, ref bool inside, ref double[] pint)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_CONTAINS_LINE_EXP_3D finds if a line is inside a triangle in 3D.
        //
        //  Discussion:
        //
        //    A line will "intersect" the plane of a triangle in 3D if
        //    * the line does not lie in the plane of the triangle
        //      (there would be infinitely many intersections), AND
        //    * the line does not lie parallel to the plane of the triangle
        //      (there are no intersections at all).
        //
        //    Therefore, if a line intersects the plane of a triangle, it does so
        //    at a single point.  We say the line is "inside" the triangle if,
        //    regarded as 2D objects, the intersection point of the line and the plane
        //    is inside the triangle.
        //
        //    The explicit form of a line in 3D is:
        //
        //      the line through the points P1, P2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Steve Marschner, Cornell University,
        //    CS465 Notes: Simple Ray-Triangle Intersection
        //
        //  Parameters:
        //
        //    Input, double T[3*3], the triangle vertices.
        //    The vertices should be given in counter clockwise order.
        //
        //    Input, double P1[3], P2[3], two points on the line.
        //
        //    Output, bool INSIDE, is TRUE if the line is inside the triangle.
        //
        //    Output, double PINT[3], the point where the line
        //    intersects the plane of the triangle.
        //
    {
        int DIM_NUM = 3;

        int i;
        int ival;
        double[] normal = new double[DIM_NUM];
        double[] normal2 = new double[DIM_NUM];
        double temp;
        double[] v1 = new double[DIM_NUM];
        double[] v2 = new double[DIM_NUM];
        //
        //  Make sure the line is not degenerate.
        //
        if (LineNS.Geometry.line_exp_is_degenerate_nd(DIM_NUM, p1, p2))
        {
            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_CONTAINS_LINE_EXP_3D - Fatal error!");
            Console.WriteLine("  The explicit line is degenerate.");
            return;
        }

        //
        //  Make sure the triangle is not degenerate.
        //
        if (triangle_is_degenerate_nd(DIM_NUM, t))
        {
            Console.WriteLine("");
            Console.WriteLine("TRIANGLE_CONTAINS_LINE_EXP_3D - Fatal error!");
            Console.WriteLine("  The triangle is degenerate.");
            return;
        }

        //
        //  Determine a unit normal vector associated with the plane of
        //  the triangle.
        //
        for (i = 0; i < DIM_NUM; i++)
        {
            v1[i] = t[i + 1 * DIM_NUM] - t[i + 0 * DIM_NUM];
        }

        for (i = 0; i < DIM_NUM; i++)
        {
            v2[i] = t[i + 2 * DIM_NUM] - t[i + 0 * DIM_NUM];
        }

        normal[0] = v1[1] * v2[2] - v1[2] * v2[1];
        normal[1] = v1[2] * v2[0] - v1[0] * v2[2];
        normal[2] = v1[0] * v2[1] - v1[1] * v2[0];

        temp = 0.0;
        for (i = 0; i < DIM_NUM; i++)
        {
            temp += Math.Pow(normal[i], 2);
        }

        temp = Math.Sqrt(temp);

        for (i = 0; i < DIM_NUM; i++)
        {
            normal[i] /= temp;
        }

        //
        //  Find the intersection of the plane and the line.
        //
        ival = Plane.Geometry.plane_normal_line_exp_int_3d(t, normal, p1, p2, ref pint);

        switch (ival)
        {
            case 0:
            {
                inside = false;
                for (i = 0; i < DIM_NUM; i++)
                {
                    pint[i] = typeMethods.r8_huge();
                }

                return;
            }
            case 2:
                inside = false;
                typeMethods.r8vec_copy(DIM_NUM, p1, ref pint);
                return;
        }

        //
        //  Now, check that all three triangles made by two vertices and
        //  the intersection point have the same "clock sense" as the
        //  triangle's normal vector.
        //
        for (i = 0; i < DIM_NUM; i++)
        {
            v1[i] = t[i + 1 * DIM_NUM] - t[i + 0 * DIM_NUM];
        }

        for (i = 0; i < DIM_NUM; i++)
        {
            v2[i] = pint[i] - t[i + 0 * DIM_NUM];
        }

        normal2[0] = v1[1] * v2[2] - v1[2] * v2[1];
        normal2[1] = v1[2] * v2[0] - v1[0] * v2[2];
        normal2[2] = v1[0] * v2[1] - v1[1] * v2[0];

        if (typeMethods.r8vec_dot_product(DIM_NUM, normal, normal2) < 0.0)
        {
            inside = false;
            return;
        }

        for (i = 0; i < DIM_NUM; i++)
        {
            v1[i] = t[i + 2 * DIM_NUM] - t[i + 1 * DIM_NUM];
        }

        for (i = 0; i < DIM_NUM; i++)
        {
            v2[i] = pint[i] - t[i + 1 * DIM_NUM];
        }

        normal2[0] = v1[1] * v2[2] - v1[2] * v2[1];
        normal2[1] = v1[2] * v2[0] - v1[0] * v2[2];
        normal2[2] = v1[0] * v2[1] - v1[1] * v2[0];

        if (typeMethods.r8vec_dot_product(DIM_NUM, normal, normal2) < 0.0)
        {
            inside = false;
            return;
        }

        for (i = 0; i < DIM_NUM; i++)
        {
            v1[i] = t[i + 0 * DIM_NUM] - t[i + 2 * DIM_NUM];
        }

        for (i = 0; i < DIM_NUM; i++)
        {
            v2[i] = pint[i] - t[i + 2 * DIM_NUM];
        }

        normal2[0] = v1[1] * v2[2] - v1[2] * v2[1];
        normal2[1] = v1[2] * v2[0] - v1[0] * v2[2];
        normal2[2] = v1[0] * v2[1] - v1[1] * v2[0];

        if (typeMethods.r8vec_dot_product(DIM_NUM, normal, normal2) < 0.0)
        {
            inside = false;
            return;
        }

        inside = true;

    }

    public static void triangle_contains_line_par_3d(double[] t, double[] p0, double[] pd,
            ref bool inside, ref double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_CONTAINS_LINE_PAR_3D: finds if a line is inside a triangle in 3D.
        //
        //  Discussion:
        //
        //    A line will "intersect" the plane of a triangle in 3D if
        //    * the line does not lie in the plane of the triangle
        //      (there would be infinitely many intersections), AND
        //    * the line does not lie parallel to the plane of the triangle
        //      (there are no intersections at all).
        //
        //    Therefore, if a line intersects the plane of a triangle, it does so
        //    at a single point.  We say the line is "inside" the triangle if,
        //    regarded as 2D objects, the intersection point of the line and the plane
        //    is inside the triangle.
        //
        //    A triangle in 3D is determined by three points:
        //
        //      T(1:3,1), T(1:3,2) and T(1:3,3).
        //
        //    The parametric form of a line in 3D is:
        //
        //      P(1:3) = P0(1:3) + PD(1:3) * T
        //
        //    We can normalize by requiring PD to have euclidean norm 1,
        //    and the first nonzero entry positive.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 February 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Adrian Bowyer, John Woodwark,
        //    A Programmer's Geometry,
        //    Butterworths, 1983, page 111.
        //
        //  Parameters:
        //
        //    Input, double T[3*3], the three points that define
        //    the triangle.
        //
        //    Input, double P0[3], PD(3], parameters that define the
        //    parametric line.
        //
        //    Output, bool *INSIDE, is TRUE if (the intersection point of)
        //    the line is inside the triangle.
        //
        //    Output, double P[3], is the point of intersection of the line
        //    and the plane of the triangle, unless they are parallel.
        //
    {
        int DIM_NUM = 3;

        double a;
        double angle_sum;
        double b;
        double c;
        double d;
        double denom;
        int dim;
        //bool intersect;
        double norm;
        double norm1;
        double norm2;
        double t_int;
        double tol = 0.00001;
        double[] v1 = new double[DIM_NUM];
        double[] v2 = new double[DIM_NUM];
        double[] v3 = new double[DIM_NUM];
        //
        //  Determine the implicit form (A,B,C,D) of the plane containing the
        //  triangle.
        //
        a = (t[1 + 1 * 3] - t[1 + 0 * 3]) * (t[2 + 2 * 3] - t[2 + 0 * 3])
            - (t[2 + 1 * 3] - t[2 + 0 * 3]) * (t[1 + 2 * 3] - t[1 + 0 * 3]);

        b = (t[2 + 1 * 3] - t[2 + 0 * 3]) * (t[0 + 2 * 3] - t[0 + 0 * 3])
            - (t[0 + 1 * 3] - t[0 + 0 * 3]) * (t[2 + 2 * 3] - t[2 + 0 * 3]);

        c = (t[0 + 1 * 3] - t[0 + 0 * 3]) * (t[1 + 2 * 3] - t[1 + 0 * 3])
            - (t[1 + 1 * 3] - t[1 + 0 * 3]) * (t[0 + 2 * 3] - t[0 + 0 * 3]);

        d = -t[0 + 1 * 3] * a - t[1 + 1 * 3] * b - t[2 + 1 * 3] * c;
        //
        //  Make sure the plane is well-defined.
        //
        norm1 = Math.Sqrt(a * a + b * b + c * c);

        switch (norm1)
        {
            case 0.0:
                Console.WriteLine("");
                Console.WriteLine("TRIANGLE_LINE_PAR_INT_3D - Fatal error!");
                Console.WriteLine("  The plane normal vector is null.");
                inside = false;
                typeMethods.r8vec_zero(DIM_NUM, ref p);
                return;
        }

        //
        //  Make sure the implicit line is well defined.
        //
        norm2 = typeMethods.r8vec_norm(DIM_NUM, pd);

        switch (norm2)
        {
            case 0.0:
                Console.WriteLine("");
                Console.WriteLine("TRIANGLE_LINE_PAR_INT_3D - Fatal error!");
                Console.WriteLine("  The line direction vector is null.");
                inside = false;
                typeMethods.r8vec_zero(DIM_NUM, ref p);
                return;
        }

        //
        //  Determine the denominator of the parameter in the
        //  implicit line definition that determines the intersection
        //  point.
        //
        denom = a * pd[0] + b * pd[1] + c * pd[2];
        //
        //  If DENOM is zero, or very small, the line and the plane may be
        //  parallel or almost so.
        //
        if (Math.Abs(denom) < tol * norm1 * norm2)
        {
            switch (a * p0[0] + b * p0[1] + c * p0[2] + d)
            {
                //  The line may actually lie in the plane.  We're not going
                //  to try to address this possibility.
                //
                case 0.0:
                    inside = false;
                    typeMethods.r8vec_copy(DIM_NUM, p0, ref p);
                    break;
                //
                default:
                    inside = false;
                    typeMethods.r8vec_zero(DIM_NUM, ref p);
                    break;
            }
        }
        //
        //  The line and plane intersect at a single point P.
        //
        else
        {
            t_int = -(a * p0[0] + b * p0[1] + c * p0[2] + d) / denom;
            for (dim = 0; dim < DIM_NUM; dim++)
            {
                p[dim] = p0[dim] + t_int * pd[dim];
            }

            //
            //  To see if P is included in the triangle, sum the angles
            //  formed by P and pairs of the vertices.  If the point is in the
            //  triangle, we get a total 360 degree view.  Otherwise, we
            //  get less than 180 degrees.
            //
            for (dim = 0; dim < DIM_NUM; dim++)
            {
                v1[dim] = t[dim + 0 * 3] - p[dim];
            }

            for (dim = 0; dim < DIM_NUM; dim++)
            {
                v2[dim] = t[dim + 1 * 3] - p[dim];
            }

            for (dim = 0; dim < DIM_NUM; dim++)
            {
                v3[dim] = t[dim + 2 * 3] - p[dim];
            }

            norm = typeMethods.r8vec_norm(DIM_NUM, v1);

            switch (norm)
            {
                case 0.0:
                    inside = true;
                    return;
            }

            for (dim = 0; dim < DIM_NUM; dim++)
            {
                v1[dim] /= norm;
            }

            norm = typeMethods.r8vec_norm(DIM_NUM, v2);

            switch (norm)
            {
                case 0.0:
                    inside = true;
                    return;
            }

            for (dim = 0; dim < DIM_NUM; dim++)
            {
                v2[dim] /= norm;
            }

            norm = typeMethods.r8vec_norm(DIM_NUM, v3);

            switch (norm)
            {
                case 0.0:
                    inside = true;
                    return;
            }

            for (dim = 0; dim < DIM_NUM; dim++)
            {
                v3[dim] /= norm;
            }

            angle_sum = typeMethods.r8_acos(typeMethods.r8vec_dot_product(DIM_NUM, v1, v2))
                        + typeMethods.r8_acos(typeMethods.r8vec_dot_product(DIM_NUM, v2, v3))
                        + typeMethods.r8_acos(typeMethods.r8vec_dot_product(DIM_NUM, v3, v1));

            if (typeMethods.r8_nint(angle_sum / Math.PI) == 2)
            {
                inside = true;
            }
            else
            {
                inside = false;
            }

        }

    }

    public static bool triangle_contains_point_2d_1(double[] t, double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_CONTAINS_POINT_2D_1 finds if a point is inside a triangle in 2D.
        //
        //  Discussion:
        //
        //    It is conventional to list the triangle vertices in counter clockwise
        //    order.  However, this routine does not require a particular order
        //    for the vertices.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double T[2*3], the triangle vertices.
        //
        //    Input, double P[2], the point to be checked.
        //
        //    Output, bool TRIANGLE_CONTAINS_POINT_2D_1, is TRUE if the points
        //    is inside the triangle or on its boundary, and FALSE otherwise.
        //
    {
        double[] c;
        int i;
        bool value;

        c = triangle_barycentric_2d(t, p);

        value = true;

        for (i = 0; i < 3; i++)
        {
            value = c[i] switch
            {
                < 0.0 => false,
                _ => value
            };
        }

        return value;
    }

    public static bool triangle_contains_point_2d_2(double[] t, double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_CONTAINS_POINT_2D_2 finds if a point is inside a triangle in 2D.
        //
        //  Discussion:
        //
        //    The routine assumes that the vertices are given in counter clockwise
        //    order.  If the triangle vertices are actually given in clockwise
        //    order, this routine will behave as though the triangle contains
        //    no points whatsoever!
        //
        //    The routine determines if P is "to the right of" each of the lines
        //    that bound the triangle.  It does this by computing the cross product
        //    of vectors from a vertex to its next vertex, and to P.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 June 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double T[2*3], the triangle vertices.
        //    The vertices should be given in counter clockwise order.
        //
        //    Input, double P[2], the point to be checked.
        //
        //    Output, bool TRIANGLE_CONTAINS_POINT_2D_2, is TRUE if P is inside
        //    the triangle or on its boundary.
        //
    {
        int j;
        int k;

        for (j = 0; j < 3; j++)
        {
            k = (j + 1) % 3;
            switch ((p[0] - t[0 + j * 2]) * (t[1 + k * 2] - t[1 + j * 2])
                    - (p[1] - t[1 + j * 2]) * (t[0 + k * 2] - t[0 + j * 2]))
            {
                case > 0.0:
                    return false;
            }
        }

        return true;

    }

    public static bool triangle_contains_point_2d_3(double[] t, double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_CONTAINS_POINT_2D_3 finds if a point is inside a triangle in 2D.
        //
        //  Discussion:
        //
        //    This routine is the same as TRIANGLE_CONTAINS_POINT_2D_2, except
        //    that it does not assume an ordering of the points.  It should
        //    work correctly whether the vertices of the triangle are listed
        //    in clockwise or counter clockwise order.
        //
        //    The routine determines if a point P is "to the right of" each of the lines
        //    that bound the triangle.  It does this by computing the cross product
        //    of vectors from a vertex to its next vertex, and to P.
        //
        //    The point is inside the triangle if it is to the right of all
        //    the lines, or to the left of all the lines.
        //
        //    This version was suggested by Paulo Ernesto of Maptek Brasil.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 June 2006
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double T[2*3], the triangle vertices.
        //
        //    Input, double P[2], the point to be checked.
        //
        //    Output, bool TRIANGLE_CONTAINS_POINT_2D_3, is TRUE if P is inside
        //    the triangle or on its boundary.
        //
    {
        double dir_new;
        double dir_old;
        int j;
        int k;

        dir_old = 0.0;

        for (j = 0; j < 3; j++)
        {
            k = (j + 1) % 3;

            dir_new = (p[0] - t[0 + j * 2]) * (t[1 + k * 2] - t[1 + j * 2])
                      - (p[1] - t[1 + j * 2]) * (t[0 + k * 2] - t[0 + j * 2]);

            switch (dir_new * dir_old)
            {
                case < 0.0:
                    return false;
            }

            if (dir_new != 0.0)
            {
                dir_old = dir_new;
            }
        }

        return true;

    }

    public static double triangle_diameter_2d(double[] t)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_DIAMETER_2D computes the diameter of a triangle in 2D.
        //
        //  Discussion:
        //
        //    The diameter of a triangle is the diameter of the smallest circle
        //    that can be drawn around the triangle.  At least two of the vertices
        //    of the triangle will intersect the circle, but not necessarily
        //    all three!
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
        //  Parameters:
        //
        //    Input, double T[2*3], the triangle vertices.
        //
        //    Output, double TRIANGLE_DIAMETER_2D, the diameter of the triangle.
        //
    {
        double a;
        double b;
        double c;
        double diam;
        //
        //  Compute the (squares of) the lengths of the sides.
        //
        a = Math.Sqrt(Math.Pow(t[0 + 1 * 2] - t[0 + 0 * 2], 2) + Math.Pow(t[1 + 1 * 2] - t[1 + 0 * 2], 2));
        b = Math.Sqrt(Math.Pow(t[0 + 2 * 2] - t[0 + 1 * 2], 2) + Math.Pow(t[1 + 2 * 2] - t[1 + 1 * 2], 2));
        c = Math.Sqrt(Math.Pow(t[0 + 0 * 2] - t[0 + 2 * 2], 2) + Math.Pow(t[1 + 0 * 2] - t[1 + 2 * 2], 2));
        switch (a)
        {
            //
            //  Take care of a zero side.
            //
            case 0.0:
                return Math.Sqrt(b);
        }

        switch (b)
        {
            case 0.0:
                return Math.Sqrt(c);
        }
        switch (c)
        {
            case 0.0:
                return Math.Sqrt(a);
        }

        //
        //  Make A the largest.
        //
        if (a < b)
        {
            typeMethods.r8_swap(ref a, ref b);
        }

        if (a < c)
        {
            typeMethods.r8_swap(ref a, ref c);
        }

        //
        //  If A is very large...
        //
        if (b + c < a)
        {
            diam = Math.Sqrt(a);
        }
        else
        {
            a = Math.Sqrt(a);
            b = Math.Sqrt(b);
            c = Math.Sqrt(c);
            diam = 2.0 * a * b * c / Math.Sqrt((a + b + c) * (-a + b + c)
                                                           * (a - b + c) * (a + b - c));
        }

        return diam;

    }

    public static double[] triangle_edge_length_2d(double[] t)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_EDGE_LENGTH_2D returns edge lengths of a triangle in 2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    17 August 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double T[2*3], the triangle vertices.
        //
        //    Output, double TRIANGLE_EDGE_LENGTH[3], the length of the edges.
        //
    {
        double[] edge_length;
        int j1;
        int j2;

        edge_length = new double[3];

        for (j1 = 0; j1 < 3; j1++)
        {
            j2 = typeMethods.i4_wrap(j1 + 1, 0, 2);
            edge_length[j1] = Math.Sqrt(Math.Pow(t[0 + j2 * 2] - t[0 + j1 * 2], 2)
                                        + Math.Pow(t[1 + j2 * 2] - t[1 + j1 * 2], 2));
        }

        return edge_length;
    }

    public static void triangle_gridpoints_2d(double[] t, int sub_num, int grid_max,
            ref int grid_num, ref double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_GRIDPOINTS_2D computes gridpoints within a triangle in 2D.
        //
        //  Discussion:
        //
        //    The gridpoints are computed by repeated halving of the triangle.
        //    The 0-th set of grid points is the vertices themselves.
        //    The first set of grid points is the midpoints of the sides.
        //    These points can be used to draw 4 triangles that make up the original
        //    triangle.  The second set of grid points is the side midpoints and centers
        //    of these four triangles.
        //
        //    SUB_NUM                     GRID_NUM
        //    -----                        -----
        //        0      1                  =  1  (centroid)
        //        1      1 + 2              =  3  (vertices)
        //        2      1 + 2 + 3          =  6
        //        3      1 + 2 + 3 + 4      = 10
        //        4      1 + 2 + 3 + 4 + 5  = 15
        //
        //    GRID_NUM is the sum of the integers from 1 to SUB_NUM+1 or
        //
        //      GRID_NUM = (SUB_NUM+1) * (SUB_NUM+2) / 2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double T[2*3], the triangle vertices.
        //
        //    Input, int SUB_NUM, the number of subdivisions.
        //
        //    Input, int GRID_MAX, the maximum number of grid points.
        //
        //    Output, int *GRID_NUM, the number of grid points returned.
        //
        //    Output, double P[2*(*GRID_NUM)], coordinates of the grid points.
        //
    {
        int i;
        int j;

        grid_num = 0;
        switch (sub_num)
        {
            //
            //  Special case, SUB_NUM = 0.
            //
            case 0:
            {
                if (grid_num + 1 <= grid_max)
                {
                    p[0 + grid_num * 2] = (t[0 + 0 * 2] + t[0 + 1 * 2] + t[0 + 2 * 2]) / 3.0;
                    p[1 + grid_num * 2] = (t[1 + 0 * 2] + t[1 + 1 * 2] + t[1 + 2 * 2]) / 3.0;
                    grid_num += 1;
                }

                return;
            }
        }

        for (i = 0; i <= sub_num; i++)
        {
            for (j = 0; j <= sub_num - i; j++)
            {
                if (grid_max <= grid_num)
                {
                    return;
                }

                p[0 + grid_num * 2] = (i * t[0 + 0 * 2]
                                       + j * t[0 + 1 * 2]
                                       + (sub_num - i - j) * t[0 + 2 * 2])
                                      / sub_num;

                p[1 + grid_num * 2] = (i * t[1 + 0 * 2]
                                       + j * t[1 + 1 * 2]
                                       + (sub_num - i - j) * t[1 + 2 * 2])
                                      / sub_num;

                grid_num += 1;
            }
        }

    }

    public static void triangle_incenter_2d(double[] t, ref double[] pc)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_INCENTER_2D computes the incenter of a triangle in 2D.
        //
        //  Discussion:
        //
        //    The incenter of a triangle is the center of the inscribed circle.
        //
        //    The inscribed circle of a triangle is the largest circle that can
        //    be drawn inside the triangle.
        //
        //    The inscribed circle is tangent to all three sides of the triangle.
        //
        //    The angle bisectors of the triangle intersect at the center of the
        //    inscribed circle.
        //
        //    In geometry, the incenter is often represented by "I".
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 July 2005
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
        //    Input, double T[2*3], the triangle vertices.
        //
        //    Output, double PC[2], the coordinates of the center of the
        //    inscribed circle.
        //
    {
        double perim;
        double s12;
        double s23;
        double s31;

        s12 = Math.Sqrt(Math.Pow(t[0 + 1 * 2] - t[0 + 0 * 2], 2)
                        + Math.Pow(t[1 + 1 * 2] - t[1 + 0 * 2], 2));
        s23 = Math.Sqrt(Math.Pow(t[0 + 2 * 2] - t[0 + 1 * 2], 2)
                        + Math.Pow(t[1 + 2 * 2] - t[1 + 1 * 2], 2));
        s31 = Math.Sqrt(Math.Pow(t[0 + 0 * 2] - t[0 + 2 * 2], 2)
                        + Math.Pow(t[1 + 0 * 2] - t[1 + 2 * 2], 2));

        perim = s12 + s23 + s31;

        switch (perim)
        {
            case 0.0:
                pc[0] = t[0 + 0 * 2];
                pc[1] = t[1 + 0 * 2];
                break;
            default:
                pc[0] = (s23 * t[0 + 0 * 2] + s31 * t[0 + 1 * 2] + s12 * t[0 + 2 * 2]) / perim;
                pc[1] = (s23 * t[1 + 0 * 2] + s31 * t[1 + 1 * 2] + s12 * t[1 + 2 * 2]) / perim;
                break;
        }

    }

    public static void triangle_incircle_2d(double[] t, ref double[] pc, ref double r)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_INCIRCLE_2D computes the inscribed circle of a triangle in 2D.
        //
        //  Discussion:
        //
        //    The inscribed circle of a triangle is the largest circle that can
        //    be drawn inside the triangle.  It is tangent to all three sides,
        //    and the lines from its center to the vertices bisect the angles
        //    made by each vertex.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 July 2005
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
        //    Input, double T[2*3], the triangle vertices.
        //
        //    Output, double PC[2], *R, the center of the inscribed circle, and its radius.
        //
    {
        double perim;
        double s12;
        double s23;
        double s31;

        s12 = Math.Sqrt(Math.Pow(t[0 + 1 * 2] - t[0 + 0 * 2], 2)
                        + Math.Pow(t[1 + 1 * 2] - t[1 + 0 * 2], 2));
        s23 = Math.Sqrt(Math.Pow(t[0 + 2 * 2] - t[0 + 1 * 2], 2)
                        + Math.Pow(t[1 + 2 * 2] - t[1 + 1 * 2], 2));
        s31 = Math.Sqrt(Math.Pow(t[0 + 0 * 2] - t[0 + 2 * 2], 2)
                        + Math.Pow(t[1 + 0 * 2] - t[1 + 2 * 2], 2));

        perim = s12 + s23 + s31;

        switch (perim)
        {
            case 0.0:
                r = 0.0;
                pc[0] = t[0 + 0 * 2];
                pc[1] = t[1 + 0 * 2];
                break;
            default:
                pc[0] = (s23 * t[0 + 0 * 2] + s31 * t[0 + 1 * 2] + s12 * t[0 + 2 * 2]) / perim;
                pc[1] = (s23 * t[1 + 0 * 2] + s31 * t[1 + 1 * 2] + s12 * t[1 + 2 * 2]) / perim;

                r = 0.5 * Math.Sqrt(
                    (-s12 + s23 + s31)
                    * (+s12 - s23 + s31)
                    * (+s12 + s23 - s31) / perim);
                break;
        }

    }

    public static double triangle_inradius_2d(double[] t)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_INRADIUS_2D computes the inradius of a triangle in 2D.
        //
        //  Discussion:
        //
        //    The inscribed circle of a triangle is the largest circle that can
        //    be drawn inside the triangle.  It is tangent to all three sides,
        //    and the lines from its center to the vertices bisect the angles
        //    made by each vertex.
        //
        //    The inradius is the radius of the inscribed circle.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 July 2005
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
        //    Input, double T[2*3], the triangle vertices.
        //
        //    Output, double TRIANGLE_INRADIUS_2D, the inradius.
        //
    {
        double perim;
        double r;
        double s12;
        double s23;
        double s31;

        s12 = Math.Sqrt(Math.Pow(t[0 + 1 * 2] - t[0 + 0 * 2], 2)
                        + Math.Pow(t[1 + 1 * 2] - t[1 + 0 * 2], 2));
        s23 = Math.Sqrt(Math.Pow(t[0 + 2 * 2] - t[0 + 1 * 2], 2)
                        + Math.Pow(t[1 + 2 * 2] - t[1 + 1 * 2], 2));
        s31 = Math.Sqrt(Math.Pow(t[0 + 0 * 2] - t[0 + 2 * 2], 2)
                        + Math.Pow(t[1 + 0 * 2] - t[1 + 2 * 2], 2));

        perim = s12 + s23 + s31;

        r = perim switch
        {
            0.0 => 0.0,
            _ => 0.5 * Math.Sqrt((-s12 + s23 + s31) * (+s12 - s23 + s31) * (+s12 + s23 - s31) / perim)
        };

        return r;

    }

    public static bool triangle_is_degenerate_nd(int dim_num, double[] t)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_IS_DEGENERATE_ND finds if a triangle is degenerate in ND.
        //
        //  Discussion:
        //
        //    A triangle in ND is described by the coordinates of its 3 vertices.
        //
        //    A triangle in ND is degenerate if any two vertices are equal.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, double T[DIM_NUM*3], the triangle vertices.
        //
        //    Output, bool TRIANGLE_IS_DEGENERATE_ND, is TRUE if the
        //    triangle is degenerate.
        //
    {
        bool value;

        value =
            typeMethods.r8vec_eq(dim_num, t, t, +0 * dim_num, +1 * dim_num) ||
            typeMethods.r8vec_eq(dim_num, t, t, +1 * dim_num, +2 * dim_num) ||
            typeMethods.r8vec_eq(dim_num, t, t, +2 * dim_num, +0 * dim_num);

        return value;
    }

    public static void triangle_lattice_layer_point_next(int[] c, ref int[] v, ref bool more)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_LATTICE_LAYER_POINT_NEXT: next triangle lattice layer point.
        //
        //  Discussion:
        //
        //    The triangle lattice layer L is bounded by the lines
        //
        //      0 <= X,
        //      0 <= Y,
        //      L - 1 < X / C[0] + Y / C[1] <= L.
        //
        //    In particular, layer L = 0 always contains the single point (0,0).
        //
        //    This function returns, one at a time, the points that lie within
        //    a given triangle lattice layer.
        //
        //    Thus, if we set C[0] = 2, C[1] = 3, then we get the following layers:
        //
        //    L = 0: (0,0)
        //    L = 1: (1,0), (2,0), (0,1), (1,1), (0,2), (0,3)
        //    L = 2: (3,0), (4,0), (2,1), (3,1), (1,2), (2,2), (1,3), (2,3),
        //           (0,4), (1,4), (0,5), (0,6).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 July 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int C[3], coefficients defining the
        //    lattice layer.  Entry C[2] contains the layer index.
        //    C[0] and C[1] should be positive, and C[2] must be nonnegative.
        //
        //    Input/output, int V[2].  On first call for a given layer,
        //    the input value of V is not important.  On a repeated call for the same
        //    layer, the input value of V should be the output value from the previous
        //    call.  On output, V contains the next lattice layer point.
        //
        //    Input/output, bool *MORE.  On input, set MORE to FALSE to indicate
        //    that this is the first call for a given layer.  Thereafter, the input
        //    value should be the output value from the previous call.  On output,
        //    MORE is TRUE if the returned value V is a new point.
        //    If the output value is FALSE, then no more points were found,
        //    and V was reset to 0, and the lattice layer has been exhausted.
        //
    {
        int c1n;
        int n = 2;
        int rhs1;
        int rhs2;
        switch (c[n])
        {
            //
            //  Treat layer C[N] = 0 specially.
            //
            case 0:
            {
                switch (more)
                {
                    case false:
                        v[0] = 0;
                        v[1] = 0;
                        more = true;
                        break;
                    default:
                        more = false;
                        break;
                }

                return;
            }
        }

        switch (more)
        {
            //
            //  Compute first point.
            //
            case false:
                v[0] = (c[n] - 1) * c[0] + 1;
                v[1] = 0;
                more = true;
                break;
            default:
            {
                c1n = typeMethods.i4vec_lcm(n, c);
                rhs1 = c1n * (c[n] - 1);
                rhs2 = c1n * c[n];

                if (c[1] * (v[0] + 1) + c[0] * v[1] <= rhs2)
                {
                    v[0] += 1;
                }
                else
                {
                    v[0] = (rhs1 - c[0] * (v[1] + 1)) / c[1];
                    v[0] = Math.Max(v[0], 0);
                    v[1] += 1;
                    if (c[1] * v[0] + c[0] * v[1] <= rhs1)
                    {
                        v[0] += 1;
                    }

                    if (c[1] * v[0] + c[0] * v[1] <= rhs2)
                    {
                    }
                    else
                    {
                        v[0] = 0;
                        v[1] = 0;
                        more = false;
                    }
                }

                break;
            }
        }
    }

    public static void triangle_lattice_point_next(int[] c, ref int[] v, ref bool more)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_LATTICE_POINT_NEXT returns the next triangle lattice point.
        //
        //  Discussion:
        //
        //    The lattice triangle is defined by the vertices:
        //
        //      (0,0), (C[2]/C[0], 0) and (0,C[2]/C[1])
        //
        //    The lattice triangle is bounded by the lines
        //
        //      0 <= X,
        //      0 <= Y
        //      X / C[0] + Y / C[1] <= C[2]
        //
        //    Lattice points are listed one at a time, starting at the origin,
        //    with X increasing first.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 July 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int C[3], coefficients defining the
        //    lattice triangle.  These should be positive.
        //
        //    Input/output, int V[2].  On first call, the input
        //    value is not important.  On a repeated call, the input value should
        //    be the output value from the previous call.  On output, V contains
        //    the next lattice point.
        //
        //    Input/output, bool *MORE.  On input, set MORE to FALSE to indicate
        //    that this is the first call for a given triangle.  Thereafter, the input
        //    value should be the output value from the previous call.  On output,
        //    MORE is TRUE if the returned value V is a new lattice point.
        //    If the output value is FALSE, then no more lattice points were found,
        //    and V was reset to 0, and the routine should not be called further
        //    for this triangle.
        //
    {
        int c1n;
        int n = 2;
        int rhs;

        switch (more)
        {
            case false:
                v[0] = 0;
                v[1] = 0;
                more = true;
                break;
            default:
            {
                c1n = typeMethods.i4vec_lcm(n, c);
                rhs = c1n * c[n];

                if (c[1] * (v[0] + 1) + c[0] * v[1] <= rhs)
                {
                    v[0] += 1;
                }
                else
                {
                    v[0] = 0;
                    if (c[1] * v[0] + c[0] * (v[1] + 1) <= rhs)
                    {
                        v[1] += 1;
                    }
                    else
                    {
                        v[1] = 0;
                        more = false;
                    }
                }

                break;
            }
        }
    }

    public static void triangle_line_imp_int_2d(double[] t, double a, double b, double c,
            ref int int_num, ref double[] pint)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_LINE_IMP_INT_2D finds where an implicit line intersects a triangle in 2D.
        //
        //  Discussion:
        //
        //    An implicit line is the set of points P satisfying
        //
        //      A * P[0] + B * P[1] + C = 0
        //
        //    where at least one of A and B is not zero.
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
        //  Parameters:
        //
        //    Input, double T[2*3], the triangle vertices.
        //
        //    Input, double A, B, C, determine the equation of the line:
        //    A*X + B*Y + C = 0.
        //
        //    Output, int *INT_NUM, the number of points of intersection
        //    of the line with the triangle.  INT_NUM may be 0, 1, 2 or 3.
        //
        //    Output, double PINT[3*3], contains the intersection points.
        //
    {
        int DIM_NUM = 2;

        double a1 = 0;
        double b1 = 0;
        double c1 = 0;
        int ival = 0;
        int n;
        double[] p = new double[DIM_NUM];
        int r;
        int s;
        double test1;
        double test2;

        n = 0;

        for (r = 0; r < 3; r++)
        {
            s = typeMethods.i4_wrap(r + 1, 0, 2);
            //
            //  Get the implicit form of the line through vertices R and R+1.
            //
            LineNS.Geometry.line_exp2imp_2d(t, t, ref a1, ref b1, ref c1, +0 + r * 2, +0 + s * 2);
            //
            //  Seek an intersection with the original line.
            //
            LineNS.Geometry.lines_imp_int_2d(a, b, c, a1, b1, c1, ref ival, ref p);
            switch (ival)
            {
                //
                //  If there is an intersection, then determine if it happens between
                //  the two vertices.
                //
                case 1:
                {
                    test1 = (p[0] - t[0 + r * 2]) * (t[0 + s * 2] - t[0 + r * 2])
                            + (p[1] - t[1 + r * 2]) * (t[1 + s * 2] - t[1 + r * 2]);

                    test2 = (t[0 + s * 2] - t[0 + r * 2]) * (t[0 + s * 2] - t[0 + r * 2])
                            + (t[1 + s * 2] - t[1 + r * 2]) * (t[1 + s * 2] - t[1 + r * 2]);

                    switch (test1)
                    {
                        case >= 0 when test1 <= test2:
                            pint[0 + n * 2] = p[0];
                            pint[1 + n * 2] = p[1];
                            n += 1;
                            break;
                    }

                    break;
                }
            }
        }

        int_num = n;

    }

    public static int triangle_orientation_2d(double[] t)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_ORIENTATION_2D determines the orientation of a triangle in 2D.
        //
        //  Discussion:
        //
        //    Three distinct non-colinear points in the plane define a circle.
        //    If the points are visited in the order (x1,y1), (x2,y2), and then
        //    (x3,y3), this motion defines a clockwise or counter clockwise
        //    rotation along the circle.
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
        //  Parameters:
        //
        //    Input, double T[2*3], the triangle vertices.
        //
        //    Output, int TRIANGLE_ORIENTATION_2D, reports if the three points lie
        //    clockwise on the circle that passes through them.  The possible
        //    return values are:
        //    0, the points are distinct, noncolinear, and lie counter clockwise
        //    on their circle.
        //    1, the points are distinct, noncolinear, and lie clockwise
        //    on their circle.
        //    2, the points are distinct and colinear.
        //    3, at least two of the points are identical.
        //
    {
        double det;
        int value = 0;

        if (typeMethods.r8vec_eq(2, t, t, +0 * 2, +1 * 2) ||
            typeMethods.r8vec_eq(2, t, t, +1 * 2, +2 * 2) ||
            typeMethods.r8vec_eq(2, t, t, +2 * 2, +0 * 2))
        {
            value = 3;
            return value;
        }

        det = (t[0 + 0 * 2] - t[0 + 2 * 2]) * (t[1 + 1 * 2] - t[1 + 2 * 2])
              - (t[0 + 1 * 2] - t[0 + 2 * 2]) * (t[1 + 0 * 2] - t[1 + 2 * 2]);

        value = det switch
        {
            0.0 => 2,
            < 0.0 => 1,
            > 0.0 => 0,
            _ => value
        };

        return value;

    }

    public static void triangle_orthocenter_2d(double[] t, ref double[] p, ref bool flag)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_ORTHOCENTER_2D computes the orthocenter of a triangle in 2D.
        //
        //  Discussion:
        //
        //    The orthocenter is defined as the intersection of the three altitudes
        //    of a triangle.
        //
        //    An altitude of a triangle is the line through a vertex of the triangle
        //    and perpendicular to the opposite side.
        //
        //    In geometry, the orthocenter of a triangle is often symbolized by "H".
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 July 2009
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
        //    Input, double T[2*3], the triangle vertices.
        //
        //    Output, double P[2], the coordinates of the orthocenter of the triangle.
        //
        //    Output, bool *FLAG, is TRUE if the value could not be computed.
        //
    {
        int ival = 0;
        double[] p23;
        double[] p31;
        //
        //  Determine a point P23 common to the line through P2 and P3 and
        //  its perpendicular through P1.
        //
        p23 = LineNS.Geometry.line_exp_perp_2d(t, t, t, ref flag, +1 * 2, +2 * 2 + 0 * 2);

        switch (flag)
        {
            case true:
                p[0] = typeMethods.r8_huge();
                p[1] = typeMethods.r8_huge();
                break;
        }

        //
        //  Determine a point P31 common to the line through P3 and P1 and
        //  its perpendicular through P2.
        //
        p31 = LineNS.Geometry.line_exp_perp_2d(t, t, t, ref flag, +2 * 2, +0 * 2, +1 * 2);

        switch (flag)
        {
            case true:
                p[0] = typeMethods.r8_huge();
                p[1] = typeMethods.r8_huge();
                break;
        }

        //
        //  Determine P, the intersection of the lines through P1 and P23, and
        //  through P2 and P31.
        //
        LineNS.Geometry.lines_exp_int_2d(t, p23, t, p31, ref ival, ref p, p1Index: +0 * 2,
            p3Index: +1 * 2);

        if (ival != 1)
        {
            p[0] = typeMethods.r8_huge();
            p[1] = typeMethods.r8_huge();
            flag = true;
        }

    }

    public static double triangle_point_dist_2d(double[] t, double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_POINT_DIST_2D: distance ( triangle, point ) in 2D.
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
        //  Parameters:
        //
        //    Input, double T[2*3], the triangle vertices.
        //
        //    Input, double P[2], the point which is to be checked.
        //
        //    Output, double TRIANGLE_POINT_DIST_2D, the distance from the point to the triangle.
        //    DIST is zero if the point lies exactly on the triangle.
        //
    {
        int DIM_NUM = 2;

        double value = 0;

        value =
            Segments.segment_point_dist_2d(t, t, p,  + 0 * DIM_NUM,  + 1 * DIM_NUM);
        value = Math.Min(value,
            Segments.segment_point_dist_2d(t, t, p,  + 1 * DIM_NUM,  + 2 * DIM_NUM));
        value = Math.Min(value,
            Segments.segment_point_dist_2d(t, t, p,  + 2 * DIM_NUM,  + 0 * DIM_NUM));

        return value;

    }

    public static double triangle_point_dist_3d(double[] t, double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_POINT_DIST_3D: distance ( triangle, point ) in 3D.
        //
        //  Discussion:
        //
        //    Thanks to Gordon Griesel for pointing out that a triangle in 3D
        //    has to have coordinates in 3D as well.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 October 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double T[3*3], the triangle vertices.
        //
        //    Input, double P[3], the point which is to be checked.
        //
        //    Output, double TRIANGLE_POINT_DIST_3D, the distance from the point
        //    to the triangle.
        //
    {
        int DIM_NUM = 3;

        double value = 0;

        value =
            Segments.segment_point_dist_3d(t, t, p,  + 0 * DIM_NUM,  + 1 * DIM_NUM);
        value = Math.Min(value,
            Segments.segment_point_dist_3d(t, t, p,  + 1 * DIM_NUM,  + 2 * DIM_NUM));
        value = Math.Min(value,
            Segments.segment_point_dist_3d(t, t, p,  + 2 * DIM_NUM,  + 0 * DIM_NUM));

        return value;

    }

    public static double triangle_point_dist_signed_2d(double[] t, double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_POINT_DIST_SIGNED_2D: signed distance ( triangle, point ) in 2D.
        //
        //  Discussion:
        //
        //    If the signed distance is:
        //    0, the point is on the boundary of the triangle;
        //    negative, the point is in the triangle;
        //    positive, the point is outside the triangle.
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
        //  Parameters:
        //
        //    Input, double T[2*3], the triangle vertices.
        //    These should be given in counter clockwise order.
        //
        //    Input, double P[2], the point which is to be checked.
        //
        //    Output, double TRIANGLE_POINT_DIST_SIGNED_2D, the signed distance from the
        //    point to the triangle.
        //
    {
        double dis12;
        double dis23;
        double dis31;
        double value = 0;
        //
        //  Compute the signed line-distances to the point.
        //
        dis12 = LineNS.Geometry.line_exp_point_dist_signed_2d(t, t, p, +0 * 2, +1 * 2);
        dis23 = LineNS.Geometry.line_exp_point_dist_signed_2d(t, t, p, +1 * 2, +2 * 2);
        dis31 = LineNS.Geometry.line_exp_point_dist_signed_2d(t, t, p, +2 * 2, +0 * 2);
        switch (dis12)
        {
            //
            //  If the point is inside the triangle, all the line-distances are negative.
            //  The largest (negative) line-distance has the smallest magnitude,
            //  and is the signed triangle-distance.
            //
            case <= 0.0 when dis23 <= 0.0 && dis31 <= 0.0:
                value = dis12;
                value = Math.Max(value, dis23);
                value = Math.Max(value, dis31);
                break;
            //
            default:
                value = Segments.segment_point_dist_2d(t, t, p,  + 0 * 2,  + 1 * 2);
                value = Math.Min(value, Segments.segment_point_dist_2d(t, t, p,  + 1 * 2,  + 2 * 2));
                value = Math.Min(value, Segments.segment_point_dist_2d(t, t, p,  + 2 * 2,  + 0 * 2));
                break;
        }

        return value;

    }

    public static void triangle_point_near_2d(double[] t, double[] p, ref double[] pn,
            ref double dist)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_POINT_NEAR_2D computes the nearest triangle point to a point in 2D.
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
        //  Parameters:
        //
        //    Input, double T[2*3], the triangle vertices.
        //
        //    Input, double P[2], the point whose nearest neighbor
        //    on the line is to be determined.
        //
        //    Output, double PN[2], the nearest point to P.
        //
        //    Output, double *DIST, the distance from the point to the triangle.
        //
    {
        int DIM_NUM = 2;

        double dist12 = 0;
        double dist23 = 0;
        double dist31 = 0;
        double tval = 0;
        double[] pn12 = new double[DIM_NUM];
        double[] pn23 = new double[DIM_NUM];
        double[] pn31 = new double[DIM_NUM];
        //
        //  Find the distance to each of the line segments that make up the edges
        //  of the triangle.
        //
        Segments.segment_point_near_2d(t, t, p, ref pn12, ref dist12, ref tval,  + 0 * DIM_NUM, + 1 * DIM_NUM);

        Segments.segment_point_near_2d(t, t, p, ref pn23, ref dist23, ref tval,  + 1 * DIM_NUM, + 2 * DIM_NUM);

        Segments.segment_point_near_2d(t, t, p, ref pn31, ref dist31, ref tval,  + 2 * DIM_NUM, + 0 * DIM_NUM);

        if (dist12 <= dist23 && dist12 <= dist31)
        {
            dist = dist12;
            typeMethods.r8vec_copy(DIM_NUM, pn12, ref pn);
        }
        else if (dist23 <= dist12 && dist23 <= dist31)
        {
            dist = dist23;
            typeMethods.r8vec_copy(DIM_NUM, pn23, ref pn);
        }
        else
        {
            dist = dist31;
            typeMethods.r8vec_copy(DIM_NUM, pn31, ref pn);
        }

    }

    public static double triangle_quality_2d(double[] t)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_QUALITY_2D: "quality" of a triangle in 2D.
        //
        //  Discussion:
        //
        //    The quality of a triangle is 2 times the ratio of the radius of the
        //    inscribed circle divided by that of the circumscribed circle.  An
        //    equilateral triangle achieves the maximum possible quality of 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    31 July 2009
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
        //    Input, double T[2*3], the triangle vertices.
        //
        //    Output, double TRIANGLE_QUALITY_2D, the quality of the triangle.
        //
    {
        int DIM_NUM = 2;

        double a;
        double b;
        double c;
        int i;
        double value = 0;
        //
        //  Compute the length of each side.
        //
        a = 0.0;
        b = 0.0;
        c = 0.0;

        for (i = 0; i < DIM_NUM; i++)
        {
            a += Math.Pow(t[i + 0 * DIM_NUM] - t[i + 1 * DIM_NUM], 2);
            b += Math.Pow(t[i + 1 * DIM_NUM] - t[i + 2 * DIM_NUM], 2);
            c += Math.Pow(t[i + 2 * DIM_NUM] - t[i + 0 * DIM_NUM], 2);
        }

        a = Math.Sqrt(a);
        b = Math.Sqrt(b);
        c = Math.Sqrt(c);

        value = (a * b * c) switch
        {
            0.0 => 0.0,
            _ => (-a + b + c) * (a - b + c) * (a + b - c) / (a * b * c)
        };

        return value;

    }

    public static int triangle_right_lattice_point_num_2d(int a, int b)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_RIGHT_LATTICE_POINT_NUM_2D: count lattice points.
        //
        //  Discussion:
        //
        //    The triangle is assumed to be a right triangle which, without loss
        //    of generality, has the coordinates:
        //
        //    ( (0,0), (a,0), (0,b) )
        //
        //    The routine returns the number of integer lattice points that appear
        //    inside the triangle or on its edges or vertices.
        //
        //    The formula for this function occurred to me (JVB) after some thought,
        //    on 06 July 2009.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 July 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int A, B, define the vertices.
        //
        //    Output, int N, the number of lattice points.
        //
    {
        int n;

        n = ((a + 1) * (b + 1) + typeMethods.i4_gcd(a, b) + 1) / 2;

        return n;
    }

    public static void triangle_sample(double[] t, int n, ref int seed, ref double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_SAMPLE returns random points in a triangle.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double T[2*3], the triangle vertices.
        //
        //    Input, int N, the number of points to sample.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double P[2*N], a random point in the triangle.
        //
    {
        int DIM_NUM = 2;

        double alpha;
        double beta;
        int j;
        double r;
        double[] p12 = new double[DIM_NUM];
        double[] p13 = new double[DIM_NUM];

        for (j = 0; j < n; j++)
        {
            r = UniformRNG.r8_uniform_01(ref seed);
            //
            //  Interpret R as a percentage of the triangle's area.
            //
            //  Imagine a line L, parallel to side 1, so that the area between
            //  vertex 1 and line L is R percent of the full triangle's area.
            //
            //  The line L will intersect sides 2 and 3 at a fraction
            //  ALPHA = Math.Sqrt ( R ) of the distance from vertex 1 to vertices 2 and 3.
            //
            alpha = Math.Sqrt(r);
            //
            //  Determine the coordinates of the points on sides 2 and 3 intersected
            //  by line L.
            //
            p12[0] = (1.0 - alpha) * t[0 + 0 * 2] + alpha * t[0 + 1 * 2];
            p12[1] = (1.0 - alpha) * t[1 + 0 * 2] + alpha * t[1 + 1 * 2];

            p13[0] = (1.0 - alpha) * t[0 + 0 * 2] + alpha * t[0 + 2 * 2];
            ;
            p13[1] = (1.0 - alpha) * t[1 + 0 * 2] + alpha * t[1 + 2 * 2];
            ;
            //
            //  Now choose, uniformly at random, a point on the line L.
            //
            beta = UniformRNG.r8_uniform_01(ref seed);

            p[0 + j * 2] = (1.0 - beta) * p12[0] + beta * p13[0];
            p[1 + j * 2] = (1.0 - beta) * p12[1] + beta * p13[1];
        }

    }

    public static int triangle_unit_lattice_point_num_2d(int s)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_UNIT_LATTICE_POINT_NUM_2D: count lattice points.
        //
        //  Discussion:
        //
        //    The triangle is assumed to be the unit triangle:
        //
        //    ( (0,0), (1,0), (0,1) )
        //
        //    or a copy of this triangle scaled by an integer S:
        //
        //    ( (0,0), (S,0), (0,S) ).
        //
        //    The routine returns the number of integer lattice points that appear
        //    inside the triangle or on its edges or vertices.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 July 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Matthias Beck, Sinai Robins,
        //    Computing the Continuous Discretely,
        //    Springer, 2006,
        //    ISBN13: 978-0387291390,
        //    LC: QA640.7.B43.
        //
        //  Parameters:
        //
        //    Input, int S, the scale factor.
        //
        //    Output, int TRIANGLE_UNIT_LATTICE_POINT_NUM_2D, the number of lattice points.
        //
    {
        int n;

        n = (s + 2) * (s + 1) / 2;

        return n;
    }

    public static void triangle_xsi_to_xy_2d(double[] t, double[] xsi, ref double[] p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_XSI_TO_XY_2D converts from barycentric to XY coordinates in 2D.
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
        //  Parameters:
        //
        //    Input, double T[2*3], the triangle vertices.
        //
        //    Input, double XSI[3], the barycentric coordinates of a point.
        //
        //    Output, double P[2], the Cartesian coordinates of the point.
        //
    {
        p[0] = xsi[0] * t[0 + 0 * 2] + xsi[1] * t[0 + 1 * 2] + xsi[2] * t[0 + 2 * 2];
        p[1] = xsi[0] * t[1 + 0 * 2] + xsi[1] * t[1 + 1 * 2] + xsi[2] * t[1 + 2 * 2];

    }

    public static void triangle_xy_to_xsi_2d(double[] t, double[] p, double[] xsi)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_XY_TO_XSI_2D converts from XY to barycentric in 2D.
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
        //  Parameters:
        //
        //    Input, double T[2*3], the triangle vertices.
        //
        //    Input, double P[2], the XY coordinates of a point.
        //
        //    Output, double XSI[3], the barycentric coordinates of the point.
        //
    {
        double det;

        det = (t[0 + 0 * 2] - t[0 + 2 * 2]) * (t[1 + 1 * 2] - t[1 + 2 * 2])
              - (t[0 + 1 * 2] - t[0 + 2 * 2]) * (t[1 + 0 * 2] - t[1 + 2 * 2]);

        xsi[0] = ((t[1 + 1 * 2] - t[1 + 2 * 2]) * (p[0] - t[0 + 2 * 2])
                  - (t[0 + 1 * 2] - t[0 + 2 * 2]) * (p[1] - t[1 + 2 * 2])) / det;

        xsi[1] = (-(t[1 + 0 * 2] - t[1 + 2 * 2]) * (p[0] - t[0 + 2 * 2])
                  + (t[0 + 0 * 2] - t[0 + 2 * 2]) * (p[1] - t[1 + 2 * 2])) / det;

        xsi[2] = 1.0 - xsi[0] - xsi[1];

    }
}