using System;
using Burkardt.Types;

namespace Burkardt.Vector
{
    public static class Geometry
    {
        public static void vector_directions_nd(int dim_num, double[] v, ref double[] angle)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    VECTOR_DIRECTIONS_ND returns the direction angles of a vector in ND.
            //
            //  Discussion:
            //
            //    Let V be the vector, and let E(I) be the I-th unit coordinate axis vector.
            //    The I-th direction angle is the angle between V and E(I), which is
            //    the angle whose cosine is equal to the direction cosine:
            //
            //      Direction_Cosine(I) = V dot E(I) / |V|.
            //
            //    If V is the null or zero vector, then the direction cosines and
            //    direction angles are undefined, and this routine simply returns
            //    zeroes.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    01 August 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the spatial dimension.
            //
            //    Input, double V[DIM_NUM], the vector.
            //
            //    Output, double ANGLE[DIM_NUM], the direction angles, in radians,
            //    that the vector V makes with the coordinate axes.
            //
        {
            int i;
            double vnorm;
            //
            //  Get the norm of the vector.
            //
            vnorm = typeMethods.r8vec_norm(dim_num, v);

            if (vnorm == 0.0)
            {
                typeMethods.r8vec_zero(dim_num, ref angle);
                return;
            }

            for (i = 0; i < dim_num; i++)
            {
                angle[i] = Math.Acos(v[i] / vnorm);
            }

        }

        public static void vector_rotate_2d(double[] v1, double angle, ref double[] v2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    VECTOR_ROTATE_2D rotates a vector around the origin in 2D.
            //
            //  Discussion:
            //
            //    To see why this formula is so, consider that the original point
            //    has the form ( R Math.Cos Theta, R Math.Sin Theta ), and the rotated point
            //    has the form ( R Math.Cos ( Theta + Angle ), R Math.Sin ( Theta + Angle ) ).
            //    Now use the addition formulas for cosine and sine to relate
            //    the new point to the old one:
            //
            //      ( X2 ) = ( Math.Cos Angle  - Math.Sin Angle ) * ( X1 )
            //      ( Y2 )   ( Math.Sin Angle    Math.Cos Angle )   ( Y1 )
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
            //    Input, double V1[2], the vector to be rotated.
            //
            //    Input, double ANGLE, the angle, in radians, of the rotation to be
            //    carried out.  A positive angle rotates the vector in the
            //    counter clockwise direction.
            //
            //    Output, double V2[2], the rotated vector.
            //
        {
            v2[0] = Math.Cos(angle) * v1[0] - Math.Sin(angle) * v1[1];
            v2[1] = Math.Sin(angle) * v1[0] + Math.Cos(angle) * v1[1];

        }

        public static void vector_rotate_3d(double[] p1, double[] pa, double angle, ref double[] p2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    VECTOR_ROTATE_3D rotates a vector around an axis vector in 3D.
            //
            //  Discussion:
            //
            //    Thanks to Cody Farnell for correcting some errors in a previous
            //    version of this routine!
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    18 May 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double P1[3], the components of the vector to be rotated.
            //
            //    Input, double PA[3], the vector about which the rotation is to
            //    be carried out.
            //
            //    Input, double ANGLE, the angle, in radians, of the rotation to be
            //    carried out.
            //
            //    Output, double P2[3], the rotated vector.
            //
        {
            int DIM_NUM = 3;

            double axis_norm;
            double dot;
            double normn;
            double[] pn = new double[DIM_NUM];
            double[] pn2;
            double[] pp = new double[DIM_NUM];
            double[] pr = new double[DIM_NUM];
            //
            //  Compute the length of the rotation axis.
            //
            axis_norm = typeMethods.r8vec_norm(DIM_NUM, pa);

            if (axis_norm == 0.0)
            {
                typeMethods.r8vec_copy(DIM_NUM, p1, ref p2);
                return;
            }

            //
            //  Compute the dot product of the vector and the (unit) rotation axis.
            //
            dot = typeMethods.r8vec_dot_product(DIM_NUM, p1, pa) / axis_norm;
            //
            //  Compute the parallel component of the vector.
            //
            pp[0] = dot * pa[0] / axis_norm;
            pp[1] = dot * pa[1] / axis_norm;
            pp[2] = dot * pa[2] / axis_norm;
            //
            //  Compute the normal component of the vector.
            //
            pn[0] = p1[0] - pp[0];
            pn[1] = p1[1] - pp[1];
            pn[2] = p1[2] - pp[2];

            normn = typeMethods.r8vec_norm(DIM_NUM, pn);

            if (normn == 0.0)
            {
                typeMethods.r8vec_copy(DIM_NUM, pp, ref p2);
                return;
            }

            vector_unit_nd(3, ref pn);
            //
            //  Compute a second vector, lying in the plane, perpendicular
            //  to P1, and forming a right-handed system...
            //
            pn2 = typeMethods.r8vec_cross_product_3d(pa, pn);

            vector_unit_nd(3, ref pn2);
            //
            //  Rotate the normal component by the angle.
            //
            pr[0] = normn * (Math.Cos(angle) * pn[0] + Math.Sin(angle) * pn2[0]);
            pr[1] = normn * (Math.Cos(angle) * pn[1] + Math.Sin(angle) * pn2[1]);
            pr[2] = normn * (Math.Cos(angle) * pn[2] + Math.Sin(angle) * pn2[2]);

            //
            //  The rotated vector is the parallel component plus the rotated
            //  component.
            //
            p2[0] = pp[0] + pr[0];
            p2[1] = pp[1] + pr[1];
            p2[2] = pp[2] + pr[2];
        }

        public static void vector_rotate_base_2d(double[] p1, double[] pb, double angle,
                ref double[] p2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    VECTOR_ROTATE_BASE_2D rotates a vector around a base point in 2D.
            //
            //  Discussion:
            //
            //    The original vector is assumed to be P1-PB, and the
            //    rotated vector is P2-PB.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    30 June 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double P1[2], the endpoint of the original vector.
            //
            //    Input, double PB[2], the location of the base point.
            //
            //    Input, double ANGLE, the angle, in radians, of the rotation to be
            //    carried out.  A positive angle rotates the vector in the
            //    counter clockwise direction.
            //
            //    Output, double P2[2], the endpoint of the rotated vector.
            //
        {
            p2[0] = pb[0] + Math.Cos(angle) * (p1[0] - pb[0])
                    - Math.Sin(angle) * (p1[1] - pb[1]);
            p2[1] = pb[1] + Math.Sin(angle) * (p1[0] - pb[0])
                          + Math.Cos(angle) * (p1[1] - pb[1]);

        }

        public static double vector_separation_2d(double[] v1, double[] v2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    VECTOR_SEPARATION_2D finds the angular separation between vectors in 2D.
            //
            //  Discussion:
            //
            //    Any two vectors lie in a plane, and are separated by a plane angle.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    17 August 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double V1[2], V2[2], the two vectors.
            //
            //    Output, double VECTOR_SEPARATION_2D, the angle between the two vectors.
            //
        {
            int DIM_NUM = 2;

            double cos_theta;
            double v1_norm;
            double v2_norm;

            v1_norm = typeMethods.r8vec_norm(DIM_NUM, v1);
            v2_norm = typeMethods.r8vec_norm(DIM_NUM, v2);

            cos_theta = typeMethods.r8vec_dot_product(DIM_NUM, v1, v2) / (v1_norm * v2_norm);

            return (typeMethods.r8_acos(cos_theta));
        }

        public static double vector_separation_3d(double[] v1, double[] v2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    VECTOR_SEPARATION_3D finds the angular separation between vectors in 3D.
            //
            //  Discussion:
            //
            //    Any two vectors lie in a plane, and are separated by a plane angle.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    17 August 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double V1[3], V2[3], the two vectors.
            //
            //    Output, double VECTOR_SEPARATION_3D, the angle between the two vectors.
            //
        {
            int DIM_NUM = 3;

            double cos_theta;
            double v1_norm;
            double v2_norm;

            v1_norm = typeMethods.r8vec_norm(DIM_NUM, v1);
            v2_norm = typeMethods.r8vec_norm(DIM_NUM, v2);

            cos_theta = typeMethods.r8vec_dot_product(DIM_NUM, v1, v2) / (v1_norm * v2_norm);

            return (typeMethods.r8_acos(cos_theta));
        }

        public static double vector_separation_nd(int dim_num, double[] v1, double[] v2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    VECTOR_SEPARATION_ND finds the angular separation between vectors in ND.
            //
            //  Discussion:
            //
            //    Any two vectors lie in a plane, and are separated by a plane angle.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 September 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the dimension of the vectors.
            //
            //    Input, double V1[DIM_NUM], V2[DIM_NUM], the two vectors.
            //
            //    Output, double VECTOR_SEPARATION_ND, the angle between the two vectors.
            //
        {
            double cos_theta;
            double v1_norm;
            double v2_norm;

            v1_norm = typeMethods.r8vec_norm(dim_num, v1);
            v2_norm = typeMethods.r8vec_norm(dim_num, v2);
            cos_theta = typeMethods.r8vec_dot_product(dim_num, v1, v2) / (v1_norm * v2_norm);

            return (typeMethods.r8_acos(cos_theta));
        }

        public static void vector_unit_nd(int dim_num, ref double[] p)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    VECTOR_UNIT_ND normalizes a vector in ND.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    29 June 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the dimension of the vector.
            //
            //    Input/output, double P[DIM_NUM], the vector to be normalized.  On output,
            //    the vector should have unit Euclidean norm.
            //    However, if the input vector has zero Euclidean norm, it is
            //    not altered.
            //
        {
            int i;
            double norm;

            norm = typeMethods.r8vec_norm(dim_num, p);

            if (norm != 0.0)
            {
                for (i = 0; i < dim_num; i++)
                {
                    p[i] = p[i] / norm;
                }
            }

        }

        public static void provec(int m, int n, double[] base_, double[] vecm, ref double[] vecn,
                ref double[] vecnm)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PROVEC projects a vector from M space into N space.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    20 July 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the dimension of the higher order space.
            //
            //    Input, int N, the dimension of the lower order space.
            //
            //    Input, double BASE[M*N].  The columns of BASE contain
            //    N vectors, each of length M, which form the basis for
            //    a space of dimension N.
            //
            //    Input, double VECM[M], is an M dimensional vector.
            //
            //    Output, double VECN[N], the projection of VECM into the
            //    lower dimensional space.  These values represent
            //    coordinates in the lower order space.
            //
            //    Output, double VECNM[M], the projection of VECM into the
            //    lower dimensional space, but using coordinates in
            //    the higher dimensional space.
            //
        {
            int i;
            int j;
            int k;
            double temp;
            //
            //  For each vector, remove all projections onto previous vectors,
            //  and then normalize.  This should result in a matrix BASE
            //  whose columns are orthonormal.
            //
            for (j = 0; j < n; j++)
            {
                for (k = 0; k < j; k++)
                {
                    temp = typeMethods.r8vec_dot_product(m, base_, base_, a1Index: +k * m, a2Index: +j * m);

                    for (i = 0; i < m; i++)
                    {
                        base_[i + j * m] = base_[i + j * m] - temp * base_[i + k * m];
                    }
                }

                temp = 0.0;
                for (i = 0; i < m; i++)
                {
                    temp = temp + Math.Pow(base_[i + j * m], 2);
                }

                temp = Math.Sqrt(temp);

                if (0.0 < temp)
                {
                    for (i = 0; i < m; i++)
                    {
                        base_[i + j * m] = base_[i + j * m] / temp;
                    }
                }
            }

            //
            //  Compute the coordinates of the projection of the vector
            //  simply by taking dot products.
            //
            for (j = 0; j < n; j++)
            {
                vecn[j] = typeMethods.r8vec_dot_product(m, vecm, base_, a2Index: +j * m);
            }

            //
            //  Compute the coordinates of the projection in terms of
            //  the original space.
            //
            for (i = 0; i < m; i++)
            {
                vecnm[i] = 0.0;
                for (j = 0; j < n; j++)
                {
                    vecnm[i] = vecnm[i] + base_[i + j * n] * vecn[j];
                }
            }

        }

    }
}