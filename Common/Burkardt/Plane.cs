using System;
using Burkardt.Types;

namespace Burkardt.PlaneNS
{
    public static class Plane
    {
        public static void plane_normal_basis_3d(double[] pp, ref double[] pn, ref double[] pq,
                ref double[] pr)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    PLANE_NORMAL_BASIS_3D finds two perpendicular vectors in a plane in 3D.
            //
            //  Discussion:
            //
            //    The normal form of a plane in 3D is:
            //
            //      PP is a point on the plane,
            //      N is a normal vector to the plane.
            //
            //    The two vectors to be computed, PQ and PR, can be regarded as
            //    the basis of a Cartesian coordinate system for points in the plane.
            //    Any point in the plane can be described in terms of the "origin"
            //    point PP plus a weighted sum of the two vectors PQ and PR:
            //
            //      P = PP + a * PQ + b * PR.
            //
            //    The vectors PQ and PR have unit length, and are perpendicular to N
            //    and to each other.
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
            //    Input, double PP[3], a point on the plane.
            //
            //    Input, double PN[3], a normal vector to the plane.  The
            //    vector must not have zero length, but it is not necessary for PN
            //    to have unit length.
            //
            //    Output, double PQ[3], a vector of unit length, perpendicular
            //    to the vector PN and the vector PR.
            //
            //    Output, double PR[3], a vector of unit length, perpendicular
            //    to the vector PN and the vector PQ.
            //
        {
            int DIM_NUM = 3;

            int i;
            double normal_norm;
            double pr_norm;
            double[] temp;
            //
            //  Compute the length of NORMAL.
            //
            normal_norm = typeMethods.r8vec_norm(DIM_NUM, pn);

            if (normal_norm == 0.0)
            {
                Console.WriteLine("");
                Console.WriteLine("PLANE_NORMAL_BASIS_3D - Fatal error!");
                Console.WriteLine("  The normal vector is 0.");
                return;
            }

            //
            //  Find a vector PQ that is normal to PN and has unit length.
            //
            temp = typeMethods.r8vec_any_normal(DIM_NUM, pn);
            typeMethods.r8vec_copy(DIM_NUM, temp, ref pq);
            //
            //  Now just take the cross product PN x PQ to get the PR vector.
            //
            temp = typeMethods.r8vec_cross_product_3d(pn, pq);

            pr_norm = typeMethods.r8vec_norm(DIM_NUM, temp);

            for (i = 0; i < DIM_NUM; i++)
            {
                pr[i] = temp[i] / pr_norm;
            }
        }
    }
}