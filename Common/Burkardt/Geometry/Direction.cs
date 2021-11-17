using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.Geometry;

public static class Direction
{
    public static double[] direction_pert_3d(double sigma, double[] vbase, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIRECTION_PERT_3D randomly perturbs a direction vector in 3D.
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
        //    Input, double SIGMA, determines the strength of the perturbation.
        //    SIGMA <= 0 results in a completely random direction.
        //    1 <= SIGMA results in VBASE.
        //    0 < SIGMA < 1 results in a perturbation from VBASE, which is
        //    large when SIGMA is near 0, and small when SIGMA is near 1.
        //
        //    Input, double VBASE[3], the base direction vector, which should have
        //    unit norm.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double DIRECTION_PERT_3D[3], the perturbed vector, which will
        //    have unit norm.
        //
    {
        int DIM_NUM = 3;

        double dphi;
        double[] p = new double[DIM_NUM];
        double phi;
        double psi;
        double theta;
        double vdot;
        double[] vran;
        double x;
        //
        //  0 <= SIGMA, just use the base vector.
        //
        vran = new double[DIM_NUM];

        switch (sigma)
        {
            case >= 1.0:
                typeMethods.r8vec_copy(DIM_NUM, vbase, ref vran);
                break;
            case <= 0.0:
                vdot = UniformRNG.r8_uniform_01(ref seed);
                vdot = 2.0 * vdot - 1.0;
                phi = Math.Acos(vdot);
                theta = UniformRNG.r8_uniform_01(ref seed);
                theta = 2.0 * Math.PI * theta;

                vran[0] = Math.Cos(theta) * Math.Sin(phi);
                vran[1] = Math.Sin(theta) * Math.Sin(phi);
                vran[2] = Math.Cos(phi);
                break;
            default:
                phi = Math.Acos(vbase[2]);
                theta = Math.Atan2(vbase[1], vbase[0]);
                //
                //  Pick VDOT, which must be between -1 and 1.  This represents
                //  the dot product of the perturbed vector with the base vector.
                //
                //  UniformRNG.r8_uniform_01 returns a uniformly random value between 0 and 1.
                //  The operations we perform on this quantity tend to bias it
                //  out towards 1, as SIGMA grows from 0 to 1.
                //
                //  VDOT, in turn, is a value between -1 and 1, which, for large
                //  SIGMA, we want biased towards 1.
                //
                x = UniformRNG.r8_uniform_01(ref seed);
                x = Math.Exp((1.0 - sigma) * Math.Log(x));
                vdot = 2.0 * x - 1.0;
                dphi = Math.Acos(vdot);
                //
                //  Now we know enough to write down a vector that is rotated DPHI
                //  from the base vector.
                //
                p[0] = Math.Cos(theta) * Math.Sin(phi + dphi);
                p[1] = Math.Sin(theta) * Math.Sin(phi + dphi);
                p[2] = Math.Cos(phi + dphi);
                //
                //  Pick a uniformly random rotation between 0 and 2 Pi around the
                //  axis of the base vector.
                //
                psi = UniformRNG.r8_uniform_01(ref seed);
                psi = 2.0 * Math.PI * psi;

                Vector.Geometry.vector_rotate_3d(p, vbase, psi, ref vran);
                break;
        }

        return vran;
    }

    public static double[] direction_uniform_2d(ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIRECTION_UNIFORM_2D picks a random direction vector in 2D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    16 June 2002
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double DIRECTION_UNIFORM_2D[2], the random direction vector,
        //    with unit norm.
        //
    {
        int DIM_NUM = 2;

        double theta;
        double[] vran;

        theta = UniformRNG.r8_uniform_01(ref seed);
        theta = 2.0 * Math.PI * theta;

        vran = new double[DIM_NUM];

        vran[0] = Math.Cos(theta);
        vran[1] = Math.Sin(theta);

        return vran;
    }

    public static double[] direction_uniform_3d(ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIRECTION_UNIFORM_3D picks a random direction vector in 3D.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double DIRECTION_UNIFORM_3D[3], the random direction vector,
        //    with unit norm.
        //
    {
        int DIM_NUM = 3;

        double phi;
        double theta;
        double vdot;
        double[] vran;
        //
        //  Pick a uniformly random VDOT, which must be between -1 and 1.
        //  This represents the dot product of the random vector with the Z unit vector.
        //
        //  This works because the surface area of the sphere between
        //  Z and Z + dZ is independent of Z.  So choosing Z uniformly chooses
        //  a patch of area uniformly.
        //
        vdot = UniformRNG.r8_uniform_01(ref seed);
        vdot = 2.0 * vdot - 1.0;
        phi = Math.Acos(vdot);
        //
        //  Pick a uniformly random rotation between 0 and 2 Pi around the
        //  axis of the Z vector.
        //
        theta = UniformRNG.r8_uniform_01(ref seed);
        theta = 2.0 * Math.PI * theta;

        vran = new double[DIM_NUM];

        vran[0] = Math.Cos(theta) * Math.Sin(phi);
        vran[1] = Math.Sin(theta) * Math.Sin(phi);
        vran[2] = Math.Cos(phi);

        return vran;
    }

    public static double[] direction_uniform_nd(int dim_num, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DIRECTION_UNIFORM_ND generates a random direction vector in ND.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the dimension of the space.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double DIRECTION_UNIFORM_ND[DIM_NUM], a random direction vector, 
        //    with unit norm.
        //
    {
        double[] a;
        typeMethods.r8vecNormalData data = new();
        //
        //  Take DIM_NUM random samples from the normal distribution.
        //
        a = typeMethods.r8vec_normal_01_new(dim_num, ref data, ref seed);
        //
        //  Normalize the vector.
        //
        Vector.Geometry.vector_unit_nd(dim_num, ref a);

        return a;
    }
}