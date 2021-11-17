using System;

namespace Burkardt.Uniform;

public static class Annulus
{
    public static double[] uniform_in_annulus(double[] pc, double r1, double r2, int n,
            ref int seed )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UNIFORM_IN_ANNULUS samples a circular annulus.
        //
        //  Discussion:
        //
        //    A circular annulus with center PC, inner radius R1 and
        //    outer radius R2, is the set of points P so that
        //
        //      R1^2 <= (P(1)-PC(1))^2 + (P(2)-PC(2))^2 <= R2^2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Peter Shirley,
        //    Nonuniform Random Point Sets Via Warping,
        //    Graphics Gems, Volume III,
        //    edited by David Kirk,
        //    AP Professional, 1992,
        //    ISBN: 0122861663,
        //    LC: T385.G6973.
        //
        //  Parameters:
        //
        //    Input, double PC[2], the center.
        //
        //    Input, double R1, R2, the inner and outer radii.
        //
        //    Input, int N, the number of points to generate.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double UNIFORM_IN_ANNULUS[2*N], sample points.
        //
    {
        const int DIM_NUM = 2;

        int j;

        double[] p = new double[DIM_NUM * n];

        for (j = 0; j < n; j++)
        {
            double u = UniformRNG.r8_uniform_01(ref seed);
            double theta = u * 2.0 * Math.PI;
            double v = UniformRNG.r8_uniform_01(ref seed);
            double r = Math.Sqrt((1.0 - v) * r1 * r1
                                 + v * r2 * r2);

            p[0 + j * 2] = pc[0] + r * Math.Cos(theta);
            p[1 + j * 2] = pc[1] + r * Math.Sin(theta);
        }

        return p;
    }

    public static double[] uniform_in_annulus_accept(double[] pc, double r1, double r2, int n,
            ref int seed )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UNIFORM_IN_ANNULUS_ACCEPT accepts points in an annulus.
        //
        //  Discussion:
        //
        //    A circular annulus with center PC, inner radius R1 and
        //    outer radius R2, is the set of points P so that
        //
        //      R1^2 <= (P(1)-PC(1))^2 + (P(2)-PC(2))^2 <= R2^2
        //
        //    The acceptance/rejection method is used.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double PC[2], the center.
        //
        //    Input, double R1, R2, the inner and outer radii.
        //
        //    Input, int N, the number of points.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double UNIFORM_IN_ANNULUS_ACCEPT[2*N], the points.
        //
    {
        const int DIM_NUM = 2;

        int j;
        double[] u = new double[DIM_NUM];

        if (r2 <= r1)
        {
            Console.WriteLine("");
            Console.WriteLine("UNIFORM_IN_ANNULUS_ACCEPT - Fatal error!");
            Console.WriteLine("  R2 <= R1.");
            return null;
        }

        double[] p = new double[DIM_NUM * n];
        //
        //  Generate points in a square of "radius" R2.
        //  Accept those points which lie inside the circle of radius R2, and outside
        //  the circle of radius R1.
        //
        for (j = 0; j < n; j++)
        {
            int i;
            for (;;)
            {
                UniformRNG.r8vec_uniform_01(DIM_NUM, ref seed, ref u);

                for (i = 0; i < DIM_NUM; i++)
                {
                    u[i] = (2.0 * u[i] - 1.0) * r2;
                }

                double r_squared = 0.0;
                for (i = 0; i < DIM_NUM; i++)
                {
                    r_squared += u[i] * u[i];
                }

                if (r1 * r1 <= r_squared && r_squared <= r2 * r2)
                {
                    break;
                }
            }

            for (i = 0; i < DIM_NUM; i++)
            {
                p[i + j * DIM_NUM] = pc[i] + u[i];
            }

        }

        return p;
    }

    public static double[] uniform_in_annulus_sector(double[] pc, double r1, double r2,
            double theta1, double theta2, int n, ref int seed )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UNIFORM_IN_ANNULUS_SECTOR samples an annular sector in 2D.
        //
        //  Discussion:
        //
        //    An annular sector with center PC, inner radius R1 and
        //    outer radius R2, and angles THETA1, THETA2, is the set of points
        //    P so that
        //
        //      R1^2 <= (P(1)-PC(1))^2 + (P(2)-PC(2))^2 <= R2^2
        //
        //    and
        //
        //      THETA1 <= THETA ( P - PC ) <= THETA2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    03 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Peter Shirley,
        //    Nonuniform Random Point Sets Via Warping,
        //    Graphics Gems, Volume III,
        //    edited by David Kirk,
        //    AP Professional, 1992,
        //    ISBN: 0122861663,
        //    LC: T385.G6973.
        //
        //  Parameters:
        //
        //    Input, double PC[2], the center.
        //
        //    Input, double R1, R2, the inner and outer radii.
        //
        //    Input, double THETA1, THETA2, the angles.
        //
        //    Input, int N, the number of points to generate.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double UNIFORM_IN_ANNULUS_SECTOR[2*N], sample points.
        //
    {
        const int DIM_NUM = 2;

        int j;

        double[] p = new double[DIM_NUM * n];

        for (j = 0; j < n; j++)
        {
            double u = UniformRNG.r8_uniform_01(ref seed);

            double theta = (1.0 - u) * theta1
                           + u * theta2;

            double v = UniformRNG.r8_uniform_01(ref seed);

            double r = Math.Sqrt((1.0 - v) * r1 * r1
                                 + v * r2 * r2);

            p[0 + j * 2] = pc[0] + r * Math.Cos(theta);
            p[1 + j * 2] = pc[1] + r * Math.Sin(theta);
        }

        return p;
    }
}