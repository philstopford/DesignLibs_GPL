using System;

namespace Burkardt.Uniform;

public static class Sector
{
    public static double[] uniform_in_sector_map ( double r1, double r2, double t1,
            double t2, int n, ref int seed )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UNIFORM_IN_SECTOR_MAP maps uniform points into a circular sector.
        //
        //  Discussion:
        //
        //    The sector lies between circles with center at 0 and radius R1 and R2,
        //    and between rays from the center at the angles T1 and T2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 August 2004
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
        //    Input, double R1, R2, the two radii.
        //
        //    Input, double T1, T2, the two angles, which should
        //    be measured in radians, with T1 < T2.
        //
        //    Input, int N, the number of points.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double UNIFORM_IN_SECTOR_MAP[2*N], the points.
        //
    {
        const int DIM_NUM = 2;

        int j;

        double[] x = new double[DIM_NUM*n];
        double[] r = new double[n];
        double[] t = new double[n];
        double[] u = new double[n];
        double[] v = new double[n];

        UniformRNG.r8vec_uniform_01 ( n, ref seed, ref u );
        UniformRNG.r8vec_uniform_01 ( n, ref seed, ref v );

        for ( j = 0; j < n; j++ )
        {
            t[j] = ( 1.0 - u[j] ) * t1 + u[j] * t2;
            r[j] = Math.Sqrt ( ( 1.0 - v[j] ) * r1 * r1 + v[j] * r2 * r2 );

            x[0+j*DIM_NUM] = r[j] * Math.Cos ( t[j] );
            x[1+j*DIM_NUM] = r[j] * Math.Sin ( t[j] );
        }

        return x;
    }
}