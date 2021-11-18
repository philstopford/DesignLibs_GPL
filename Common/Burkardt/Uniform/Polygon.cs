﻿using Burkardt.Types;

namespace Burkardt.Uniform;

public static class Polygon
{
    public static double[] uniform_in_polygon_map(int nv, double[] v, int n, ref int seed )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UNIFORM_IN_POLYGON_MAP maps uniform points into a polygon.
        //
        //  Discussion:
        //
        //    If the polygon is regular, or convex, or at least star-shaped,
        //    this routine will work.
        //
        //    This routine assumes that all points between the centroid and
        //    any point on the boundary lie within the polygon.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 August 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NV, the number of vertices.
        //
        //    Input, double V[2*NV], the vertices of the polygon, listed in
        //    clockwise or counterclockwise order.
        //
        //    Input, int N, the number of points.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, double UNIFORM_IN_POLYGON_MAP[2*N], the points.
        //
    {
        const int DIM_NUM = 2;

        int i;
        int j;
        double[] r = new double[DIM_NUM];
        double[] t = new double[2 * 3];
        int triangle = 0;

        double[] area = new double[nv];
        double[] x = new double[DIM_NUM * n];
        //
        //  Find the centroid.
        //
        double[] centroid = typeMethods.polygon_centroid_2d(nv, v);
        //
        //  Determine the areas of each triangle.
        //
        double total = 0.0;
        for (i = 0; i < nv; i++)
        {
            int ip1;
            if (i < nv - 1)
            {
                ip1 = i + 1;
            }
            else
            {
                ip1 = 0;
            }

            t[0 + 0 * 2] = v[0 + i * 2];
            t[1 + 0 * 2] = v[1 + i * 2];

            t[0 + 1 * 2] = v[0 + ip1 * 2];
            t[1 + 1 * 2] = v[1 + ip1 * 2];

            t[0 + 2 * 2] = centroid[0];
            t[1 + 2 * 2] = centroid[1];

            area[i] = typeMethods.triangle_area_2d(t);

            total += area[i];
        }

        //
        //  Normalize the areas.
        //
        for (i = 0; i < nv; i++)
        {
            area[i] /= total;
        }

        //
        //  Replace each area by the sum of itself and all previous ones.
        //
        for (i = 1; i < nv; i++)
        {
            area[i] += area[i - 1];
        }

        for (j = 0; j < n; j++)
        {
            //
            //  Choose a triangle T at random, based on areas.
            //
            double u = UniformRNG.r8_uniform_01(ref seed);

            int k;
            for (k = 0; k < nv; k++)
            {
                if (!(u <= area[k]))
                {
                    continue;
                }

                triangle = k;
                break;
            }

            //
            //  Now choose a point at random in the triangle.
            //
            int triangle_p1;
            if (triangle < nv - 1)
            {
                triangle_p1 = triangle + 1;
            }
            else
            {
                triangle_p1 = 0;
            }

            UniformRNG.r8vec_uniform_01(DIM_NUM, ref seed, ref r);

            switch (r[0] + r[1])
            {
                case > 1.0:
                    r[0] = 1.0 - r[0];
                    r[1] = 1.0 - r[1];
                    break;
            }

            for (i = 0; i < DIM_NUM; i++)
            {
                x[i + j * DIM_NUM] = (1.0 - r[0] - r[1]) * v[i + DIM_NUM * triangle]
                                     + r[0] * v[i + DIM_NUM * triangle_p1]
                                     + r[1] * centroid[i];
            }

        }

        return x;
    }
}