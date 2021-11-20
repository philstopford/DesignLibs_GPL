using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.Polygon;

public static class MonteCarlo
{
    public static double polygon_area(int nv, double[] v)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_AREA determines the area of a polygon.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NV, the number of vertices of the polygon.
        //
        //    Input, double V[2*N], the vertex coordinates.
        //
        //    Output, double POLYGON_AREA, the area of the polygon.
        //
    {
        int[] e = new int[2];

        e[0] = 0;
        e[1] = 0;

        double area = polygon_monomial_integral(nv, v, e);

        return area;
    }

    public static double polygon_monomial_integral(int nv, double[] v, int[] e)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_INTEGRAL integrates a monomial over a polygon.
        //
        //  Discussion:
        //
        //    Nu(P,Q) = Integral ( x, y in polygon ) x^e(1) y^e(2) dx dy
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    09 May 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Carsten Steger,
        //    On the calculation of arbitrary moments of polygons,
        //    Technical Report FGBV-96-05,
        //    Forschungsgruppe Bildverstehen, Informatik IX,
        //    Technische Universitaet Muenchen, October 1996.
        //
        //  Parameters:
        //
        //    Input, int NV, the number of vertices of the polygon.
        //
        //    Input, double V[2*N], the vertex coordinates.
        //
        //    Input, int E[2], the exponents of the monomial.
        //
        //    Output, double POLYGON_MONOMIAL_INTEGRAL, the unnormalized moment Nu(P,Q).
        //
    {
        int i;

        int p = e[0];
        int q = e[1];

        double nu_pq = 0.0;

        double xj = v[0 + (nv - 1) * 2];
        double yj = v[1 + (nv - 1) * 2];

        for (i = 0; i < nv; i++)
        {
            double xi = v[0 + i * 2];
            double yi = v[1 + i * 2];

            double s_pq = 0.0;

            int k;
            for (k = 0; k <= p; k++)
            {
                int l;
                for (l = 0; l <= q; l++)
                {
                    s_pq += typeMethods.r8_choose(k + l, l) * typeMethods.r8_choose(p + q - k - l, q - l)
                                                            * Math.Pow(xi, k) * Math.Pow(xj, p - k)
                                                            * Math.Pow(yi, l) * Math.Pow(yj, q - l);
                }
            }

            nu_pq += (xj * yi - xi * yj) * s_pq;

            xj = xi;
            yj = yi;
        }

        nu_pq = nu_pq / (p + q + 2)
                      / (p + q + 1)
                      / typeMethods.r8_choose(p + q, p);

        return nu_pq;
    }

    public static double[] polygon_sample(int nv, double[] v, int n, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_SAMPLE uniformly samples a polygon.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 May 2014
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
        //    counterclockwise order.
        //
        //    Input, int N, the number of points to create.
        //
        //    Input/output, int[] SEED, a seed for the random
        //    number generator.
        //
        //    Output, double POLYGON_SAMPLE[2*N], the points.
        //
    {
        int i;
        int j;
        //
        //  Triangulate the polygon.
        //
        double[] x = new double[nv];
        double[] y = new double[nv];
        for (i = 0; i < nv; i++)
        {
            x[i] = v[0 + i * 2];
            y[i] = v[1 + i * 2];
        }

        int[] triangles = polygon_triangulate(nv, x, y);
        //
        //  Determine the areas of each triangle.
        //
        double[] area_triangle = new double[nv - 2];

        for (i = 0; i < nv - 2; i++)
        {
            area_triangle[i] = Properties.triangle_area(
                v[0 + triangles[0 + i * 3] * 2], v[1 + triangles[0 + i * 3] * 2],
                v[0 + triangles[1 + i * 3] * 2], v[1 + triangles[1 + i * 3] * 2],
                v[0 + triangles[2 + i * 3] * 2], v[1 + triangles[2 + i * 3] * 2]);
        }

        //
        //  Normalize the areas.
        //
        double area_polygon = typeMethods.r8vec_sum(nv - 2, area_triangle);

        double[] area_relative = new double[nv - 2];

        for (i = 0; i < nv - 2; i++)
        {
            area_relative[i] = area_triangle[i] / area_polygon;
        }

        //
        //  Replace each area by the sum of itself and all previous ones.
        //
        double[] area_cumulative = new double[nv - 2];

        area_cumulative[0] = area_relative[0];
        for (i = 1; i < nv - 2; i++)
        {
            area_cumulative[i] = area_relative[i] + area_cumulative[i - 1];
        }

        double[] s = new double[2 * n];

        for (j = 0; j < n; j++)
        {
            //
            //  Choose triangle I at random, based on areas.
            //
            double area_percent = UniformRNG.r8_uniform_01(ref seed);

            int k;
            for (k = 0; k < nv - 2; k++)
            {
                i = k;

                if (area_percent <= area_cumulative[k])
                {
                    break;
                }
            }

            //
            //  Now choose a point at random in triangle I.
            //
            double[] r = UniformRNG.r8vec_uniform_01_new(2, ref seed);

            switch (r[0] + r[1])
            {
                case > 1.0:
                    r[0] = 1.0 - r[0];
                    r[1] = 1.0 - r[1];
                    break;
            }

            s[0 + j * 2] = (1.0 - r[0] - r[1]) * v[0 + triangles[0 + i * 3] * 2]
                           + r[0] * v[0 + triangles[1 + i * 3] * 2]
                           + r[1] * v[0 + triangles[2 + i * 3] * 2];

            s[1 + j * 2] = (1.0 - r[0] - r[1]) * v[1 + triangles[0 + i * 3] * 2]
                           + r[0] * v[1 + triangles[1 + i * 3] * 2]
                           + r[1] * v[1 + triangles[2 + i * 3] * 2];
        }

        return s;
    }

    public static int[] polygon_triangulate(int n, double[] x, double[] y)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYGON_TRIANGULATE determines a triangulation of a polygon.
        //
        //  Discussion:
        //
        //    There are N-3 triangles in the triangulation.
        //
        //    For the first N-2 triangles, the first edge listed is always an
        //    internal diagonal.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 May 2014
        //
        //  Author:
        //
        //    Original C version by Joseph ORourke.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Joseph ORourke,
        //    Computational Geometry in C,
        //    Cambridge, 1998,
        //    ISBN: 0521649765,
        //    LC: QA448.D38.
        //
        //  Parameters:
        //
        //    Input, int N, the number of vertices.
        //
        //    Input, double X[N], Y[N], the coordinates of each vertex.
        //
        //    Output, int TRIANGLES[3*(N-2)], the triangles of the 
        //    triangulation.
        //
    {
        int i1;
        int i3;
        int node;
        switch (n)
        {
            //
            //  We must have at least 3 vertices.
            //
            case < 3:
                Console.WriteLine("");
                Console.WriteLine("POLYGON_TRIANGULATE - Fatal error!");
                Console.WriteLine("  N < 3.");
                return null;
        }

        //
        //  Consecutive vertices cannot be equal.
        //
        int node_m1 = n - 1;
        for (node = 0; node < n; node++)
        {
            if (Math.Abs(x[node_m1] - x[node]) <= double.Epsilon && Math.Abs(y[node_m1] - y[node]) <= double.Epsilon)
            {
                Console.WriteLine("");
                Console.WriteLine("POLYGON_TRIANGULATE - Fatal error!");
                Console.WriteLine("  Two consecutive nodes are identical.");
                return null;
            }

            node_m1 = node;
        }

        //
        //  Area must be positive.
        //
        double area = 0.0;
        for (node = 0; node < n - 2; node++)
        {
            area += 0.5 *
                    (
                        (x[node + 1] - x[node]) * (y[node + 2] - y[node])
                        - (x[node + 2] - x[node]) * (y[node + 1] - y[node])
                    );
        }

        switch (area)
        {
            case <= 0.0:
                Console.WriteLine("");
                Console.WriteLine("POLYGON_TRIANGULATE - Fatal error!");
                Console.WriteLine("  Polygon has zero or negative area.");
                return null;
        }

        int[] triangles = new int[3 * (n - 2)];
        //
        //  PREV and NEXT point to the previous and next nodes.
        //
        int[] prev = new int[n];
        int[] next = new int[n];

        int i = 0;
        prev[i] = n - 1;
        next[i] = 1;

        for (i = 1; i < n - 1; i++)
        {
            prev[i] = i - 1;
            next[i] = i + 1;
        }

        i = n - 1;
        prev[i] = i - 1;
        next[i] = 0;
        //
        //  EAR indicates whether the node and its immediate neighbors form an ear
        //  that can be sliced off immediately.
        //
        bool[] ear = new bool[n];
        for (i = 0; i < n; i++)
        {
            ear[i] = Properties.diagonal(prev[i], next[i], n, prev, next, x, y);
        }

        int triangle_num = 0;

        int i2 = 0;

        while (triangle_num < n - 3)
        {
            switch (ear[i2])
            {
                //
                //  If I2 is an ear, gather information necessary to carry out
                //  the slicing operation and subsequent "healing".
                //
                case true:
                    i3 = next[i2];
                    int i4 = next[i3];
                    i1 = prev[i2];
                    int i0 = prev[i1];
                    //
                    //  Make vertex I2 disappear.
                    //
                    next[i1] = i3;
                    prev[i3] = i1;
                    //
                    //  Update the earity of I1 and I3, because I2 disappeared.
                    //
                    ear[i1] = Properties.diagonal(i0, i3, n, prev, next, x, y);
                    ear[i3] = Properties.diagonal(i1, i4, n, prev, next, x, y);
                    //
                    //  Add the diagonal [I3, I1, I2] to the list.
                    //
                    triangles[0 + triangle_num * 3] = i3;
                    triangles[1 + triangle_num * 3] = i1;
                    triangles[2 + triangle_num * 3] = i2;
                    triangle_num += 1;
                    break;
            }

            //
            //  Try the next vertex.
            //
            i2 = next[i2];
        }

        //
        //  The last triangle is formed from the three remaining vertices.
        //
        i3 = next[i2];
        i1 = prev[i2];

        triangles[0 + triangle_num * 3] = i3;
        triangles[1 + triangle_num * 3] = i1;
        triangles[2 + triangle_num * 3] = i2;
        triangle_num += 1;

        return triangles;
    }
}