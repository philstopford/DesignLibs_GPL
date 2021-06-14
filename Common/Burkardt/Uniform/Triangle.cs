using System;

namespace Burkardt.Uniform
{
    public static class Triangle
    {
        public static double[] uniform_in_triangle_map1(double[] v1, double[] v2, double[] v3,
                int n, ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    UNIFORM_IN_TRIANGLE_MAP1 maps uniform points into a triangle.
            //
            //  Discussion:
            //
            //    This routine uses Turk's rule 1.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    18 August 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Greg Turk,
            //    Generating Random Points in a Triangle,
            //    in Graphics Gems,
            //    edited by Andrew Glassner,
            //    AP Professional, 1990, pages 24-28.
            //
            //  Parameters:
            //
            //    Input, double V1[2], V2[2], V3[2], the vertices.
            //
            //    Input, int N, the number of points.
            //
            //    Input/output, int &SEED, a seed for the random number generator.
            //
            //    Output, double UNIFORM_IN_TRIANGLE_MAP1[2*N], the points.
            //
        {
            int DIM_NUM = 2;

            double a;
            double b;
            double c;
            int i;
            int j;
            double[] r = new double[DIM_NUM];
            double[] x;

            x = new double[DIM_NUM * n];

            for (j = 0; j < n; j++)
            {
                UniformRNG.r8vec_uniform_01(DIM_NUM, ref seed, ref r);

                r[1] = Math.Sqrt(r[1]);

                a = 1.0 - r[1];
                b = (1.0 - r[0]) * r[1];
                c = r[0] * r[1];

                for (i = 0; i < DIM_NUM; i++)
                {
                    x[i + j * DIM_NUM] = a * v1[i] + b * v2[i] + c * v3[i];
                }
            }

            return x;
        }

        public static double[] uniform_in_triangle_map2(double[] v1, double[] v2, double[] v3,
                int n, ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    UNIFORM_IN_TRIANGLE_MAP2 maps uniform points into a triangle.
            //
            //  Discussion:
            //
            //    The triangle is defined by three vertices.
            //
            //    This routine uses Turk's rule 2.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    18 August 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Greg Turk,
            //    Generating Random Points in a Triangle,
            //    in Graphics Gems,
            //    edited by Andrew Glassner,
            //    AP Professional, 1990, pages 24-28.
            //
            //  Parameters:
            //
            //    Input, double V1[2], V2[2], V3[2], the vertices.
            //
            //    Input, int N, the number of points.
            //
            //    Input/output, int &SEED, a seed for the random number generator.
            //
            //    Output, double UNIFORM_IN_TRIANGLE_MAP2[2*N], the points.
            //
        {
            int DIM_NUM = 2;

            int i;
            int j;
            double[] r = new double[DIM_NUM];
            double[] x;

            x = new double[DIM_NUM * n];

            for (j = 0; j < n; j++)
            {
                UniformRNG.r8vec_uniform_01(DIM_NUM, ref seed, ref r);

                if (1.0 < r[0] + r[1])
                {
                    r[0] = 1.0 - r[0];
                    r[1] = 1.0 - r[1];
                }

                for (i = 0; i < DIM_NUM; i++)
                {
                    x[i + j * DIM_NUM] = (1.0 - r[0] - r[1]) * v1[i]
                                         + r[0] * v2[i]
                                         + r[1] * v3[i];
                }
            }

            return x;
        }

        public static double[] uniform_in_triangle01_map(int n, ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    UNIFORM_IN_TRIANGLE01_MAP maps uniform points into the unit triangle.
            //
            //  Discussion:
            //
            //    The triangle is defined by the three vertices (1,0), (0,1) and (0,0).
            //    Because this is a right triangle, it is easy to generate sample points.
            //    In the case of a general triangle, more care must be taken.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    19 August 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of points.
            //
            //    Input/output, int &SEED, a seed for the random number geneator.
            //
            //    Output, double UNIFORM_IN_TRIANGLE01_MAP[2*N], the points.
            //
        {
            int DIM_NUM = 2;

            int i;
            int j;
            double[] r = new double[DIM_NUM];
            double total;
            double[] x;

            x = new double[DIM_NUM * n];

            for (j = 0; j < n; j++)
            {
                UniformRNG.r8vec_uniform_01(DIM_NUM, ref seed, ref r);

                total = 0.0;
                for (i = 0; i < DIM_NUM; i++)
                {
                    total = total + r[i];
                }

                if (1.0 < total)
                {
                    for (i = 0; i < DIM_NUM; i++)
                    {
                        r[i] = 1.0 - r[i];
                    }
                }

                for (i = 0; i < DIM_NUM; i++)
                {
                    x[i + j * DIM_NUM] = r[i];
                }
            }

            return x;
        }

        public static double[] uniform_on_triangle(int n, double[] v, ref int seed )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UNIFORM_ON_TRIANGLE maps uniform points onto the boundary of a triangle.
        //
        //  Discussion:
        //
        //    The triangle is defined by the three vertices V1, V2, V3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    07 May 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int  N, the number of points.
        //
        //    Input, double V[2*3], the vertices of the triangle.
        //
        //    Input/output, int &SEED, a seed for the random
        //    number generator.
        //
        //    Output, double UNIFORM_ON_TRIANGLE[2*N], the points.
        //
        {
            int j;
            double l1;
            double l2;
            double l3;
            int m = 2;
            double r;
            double s;
            double t;
            double[] x;

            l1 = Math.Sqrt(Math.Pow(v[0 + 1 * 2] - v[0 + 0 * 2], 2)
                      + Math.Pow(v[1 + 1 * 2] - v[1 + 0 * 2], 2));

            l2 = Math.Sqrt(Math.Pow(v[0 + 2 * 2] - v[0 + 1 * 2], 2)
                           + Math.Pow(v[1 + 2 * 2] - v[1 + 1 * 2], 2));

            l3 = Math.Sqrt(Math.Pow(v[0 + 0 * 2] - v[0 + 2 * 2], 2)
                           + Math.Pow(v[1 + 0 * 2] - v[1 + 2 * 2], 2));

            x = new double[m * n];

            for (j = 0; j < n; j++)
            {
                //
                //  R can be regarded as the distance of the point on the perimeter,
                //  as measured from the origin, along the perimeter.
                //
                r = (l1 + l2 + l3) * UniformRNG.r8_uniform_01(ref seed);
                //
                //  Case 1: between V1 and V2.
                //
                if (r <= l1)
                {
                    s = (l1 - r) / l1;
                    t = r / l1;
                    x[0 + j * 2] = s * v[0 + 0 * 2] + t * v[0 + 1 * 2];
                    x[1 + j * 2] = s * v[1 + 0 * 2] + t * v[1 + 1 * 2];
                }
                //
                //  Case 2: between V2 and V3.
                //
                else if (r <= l1 + l2)
                {
                    s = (l2 - r + l1) / l2;
                    t = (r - l1) / l2;
                    x[0 + j * 2] = s * v[0 + 1 * 2] + t * v[0 + 2 * 2];
                    x[1 + j * 2] = s * v[1 + 1 * 2] + t * v[1 + 2 * 2];
                }
                //
                //  Case 3: between V3 and V1.
                //
                else
                {
                    s = (l3 - r + l1 + l2) / l3;
                    t = (r - l1 - l2) / l3;
                    x[0 + j * 2] = s * v[0 + 2 * 2] + t * v[0 + 0 * 2];
                    x[1 + j * 2] = s * v[1 + 2 * 2] + t * v[1 + 0 * 2];
                }
            }

            return x;
        }

        public static double[] uniform_on_triangle01(int n, ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    UNIFORM_ON_TRIANGLE01 maps uniform points onto the unit triangle.
            //
            //  Discussion:
            //
            //    The unit triangle is defined by the three vertices (1,0), (0,1) and (0,0).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 May 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int N, the number of points.
            //
            //    Input/output, int SEED, a seed for the random
            //    number generator.
            //
            //    Output, double X[2*N], the points.
            //
        {
            double a;
            double b;
            int j;
            int m = 2;
            double r;
            double s;
            double[] x;

            s = Math.Sqrt(2.0);

            a = 1.0 / (2.0 + s);
            b = (1.0 + s) / (2.0 + s);

            x = new double[m * n];

            for (j = 0; j < n; j++)
            {
                //
                //  R can be regarded as the distance of the point on the perimeter,
                //  as measured from the origin, along the perimeter.
                //
                r = (2.0 + s) * UniformRNG.r8_uniform_01(ref seed);
                //
                //  Case 1: between (0,0) and (1,0).
                //
                if (r <= a)
                {
                    x[0 + j * 2] = 0.0;
                    x[1 + j * 2] = r;
                }
                //
                //  Case 2: between (1,0) and (0,1).
                //
                else if (r <= b)
                {
                    x[0 + j * 2] = 1.0 - (r - a) * Math.Sqrt(2.0) / 2.0;
                    x[1 + j * 2] = 0.0 + (r - a) * Math.Sqrt(2.0) / 2.0;
                }
                //
                //  Case 3: between (0,1) and (0,0).
                //
                else
                {
                    x[0 + j * 2] = 0.0;
                    x[1 + j * 2] = 1.0 - (r - b);
                }
            }

            return x;
        }

    }
}