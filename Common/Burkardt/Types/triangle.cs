using System;
using Burkardt.Uniform;

namespace Burkardt.Types
{
    public static partial class typeMethods
    {
        public static void triangle_reference_sample(int n, ref int seed, ref double[] p)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TRIANGLE_REFERENCE_SAMPLE returns random points in the reference triangle.
            //
            //  Diagram:
            //
            //       3
            //    s  |.
            //    i  | .
            //    d  |  .
            //    e  |   .  side 2
            //       |    .
            //    3  |     .
            //       |      .
            //       1-------2
            //
            //         side 1
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    07 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, integer N, the number of points to sample.
            //
            //    Input/output, int *SEED, a seed for the random number generator.
            //
            //    Output, double P[2*N], a random point in the triangle.
            //
        {
            int DIM_NUM = 2;

            double alpha;
            double beta;
            int j;
            double r;

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
                //  ALPHA = SQRT ( R ) of the distance from vertex 1 to vertices 2 and 3.
                //
                alpha = Math.Sqrt(r);
                //
                //  Now choose, uniformly at random, a point on the line L.
                //
                beta = UniformRNG.r8_uniform_01(ref seed);

                p[0 + j * 2] = (1.0 - beta) * alpha;
                p[1 + j * 2] = beta * alpha;
            }
        }

        public static void triangle_sample(double[] t, int n, ref int seed, ref double[] p, int pIndex = 0)

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
            //    06 December 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double T[2*3], the triangle vertices.
            //
            //    Input, integer N, the number of points to sample.
            //
            //    Input/output, int *SEED, a seed for the random number generator.
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
                //  ALPHA = SQRT ( R ) of the distance from vertex 1 to vertices 2 and 3.
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

                p[pIndex + (0 + j * 2)] = (1.0 - beta) * p12[0] + beta * p13[0];
                p[pIndex + (1 + j * 2)] = (1.0 - beta) * p12[1] + beta * p13[1];
            }
        }

        public static double[] triangle_angles_2d_new ( double[] t )

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
            double[] angle;
            double a;
            double b;
            double c;
            double pi = 3.141592653589793;

            angle = new double[3];

            a = Math.Sqrt ( Math.Pow ( t[0+1*2] - t[0+0*2], 2 )
                       + Math.Pow ( t[1+1*2] - t[1+0*2], 2 ) );

            b = Math.Sqrt ( Math.Pow ( t[0+2*2] - t[0+1*2], 2 )
                            + Math.Pow ( t[1+2*2] - t[1+1*2], 2 ) );

            c = Math.Sqrt ( Math.Pow ( t[0+0*2] - t[0+2*2], 2 )
                            + Math.Pow ( t[1+0*2] - t[1+2*2], 2 ) );
            //
            //  Take care of a ridiculous special case.
            //
            if ( a == 0.0 && b == 0.0 && c == 0.0 )
            {
                angle[0] = 2.0 * pi / 3.0;
                angle[1] = 2.0 * pi / 3.0;
                angle[2] = 2.0 * pi / 3.0;
                return angle;
            }

            if ( c == 0.0 || a == 0.0 )
            {
                angle[0] = pi;
            }
            else
            {
                angle[0] = Helpers.arc_cosine ( ( c * c + a * a - b * b ) / ( 2.0 * c * a ) );
            }

            if ( a == 0.0 || b == 0.0 )
            {
                angle[1] = pi;
            }
            else
            {
                angle[1] = Helpers.arc_cosine ( ( a * a + b * b - c * c ) / ( 2.0 * a * b ) );
            }

            if ( b == 0.0 || c == 0.0 )
            {
                angle[2] = pi;
            }
            else
            {
                angle[2] = Helpers.arc_cosine ( ( b * b + c * c - a * a ) / ( 2.0 * b * c ) );
            }

            return angle;
        }
        
        public static double triangle_area_2d ( double[] t )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TRIANGLE_AREA_2D computes the area of a triangle in 2D.
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
        //    Input, double T[2*3], the vertices of the triangle.
        //
        //    Output, double TRIANGLE_AREA_2D, the area of the triangle.  AREA will
        //    be nonnegative.
        //
        {
            double area = Math.Abs ( 0.5 * (
                t[0+0*2] * ( t[1+2*2] - t[1+1*2] ) +
                t[0+1*2] * ( t[1+0*2] - t[1+2*2] ) +
                t[0+2*2] * ( t[1+1*2] - t[1+0*2] ) ) );

            return area;
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
            //    Output, double *TRIANGLE_CIRCUMCENTER_2D[2], the circumcenter of
            //    the triangle.
            //
        {
            int DIM_NUM = 2;

            double asq;
            double bot;
            double[] center;
            double csq;
            double top1;
            double top2;

            center = new double[DIM_NUM];

            asq = (t[0 + 1 * 2] - t[0 + 0 * 2]) * (t[0 + 1 * 2] - t[0 + 0 * 2])
                  + (t[1 + 1 * 2] - t[1 + 0 * 2]) * (t[1 + 1 * 2] - t[1 + 0 * 2]);

            csq = (t[0 + 2 * 2] - t[0 + 0 * 2]) * (t[0 + 2 * 2] - t[0 + 0 * 2])
                  + (t[1 + 2 * 2] - t[1 + 0 * 2]) * (t[1 + 2 * 2] - t[1 + 0 * 2]);

            top1 = (t[1 + 1 * 2] - t[1 + 0 * 2]) * csq - (t[1 + 2 * 2] - t[1 + 0 * 2]) * asq;
            top2 = -(t[0 + 1 * 2] - t[0 + 0 * 2]) * csq + (t[0 + 2 * 2] - t[0 + 0 * 2]) * asq;

            bot = (t[1 + 1 * 2] - t[1 + 0 * 2]) * (t[0 + 2 * 2] - t[0 + 0 * 2])
                  - (t[1 + 2 * 2] - t[1 + 0 * 2]) * (t[0 + 1 * 2] - t[0 + 0 * 2]);

            center[0] = t[0 + 0 * 2] + 0.5 * top1 / bot;
            center[1] = t[1 + 0 * 2] + 0.5 * top2 / bot;

            return center;
        }
    }
}
