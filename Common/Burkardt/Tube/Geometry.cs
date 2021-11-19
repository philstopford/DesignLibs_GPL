using System;
using System.Linq;
using Burkardt.Geometry;
using Burkardt.Types;

namespace Burkardt.Tube;

public static class Geometry
{
    public static void tube_2d(double dist, int n, double[] p, ref double[] p1, ref double[] p2)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TUBE_2D constructs a "tube" of given width around a path in 2D.
        //
        //  Discussion:
        //
        //    The routine is given a sequence of N points, and a distance DIST.
        //
        //    It returns the coordinates of the corners of the top and bottom
        //    of a tube of width 2*DIST, which envelopes the line connecting
        //    the points.
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
        //    Input, double DIST, the radius of the tube.
        //
        //    Input, int N, the number of points defining the line.
        //    N must be at least 2.
        //
        //    Input, double P[2*N], the points which comprise the broken
        //    line which is to be surrounded by the tube.  Points should
        //    not be immediately repeated.
        //
        //    Output, double P1[N], P2[N], the points P1
        //    form one side of the tube, and P2 the other.
        //
    {
        const int DIM_NUM = 2;

        int j;
        double[] p4 = new double[DIM_NUM];
        double[] p5 = new double[DIM_NUM];
        switch (n)
        {
            //
            //  Check that N is at least 3.
            //
            case < 3:
                Console.WriteLine("");
                Console.WriteLine("TUBE_2D - Fatal error!");
                Console.WriteLine("  N must be at least 3");
                Console.WriteLine("  but your input value was N = " + n + "");
                return;
        }

        //
        //  Check that consecutive points are distinct.
        //
        for (j = 0; j < n - 1; j++)
        {
            if (!typeMethods.r8vec_eq(DIM_NUM, p, p, startIndexA1: +j * DIM_NUM, startIndexA2: +(j + 1) * DIM_NUM))
            {
                continue;
            }

            Console.WriteLine("");
            Console.WriteLine("TUBE_2D - Fatal error!");
            Console.WriteLine("  P[1:2,J] = P[1:2,J+1] for J = " + j + "");
            Console.WriteLine("  P[0,J] = " + p[0 + j * DIM_NUM] + "");
            Console.WriteLine("  P[1,J] = " + p[1 + j * DIM_NUM] + "");
            return;
        }

        for (j = 1; j <= n; j++)
        {
            double[] pim1 = j switch
            {
                1 => p.Skip(+(j - 1) * DIM_NUM).ToArray(),
                _ => p.Skip(+(j - 2) * DIM_NUM).ToArray()
            };

            double[] pi = p.Skip(+(j - 1) * DIM_NUM).ToArray();

            double[] pip1 = j < n ? p.Skip(+j * DIM_NUM).ToArray() : p.Skip(+(j - 1) * DIM_NUM).ToArray();

            Angle.angle_box_2d(dist, pim1, pi, pip1, ref p4, ref p5);

            p1[0 + (j - 1) * DIM_NUM] = p4[0];
            p1[1 + (j - 1) * DIM_NUM] = p4[1];
            p2[0 + (j - 1) * DIM_NUM] = p5[0];
            p2[1 + (j - 1) * DIM_NUM] = p5[1];
            double temp;
            switch (j)
            {
                //
                //  On the first and last steps, translate the corner points DIST units
                //  along the line, to make an extra buffer.
                //
                case 1:
                    temp = Math.Sqrt(Math.Pow(p[0 + 1 * DIM_NUM] - p[0 + 0 * DIM_NUM], 2)
                                     + Math.Pow(p[1 + 1 * DIM_NUM] - p[1 + 0 * DIM_NUM], 2));

                    p1[0 + 0 * DIM_NUM] -= dist * (p[0 + 1 * DIM_NUM] - p[0 + 0 * DIM_NUM]) / temp;
                    p1[1 + 0 * DIM_NUM] -= dist * (p[1 + 1 * DIM_NUM] - p[1 + 0 * DIM_NUM]) / temp;
                    p2[0 + 0 * DIM_NUM] -= dist * (p[0 + 1 * DIM_NUM] - p[0 + 0 * DIM_NUM]) / temp;
                    p2[1 + 0 * DIM_NUM] -= dist * (p[1 + 1 * DIM_NUM] - p[1 + 0 * DIM_NUM]) / temp;
                    break;
                default:
                {
                    if (j == n)
                    {
                        temp = Math.Sqrt(Math.Pow(p[0 + (n - 1) * DIM_NUM] - p[0 + (n - 2) * DIM_NUM], 2)
                                         + Math.Pow(p[1 + (n - 1) * DIM_NUM] - p[1 + (n - 2) * DIM_NUM], 2));

                        p1[0 + (n - 1) * DIM_NUM] += dist * (p[0 + (n - 1) * DIM_NUM] - p[0 + (n - 2) * DIM_NUM]) / temp;
                        p1[1 + (n - 1) * DIM_NUM] += dist * (p[1 + (n - 1) * DIM_NUM] - p[1 + (n - 2) * DIM_NUM]) / temp;
                        p2[0 + (n - 1) * DIM_NUM] += dist * (p[0 + (n - 1) * DIM_NUM] - p[0 + (n - 2) * DIM_NUM]) / temp;
                        p2[1 + (n - 1) * DIM_NUM] += dist * (p[1 + (n - 1) * DIM_NUM] - p[1 + (n - 2) * DIM_NUM]) / temp;
                    }

                    break;
                }
            }

            switch (j)
            {
                //
                //  The new points may need to be swapped.
                //
                //  Compute the signed distance from the points to the line.
                //
                case > 1:
                {
                    double a = p[1 + (j - 2) * DIM_NUM] - p[1 + (j - 1) * DIM_NUM];
                    double b = p[0 + (j - 1) * DIM_NUM] - p[0 + (j - 2) * DIM_NUM];
                    double c = p[0 + (j - 2) * DIM_NUM] * p[1 + (j - 1) * DIM_NUM]
                               - p[0 + (j - 1) * DIM_NUM] * p[1 + (j - 2) * DIM_NUM];

                    double dis1 = (a * p1[0 + (j - 2) * DIM_NUM] + b * p1[1 + (j - 2) * DIM_NUM] + c)
                                  / Math.Sqrt(a * a + b * b);

                    double dis2 = (a * p1[0 + (j - 1) * DIM_NUM] + b * p1[1 + (j - 1) * DIM_NUM] + c)
                                  / Math.Sqrt(a * a + b * b);

                    if (Math.Abs(typeMethods.r8_sign(dis1) - typeMethods.r8_sign(dis2)) > double.Epsilon)
                    {
                        temp = p1[0 + (j - 1) * DIM_NUM];
                        p1[0 + (j - 1) * DIM_NUM] = p2[0 + (j - 1) * DIM_NUM];
                        p2[0 + (j - 1) * DIM_NUM] = temp;
                        temp = p1[1 + (j - 1) * DIM_NUM];
                        p1[1 + (j - 1) * DIM_NUM] = p2[1 + (j - 1) * DIM_NUM];
                        p2[1 + (j - 1) * DIM_NUM] = temp;
                    }

                    break;
                }
            }
        }

    }
}