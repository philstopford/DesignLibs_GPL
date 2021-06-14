using System;

namespace Burkardt.RandomNS
{
    public static class Scale
    {
        public static void scale_from_simplex01(int dim_num, int n, double[] t, ref double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SCALE_FROM_SIMPLEX01 rescales data from a unit to non-unit simplex.
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
        //    Reuven Rubinstein,
        //    Monte Carlo Optimization, Simulation, and Sensitivity
        //    of Queueing Networks,
        //    Krieger, 1992,
        //    ISBN: 0894647644,
        //    LC: QA298.R79.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the dimension of the space.
        //
        //    Input, int N, the number of points.
        //
        //    Input, double T[DIM_NUM*(DIM_NUM+1)], the coordinates of the DIM_NUM+1
        //    points that define the simplex.  T[0:DIM_NUM-1,0] corresponds to the
        //    origin, and T[0:DIM_NUM-1,J] will be the image of the J-th unit 
        //    coordinate vector.
        //
        //    Input/output, double X[DIM_NUM*N], the data to be modified.
        //
        {
            double[] a;
            int i;
            int j;
            double[] v;

            a = new double[dim_num * dim_num];
            v = new double[dim_num];

            for (j = 0; j < dim_num; j++)
            {
                for (i = 0; i < dim_num; i++)
                {
                    a[i + j * dim_num] = t[i + j * (dim_num + 1)] - t[i + 0 * (dim_num + 1)];
                }
            }

            for (j = 0; j < n; j++)
            {

                for (i = 0; i < dim_num; i++)
                {
                    v[i] = x[i + j * dim_num];
                }

                for (i = 0; i < dim_num; i++)
                {
                    x[i + j * dim_num] = t[i + 0 * (dim_num + 1)];
                    for (j = 0; j < n; j++)
                    {
                        x[i + j * dim_num] = x[i + j * dim_num] + a[i + j * dim_num] * v[j];
                    }
                }

            }

            return;
        }

        public static void scale_to_ball01(int dim_num, int n, ref double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SCALE_TO_BALL01 translates and rescales data to fit within the unit ball.
            //
            //  Discussion:
            //
            //    Completely arbitrary input data is given.
            //
            //    The average of the data is computed, and taken as the coordinates
            //    of the center C of a sphere.  The radius R of that sphere is the
            //    distance from the center to the furthest point in the data set.
            //
            //    Then each point is transformed to the ball of center 0 and radius
            //    1 by subtracting C and dividing by R:
            //
            //      X(1:DIM_NUM,J) -> ( X(1:DIM_NUM,J) - C(1:DIM_NUM) ) / R
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    17 August 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the dimension of the space.
            //
            //    Input, int N, the number of points.
            //
            //    Input/output, double X[DIM_NUM*N], the data to be modified.
            //
        {
            int i;
            int j;
            double r;
            double scale = 0;
            double[] xave;
            //
            //  Determine the center.
            //
            xave = new double[dim_num];

            for (i = 0; i < dim_num; i++)
            {
                xave[i] = 0.0;
                for (j = 0; j < n; j++)
                {
                    xave[i] = xave[i] + x[i + j * dim_num];
                }

                xave[i] = xave[i] / (double) n;
            }

            //
            //  Determine the maximum distance of any point from the center.
            //
            for (j = 0; j < n; j++)
            {
                r = 0.0;
                for (i = 0; i < dim_num; i++)
                {
                    r = r + Math.Pow(x[i + j * dim_num] - xave[i], 2);
                }

                if (scale < r)
                {
                    scale = r;
                }
            }

            scale = Math.Sqrt(scale);

            if (0.0 < scale)
            {
                for (j = 0; j < n; j++)
                {
                    for (i = 0; i < dim_num; i++)
                    {
                        x[i + j * dim_num] = (x[i + j * dim_num] - xave[i]) / scale;
                    }
                }
            }
            else
            {
                for (j = 0; j < n; j++)
                {
                    for (i = 0; i < dim_num; i++)
                    {
                        x[i + j * dim_num] = 0.0;
                    }
                }
            }
        }

        public static void scale_to_block01(int dim_num, int n, ref double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SCALE_TO_BLOCK01 translates and rescales data to fit in the unit block.
            //
            //  Discussion:
            //
            //    The minimum and maximum coordinate values M1(I) and M2(I) are
            //    determined, and the maximum of M2(I) - M1(I) is used to scale
            //    all the coordinates by the same factor.
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
            //  Parameters:
            //
            //    Input, int DIM_NUM, the dimension of the space.
            //
            //    Input, int N, the number of points.
            //
            //    Input/output, double X[DIM_NUM*N], the data to be modified.
            //
        {
            int i;
            int j;
            double[] xmax;
            double[] xmin;
            double xrange;
            double xrange2;

            xmax = new double[dim_num];
            xmin = new double[dim_num];
            //
            //  Determine the extremes in each dimension.
            //
            xrange = 0.0;
            for (i = 0; i < dim_num; i++)
            {
                xmin[i] = x[i + 0 * dim_num];
                xmax[i] = x[i + 0 * dim_num];
                for (j = 1; j < n; j++)
                {
                    xmin[i] = Math.Min(xmin[i], x[i + j * dim_num]);
                    xmax[i] = Math.Max(xmax[i], x[i + j * dim_num]);
                }

                xrange = Math.Max(xrange, xmax[i] - xmin[i]);
            }

            //
            //  Extend all the extremes so that the range is the same in each dimension.
            //
            for (i = 0; i < dim_num; i++)
            {
                xrange2 = xrange - (xmax[i] - xmin[i]);
                xmax[i] = xmax[i] + 0.5 * xrange2;
                xmin[i] = xmin[i] - 0.5 * xrange2;
            }

            //
            //  Now map the data to [0,1], using a single dilation factor
            //  for all dimensions.
            //
            if (0.0 == xrange)
            {
                for (j = 0; j < n; j++)
                {
                    for (i = 0; i < dim_num; i++)
                    {
                        x[i + j * dim_num] = 0.5;
                    }
                }
            }
            else
            {
                for (j = 0; j < n; j++)
                {
                    for (i = 0; i < dim_num; i++)
                    {
                        x[i + j * dim_num] = (x[i + j * dim_num] - xmin[i]) / xrange;
                    }
                }
            }
        }

        public static void scale_to_cube01(int dim_num, int n, ref double[] x)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    SCALE_TO_CUBE01 translates and rescales data to the unit hypercube.
            //
            //  Discussion:
            //
            //    In each coordinate dimension I, the minimum and maximum coordinate
            //    values M1(I) and M2(I) are determined.
            //
            //    Then, in each coordinate, the points are rescaled as
            //
            //      X(I) -> ( X(I) - M1(I) ) / ( M2(I) - M1(I) ).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    17 August 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the dimension of the space.
            //
            //    Input, int N, the number of points.
            //
            //    Input/output, double X[DIM_NUM*N], the data to be modified.
            //
        {
            int i;
            int j;
            double xmax;
            double xmin;

            for (i = 0; i < dim_num; i++)
            {
                xmin = x[i + 0 * dim_num];
                xmax = x[i + 0 * dim_num];
                for (j = 1; j < n; j++)
                {
                    if (x[i + j * dim_num] < xmin)
                    {
                        xmin = x[i + j * dim_num];
                    }

                    if (xmax < x[i + j * dim_num])
                    {
                        xmax = x[i + j * dim_num];
                    }
                }

                if (0.0 < xmax - xmin)
                {
                    for (j = 0; j < n; j++)
                    {
                        x[i + j * dim_num] = (x[i + j * dim_num] - xmin) / (xmax - xmin);
                    }
                }
                else
                {
                    for (j = 0; j < n; j++)
                    {
                        x[i + j * dim_num] = 0.0;
                    }
                }
            }
        }
    }
}