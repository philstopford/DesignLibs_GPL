using System;
using Burkardt.Types;

namespace Burkardt.IHSNS
{
    public class covariance
    {
        public double average, std, covc;

        public covariance(int dim_num, int n, int[] x)
        {
            const double r8_huge = 1.0E+30;
            //
            //  Find the minimum distance for each point.
            //
            double[] mindist = new double[n];

            for (int i = 0; i < n; i++ )
            {
                mindist[i] = r8_huge;
                for (int j = 0; j < n; j++ )
                {
                    if ( i != j )
                    {
                        double dist = 0.0;
                        for (int k = 0; k < dim_num; k++ )
                        {
                            dist = dist + ( ( double )
                                ( ( x[k+i*dim_num] - x[k+j*dim_num] )
                                  * ( x[k+i*dim_num] - x[k+j*dim_num] ) ) );
                        }
                        dist = Math.Sqrt ( dist );
                        if ( dist < mindist[i] )
                        {
                            mindist[i] = dist;
                        }
                    }
                }
            }
            //
            //  Find the average minimum distance.
            //
            average = average_( n, mindist );
            //
            //  Compute the standard deviation of the distances.
            //
            std = std_( n, mindist );
            //
            //  Compute the covariance.
            //
            covc = std / average;
            
        }

        double average_(int n, double[] a)
        {
            foreach (double t in a)
            {
                average += t;
            }

            average /= a.Length;

            return average;
        }

        double std_(int n, double[] a)
        {
            if ( n < 2 )
            {
                std = 0.0;
            }
            else
            {
                average = average_( n, a );

                std = 0.0;
                for (int i = 0; i < n; i++ )
                {
                    std = std + ( a[i] - average ) * ( a[i] - average );
                }
                std = Math.Sqrt ( std / ( ( double ) ( n - 1 ) ) );

            }
            return std;
        }
    }


    public static class IHS
    {
        public static int[] ihs ( int dim_num, int n, int d, ref int seed )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    IHS implements the improved distributed hypercube sampling algorithm.
        //
        //  Discussion:
        //
        //    N Points in a DIM_NUM dimensional Latin hypercube are to be selected.
        //
        //    Each of the DIM_NUM coordinate dimensions is discretized to the values
        //    1 through N.  The points are to be chosen in such a way that
        //    no two points have any coordinate value in common.  This is
        //    a standard Latin hypercube requirement, and there are many
        //    solutions.
        //
        //    This algorithm differs in that it tries to pick a solution
        //    which has the property that the points are "spread out"
        //    as evenly as possible.  It does this by determining an optimal
        //    even spacing, and using the duplication factor D to allow it
        //    to choose the best of the various options available to it.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 April 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Brian Beachkofski, Ramana Grandhi,
        //    Improved Distributed Hypercube Sampling,
        //    American Institute of Aeronautics and Astronautics Paper 2002-1274.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int N, the number of points to be generated.
        //
        //    Input, int D, the duplication factor.  This must
        //    be at least 1.  A value of 5 is reasonable.
        //
        //    Input/output, int &SEED, a seed for the random number generator.
        //
        //    Output, int IHS[DIM_NUM*N], the points.
        //
        {
           
            int[] avail = new int [ dim_num * n ];
            int[] list = new int [ d * n ];
            int[] point = new int [ dim_num * d * n ];
            int[] x = new int[dim_num*n];

            double opt = ( ( double ) n ) /
                         Math.Pow ( ( double ) n, ( double ) ( 1.0 / ( double ) dim_num ) );
            //
            //  Pick the first point.
            //
            for (int i = 0; i < dim_num; i++ )
            {
                x[i+(n-1)*dim_num] = Uniform.UniformRNG.i4_uniform_ab( 1, n, ref seed );
            }
            //
            //  Initialize AVAIL,
            //  and set an entry in a random row of each column of AVAIL to N.
            //
            for (int j = 0; j < n; j++ )
            {
                for (int i = 0; i < dim_num; i++ )
                {
                    avail[i+j*dim_num] = j + 1;
                }
            }

            for (int i = 0; i < dim_num; i++ )
            {
                avail[i+(x[i+(n-1)*dim_num]-1)*dim_num] = n;
            }
            //
            //  Main loop:
            //  Assign a value to X(1:M,COUNT) for COUNT = N-1 down to 2.
            //
            for (int count = n - 1; 2 <= count; count-- )
            {
                //
                //  Generate valid points.
                //
                for (int i = 0; i < dim_num; i++ )
                {
                    for (int k = 0; k < d; k++ )
                    {
                        for (int j = 0; j < count; j++ )
                        {
                            list[j+k*count] = avail[i+j*dim_num];
                        }
                    }

                    for (int k = count*d - 1; 0 <= k; k-- )
                    {
                        int point_index = Uniform.UniformRNG.i4_uniform_ab( 0, k, ref seed );
                        point[i+k*dim_num] = list[point_index];
                        list[point_index] = list[k];
                    }
                }
                //
                //  For each candidate, determine the distance to all the
                //  points that have already been selected, and save the minimum value.
                //
                double min_all = typeMethods.r8_huge();
                int best = 0;

                for (int k = 0; k < d * count; k++ )
                {
                    double min_can = typeMethods.r8_huge();

                    for (int j = count; j < n; j++ )
                    {

                        double dist = 0.0;
                        for (int i = 0; i < dim_num; i++ )
                        {
                            dist = dist + ( point[i+k*dim_num] - x[i+j*dim_num] )
                                      * ( point[i+k*dim_num] - x[i+j*dim_num] );
                        }
                        dist = Math.Sqrt ( dist );

                        if ( dist < min_can )
                        {
                            min_can = dist;
                        }
                    }

                    if ( Math.Abs ( min_can - opt ) < min_all )
                    {
                        min_all = Math.Abs ( min_can - opt );
                        best = k;
                    }

                }

                for (int i = 0; i < dim_num; i++ )
                {
                    x[i+(count-1)*dim_num] = point[i+best*dim_num];
                }
                //
                //  Having chosen X(*,COUNT), update AVAIL.
                //
                for (int i = 0; i < dim_num; i++ )
                {
                    for (int j = 0; j < n; j++ )
                    {
                        if ( avail[i+j*dim_num] == x[i+(count-1)*dim_num] )
                        {
                            avail[i+j*dim_num] = avail[i+(count-1)*dim_num];
                        }
                    }
                }
            }
            //
            //  For the last point, there's only one choice.
            //
            for (int i = 0; i < dim_num; i++ )
            {
                x[i+0*dim_num] = avail[i+0*dim_num];
            }

            return x;
        }        
        
    }
}