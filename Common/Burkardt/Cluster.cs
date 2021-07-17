using System;
using Burkardt.Sampling;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt
{
    public static class Cluster
    {
        public static double cluster_energy(ref RegionData data, int dim_num, int n, double[] cell_generator,
                int sample_num_cvt, int sample_function_cvt, ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CLUSTER_ENERGY returns the energy of a dataset.
            //
            //  Discussion:
            //
            //    The energy is the integral of the square of the distance from each point
            //    in the region to its nearest generator.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    08 September 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the spatial dimension.
            //
            //    Input, int N, the number of generators.
            //
            //    Input, double CELL_GENERATOR[DIM_NUM*N], the coordinates of the points.
            //
            //    Input, int SAMPLE_NUM_CVT, the number of sample points to use.
            //
            //    Input, int SAMPLE_FUNCTION_CVT, specifies how the sampling is done.
            //    -1, 'RANDOM', using C++ RANDOM function;
            //     0, 'UNIFORM', using a simple uniform RNG;
            //     1, 'HALTON', from a Halton sequence;
            //     2, 'GRID', points from a grid;
            //     3, 'USER', call "user" routine.
            //
            //    Input/output, int *SEED, a seed for the random number generator.
            //
            //    Output, double CLUSTER_ENERGY, the estimated energy.
            //
        {
            double energy;
            int i;
            int j;
            int nearest;
            bool reset;
            double[] x;

            x = new double [dim_num];

            energy = 0.0;
            reset = true;

            for (j = 0; j < sample_num_cvt; j++)
            {
                //
                //  Generate a sampling point X.
                //
                Region.region_sampler(ref data, dim_num, 1, sample_num_cvt, x, sample_function_cvt,
                    reset, ref seed);

                reset = false;
                //
                //  Find the nearest cell generator.
                //
                nearest = find_closest(dim_num, n, x, cell_generator);

                for (i = 0; i < dim_num; i++)
                {
                    energy = energy
                             + Math.Pow(x[i] - cell_generator[i + nearest * dim_num], 2);
                }
            }

            //
            //  Add the contribution to the energy.
            //
            energy = energy / (double) (sample_num_cvt);

            return energy;
        }

        public static double[] cluster_energy_compute(int dim_num, int point_num, int cluster_num,
                double[] point, int[] cluster, double[] cluster_center)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CLUSTER_ENERGY_COMPUTE computes the energy of the clusters.
            //
            //  Discussion:
            //
            //    The cluster energy is defined as the sum of the distance
            //    squared from each point to its cluster center.  It is the goal
            //    of the H-means and K-means algorithms to find, for a fixed number
            //    of clusters, a clustering that minimizes this energy
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    07 October 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the number of spatial dimensions.
            //
            //    Input, int POINT_NUM, the number of data points.
            //
            //    Input, int CLUSTER_NUM, the number of clusters.
            //
            //    Input, double POINT[DIM_NUM*POINT_NUM], the data points.
            //
            //    Input, int CLUSTER[POINT_NUM], the cluster to which each
            //    data point belongs.  These values are 0-based.
            //
            //    Input, double CLUSTER_CENTER[DIM_NUM*CLUSTER_NUM], the 
            //    centers associated with the minimal energy clustering.
            //
            //    Output, double CLUSTER_ENERGY_COMPUTE[CLUSTER_NUM], the energy
            //    associated with each cluster.
            //
        {
            double[] cluster_energy;
            int i = 0;
            int j;
            int k;
            double point_energy;

            cluster_energy = typeMethods.r8vec_zero_new(cluster_num);

            for (j = 0; j < point_num; j++)
            {
                k = cluster[i];
                point_energy = 0.0;
                for (i = 0; i < dim_num; i++)
                {
                    point_energy = point_energy
                                   + Math.Pow(point[i + j * dim_num] - cluster_center[i + k * dim_num], 2);
                }

                cluster_energy[k] = cluster_energy[k] + point_energy;
            }

            return cluster_energy;
        }

        public static double[] cluster_initialize_1(int dim_num, int point_num, int cluster_num,
                double[] point)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CLUSTER_INITIALIZE_1 initializes the clusters to data points.
            //
            //  Discussion:
            //
            //    The cluster centers are simply chosen to be the first data points.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    08 October 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the number of spatial dimensions.
            //
            //    Input, int POINT_NUM, the number of points.
            //
            //    Input, int CLUSTER_NUM, the number of clusters.
            //
            //    Input, double POINT[DIM_NUM*POINT_NUM], the coordinates 
            //    of the points.
            //
            //    Output, double CLUSTER_INITIALIZE_1[DIM_NUM*CLUSTER_NUM],
            //    the coordinates of the cluster centers.
            //
        {
            double[] cluster_center;
            int i;
            int j;

            cluster_center = new double[dim_num * cluster_num];

            for (j = 0; j < cluster_num; j++)
            {
                for (i = 0; i < dim_num; i++)
                {
                    cluster_center[i + j * dim_num] = point[i + j * dim_num];
                }
            }

            return cluster_center;
        }

        public static double[] cluster_initialize_2(int dim_num, int point_num, int cluster_num,
                double[] point, ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CLUSTER_INITIALIZE_2 initializes the cluster centers to random values.
            //
            //  Discussion:
            //
            //    In this case, the hyperbox containing the data is computed.
            //
            //    Then the cluster centers are chosen uniformly at random within
            //    this hyperbox.
            //
            //    Of course, if the data is not smoothly distributed throughout
            //    the box, many cluster centers will be isolated.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    08 October 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the number of spatial dimensions.
            //
            //    Input, int POINT_NUM, the number of points.
            //
            //    Input, int CLUSTER_NUM, the number of clusters.
            //
            //    Input, double POINT[DIM_NUM*POINT_NUM], the coordinates 
            //    of the points.
            //
            //    Input/output, int *SEED, a seed for the random 
            //    number generator.
            //
            //    Output, double CLUSTER_CENTER[DIM_NUM*CLUSTER_NUM], 
            //    the coordinates of the cluster centers.
            //
        {
            double[] cluster_center;
            int i;
            int j;
            double[] r;
            double[] r_max;
            double[] r_min;

            cluster_center = new double[dim_num * cluster_num];

            r = new double[dim_num];
            r_min = new double[dim_num];
            r_max = new double[dim_num];

            j = 0;
            for (i = 0; i < dim_num; i++)
            {
                r_max[i] = point[i + j * dim_num];
                r_min[i] = point[i + j * dim_num];
            }

            for (j = 1; j < point_num; j++)
            {
                for (i = 0; i < dim_num; i++)
                {
                    r_max[i] = Math.Max(r_max[i], point[i + j * dim_num]);
                    r_min[i] = Math.Min(r_min[i], point[i + j * dim_num]);
                }
            }

            for (j = 0; j < cluster_num; j++)
            {
                UniformRNG.r8vec_uniform_01(dim_num, ref seed, ref r);
                for (i = 0; i < dim_num; i++)
                {
                    cluster_center[i + j * dim_num] = (1.0 - r[i]) * r_min[i] + r[i] * r_max[i];
                }
            }

            return cluster_center;
        }

        public static double[] cluster_initialize_3(int dim_num, int point_num, int cluster_num,
                double[] point, ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //   CLUSTER_INITIALIZE_3 initializes the cluster centers to random values.
            //
            //  Discussion:
            //
            //    In this case, each point is randomly assigned to a cluster, and
            //    the cluster centers are then computed as the centroids of the points 
            //    in the cluster.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    08 October 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the number of spatial dimensions.
            //
            //    Input, int POINT_NUM, the number of points.
            //
            //    Input, int CLUSTER_NUM, the number of clusters.
            //
            //    Input, double POINT[DIM_NUM*POINT_NUM], the coordinates 
            //    of the points.
            //
            //    Input/output, int SEED, a seed for the random 
            //    number generator.
            //
            //    Output, double CLUSTER_INITIALIZE_3[DIM_NUM*CLUSTER_NUM], 
            //    the coordinates of the cluster centers.
            //
        {
            double[] cluster_center;
            int[] cluster_population;
            int i;
            int j;
            int k;
            //
            //  Assign one point to each cluster center.
            //
            cluster_center = new double[dim_num * cluster_num];

            for (k = 0; k < cluster_num; k++)
            {
                for (i = 0; i < dim_num; i++)
                {
                    cluster_center[i + k * dim_num] = point[i + k * dim_num];
                }
            }

            cluster_population = new int[cluster_num];

            for (k = 0; k < cluster_num; k++)
            {
                cluster_population[k] = 1;
            }

            //
            //  The rest of the points get assigned randomly.
            //
            for (j = cluster_num; j < point_num; j++)
            {
                k = UniformRNG.i4_uniform(1, cluster_num, ref seed);
                for (i = 0; i < dim_num; i++)
                {
                    cluster_center[i + k * dim_num] = cluster_center[i + k * dim_num]
                                                      + point[i + j * dim_num];
                }

                cluster_population[k] = cluster_population[k] + 1;
            }

            //
            //  Now average the points to get the centroid.
            //
            for (k = 0; k < cluster_num; k++)
            {
                if (cluster_population[k] != 0)
                {
                    for (i = 0; i < dim_num; i++)
                    {
                        cluster_center[i + k * dim_num] = cluster_center[i + k * dim_num] /
                                                          (double) (cluster_population[k]);
                    }
                }
            }

            return cluster_center;
        }

        public static double[] cluster_initialize_4(int dim_num, int point_num, int cluster_num,
                double[] point, ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //   CLUSTER_INITIALIZE_4 initializes the cluster centers to random values.
            //
            //  Discussion:
            //
            //    In this case, each data point is divided randomly among the
            //    the cluster centers.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    08 October 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the number of spatial dimensions.
            //
            //    Input, int POINT_NUM, the number of points.
            //
            //    Input, int CLUSTER_NUM, the number of clusters.
            //
            //    Input, double POINT[DIM_NUM*POINT_NUM], the coordinates 
            //    of the points.
            //
            //    Input/output, int SEED, a seed for the random 
            //    number generator.
            //
            //    Output, double CLUSTER_INITIALIZE_4[DIM_NUM*CLUSTER_NUM], 
            //    the coordinates of the cluster centers.
            //
        {
            double[] cluster_center;
            double[] cluster_factor;
            double[] cluster_weight;
            double divisor;
            int i;
            int j;
            int k;

            cluster_center = typeMethods.r8vec_zero_new(dim_num * cluster_num);

            cluster_factor = new double[cluster_num];

            cluster_weight = typeMethods.r8vec_zero_new(cluster_num);

            for (j = 0; j < point_num; j++)
            {
                UniformRNG.r8vec_uniform_01(cluster_num, ref seed, ref cluster_factor);

                divisor = typeMethods.r8vec_sum(cluster_num, cluster_factor);

                for (k = 0; k < cluster_num; k++)
                {
                    cluster_factor[k] = cluster_factor[k] / divisor;
                }

                for (k = 0; k < cluster_num; k++)
                {
                    for (i = 0; i < dim_num; i++)
                    {
                        cluster_center[i + k * dim_num] = cluster_center[i + k * dim_num]
                                                          + cluster_factor[k] * point[i + j * dim_num];
                    }
                }

                for (k = 0; k < cluster_num; k++)
                {
                    cluster_weight[k] = cluster_weight[k] + cluster_factor[k];
                }
            }

            //
            //  Now normalize,  so that each cluster center is now a convex 
            //  combination of the points.
            //
            for (k = 0; k < cluster_num; k++)
            {
                for (i = 0; i < dim_num; i++)
                {
                    cluster_center[i + k * dim_num] = cluster_center[i + k * dim_num]
                                                      / cluster_weight[k];
                }
            }

            return cluster_center;
        }

        public static double[] cluster_initialize_5(int dim_num, int point_num, int cluster_num,
                double[] point, ref int seed)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //   CLUSTER_INITIALIZE_5 initializes the cluster centers to random values.
            //
            //  Discussion:
            //
            //    In this case, each cluster center is a random convex combination 
            //    of the data points.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    08 October 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the number of spatial dimensions.
            //
            //    Input, int POINT_NUM, the number of points.
            //
            //    Input, int CLUSTER_NUM, the number of clusters.
            //
            //    Input, double POINT[DIM_NUM*POINT_NUM], the coordinates 
            //    of the points.
            //
            //    Input/output, int *SEED, a seed for the random 
            //    number generator.
            //
            //    Output, double CLUSTER_CENTER[DIM_NUM*CLUSTER_NUM],
            //    the coordinates of the cluster centers.
            //
        {
            double[] cluster_center;
            double column_sum;
            double[] factor;
            int j;
            int k;
            //
            //  Get a PxC block of random factors.
            //
            factor = UniformRNG.r8mat_uniform_01_new(point_num, cluster_num, ref seed);
            //
            //  Make each column of factors have unit sum.
            //
            for (k = 0; k < cluster_num; k++)
            {
                column_sum = 0.0;
                for (j = 0; j < point_num; j++)
                {
                    column_sum = column_sum + factor[j + k * point_num];
                }

                for (j = 0; j < point_num; j++)
                {
                    factor[j + k * point_num] = factor[j + k * point_num] / column_sum;
                }
            }

            //
            //  Set centers = points * factors.
            //
            cluster_center = typeMethods.r8mat_mm_new(dim_num, point_num, cluster_num, point,
                factor);

            return cluster_center;
        }

        public static void cluster_print_summary(int point_num, int cluster_num,
                int[] cluster_population, double[] cluster_energy, double[] cluster_variance)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //   CLUSTER_PRINT_SUMMARY prints a summary of data about a clustering.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    08 October 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int POINT_NUM, the number of points.
            //
            //    Input, int CLUSTER_NUM, the number of clusters.
            //
            //    Input, int CLUSTER_POPULATION[CLUSTER_NUM], the number of
            //    points assigned to each cluster.
            //
            //    Input, double CLUSTER_ENERGY[CLUSTER_NUM], the energy of 
            //    the clusters.
            //
            //    Input, double CLUSTER_VARIANCE[CLUSTER_NUM], the variance of 
            //    the clusters.
            //
        {
            double ce;
            int cep;
            double ce_total;
            int cp;
            int cpp;
            double cv;
            int k;

            ce_total = typeMethods.r8vec_sum(cluster_num, cluster_energy);

            Console.WriteLine("");
            Console.WriteLine("  Clustering statistics:");
            Console.WriteLine("");
            Console.WriteLine("    Number of clusters is " + cluster_num + "");
            Console.WriteLine("    Number of points is   " + point_num + "");
            Console.WriteLine("    Total energy is       " + ce_total + "");
            Console.WriteLine("");
            Console.WriteLine("    Cluster   Population        Energy          Variance");
            Console.WriteLine("    -------  -----------  -----------------  --------------");
            Console.WriteLine("                  #    %     value        %");
            Console.WriteLine("");

            for (k = 0; k < cluster_num; k++)
            {
                cp = cluster_population[k];
                cpp = (int) ((double) (100 * cp) / (double) (point_num));
                ce = cluster_energy[k];
                cep = (int) ((ce * 100.0) / ce_total);
                cv = cluster_variance[k];
                Console.WriteLine("  " + k.ToString().PadLeft(7)
                                       + "  " + cp.ToString().PadLeft(8)
                                       + "  " + cpp.ToString().PadLeft(3)
                                       + "  " + ce.ToString().PadLeft(12)
                                       + "  " + cep.ToString().PadLeft(3)
                                       + "  " + cv.ToString().PadLeft(12) + "");
            }

            cp = typeMethods.i4vec_sum(cluster_num, cluster_population);
            cpp = 100;
            ce = typeMethods.r8vec_sum(cluster_num, cluster_energy);
            cep = 100;
            cv = typeMethods.r8vec_i4vec_dot_product(cluster_num, cluster_variance,
                cluster_population) / cp;

            Console.WriteLine("");
            Console.WriteLine("  " + "  Total"
                                   + "  " + cp.ToString().PadLeft(8)
                                   + "  " + cpp.ToString().PadLeft(3)
                                   + "  " + ce.ToString().PadLeft(12)
                                   + "  " + cep.ToString().PadLeft(3)
                                   + "  " + cv.ToString().PadLeft(12) + "");

            return;
        }

        public static double[] cluster_variance_compute(int dim_num, int point_num, int cluster_num,
                double[] point, int[] cluster, double[] cluster_center)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //   CLUSTER_VARIANCE_COMPUTE computes the variance of the clusters.
            //
            //  Discussion:
            //
            //    The cluster variance (from the cluster center) is the average of the 
            //    sum of the squares of the distances of each point in the cluster to the 
            //    cluster center.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    06 October 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the number of spatial dimensions.
            //
            //    Input, int POINT_NUM, the number of data points.
            //
            //    Input, int CLUSTER_NUM, the number of clusters.
            //
            //    Input, double POINT[DIM_NUM*POINT_NUM], the data points.
            //
            //    Input, int CLUSTER[POINT_NUM], the cluster to which each
            //    data point belongs.
            //
            //    Input, double CLUSTER_CENTER[DIM_NUM*CLUSTER_NUM], the 
            //    centers associated with the minimal energy clustering.
            //
            //    Output, double CLUSTER_VARIANCE_COMPUTE[CLUSTER_NUM], the variance
            //    associated with each cluster.
            //
        {
            int[] cluster_population;
            double[] cluster_variance;
            int i;
            int j;
            int k;
            double point_variance;

            cluster_population = typeMethods.i4vec_zero_new(cluster_num);
            cluster_variance = typeMethods.r8vec_zero_new(cluster_num);

            for (j = 0; j < point_num; j++)
            {
                k = cluster[j];

                point_variance = 0.0;
                for (i = 0; i < dim_num; i++)
                {
                    point_variance = point_variance +
                                     Math.Pow(point[i + j * dim_num] - cluster_center[i + k * dim_num], 2);
                }

                cluster_variance[k] = cluster_variance[k] + point_variance;
                cluster_population[k] = cluster_population[k] + 1;
            }

            for (k = 0; k < cluster_num; k++)
            {
                if (0 < cluster_population[k])
                {
                    cluster_variance[k] = cluster_variance[k] / cluster_population[k];
                }
            }

            return cluster_variance;
        }

        public static int find_closest(int m, int n, double[] x, double[] generator)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FIND_CLOSEST finds the Voronoi cell generator closest to a point X.
            //
            //  Discussion:
            //
            //    This routine finds the closest Voronoi cell generator by checking every
            //    one.  For problems with many cells, this process can take the bulk
            //    of the CPU time.  Other approaches, which group the cell generators into
            //    bins, can run faster by a large factor.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    24 September 2006
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the spatial dimension.
            //
            //    Input, int N, the number of cell generators.
            //
            //    Input, double X[M], the point to be checked.
            //
            //    Input, double GENERATOR[M*N], the cell generators.
            //
            //    Output, int FIND_CLOSEST, the index of the nearest cell generators.
            //
        {
            double dist_min;
            double dist;
            int i;
            int j;
            int nearest;

            nearest = 0;
            dist_min = 0.0;

            for (j = 0; j < n; j++)
            {
                dist = 0.0;
                for (i = 0; i < m; i++)
                {
                    dist = dist + Math.Pow(x[i] - generator[i + j * m], 2);
                }

                if (j == 0 || dist < dist_min)
                {
                    dist_min = dist;
                    nearest = j;
                }

            }

            return nearest;
        }

    }
}