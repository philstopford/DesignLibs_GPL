using System;
using Burkardt.Types;

namespace Burkardt.Means;

public static class KMeans
{
    public static void kmeans_01(int dim_num, int point_num, int cluster_num, int it_max,
            ref int it_num, double[] point, ref int[] cluster, ref double[] cluster_center,
            ref int[] cluster_population, ref double[] cluster_energy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //   KMEANS_01 applies the K-Means algorithm.
        //
        //  Discussion:
        //
        //    Given a matrix of POINT_NUM observations on DIM_NUM variables, the
        //    observations are to be allocated to CLUSTER_NUM clusters in such 
        //    a way that the within-cluster sum of squares is minimized.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 October 2011
        //
        //  Author:
        //
        //    Original FORTRAN77 version by David Sparks.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    David Sparks,
        //    Algorithm AS 58: 
        //    Euclidean Cluster Analysis,
        //    Applied Statistics,
        //    Volume 22, Number 1, 1973, pages 126-130.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the number of spatial dimensions.
        //
        //    Input, int POINT_NUM, the number of points.
        //
        //    Input, int CLUSTER_NUM, the number of clusters.
        //
        //    Input, int IT_MAX, the maximum number of iterations.
        //
        //    Output, int IT_NUM, the number of iterations taken.
        //
        //    Input, double POINT[DIM_NUM*POINT_NUM], the points.
        //
        //    Output, int CLUSTER[POINT_NUM], indicates which cluster
        //    each point belongs to.
        //
        //    Input/output, double CLUSTER_CENTER[DIM_NUM*CLUSTER_NUM],
        //    the cluster centers.
        //
        //    Output, int CLUSTER_POPULATION[CLUSTER_NUM], the number 
        //    of points in each cluster.
        //
        //    Output, double CLUSTER_ENERGY[CLUSTER_NUM], the 
        //    cluster energies.
        //
    {
        double dc;
        double de;
        double[] f;
        int i;
        int il;
        int ir;
        int j;
        int j2;
        int k;
        double point_energy;
        double point_energy_min;
        int swap;

        it_num = 0;
        switch (cluster_num)
        {
            //
            //  Idiot checks.
            //
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("KMEANS_01 - Fatal error!");
                Console.WriteLine("  CLUSTER_NUM < 1.");
                return;
        }

        switch (dim_num)
        {
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("KMEANS_01 - Fatal error!");
                Console.WriteLine("  DIM_NUM < 1.");
                return;
        }

        switch (point_num)
        {
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("KMEANS_01 - Fatal error!");
                Console.WriteLine("  POINT_NUM < 1.");
                return;
        }

        //
        //  For each observation, calculate the distance from each cluster
        //  center, and assign to the nearest.
        //
        for (j = 0; j < point_num; j++)
        {
            point_energy_min = typeMethods.r8_huge();
            cluster[j] = -1;

            for (k = 0; k < cluster_num; k++)
            {
                point_energy = 0.0;
                for (i = 0; i < dim_num; i++)
                {
                    point_energy += Math.Pow(point[i + j * dim_num] - cluster_center[i + k * dim_num], 2);
                }

                if (point_energy < point_energy_min)
                {
                    point_energy_min = point_energy;
                    cluster[j] = k;
                }
            }
        }

        //
        //  Determine the cluster population counts.
        //
        typeMethods.i4vec_zero(cluster_num, ref cluster_population);

        for (j = 0; j < point_num; j++)
        {
            k = cluster[j];
            cluster_population[k] += 1;
        }

        //
        //  Calculate the mean and sum of squares for each cluster.
        //
        typeMethods.r8vec_zero(dim_num * cluster_num, ref cluster_center);

        for (j = 0; j < point_num; j++)
        {
            k = cluster[j];
            for (i = 0; i < dim_num; i++)
            {
                cluster_center[i + k * dim_num] += point[i + j * dim_num];
            }
        }

        for (k = 0; k < cluster_num; k++)
        {
            switch (cluster_population[k])
            {
                case > 0:
                {
                    for (i = 0; i < dim_num; i++)
                    {
                        cluster_center[i + k * dim_num] /= cluster_population[k];
                    }

                    break;
                }
            }
        }

        //
        //  Set the point energies.
        //
        f = typeMethods.r8vec_zero_new(point_num);

        for (j = 0; j < point_num; j++)
        {
            k = cluster[j];
            for (i = 0; i < dim_num; i++)
            {
                f[j] += Math.Pow(point[i + j * dim_num] - cluster_center[i + k * dim_num], 2);
            }
        }

        //
        //  Set the cluster energies.
        //
        typeMethods.r8vec_zero(cluster_num, ref cluster_energy);

        for (j = 0; j < point_num; j++)
        {
            k = cluster[j];
            cluster_energy[k] += f[j];
        }

        //
        //  Adjust the point energies by a weight factor.
        //
        for (j = 0; j < point_num; j++)
        {
            k = cluster[j];
            f[j] = cluster_population[k] switch
            {
                > 1 => f[j] * cluster_population[k] / (cluster_population[k] - 1),
                _ => f[j]
            };
        }

        //
        //  Examine each observation in turn to see if it should be
        //  reassigned to a different cluster.
        //
        it_num = 0;

        while (it_num < it_max)
        {
            it_num += 1;

            swap = 0;

            for (j = 0; j < point_num; j++)
            {
                il = cluster[j];
                ir = il;

                switch (cluster_population[il])
                {
                    case <= 1:
                        continue;
                }

                dc = f[j];

                for (k = 0; k < cluster_num; k++)
                {
                    if (k != il)
                    {
                        de = 0.0;
                        for (i = 0; i < dim_num; i++)
                        {
                            de += Math.Pow(point[i + j * dim_num] - cluster_center[i + k * dim_num], 2);
                        }

                        de = de * cluster_population[k]
                             / (cluster_population[k] + 1);

                        if (de < dc)
                        {
                            dc = de;
                            ir = k;
                        }
                    }
                }

                //
                //  If the lowest value was obtained by staying in the current cluster,
                //  then cycle.
                //
                if (ir == il)
                {
                    continue;
                }

                //
                //  Reassign the point from cluster IL to cluster IR.
                //
                for (i = 0; i < dim_num; i++)
                {
                    cluster_center[i + il * dim_num] = (cluster_center[i + il * dim_num]
                                                        * cluster_population[il] -
                                                        point[i + j * dim_num])
                                                       / (cluster_population[il] - 1);

                    cluster_center[i + ir * dim_num] = (cluster_center[i + ir * dim_num]
                                                        * cluster_population[ir] +
                                                        point[i + j * dim_num])
                                                       / (cluster_population[ir] + 1);
                }

                cluster_energy[il] -= f[j];
                cluster_energy[ir] += dc;
                cluster_population[ir] += 1;
                cluster_population[il] -= 1;

                cluster[j] = ir;
                //
                //  Adjust the value of F for points in clusters IL and IR.
                //  
                for (j2 = 0; j2 < point_num; j2++)
                {
                    k = cluster[j2];

                    if (k == il || k == ir)
                    {
                        f[j2] = 0.0;
                        for (i = 0; i < dim_num; i++)
                        {
                            f[j2] += Math.Pow(point[i + j2 * dim_num] - cluster_center[i + k * dim_num], 2);
                        }

                        f[j2] = cluster_population[k] switch
                        {
                            > 1 => f[j2] * cluster_population[k] / (cluster_population[k] - 1),
                            _ => f[j2]
                        };
                    }
                }

                swap += 1;
            }

            //
            //  Exit if no reassignments were made during this iteration.
            //
            if (swap == 0)
            {
                break;
            }
        }

        //
        //  Compute the cluster energies.
        //
        typeMethods.r8vec_zero(cluster_num, ref cluster_energy);

        for (j = 0; j < point_num; j++)
        {
            k = cluster[j];

            point_energy = 0.0;
            for (i = 0; i < dim_num; i++)
            {
                point_energy += Math.Pow(point[i + j * dim_num] - cluster_center[i + k * dim_num], 2);
            }

            cluster_energy[k] += point_energy;
        }
    }

    public static void kmeans_02(int dim_num, int point_num, int cluster_num, int it_max,
            ref int it_num, double[] point, ref int[] cluster, ref double[] cluster_center,
            ref int[] cluster_population, ref double[] cluster_energy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //   KMEANS_02 applies the K-Means algorithm.
        //
        //  Discussion:
        //
        //    The routine attempts to divide POINT_NUM points in 
        //    DIM_NUM-dimensional space into CLUSTER_NUM clusters so that the within 
        //    cluster sum of squares is minimized.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 October 2011
        //
        //  Author:
        //
        //    Original FORTRAN77 by John Hartigan, Manchek Wong.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    John Hartigan, Manchek Wong,
        //    Algorithm AS 136:
        //    A K-Means Clustering Algorithm,
        //    Applied Statistics,
        //    Volume 28, Number 1, 1979, pages 100-108.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the number of spatial dimensions.
        //
        //    Input, int POINT_NUM, the number of points.
        //
        //    Input, int CLUSTER_NUM, the number of clusters.
        //
        //    Input, int IT_MAX, the maximum number of iterations.
        //
        //    Output, int &IT_NUM, the number of iterations taken.
        //
        //    Input, double POINT[DIM_NUM*POINT_NUM], the coordinates 
        //    of the points.
        //
        //    Output, int CLUSTER[POINT_NUM], the cluster each 
        //    point belongs to.
        //
        //    Input/output, double CLUSTER_CENTER[DIM_NUM*CLUSTER_NUM],
        //    the cluster centers.
        //
        //    Output, int CLUSTER_POPULATION[CLUSTER_NUM], the number 
        //    of points in each cluster.
        //
        //    Output, double CLUSTER_ENERGY[CLUSTER_NUM], the 
        //    within-cluster sum of squares.
        //
    {
        double[] an1;
        double[] an2;
        int[] cluster2;
        double[] d;
        double db;
        double[] dt = new double[2];
        int i;
        int ifault;
        int il;
        int indx;
        int[] itran;
        int j;
        int k;
        int[] live;
        int[] ncp;
        double point_energy;
        double temp;

        it_num = 0;

        switch (cluster_num)
        {
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("KMEANS_02 - Fatal error!");
                Console.WriteLine("  CLUSTER_NUM < 1.");
                return;
        }

        if (point_num <= cluster_num)
        {
            Console.WriteLine("");
            Console.WriteLine("KMEANS_02 - Fatal error!");
            Console.WriteLine("  POINT_NUM <= CLUSTER_NUM.");
            return;
        }

        switch (dim_num)
        {
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("KMEANS_02 - Fatal error!");
                Console.WriteLine("  DIM_NUM < 1.");
                return;
        }

        switch (point_num)
        {
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("KMEANS_02 - Fatal error!");
                Console.WriteLine("  POINT_NUM < 1.");
                return;
        }

        //
        //  For each point I, find its two closest centers, CLUSTER(I) and CLUSTER2(I).
        //  Assign it to CLUSTER(I).
        //
        cluster2 = new int[point_num];

        for (j = 0; j < point_num; j++)
        {
            cluster[j] = 0;
            cluster2[j] = 1;

            for (il = 0; il < 2; il++)
            {
                dt[il] = 0.0;
                for (i = 0; i < dim_num; i++)
                {
                    dt[il] += Math.Pow(point[i + j * dim_num] - cluster_center[i + il * dim_num], 2);
                }
            }

            if (dt[1] < dt[0])
            {
                cluster[j] = 1;
                cluster2[j] = 0;
                temp = dt[0];
                dt[0] = dt[1];
                dt[1] = temp;
            }

            for (k = 2; k < cluster_num; k++)
            {
                db = 0.0;
                for (i = 0; i < dim_num; i++)
                {
                    db += Math.Pow(point[i + j * dim_num] - cluster_center[i + k * dim_num], 2);
                }

                if (db < dt[0])
                {
                    dt[1] = dt[0];
                    cluster2[j] = cluster[j];
                    dt[0] = db;
                    cluster[j] = k;
                }
                else if (db < dt[1])
                {
                    dt[1] = db;
                    cluster2[j] = k;
                }
            }
        }

        //
        //  Update cluster centers to be the average of points contained
        //  within them.
        //
        typeMethods.i4vec_zero(cluster_num, ref cluster_population);
        typeMethods.r8vec_zero(dim_num * cluster_num, ref cluster_center);

        for (j = 0; j < point_num; j++)
        {
            k = cluster[j];
            cluster_population[k] += 1;
            for (i = 0; i < dim_num; i++)
            {
                cluster_center[i + k * dim_num] += point[i + j * dim_num];
            }
        }

        //
        //  Check to see if there is any empty cluster.
        //
        an1 = new double[cluster_num];
        an2 = new double[cluster_num];
        itran = new int[cluster_num];
        ncp = new int[cluster_num];

        for (k = 0; k < cluster_num; k++)
        {
            switch (cluster_population[k])
            {
                case 0:
                    Console.WriteLine("");
                    Console.WriteLine("KMEANS_02 - Fatal error!");
                    Console.WriteLine("  Cluster " + k + " has zero population.");
                    return;
            }

            for (i = 0; i < dim_num; i++)
            {
                cluster_center[i + k * dim_num] /= cluster_population[k];
            }

            //
            //  Initialize AN1, AN2, ITRAN and NCP
            //  AN1(K) = CLUSTER_POPULATION(K) / (CLUSTER_POPULATION(K) - 1)
            //  AN2(K) = CLUSTER_POPULATION(K) / (CLUSTER_POPULATION(K) + 1)
            //  ITRAN(K) = 1 if cluster K is updated in the quick-transfer stage,
            //           = 0 otherwise
            //  In the optimal-transfer stage, NCP(K) stores the step at which
            //  cluster K is last updated.
            //  In the quick-transfer stage, NCP(K) stores the step at which
            //  cluster K is last updated plus POINT_NUM.
            //
            an2[k] = cluster_population[k]
                     / (double) (cluster_population[k] + 1);

            an1[k] = cluster_population[k] switch
            {
                > 1 => cluster_population[k] / (double) (cluster_population[k] - 1),
                _ => typeMethods.r8_huge()
            };

            itran[k] = 1;
            ncp[k] = -1;
        }

        indx = 0;
        ifault = 2;
        it_num = 0;

        d = new double[point_num];
        live = new int[cluster_num];

        while (it_num < it_max)
        {
            it_num += 1;
            //
            //  In this stage, there is only one pass through the data.   Each
            //  point is re-allocated, if necessary, to the cluster that will
            //  induce the maximum reduction in within-cluster sum of squares.
            //
            kmeans_02_optra(dim_num, point_num, cluster_num, point,
                ref cluster_center, ref cluster, ref cluster2, ref cluster_population, ref an1, ref an2,
                ref ncp, ref d, ref itran, ref live, ref indx);
            //
            //  Stop if no transfer took place in the last POINT_NUM optimal transfer steps.
            //
            if (indx == point_num)
            {
                ifault = 0;
                break;
            }

            //
            //  Each point is tested in turn to see if it should be re-allocated
            //  to the cluster to which it is most likely to be transferred,
            //  CLUSTER2(I), from its present cluster, CLUSTER(I).   Loop through the
            //  data until no further change is to take place.
            //
            kmeans_02_qtran(dim_num, point_num, cluster_num, point,
                ref cluster_center, ref cluster, ref cluster2, ref cluster_population, ref an1, ref an2,
                ref ncp, ref d, ref itran, ref indx);
            //
            //  If there are only two clusters, there is no need to re-enter the
            //  optimal transfer stage.
            //
            if (cluster_num == 2)
            {
                ifault = 0;
                break;
            }

            //
            //  NCP has to be set to 0 before entering OPTRA.
            //
            typeMethods.i4vec_zero(cluster_num, ref ncp);
        }

        switch (ifault)
        {
            case 2:
                Console.WriteLine("");
                Console.WriteLine("KMEANS_02 - Warning!");
                Console.WriteLine("  Maximum number of iterations reached");
                Console.WriteLine("  without convergence.");
                break;
        }

        //
        //  Compute the within-cluster sum of squares for each cluster.
        //
        typeMethods.r8vec_zero(dim_num * cluster_num, ref cluster_center);

        for (j = 0; j < point_num; j++)
        {
            k = cluster[j];
            for (i = 0; i < dim_num; i++)
            {
                cluster_center[i + k * dim_num] += point[i + j * dim_num];
            }
        }

        for (k = 0; k < cluster_num; k++)
        {
            for (i = 0; i < dim_num; i++)
            {
                cluster_center[i + k * dim_num] /= cluster_population[k];
            }
        }

        //
        //  Compute the cluster energies.
        //
        typeMethods.r8vec_zero(cluster_num, ref cluster_energy);

        for (j = 0; j < point_num; j++)
        {
            k = cluster[j];

            point_energy = 0.0;
            for (i = 0; i < dim_num; i++)
            {
                point_energy += Math.Pow(point[i + j * dim_num] - cluster_center[i + k * dim_num], 2);
            }

            cluster_energy[k] += point_energy;
        }
    }

    public static void kmeans_02_optra(int dim_num, int point_num, int cluster_num,
            double[] point, ref double[] cluster_center, ref int[] cluster, ref int[] cluster2,
            ref int[] cluster_population, ref double[] an1, ref double[] an2, ref int[] ncp,
            ref double[] d, ref int[] itran, ref int[] live, ref int indx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //   KMEANS_02_OPTRA carries out the optimal transfer stage.
        //
        //  Discussion:
        //
        //    Each point is re-allocated, if necessary, to the cluster that
        //    will induce a maximum reduction in the within-cluster sum of
        //    squares.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 October 2011
        //
        //  Author:
        //
        //    Original FORTRAN77 by John Hartigan, Manchek Wong.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    John Hartigan, Manchek Wong,
        //    Algorithm AS 136:
        //    A K-Means Clustering Algorithm,
        //    Applied Statistics,
        //    Volume 28, Number 1, 1979, pages 100-108.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the number of spatial dimensions.
        //
        //    Input, int POINT_NUM, the number of points.
        //
        //    Input, int CLUSTER_NUM, the number of clusters.
        //
        //    Input, double POINT[DIM_NUM*POINT_NUM], the coordinates of 
        //    the points.
        //
        //    Input/output, double CLUSTER_CENTER[DIM_NUM*CLUSTER_NUM],
        //    the cluster centers.
        //
        //    Input/output, int CLUSTER[POINT_NUM], the cluster 
        //    each point belongs to.
        //
        //    Input/output, int CLUSTER2[POINT_NUM], the cluster 
        //    to which each point is most likely to be transferred to.
        //
        //    Input/output, int CLUSTER_POPULATION[CLUSTER_NUM], 
        //    the number of points in each cluster.
        //
        //    Input/output, double AN1[CLUSTER_NUM], 
        //    CLUSTER_POPULATION(L) / (CLUSTER_POPULATION(L) - 1)
        //
        //    Input/output, double AN2[CLUSTER_NUM], 
        //    CLUSTER_POPULATION(L) / (CLUSTER_POPULATION(L) + 1)
        //
        //    Input/output, int NCP[CLUSTER_NUM], ?
        //
        //    Input/output, double D[POINT_NUM], ?
        //
        //    Input/output, int ITRAN[CLUSTER_NUM], 
        //    1 if cluster L is updated in the quick-transfer stage,
        //    0 otherwise.  Reset to zero on output.
        //
        //    Input/output, int LIVE[CLUSTER_NUM], ?
        //
        //    Input/output, int INDX, ?
        //
    {
        double al1;
        double al2;
        double alt;
        double alw;
        double dc;
        int i;
        int j;
        int k;
        int l;
        int l1;
        int l2;
        int ll;
        double r2;
        double rr;
        //
        //  If cluster L is updated in the last quick-transfer stage, it
        //  belongs to the live set throughout this stage.   Otherwise, at
        //  each step, it is not in the live set if it has not been updated
        //  in the last POINT_NUM optimal transfer steps.
        //
        for (k = 0; k < cluster_num; k++)
        {
            live[k] = itran[k] switch
            {
                1 => point_num + 1,
                _ => live[k]
            };
        }

        for (j = 0; j < point_num; j++)
        {
            indx += 1;
            l1 = cluster[j];
            l2 = cluster2[j];
            ll = l2;
            switch (cluster_population[l1])
            {
                //
                //  If point J is the only member of cluster L1, no transfer.
                //
                case > 1:
                {
                    //
                    //  If L1 has been updated in this stage, re-compute D(I).
                    //
                    if (ncp[l1] != 0)
                    {
                        d[j] = 0.0;
                        for (i = 0; i < dim_num; i++)
                        {
                            d[j] += Math.Pow(point[i + j * dim_num] - cluster_center[i + l1 * dim_num], 2);
                        }

                        d[j] = an1[l1] * d[j];
                    }

                    //
                    //  Find the cluster with minimum R2.
                    //
                    r2 = 0.0;
                    for (i = 0; i < dim_num; i++)
                    {
                        r2 += Math.Pow(point[i + j * dim_num] - cluster_center[i + l2 * dim_num], 2);
                    }

                    r2 = an2[l2] * r2;

                    for (l = 0; l < cluster_num; l++)
                    {
                        //
                        //  If LIVE(L1) <= J, then L1 is not in the live set.   If this is
                        //  true, we only need to consider clusters that are in the live set
                        //  for possible transfer of point J.   
                        //
                        //  Otherwise, we need to consider all possible clusters.
                        //
                        if ((j < live[l1] || j < live[l]) && l != l1 && l != ll)
                        {
                            rr = r2 / an2[l];

                            dc = 0.0;
                            for (i = 0; i < dim_num; i++)
                            {
                                dc += Math.Pow(point[i + j * dim_num] - cluster_center[i + l * dim_num], 2);
                            }

                            if (dc < rr)
                            {
                                r2 = dc * an2[l];
                                l2 = l;
                            }
                        }
                    }

                    //
                    //  If no transfer is necessary, L2 is the new CLUSTER2(J).
                    // 
                    if (d[j] <= r2)
                    {
                        cluster2[j] = l2;
                    }
                    //
                    //  Update cluster centers, LIVE, NCP, AN1 and AN2 for clusters L1 and
                    //  L2, and update CLUSTER(J) and CLUSTER2(J).
                    //
                    else
                    {
                        indx = 0;
                        live[l1] = point_num + j;
                        live[l2] = point_num + j;
                        ncp[l1] = j;
                        ncp[l2] = j;
                        al1 = cluster_population[l1];
                        alw = al1 - 1.0;
                        al2 = cluster_population[l2];
                        alt = al2 + 1.0;

                        for (i = 0; i < dim_num; i++)
                        {
                            cluster_center[i + l1 * dim_num] = (cluster_center[i + l1 * dim_num] * al1
                                                                - point[i + j * dim_num]) / alw;

                            cluster_center[i + l2 * dim_num] = (cluster_center[i + l2 * dim_num] * al2
                                                                + point[i + j * dim_num]) / alt;
                        }

                        cluster_population[l1] -= 1;
                        cluster_population[l2] += 1;
                        an2[l1] = alw / al1;

                        an1[l1] = alw switch
                        {
                            > 1.0 => alw / (alw - 1.0),
                            _ => typeMethods.r8_huge()
                        };

                        an1[l2] = alt / al2;
                        an2[l2] = alt / (alt + 1.0);
                        cluster[j] = l2;
                        cluster2[j] = l1;
                    }

                    break;
                }
            }

            if (indx == point_num)
            {
                return;
            }

        }

        //
        //  ITRAN(L) = 0 before entering QTRAN.
        //
        typeMethods.i4vec_zero(cluster_num, ref itran);
        //
        //  LIVE(L) has to be decreased by POINT_NUM before re-entering OPTRA.
        //
        for (k = 0; k < cluster_num; k++)
        {
            live[k] -= point_num;
        }
    }

    public static void kmeans_02_qtran(int dim_num, int point_num, int cluster_num,
            double[] point, ref double[] cluster_center, ref int[] cluster, ref int[] cluster2,
            ref int[] cluster_population, ref double[] an1, ref double[] an2, ref int[] ncp, ref double[] d,
            ref int[] itran, ref int indx)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //   KMEANS_02_QTRAN carries out the quick transfer stage.
        //
        //  Discussion:
        //
        //    For each point I, CLUSTER(I) and CLUSTER2(I) are switched, if necessary, 
        //    to reduce within-cluster sum of squares.  The cluster centers are
        //    updated after each step.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 October 2011
        //
        //  Author:
        //
        //    Original FORTRAN77 by John Hartigan, Manchek Wong.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    John Hartigan, Manchek Wong,
        //    Algorithm AS 136:
        //    A K-Means Clustering Algorithm,
        //    Applied Statistics,
        //    Volume 28, Number 1, 1979, pages 100-108.
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
        //    Input/output, double CLUSTER_CENTER[DIM_NUM*CLUSTER_NUM],
        //    the cluster centers.
        //
        //    Input/output, int CLUSTER[POINT_NUM], the cluster 
        //    each point belongs to.
        //
        //    Input/output, int CLUSTER2[POINT_NUM], the cluster to 
        //    which each point is most likely to be transferred to.
        //
        //    Input/output, int CLUSTER_POPULATION[CLUSTER_NUM], 
        //    the number of points in each cluster.
        //
        //    Input/output, double AN1[CLUSTER_NUM], 
        //    CLUSTER_POPULATION(L) / (CLUSTER_POPULATION(L) - 1).
        //
        //    Input/output, double AN2[CLUSTER_NUM], 
        //    CLUSTER_POPULATION(L) / (CLUSTER_POPULATION(L) + 1).
        //
        //    Input/output, int NCP[CLUSTER_NUM], ?
        //
        //    Input/output, double D[POINT_NUM], ?
        //
        //    Input/output, int ITRAN[CLUSTER_NUM], 
        //    1 if cluster L is updated in the quick-transfer stage,
        //    0 otherwise.
        //
        //    Input/output, int INDX, is set to 0 if any 
        //    updating occurs.
        //
    {
        double al1;
        double al2;
        double alt;
        double alw;
        int count;
        double dd;
        int i;
        int j;
        int l1;
        int l2;
        double r2;
        int step;
        //
        //  In the optimal transfer stage, NCP(L) indicates the step at which
        //  cluster L is last updated.   In the quick transfer stage, NCP(L)
        //  is equal to the step at which cluster L is last updated plus POINT_NUM.
        //
        count = 0;
        step = 0;

        for (;;)
        {
            for (j = 0; j < point_num; j++)
            {
                count += 1;
                step += 1;
                l1 = cluster[j];
                l2 = cluster2[j];
                switch (cluster_population[l1])
                {
                    //
                    //  If point I is the only member of cluster L1, no transfer.
                    //
                    case > 1:
                    {
                        //
                        //  If NCP(L1) < STEP, no need to re-compute distance from point I to
                        //  cluster L1.   Note that if cluster L1 is last updated exactly POINT_NUM
                        //  steps ago, we still need to compute the distance from point I to
                        //  cluster L1.
                        //
                        if (step <= ncp[l1])
                        {
                            d[j] = 0.0;
                            for (i = 0; i < dim_num; i++)
                            {
                                d[j] += Math.Pow(point[i + j * dim_num] - cluster_center[i + l1 * dim_num], 2);
                            }

                            d[j] = an1[l1] * d[j];
                        }

                        //
                        //  If STEP >= both NCP(L1) and NCP(L2) there will be no transfer of
                        //  point I at this step.
                        //
                        if (step < ncp[l1] || step < ncp[l2])
                        {
                            r2 = d[j] / an2[l2];

                            dd = 0.0;
                            for (i = 0; i < dim_num; i++)
                            {
                                dd += Math.Pow(point[i + j * dim_num] - cluster_center[i + l2 * dim_num], 2);
                            }

                            //
                            //  Update cluster centers, NCP, CLUSTER_POPULATION, ITRAN, AN1 and AN2 
                            //  for clusters L1 and L2.   Also update CLUSTER(J) and CLUSTER2(J).   
                            //
                            //  Note that if any updating occurs in this stage, INDX is set back to 0.
                            //
                            if (dd < r2)
                            {
                                count = 0;
                                indx = 0;
                                itran[l1] = 1;
                                itran[l2] = 1;
                                ncp[l1] = step + point_num;
                                ncp[l2] = step + point_num;
                                al1 = cluster_population[l1];
                                alw = al1 - 1.0;
                                al2 = cluster_population[l2];
                                alt = al2 + 1.0;

                                for (i = 0; i < dim_num; i++)
                                {
                                    cluster_center[i + l1 * dim_num] = (cluster_center[i + l1 * dim_num] * al1
                                                                        - point[i + j * dim_num]) / alw;

                                    cluster_center[i + l2 * dim_num] = (cluster_center[i + l2 * dim_num] * al2
                                                                        + point[i + j * dim_num]) / alt;
                                }

                                cluster_population[l1] -= 1;
                                cluster_population[l2] += 1;
                                an2[l1] = alw / al1;

                                an1[l1] = alw switch
                                {
                                    > 1.0 => alw / (alw - 1.0),
                                    _ => typeMethods.r8_huge()
                                };

                                an1[l2] = alt / al2;
                                an2[l2] = alt / (alt + 1.0);
                                cluster[j] = l2;
                                cluster2[j] = l1;
                            }
                        }

                        break;
                    }
                }

                //
                //  If no re-allocation took place in the last POINT_NUM steps, return.
                //
                if (count == point_num)
                {
                    return;
                }
            }
        }

        return;
    }

    public static void kmeans_03(int dim_num, int point_num, int cluster_num, int it_max,
            ref int it_num, double[] point, ref int[] cluster, ref double[] cluster_center,
            ref int[] cluster_population, ref double[] cluster_energy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //   KMEANS_03 applies the K-Means algorithm.
        //
        //  Discussion:
        //
        //    It is possible for a straightforward K-Means algorithm to
        //    halt at a non-optimal partition of the points.  This routine
        //    tries to improve the input partition if possible.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Wendy Martinez, Angel Martinez,
        //    Computational Statistics Handbook with MATLAB,
        //    pages 373-376,
        //    Chapman and Hall / CRC, 2002.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the number of spatial dimensions.
        //
        //    Input, int POINT_NUM, the number of data points.
        //
        //    Input, int CLUSTER_NUM, the number of clusters.
        //
        //    Input, int IT_MAX, the maximum number of iterations.
        //
        //    Output, int IT_NUM, the number of iterations taken.
        //
        //    Input, double POINT[DIM_NUM*POINT_NUM], the data points.
        //
        //    Output, int CLUSTER[POINT_NUM], the cluster to which
        //    each point belongs.
        //
        //    Input/output, double CLUSTER_CENTER[DIM_NUM*CLUSTER_NUM], the 
        //    centers associated with the clustering.  On output, these may 
        //    have been altered.
        //
        //    Output, int CLUSTER_POPULATION[CLUSTER_NUM], the number
        //    of points in each cluster.
        //
        //    Output, double CLUSTER_ENERGY[CLUSTER_NUM], the energy of 
        //    the clusters.
        //
    {
        int ci;
        int cj;
        double[] distsq;
        int i;
        int j;
        int k;
        double point_energy;
        double point_energy_min;
        int swap;
        switch (cluster_num)
        {
            //
            //  Check the input.
            //
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("KMEANS_03 - Fatal error!");
                Console.WriteLine("  CLUSTER_NUM < 1.");
                return;
        }

        switch (dim_num)
        {
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("KMEANS_03 - Fatal error!");
                Console.WriteLine("  DIM_NUM < 1.");
                return;
        }

        switch (point_num)
        {
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("KMEANS_03 - Fatal error!");
                Console.WriteLine("  POINT_NUM < 1.");
                return;
        }

        switch (it_max)
        {
            case < 0:
                Console.WriteLine("");
                Console.WriteLine("KMEANS_03 - Fatal error!");
                Console.WriteLine("  IT_MAX < 0.");
                return;
        }

        //
        //  Assign each point to the nearest cluster center.
        //
        for (j = 0; j < point_num; j++)
        {
            point_energy_min = typeMethods.r8_huge();
            cluster[j] = -1;

            for (k = 0; k < cluster_num; k++)
            {
                point_energy = 0.0;
                for (i = 0; i < dim_num; i++)
                {
                    point_energy += Math.Pow(point[i + j * dim_num] - cluster_center[i + k * dim_num], 2);
                }

                if (point_energy < point_energy_min)
                {
                    point_energy_min = point_energy;
                    cluster[j] = k;
                }
            }
        }

        //
        //  Determine the cluster populations.
        //
        typeMethods.i4vec_zero(cluster_num, ref cluster_population);

        for (j = 0; j < point_num; j++)
        {
            k = cluster[j];
            cluster_population[k] += 1;
        }

        //
        //  Average the points in each cluster to get a new cluster center.
        //
        typeMethods.r8vec_zero(dim_num * cluster_num, ref cluster_center);

        for (j = 0; j < point_num; j++)
        {
            k = cluster[j];
            for (i = 0; i < dim_num; i++)
            {
                cluster_center[i + k * dim_num] += point[i + j * dim_num];
            }
        }

        for (k = 0; k < cluster_num; k++)
        {
            for (i = 0; i < dim_num; i++)
            {
                cluster_center[i + k * dim_num] /= cluster_population[k];
            }
        }

        //
        //  Carry out the iteration.
        //
        it_num = 0;
        distsq = new double[cluster_num];

        while (it_num < it_max)
        {
            it_num += 1;

            swap = 0;

            for (j = 0; j < point_num; j++)
            {
                ci = cluster[j];

                switch (cluster_population[ci])
                {
                    case <= 1:
                        continue;
                }

                for (cj = 0; cj < cluster_num; cj++)
                {
                    if (cj == ci)
                    {
                        distsq[cj] = 0.0;
                        for (i = 0; i < dim_num; i++)
                        {
                            distsq[cj] += Math.Pow(point[i + j * dim_num] - cluster_center[i + cj * dim_num], 2);
                        }

                        distsq[cj] = distsq[cj] * cluster_population[cj]
                                     / (cluster_population[cj] - 1);
                    }
                    else
                    {
                        switch (cluster_population[cj])
                        {
                            case 0:
                            {
                                for (i = 0; i < dim_num; i++)
                                {
                                    cluster_center[i + cj * dim_num] = point[i + j * dim_num];
                                }

                                distsq[cj] = 0.0;
                                break;
                            }
                            default:
                            {
                                distsq[cj] = 0.0;
                                for (i = 0; i < dim_num; i++)
                                {
                                    distsq[cj] += Math.Pow(point[i + j * dim_num] - cluster_center[i + cj * dim_num], 2);
                                }

                                distsq[cj] = distsq[cj] * cluster_population[cj]
                                             / (cluster_population[cj] + 1);
                                break;
                            }
                        }
                    }
                }

                //
                //  Find the index of the minimum value of DISTSQ.
                //
                k = typeMethods.r8vec_min_index(cluster_num, distsq);
                //
                //  If that is not the cluster to which point I now belongs, move it there.
                //
                if (k == ci)
                {
                    continue;
                }

                cj = k;

                for (i = 0; i < dim_num; i++)
                {
                    cluster_center[i + ci * dim_num] = (cluster_population[ci]
                                                        * cluster_center[i + ci * dim_num] -
                                                        point[i + j * dim_num])
                                                       / (cluster_population[ci] - 1);

                    cluster_center[i + cj * dim_num] = (cluster_population[cj]
                                                        * cluster_center[i + cj * dim_num] +
                                                        point[i + j * dim_num])
                                                       / (cluster_population[cj] + 1);
                }

                cluster_population[ci] -= 1;
                cluster_population[cj] += 1;

                cluster[j] = cj;
                swap += 1;
            }

            //
            //  Exit if no reassignments were made during this iteration.
            //
            if (swap == 0)
            {
                break;
            }
        }

        //
        //  Compute the cluster energies.
        //
        typeMethods.r8vec_zero(cluster_num, ref cluster_energy);

        for (j = 0; j < point_num; j++)
        {
            k = cluster[j];

            point_energy = 0.0;
            for (i = 0; i < dim_num; i++)
            {
                point_energy += Math.Pow(point[i + j * dim_num] - cluster_center[i + k * dim_num], 2);
            }

            cluster_energy[k] += point_energy;
        }
    }

    public static void kmeans_w_01(int dim_num, int point_num, int cluster_num, int it_max,
            ref int it_num, double[] point, double[] weight, ref int[] cluster,
            ref double[] cluster_center, ref int[] cluster_population, ref double[] cluster_energy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //   KMEANS_W_01 applies the weighted K-Means algorithm.
        //
        //  Discussion:
        //
        //    The input data for the weight K-Means problem includes:
        //    * a set of N data points X in M dimensions, 
        //    * a set of N nonnegative weights W,
        //    * a desired number of clusters K.
        //    * an initial set of cluster centers Z,
        //    * an (optional) initial set of cluster assignments.
        //
        //    The goal is to determine K points Z, called cluster centers, and
        //    to assign each point X(I) to some cluster Z(J), so that we minimize
        //    the weighted standard deviation of the distance of each data point
        //    to the center of its cluster.  Writing J = CLUSTER(I) to
        //    indicate the index of the nearest cluster center Z(J) to the 
        //    point X(I), the quantity we are trying to minimize is the sum
        //    of the weighted cluster energies E(J), where:
        //
        //      E(J) = Sum ( 1 <= I <= N ) W(I) * || X(I) - Z(J) ||^2
        //
        //    Here, we assume that we are using the Euclidean norm, so that
        //    
        //      || X(I) - Z(J) ||^2 = Sum ( 1 <= K <= M )
        //         ( X(I)(K) - Z(J)(K) )^2
        //
        //    In this notation, X(I)(K) is the K-th spatial component of the
        //    I-th data point.
        //
        //    Note that this routine should give the same results as KMEANS_01
        //    in any case in which all the entries of the WEIGHT vector are equal.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    David Sparks,
        //    Algorithm AS 58: Euclidean Cluster Analysis,
        //    Applied Statistics,
        //    Volume 22, Number 1, 1973, pages 126-130.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the number of spatial dimensions.
        //
        //    Input, int POINT_NUM, the number of points.
        //
        //    Input, int CLUSTER_NUM, the number of clusters.
        //
        //    Input, int IT_MAX, the maximum number of iterations.
        //
        //    Output, int &IT_NUM, the number of iterations taken.
        //
        //    Input, double POINT[DIM_NUM*POINT_NUM], the points.
        //
        //    Input, double WEIGHT[POINT_NUM], the weights.
        //
        //    Output, int CLUSTER[POINT_NUM], indicates which cluster
        //    each point belongs to.
        //
        //    Input/output, double CLUSTER_CENTER[DIM_NUM*CLUSTER_NUM],
        //    the cluster centers.
        //
        //    Output, int CLUSTER_POPULATION[CLUSTER_NUM], the number 
        //    of points in each cluster.
        //
        //    Output, double CLUSTER_ENERGY[CLUSTER_NUM], the 
        //    cluster energies.
        //
    {
        double[] cluster_weight;
        double dc;
        double de;
        double[] f;
        int i;
        int il;
        int ir;
        int j;
        int k;
        double point_energy;
        double point_energy_min;
        int swap;

        it_num = 0;
        switch (cluster_num)
        {
            //
            //  Idiot checks.
            //
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("KMEANS_W_01 - Fatal error!");
                Console.WriteLine("  CLUSTER_NUM < 1.");
                return;
        }

        switch (dim_num)
        {
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("KMEANS_W_01 - Fatal error!");
                Console.WriteLine("  DIM_NUM < 1.");
                return;
        }

        switch (point_num)
        {
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("KMEANS_W_01 - Fatal error!");
                Console.WriteLine("  POINT_NUM < 1.");
                return;
        }

        switch (it_max)
        {
            case < 0:
                Console.WriteLine("");
                Console.WriteLine("KMEANS_W_01 - Fatal error!");
                Console.WriteLine("  IT_MAX < 0.");
                return;
        }

        if (typeMethods.r8vec_any_negative(point_num, weight))
        {
            Console.WriteLine("");
            Console.WriteLine("KMEANS_W_01 - Fatal error!");
            Console.WriteLine("  Some weight entry is negative.");
            return;
        }

        if (typeMethods.r8vec_all_nonpositive(point_num, weight))
        {
            Console.WriteLine("");
            Console.WriteLine("KMEANS_W_01 - Fatal error!");
            Console.WriteLine("  No weight entry is positive.");
            return;
        }

        //
        //  Assign each point to the nearest cluster.
        //
        for (j = 0; j < point_num; j++)
        {
            point_energy_min = typeMethods.r8_huge();
            cluster[j] = -1;

            for (k = 0; k < cluster_num; k++)
            {
                point_energy = 0.0;
                for (i = 0; i < dim_num; i++)
                {
                    point_energy += Math.Pow(point[i + j * dim_num] - cluster_center[i + k * dim_num], 2);
                }

                if (point_energy < point_energy_min)
                {
                    point_energy_min = point_energy;
                    cluster[j] = k;
                }
            }
        }

        //
        //  Determine the cluster populations and weights.
        //
        typeMethods.i4vec_zero(cluster_num, ref cluster_population);
        cluster_weight = typeMethods.r8vec_zero_new(cluster_num);

        for (j = 0; j < point_num; j++)
        {
            k = cluster[j];
            cluster_population[k] += 1;
            cluster_weight[k] += weight[j];
        }

        //
        //  Calculate the mean and sum of squares for each cluster.
        //
        typeMethods.r8vec_zero(dim_num * cluster_num, ref cluster_center);

        for (j = 0; j < point_num; j++)
        {
            k = cluster[j];
            for (i = 0; i < dim_num; i++)
            {
                cluster_center[i + k * dim_num] += weight[j] * point[i + j * dim_num];
            }
        }

        for (k = 0; k < cluster_num; k++)
        {
            switch (cluster_weight[k])
            {
                case > 0.0:
                {
                    for (i = 0; i < dim_num; i++)
                    {
                        cluster_center[i + k * dim_num] /= cluster_weight[k];
                    }

                    break;
                }
            }
        }

        //
        //  Set the point energies.
        //
        f = typeMethods.r8vec_zero_new(point_num);

        for (j = 0; j < point_num; j++)
        {
            k = cluster[j];
            for (i = 0; i < dim_num; i++)
            {
                f[j] += Math.Pow(point[i + j * dim_num] - cluster_center[i + k * dim_num], 2);
            }
        }

        //
        //  Set the cluster energies.
        //
        typeMethods.r8vec_zero(cluster_num, ref cluster_energy);

        for (j = 0; j < point_num; j++)
        {
            k = cluster[j];
            cluster_energy[k] += weight[j] * f[j];
        }

        //
        //  Adjust the point energies by a weight factor.
        //
        for (j = 0; j < point_num; j++)
        {
            k = cluster[j];
            if (weight[j] < cluster_weight[k])
            {
                f[j] = f[j] * cluster_weight[k] / (cluster_weight[k] - weight[j]);
            }
        }

        //
        //  Examine each observation in turn to see if it should be
        //  reassigned to a different cluster.
        //
        it_num = 0;

        while (it_num < it_max)
        {
            it_num += 1;

            swap = 0;

            for (j = 0; j < point_num; j++)
            {
                il = cluster[j];
                ir = il;

                switch (cluster_population[il])
                {
                    case <= 1:
                        continue;
                }

                dc = f[j];

                for (k = 0; k < cluster_num; k++)
                {
                    if (k != il)
                    {
                        de = 0.0;
                        for (i = 0; i < dim_num; i++)
                        {
                            de += Math.Pow(point[i + j * dim_num] - cluster_center[i + k * dim_num], 2)
                                * cluster_weight[k] / (cluster_weight[k] + weight[j]);
                        }

                        if (de < dc)
                        {
                            dc = de;
                            ir = k;
                        }
                    }
                }

                //
                //  If the lowest value was obtained by staying in the current cluster,
                //  then cycle.
                //
                if (ir == il)
                {
                    continue;
                }

                //
                //  Reassign the point from cluster IL to cluster IR.
                //
                for (i = 0; i < dim_num; i++)
                {
                    cluster_center[i + il * dim_num] =
                        (cluster_weight[il] * cluster_center[i + il * dim_num]
                         - weight[j] * point[i + j * dim_num])
                        / (cluster_weight[il] - weight[j]);

                    cluster_center[i + ir * dim_num] =
                        (cluster_weight[ir] * cluster_center[i + ir * dim_num]
                         + weight[j] * point[i + j * dim_num])
                        / (cluster_weight[ir] + weight[j]);
                }

                cluster_weight[il] -= weight[j];
                cluster_weight[ir] += weight[j];

                cluster_energy[il] -= f[j];
                cluster_energy[ir] += dc;

                cluster_population[ir] += 1;
                cluster_population[il] -= 1;

                cluster[j] = ir;
                //
                //  Adjust the value of F for all points in clusters IL and IR.
                //
                for (j = 0; j < point_num; j++)
                {
                    k = cluster[j];

                    if (k == il || k == ir)
                    {
                        f[j] = 0.0;
                        for (i = 0; i < dim_num; i++)
                        {
                            f[j] += Math.Pow(point[i + j * dim_num] - cluster_center[i + k * dim_num], 2);
                        }

                        if (weight[j] < cluster_weight[k])
                        {
                            f[j] = f[j] * cluster_weight[k] / (cluster_weight[k] - weight[j]);
                        }
                    }
                }

                swap += 1;
            }

            //
            //  Exit if no reassignments were made during this iteration.
            //
            if (swap == 0)
            {
                break;
            }
        }

        //
        //  Compute the energy based on the final value of the cluster centers.
        //
        typeMethods.r8vec_zero(cluster_num, ref cluster_energy);

        for (j = 0; j < point_num; j++)
        {
            k = cluster[j];

            point_energy = 0.0;
            for (i = 0; i < dim_num; i++)
            {
                point_energy += Math.Pow(point[i + j * dim_num] - cluster_center[i + k * dim_num], 2);
            }

            cluster_energy[k] += weight[j] * point_energy;
        }
    }

    public static void kmeans_w_03(int dim_num, int point_num, int cluster_num, int it_max,
            ref int it_num, double[] point, double[] weight, ref int[] cluster,
            ref double[] cluster_center, ref int[] cluster_population, ref double[] cluster_energy)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //   KMEANS_W_03 applies the weighted K-Means algorithm.
        //
        //  Discussion:
        //
        //    The input data for the weight K-Means problem includes:
        //    * a set of N data points X in M dimensions, 
        //    * a set of N nonnegative weights W,
        //    * a desired number of clusters K.
        //    * an initial set of cluster centers Z,
        //    * an (optional) initial set of cluster assignments.
        //
        //    The goal is to determine K points Z, called cluster centers, and
        //    to assign each point X(I) to some cluster Z(J), so that we minimize
        //    the weighted standard deviation of the distance of each data point
        //    to the center of its cluster.  Writing J = CLUSTER(I) to
        //    indicate the index of the nearest cluster center Z(J) to the 
        //    point X(I), the quantity we are trying to minimize is the sum
        //    of the weighted cluster energies E(J), where:
        //
        //      E(J) = Sum ( 1 <= I <= N ) W(I) * || X(I) - Z(J) ||^2
        //
        //    Here, we assume that we are using the Euclidean norm, so that
        //    
        //      || X(I) - Z(J) ||^2 = Sum ( 1 <= K <= M )
        //        ( X(I)(K) - Z(J)(K) )^2
        //
        //    In this notation, X(I)(K) is the K-th spatial component of the
        //    I-th data point.
        //
        //    Note that this routine should give the same results as KMEANS_01
        //    in any case in which all the entries of the WEIGHT vector are equal.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Wendy Martinez, Angel Martinez,
        //    Computational Statistics Handbook with MATLAB,
        //    pages 373-376,
        //    Chapman and Hall / CRC, 2002.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the number of spatial dimensions.
        //
        //    Input, int POINT_NUM, the number of data points.
        //
        //    Input, int CLUSTER_NUM, the number of clusters.
        //
        //    Input, int IT_MAX, the maximum number of iterations.
        //
        //    Output, int IT_NUM, the number of iterations taken.
        //
        //    Input, double POINT[DIM_NUM*POINT_NUM], the data points.
        //
        //    Input, double WEIGHT[POINT_NUM], the weights.
        //
        //    Input/output, int CLUSTER[POINT_NUM], the cluster 
        //    to which each point belongs.  On output, these may have been altered.
        //
        //    Input/output, double CLUSTER_CENTER[DIM_NUM*CLUSTER_NUM], the
        //    centers associated with the clustering.  On output, these may
        //    have been altered.
        //
        //    Output, int CLUSTER_POPULATION[CLUSTER_NUM], the number
        //    of points in each cluster.
        //
        //    Output, double CLUSTER_ENERGY[CLUSTER_NUM], the energy of
        //    the clusters.
        //
    {
        int ci;
        int cj;
        double[] cluster_weight;
        double[] distsq;
        int i;
        int j;
        int k;
        double point_energy;
        double point_energy_min;
        int swap;
        switch (cluster_num)
        {
            //
            //  Check the input.
            //
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("KMEANS_W_03 - Fatal error!");
                Console.WriteLine("  CLUSTER_NUM < 1.");
                return;
        }

        switch (dim_num)
        {
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("KMEANS_W_03 - Fatal error!");
                Console.WriteLine("  DIM_NUM < 1.");
                return;
        }

        switch (point_num)
        {
            case < 1:
                Console.WriteLine("");
                Console.WriteLine("KMEANS_W_03 - Fatal error!");
                Console.WriteLine("  POINT_NUM < 1.");
                return;
        }

        switch (it_max)
        {
            case < 0:
                Console.WriteLine("");
                Console.WriteLine("KMEANS_W_03 - Fatal error!");
                Console.WriteLine("  IT_MAX < 0.");
                return;
        }

        if (typeMethods.r8vec_any_negative(point_num, weight))
        {
            Console.WriteLine("");
            Console.WriteLine("KMEANS_W_03 - Fatal error!");
            Console.WriteLine("  Some weight entry is negative.");
            return;
        }

        if (typeMethods.r8vec_all_nonpositive(point_num, weight))
        {
            Console.WriteLine("");
            Console.WriteLine("KMEANS_W_03 - Fatal error!");
            Console.WriteLine("  No weight entry is positive.");
            return;
        }

        //
        //  Assign each observation to the nearest cluster center.
        //
        for (j = 0; j < point_num; j++)
        {
            point_energy_min = typeMethods.r8_huge();
            cluster[j] = -1;

            for (k = 0; k < cluster_num; k++)
            {
                point_energy = 0.0;
                for (i = 0; i < dim_num; i++)
                {
                    point_energy += Math.Pow(point[i + j * dim_num] - cluster_center[i + k * dim_num], 2);
                }

                if (point_energy < point_energy_min)
                {
                    point_energy_min = point_energy;
                    cluster[j] = k;
                }
            }
        }

        //
        //  Determine the cluster populations and weights.
        //
        typeMethods.i4vec_zero(cluster_num, ref cluster_population);
        cluster_weight = typeMethods.r8vec_zero_new(cluster_num);

        for (j = 0; j < point_num; j++)
        {
            k = cluster[j];
            cluster_population[k] += 1;
            cluster_weight[k] += weight[j];
        }

        //
        //  Average the points in each cluster to get a new cluster center.
        //
        typeMethods.r8vec_zero(dim_num * cluster_num, ref cluster_center);

        for (j = 0; j < point_num; j++)
        {
            k = cluster[j];
            for (i = 0; i < dim_num; i++)
            {
                cluster_center[i + k * dim_num] += weight[j] * point[i + j * dim_num];
            }
        }

        for (k = 0; k < cluster_num; k++)
        {
            if (cluster_weight[k] != 0.0)
            {
                for (i = 0; i < dim_num; i++)
                {
                    cluster_center[i + k * dim_num] /= cluster_weight[k];
                }
            }
        }

        //
        //  Carry out the iteration.
        //
        it_num = 0;
        distsq = new double[cluster_num];

        while (it_num < it_max)
        {
            it_num += 1;

            swap = 0;

            for (j = 0; j < point_num; j++)
            {
                ci = cluster[j];

                switch (cluster_population[ci])
                {
                    case <= 1:
                        continue;
                }

                for (cj = 0; cj < cluster_num; cj++)
                {
                    if (cj == ci)
                    {
                        distsq[cj] = 0.0;
                        for (i = 0; i < dim_num; i++)
                        {
                            distsq[cj] += Math.Pow(point[i + j * dim_num] - cluster_center[i + cj * dim_num], 2);
                        }

                        distsq[cj] = distsq[cj] * cluster_weight[cj]
                                     / (cluster_weight[cj] - weight[j]);
                    }
                    else
                    {
                        switch (cluster_population[cj])
                        {
                            case 0:
                            {
                                for (i = 0; i < dim_num; i++)
                                {
                                    cluster_center[i + cj * dim_num] = point[i + j * dim_num];
                                }

                                distsq[cj] = 0.0;
                                break;
                            }
                            default:
                            {
                                distsq[cj] = 0.0;
                                for (i = 0; i < dim_num; i++)
                                {
                                    distsq[cj] += Math.Pow(point[i + j * dim_num] - cluster_center[i + cj * dim_num], 2);
                                }

                                distsq[cj] = distsq[cj] * cluster_weight[cj]
                                             / (cluster_weight[cj] + weight[j]);
                                break;
                            }
                        }
                    }
                }

                //
                //  Find the index of the minimum value of DISTSQ.
                //
                k = typeMethods.r8vec_min_index(cluster_num, distsq);
                //
                //  If that is not the cluster to which point I now belongs, move it there.
                //
                if (k == ci)
                {
                    continue;
                }

                cj = k;

                for (i = 0; i < dim_num; i++)
                {
                    cluster_center[i + ci * dim_num] =
                        (cluster_weight[ci] * cluster_center[i + ci * dim_num]
                         - weight[j] * point[i + j * dim_num])
                        / (cluster_weight[ci] - weight[j]);

                    cluster_center[i + cj * dim_num] =
                        (cluster_weight[cj] * cluster_center[i + cj * dim_num]
                         + weight[j] * point[i + j * dim_num])
                        / (cluster_weight[cj] + weight[j]);
                }

                cluster_population[ci] -= 1;
                cluster_population[cj] += 1;

                cluster_weight[ci] -= weight[j];
                cluster_weight[cj] += weight[j];

                cluster[j] = cj;

                swap += 1;
            }

            //
            //  Exit if no reassignments were made during this iteration.
            //
            if (swap == 0)
            {
                break;
            }
        }

        //
        //  Compute the energy based on the final value of the cluster centers.
        //
        typeMethods.r8vec_zero(cluster_num, ref cluster_energy);

        for (j = 0; j < point_num; j++)
        {
            k = cluster[j];

            point_energy = 0.0;
            for (i = 0; i < dim_num; i++)
            {
                point_energy += Math.Pow(point[i + j * dim_num] - cluster_center[i + k * dim_num], 2);
            }

            cluster_energy[k] += weight[j] * point_energy;
        }
    }
}