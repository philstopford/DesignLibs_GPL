﻿using System;
using Burkardt.Types;

namespace Burkardt.AppliedStatistics;

public static partial class Algorithms
{
    public static void kmns ( double[] a, int m, int n, ref double[] c, int k, ref int[] ic1, ref int[] nc,
            int iter, ref double[] wss, ref int ifault )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    KMNS carries out the K-means algorithm.
        //
        //  Discussion:
        //
        //    This routine attempts to divide M points in N-dimensional space into 
        //    K clusters so that the within cluster sum of squares is minimized.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 February 2008
        //
        //  Author:
        //
        //    Original FORTRAN77 version by John Hartigan, Manchek Wong.
        //    C++ version by John Burkardt
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
        //    Input, double A(M,N), the points.
        //
        //    Input, int M, the number of points.
        //
        //    Input, int N, the number of spatial dimensions.
        //
        //    Input/output, double C(K,N), the cluster centers.
        //
        //    Input, int K, the number of clusters.
        //
        //    Output, int IC1(M), the cluster to which each point 
        //    is assigned.
        //
        //    Output, int NC(K), the number of points in each cluster.
        //
        //    Input, int ITER, the maximum number of iterations allowed.
        //
        //    Output, double WSS(K), the within-cluster sum of squares
        //    of each cluster.
        //
        //    Output, int *IFAULT, error indicator.
        //    0, no error was detected.
        //    1, at least one cluster is empty after the initial assignment.  A better
        //       set of initial cluster centers is needed.
        //    2, the allowed maximum number off iterations was exceeded.
        //    3, K is less than or equal to 1, or greater than or equal to M.
        //
    {
        double[] dt = new double[2];

        ifault = 0;

        if (k <= 1 || m <= k)
        {
            ifault = 3;
            return;
        }

        int[] ic2 = new int[m];
        double[] an1 = new double[k];
        double[] an2 = new double[k];
        int[] ncp = new int[k];
        double[] d = new double[m];
        int[] itran = new int[k];
        int[] live = new int[k];
        //
        //  For each point I, find its two closest centers, IC1(I) and
        //  IC2(I).  Assign the point to IC1(I).
        //
        for (int i = 1; i <= m; i++)
        {
            ic1[i - 1] = 1;
            ic2[i - 1] = 2;

            for (int il = 1; il <= 2; il++)
            {
                dt[il - 1] = 0.0;
                for (int j = 1; j <= n; j++)
                {
                    double da = a[i - 1 + (j - 1) * m] - c[il - 1 + (j - 1) * k];
                    dt[il - 1] += da * da;
                }
            }

            if (dt[1] < dt[0])
            {
                ic1[i - 1] = 2;
                ic2[i - 1] = 1;
                (dt[0], dt[1]) = (dt[1], dt[0]);
            }

            for (int l = 3; l <= k; l++)
            {
                double db = 0.0;
                for (int j = 1; j <= n; j++)
                {
                    double dc = a[i - 1 + (j - 1) * m] - c[l - 1 + (j - 1) * k];
                    db += dc * dc;
                }

                if (!(db < dt[1]))
                {
                    continue;
                }

                if (dt[0] <= db)
                {
                    dt[1] = db;
                    ic2[i - 1] = l;
                }
                else
                {
                    dt[1] = dt[0];
                    ic2[i - 1] = ic1[i - 1];
                    dt[0] = db;
                    ic1[i - 1] = l;
                }
            }
        }

        //
        //  Update cluster centers to be the average of points contained within them.
        //
        for (int l = 1; l <= k; l++)
        {
            nc[l - 1] = 0;
            for (int j = 1; j <= n; j++)
            {
                c[l - 1 + (j - 1) * k] = 0.0;
            }
        }

        for (int i = 1; i <= m; i++)
        {
            int l = ic1[i - 1];
            nc[l - 1] += 1;
            for (int j = 1; j <= n; j++)
            {
                c[l - 1 + (j - 1) * k] += a[i - 1 + (j - 1) * m];
            }
        }

        //
        //  Check to see if there is any empty cluster at this stage.
        //
        ifault = 1;

        for (int l = 1; l <= k; l++)
        {
            switch (nc[l - 1])
            {
                case 0:
                    ifault = 1;
                    return;
            }
        }

        ifault = 0;

        for (int l = 1; l <= k; l++)
        {
            double aa = nc[l - 1];

            for (int j = 1; j <= n; j++)
            {
                c[l - 1 + (j - 1) * k] /= aa;
            }

            //
            //  Initialize AN1, AN2, ITRAN and NCP.
            //
            //  AN1(L) = NC(L) / (NC(L) - 1)
            //  AN2(L) = NC(L) / (NC(L) + 1)
            //  ITRAN(L) = 1 if cluster L is updated in the quick-transfer stage,
            //           = 0 otherwise
            //
            //  In the optimal-transfer stage, NCP(L) stores the step at which
            //  cluster L is last updated.
            //
            //  In the quick-transfer stage, NCP(L) stores the step at which
            //  cluster L is last updated plus M.
            //
            an2[l - 1] = aa / (aa + 1.0);

            an1[l - 1] = aa switch
            {
                > 1.0 => aa / (aa - 1.0),
                _ => typeMethods.r8_huge()
            };

            itran[l - 1] = 1;
            ncp[l - 1] = -1;
        }

        int indx = 0;
        ifault = 2;

        for (int ij = 1; ij <= iter; ij++)
        {
            //
            //  In this stage, there is only one pass through the data.   Each
            //  point is re-allocated, if necessary, to the cluster that will
            //  induce the maximum reduction in within-cluster sum of squares.
            //
            optra(a, m, n, ref c, k, ref ic1, ref ic2, ref nc, ref an1, ref an2, ref ncp, ref d, ref itran, ref live, ref indx);
            //
            //  Stop if no transfer took place in the last M optimal transfer steps.
            //
            if (indx == m)
            {
                ifault = 0;
                break;
            }

            //
            //  Each point is tested in turn to see if it should be re-allocated
            //  to the cluster to which it is most likely to be transferred,
            //  IC2(I), from its present cluster, IC1(I).   Loop through the
            //  data until no further change is to take place.
            //
            qtran(a, m, n, ref c, k, ref ic1, ref ic2, ref nc, ref an1, ref an2, ref ncp, ref d, ref itran, ref indx);
            //
            //  If there are only two clusters, there is no need to re-enter the
            //  optimal transfer stage.
            //
            if (k == 2)
            {
                ifault = 0;
                break;
            }

            //
            //  NCP has to be set to 0 before entering OPTRA.
            //
            for (int l = 1; l <= k; l++)
            {
                ncp[l - 1] = 0;
            }

        }

        switch (ifault)
        {
            //
            //  If the maximum number of iterations was taken without convergence,
            //  IFAULT is 2 now.  This may indicate unforeseen looping.
            //
            case 2:
                Console.WriteLine("");
                Console.WriteLine("KMNS - Warning!");
                Console.WriteLine("  Maximum number of iterations reached");
                Console.WriteLine("  without convergence.");
                break;
        }

        //
        //  Compute the within-cluster sum of squares for each cluster.
        //
        for (int l = 1; l <= k; l++)
        {
            wss[l - 1] = 0.0;
            for (int j = 1; j <= n; j++)
            {
                c[l - 1 + (j - 1) * k] = 0.0;
            }
        }

        for (int i = 1; i <= m; i++)
        {
            int ii = ic1[i - 1];
            for (int j = 1; j <= n; j++)
            {
                c[ii - 1 + (j - 1) * k] += a[i - 1 + (j - 1) * m];
            }
        }

        for (int j = 1; j <= n; j++)
        {
            for (int l = 1; l <= k; l++)
            {
                c[l - 1 + (j - 1) * k] /= nc[l - 1];
            }

            for (int i = 1; i <= m; i++)
            {
                int ii = ic1[i - 1];
                double da = a[i - 1 + (j - 1) * m] - c[ii - 1 + (j - 1) * k];
                wss[ii - 1] += da * da;
            }
        }

    }

    public static void optra(double[] a, int m, int n, ref double[] c, int k, ref int[] ic1,
            ref int[] ic2, ref int[] nc, ref double[] an1, ref double[] an2, ref int[] ncp, ref double[] d,
            ref int[] itran, ref int[] live, ref int indx )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    OPTRA carries out the optimal transfer stage.
        //
        //  Discussion:
        //
        //    This is the optimal transfer stage.
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
        //    15 February 2008
        //
        //  Author:
        //
        //    Original FORTRAN77 version by John Hartigan, Manchek Wong.
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
        //    Input, double A(M,N), the points.
        //
        //    Input, int M, the number of points.
        //
        //    Input, int N, the number of spatial dimensions.
        //
        //    Input/output, double C(K,N), the cluster centers.
        //
        //    Input, int K, the number of clusters.
        //
        //    Input/output, int IC1(M), the cluster to which each 
        //    point is assigned.
        //
        //    Input/output, int IC2(M), used to store the cluster 
        //    which each point is most likely to be transferred to at each step.
        //
        //    Input/output, int NC(K), the number of points in 
        //    each cluster.
        //
        //    Input/output, double AN1(K).
        //
        //    Input/output, double AN2(K).
        //
        //    Input/output, int NCP(K).
        //
        //    Input/output, double D(M).
        //
        //    Input/output, int ITRAN(K).
        //
        //    Input/output, int LIVE(K).
        //
        //    Input/output, int *INDX, the number of steps since a 
        //    transfer took place.
        //
    {
        //
        //  If cluster L is updated in the last quick-transfer stage, it
        //  belongs to the live set throughout this stage.   Otherwise, at
        //  each step, it is not in the live set if it has not been updated
        //  in the last M optimal transfer steps.
        //
        for (int l = 1; l <= k; l++)
        {
            live[l - 1] = itran[l - 1] switch
            {
                1 => m + 1,
                _ => live[l - 1]
            };
        }

        for (int i = 1; i <= m; i++)
        {
            indx += 1;
            int l1 = ic1[i - 1];
            int l2 = ic2[i - 1];
            int ll = l2;
            switch (nc[l1 - 1])
            {
                //
                //  If point I is the only member of cluster L1, no transfer.
                //
                case > 1:
                {
                    //
                    //  If L1 has not yet been updated in this stage, no need to
                    //  re-compute D(I).
                    //
                    if (ncp[l1 - 1] != 0)
                    {
                        double de = 0.0;
                        for (int j = 1; j <= n; j++)
                        {
                            double df = a[i - 1 + (j - 1) * m] - c[l1 - 1 + (j - 1) * k];
                            de += df * df;
                        }

                        d[i - 1] = de * an1[l1 - 1];
                    }

                    //
                    //  Find the cluster with minimum R2.
                    //
                    double da = 0.0;
                    for (int j = 1; j <= n; j++)
                    {
                        double db = a[i - 1 + (j - 1) * m] - c[l2 - 1 + (j - 1) * k];
                        da += db * db;
                    }

                    double r2 = da * an2[l2 - 1];

                    for (int l = 1; l <= k; l++)
                    {
                        //
                        //  If LIVE(L1) <= I, then L1 is not in the live set.   If this is
                        //  true, we only need to consider clusters that are in the live set
                        //  for possible transfer of point I.   Otherwise, we need to consider
                        //  all possible clusters.
                        //
                        if ((i >= live[l1 - 1] && i >= live[l2 - 1]) || l == l1 || l == ll)
                        {
                            continue;
                        }

                        double rr = r2 / an2[l - 1];

                        double dc = 0.0;
                        for (int j = 1; j <= n; j++)
                        {
                            double dd = a[i - 1 + (j - 1) * m] - c[l - 1 + (j - 1) * k];
                            dc += dd * dd;
                        }

                        if (!(dc < rr))
                        {
                            continue;
                        }

                        r2 = dc * an2[l - 1];
                        l2 = l;
                    }

                    //
                    //  If no transfer is necessary, L2 is the new IC2(I).
                    //
                    if (d[i - 1] <= r2)
                    {
                        ic2[i - 1] = l2;
                    }
                    //
                    //  Update cluster centers, LIVE, NCP, AN1 and AN2 for clusters L1 and
                    //  L2, and update IC1(I) and IC2(I).
                    //
                    else
                    {
                        indx = 0;
                        live[l1 - 1] = m + i;
                        live[l2 - 1] = m + i;
                        ncp[l1 - 1] = i;
                        ncp[l2 - 1] = i;
                        double al1 = nc[l1 - 1];
                        double alw = al1 - 1.0;
                        double al2 = nc[l2 - 1];
                        double alt = al2 + 1.0;
                        for (int j = 1; j <= n; j++)
                        {
                            c[l1 - 1 + (j - 1) * k] = (c[l1 - 1 + (j - 1) * k] * al1 - a[i - 1 + (j - 1) * m]) / alw;
                            c[l2 - 1 + (j - 1) * k] = (c[l2 - 1 + (j - 1) * k] * al2 + a[i - 1 + (j - 1) * m]) / alt;
                        }

                        nc[l1 - 1] -= 1;
                        nc[l2 - 1] += 1;
                        an2[l1 - 1] = alw / al1;
                        an1[l1 - 1] = alw switch
                        {
                            > 1.0 => alw / (alw - 1.0),
                            _ => typeMethods.r8_huge()
                        };

                        an1[l2 - 1] = alt / al2;
                        an2[l2 - 1] = alt / (alt + 1.0);
                        ic1[i - 1] = l2;
                        ic2[i - 1] = l1;
                    }

                    break;
                }
            }

            if (indx == m)
            {
                return;
            }
        }

        //
        //  ITRAN(L) = 0 before entering QTRAN.   Also, LIVE(L) has to be
        //  decreased by M before re-entering OPTRA.
        //
        for (int l = 1; l <= k; l++)
        {
            itran[l - 1] = 0;
            live[l - 1] -= m;
        }
    }

    public static void qtran(double[] a, int m, int n, ref double[] c, int k, ref int[] ic1,
            ref int[] ic2, ref int[] nc, ref double[] an1, ref double[] an2, ref int[] ncp, ref double[] d,
            ref int[] itran, ref int indx )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    QTRAN carries out the quick transfer stage.
        //
        //  Discussion:
        //
        //    This is the quick transfer stage.
        //
        //    IC1(I) is the cluster which point I belongs to.
        //    IC2(I) is the cluster which point I is most likely to be
        //    transferred to.
        //
        //    For each point I, IC1(I) and IC2(I) are switched, if necessary, to
        //    reduce within-cluster sum of squares.  The cluster centers are
        //    updated after each step.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    15 February 2008
        //
        //  Author:
        //
        //    Original FORTRAN77 version by John Hartigan, Manchek Wong.
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
        //    Input, double A(M,N), the points.
        //
        //    Input, int M, the number of points.
        //
        //    Input, int N, the number of spatial dimensions.
        //
        //    Input/output, double C(K,N), the cluster centers.
        //
        //    Input, int K, the number of clusters.
        //
        //    Input/output, int IC1(M), the cluster to which each 
        //    point is assigned.
        //
        //    Input/output, int IC2(M), used to store the cluster 
        //    which each point is most likely to be transferred to at each step.
        //
        //    Input/output, int NC(K), the number of points in 
        //    each cluster.
        //
        //    Input/output, double AN1(K).
        //
        //    Input/output, double AN2(K).
        //
        //    Input/output, int NCP(K).
        //
        //    Input/output, double D(M).
        //
        //    Input/output, int ITRAN(K).
        //
        //    Input/output, int INDX, counts the number of steps 
        //    since the last transfer.
        //
    {
        //
        //  In the optimal transfer stage, NCP(L) indicates the step at which
        //  cluster L is last updated.   In the quick transfer stage, NCP(L)
        //  is equal to the step at which cluster L is last updated plus M.
        //
        int icoun = 0;
        int istep = 0;

        for (;;)
        {
            int i;
            for (i = 1; i <= m; i++)
            {
                icoun += 1;
                istep += 1;
                int l1 = ic1[i - 1];
                int l2 = ic2[i - 1];
                switch (nc[l1 - 1])
                {
                    //
                    //  If point I is the only member of cluster L1, no transfer.
                    //
                    case > 1:
                    {
                        //
                        //  If NCP(L1) < ISTEP, no need to re-compute distance from point I to
                        //  cluster L1.   Note that if cluster L1 is last updated exactly M
                        //  steps ago, we still need to compute the distance from point I to
                        //  cluster L1.
                        //
                        int j;
                        if (istep <= ncp[l1 - 1])
                        {
                            double da = 0.0;
                            for (j = 1; j <= n; j++)
                            {
                                double db = a[i - 1 + (j - 1) * m] - c[l1 - 1 + (j - 1) * k];
                                da += db * db;
                            }

                            d[i - 1] = da * an1[l1 - 1];
                        }

                        //
                        //  If NCP(L1) <= ISTEP and NCP(L2) <= ISTEP, there will be no transfer of
                        //  point I at this step.
                        //
                        if (istep < ncp[l1 - 1] || istep < ncp[l2 - 1])
                        {
                            double r2 = d[i - 1] / an2[l2 - 1];

                            double dd = 0.0;
                            for (j = 1; j <= n; j++)
                            {
                                double de = a[i - 1 + (j - 1) * m] - c[l2 - 1 + (j - 1) * k];
                                dd += de * de;
                            }

                            //
                            //  Update cluster centers, NCP, NC, ITRAN, AN1 and AN2 for clusters
                            //  L1 and L2.   Also update IC1(I) and IC2(I).   Note that if any
                            //  updating occurs in this stage, INDX is set back to 0.
                            //
                            if (dd < r2)
                            {
                                icoun = 0;
                                indx = 0;
                                itran[l1 - 1] = 1;
                                itran[l2 - 1] = 1;
                                ncp[l1 - 1] = istep + m;
                                ncp[l2 - 1] = istep + m;
                                double al1 = nc[l1 - 1];
                                double alw = al1 - 1.0;
                                double al2 = nc[l2 - 1];
                                double alt = al2 + 1.0;
                                for (j = 1; j <= n; j++)
                                {
                                    c[l1 - 1 + (j - 1) * k] =
                                        (c[l1 - 1 + (j - 1) * k] * al1 - a[i - 1 + (j - 1) * m]) / alw;
                                    c[l2 - 1 + (j - 1) * k] =
                                        (c[l2 - 1 + (j - 1) * k] * al2 + a[i - 1 + (j - 1) * m]) / alt;
                                }

                                nc[l1 - 1] -= 1;
                                nc[l2 - 1] += 1;
                                an2[l1 - 1] = alw / al1;
                                an1[l1 - 1] = alw switch
                                {
                                    > 1.0 => alw / (alw - 1.0),
                                    _ => typeMethods.r8_huge()
                                };

                                an1[l2 - 1] = alt / al2;
                                an2[l2 - 1] = alt / (alt + 1.0);
                                ic1[i - 1] = l2;
                                ic2[i - 1] = l1;
                            }
                        }

                        break;
                    }
                }

                //
                //  If no re-allocation took place in the last M steps, return.
                //
                if (icoun == m)
                {
                    return;
                }
            }
        }
    }
}