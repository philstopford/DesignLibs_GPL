namespace Burkardt.AppliedStatistics
{
    public static partial class Algorithms
    {
        public static void clustr(double[] x, ref double[] d, ref double[] dev, ref int[] b, double[] f,
        ref int[] e, int observations, int variables, int clusters, int minobserv, int maxclusters )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CLUSTR uses the K-means algorithm to cluster data.
        //
        //  Discussion:
        //
        //    Given a matrix of I observations on J variables, the
        //    observations are allocated to N clusters in such a way that the
        //    within-cluster sum of squares is minimised.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    04 February 2008
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
        //    Input, double X[I*J], the observed data.
        //
        //    Input/output, double D[K*J], the cluster centers.
        //    On input, the user has chosen these.  On output, they have been
        //    updated.
        //
        //    Output, double DEV[K], the sums of squared deviations
        //    of observations from their cluster centers.
        //
        //    Output, int B[I], indicates the cluster to which
        //    each observation has been assigned.
        //
        //    Workspace, double F[I].
        //
        //    Output, int E[K], the number of observations assigned
        //    to each cluster.
        //
        //    Input, int I, the number of observations.
        //
        //    Input, int J, the number of variables.
        //
        //    Input, int N, the number of clusters.
        //
        //    Input, int NZ, the minimum number of observations
        //    which any cluster is allowed to have.
        //
        //    Input, int K, the maximum number of clusters.
        //
        {
            double big = 1.0E+10;

            for (int i = 1; i <= clusters; i++)
            {
                e[i - 1] = 0;
            }

            //
            //  For each observation, calculate the distance from each cluster
            //  center, and assign to the nearest.
            //
            for (int i = 1; i <= observations; i++)
            {
                f[i - 1] = 0.0;
                double da = big;

                for (int j = 1; j <= clusters; j++)
                {
                    double db = 0.0;
                    for (int k = 1; k <= variables; k++)
                    {
                        double dc = x[i - 1 + (k - 1) * observations] - d[j - 1 + (k - 1) * maxclusters];
                        db = db + dc * dc;
                    }

                    if (db < da)
                    {
                        da = db;
                        b[i - 1] = j;
                    }
                }

                int ig = b[i - 1];
                e[ig - 1] = e[ig - 1] + 1;
            }

            //
            //  Calculate the mean and sum of squares for each cluster.
            //
            for (int i = 1; i <= clusters; i++)
            {
                dev[i - 1] = 0.0;
                for (int j = 1; j <= variables; j++)
                {
                    d[i - 1 + (j - 1) * maxclusters] = 0.0;
                }
            }

            for (int i = 1; i <= observations; i++)
            {
                int ig = b[i - 1];
                for (int j = 1; j <= variables; j++)
                {
                    d[ig - 1 + (j - 1) * maxclusters] = d[ig - 1 + (j - 1) * maxclusters] + x[i - 1 + (j - 1) * observations];
                }
            }

            for (int i = 1; i <= variables; i++)
            {
                for (int j = 1; j <= clusters; j++)
                {
                    d[j - 1 + (i - 1) * maxclusters] = d[j - 1 + (i - 1) * maxclusters] / (double) e[j - 1];
                }
            }

            for (int i = 1; i <= variables; i++)
            {
                for (int j = 1; j <= observations; j++)
                {
                    int il = b[j - 1];
                    double da = x[j - 1 + (i - 1) * observations] - d[il - 1 + (i - 1) * maxclusters];
                    double db = da * da;
                    f[j - 1] = f[j - 1] + db;
                    dev[il - 1] = dev[il - 1] + db;
                }
            }

            for (int i = 1; i <= observations; i++)
            {
                int il = b[i - 1];
                double fl = e[il - 1];
                if (2.0 <= fl)
                {
                    f[i - 1] = f[i - 1] * fl / (fl - 1.0);
                }
            }

            //
            //  Examine each observation in turn to see if it should be
            //  reassigned to a different cluster.
            //
            for (;;)
            {
                int iw = 0;

                for (int i = 1; i <= observations; i++)
                {
                    int il = b[i - 1];
                    int ir = il;
                    //
                    //  If the number of cluster points is less than or equal to the
                    //  specified minimum, NZ, then bypass this iteration.
                    //
                    if (minobserv < e[il - 1])
                    {
                        double fl = e[il - 1];
                        double dc = f[i - 1];

                        for (int  j = 1; j <= clusters; j++ )
                        {
                            if ( j != il )
                            {
                                double fm = e[j -1];
                                fm = fm / (fm + 1.0);

                                double de = 0.0;
                                for (int k = 1; k <= variables; k++)
                                {
                                    double da = x[i - 1 + (k - 1) * observations] - d[j -1 + (k - 1) * maxclusters];
                                    de = de + da * da * fm;
                                }

                                if (de < dc)
                                {
                                    dc = de;
                                    ir = j;
                                }
                            }
                        }
                        //
                        //  Reassignment is made here if necessary.
                        //
                        if (ir != il)
                        {
                            double fq = e[ir - 1];
                            dev[il - 1] = dev[il - 1] - f[i - 1];
                            dev[ir - 1] = dev[ir - 1] + dc;
                            e[ir - 1] = e[ir - 1] + 1;
                            e[il - 1] = e[il - 1] - 1;

                            for (int j = 1; j <= variables; j++ )
                            {
                                d[il - 1 + ( j -1)*maxclusters] = (d[il - 1 + ( j -1) * maxclusters]
                                *fl - x[i - 1 + ( j -1)*observations] ) / (fl - 1.0);
                                d[ir - 1 + ( j -1)*maxclusters] = (d[ir - 1 + ( j -1) * maxclusters]
                                *fq + x[i - 1 + ( j -1)*observations] ) / (fq + 1.0);
                            }

                            b[i - 1] = ir;

                            for (int j = 1; j <= observations; j++)
                            {
                                int ij = b[j - 1];

                                if (ij == il || ij == ir)
                                {
                                    f[j - 1] = 0.0;
                                    for (int k = 1; k <= variables; k++)
                                    {
                                        double da = x[j - 1 + (k - 1) * observations] - d[ij - 1 + (k - 1) * maxclusters];
                                        f[j - 1] = f[j - 1] + da * da;
                                    }

                                    fl = e[ij - 1];
                                    f[j - 1] = f[j - 1] * fl / (fl - 1.0);
                                }
                            }

                            iw = iw + 1;
                        }
                    }
                }

                //
                //  If any reassignments were made on this pass, then do another pass.
                //
                if (iw == 0)
                {
                    break;
                }
            }
        }


    }
}