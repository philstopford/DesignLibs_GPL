using System;
using Burkardt.Sequence;
using Burkardt.Types;
using Burkardt.Uniform;
using entropyRNG;

namespace Burkardt.Tessellation
{
    public class CVTHaltonData
    {
        public BTupleData tupledata = new BTupleData();
        public int[] halton_base;
        public int[] halton_leap;
        public int[] halton_seed;
        public int ngrid ;
        public int rank;
        public int[] tuple;

        public CVTHaltonData(int dim_num)
        {
            tupledata.base_ = new int[dim_num];
        }
    }
    
    public static class CentroidalVoronoi
    {
        public static void cvt(ref CVTHaltonData data, int dim_num, int n, int batch, int init, int sample, int sample_num,
            int it_max, int it_fixed, ref int seed, ref double[] r, ref int it_num, ref double it_diff,
        ref double energy )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CVT computes a Centroidal Voronoi Tessellation.
        //
        //  Discussion:
        //
        //    This routine initializes the data, and carries out the
        //    CVT iteration.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    23 June 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Qiang Du, Vance Faber, and Max Gunzburger,
        //    Centroidal Voronoi Tessellations: Applications and Algorithms,
        //    SIAM Review, Volume 41, 1999, pages 637-676.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int N, the number of Voronoi cells.
        //
        //    Input, int BATCH, sets the maximum number of sample points
        //    generated at one time.  It is inefficient to generate the sample
        //    points 1 at a time, but memory intensive to generate them all
        //    at once.  You might set BATCH to min ( SAMPLE_NUM, 10000 ), for instance.
        //    BATCH must be at least 1.
        //
        //    Input, int INIT, specifies how the points are to be initialized.
        //    -1, 'RANDOM', using C++ RANDOM function;
        //     0, 'UNIFORM', using a simple uniform RNG;
        //     1, 'HALTON', from a Halton sequence;
        //     2, 'GRID', points from a grid;
        //     3, 'USER', call "user" routine;
        //     4, points are already initialized on input.
        //
        //    Input, int SAMPLE, specifies how the sampling is done.
        //    -1, 'RANDOM', using C++ RANDOM function;
        //     0, 'UNIFORM', using a simple uniform RNG;
        //     1, 'HALTON', from a Halton sequence;
        //     2, 'GRID', points from a grid;
        //     3, 'USER', call "user" routine.
        //
        //    Input, int SAMPLE_NUM, the number of sample points.
        //
        //    Input, int IT_MAX, the maximum number of iterations.
        //
        //    Input, int IT_FIXED, the maximum number of iterations to take
        //    with a fixed set of sample points.
        //
        //    Input/output, int *SEED, the random number seed.
        //
        //    Input/output, double R[DIM_NUM*N], the approximate CVT points.
        //    If INIT = 4 on input, then it is assumed that these values have been
        //    initialized.  On output, the CVT iteration has been applied to improve
        //    the value of the points.
        //
        //    Output, int *IT_NUM, the number of iterations taken.  Generally,
        //    this will be equal to IT_MAX, unless the iteration tolerance was
        //    satisfied early.
        //
        //    Output, double *IT_DIFF, the L2 norm of the difference
        //    between the iterates.
        //
        //    Output, double *ENERGY,  the discrete "energy", divided
        //    by the number of sample points.
        //
        {
            bool DEBUG = true;
            bool initialize;
            int seed_base = 0;
            int seed_init = 0;

            if (batch < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("CVT - Fatal error!");
                Console.WriteLine("  The input value BATCH < 1.");
                return;
            }

            if (seed <= 0)
            {
                Console.WriteLine("");
                Console.WriteLine("CVT - Fatal error!");
                Console.WriteLine("  The input value SEED <= 0.");
                return;
            }

            if (DEBUG)
            {
                Console.WriteLine("");
                Console.WriteLine("  Step       SEED          L2-Change        Energy");
                Console.WriteLine("");
            }

            it_num = 0;
            it_diff = 0.0;
            energy = 0.0;
            seed_init = seed;
            //
            //  Initialize the data, unless the user has already done that.
            //
            if (init != 4)
            {
                initialize = true;
                cvt_sample(ref data, dim_num, n, n, init, initialize, ref seed, r);
            }

            if (DEBUG)
            {
                Console.WriteLine("  "
                     + it_num.ToString().PadLeft(4) + "  "
                     + seed_init.ToString().PadLeft(12) + "");
            }

            //
            //  If the initialization and sampling steps use the same random number
            //  scheme, then the sampling scheme does not have to be initialized.
            //
            if (init == sample)
            {
                initialize = false;
            }
            else
            {
                initialize = true;
            }

            //
            //  Carry out the iteration.
            //
            while (it_num < it_max)
            {
                //
                //  If it's time to update the seed, save its current value
                //  as the starting value for all iterations in this cycle.
                //  If it's not time to update the seed, restore it to its initial
                //  value for this cycle.
                //
                if (((it_num) % it_fixed) == 0)
                {
                    seed_base = seed;
                }
                else
                {
                    seed = seed_base;
                }

                it_num = it_num + 1;
                seed_init = seed;

                cvt_iterate(ref data, dim_num, n, batch, sample, initialize, sample_num, ref seed,
                    r, ref it_diff, ref energy);

                initialize = false;

                if (DEBUG)
                {
                    Console.WriteLine("  "
                         + it_num.ToString().PadLeft(4) + "  "
                         + seed_init.ToString().PadLeft(12) + "  "
                         + it_diff.ToString().PadLeft(14) + "  "
                         + energy.ToString().PadLeft(14) + "");
                }
            }
        }

        public static double cvt_energy(ref CVTHaltonData data, int dim_num, int n, int batch, int sample, bool initialize,
                int sample_num, ref int seed, double[] r)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CVT_ENERGY computes the CVT energy of a dataset.
            //
            //  Discussion:
            //
            //    For a given number of generators, a CVT is a minimizer (or at least
            //    a local minimizer) of the CVT energy.  During a CVT iteration,
            //    it should generally be the case that the CVT energy decreases from
            //    step to step, and that perturbations or adjustments of an
            //    approximate CVT will almost always have higher CVT energy.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    03 December 2004
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
            //    Input, int BATCH, the maximum number of sample points to generate
            //    at one time.
            //
            //    Input, int SAMPLE, specifies how the sampling is done.
            //    -1, 'RANDOM', using C++ RANDOM function;
            //     0, 'UNIFORM', using a simple uniform RNG;
            //     1, 'HALTON', from a Halton sequence;
            //     2, 'GRID', points from a grid;
            //     3, 'USER', call "user" routine.
            //
            //    Input, bool INITIALIZE, is TRUE if the pseudorandom process 
            //    should be reinitialized.
            //
            //    Input, int SAMPLE_NUM, the number of sample points to use.
            //
            //    Input/output, int *SEED, a seed for the random number generator.
            //
            //    Input, double R[DIM_NUM*N], the coordinates of the points.
            //
            //    Output, double CVT_ENERGY, the estimated CVT energy.
            //
        {
            double energy;
            int get;
            int have;
            int i;
            int j;
            int[] nearest;
            double[] s;

            nearest = new int[batch];
            s = new double [dim_num * batch];

            have = 0;
            energy = 0.0;
            
            while (have < sample_num)
            {
                get = Math.Min(sample_num - have, batch);

                cvt_sample(ref data, dim_num, sample_num, get, sample, initialize, ref seed, s);

                have = have + get;

                find_closest(dim_num, n, get, s, r, ref nearest);

                for (j = 0; j < get; j++)
                {
                    for (i = 0; i < dim_num; i++)
                    {
                        energy = energy + (s[i + j * dim_num] - r[i + nearest[j] * dim_num])
                            * (s[i + j * dim_num] - r[i + nearest[j] * dim_num]);
                    }
                }

            }

            energy = energy / (double) (sample_num);

            return energy;
        }

        public static void cvt_iterate(ref CVTHaltonData data, int dim_num, int n, int batch, int sample, bool initialize,
            int sample_num, ref int seed, double[] r, ref double it_diff, ref double energy )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CVT_ITERATE takes one step of the CVT iteration.
        //
        //  Discussion:
        //
        //    The routine is given a set of points, called "generators", which
        //    define a tessellation of the region into Voronoi cells.  Each point
        //    defines a cell.  Each cell, in turn, has a centroid, but it is
        //    unlikely that the centroid and the generator coincide.
        //
        //    Each time this CVT iteration is carried out, an attempt is made
        //    to modify the generators in such a way that they are closer and
        //    closer to being the centroids of the Voronoi cells they generate.
        //
        //    A large number of sample points are generated, and the nearest generator
        //    is determined.  A count is kept of how many points were nearest to each
        //    generator.  Once the sampling is completed, the location of all the
        //    generators is adjusted.  This step should decrease the discrepancy
        //    between the generators and the centroids.
        //
        //    The centroidal Voronoi tessellation minimizes the "energy",
        //    defined to be the integral, over the region, of the square of
        //    the distance between each point in the region and its nearest generator.
        //    The sampling technique supplies a discrete estimate of this
        //    energy.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    20 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Qiang Du, Vance Faber, and Max Gunzburger,
        //    Centroidal Voronoi Tessellations: Applications and Algorithms,
        //    SIAM Review, Volume 41, 1999, pages 637-676.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int N, the number of Voronoi cells.
        //
        //    Input, int BATCH, sets the maximum number of sample points
        //    generated at one time.  It is inefficient to generate the sample
        //    points 1 at a time, but memory intensive to generate them all
        //    at once.  You might set BATCH to min ( SAMPLE_NUM, 10000 ), for instance.
        //    BATCH must be at least 1.
        //
        //    Input, int SAMPLE, specifies how the sampling is done.
        //    -1, 'RANDOM', using C++ RANDOM function;
        //     0, 'UNIFORM', using a simple uniform RNG;
        //     1, 'HALTON', from a Halton sequence;
        //     2, 'GRID', points from a grid;
        //     3, 'USER', call "user" routine.
        //
        //    Input, bool INITIALIZE, is TRUE if the SEED must be reset to SEED_INIT
        //    before computation.  Also, the pseudorandom process may need to be
        //    reinitialized.
        //
        //    Input, int SAMPLE_NUM, the number of sample points.
        //
        //    Input/output, int *SEED, the random number seed.
        //
        //    Input/output, double R[DIM_NUM*N], the Voronoi
        //    cell generators.  On output, these have been modified
        //
        //    Output, double *IT_DIFF, the L2 norm of the difference
        //    between the iterates.
        //
        //    Output, double *ENERGY,  the discrete "energy", divided
        //    by the number of sample points.
        //
        {
            int[] count;
            int get;
            int have;
            int i;
            int j;
            int j2;
            int[] nearest;
            double[] r2;
            double[] s;
            double term;
            //
            //  Take each generator as the first sample point for its region.
            //  This can slightly slow the convergence, but it simplifies the
            //  algorithm by guaranteeing that no region is completely missed
            //  by the sampling.
            //
            energy = 0.0;
            r2 = new double[dim_num * n];
            count = new int[n];
            nearest = new int[sample_num];
            s = new double[dim_num * sample_num];

            for (j = 0; j < n; j++)
            {
                for (i = 0; i < dim_num; i++)
                {
                    r2[i + j * dim_num] = r[i + j * dim_num];
                }
            }

            for (j = 0; j < n; j++)
            {
                count[j] = 1;
            }

            //
            //  Generate the sampling points S.
            //
            have = 0;
            
            while (have < sample_num)
            {
                get = Math.Min(sample_num - have, batch);

                cvt_sample(ref data, dim_num, sample_num, get, sample, initialize, ref seed, s);

                initialize = false;
                have = have + get;
                //
                //  Find the index N of the nearest cell generator to each sample point S.
                //
                find_closest(dim_num, n, get, s, r, ref nearest);
                //
                //  Add S to the centroid associated with generator N.
                //
                for (j = 0; j < get; j++)
                {
                    j2 = nearest[j];
                    for (i = 0; i < dim_num; i++)
                    {
                        r2[i + j2 * dim_num] = r2[i + j2 * dim_num] + s[i + j * dim_num];
                    }

                    for (i = 0; i < dim_num; i++)
                    {
                        energy = energy + Math.Pow(r[i + j2 * dim_num] - s[i + j * dim_num], 2);
                    }

                    count[j2] = count[j2] + 1;
                }
            }

            //
            //  Estimate the centroids.
            //
            for (j = 0; j < n; j++)
            {
                for (i = 0; i < dim_num; i++)
                {
                    r2[i + j * dim_num] = r2[i + j * dim_num] / (double) (count[j]);
                }
            }

            //
            //  Determine the sum of the distances between generators and centroids.
            //
            it_diff = 0.0;

            for (j = 0; j < n; j++)
            {
                term = 0.0;
                for (i = 0; i < dim_num; i++)
                {
                    term = term + (r2[i + j * dim_num] - r[i + j * dim_num])
                        * (r2[i + j * dim_num] - r[i + j * dim_num]);
                }

                it_diff = it_diff + Math.Sqrt(term);
            }

            //
            //  Replace the generators by the centroids.
            //
            for (j = 0; j < n; j++)
            {
                for (i = 0; i < dim_num; i++)
                {
                    r[i + j * dim_num] = r2[i + j * dim_num];
                }
            }

            //
            //  Normalize the discrete energy estimate.
            //
            energy = energy / sample_num;
        }

        public static void cvt_sample(ref CVTHaltonData data, int dim_num, int n, int n_now, int sample, bool initialize,
                ref int seed, double[] r)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CVT_SAMPLE returns sample points.
            //
            //  Discussion:
            //
            //    N sample points are to be taken from the unit box of dimension DIM_NUM.
            //
            //    These sample points are usually created by a pseudorandom process
            //    for which the points are essentially indexed by a quantity called
            //    SEED.  To get N sample points, we generate values with indices
            //    SEED through SEED+N-1.
            //
            //    It may not be practical to generate all the sample points in a 
            //    single call.  For that reason, the routine allows the user to
            //    request a total of N points, but to require that only N_NOW be
            //    generated now (on this call).  
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    23 June 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the spatial dimension.
            //
            //    Input, int N, the number of Voronoi cells.
            //
            //    Input, int N_NOW, the number of sample points to be generated
            //    on this call.  N_NOW must be at least 1.
            //
            //    Input, int SAMPLE, specifies how the sampling is done.
            //    -1, 'RANDOM', using C++ RANDOM function;
            //     0, 'UNIFORM', using a simple uniform RNG;
            //     1, 'HALTON', from a Halton sequence;
            //     2, 'GRID', points from a grid;
            //     3, 'USER', call "user" routine.
            //
            //    Input, bool INITIALIZE, is TRUE if the pseudorandom process should be
            //    reinitialized.
            //
            //    Input/output, int *SEED, the random number seed.
            //
            //    Output, double R[DIM_NUM*N_NOW], the sample points.
            //
        {
            double exponent;
            int halton_step;
            int i;
            int j;
            int rank_max;

            if (n_now < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("CVT_SAMPLE - Fatal error!");
                Console.WriteLine("  N_NOW < 1.");
                return;
            }

            if (sample == -1)
            {
                /*
                if (initialize)
                {
                    random_initialize(ref seed);
                }
                */
                for (j = 0; j < n_now; j++)
                {
                    for (i = 0; i < dim_num; i++)
                    {
                        r[i + j * dim_num] = RNG.nextdouble();//(double) random() / (double) RAND_MAX;
                    }
                }

                seed = (seed) + n_now * dim_num;
            }
            else if (sample == 0)
            {
                r = UniformRNG.r8mat_uniform_01(dim_num, n_now, ref seed);
            }
            else if (sample == 1)
            {
                data.halton_seed = new int[dim_num];
                data.halton_leap = new int[dim_num];
                data.halton_base = new int[dim_num];

                halton_step = seed;

                for (i = 0; i < dim_num; i++)
                {
                    data.halton_seed[i] = 0;
                }

                for (i = 0; i < dim_num; i++)
                {
                    data.halton_leap[i] = 1;
                }

                for (i = 0; i < dim_num; i++)
                {
                    data.halton_base[i] = Prime.prime(i + 1);
                }

                Halton.i4_to_halton_sequence(dim_num, n_now, halton_step, data.halton_seed,
                    data.halton_leap, data.halton_base, ref r);

                seed = seed + n_now;
            }
            else if (sample == 2)
            {
                exponent = 1.0 / (double) (dim_num);
                data.ngrid = (int) Math.Pow((double) n, exponent);
                rank_max = (int) Math.Pow((double) data.ngrid, (double) dim_num);
                data.tuple = new int[dim_num];

                if (rank_max < n)
                {
                    data.ngrid = data.ngrid + 1;
                    rank_max = (int) Math.Pow((double) data.ngrid, (double) dim_num);
                }

                if (initialize)
                {
                    data.rank = -1;
                    BTuple.tuple_next_fast(ref data.tupledata, data.ngrid, dim_num, data.rank, ref data.tuple);
                }

                data.rank = (seed) % rank_max;

                for (j = 0; j < n_now; j++)
                {
                    BTuple.tuple_next_fast(ref data.tupledata, data.ngrid, dim_num, data.rank, ref data.tuple);
                    data.rank = data.rank + 1;
                    data.rank = data.rank % rank_max;
                    for (i = 0; i < dim_num; i++)
                    {
                        r[i + j * dim_num] = (2.0 * data.tuple[i] - 1) / (2.0 * data.ngrid);
                    }
                }

                seed = seed + n_now;
            }
            else if (sample == 3)
            {
                user(dim_num, n_now, ref seed, ref r);
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("CVT_SAMPLE - Fatal error!");
                Console.WriteLine("  The value of SAMPLE = " + sample + " is illegal.");
            }
        }
        
        public static void find_closest ( int dim_num, int n, int sample_num, double[] s, double[] r,
        ref int[] nearest )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FIND_CLOSEST finds the nearest R point to each S point.
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
        //    21 October 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int N, the number of cell generators.
        //
        //    Input, int SAMPLE_NUM, the number of sample points.
        //
        //    Input, double S[DIM_NUM*SAMPLE_NUM], the points to be checked.
        //
        //    Input, double R[DIM_NUM*N], the cell generators.
        //
        //    Output, int NEAREST[SAMPLE_NUM], the (0-based) index of the nearest 
        //    cell generator.
        //
        {
            double dist_sq_min;
            double dist_sq;
            int i;
            int jr;
            int js;

            for ( js = 0; js < sample_num; js++ )
            {
                dist_sq_min = typeMethods.r8_huge();
                nearest[js] = -1;

                for ( jr = 0; jr < n; jr++ )
                {
                    dist_sq = 0.0;
                    for ( i = 0; i < dim_num; i++ )
                    {
                        dist_sq = dist_sq + ( s[i+js*dim_num] - r[i+jr*dim_num] ) 
                            * ( s[i+js*dim_num] - r[i+jr*dim_num] );
                    }

                    if ( jr == 0 || dist_sq < dist_sq_min )
                    {
                        dist_sq_min = dist_sq;
                        nearest[js] = jr;
                    }
                }
            }
        }
        
        public static void user ( int dim_num, int n, ref int seed, ref double[] r )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    USER samples points in a user-specified region with given density.
            //
            //  Discussion:
            //
            //    This routine can be used to 
            //
            //    * specify an interesting initial configuration for the data,
            //      by specifing that USER be used for initialization (INIT = 3);
            //
            //    * specify the shape of the computational region, by specifying
            //      that sample points are to be generated by this routine, 
            //      (SAMPLE = 3) and then returning sample points uniformly at random.
            //
            //    * specify the distribution or density function, by specifying
            //      that sample points are to be generated by this routine, 
            //      (SAMPLE = 3 ) and then returning sample points according to a 
            //      given probability density function.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    23 June 2005
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, integer DIM_NUM, the spatial dimension.
            //
            //    Input, integer N, the number of sample points desired.
            //
            //    Input/output, int *SEED, the "seed" value.  On output, SEED has 
            //    been updated.
            //
            //    Output, double R[DIM_NUM*N], the sample values.
            //
        {
            double angle;
            int j;
            double radius;

            for ( j = 0; j < n; j++ )
            {
                angle = 2.0 * Math.PI * RNG.nextdouble();
                radius = Math.Sqrt ( RNG.nextdouble() );
                r[0+j*2] = radius * Math.Cos ( angle );
                r[1+j*2] = radius * Math.Sin ( angle );
            }
        }
        
    }
}