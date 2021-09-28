using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.Sampling;

namespace Burkardt.Tessellation
{
    public class LCVData
    {
        public RegionData data;
        public LCVData()
        {
            init();
        }
        public LCVData(int M)
        {
            init();
            data.tdata.base_ = new int[M];
            for (int i = 0; i < M; i++)
            {
                data.tdata.base_[i] = 1;
            }
        }

        void init()
        {
            data = new RegionData();
        }
    }
    
    public static class LatinCentroidalVoronoi
    {
        public static void cvt(ref LCVData data, int m, int n, int sample_function_init,
                int sample_function_cvt, int sample_num_cvt, int maxit, ref int seed,
                ref double[] generator)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    CVT computes a Centroidal Voronoi Tessellation.
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
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 May 2003
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int M, the spatial dimension.
            //
            //    Input, int N, the number of Voronoi cells.
            //
            //    Input, int SAMPLE_FUNCTION_INIT, generator initialization function:
            //    -1, initializing function is RANDOM (C++ STDLIB library function);
            //    0, initializing function is UNIFORM;
            //    1, initializing function is HALTON;
            //    2, initializing function is GRID;
            //    3, initial values are set by the user.
            //
            //    Input, int SAMPLE_FUNCTION_CVT, region sampling function:
            //    -1, sampling function is RANDOM (C++ STDLIB library function);
            //    0, sampling function is UNIFORM;
            //    1, sampling function is HALTON;
            //    2, sampling function is GRID;
            //
            //    Input, int SAMPLE_NUM_CVT, the number of sample points.
            //
            //    Input, int MAXIT, the maximum number of correction iterations
            //    used in the Voronoi calculation.
            //
            //    Input/output, int *SEED, the random number seed.
            //
            //    Input/output, double GENERATOR[M*N], the Voronoi cell generators.
            //    On input, if SAMPLE_FUNCTION_INIT = 3, the user has initialized these.
            //    On output, the values have been through the CVT iteration.
            //
        {
            double change_l2 = 0;
            int it;
            bool reset;
            bool verbose = true;

            //
            //  Initialize the generators.
            //
            if (sample_function_init != 3)
            {
                reset = true;

                Region.region_sampler(ref data.data, m, n, n, generator, sample_function_init, reset, ref seed);
            }

            //
            //  Carry out the iteration.
            //
            if (verbose)
            {
                Console.WriteLine("");
                Console.WriteLine("  STEP  L2 Change");
                Console.WriteLine("");
            }

            for (it = 1; it <= maxit; it++)
            {
                cvt_iteration(ref data, m, n, ref generator, sample_num_cvt, sample_function_cvt,
                    ref seed, ref change_l2);

                if (verbose)
                {
                    Console.WriteLine(it.ToString().PadLeft(4) + "  "
                        + change_l2.ToString().PadLeft(12) + "");
                }

            }
        }

        public static void cvt_iteration(ref LCVData data, int m, int n, ref double[] generator, int sample_num_cvt,
        int sample_function_cvt, ref int seed, ref double change_l2 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CVT_ITERATION takes one step of the CVT iteration.
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
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 May 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the spatial dimension.
        //
        //    Input, int N, the number of Voronoi cells.
        //
        //    Input/output, double GENERATOR[M*N], the Voronoi
        //    cell generators.  On output, these have been modified
        //
        //    Input, int SAMPLE_NUM_CVT, the number of sample points.
        //
        //    Input, int SAMPLE_FUNCTION_CVT, region sampling function:
        //    -1, sampling function is RANDOM (C++ STDLIB library function);
        //    0, sampling function is UNIFORM;
        //    1, sampling function is HALTON;
        //    2, sampling function is GRID;
        //
        //    Input/output, int *SEED, the random number seed.
        //
        //    Output, double *CHANGE_L2, the L2 norm of the difference between
        //    the input and output data.
        //
        {
            double[] generator2;
            int[] count;
            int i;
            int j;
            int k;
            int nearest;
            bool reset;
            double[] x;
            //
            generator2 = new double[m * n];
            for (k = 0; k < m * n; k++)
            {
                generator2[k] = 0.0;
            }

            count = new int[n];
            for (i = 0; i < n; i++)
            {
                count[i] = 0;
            }

            x = new double[m];

            reset = true;
            //
            //  If we are using the C++ random number generator, then initialize using the current seed.
            //  (Currently, if we are using RANDOM for both the initializing and sampling, we make this
            //  call twice, which is inefficient and possibly misleading.)
            //
            if (sample_function_cvt == -1)
            {
                seed = entropyRNG.RNG.nextint();
            }

            for (j = 0; j < sample_num_cvt; j++)
            {
                //
                //  Generate a sampling point X.
                //
                Region.region_sampler(ref data.data, m, 1, sample_num_cvt, x, sample_function_cvt, reset,
                    ref seed);

                reset = false;
                //
                //  Find the nearest cell generator.
                //
                nearest = find_closest(m, n, x, generator);
                //
                //  Add X to the averaging data for GENERATOR(*,NEAREST).
                //
                for (i = 0; i < m; i++)
                {
                    generator2[nearest * m + i] = generator2[nearest * m + i] + x[i];
                }

                count[nearest] = count[nearest] + 1;
            }

            //
            //  Compute the new generators.
            //
            for (j = 0; j < n; j++)
            {
                if (count[j] != 0)
                {
                    for (i = 0; i < m; i++)
                    {
                        generator2[j * m + i] = generator2[j * m + i] / ((double) count[j]);
                    }
                }
            }

            //
            //  Determine the change.
            //
            change_l2 = 0.0;
            for (k = 0; k < m * n; k++)
            {
                change_l2 = change_l2 + Math.Pow((generator2[k] - generator[k]), 2);
            }

            change_l2 = Math.Sqrt(change_l2);
            //
            //  Update.
            //
            for (k = 0; k < m * n; k++)
            {
                generator[k] = generator2[k];
            }
        }

        public static void cvt_write(int dim_num, int n, int batch, int seed_init, int seed,
            string init_string, int it_max, int it_fixed, int it_num,
            double it_diff, double energy, string sample_string, int sample_num,
            double[] r, string file_out_name, bool comment )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CVT_WRITE writes a CVT dataset to a file.
        //
        //  Discussion:
        //
        //    The initial lines of the file are comments, which begin with a
        //    "#" character.
        //
        //    Thereafter, each line of the file contains the M-dimensional
        //    components of the next entry of the dataset.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 July 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int N, the number of points.
        //
        //    Input, int BATCH, sets the maximum number of sample points
        //    generated at one time.  It is inefficient to generate the sample
        //    points 1 at a time, but memory intensive to generate them all
        //    at once.  You might set BATCH to min ( SAMPLE_NUM, 10000 ), for instance.
        //
        //    Input, int SEED_INIT, the initial random number seed.
        //
        //    Input, int SEED, the current random number seed.
        //
        //    Input, char *INIT_STRING, specifies how the initial
        //    generators are chosen:
        //    filename, by reading data from a file;
        //    'GRID', picking points from a grid;
        //    'HALTON', from a Halton sequence;
        //    'RANDOM', using the C++ RANDOM function;
        //    'UNIFORM', using a simple uniform RNG;
        //    'USER', call "user" routine.
        //
        //    Input, int IT_MAX, the maximum number of iterations allowed.
        //
        //    Input, int IT_FIXED, the number of iterations to take with a
        //    fixed set of sample points.
        //
        //    Input, int IT_NUM, the actual number of iterations taken.
        //
        //    Input, double IT_DIFF, the L2 norm of the change
        //    in the CVT coordinates on the last iteration.
        //
        //    Input, double *ENERGY,  the discrete "energy", divided
        //    by the number of sample points.
        //
        //    Input, char *SAMPLE_STRING, specifies how the region is sampled:
        //    'GRID', picking points from a grid;
        //    'HALTON', from a Halton sequence;
        //    'RANDOM', using the C++ RANDOM function;
        //    'UNIFORM', using a simple uniform RNG;
        //    'USER', call "user" routine.
        //
        //    Input, int SAMPLE_NUM, the number of sampling points used on
        //    each iteration.
        //
        //    Input, double R(DIM_NUM,N), the points.
        //
        //    Input, char *FILE_OUT_NAME, the name of the output file.
        //
        //    Input, bool COMMENT, is true if comments may be included in the file.
        //
        {
            List<string> file_out = new List<string>();
            int i;
            int j;
            string s;
            
            s = DateTime.Now.ToString();

            if (comment)
            {
                file_out.Add("#  " + file_out_name + "");
                file_out.Add("#  created by routine CVT_WRITE.C" + "");
                file_out.Add("#  at " + s + "");
                file_out.Add("#");

                file_out.Add("#  Dimension DIM_NUM =        " + dim_num + "");
                file_out.Add("#  Number of points N =       " + n + "");
                file_out.Add("#  Initial SEED_INIT =        " + seed_init + "");
                file_out.Add("#  Current SEED =             " + seed + "");
                file_out.Add("#  INIT =                    \"" + init_string + "\".");
                file_out.Add("#  Max iterations IT_MAX =    " + it_max + "");
                file_out.Add("#  IT_FIXED (fixed samples) = " + it_fixed + "");
                file_out.Add("#  Iterations IT_NUM =        " + it_num + "");
                file_out.Add("#  Difference IT_DIFF =       " + it_diff + "");
                file_out.Add("#  CVT ENERGY =               " + energy + "");
                file_out.Add("#  SAMPLE =                  \"" + sample_string + "\".");
                file_out.Add("#  Samples SAMPLE_NUM =       " + sample_num + "");
                file_out.Add("#  Sampling BATCH size =      " + batch + "");
                file_out.Add("#  EPSILON (unit roundoff) =  " + double.Epsilon + "");
                file_out.Add("#");
            }

            for (j = 0; j < n; j++)
            {
                string cout = "";
                for (i = 0; i < dim_num; i++)
                {
                    cout += r[i + j * dim_num].ToString().PadLeft(10) + "  ";
                }

                file_out.Add(cout);
            }

            try
            {
                File.WriteAllLines(file_out_name, file_out);
            }
            catch
            {
                Console.WriteLine("");
                Console.WriteLine("CVT_WRITE - Fatal error!");
                Console.WriteLine("  Could not open the output file.");
            }
        }
        
        public static int find_closest ( int m, int n, double[] x, double[] generator )

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

            for ( j = 0; j < n; j++ )
            {
                dist = 0.0;
                for ( i = 0; i < m; i++ )
                {
                    dist = dist + Math.Pow ( x[i] - generator[i+j*m], 2 );
                }

                if ( j == 0 || dist < dist_min )
                {
                    dist_min = dist;
                    nearest = j;
                }

            }

            return nearest;
        }
    }
}