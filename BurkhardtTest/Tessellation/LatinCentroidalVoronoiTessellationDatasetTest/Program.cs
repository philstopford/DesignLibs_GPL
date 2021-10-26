using System;
using Burkardt.Sampling;
using Burkardt.SolveNS;
using Burkardt.Tessellation;
using Burkardt.Types;

namespace LatinCentroidalVoronoiTessellationDatasetTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for LCVT_DATASET.
            //
            //  Discussion:
            //
            //    LCVT_DATASET computes a Latinized CVT dataset and writes it to a file.
            //
            //    This program is meant to be used interactively.  It's also
            //    possible to prepare a simple input file beforehand and use it
            //    in batch mode.
            //
            //    The program requests input values from the user:
            //
            //    * DIM_NUM, the spatial dimension;
            //    * N, the number of points to generate;
            //    * SEED_INIT, a seed to use for random number generation;
            //    * INIT, initialize the points:
            //      ** file, by reading data from file;
            //      ** GRID, picking points from a grid;
            //      ** HALTON, from a Halton sequence;
            //      ** RANDOM, using C++ RANDOM function;
            //      ** UNIFORM, using a simple uniform RNG;
            //      ** USER, call the "user" routine;
            //    * CVT_IT_NUM, the maximum number of iterations;
            //    * SAMPLE, how to conduct the sampling:
            //      ** GRID, picking points from a grid;
            //      ** HALTON, from a Halton sequence;
            //      ** RANDOM, using C++ RANDOM function;
            //      ** UNIFORM, using a simple uniform RNG;
            //      ** USER, call the "user" routine.
            //    * SAMPLE_NUM, the number of sampling points;
            //    * BATCH, the number of sampling points to create at one time.
            //    * LAT_IT_NUM, the maximum number of iterations;
            //    * OUTPUT, a file in which to store the data.
            //
            //    To indicate that no further computations are desired, it is
            //    enough to input a nonsensical value, such as -1.
            //
            //  Modified:
            //
            //    12 April 2007
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    double CELL_GENERATOR(M,N), the Voronoi cell generators
            //    of the Voronoi tessellation, as approximated by the CVT algorithm.  This
            //    is the output quantity of most interest.
            //
            //    int CVT_IT, the number of iterations used in the Centroidal
            //    Voronoi Tesselation calculation.  The default value is 10.
            //
            //    int LATIN_IT, the number of Latin hypercube iterations to carry out.
            //    This defaults to 5.
            //
            //    int M, the spatial dimension.
            //
            //    int N, the number of Voronoi cells to generate.
            //
            //    int SAMPLE_FUNCTION_CVT, specifies how the region is sampled:
            //    -1, the sampling function is C++ RANDOM_NUMBER,
            //    0, the sampling function is UNIFORM,
            //    1, the sampling function is HALTON,
            //    2, the sampling function is GRID.
            //
            //    int SAMPLE_FUNCTION_INIT, specifies how the initial
            //    generators are chosen:
            //    -1, the initialization function is C++ RANDOM,
            //    0, the initialization function is UNIFORM,
            //    1, the initialization function is HALTON,
            //    2, the initialization function is GRID,
            //    3, the initial values are read in from a file.
            //
            //    int SAMPLE_NUM, the number of sampling points used on
            //    each CVT iteration.  A typical value is 5000 * N.
            //
            //    int SEED, determines how to initialize the random number routine.
            //    If SEED is zero, then RANDOM_INITIALIZE will make up a seed
            //    from the current real time clock reading.
            //    If SEED is nonzero, then a reproducible sequence of random numbers
            //    defined by SEED will be chosen.
            //    By default, SEED initially has a value chosen by RANDOM_INITIALIZE,
            //    but the user can reset SEED at any time.
            //
        {
            int batch;
            double cvt_energy = 0;
            double cvt_it_diff = 0;
            int cvt_it;
            int cvt_it_num;
            bool debug = true;
            int dim_num;
            string input_file_name = "";
            int init;
            string init_string;
            double lat_energy = 0;
            int lat_it;
            int lat_it_num;
            int n;
            int n_total;
            string output_file_name = "";
            double[] r;
            bool reset;
            int sample;
            int sample_num;
            string sample_string;
            int seed;
            int seed_init;
            RegionData rdata = new RegionData();
            LCVData ldata = new LCVData();
            //
            //  Print introduction and options.
            //

            Console.WriteLine("");
            Console.WriteLine("LCVT_DATASET");
            Console.WriteLine("  C++ version");
            Console.WriteLine("  Create a \"Latinized\" CVT datasets.");
            Console.WriteLine("");
            Console.WriteLine("  This program is meant to be used interactively.");
            Console.WriteLine("  It is also possible to prepare a simple input ");
            Console.WriteLine("  file beforehand and use it in batch mode.");
            Console.WriteLine("");
            Console.WriteLine("  The program requests input values from the user:");
            Console.WriteLine("");
            Console.WriteLine("  * DIM_NUM, the spatial dimension,");
            Console.WriteLine("  * N, the number of points to generate,");
            Console.WriteLine("  * SEED_INIT, a seed to use for random number generation,");
            Console.WriteLine("  * INIT, initialize the points:");
            Console.WriteLine("    ** file, read data from a file;");
            Console.WriteLine("    ** GRID, by picking points from a grid;");
            Console.WriteLine("    ** HALTON, from a Halton sequence;");
            Console.WriteLine("    ** RANDOM, using C++ RANDOM function;");
            Console.WriteLine("    ** UNIFORM, using a simple uniform RNG;");
            Console.WriteLine("    ** USER, call the \"user\" routine;");
            Console.WriteLine("  * CVT_IT_NUM, the maximum number of CVT iterations.");
            Console.WriteLine("  * SAMPLE, how to conduct the sampling.");
            Console.WriteLine("    ** GRID, by picking points from a grid;");
            Console.WriteLine("    ** HALTON, from a Halton sequence;");
            Console.WriteLine("    ** RANDOM, using C++ RANDOM function;");
            Console.WriteLine("    ** UNIFORM, using a simple uniform RNG;");
            Console.WriteLine("    ** USER, call the \"user\" routine;");
            Console.WriteLine("  * SAMPLE_NUM, the number of sample points.");
            Console.WriteLine("  * BATCH, number of sample points to create at one time.");
            Console.WriteLine("  * LAT_IT_NUM, the number of Latinizing iterations.");
            Console.WriteLine("  * OUTPUT, a file in which to store the data.");
            Console.WriteLine("");
            Console.WriteLine("  To indicate that no further computations are");
            Console.WriteLine("  desired, it is enough to input a nonsensical value,");
            Console.WriteLine("  such as -1.");

            Console.WriteLine("  *");
            Console.WriteLine(" *");
            Console.WriteLine("*  Ready to generate a new dataset:");
            Console.WriteLine(" *");
            Console.WriteLine("  *");

            Console.WriteLine("  Enter DIM_NUM, the spatial dimension:");
            Console.WriteLine("  (Try \"2\" if you do not have a preference.)");
            Console.WriteLine("  (0 or any negative value terminates execution).");

            try
            {
                dim_num = Convert.ToInt32(Console.ReadLine());
            }
            catch (Exception e)
            {
                Console.WriteLine("");
                Console.WriteLine("LCVT_DATASET - Warning!");
                Console.WriteLine("  Terminating abnormally because of an I/O error");
                Console.WriteLine("  while expecting input for DIM_NUM.");
                return;
            }

            Console.WriteLine("  User input DIM_NUM = " + dim_num + "");

            if (dim_num < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("LCVT_DATASET");
                Console.WriteLine("  The input value of DIM_NUM = " + dim_num + "");
                Console.WriteLine("  is interpreted as a request for termination.");
                Console.WriteLine("  Normal end of execution.");
                return;
            }

            Console.WriteLine("");
            Console.WriteLine("  Enter N, the number of points to generate:");
            Console.WriteLine("  (Try \"25\" if you do not have a preference.)");
            Console.WriteLine("  (0 or any negative value terminates execution).");

            try
            {
                n = Convert.ToInt32(Console.ReadLine());
            }
            catch (Exception e)
            {
                Console.WriteLine("");
                Console.WriteLine("LCVT_DATASET - Warning!");
                Console.WriteLine("  Terminating abnormally because of an I/O error");
                Console.WriteLine("  while expecting input for N.");
                return;
            }

            Console.WriteLine("  User input N = " + n + "");

            if (n < 1)
            {
                Console.WriteLine("");
                Console.WriteLine("LCVT_DATASET");
                Console.WriteLine("  The input value of N = " + n + "");
                Console.WriteLine("  is interpreted as a request for termination.");
                Console.WriteLine("  Normal end of execution.");
                return;
            }

            Console.WriteLine("");
            Console.WriteLine("  Enter SEED_INIT, a seed for the random number generator:");
            Console.WriteLine("  (Try \"123456789\" if you do not have a preference.)");
            Console.WriteLine("  (Any negative value terminates execution).");

            try
            {
                seed_init = Convert.ToInt32(Console.ReadLine());
            }
            catch (Exception e)
            {
                Console.WriteLine("");
                Console.WriteLine("LCVT_DATASET - Warning!");
                Console.WriteLine("  Terminating abnormally because of an I/O error");
                Console.WriteLine("  while expecting input for SEED_INIT.");
                return;
            }

            Console.WriteLine("  User input SEED_INIT = " + seed_init + "");

            if (seed_init < 0)
            {
                Console.WriteLine("");
                Console.WriteLine("LCVT_DATASET");
                Console.WriteLine("  The input value of SEED_INIT = " + seed_init + "");
                Console.WriteLine("  is interpreted as a request for termination.");
                Console.WriteLine("  Normal end of execution.");
                return;
            }

            seed = seed_init;

            Console.WriteLine("");
            Console.WriteLine("  INIT is the method of initializing the data:");
            Console.WriteLine("");
            Console.WriteLine("  file     read data from a file;");
            Console.WriteLine("  GRID     by picking points from a grid;");
            Console.WriteLine("  HALTON   from a Halton sequence;");
            Console.WriteLine("  RANDOM   using C++ RANDOM function;");
            Console.WriteLine("  UNIFORM  using a simple uniform RNG;");
            Console.WriteLine("  USER     call the \"user\" routine;");
            Console.WriteLine("");
            Console.WriteLine("  (Try \"RANDOM\" if you do not have a preference.)");
            Console.WriteLine("  (A blank value terminates execution).");
            Console.WriteLine("");
            Console.WriteLine("  Enter INIT:");

            try
            {
                init_string = Console.ReadLine();
            }
            catch
            {
                Console.WriteLine("");
                Console.WriteLine("LCVT_DATASET - Warning!");
                Console.WriteLine("  Terminating abnormally because of an I/O error");
                Console.WriteLine("  while expecting input for INIT");
                return;
            }

            if (typeMethods.s_eqi(init_string, "RANDOM"))
            {
                init = -1;
                Console.WriteLine("  User input INIT = \"RANDOM\"");
            }
            else if (typeMethods.s_eqi(init_string, "UNIFORM"))
            {
                init = 0;
                Console.WriteLine("  User input INIT = \"UNIFORM\"");
            }
            else if (typeMethods.s_eqi(init_string, "HALTON"))
            {
                init = 1;
                Console.WriteLine("  User input INIT = \"HALTON\".");
            }
            else if (typeMethods.s_eqi(init_string, "GRID"))
            {
                init = 2;
                Console.WriteLine("  User input INIT = \"GRID\".");
            }
            else if (typeMethods.s_eqi(init_string, "USER"))
            {
                init = 3;
                Console.WriteLine("  User input INIT = \"USER\".");
            }
            else if (0 < typeMethods.s_len_trim(init_string))
            {
                init = 4;
                Console.WriteLine("  User input INIT = FILE_NAME = \"" + init_string + "\".");
                input_file_name = init_string;
            }
            else
            {
                Console.WriteLine();
                Console.WriteLine("LCVT_DATASET");
                Console.WriteLine("  The input value of INIT ");
                Console.WriteLine("  is interpreted as a request for termination.");
                Console.WriteLine("  Normal end of execution.");
                return;
            }

            Console.WriteLine("");
            Console.WriteLine("  CVT_IT_NUM is the number of CVT iterations.");
            Console.WriteLine("");
            Console.WriteLine("  A CVT iteration carries out the following steps:");
            Console.WriteLine("  * the Voronoi region associated with each");
            Console.WriteLine("    generator is estimated by sampling;");
            Console.WriteLine("  * the centroid of each Voronoi region is estimated.");
            Console.WriteLine("  * the generator is replaced by the centroid.");
            Console.WriteLine("");
            Console.WriteLine("  If \"enough\" sampling points are used,");
            Console.WriteLine("  and \"enough\" iterations are taken, this process");
            Console.WriteLine("  will converge");
            Console.WriteLine("");
            Console.WriteLine("  (Try \"50\" if you do not have a preference.)");
            Console.WriteLine("  (A negative value terminates execution).");
            Console.WriteLine("");
            Console.WriteLine("  Enter CVT_IT_NUM:");

            try
            {
                cvt_it_num = Convert.ToInt32(Console.ReadLine());

            }
            catch
            {
                Console.WriteLine("");
                Console.WriteLine("LCVT_DATASET - Warning!");
                Console.WriteLine("  Terminating abnormally because of an I/O error");
                Console.WriteLine("  while expecting input for CVT_IT_NUM.");
                return;
            }

            Console.WriteLine("  User input CVT_IT_NUM = " + cvt_it_num + "");

            if (cvt_it_num < 0)
            {
                Console.WriteLine("");
                Console.WriteLine("LCVT_DATASET");
                Console.WriteLine("  The input value of CVT_IT_NUM = " + cvt_it_num + "");
                Console.WriteLine("  is interpreted as a request for termination.");
                Console.WriteLine("  Normal end of execution.");
                return;
            }

            Console.WriteLine("");
            Console.WriteLine("  SAMPLE is the method of sampling the region:");
            Console.WriteLine("");
            Console.WriteLine("  GRID     by picking points from a grid;");
            Console.WriteLine("  HALTON   from a Halton sequence;");
            Console.WriteLine("  RANDOM   using C++ RANDOM function;");
            Console.WriteLine("  UNIFORM  using a simple uniform RNG;");
            Console.WriteLine("  USER     call the \"user\" routine;");
            Console.WriteLine("");
            Console.WriteLine("  (Try \"RANDOM\" if you do not have a preference.)");
            Console.WriteLine("  (A blank value terminates execution).");
            Console.WriteLine("");
            Console.WriteLine("  Enter SAMPLE:");

            try
            {
                sample_string = Console.ReadLine();
            }
            catch
            {
                Console.WriteLine("");
                Console.WriteLine("LCVT_DATASET - Warning!");
                Console.WriteLine("  Terminating abnormally because of an I/O error");
                Console.WriteLine("  while expecting input for SAMPLE.");
                return;
            }

            if (typeMethods.s_eqi(sample_string, "RANDOM"))
            {
                Console.WriteLine("  User input INIT = \"RANDOM\".");
                sample = -1;
            }
            else if (typeMethods.s_eqi(sample_string, "UNIFORM"))
            {
                Console.WriteLine("  User input INIT = \"UNIFORM\".");
                sample = 0;
            }
            else if (typeMethods.s_eqi(sample_string, "HALTON"))
            {
                Console.WriteLine("  User input INIT = \"HALTON\".");
                sample = 1;
            }
            else if (typeMethods.s_eqi(sample_string, "GRID"))
            {
                Console.WriteLine("  User input INIT = \"GRID\".");
                sample = 2;
            }
            else if (typeMethods.s_eqi(sample_string, "USER"))
            {
                Console.WriteLine("  User input INIT = \"USER\".");
                sample = 3;
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("LCVT_DATASET");
                Console.WriteLine("  The input value of SAMPLE ");
                Console.WriteLine("  is interpreted as a request for termination.");
                Console.WriteLine("  Normal end of execution.");
                return;
            }

            Console.WriteLine("");
            Console.WriteLine("  SAMPLE_NUM is the number of sample points for CVT.");
            Console.WriteLine("");
            Console.WriteLine("  The Voronoi regions will be explored by generating");
            Console.WriteLine("  SAMPLE_NUM points.  For each sample point, the");
            Console.WriteLine("  nearest generator is found.  Using more points");
            Console.WriteLine("  gives a better estimate of these regions.");
            Console.WriteLine("");
            Console.WriteLine("  SAMPLE_NUM should be much larger than N, the");
            Console.WriteLine("  number of generators.");
            Console.WriteLine("");
            Console.WriteLine("  (Try \"10000\" if you do not have a preference.)");
            Console.WriteLine("  (A zero or negative value terminates execution.)");
            Console.WriteLine("");
            Console.WriteLine("  Enter SAMPLE_NUM:");

            try
            {
                sample_num = Convert.ToInt32(Console.ReadLine());
            }
            catch
            {
                Console.WriteLine("");
                Console.WriteLine("LCVT_DATASET - Warning!");
                Console.WriteLine("  Terminating abnormally because of an I/O error");
                Console.WriteLine("  while expecting input for SAMPLE_NUM.");
                return;
            }

            Console.WriteLine("  User input SAMPLE_NUM = " + sample_num + "");

            if (sample_num <= 0)
            {
                Console.WriteLine("");
                Console.WriteLine("LCVT_DATASET");
                Console.WriteLine("  The input value of SAMPLE_NUM = " + sample_num + "");
                Console.WriteLine("  is interpreted as a request for termination.");
                Console.WriteLine("  Normal end of execution.");
                return;
            }

            Console.WriteLine("");
            Console.WriteLine("  BATCH is the number of sample points to create");
            Console.WriteLine("  at one time.");
            Console.WriteLine("");
            Console.WriteLine("  BATCH should be between 1 and SAMPLE_NUM.");
            Console.WriteLine("");
            Console.WriteLine("  It is FASTER to set BATCH to SAMPLE_NUM;");
            Console.WriteLine("  setting BATCH to 1 requires the least memory.");
            Console.WriteLine("");
            Console.WriteLine("  (Try \"" + Math.Min(sample_num, 1000) +
                              "\" if you do not have a preference.)");
            Console.WriteLine("  (A zero or negative value terminates execution.)");
            Console.WriteLine("");
            Console.WriteLine("  Enter BATCH:");

            try
            {
                batch = Convert.ToInt32(Console.ReadLine());
            }
            catch
            {
                Console.WriteLine("");
                Console.WriteLine("LCVT_DATASET - Warning!");
                Console.WriteLine("  Terminating abnormally because of an I/O error");
                Console.WriteLine("  while expecting input for SAMPLE_NUM.");
                return;
            }

            Console.WriteLine("  User input BATCH = " + batch + "");

            if (batch <= 0)
            {
                Console.WriteLine("");
                Console.WriteLine("LCVT_DATASET");
                Console.WriteLine("  The input value of BATCH = " + batch + "");
                Console.WriteLine("  is interpreted as a request for termination.");
                Console.WriteLine("  Normal end of execution.");
                return;
            }

            Console.WriteLine("");
            Console.WriteLine("  LAT_IT_NUM is the number of Latinizing iterations.");
            Console.WriteLine("");
            Console.WriteLine("  Each step of the latinizing iteration begins");
            Console.WriteLine("  by carrying out CVT_IT_NUM steps of CVT iteration,");
            Console.WriteLine("  after which the data is \"latinized\".");
            Console.WriteLine("");
            Console.WriteLine("  Often, one latinizing step is enough.");
            Console.WriteLine("");
            Console.WriteLine("  In some cases, it may be worth while to carry");
            Console.WriteLine("  out several latinizing steps; that is, the");
            Console.WriteLine("  Latinized data is smoothed by another series");
            Console.WriteLine("  of CVT steps, then latinized, and so on.");
            Console.WriteLine("");
            Console.WriteLine("  (Try \"1\" if you do not have a preference.)");
            Console.WriteLine("  (A negative value terminates execution).");
            Console.WriteLine("");
            Console.WriteLine("  Enter LAT_IT_NUM:");

            try
            {
                lat_it_num = Convert.ToInt32(Console.ReadLine());
            }
            catch
            {
                Console.WriteLine("");
                Console.WriteLine("LCVT_DATASET - Warning!");
                Console.WriteLine("  Terminating abnormally because of an I/O error");
                Console.WriteLine("  while expecting input for LAT_IT_NUM.");
                return;
            }

            Console.WriteLine("  User input LAT_IT_NUM = " + lat_it_num + "");

            if (cvt_it_num < 0)
            {
                Console.WriteLine("");
                Console.WriteLine("LCVT_DATASET");
                Console.WriteLine("  The input value of LAT_IT_NUM = " + lat_it_num + "");
                Console.WriteLine("  is interpreted as a request for termination.");
                Console.WriteLine("  Normal end of execution.");
                return;
            }

            Console.WriteLine("");
            Console.WriteLine("  OUTPUT is the name of a file into which");
            Console.WriteLine("  the computed data may be stored.");
            Console.WriteLine("");
            Console.WriteLine("  (Try \"lcvt.txt\" if you do not have a preference.)");
            Console.WriteLine("  (A blank value terminates execution).");
            Console.WriteLine("");
            Console.WriteLine("  Enter OUTPUT:");

            try
            {
                output_file_name = Console.ReadLine();
            }
            catch
            {
                Console.WriteLine("");
                Console.WriteLine("LCVT_DATASET - Warning!");
                Console.WriteLine("  Terminating abnormally because of an I/O error");
                Console.WriteLine("  while expecting input for OUTPUT.");
                return;
            }

            Console.WriteLine("  User input OUTPUT = \"" + output_file_name + "\".");

            if (typeMethods.s_len_trim(output_file_name) <= 0)
            {
                Console.WriteLine("");
                Console.WriteLine("LCVT_DATASET");
                Console.WriteLine("  The input value of OUTPUT ");
                Console.WriteLine("  is interpreted as a request for termination.");
                Console.WriteLine("  Normal end of execution.");
                return;
            }

            //
            //  Initialize the data.
            //
            if (init == 4)
            {
                r = typeMethods.r8table_data_read(input_file_name, dim_num, n);
            }
            else
            {
                n_total = n;
                r = new double[dim_num * n];
                reset = true;

                Region.region_sampler(ref rdata, dim_num, n, n_total, ref r, init, reset, ref seed);
            }

            if (debug)
            {
                Console.WriteLine("");
                Console.WriteLine("  Latin IT      CVT Energy    Latin Energy");
                Console.WriteLine("");
            }

            for (lat_it = 1; lat_it <= lat_it_num; lat_it++)
            {
                if (debug)
                {
                    Console.WriteLine("");
                    Console.WriteLine("    CVT IT  Change");
                    Console.WriteLine("");
                }

                for (cvt_it = 1; cvt_it <= cvt_it_num; cvt_it++)
                {
                    LatinCentroidalVoronoi.cvt_iteration(ref ldata, dim_num, n, ref r, sample_num, sample, ref seed,
                        ref cvt_it_diff);

                    if (debug)
                    {
                        Console.WriteLine("  " + cvt_it.ToString().PadLeft(8)
                                               + "  " + cvt_it_diff.ToString().PadLeft(14) + "");
                    }

                }

                if (debug)
                {
                    Console.WriteLine("");
                }

                cvt_energy = Cluster.cluster_energy(ref rdata, dim_num, n, r, sample_num, sample, ref seed);

                typeMethods.r8mat_latinize(dim_num, n, ref r);

                lat_energy = Cluster.cluster_energy(ref rdata, dim_num, n, r, sample_num, sample, ref seed);

                Console.WriteLine("  " + lat_it.ToString().PadLeft(8)
                                       + "  " + cvt_energy.ToString().PadLeft(14)
                                       + "  " + lat_energy.ToString().PadLeft(14) + "");
            }

            //
            //  Write the data to a file.
            //
            LatinCentroidalVoronoi.lcvt_write(dim_num, n, seed_init, init, input_file_name, sample,
                sample_num, cvt_it_num, cvt_energy, lat_it_num, lat_energy,
                r, output_file_name);

            Console.WriteLine("");
            Console.WriteLine("  The data was written to the file \""
                              + output_file_name + "\".");


            Console.WriteLine("");
            Console.WriteLine("LCVT_DATASET:");
            Console.WriteLine("  Normal end of execution.");

            Console.WriteLine("");
        }
    }
}