using System;
using Burkardt.Function;
using Burkardt.Sequence;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.Sampling;

public class RegionData
{
    public BTupleData tdata = new();
    public int[] halton_base;
    public int halton_seed = 1;
    public int ngrid;
    public int rank;
    public int[] tuple;
}
    
public static class Region
{
    public static void region_sampler(ref RegionData data, int m, int n, int n_total, ref double[] x,
            int sample_function, bool reset, ref int seed )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    REGION_SAMPLER returns a sample point in the physical region.
        //
        //  Discussion:
        //
        //    This routine original interfaced with a lower routine called
        //    TEST_REGION, which tested whether the points generated in the
        //    bounding box were actually inside a possibly smaller physical
        //    region of interest.  It's been a long time since that option
        //    was actually used, so it's been dropped.
        //
        //    A point is chosen in the bounding box, either by a uniform random
        //    number generator, or from a vector Halton sequence.
        //
        //    The entries of the local vector HALTON_BASE should be distinct primes.
        //    Right now, we're assuming M is no greater than 3.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    04 March 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the spatial dimension.
        //
        //    Input, int N, the number of points to generate now.
        //
        //    Input, int N_TOTAL, the total number of points to generate.
        //
        //    Output, double X[M*N], the sample points.
        //
        //    Input, int SAMPLE_FUNCTION, region sampling function:
        //    -1, sampling function is RANDOM (C++ STDLIB library function);
        //    0, sampling function is UNIFORM;
        //    1, sampling function is HALTON;
        //    2, sampling function is GRID;
        //    3, sample points are generated elsewhere, and this routine is skipped.
        //
        //    Input, bool RESET, if true, then this is the first call for a particular
        //    calculation, and initialization should be taken care of.
        //
        //    Input/output, int *SEED, the random number seed.
        //
    {
        double exponent;
        int i;
        int j;
        int k;
        switch (sample_function)
        {
            //
            case -1:
            {
                for (k = 0; k < m * n; k++)
                {
                    x[k] = entropyRNG.RNG.nextdouble();
                }

                break;
            }
            case 0:
            {
                for (k = 0; k < m * n; k++)
                {
                    x[k] = UniformRNG.r8_uniform_01(ref seed);
                }

                break;
            }
            case 1:
            {
                switch (reset)
                {
                    case true:
                    {
                        data.halton_seed = 1;
                        reset = false;

                        data.halton_base = new int[m];
                        for (i = 0; i < m; i++)
                        {
                            data.halton_base[i] = Prime.prime(i + 1);
                        }

                        break;
                    }
                }

                //
                //  The unusual syntax X+J*M essentially means pass the address of the beginning
                //  of the J-th vector of length M in X.
                //
                for (j = 0; j < n; j++)
                {
                    Halton.i4_to_halton(data.halton_seed, data.halton_base, m, ref x, rIndex: + j * m);
                    data.halton_seed += 1;
                }

                break;
            }
            case 2:
            {
                switch (reset)
                {
                    case true:
                    {
                        data.rank = 0;
                        exponent = 1.0 / m;

                        data.ngrid = (int) Math.Pow(n_total, exponent);

                        if (Math.Pow(data.ngrid, m) < n_total)
                        {
                            data.ngrid += 1;
                        }
                    
                        data.tuple = new int[m];
                        reset = false;
                        break;
                    }
                }

                for (j = 0; j < n; j++)
                {
                    BTuple.tuple_next_fast(ref data.tdata, data.ngrid, m, data.rank, ref data.tuple);
                    data.rank += 1;
                    for (i = 0; i < m; i++)
                    {
                        x[j * m + i] = (2 * data.tuple[i] - 1)
                                       / (double) (2 * data.ngrid);
                    }
                }

                break;
            }
            case 3:
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("REGION_SAMPLER - Fatal error!");
                Console.WriteLine("  Illegal SAMPLE_FUNCTION = " + sample_function + "");
                break;
        }
    }
}