using System;
using Burkardt.Function;
using Burkardt.PolynomialNS;
using Burkardt.Types;

namespace Burkardt.RandomNS
{
    public class FaureData
    {
        public int[] coef { get; set; }
        public int[] ytemp { get; set; }
        public int qs { get; set; }
        public int hisum_save { get; set; }

        public FaureData()
        {
            qs = -1;
            hisum_save = -1;
        }
    }

    public static class Faure
    {
        public static void faure(ref FaureData data, int dim_num, ref int seed, ref double[] quasi, int quasiIndex = 0)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FAURE generates a new quasirandom Faure vector with each call.
            //
            //  Discussion:
            //
            //    This routine implements the Faure method for computing
            //    quasirandom numbers.  It is a merging and adaptation of
            //    Bennett Fox's routines INFAUR and GOFAUR from ACM TOMS Algorithm 647.
            //
            //    Michael Baudin suggested a change so that whenever HISUM is altered,
            //    the binomial table is recomputed, which allows a user to go back or
            //    forth in the sequence.  08 December 2009.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    08 December 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Henri Faure,
            //    Discrepance de suites associees a un systeme de numeration
            //    (en dimension s),
            //    Acta Arithmetica,
            //    Volume 41, 1982, pages 337-351.
            //
            //    Bennett Fox,
            //    Algorithm 647:
            //    Implementation and Relative Efficiency of Quasirandom 
            //    Sequence Generators,
            //    ACM Transactions on Mathematical Software,
            //    Volume 12, Number 4, December 1986, pages 362-376.
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the spatial dimension, which should be
            //    at least 2.
            //
            //    Input/output, int *SEED, the seed, which can be used to index
            //    the values.  On first call, set the input value of SEED to be 0
            //    or negative.  The routine will automatically initialize data,
            //    and set SEED to a new value.  Thereafter, to compute successive
            //    entries of the sequence, simply call again without changing
            //    SEED.  On the first call, if SEED is negative, it will be set
            //    to a positive value that "skips over" an early part of the sequence
            //    (This is recommended for better results).
            //
            //    Output, double QUASI[DIM_NUM], the next quasirandom vector.
            //
        {
            int hisum;
            int i;
            int j;
            int k;
            int ktemp;
            int ltemp;
            int mtemp;
            double r;
            int ztemp;
            //
            //  Initialization required or requested?
            //
            if (data.qs <= 0 || seed <= 0)
            {
                data.qs = Prime.prime_ge(dim_num);

                if (data.qs < 1)
                {
                    Console.WriteLine("");
                    Console.WriteLine("FAURE - Fatal error!");
                    Console.WriteLine("  PRIME_GE failed.");
                    return;
                }

                data.hisum_save = -1;
            }

            //
            //  If SEED < 0, reset for recommended initial skip.
            //
            if (seed < 0)
            {
                hisum = 3;
                seed = (int)Math.Pow(data.qs, hisum + 1) - 1;
            }
            else if (seed == 0)
            {
                hisum = 0;
            }
            else
            {
                hisum = (int)Math.Log(seed, data.qs);
            }

            //
            //  Is it necessary to recompute the coefficient table?
            //
            if (data.hisum_save != hisum)
            {
                data.hisum_save = hisum;

                data.coef = Coefficients.binomial_table(data.qs, hisum, hisum);

                data.ytemp = new int[hisum + 1];
            }

            //
            //  Find QUASI(1) using the method of Faure.
            //
            //  SEED has a representation in base QS of the form: 
            //
            //    Sum ( 0 <= J <= HISUM ) YTEMP(J) * QS**J
            //
            //  We now compute the YTEMP(J)'s.
            //
            ktemp = (int)Math.Pow(data.qs, hisum + 1);
            ltemp = seed;

            for (i = hisum; 0 <= i; i--)
            {
                ktemp = ktemp / data.qs;
                mtemp = ltemp % ktemp;
                data.ytemp[i] = (ltemp - mtemp) / ktemp;
                ltemp = mtemp;
            }

            //
            //  QUASI(K) has the form
            //
            //    Sum ( 0 <= J <= HISUM ) YTEMP(J) / QS**(J+1)
            //
            //  Compute QUASI(1) using nested multiplication.
            //
            r = ((double) data.ytemp[hisum]);
            for (i = hisum - 1; 0 <= i; i--)
            {
                r = ((double) data.ytemp[i]) + r / ((double) data.qs);
            }

            quasi[quasiIndex + 0] = r / ((double) data.qs);
            //
            //  Find components QUASI(2:DIM_NUM) using the Faure method.
            //
            for (k = 1; k < dim_num; k++)
            {
                quasi[quasiIndex + k] = 0.0;
                r = 1.0 / ((double) data.qs);

                for (j = 0; j <= hisum; j++)
                {
                    ztemp = 0;
                    for (i = j; i <= hisum; i++)
                    {
                        ztemp = ztemp + data.ytemp[i] * data.coef[i + j * (hisum + 1)];
                    }

                    //
                    //  New YTEMP(J) is:
                    //
                    //    Sum ( J <= I <= HISUM ) ( old ytemp(i) * binom(i,j) ) mod QS.
                    //
                    data.ytemp[j] = ztemp % data.qs;
                    quasi[quasiIndex + k] = quasi[quasiIndex + k] + ((double) data.ytemp[j]) * r;
                    r = r / ((double) data.qs);
                }
            }

            //
            //  Update SEED.
            //
            seed = seed + 1;
        }

        public static double[] faure_generate(int dim_num, int n, int skip)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    FAURE_GENERATE generates a Faure dataset.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    11 December 2009
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int DIM_NUM, the spatial dimension.
            //
            //    Input, int N, the number of points to generate.
            //
            //    Input, int SKIP, the number of initial points to skip.
            //
            //    Output, double FAURE_GENERATE[DIM_NUM*N], the points.
            //
        {
            int j;
            double[] r;
            int seed;

            r = new double[dim_num * n];

            seed = skip;
            FaureData data = new FaureData();
            for (j = 0; j < n; j++)
            {
                faure(ref data, dim_num, ref seed, ref r, + j * dim_num);
            }

            return r;
        }
    }
}