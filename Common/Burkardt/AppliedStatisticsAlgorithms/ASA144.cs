﻿using Burkardt.Uniform;

namespace Burkardt.AppliedStatistics
{
    public static partial class Algorithms
    {
        public static void rcont(int nrow, int ncol, int[] nrowt, int[] ncolt, ref int[] nsubt,
                            ref int[] matrix, ref bool key, ref int ifault )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RCONT generates a random two-way table with given marginal totals.
        //
        //  Discussion:
        //
        //    Each time the program is called, another table will be randomly
        //    generated.
        //
        //    Note that it should be the case that the sum of the row totals
        //    is equal to the sum of the column totals.  However, this program
        //    does not check for that condition.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    28 January 2008
        //
        //  Author:
        //
        //    Original FORTRAN77 version by James Boyett.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    James Boyett,
        //    Algorithm AS 144:
        //    Random R x C Tables with Given Row and Column Totals,
        //    Applied Statistics,
        //    Volume 28, Number 3, pages 329-332, 1979.
        //
        //  Parameters:
        //
        //    Input, int NROW, the number of rows in the observed matrix.
        //
        //    Input, int NCOL, the number of columns in the observed matrix.
        //
        //    Input, int NROWT[NROW], the row totals of the observed matrix.
        //
        //    Input, int NCOLT[NCOL], the column totals of the observed matrix.
        //
        //    Input/output, int NSUBT[NCOL], used by RCONT for partial column sums.
        //    Must not be changed by the calling program.
        //
        //    Output, int MATRIX[NROW*NCOL], the randomly generated matrix.
        //
        //    Input/output, bool *KEY, should be set to FALSE by the user before
        //    the initial call.  RCONT will reset it to TRUE, and it should be left
        //    at that value for subsequent calls in which the same values of NROW,
        //    NCOL, NROWT and NCOLT are being used.
        //
        //    Output, int *IFAULT, fault indicator.
        //    0, no error occured.
        //    1, NROW <= 0.
        //    2, NCOL <= 1.
        //    3, some entry of NROWT is less than 0.
        //    4, some entry of NCOLT is less than 0.
        //
        {
            int ntotal = 0;
            int[] nvect = new int[1];
            int seed =  0;

            ifault = 0;

            if (!(key))
            {
                //
                //  Set KEY for subsequent calls.
                //
                key = true;
                seed = 123456789;
                //
                //  Check for faults and prepare for future calls.
                //
                if (nrow <= 0)
                {
                    ifault = 1;
                    return;
                }

                if (ncol <= 1)
                {
                    ifault = 2;
                    return;
                }

                for (int i = 0; i < nrow; i++)
                {
                    if (nrowt[i] <= 0)
                    {
                        ifault = 3;
                        return;
                    }
                }

                if (ncolt[0] <= 0)
                {
                    ifault = 4;
                    return;
                }

                nsubt[0] = ncolt[0];

                for (int j = 1; j < ncol; j++)
                {
                    if (ncolt[j] <= 0)
                    {
                        ifault = 4;
                        return;
                    }

                    nsubt[j] = nsubt[j - 1] + ncolt[j];
                }

                ntotal = nsubt[ncol - 1];

                nvect = new int[ntotal];
                //
                //  Initialize vector to be permuted.
                //
                for (int i = 0; i < ntotal; i++)
                {
                    nvect[i] = i + 1;
                }
            }

            //
            //  Initialize vector to be permuted.
            //
            int[] nnvect = new int[ntotal];

            for (int i = 0; i < ntotal; i++)
            {
                nnvect[i] = nvect[i];
            }

            //
            //  Permute vector.
            //
            int ntemp = ntotal;

            for (int i = 0; i < ntotal; i++)
            {
                int noct = (int) (UniformRNG.r8_uniform_01(ref seed) * (double) (ntemp) + 1.0);
                nvect[i] = nnvect[noct - 1];
                nnvect[noct - 1] = nnvect[ntemp - 1];
                ntemp = ntemp - 1;
            }

            //
            //  Construct random matrix.
            //
            for (int j = 0; j < ncol; j++)
            {
                for (int i = 0; i < nrow; i++)
                {
                    matrix[i + j * nrow] = 0;
                }
            }

            int ii = 0;

            for (int i = 0; i < nrow; i++)
            {
                int limit = nrowt[i];

                for (int k = 0; k < limit; k++)
                {
                    for (int j = 0; j < ncol; j++)
                    {
                        if (nvect[ii] <= nsubt[j])
                        {
                            ii = ii + 1;
                            matrix[i + j * nrow] = matrix[i + j * nrow] + 1;
                            break;
                        }
                    }
                }
            }
        }
    }
}