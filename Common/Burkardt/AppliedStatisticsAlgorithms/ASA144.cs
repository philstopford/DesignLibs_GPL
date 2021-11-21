using Burkardt.Uniform;

namespace Burkardt.AppliedStatistics;

public static partial class Algorithms
{
    public class RContData
    {
        public int ntotal;
        public int[] nvect = new int[1];
        public int seed;
    }
    public static void rcont(ref RContData data, int nrow, int ncol, int[] nrowt, int[] ncolt, ref int[] nsubt,
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

        ifault = 0;

        switch (key)
        {
            case false:
            {
                //
                //  Set KEY for subsequent calls.
                //
                key = true;
                data.seed = 123456789;
                switch (nrow)
                {
                    //
                    //  Check for faults and prepare for future calls.
                    //
                    case <= 0:
                        ifault = 1;
                        return;
                }

                switch (ncol)
                {
                    case <= 1:
                        ifault = 2;
                        return;
                }

                for (int i = 0; i < nrow; i++)
                {
                    switch (nrowt[i])
                    {
                        case <= 0:
                            ifault = 3;
                            return;
                    }
                }

                switch (ncolt[0])
                {
                    case <= 0:
                        ifault = 4;
                        return;
                }

                nsubt[0] = ncolt[0];

                for (int j = 1; j < ncol; j++)
                {
                    switch (ncolt[j])
                    {
                        case <= 0:
                            ifault = 4;
                            return;
                        default:
                            nsubt[j] = nsubt[j - 1] + ncolt[j];
                            break;
                    }
                }

                data.ntotal = nsubt[ncol - 1];

                data.nvect = new int[data.ntotal];
                //
                //  Initialize vector to be permuted.
                //
                for (int i = 0; i < data.ntotal; i++)
                {
                    data.nvect[i] = i + 1;
                }

                break;
            }
        }

        //
        //  Initialize vector to be permuted.
        //
        int[] nnvect = new int[data.ntotal];

        for (int i = 0; i < data.ntotal; i++)
        {
            nnvect[i] = data.nvect[i];
        }

        //
        //  Permute vector.
        //
        int ntemp = data.ntotal;

        for (int i = 0; i < data.ntotal; i++)
        {
            int noct = (int) (UniformRNG.r8_uniform_01(ref data.seed) * ntemp + 1.0);
            data.nvect[i] = nnvect[noct - 1];
            nnvect[noct - 1] = nnvect[ntemp - 1];
            ntemp -= 1;
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
                    if (data.nvect[ii] > nsubt[j])
                    {
                        continue;
                    }

                    ii += 1;
                    matrix[i + j * nrow] += 1;
                    break;
                }
            }
        }
    }
}