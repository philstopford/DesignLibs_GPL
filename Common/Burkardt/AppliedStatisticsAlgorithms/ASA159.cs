﻿using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.AppliedStatistics
{
    public static partial class Algorithms
    {

        public static void rcont2(int nrow, int ncol, int[] nrowt, int[] ncolt, ref bool key,
        ref int seed, ref int[] matrix, ref int ierror )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RCONT2 constructs a random two-way contingency table with given sums.
        //
        //  Discussion:
        //
        //    It is possible to specify row and column sum vectors which
        //    correspond to no table at all.  As far as I can see, this routine does
        //    not detect such a case.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    10 March 2009
        //
        //  Author:
        //
        //    Original FORTRAN77 version by WM Patefield.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    WM Patefield,
        //    Algorithm AS 159:
        //    An Efficient Method of Generating RXC Tables with
        //    Given Row and Column Totals,
        //    Applied Statistics,
        //    Volume 30, Number 1, 1981, pages 91-97.
        //
        //  Parameters:
        //
        //    Input, int NROW, NCOL, the number of rows and columns 
        //    in the table.  NROW and NCOL must each be at least 2.
        //
        //    Input, int NROWT[NROW], NCOLT[NCOL], the row and column 
        //    sums.  Each entry must be positive.
        //
        //    Input/output, bool *KEY, a flag that indicates whether data has
        //    been initialized for this problem.  Set KEY = .FALSE. before the first
        //    call.
        //
        //    Input/output, int *SEED, a seed for the random number generator.
        //
        //    Output, int MATRIX[NROW*NCOL], the matrix.
        //
        //    Output, int *IERROR, an error flag, which is returned 
        //    as 0 if no error occurred.
        //
        {
            bool done1 = false;
            bool done2 = false;
            double[] fact = { };
            int i;
            int ia;
            int iap;
            int ib = 0;
            int ic;
            int id;
            int idp;
            int ie;
            int igp;
            int ihp;
            int ii;
            int iip;
            int j;
            int jc;
            int[] jwork;
            int l;
            bool lsm;
            bool lsp;
            int m;
            int nll;
            int nlm = 0;
            int nlmp;
            int nrowtl;
            int ntotal = 0;
            double r;
            double sumprb;
            double x;
            double y;

           ierror = 0;
//
//  On user's signal, set up the factorial table.
//
            if (!(key))
            {

                key = true;

                if (nrow <= 1)
                {
                    Console.WriteLine("");
                    Console.WriteLine("RCONT - Fatal error!");
                    Console.WriteLine("  Input number of rows is less than 2.");
                    ierror = 1;
                    return;
                }

                if (ncol <= 1)
                {
                    Console.WriteLine("");
                    Console.WriteLine("RCONT - Fatal error!");
                    Console.WriteLine("  The number of columns is less than 2.");
                    ierror = 2;
                    return;
                }

                for (i = 0; i < nrow; i++)
                {
                    if (nrowt[i] <= 0)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("RCONT - Fatal error!");
                        Console.WriteLine("  An entry in the row sum vector is not positive.");
                        ierror = 3;
                        return;
                    }
                }

                for (j = 0; j < ncol; j++)
                {
                    if (ncolt[j] <= 0)
                    {
                        Console.WriteLine("");
                        Console.WriteLine("RCONT - Fatal error!");
                        Console.WriteLine("  An entry in the column sum vector is not positive.");
                        ierror = 4;
                        return;
                    }
                }

                if (typeMethods.i4vec_sum(ncol, ncolt) != typeMethods.i4vec_sum(nrow, nrowt))
                {
                    Console.WriteLine("");
                    Console.WriteLine("RCONT - Fatal error!");
                    Console.WriteLine("  The row and column sum vectors do not have the same sum.");
                    ierror = 6;
                    return;
                }

                ntotal = typeMethods.i4vec_sum(ncol, ncolt);

                fact = new double[ntotal + 1];
//
//  Calculate log-factorials.
//
                x = 0.0;
                fact[0] = 0.0;
                for (i = 1; i <= ntotal; i++)
                {
                    x = x + Math.Log((double) (i));
                    fact[i] = x;
                }

            }

//
//  Construct a random matrix.
//
            jwork = new int[ncol];

            for (i = 0; i < ncol - 1; i++)
            {
                jwork[i] = ncolt[i];
            }

            jc = ntotal;

            for (l = 0; l < nrow - 1; l++)
            {
                nrowtl = nrowt[l];
                ia = nrowtl;
                ic = jc;
                jc = jc - nrowtl;

                for (m = 0; m < ncol - 1; m++)
                {
                    id = jwork[m];
                    ie = ic;
                    ic = ic - id;
                    ib = ie - ia;
                    ii = ib - id;
//
//  Test for zero entries in matrix.
//
                    if (ie == 0)
                    {
                        ia = 0;
                        for (j = m; j < ncol; j++)
                        {
                            matrix[l + j * nrow] = 0;
                        }

                        break;
                    }

//
//  Generate a pseudo-random number.
//
                    r = UniformRNG.r8_uniform_01(ref seed);
//
//  Compute the conditional expected value of MATRIX(L,M).
//
                    done1 = false;

                    for (;;)
                    {
                        if (fact.Length == 0)
                        {
                            break;
                        }
                        nlm = (int) ((double) (ia * id) / (double) (ie) + 0.5);
                        iap = ia + 1;
                        idp = id + 1;
                        igp = idp - nlm;
                        ihp = iap - nlm;
                        nlmp = nlm + 1;
                        iip = ii + nlmp;
                        x = Math.Exp(fact[iap - 1] + fact[ib] + fact[ic] + fact[idp - 1] -
                                fact[ie] - fact[nlmp - 1] - fact[igp - 1] - fact[ihp - 1] - fact[iip - 1]);

                        if (r <= x)
                        {
                            break;
                        }

                        sumprb = x;
                        y = x;
                        nll = nlm;
                        lsp = false;
                        lsm = false;
//
//  Increment entry in row L, column M.
//
                        while (!lsp)
                        {
                            j = (id - nlm) * (ia - nlm);

                            if (j == 0)
                            {
                                lsp = true;
                            }
                            else
                            {
                                nlm = nlm + 1;
                                x = x * (double) (j) / (double) (nlm * (ii + nlm));
                                sumprb = sumprb + x;

                                if (r <= sumprb)
                                {
                                    done1 = true;
                                    break;
                                }
                            }

                            done2 = false;

                            while (!lsm)
                            {
//
//  Decrement the entry in row L, column M.
//
                                j = nll * (ii + nll);

                                if (j == 0)
                                {
                                    lsm = true;
                                    break;
                                }

                                nll = nll - 1;
                                y = y * (double) (j) / (double) ((id - nll) * (ia - nll));
                                sumprb = sumprb + y;

                                if (r <= sumprb)
                                {
                                    nlm = nll;
                                    done2 = true;
                                    break;
                                }

                                if (!lsp)
                                {
                                    break;
                                }

                            }

                            if (done2)
                            {
                                break;
                            }

                        }

                        if (done1)
                        {
                            break;
                        }

                        if (done2)
                        {
                            break;
                        }

                        r = UniformRNG.r8_uniform_01(ref seed);
                        r = sumprb * r;

                    }

                    matrix[l + m * nrow] = nlm;
                    ia = ia - nlm;
                    jwork[m] = jwork[m] - nlm;

                }

                matrix[l + (ncol - 1) * nrow] = ia;
            }

//
//  Compute the last row.
//
            for (j = 0; j < ncol - 1; j++)
            {
                matrix[nrow - 1 + j * nrow] = jwork[j];
            }

            matrix[nrow - 1 + (ncol - 1) * nrow] = ib - matrix[nrow - 1 + (ncol - 2) * nrow];

        }

    }
}