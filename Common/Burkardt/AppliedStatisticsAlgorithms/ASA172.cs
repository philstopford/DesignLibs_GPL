﻿using System;

namespace Burkardt.AppliedStatistics;

public static partial class Algorithms
{
    public static void revers(ref int[] ivec, int kdim )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    REVERS reorders the subscript vector, if required.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    27 July 2008
        //
        //  Author:
        //
        //    Original FORTRAN77 version by M O'Flaherty, G MacKenzie.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    M O'Flaherty, G MacKenzie,
        //    Algorithm AS 172:
        //    Direct Simulation of Nested Fortran DO-LOOPS,
        //    Applied Statistics,
        //    Volume 31, Number 1, 1982, pages 71-74.
        //
        //  Parameters:
        //
        //    Input/output, int IVEC[KDIM], the subscript vector.
        //
        //    Input, int KDIM, the dimension of the subscript vector.
        //
    {
        for (int i = 0; i < kdim / 2; i++)
        {
            (ivec[i], ivec[kdim - 1 - i]) = (ivec[kdim - 1 - i], ivec[i]);
        }
    }

    public static int simdo(bool qind, bool qfor, int[] iprod, int kdim, ref int jsub, ref int[] ivec )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SIMDO generates multi-indices, simulating nested DO-loops.
        //
        //  Discussion:
        //
        //    The loops are assumed to be nested to a depth of K.
        //
        //    The R-th loop is assumed to have upper limit N(R) and increment Inc(R).
        //
        //    The total number of executions of the innermost loop is 
        //
        //      N = product ( 1 <= R <= K ) N(R).
        //
        //    Let these executions be indexed by the single integer J, which
        //    we call the index subscript.
        //
        //    Each value of J corresponds to a particular set of loop indices,
        //    which we call the subscript vector I(J).
        //
        //    This routine can start with J and find I(J), or determine
        //    J from I(J).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        // 
        //  Modified:
        //
        //    27 July 2008
        //
        //  Author:
        //
        //    Original FORTRAN77 version by M O'Flaherty, G MacKenzie.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    M O'Flaherty, G MacKenzie,
        //    Algorithm AS 172:
        //    Direct Simulation of Nested Fortran DO-LOOPS,
        //    Applied Statistics,
        //    Volume 31, Number 1, 1982, pages 71-74.
        //
        //  Parameters:
        //
        //    Input, bool QIND.
        //    TRUE to convert an index subscript J to the subscript vector I(J).
        //    FALSE to convert the subscript vector I(J) to the index subscript J.
        //
        //    Input, bool QFOR,
        //    TRUE if conversion is required in standard Fortran subscripting order,
        //    FALSE otherwise.
        //
        //    Input, int IPROD[KDIM], contains the partial products.
        //    If QFOR is FALSE, then
        //      IPROD(S) = product ( 1 <= R <= S ) N(R).
        //    If QFOR is TRUE, then
        //      IPROD(S) = product ( 1 <= R <= S ) N(KDIM+1-R).
        //
        //    Input, int KDIM, the nesting depth of the loops.
        //
        //    Input/output, int *JSUB.
        //    If QIND is TRUE, then JSUB is an input quantity, an index subscript J
        //    to be converted into the subscript vector I(J).
        //    If QIND is FALSE, then JSUB is an output quantity, the index subscript J
        //    corresponding to the subscript vector I(J).
        //
        //    Input/output, int IVEC[KDIM].
        //    if QIND is TRUE, then IVEC is an output quantity, the subscript vector I(J)
        //    corresponding to the index subscript J.
        //    If QIND is FALSE, then IVEC is an input quantity, a subscript vector I(J)
        //    for which the corresponding index subscript J is to be computed.
        //
        //    Output, int SIMDO, error flag.
        //    0, no error was detected.
        //    1, if QIND is TRUE, and the input value of JSUB exceeds IPROD(KDIM).
        //    2, if QIND is FALSE, and IVEC contains an illegal component.
        //
    {
        int i;

        int ifault = 0;
        switch (qind)
        {
            //
            //  Index subscript to subscript vector conversion.
            //
            case true when iprod[kdim - 1] < jsub:
                ifault = 1;
                Console.WriteLine("");
                Console.WriteLine("SIMDO - Fatal error!");
                Console.WriteLine("  JSUB is out of bounds.");
                return ifault;
            case true:
            {
                int itempv = jsub - 1;

                for (i = 0; i < kdim - 1; i++)
                {
                    int ik = kdim - 2 - i;
                    ivec[i] = itempv / iprod[ik];
                    itempv -= iprod[ik] * ivec[i];
                    ivec[i] += 1;
                }

                ivec[kdim - 1] = itempv + 1;
                switch (qfor)
                {
                    case true:
                        revers(ref ivec, kdim);
                        break;
                }

                break;
            }
            //
            default:
            {
                switch (qfor)
                {
                    case false:
                        revers(ref ivec, kdim);
                        break;
                }

                if (iprod[0] < ivec[0])
                {
                    ifault = 2;
                    Console.WriteLine("");
                    Console.WriteLine("SIMDO - Fatal error!");
                    Console.WriteLine("  An entry of IVEC is out of bounds.");
                    return ifault;
                }

                for (i = 1; i < kdim; i++)
                {
                    if (iprod[i] / iprod[i - 1] >= ivec[i])
                    {
                        continue;
                    }

                    ifault = 2;
                    Console.WriteLine("");
                    Console.WriteLine("SIMDO - Fatal error!");
                    Console.WriteLine("  An entry of IVEC is out of bounds.");
                    return ifault;
                }

                jsub = ivec[0];
                for (i = 1; i < kdim; i++)
                {
                    jsub += (ivec[i] - 1) * iprod[i - 1];
                }

                switch (qfor)
                {
                    //
                    //  As a courtesy to the caller, UNREVERSE the IVEC vector
                    //  if you reversed it.
                    //
                    case false:
                        revers(ref ivec, kdim);
                        break;
                }

                break;
            }
        }

        return ifault;
    }
}