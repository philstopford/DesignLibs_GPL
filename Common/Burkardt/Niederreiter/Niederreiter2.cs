﻿using System;

namespace Burkardt.NiederreiterNS;

public static class Niederreiter2
{
    private const int MAXDEG = 50;
    private const int DIM_MAX = 20;
    private const int NBITS = 31;

    public static void calcc2(ref NiederReiter2CalcData data, int dim_num, ref int[,] cj)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CALCC2 computes values of the constants C(I,J,R).
        //
        //  Discussion:
        //
        //    This program calculates the values of the constants C(I,J,R).
        //
        //    As far as possible, Niederreiter's notation is used.
        //
        //    For each value of I, we first calculate all the corresponding
        //    values of C.  These are held in the array CI.  All these
        //    values are either 0 or 1.  
        //
        //    Next we pack the values into the
        //    array CJ, in such a way that CJ(I,R) holds the values of C
        //    for the indicated values of I and R and for every value of
        //    J from 1 to NBITS.  The most significant bit of CJ(I,R)
        //    (not counting the sign bit) is C(I,1,R) and the least
        //    significant bit is C(I,NBITS,R).
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 March 2003
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Paul Bratley, Bennett Fox, Harald Niederreiter.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    R Lidl, Harald Niederreiter, 
        //    Finite Fields,
        //    Cambridge University Press, 1984, page 553.
        //
        //    Harald Niederreiter,
        //    Low-discrepancy and low-dispersion sequences,
        //    Journal of Number Theory,
        //    Volume 30, 1988, pages 51-70.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the dimension of the sequence to be generated.
        //
        //    Output, int CJ[DIM_MAX][NBITS], the packed values of 
        //    Niederreiter's C(I,J,R)
        //
        //  Local Parameters:
        //
        //    Local, int MAXE; we need DIM_MAX irreducible polynomials over Z2.
        //    MAXE is the highest degree among these.
        //
        //    Local, int MAXV, the maximum possible index used in V.
        //
    {
        const int MAXE = 6;

        int[,] add = new int [2,2];
        int[] b = new int[MAXDEG + 1];
        int[,] ci = new int[NBITS,NBITS];
        int i;
        int[,] irred =
            {
                {
                    0,1,0,0,0,0,0
                },
                {
                    1,1,0,0,0,0,0
                },
                {
                    1,1,1,0,0,0,0
                },
                {
                    1,1,0,1,0,0,0
                },
                {
                    1,0,1,1,0,0,0
                },
                {
                    1,1,0,0,1,0,0
                },
                {
                    1,0,0,1,1,0,0
                },
                {
                    1,1,1,1,1,0,0
                },
                {
                    1,0,1,0,0,1,0
                },
                {
                    1,0,0,1,0,1,0
                },
                {
                    1,1,1,1,0,1,0
                },
                {
                    1,1,1,0,1,1,0
                },
                {
                    1,1,0,1,1,1,0
                },
                {
                    1,0,1,1,1,1,0
                },
                {
                    1,1,0,0,0,0,1
                },
                {
                    1,0,0,1,0,0,1
                },
                {
                    1,1,1,0,1,0,1
                },
                {
                    1,1,0,1,1,0,1
                },
                {
                    1,0,0,0,0,1,1
                },
                {
                    1,1,1,0,0,1,1
                }
            }
            ;
        int[] irred_deg = 
            {
                1, 1, 2, 3, 3, 4, 4, 4, 5, 5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6
            }
            ;
        int maxv = NBITS + MAXE;
        int[,] mul = new int[2,2];
        int[] px = new int[MAXDEG + 1];
        int[,] sub = new int[2,2];
        int[] v = new int[NBITS + MAXE + 1];
        //
        //  Prepare to work in Z2.
        //
        setfld2(add, mul, sub);

        for (i = 0; i < dim_num; i++)
        {
            //
            //  For each dimension, we need to calculate powers of an
            //  appropriate irreducible polynomial:  see Niederreiter
            //  page 65, just below equation (19).
            //
            //  Copy the appropriate irreducible polynomial into PX,
            //  and its degree into E.  Set polynomial B = PX ** 0 = 1.
            //  M is the degree of B.  Subsequently B will hold higher
            //  powers of PX.
            //
            int e = irred_deg[i];

            int px_deg = irred_deg[i];

            int j;
            for (j = 0; j <= px_deg; j++)
            {
                px[j] = irred[i,j];
            }

            int b_deg = 0;
            b[0] = 1;
            //
            //  Niederreiter (page 56, after equation (7), defines two
            //  variables Q and U.  We do not need Q explicitly, but we do need U.
            //
            int u = 0;

            int r;
            for (j = 0; j < NBITS; j++)
            {
                switch (u)
                {
                    //
                    //  If U = 0, we need to set B to the next power of PX
                    //  and recalculate V.  This is done by subroutine CALCV.
                    //
                    case 0:
                        calcv2(ref data, maxv, px_deg, px, add, mul, sub, ref b_deg, b, v);
                        break;
                }

                //
                //  Now C is obtained from V.  Niederreiter obtains A from V (page 65, 
                //  near the bottom), and then gets C from A (page 56, equation (7)).  
                //  However this can be done in one step.  Here CI(J,R) corresponds to
                //  Niederreiter's C(I,J,R).
                //
                for (r = 0; r < NBITS; r++)
                {
                    ci[j,r] = v[r + u];
                }

                //
                //  Increment U.  
                //
                //  If U = E, then U = 0 and in Niederreiter's
                //  paper Q = Q + 1.  Here, however, Q is not used explicitly.
                //
                u += 1;
                if (u == e)
                {
                    u = 0;
                }

            }

            //
            //  The array CI now holds the values of C(I,J,R) for this value
            //  of I.  We pack them into array CJ so that CJ(I,R) holds all
            //  the values of C(I,J,R) for J from 1 to NBITS.
            //
            for (r = 0; r < NBITS; r++)
            {
                int term = 0;
                for (j = 0; j < NBITS; j++)
                {
                    term = 2 * term + ci[j,r];
                }

                cj[i,r] = term;
            }

        }
    }

    public class NiederReiter2CalcData
    {
        public const int arbit = 1;
        public const int nonzer = 1;
    }
        
    public static void calcv2(ref NiederReiter2CalcData data, int maxv, int px_deg, int[] px, int[,] add,
            int[,] mul, int[,] sub, ref int b_deg, int[] b,
            int[] v )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CALCV2 calculates the value of the constants V(J,R).
        //
        //  Discussion:
        //
        //    This program calculates the values of the constants V(J,R) as
        //    described in the reference (BFN) section 3.3.  It is called from CALCC2.  
        //
        //    Polynomials stored as arrays have the coefficient of degree N 
        //    in POLY(N).  
        //
        //    A polynomial which is identically 0 is given degree -1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 March 2003
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Paul Bratley, Bennett Fox, Harald Niederreiter.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Paul Bratley, Bennett Fox, Harald Niederreiter,
        //    Algorithm 738: 
        //    Programs to Generate Niederreiter's Low-Discrepancy Sequences,
        //    ACM Transactions on Mathematical Software,
        //    Volume 20, Number 4, pages 494-495, 1994.
        //
        //  Parameters:
        //
        //    Input, int MAXV, the dimension of the array V.
        //
        //    Input, int PX_DEG, the degree of PX.
        //
        //    Input, int PX[MAXDEG+1], the appropriate irreducible polynomial 
        //    for the dimension currently being considered.  
        //
        //    Input, int ADD[2][2], MUL[2][2], SUB[2][2], the addition, multiplication, 
        //    and subtraction tables, mod 2.
        //
        //    Input/output, int *B_DEG, the degree of the polynomial B.
        //
        //    Input/output, int B[MAXDEG+1].  On input, B is the polynomial 
        //    defined in section 2.3 of BFN.  The degree of B implicitly defines 
        //    the parameter J of section 3.3, by degree(B) = E*(J-1).  On output,
        //    B has been multiplied by PX, so its degree is now E * J.
        //
        //    Output, int V[MAXV+1], the computed V array.
        //
        //  Local Parameters:
        //
        //    Local, int ARBIT, indicates where the user can place
        //    an arbitrary element of the field of order 2.  This means 
        //    0 <= ARBIT < 2.  
        //
        //    Local, int BIGM, is the M used in section 3.3.
        //    It differs from the [little] m used in section 2.3,
        //    denoted here by M.
        //
        //    Local, int NONZER, shows where the user must put an arbitrary 
        //    non-zero element of the field.  For the code, this means 
        //    0 < NONZER < 2.
        //
    {
        int[] h = new int[MAXDEG + 1];
        int i;
        int r;
        int term;
        //
        //  The polynomial H is PX**(J-1), which is the value of B on arrival.
        //
        //  In section 3.3, the values of Hi are defined with a minus sign:
        //  don't forget this if you use them later!
        //
        int h_deg = b_deg;

        for (i = 0; i <= h_deg; i++)
        {
            h[i] = b[i];
        }

        //
        //  Multiply B by PX so B becomes PX**J.
        //  In section 2.3, the values of Bi are defined with a minus sign:
        //  don't forget this if you use them later!
        //
        int pb_deg = b_deg;

        plymul2(add, mul, px_deg, px, pb_deg, b, ref pb_deg, b);

        b_deg = pb_deg;
        int m = b_deg;
        //
        //  Now choose a value of Kj as defined in section 3.3.
        //  We must have 0 <= Kj < E*J = M.
        //  The limit condition on Kj does not seem very relevant
        //  in this program.
        //
        //
        //  Choose values of V in accordance with the conditions in section 3.3.
        //
        for (r = 0; r < h_deg; r++)
        {
            v[r] = 0;
        }

        v[h_deg] = 1;

        if (h_deg < h_deg)
        {
            term = sub[0,h[h_deg]];

            for (r = h_deg + 1; r <= h_deg - 1; r++)
            {
                v[r] = NiederReiter2CalcData.arbit;
                //
                //  Check the condition of section 3.3,
                //  remembering that the H's have the opposite sign.
                //
                term = sub[term,mul[h[r],v[r]]];

            }

            //
            //  Now V(BIGM) is anything but TERM.
            //
            v[h_deg] = add[NiederReiter2CalcData.nonzer,term];

            for (r = h_deg + 1; r <= m - 1; r++)
            {
                v[r] = NiederReiter2CalcData.arbit;
            }
        }
        else
        {
            for (r = h_deg + 1; r <= m - 1; r++)
            {
                v[r] = NiederReiter2CalcData.arbit;
            }

        }

        //
        //  Calculate the remaining V's using the recursion of section 2.3,
        //  remembering that the B's have the opposite sign.
        //
        for (r = 0; r <= maxv - m; r++)
        {
            term = 0;
            for (i = 0; i <= m - 1; i++)
            {
                term = sub[term,mul[b[i],v[r + i]]];
            }

            v[r + m] = term;
        }
    }

    public class Niederreiter2Data
    {
        public int seed_save;
        public int dim_save;
            
        public int[,] cj = new int[DIM_MAX, NBITS];
        public readonly int[] nextq = new int[DIM_MAX];
        public readonly double RECIP = 1.0 / (1 << NBITS);

        public NiederReiter2CalcData calcdata = new();

    }
    public static void niederreiter2(ref Niederreiter2Data data, int dim_num, ref int seed, ref double[] quasi, int index = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NIEDERREITER2 returns an element of the Niederreiter sequence base 2.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 March 2003
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Paul Bratley, Bennett Fox, Harald Niederreiter.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Harald Niederreiter,
        //    Low-discrepancy and low-dispersion sequences,
        //    Journal of Number Theory,
        //    Volume 30, 1988, pages 51-70.
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the dimension of the sequence to be generated.
        //
        //    Input/output, int *SEED, the index of the element entry to
        //    compute.  On output, SEED is typically reset by this routine
        //    to SEED+1.
        //
        //    Output, double QUASI[DIM_NUM], the next quasirandom vector.
        //
        //  Local Parameters:
        //
        //    Local, int CJ(DIM_MAX,0:NBITS-1), the packed values of 
        //    Niederreiter's C(I,J,R).
        //
        //    Local, int DIM_SAVE, the spatial dimension of the sequence
        //    as specified on an initialization call.
        //
        //    Local, int COUNT, the index of the current item in the sequence,
        //    expressed as an array of bits.  COUNT(R) is the same as Niederreiter's
        //    AR(N) (page 54) except that N is implicit.
        //
        //    Local, int NEXTQ[DIM_MAX], the numerators of the next item in the
        //    series.  These are like Niederreiter's XI(N) (page 54) except that
        //    N is implicit, and the NEXTQ are integers.  To obtain
        //    the values of XI(N), multiply by RECIP.
        //
    {
        int i;
        int r;
        //
        //  Initialization.
        //
        if (data.dim_save < 1 || dim_num != data.dim_save || seed <= 0)
        {
            if (dim_num <= 0 || DIM_MAX < dim_num)
            {
                Console.WriteLine("");
                Console.WriteLine("NIEDERREITER2 - Fatal error!");
                Console.WriteLine("  Bad spatial dimension.");
                return;
            }

            data.dim_save = dim_num;

            seed = seed switch
            {
                < 0 => 0,
                _ => seed
            };

            data.seed_save = seed;
            //
            //  Calculate the C array.
            //
            calcc2(ref data.calcdata, data.dim_save, ref data.cj);
        }

        //
        //  Set up NEXTQ appropriately, depending on the Gray code of SEED.
        //
        //  You can do this every time, starting NEXTQ back at 0,
        //  or you can do it once, and then carry the value of NEXTQ
        //  around from the previous computation.
        //
        if (seed != data.seed_save + 1)
        {
            int gray = seed ^ (seed / 2);

            for (i = 0; i < data.dim_save; i++)
            {
                data.nextq[i] = 0;
            }

            r = 0;

            while (gray != 0)
            {
                if (gray % 2 != 0)
                {
                    for (i = 0; i < data.dim_save; i++)
                    {
                        data.nextq[i] ^= data.cj[i,r];
                    }
                }

                gray /= 2;
                r += 1;
            }
        }

        //
        //  Multiply the numerators in NEXTQ by RECIP to get the next
        //  quasi-random vector.
        //
        for (i = 0; i < data.dim_save; i++)
        {
            quasi[index + i] = data.nextq[i] * data.RECIP;
        }

        //
        //  Find the position of the right-hand zero in SEED.  This
        //  is the bit that changes in the Gray-code representation as
        //  we go from SEED to SEED+1.
        //
        r = 0;
        i = seed;

        while (i % 2 != 0)
        {
            r += 1;
            i /= 2;
        }

        //
        //  Check that we have not passed 2^NBITS calls.
        //
        if (NBITS <= r)
        {
            Console.WriteLine("");
            Console.WriteLine("NIEDERREITER2 - Fatal error!");
            Console.WriteLine("  Too many calls!");
            return;
        }

        //
        //  Compute the new numerators in vector NEXTQ.
        //
        for (i = 0; i < data.dim_save; i++)
        {
            data.nextq[i] ^= data.cj[i,r];
        }

        data.seed_save = seed;
        seed += 1;

    }

    public static double[] niederreiter2_generate(ref Niederreiter2Data data, int dim_num, int n, ref int seed)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NIEDERREITER2_GENERATE generates a set of Niederreiter values.
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
        //    Input, int N, the number of points desired.
        //
        //    Input/output, int *SEED, a seed for the random 
        //    number generator.
        //
        //    Output, double R[DIM_NUM*N], the points.
        //
    {
        int j;

        double[] r = new double[dim_num * n];

        for (j = 0; j < n; j++)
        {
            niederreiter2(ref data, dim_num, ref seed, ref r, index: + j * dim_num);
        }

        return r;
    }

    public static void plymul2(int[,] add, int[,] mul, int pa_deg,
            int[] pa, int pb_deg, int[] pb,
            ref int pc_deg, int[] pc )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PLYMUL2 multiplies two polynomials in the field of order 2
        //
        //  Discussion:
        //
        //    Polynomials stored as arrays have the coefficient of degree N in 
        //    POLY(N), and the degree of the polynomial in POLY(-1).  
        //
        //    A polynomial which is identically 0 is given degree -1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 March 2003
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Paul Bratley, Bennett Fox, Harald Niederreiter.
        //    C++ version by John Burkardt.
        //
        //  Parameters:
        //
        //    Input, int ADD[2][2], MUL[2][2], 
        //    the addition and multiplication tables, mod 2.
        //
        //    Input, int PA_DEG, the degree of PA.
        //
        //    Input, int PA[MAXDEG+1], the first polynomial factor.
        //
        //    Input, int PB_DEG, the degree of PB.
        //
        //    Input, int PB[MAXDEG+1], the second polynomial factor.
        //
        //    Output, int *PC_DEG, the degree of the product.
        //
        //    Output, int PC[MAXDEG+1], the product polynomial.
        //
    {
        int i;
        int jlo = 0;
        int[] pt = new int[MAXDEG + 1];

        if (pa_deg == -1 || pb_deg == -1)
        {
            pc_deg = -1;
        }
        else
        {
            pc_deg = pa_deg + pb_deg;
        }

        if (MAXDEG < pc_deg)
        {
            Console.WriteLine("");
            Console.WriteLine("PLYMUL2 - Fatal error!");
            Console.WriteLine("  Degree of the product exceeds MAXDEG.");
            return;
        }

        for (i = 0; i <= pc_deg; i++)
        {
            jlo = jlo switch
            {
                < 0 => 0,
                _ => i - pa_deg
            };

            int jhi = pb_deg;
            if (i < jhi)
            {
                jhi = i;
            }

            int term = 0;

            int j;
            for (j = jlo; j <= jhi; j++)
            {
                term = add[term,mul[pa[i - j],pb[j]]];
            }

            pt[i] = term;
        }

        for (i = 0; i <= pc_deg; i++)
        {
            pc[i] = pt[i];
        }

        for (i = pc_deg + 1; i <= MAXDEG; i++)
        {
            pc[i] = 0;
        }

    }

    public static void setfld2(int[,] add, int[,] mul, int[,] sub )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SETFLD2 sets up arithmetic tables for the finite field of order 2.
        //
        //  Discussion:
        //
        //    SETFLD2 sets up addition, multiplication, and subtraction tables 
        //    for the finite field of order QIN.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    29 March 2003
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Paul Bratley, Bennett Fox, Harald Niederreiter.
        //    C++ version by John Burkardt.
        //
        //  Parameters:
        //
        //    Input, int ADD[2][2], MUL[2][2], SUB[2][2], the addition, multiplication, 
        //    and subtraction tables, mod 2.
        //
    {
        int i;
        int j;
        const int p = 2;
        const int q = 2;
        //
        for (i = 0; i < q; i++)
        {
            for (j = 0; j < q; j++)
            {
                add[i,j] = (i + j) % p;
                mul[i,j] = i * j % p;
            }
        }

        //
        //  Use the addition table to set the subtraction table.
        //
        for (i = 0; i < q; i++)
        {
            for (j = 0; j < q; j++)
            {
                sub[add[i,j],i] = j;
            }
        }
    }
}