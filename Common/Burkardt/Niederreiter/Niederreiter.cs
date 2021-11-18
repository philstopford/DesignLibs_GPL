using System;
using System.Collections.Generic;
using System.Globalization;
using System.IO;
using Burkardt.Types;

namespace Burkardt.NiederreiterNS;

public static class Niederreiter
{
    //
    //  GLOBAL DATA "/FIELD/"
    //
    //    The following GLOBAL data, used by many functions,
    //    gives the order Q of a field, its characteristic P, and its
    //    addition, multiplication, and subtraction tables.
    //
    //    Global, int DEG_MAX, the maximum degree of the polynomials
    //    to be considered.
    //
    //    Global, int P, the characteristic of the field.
    //
    //    Global, int Q, the order of the field.
    //
    //    Global, int Q_MAX, the order of the largest field to
    //    be handled.
    //
    //    Global, int ADD[Q_MAX][Q_MAX], the field addition table. 
    //
    //    Global, int MUL[Q_MAX][Q_MAX], the field multiplication table. 
    //
    //    Global, int SUB[Q_MAX][Q_MAX], the field subtraction table.
    //
    private static int DEG_MAX = 50;
    private static int P;
    private static int Q;
    private static int Q_MAX = 50;

    private static int[,] add = new int[Q_MAX, Q_MAX];
    private static int[,] mul = new int[Q_MAX, Q_MAX];

    private static int[,] sub = new int[Q_MAX, Q_MAX];

    //
    //  GLOBAL DATA "/COMM/"
    //
    //    Global, int DIM_MAX, the maximum dimension that will 
    //    be used.
    //
    //    Global, int FIG_MAX, the maximum number of base-Q digits 
    //    we can handle.  BASE^FIG_MAX - 1 must be representable in C++.
    //    For base 2, this implies that FIG_MAX could be as high as 31.
    //    In the original version of the program, FIG_MAX was set to 20.
    //
    //    Global, int C[DIM_MAX,FIG_MAX,0:FIG_MAX-1], the values of 
    //    Niederreiter's C(I,J,R).
    //
    //    Global, int COUNT[0:FIG_MAX-1], the index of the current item 
    //    in the sequence, expressed as an array of base-Q digits.  COUNT(R)
    //    is the same as Niederreiter's AR(N) (page 54) except that N is implicit.
    //
    //    Global, int D[DIM_MAX][FIG_MAX].
    //
    //    Global, int DIMEN, the dimension of the sequence to be generated.
    //
    //    Global, int NEXTQ[DIM_MAX], the numerators of the next item in 
    //    the series.  These are like Niederreiter's XI(N) (page 54) except that
    //    N is implicit, and the NEXTQ are integers.  To obtain the values of 
    //    XI(N), multiply by RECIP.
    //
    //    Global, int NFIGS, the number of base Q digits we are using.
    //
    //    Global, int QPOW[FIG_MAX], to speed things up a bit. 
    //    QPOW(I) = Q ** (NFIGS-I).
    //
    //    Global, double RECIP = 1.0 / Q^NFIGS.
    //
    private static int DIM_MAX = 50;

    //const int FIG_MAX = 20;
    private static int FIG_MAX = 31;

    private static int[,,] C = new int [DIM_MAX, FIG_MAX, FIG_MAX];
    private static int[] COUNT = new int[FIG_MAX];
    private static int[,] D = new int[DIM_MAX, FIG_MAX];
    private static int DIMEN;
    private static int[] NEXTQ = new int[DIM_MAX];
    private static int NFIGS;
    private static int[] QPOW = new int[FIG_MAX];
    private static double RECIP;

    public static void calcc()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CALCC calculates the value of the constants C(I,J,R).
        //
        //  Discussion:
        //
        //    This routine calculates the values of the constants C(I,J,R).
        //    As far as possible, we use Niederreiter's notation.
        //    We calculate the values of C for each I in turn.
        //    When all the values of C have been calculated, we return
        //    this array to the calling program.
        //
        //    Irreducible polynomials are read from file "gfplys.txt"
        //    This file should have been created earlier by running the
        //    GFPLYS program.
        //
        //    Polynomials stored as arrays have the coefficient of degree n 
        //    in POLY(N), and the degree of the polynomial in POLY(-1).  
        //    The parameter DEG is just to remind us of this last fact.  
        //    A polynomial which is identically 0 is given degree -1.
        //
        //    Thanks to Michael Baudin for pointing out that MAXE should
        //    be increased from 5 to 7, since one of the irreducible polynomials
        //    that must be stored in PX has degree 7, 07 June 2010.
        //
        //    In fact, the size of MAXE depends on the highest degree polynomial
        //    computed by GFPLYS, which in turn depends in part on the 
        //    value NPOL in that program.  To allow DIM_MAX = 50, we increased
        //    NPOL to 50 in GFPLYS and here, and hence MAXE to 8, JVB, 07 June 2010.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    07 June 2010
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Paul Bratley, Bennett Fox, 
        //    Harald Niederreiter.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Paul Bratley, Bennett Fox, Harald Niederreiter,
        //    Algorithm 738: 
        //    Programs to Generate NiederreiteSoftware,
        //    Volume 20, Number 4, 1994, pages 494-495.
        //
        //    Harald Niederreiter,
        //    Low-discrepancy and low-dispersion sequences,
        //    Journal of Number Theory,
        //    Volume 30, 1988, pages 51-70.
        //
        //  Local Parameters:
        //
        //    Local, int MAXE; we need DIM_MAX irreducible polynomials over GF(Q).
        //    MAXE is the highest degree among these.
        //
        //    Local, int V_MAX, the maximum index used in V.
        //
        //    Local, int NPOLS, the number of precalculated irreducible polynomials.
        //
    {
        //const int maxe = 5;
        int maxe = 8;

        int v_max = FIG_MAX + maxe;

        int[] b = new int[DEG_MAX + 2];
        int e;
        string[] input;
        string input_filename = "gfplys.txt";
        int i;
        int j;
        int k;
        //const int npols = 25;
        const int npols = 50;
        int[] px = new int[maxe + 2];
        int r;
        int u;
        int[] v = new int[v_max + 1];
        //
        //  Read the irreducible polynomials.
        //
        input = File.ReadAllLines(input_filename);

        int index = 0;

        while (true)
        {
            i = Convert.ToInt32(input[index].Trim());
            index++;

            if (i == Q)
            {
                break;
            }

            for (j = 1; j <= npols; j++)
            {
                string[] tokens = Helpers.splitStringByWhitespace(input[index]);
                e = Convert.ToInt32(tokens[0]);
                index++;
                for (k = 0; k <= e; k++)
                {
                    px[k + 1] = Convert.ToInt32(tokens[k+1]);
                }
            }
        }

        for (i = 0; i < DIMEN; i++)
        {
            //
            //  For each dimension, we need to calculate powers of an
            //  appropriate irreducible polynomial.  See Niederreiter
            //  page 65, just below equation (19).
            //
            //  Read the appropriate irreducible polynomial into PX,
            //  and its degree into E.  Set polynomial B = PX^0 = 1.
            //  M is the degree of B.  Subsequently B will hold higher
            //  powers of PX.
            //
            //  The polynomial PX is stored in 'gfplys.txt' in the format
            //
            //    n  a0  a1  a2  ... an
            //
            //  where n is the degree of the polynomial and the ai are
            //  its coefficients.
            //
            string[] tokens = Helpers.splitStringByWhitespace(input[index]);
            e = Convert.ToInt32(tokens[0]);
            index++;
            for (k = 0; k <= e; k++)
            {
                px[k + 1] = Convert.ToInt32(tokens[k+1]);
            }

            px[0] = e;
            b[0] = 0;
            b[1] = 1;
            //
            //  Niederreiter (page 56, after equation (7), defines two variables 
            //  Q and U.  We do not need Q explicitly, but we do need U.
            //
            u = 0;

            for (j = 0; j < NFIGS; j++)
            {
                switch (u)
                {
                    //
                    //  If U = 0, we need to set B to the next power of PX
                    //  and recalculate V.  This is done by subroutine CALCV.
                    //
                    case 0:
                        calcv(px, ref b, v, v_max);
                        break;
                }

                //
                //  Now C is obtained from V.  Neiderreiter obtains A from V 
                //  (page 65, near the bottom), and then gets C from A (page 56,
                //  equation (7)).  However this can be done in one step.
                //
                for (r = 0; r < NFIGS; r++)
                {
                    C[i, j, r] = v[r + u];
                }

                //
                //  Increment U.  If U = E, then U = 0 and in Niederreiter's
                //  paper Q = Q + 1.  Here, however, Q is not used explicitly.
                //
                u += 1;
                if (u == e)
                {
                    u = 0;
                }
            }
        }
    }

    public static void calcv(int[] px, ref int[] b, int[] v, int v_max)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    CALCV calculates the constants V(J,R).
        //
        //  Discussion:
        //
        //    This program calculates the values of the constants V(J,R) as
        //    described in Bratley, Fox and Niederreiter, section 3.3.  It 
        //    is called from either CALCC or CALCC2.  The values transmitted 
        //    through common /FIELD/ determine which field we are working in.
        //
        //    Polynomials stored as arrays have the coefficient of degree n 
        //    in POLY(N), and the degree of the polynomial in POLY(-1).  The 
        //    parameter DEG is just to remind us of this last fact.  A polynomial 
        //    which is identically 0 is given degree -1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 September 2007
        //
        //  Author:
        //
        //    Orginal FORTRAN77 version by Paul Bratley, Bennett Fox, 
        //    Harald Niederreiter.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Paul Bratley, Bennett Fox, Harald Niederreiter,
        //    Algorithm 738: 
        //    Programs to Generate Niederreiter's Low-Discrepancy Sequences,
        //    ACM Transactions on Mathematical Software,
        //    Volume 20, Number 4, 1994, pages 494-495.
        //
        //  Parameters:
        //
        //    Input, int PX[MAXDEG+2], the appropriate irreducible polynomial 
        //    for the dimension currently being considered.  The degree of PX will 
        //    be called E.
        //
        //    Input/output, int B[DEG_MAX+2].  On input, B is the polynomial 
        //    defined in section 2.3 of BFN.  The degree of B implicitly defines 
        //    the parameter J of section 3.3, by degree(B) = E*(J-1).  On output, 
        //    B has been multiplied by PX, so its degree is now E*J.
        //
        //    Input, int V[V_MAX+1], contains the values required.
        //
        //    Input, int V_MAX, the dimension of the array V.
        //
        //  Local Parameters:
        //
        //    Local, int ARBIT, indicates where the user can place
        //    an arbitrary element of the field of order Q.  For the code,
        //    this means 0 <= ARBIT < Q.  Within these limits, the user can 
        //    do what he likes.  ARBIT could be declared as a function 
        //    that returned a different arbitrary value each time it is referenced.
        //
        //    Local, int BIGM, is the M used in section 3.3.  It differs from 
        //    the [little] m used in section 2.3, denoted here by M.
        //
        //    Local, int NONZER shows where the user must put an arbitrary 
        //    non-zero element of the same field.  For the code, this means 
        //    0 < NONZER < Q.  Within these limits, the user can do what he likes.  
        //
    {
        int arbit = 1;
        int[] b2;
        int bigm;
        int[] h = new int[DEG_MAX + 2];
        int i;
        int kj;
        int m;
        int nonzer = 1;
        int r;
        int term;

        //e = px[0];
        //
        //  The polynomial H is PX^(J-1), which is the value of B on arrival.
        //
        //  In section 3.3, the values of Hi are defined with a minus sign:
        //  don't forget this if you use them later//
        //
        for (i = 0; i < DEG_MAX + 2; i++)
        {
            h[i] = b[i];
        }

        bigm = h[0];
        //
        //  Now multiply B by PX so B becomes PX^J.
        //
        //  In section 2.3, the values of Bi are defined with a minus sign:
        //  don't forget this if you use them later!
        //
        b2 = plymul(px, b);

        for (i = 0; i < DEG_MAX + 2; i++)
        {
            b[i] = b2[i];
        }

        m = b[0];
        //
        //  Now choose a value of Kj as defined in section 3.3.
        //  We must have 0 <= Kj < E*J = M.
        //  The limit condition on Kj does not seem very relevant
        //  in this program.
        //
        kj = bigm;
        //
        //  Now choose values of V in accordance with the conditions in
        //  section 3.3
        //
        for (i = 0; i < kj; i++)
        {
            v[i] = 0;
        }

        v[kj] = 1;

        if (kj < bigm)
        {
            term = sub[0, h[kj + 1]];

            for (r = kj + 1; r <= bigm - 1; r++)
            {
                v[r] = arbit;
                //
                //  Check the condition of section 3.3,
                //  remembering that the H's have the opposite sign.
                //
                term = sub[term, mul[h[r + 1], v[r]]];
            }

            //
            //  Now V(BIGM) is anything but TERM.
            //
            v[bigm] = add[nonzer, term];
            for (i = bigm + 1; i <= m - 1; i++)
            {
                v[i] = arbit;
            }
        }
        else
        {
            for (i = kj + 1; i <= m - 1; i++)
            {
                v[i] = arbit;
            }
        }

        //
        //  Calculate the remaining V's using the recursion of section 2.3,
        //  remembering that the B's have the opposite sign.
        //
        for (r = 0; r <= v_max - m; r++)
        {
            term = 0;
            for (i = 0; i <= m - 1; i++)
            {
                term = sub[term, mul[b[i + 1], v[r + i]]];
            }

            v[r + m] = term;
        }

    }

    public static void golo(double[] quasi, int index = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GOLO generates a new quasi-random vector on each call.
        //
        //  Discussion:
        //
        //    Before the first call to this routine, a call must be made
        //    to subroutine INLO to carry out some initializations.
        //
        //    Polynomials stored as arrays have the coefficient of degree n 
        //    in POLY(N), and the degree of the polynomial in POLY(-1).  
        //    The parameter DEG is just to remind us of this last fact.  
        //    A polynomial which is identically 0 is given degree -1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 September 2007
        //
        //  Author:
        //
        //    Orginal FORTRAN77 version by Paul Bratley, Bennett Fox, 
        //    Harald Niederreiter.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Paul Bratley, Bennett Fox, Harald Niederreiter,
        //    Algorithm 738: 
        //    Programs to Generate Niederreiter's Low-Discrepancy Sequences,
        //    ACM Transactions on Mathematical Software,
        //    Volume 20, Number 4, 1994, pages 494-495.
        //
        //    Harald Niederreiter,
        //    Low-discrepancy and low-dispersion sequences,
        //    Journal of Number Theory,
        //    Volume 30, 1988, pages 51-70.
        //
        //  Parameters:
        //
        //    Output, double QUASI[], the next vector in the sequence.
        //
    {
        int diff;
        int i;
        int j;
        int nq;
        int oldcnt;
        int r;
        //
        //  Multiply the numerators in NEXTQ by RECIP to get the next
        //  quasi-random vector.
        //
        for (i = 0; i < DIMEN; i++)
        {
            quasi[index + i] = NEXTQ[i] * RECIP;
        }

        //
        //  Update COUNT, treated as a base-Q integer.  Instead of
        //  recalculating the values of D from scratch, we update
        //  them for each digit of COUNT which changes.  In terms of
        //  Niederreiter page 54, NEXTQ(I) corresponds to XI(N), with
        //  N being implicit, and D(I,J) corresponds to XI(N,J), again
        //  with N implicit.  Finally COUNT(R) corresponds to AR(N).
        //
        r = 0;

        for (;;)
        {
            if (NFIGS <= r)
            {
                Console.WriteLine("");
                Console.WriteLine("GOLO - Fatal error!");
                Console.WriteLine("  Too many calls!");
                return;
            }

            oldcnt = COUNT[r];

            if (COUNT[r] < Q - 1)
            {
                COUNT[r] += 1;
            }
            else
            {
                COUNT[r] = 0;
            }

            diff = sub[COUNT[r], oldcnt];
            //
            //  Digit R has just changed.  DIFF says how much it changed
            //  by.  We use this to update the values of array D.
            //
            for (i = 0; i < DIMEN; i++)
            {
                for (j = 0; j < NFIGS; j++)
                {
                    D[i, j] = add[D[i, j], mul[C[i, j, r], diff]];
                }
            }

            //
            //  If COUNT(R) is now zero, we must propagate the carry.
            //
            if (COUNT[r] != 0)
            {
                break;
            }

            r += 1;
        }

        //
        //  Now use the updated values of D to calculate NEXTQ.
        //  Array QPOW helps to speed things up a little:
        //  QPOW(J) is Q^(NFIGS-J).
        //
        for (i = 0; i < DIMEN; i++)
        {
            nq = 0;
            for (j = 0; j < NFIGS; j++)
            {
                nq += D[i, j] * QPOW[j];
            }

            NEXTQ[i] = nq;
        }

    }

    public static void inlo(int dim, int base_, int skip)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INLO calculates the values of C(I,J,R).
        //
        //  Discussion:
        //
        //    This subroutine calculates the values of Niederreiter's
        //    C(I,J,R) and performs other initialization necessary
        //    before calling GOLO.
        //
        //    Polynomials stored as arrays have the coefficient of degree n 
        //    in POLY(N), and the degree of the polynomial in POLY(-1).  
        //    The parameter DEG is just to remind us of this last fact.  
        //    A polynomial which is identically 0 is given degree -1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 September 2007
        //
        //  Author:
        //
        //    Orginal FORTRAN77 version by Paul Bratley, Bennett Fox, 
        //    Harald Niederreiter.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Paul Bratley, Bennett Fox, Harald Niederreiter,
        //    Algorithm 738: 
        //    Programs to Generate Niederreiter's Low-Discrepancy Sequences,
        //    ACM Transactions on Mathematical Software,
        //    Volume 20, Number 4, 1994, pages 494-495.
        //
        //    Harald Niederreiter,
        //    Low-discrepancy and low-dispersion sequences,
        //    Journal of Number Theory,
        //    Volume 30, 1988, pages 51-70.
        //
        //  Parameters:
        //
        //    Input, int DIM, the dimension of the sequence to be generated.
        //    The value of DIM is copied into DIMEN in the common block.
        //
        //    Input, int BASE, the prime or prime-power base to be used.
        //
        //    Input, int SKIP, the number of values to throw away at the 
        //    beginning of the sequence.
        //
        //  Local Parameters:
        //
        //    Local, int NBITS, the number of bits in a fixed-point integer, not
        //    counting the sign.
        //
    {
        int i;
        int j;
        int nbits = 31;
        int nq;
        int r;
        double temp;

        DIMEN = dim;

        if (DIMEN <= 0 || DIM_MAX < DIMEN)
        {
            Console.WriteLine("");
            Console.WriteLine("INLO - Fatal error!");
            Console.WriteLine("  Bad spatial dimension.");
            return;
        }

        if (typeMethods.i4_characteristic(base_) == 0)
        {
            Console.WriteLine("");
            Console.WriteLine("INLO - Fatal error!");
            Console.WriteLine("  Base not prime power or out of range.");
            return;
        }

        setfld(base_);
        //
        //  Calculate how many figures to use in base Q = BASE
        //
        temp = Math.Log(Math.Pow(2.0, nbits) - 1.0) / Math.Log(Q);

        NFIGS = Math.Min(FIG_MAX, (int) temp);
        //
        //  Calculate the C array.
        //
        calcc();
        //
        //  Set RECIP.
        //
        RECIP = 1.0 / Math.Pow(Q, NFIGS);
        //
        //  Set QPOW(I) = Q^(NFIGS-I).
        //
        QPOW[NFIGS - 1] = 1;
        for (i = NFIGS - 1; 1 <= i; i--)
        {
            QPOW[i - 1] = Q * QPOW[i];
        }

        //
        //  Initialize COUNT.
        //
        i = skip;

        for (r = 0; r < NFIGS; r++)
        {
            COUNT[r] = i % Q;
            i /= Q;
        }

        if (i != 0)
        {
            Console.WriteLine("");
            Console.WriteLine("INLO - Fatal error!");
            Console.WriteLine("  SKIP is too long!");
            return;
        }

        //
        //  Initialize D.
        //
        for (i = 0; i < DIMEN; i++)
        {
            for (j = 0; j < NFIGS; j++)
            {
                D[i, j] = 0;
            }
        }

        for (r = 0; r < NFIGS; r++)
        {
            if (COUNT[r] != 0)
            {
                for (i = 0; i < DIMEN; i++)
                {
                    for (j = 0; j < NFIGS; j++)
                    {
                        D[i, j] = add[D[i, j], mul[C[i, j, r], COUNT[r]]];
                    }
                }
            }
        }

        //
        //  Initialize NEXTQ.
        //
        for (i = 0; i < DIMEN; i++)
        {
            nq = 0;
            for (j = 0; j < NFIGS; j++)
            {
                nq += D[i, j] * QPOW[j];
            }

            NEXTQ[i] = nq;
        }

    }

    public class NiederReiterData
    {
        public int dim_num_save = -1;

    }

    public static void niederreiter(ref NiederReiterData data, int dim_num, int base_, ref int seed, ref double[] r, int index = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NIEDERREITER returns an element of a Niederreiter sequence for base BASE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 September 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int BASE, the base to use for the Niederreiter sequence.
        //    The base should be a prime, or a power of a prime.
        //
        //    Input/output, int *SEED, a seed for the random number generator.
        //
        //    Output, double R[DIM_NUM], the element of the sequence.
        //
    {
        int skip;

        if (data.dim_num_save < 1 || dim_num != data.dim_num_save || seed <= 0)
        {
            skip = 1;

            inlo(dim_num, base_, skip);

            data.dim_num_save = dim_num;
        }

        golo(r, index);

        seed += 1;

    }

    public static void niederreiter_generate(ref NiederReiterData data, int dim_num, int n, int base_, ref int seed,
            ref double[] r)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NIEDERREITER_GENERATE generates a set of Niederreiter values.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    11 September 2007
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
        //    Input, int BASE, the base to use for the Niederreiter sequence.
        //    The base should be a prime, or a power of a prime.
        //
        //    Input/output, int *SEED, a seed for the random number generator.
        //
        //    Output, double R[DIM_NUM*N], the points.
        //
    {
        int j;

        for (j = 0; j < n; j++)
        {
            niederreiter(ref data, dim_num, base_, ref seed, ref r, index: +j * dim_num);
        }
    }

    public static void niederreiter_write(int dim_num, int n, int base_, int skip, double[] r,
            string output_filename)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NIEDERREITER_WRITE writes a set of Niederreiter values to a file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    12 September 2007
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
        //    Input, int BASE, the base.
        //
        //    Input, int SKIP, the number of initial points skipped.
        //
        //    Input, double R[DIM_NUM*N], the points.
        //
        //    Input, char *OUTPUT_FILENAME, the name of the
        //    file to which the output should be written.
        //
    {
        int dim;
        int j;
        List<string> output = new()
        {
            "#  " + output_filename + "",
            "#  created by NIEDERREITER.C.",
            "#",
            "#  Spatial dimension DIM_NUM = " + dim_num + "",
            "#  Number of points N = " + n + "",
            "#  EPSILON (unit roundoff) = " + typeMethods.r8_epsilon() + "",
            "#  Base: " + base_ + "",
            "#  Initial values skipped = " + skip + "",
            "#"
        };

        for (j = 0; j < n; j++)
        {
            string cout = "";
            for (dim = 0; dim < dim_num; dim++)
            {
                cout += "  " + r[dim + j * dim_num].ToString(CultureInfo.InvariantCulture).PadLeft(10);
            }

            output.Add(cout);
        }

        File.WriteAllLines(output_filename, output);

    }

    public static int[] plymul(int[] pa, int[] pb)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PLYMUL multiplies one polynomial by another.
        //
        //  Discussion:
        //
        //    Polynomial coefficients are elements of the field of order Q.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 September 2007
        //
        //  Author:
        //
        //    Orginal FORTRAN77 version by Paul Bratley, Bennet Fox, 
        //    Harald Niederreiter.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Paul Bratley, Bennett Fox, Harald Niederreiter,
        //    Algorithm 738: 
        //    Programs to Generate Niederreiter's Low-Discrepancy Sequences,
        //    ACM Transactions on Mathematical Software,
        //    Volume 20, Number 4, 1994, pages 494-495.
        //
        //  Parameters:
        //
        //    Input, int PA[DEG_MAX+2], the first polynomial.
        //
        //    Input, int PB[DEG_MAX+2], the second polynomial.
        //
        //    Output, int PLYMUL[DEG_MAX+2], the product polynomial.
        //
    {
        int dega;
        int degb;
        int degc;
        int i;
        int j;
        int[] pc;
        int term;

        pc = new int[DEG_MAX + 2];

        dega = pa[0];
        degb = pb[0];

        if (dega == -1 || degb == -1)
        {
            degc = -1;
        }
        else
        {
            degc = dega + degb;
        }

        if (DEG_MAX < degc)
        {
            Console.WriteLine("");
            Console.WriteLine("PLYMUL - Fatal error!");
            Console.WriteLine("  The degree of the product exceeds DEG_MAX.");
            return null;
        }

        for (i = 0; i <= degc; i++)
        {
            term = 0;
            for (j = Math.Max(0, i - dega); j <= Math.Min(degb, i); j++)
            {
                term = add[term, mul[pa[i - j + 1], pb[j + 1]]];
            }

            pc[i + 1] = term;
        }

        pc[0] = degc;

        for (i = degc + 1; i <= DEG_MAX; i++)
        {
            pc[i + 1] = 0;
        }

        return pc;
    }

    public static void setfld(int q_init)

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    SETFLD sets up the arithmetic tables for a finite field.
        //
        //  Discussion:
        //
        //    This subroutine sets up addition, multiplication, and
        //    subtraction tables for the finite field of order QIN.
        //
        //    A polynomial with coefficients A(*) in the field of order Q
        //    can also be stored in an integer I, with
        //
        //      I = AN*Q^N + ... + A0.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 September 2007
        //
        //  Author:
        //
        //    Orginal FORTRAN77 version by Paul Bratley, Bennet Fox, 
        //    Harald Niederreiter.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Paul Bratley, Bennett Fox, Harald Niederreiter,
        //    Algorithm 738: 
        //    Programs to Generate Niederreiter's Low-Discrepancy Sequences,
        //    ACM Transactions on Mathematical Software,
        //    Volume 20, Number 4, 1994, pages 494-495.
        //
        //  Parameters:
        //
        //    Input, int Q_INIT, the order of the field.
        //
    {
        int i;
        int j;

        if (q_init <= 1 || Q_MAX < q_init)
        {
            Console.WriteLine("");
            Console.WriteLine("SETFLD - Fatal error!");
            Console.WriteLine("  Bad value of Q = " + q_init + "");
            return;
        }

        Q = q_init;
        P = typeMethods.i4_characteristic(Q);

        switch (P)
        {
            case 0:
                Console.WriteLine("");
                Console.WriteLine("SETFLD - Fatal error!");
                Console.WriteLine("  There is no field of order Q = " + Q + "");
                return;
        }

        //
        //  Set up to handle a field of prime or prime-power order.
        //  Calculate the addition and multiplication tables.
        //
        for (i = 0; i < Q; i++)
        {
            for (j = 0; j < Q; j++)
            {
                add[i,j] = (i + j) % P;
            }
        }

        for (i = 0; i < Q; i++)
        {
            for (j = 0; j < Q; j++)
            {
                mul[i,j] = i * j % P;
            }
        }

        //
        //  Use the addition table to set the subtraction table.
        //
        for (i = 0; i < Q; i++)
        {
            for (j = 0; j < Q; j++)
            {
                sub[add[i,j],i] = j;
            }
        }
    }
}