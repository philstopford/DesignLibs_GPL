namespace Burkardt.CDFLib;

public static partial class CDF
{
    public static int ipmpar(int i)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    IPMPAR returns integer machine constants.
        //
        //  Discussion:
        //
        //    Input arguments 1 through 3 are queries about integer arithmetic.
        //    We assume integers are represented in the N-digit, base-A form
        //
        //      sign * ( X(N-1)*A**(N-1) + ... + X(1)*A + X(0) )
        //
        //    where 0 <= X(0:N-1) < A.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    14 February 2021
        //
        //    Then:
        //
        //      IPMPAR(1) = A, the base of integer arithmetic;
        //      IPMPAR(2) = N, the number of base A digits;
        //      IPMPAR(3) = A^N - 1, the largest magnitude.
        //
        //    It is assumed that the single and double precision floating
        //    point arithmetics have the same base, say B, and that the
        //    nonzero numbers are represented in the form
        //
        //      sign * (B^E) * (X(1)/B + ... + X(M)/B^M)
        //
        //    where X(1:M) is one of { 0, 1,..., B-1 }, and 1 <= X(1) and
        //    EMIN <= E <= EMAX.
        //
        //    Input argument 4 is a query about the base of real arithmetic:
        //
        //      IPMPAR(4) = B, the base of single and double precision arithmetic.
        //
        //    Input arguments 5 through 7 are queries about single precision
        //    floating point arithmetic:
        //
        //     IPMPAR(5) = M, the number of base B digits for single precision.
        //     IPMPAR(6) = EMIN, the smallest exponent E for single precision.
        //     IPMPAR(7) = EMAX, the largest exponent E for single precision.
        //
        //    Input arguments 8 through 10 are queries about double precision
        //    floating point arithmetic:
        //
        //     IPMPAR(8) = M, the number of base B digits for double precision.
        //     IPMPAR(9) = EMIN, the smallest exponent E for double precision.
        //     IPMPAR(10) = EMAX, the largest exponent E for double precision.
        //
        //  Author:
        //
        //    Barry Brown, James Lovato, Kathy Russell.
        //
        //  Reference:
        //
        //    Phyllis Fox, Andrew Hall, and Norman Schryer,
        //    Algorithm 528,
        //    Framework for a Portable FORTRAN Subroutine Library,
        //    ACM Transactions on Mathematical Software,
        //    Volume 4, 1978, pages 176-188.
        //
        //  Parameters:
        //
        //    Input, int *I, the index of the desired constant.
        //
        //    Output, int IPMPAR, the value of the desired constant.
        //
    {
        int[] imach = new int[11];
        //     MACHINE CONSTANTS FOR AMDAHL MACHINES.
        //
        //   imach[1] = 2;
        //   imach[2] = 31;
        //   imach[3] = 2147483647;
        //   imach[4] = 16;
        //   imach[5] = 6;
        //   imach[6] = -64;
        //   imach[7] = 63;
        //   imach[8] = 14;
        //   imach[9] = -64;
        //   imach[10] = 63;
        //
        //     MACHINE CONSTANTS FOR THE AT&T 3B SERIES, AT&T
        //       PC 7300, AND AT&T 6300.
        //
        //   imach[1] = 2;
        //   imach[2] = 31;
        //   imach[3] = 2147483647;
        //   imach[4] = 2;
        //   imach[5] = 24;
        //   imach[6] = -125;
        //   imach[7] = 128;
        //   imach[8] = 53;
        //   imach[9] = -1021;
        //   imach[10] = 1024;
        //
        //     MACHINE CONSTANTS FOR THE BURROUGHS 1700 SYSTEM.
        //
        //   imach[1] = 2;
        //   imach[2] = 33;
        //   imach[3] = 8589934591;
        //   imach[4] = 2;
        //   imach[5] = 24;
        //   imach[6] = -256;
        //   imach[7] = 255;
        //   imach[8] = 60;
        //   imach[9] = -256;
        //   imach[10] = 255;
        //
        //     MACHINE CONSTANTS FOR THE BURROUGHS 5700 SYSTEM.
        //
        //   imach[1] = 2;
        //   imach[2] = 39;
        //   imach[3] = 549755813887;
        //   imach[4] = 8;
        //   imach[5] = 13;
        //   imach[6] = -50;
        //   imach[7] = 76;
        //   imach[8] = 26;
        //   imach[9] = -50;
        //   imach[10] = 76;
        //
        //     MACHINE CONSTANTS FOR THE BURROUGHS 6700/7700 SYSTEMS.
        //
        //   imach[1] = 2;
        //   imach[2] = 39;
        //   imach[3] = 549755813887;
        //   imach[4] = 8;
        //   imach[5] = 13;
        //   imach[6] = -50;
        //   imach[7] = 76;
        //   imach[8] = 26;
        //   imach[9] = -32754;
        //   imach[10] = 32780;
        //
        //     MACHINE CONSTANTS FOR THE CDC 6000/7000 SERIES
        //       60 BIT ARITHMETIC, AND THE CDC CYBER 995 64 BIT
        //       ARITHMETIC (NOS OPERATING SYSTEM).
        //
        //   imach[1] = 2;
        //   imach[2] = 48;
        //   imach[3] = 281474976710655;
        //   imach[4] = 2;
        //   imach[5] = 48;
        //   imach[6] = -974;
        //   imach[7] = 1070;
        //   imach[8] = 95;
        //   imach[9] = -926;
        //   imach[10] = 1070;
        //
        //     MACHINE CONSTANTS FOR THE CDC CYBER 995 64 BIT
        //       ARITHMETIC (NOS/VE OPERATING SYSTEM).
        //
        //   imach[1] = 2;
        //   imach[2] = 63;
        //   imach[3] = 9223372036854775807;
        //   imach[4] = 2;
        //   imach[5] = 48;
        //   imach[6] = -4096;
        //   imach[7] = 4095;
        //   imach[8] = 96;
        //   imach[9] = -4096;
        //   imach[10] = 4095;
        //
        //     MACHINE CONSTANTS FOR THE CRAY 1, XMP, 2, AND 3.
        //
        //   imach[1] = 2;
        //   imach[2] = 63;
        //   imach[3] = 9223372036854775807;
        //   imach[4] = 2;
        //   imach[5] = 47;
        //   imach[6] = -8189;
        //   imach[7] = 8190;
        //   imach[8] = 94;
        //   imach[9] = -8099;
        //   imach[10] = 8190;
        //
        //     MACHINE CONSTANTS FOR THE DATA GENERAL ECLIPSE S/200.
        //
        //   imach[1] = 2;
        //   imach[2] = 15;
        //   imach[3] = 32767;
        //   imach[4] = 16;
        //   imach[5] = 6;
        //   imach[6] = -64;
        //   imach[7] = 63;
        //   imach[8] = 14;
        //   imach[9] = -64;
        //   imach[10] = 63;
        //
        //     MACHINE CONSTANTS FOR THE IBM 360/370 SERIES,
        //       THE ICL 2900, THE ITEL AS/6, THE XEROX SIGMA
        //       5/7/9 AND THE SEL SYSTEMS 85/86.
        //
        //   imach[1] = 2;
        //   imach[2] = 31;
        //   imach[3] = 2147483647;
        //   imach[4] = 16;
        //   imach[5] = 6;
        //   imach[6] = -64;
        //   imach[7] = 63;
        //   imach[8] = 14;
        //   imach[9] = -64;
        //   imach[10] = 63;
        //
        //     MACHINE CONSTANTS FOR THE IBM PC.
        //
        //   imach[1] = 2;
        //   imach[2] = 31;
        //   imach[3] = 2147483647;
        //   imach[4] = 2;
        //   imach[5] = 24;
        //   imach[6] = -125;
        //   imach[7] = 128;
        //   imach[8] = 53;
        //   imach[9] = -1021;
        //   imach[10] = 1024;
        //
        //     MACHINE CONSTANTS FOR THE MACINTOSH II - ABSOFT
        //       MACFORTRAN II.
        //
        //   imach[1] = 2;
        //   imach[2] = 31;
        //   imach[3] = 2147483647;
        //   imach[4] = 2;
        //   imach[5] = 24;
        //   imach[6] = -125;
        //   imach[7] = 128;
        //   imach[8] = 53;
        //   imach[9] = -1021;
        //   imach[10] = 1024;
        //
        //     MACHINE CONSTANTS FOR THE MICROVAX - VMS FORTRAN.
        //
        //   imach[1] = 2;
        //   imach[2] = 31;
        //   imach[3] = 2147483647;
        //   imach[4] = 2;
        //   imach[5] = 24;
        //   imach[6] = -127;
        //   imach[7] = 127;
        //   imach[8] = 56;
        //   imach[9] = -127;
        //   imach[10] = 127;
        //
        //     MACHINE CONSTANTS FOR THE PDP-10 (KA PROCESSOR).
        //
        //   imach[1] = 2;
        //   imach[2] = 35;
        //   imach[3] = 34359738367;
        //   imach[4] = 2;
        //   imach[5] = 27;
        //   imach[6] = -128;
        //   imach[7] = 127;
        //   imach[8] = 54;
        //   imach[9] = -101;
        //   imach[10] = 127;
        //
        //     MACHINE CONSTANTS FOR THE PDP-10 (KI PROCESSOR).
        //
        //   imach[1] = 2;
        //   imach[2] = 35;
        //   imach[3] = 34359738367;
        //   imach[4] = 2;
        //   imach[5] = 27;
        //   imach[6] = -128;
        //   imach[7] = 127;
        //   imach[8] = 62;
        //   imach[9] = -128;
        //   imach[10] = 127;
        //
        //     MACHINE CONSTANTS FOR THE PDP-11 FORTRAN SUPPORTING
        //       32-BIT INTEGER ARITHMETIC.
        //
        //   imach[1] = 2;
        //   imach[2] = 31;
        //   imach[3] = 2147483647;
        //   imach[4] = 2;
        //   imach[5] = 24;
        //   imach[6] = -127;
        //   imach[7] = 127;
        //   imach[8] = 56;
        //   imach[9] = -127;
        //   imach[10] = 127;
        //
        //     MACHINE CONSTANTS FOR THE SEQUENT BALANCE 8000.
        //
        //   imach[1] = 2;
        //   imach[2] = 31;
        //   imach[3] = 2147483647;
        //   imach[4] = 2;
        //   imach[5] = 24;
        //   imach[6] = -125;
        //   imach[7] = 128;
        //   imach[8] = 53;
        //   imach[9] = -1021;
        //   imach[10] = 1024;
        //
        //     MACHINE CONSTANTS FOR THE SILICON GRAPHICS IRIS-4D
        //       SERIES (MIPS R3000 PROCESSOR).
        //
        //   imach[1] = 2;
        //   imach[2] = 31;
        //   imach[3] = 2147483647;
        //   imach[4] = 2;
        //   imach[5] = 24;
        //   imach[6] = -125;
        //   imach[7] = 128;
        //   imach[8] = 53;
        //   imach[9] = -1021;
        //   imach[10] = 1024;
        //
        //     MACHINE CONSTANTS FOR IEEE ARITHMETIC MACHINES, SUCH AS THE AT&T
        //       3B SERIES, MOTOROLA 68000 BASED MACHINES (E.G. SUN 3 AND AT&T
        //       PC 7300), AND 8087 BASED MICROS (E.G. IBM PC AND AT&T 6300).

        imach[1] = 2;
        imach[2] = 31;
        imach[3] = 2147483647;
        imach[4] = 2;
        imach[5] = 24;
        imach[6] = -125;
        imach[7] = 128;
        imach[8] = 53;
        imach[9] = -1021;
        imach[10] = 1024;

        int ipmpar = imach[i];
        return ipmpar;
    }
}