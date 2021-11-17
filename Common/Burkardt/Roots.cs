using Burkardt.PolynomialNS;

namespace Burkardt.RootsNS;

public static class Roots
{
    public static void roots_to_dif(int nroots, double[] roots, ref int ntab, ref double[] xtab,
            ref double[] diftab )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ROOTS_TO_DIF sets a divided difference table for a polynomial from its roots.
        //
        //  Discussion:
        //
        //    This turns out to be a simple task, because of two facts:
        //
        //    * The divided difference polynomial of one smaller degree which
        //      passes through the values ( ROOT(I), 0 ) is the zero polynomial,
        //      and hence has a zero divided difference table.
        //
        //    * We want a polynomial of one degree higher, but we don't want it
        //      to pass through an addditional point.  Instead, we specify that
        //      the polynomial is MONIC.  This means that the divided difference
        //      table is almost the same as for the zero polynomial, except that
        //      there is one more pair of entries, an arbitrary X value, and
        //      a Y value of 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 August 2002
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NROOTS, is the number of roots.
        //
        //    Input, double ROOTS[NROOTS], the roots of the polynomial.
        //
        //    Output, int *NTAB, the size of the divided difference table.
        //    This is NROOTS+1
        //
        //    Output, double XTAB[NTAB], the abscissas of the table.
        //
        //    Output, double DIFTAB[NTAB], the divided difference table.
        //
    {
        int i;

        ntab = nroots + 1;
        //
        //  Build the appropriate difference table for the polynomial
        //  through ( ROOTS(I), 0 ) of degree NTAB-2.
        //
        for (i = 0; i < ntab - 1; i++)
        {
            diftab[i] = 0.0;
        }

        for (i = 0; i < ntab - 1; i++)
        {
            xtab[i] = roots[i];
        }

        //
        //  Append the extra data to make a monic polynomial of degree NTAB-1
        //  which is zero at the NTAB-1 roots.
        //
        xtab[ntab - 1] = 0.0;
        diftab[ntab - 1] = 1.0;
    }

    public static void roots_to_r8poly(int nroots, double[] roots, ref int nc, ref double[] c )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ROOTS_TO_R8POLY converts polynomial roots to polynomial coefficients.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    06 September 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int NROOTS, the number of roots specified.
        //
        //    Input, double ROOTS[NROOTS], the roots.
        //
        //    Output, integer *NC, the order of the polynomial, which will be NROOTS+1.
        //
        //    Output, double C[*NC], the coefficients of the polynomial.
        //
    {
        int i;

        nc = nroots + 1;
        //
        //  Initialize C to (0, 0, ..., 0, 1).
        //  Essentially, we are setting up a divided difference table.
        //
        double[] xtab = new double[nroots + 1];
        for (i = 0; i < nc - 1; i++)
        {
            xtab[i] = roots[i];
        }

        xtab[nc - 1] = 0.0;

        for (i = 0; i < nc - 1; i++)
        {
            c[i] = 0.0;
        }

        c[nc - 1] = 1.0;
        //
        //  Convert to standard polynomial form by shifting the abscissas
        //  of the divided difference table to 0.
        //
        Dif.dif_shift_zero(nc, ref xtab, ref c);
    }
}