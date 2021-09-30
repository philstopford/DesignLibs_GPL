using System;
using Burkardt.MatrixNS;
using Burkardt.Types;

namespace Burkardt.PolynomialNS
{
    public static class Dif
    {
        public static void dif_deriv(int nd, double[] xd, double[] yd, ref int ndp, double[] xdp,
                ref double[] ydp)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DIF_DERIV computes the derivative of a polynomial in divided difference form.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 June 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Carl deBoor,
            //    A Practical Guide to Splines,
            //    Springer, 2001,
            //    ISBN: 0387953663,
            //    LC: QA1.A647.v27.
            //
            //  Parameters:
            //
            //    Input, int ND, the size of the input table.
            //
            //    Input, double XD[ND], the abscissas for the divided
            //    difference table.
            //
            //    Input, double YD[ND], the divided difference table.
            //
            //    Output, int *NDP, the size of the output table, which is ND-1.
            //
            //    Input, double XDP[NP], the abscissas for the divided
            //    difference table for the derivative.
            //
            //    Output, double YDP[NDP], the divided difference
            //    table for the derivative.
            //
        {
            int i;
            double[] xd_temp;
            double[] yd_temp;
            //
            //  Using a temporary copy of the difference table, shift the
            //  abscissas to zero.
            //
            xd_temp = new double[nd];
            yd_temp = new double[nd];

            for (i = 0; i < nd; i++)
            {
                xd_temp[i] = xd[i];
            }

            for (i = 0; i < nd; i++)
            {
                yd_temp[i] = yd[i];
            }

            dif_shift_zero(nd, ref xd_temp, ref yd_temp);
            //
            //  Construct the derivative.
            //
            ndp = nd - 1;

            for (i = 0; i < ndp; i++)
            {
                xdp[i] = 0.0;
            }

            for (i = 0; i < ndp; i++)
            {
                ydp[i] = (double) (i + 1) * yd_temp[i + 1];
            }
        }

        public static void dif_antideriv(int ntab, double[] xtab, double[] diftab, ref int ntab2,
                double[] xtab2, ref double[] diftab2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DIF_ANTIDERIV integrates a polynomial in divided difference form.
            //
            //  Discussion:
            //
            //    This routine uses the divided difference representation (XTAB, DIFTAB)
            //    of a polynomial to compute the divided difference representation
            //    (XTAB, ANTTAB) of the antiderivative of the polynomial.
            //
            //    The antiderivative of a polynomial P(X) is any polynomial Q(X)
            //    with the property that d/dX Q(X) = P(X).
            //
            //    This routine chooses the antiderivative whose constant term is zero.
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
            //    Input, int NTAB, the size of the input table.
            //
            //    Input, double XTAB[NTAB], the abscissas for the divided
            //    difference table.
            //
            //    Input, double DIFTAB[NTAB], the divided difference table.
            //
            //    Output, int *NTAB2, the size of the output table, which is NTAB+1.
            //
            //    Input, double XTAB2[NTAB2], the abscissas for the divided
            //    difference table for the antiderivative.
            //
            //    Output, double DIFTAB2[NTAB2], the divided difference
            //    table for the antiderivative.
            //
        {
            int i;
            double[] xtab1;
            double[] diftab1;
            //
            //  Using a temporary copy of the difference table, shift the
            //  abscissas to zero.
            //
            xtab1 = new double[ntab];
            diftab1 = new double[ntab];

            for (i = 0; i < ntab; i++)
            {
                xtab1[i] = xtab[i];
            }

            for (i = 0; i < ntab; i++)
            {
                diftab1[i] = diftab[i];
            }

            dif_shift_zero(ntab, ref xtab1, ref diftab1);

            Console.WriteLine("");
            //
            //  Append a final zero to XTAB.
            //
            ntab2 = ntab + 1;
            for (i = 0; i < ntab2; i++)
            {
                xtab2[i] = 0.0;
            }

            //
            //  Get the antiderivative of the standard form polynomial.
            //
            typeMethods.r8poly_ant_cof(ntab, diftab1, ref diftab2);
        }

        public static void dif_append(int ntab, double[] xtab, double[] diftab, double xval,
                double yval, ref int ntab2, ref double[] xtab2, ref double[] diftab2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DIF_APPEND adds a pair of data values to a divided difference table.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 September 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int NTAB, the size of the difference table.
            //
            //    Input, double XTAB[NTAB], the abscissas of the table.
            //
            //    Input, double DIFTAB[NTAB], the difference table.
            //
            //    Input, double XVAL, the data abscissa to be added to the table.
            //
            //    Input, double YVAL, the data value to be added to the table.
            //
            //    Output, int *NTAB2, the updated size of the difference table.
            //
            //    Output, double XTAB2[*NTAB2], the updated abscissas.
            //
            //    Output, double DIFTAB2[*NTAB2], the updated difference table.
            //
        {
            int i;

            ntab2 = ntab + 1;
            //
            //  Move the original data up one index.
            //
            for (i = ntab2 - 1; 1 <= i; i--)
            {
                diftab2[i] = diftab[i - 1];
                xtab2[i] = xtab[i - 1];
            }

            //
            //  Recompute the data.
            //
            xtab2[0] = xval;
            diftab2[0] = yval;

            for (i = 1; i < ntab2; i++)
            {
                diftab2[i] = (diftab2[i] - diftab2[i - 1]) / (xtab2[i] - xtab2[0]);
            }

            return;
        }

        public static void dif_basis(int ntab, double[] xtab, ref double[] diftab)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DIF_BASIS: all Lagrange basis polynomials in divided difference form.
            //
            //  Discussion:
            //
            //    The I-th Lagrange basis polynomial for a set of NTAB X values XTAB,
            //    L(I,NTAB,XTAB)(X) is a polynomial of order NTAB-1 which is zero at
            //    XTAB(J) for J not equal to I, and 1 when J is equal to I.
            //
            //    The Lagrange basis polynomials have the property that the interpolating
            //    polynomial through a set of NTAB data points (XTAB,YTAB) may be
            //    represented as
            //
            //      P(X) = Sum ( 1 <= I <= N ) YTAB(I) * L(I,NTAB,XTAB)(X)
            //
            //    Higher order interpolation at selected points may be accomplished
            //    using repeated X values, and scaled derivative values.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 September 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Carl deBoor,
            //    A Practical Guide to Splines,
            //    Springer, 2001,
            //    ISBN: 0387953663,
            //    LC: QA1.A647.v27.
            //
            //  Parameters:
            //
            //    Input, int NTAB, the number of X data points XTAB, and the number of
            //    basis polynomials to compute.
            //
            //    Input, double XTAB[NTAB], the X values upon which the Lagrange basis
            //    polynomials are to be based.
            //
            //    Output, double DIFTAB[NTAB*NTAB], points to a list of NTAB * NTAB values,
            //    the set of divided difference tables, stored as consecutive rows.
            //    Logical row I of DIFTAB contains the table for the I-th Lagrange basis
            //    polynomial.
            //
        {
            int i;
            int j;
            double[] pointer1;
            double[] pointer2;


            pointer1 = diftab;
            int p1index = 0;

            for (i = 0; i < ntab; i++)
            {
                pointer2 = pointer1;

                for (j = 0; j < ntab; j++)
                {
                    if (j == i)
                    {
                        pointer1[p1index] = 1.0;
                    }
                    else
                    {
                        pointer1[p1index] = 0.0;
                    }

                    p1index++;
                }

                Data.data_to_dif(ntab, xtab, pointer2, ref pointer2);
            }
        }

        public static void dif_basis_deriv(int nd, double[] xd, ref double[] xdp, ref double[] ddp)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DIF_BASIS_DERIV: Lagrange basis derivative difference tables.
            //
            //  Discussion:
            //
            //    Given ND points XD, a Lagrange basis polynomial L(J)(X) is associated
            //    with each point XD(J).
            //
            //    This function computes a table DDP(*,*) whose J-th column contains
            //    the difference table for the first derivative of L(J)(X).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    01 June 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Carl deBoor,
            //    A Practical Guide to Splines,
            //    Springer, 2001,
            //    ISBN: 0387953663,
            //    LC: QA1.A647.v27.
            //
            //  Parameters:
            //
            //    Input, int ND, the number of data points.
            //
            //    Input, double XD[ND], the X values upon which the 
            //    Lagrange basis polynomials are to be based.
            //
            //    Output, double XDP[ND-1], the X values upon with
            //    the derivative difference table is based.  In fact, these are
            //    all 0.
            //
            //    Output, double DDP[(ND-1)*ND], the divided difference 
            //    tables for all the Lagrange basis polynomials.  Column J of DDP
            //    contains the table for basis polynomial associated with XD(J).
            //
        {
            double[] dd;
            int i;
            int j;
            double[] yd;
            //
            //  Process the vectors one column at a time.
            //
            dd = new double[nd];
            yd = new double[nd];

            for (j = 0; j < nd; j++)
            {
                //
                //  Set the data.
                //
                for (i = 0; i < nd; i++)
                {
                    yd[i] = 0.0;
                }

                yd[j] = 1.0;
                //
                //  Compute the divided difference table.
                //
                Data.data_to_dif(nd, xd, yd, ref dd);
                //
                //  Compute the divided difference table for the derivative.
                //
                dif_deriv_table(nd, xd, dd, ref xdp, ddp, ydpIndex: +j * (nd - 1));
            }
        }

        public static void dif_basis_derivk(int nd, double[] xd, int k, ref double[] xdp, ref double[] ddp)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DIF_BASIS_DERIVK: Lagrange basis K-th derivative difference tables.
            //
            //  Discussion:
            //
            //    Given ND points XD, a Lagrange basis polynomial L(J)(X) is associated
            //    with each point XD(J).
            //
            //    This function computes a table DDP(*,*) whose J-th column contains
            //    the difference table for the K-th derivative of L(J)(X).
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    04 June 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Carl deBoor,
            //    A Practical Guide to Splines,
            //    Springer, 2001,
            //    ISBN: 0387953663,
            //    LC: QA1.A647.v27.
            //
            //  Parameters:
            //
            //    Input, int ND, the number of data points.
            //
            //    Input, double XD[ND], the X values upon which the 
            //    Lagrange basis polynomials are to be based.
            //
            //    Input, int K, the index of the derivative.
            //
            //    Output, double XDP[ND-1], the X values upon with
            //    the derivative difference table is based.  In fact, these are
            //    all 0.
            //
            //    Output, double DDP[(ND-1)*ND], the divided difference 
            //    tables for all the Lagrange basis polynomials.  Column J of DDP
            //    contains the table for basis polynomial associated with XD(J).
            //
        {
            double[] dd;
            int i;
            int j;
            double[] yd;
            //
            //  Process the vectors one column at a time.
            //
            dd = new double[nd];
            yd = new double[nd];

            for (j = 0; j < nd; j++)
            {
                //
                //  Set the data.
                //
                for (i = 0; i < nd; i++)
                {
                    yd[i] = 0.0;
                }

                yd[j] = 1.0;
                //
                //  Compute the divided difference table.
                //
                Data.data_to_dif(nd, xd, yd, ref dd);
                //
                //  Compute the divided difference table for the derivative.
                //
                dif_derivk_table(nd, xd, dd, k, xdp, ref ddp, ddkIndex: +j * (nd - k));
            }
        }

        public static void dif_basis_i(int ival, int ntab, double[] xtab, ref double[] diftab)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DIF_BASIS_I: I-th Lagrange basis polynomial in divided difference form.
            //
            //  Discussion:
            //
            //    The I-th Lagrange basis polynomial for a set of NTAB X values XTAB,
            //    L(I,NTAB,XTAB)(X) is a polynomial of order NTAB-1 which is zero at
            //    XTAB(J) for J not equal to I, and 1 when J is equal to I.
            //
            //    The Lagrange basis polynomials have the property that the interpolating
            //    polynomial through a set of NTAB data points (XTAB,YTAB) may be
            //    represented as
            //
            //      P(X) = Sum ( 1 <= I <= N ) YTAB(I) * L(I,NTAB,XTAB)(X)
            //
            //    Higher order interpolation at selected points may be accomplished
            //    using repeated X values, and scaled derivative values.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 September 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Carl deBoor,
            //    A Practical Guide to Splines,
            //    Springer, 2001,
            //    ISBN: 0387953663,
            //    LC: QA1.A647.v27.
            //
            //  Parameters:
            //
            //    Input, int IVAL, the index of the desired Lagrange basis polynomial.
            //    IVAL should be between 1 and NTAB.
            //
            //    Input, int NTAB, the number of data points XTAB.
            //
            //    Input, double XTAB[NTAB], the X values upon which the Lagrange basis
            //    polynomial is to be based.
            //
            //    Output, double DIFTAB[NTAB], the divided difference table for the IVAL-th
            //    Lagrange basis polynomial.
            //
        {
            int i;
            //
            //  Check IVAL.
            //
            if (ival < 1 || ntab < ival)
            {
                Console.WriteLine("");
                Console.WriteLine("DIF_BASIS_I - Fatal error!");
                Console.WriteLine("  IVAL must be between 1 and " + ntab + ".");
                Console.WriteLine("  but your value is " + ival + "");
                return;
            }

            //
            //  Initialize DIFTAB to Delta(I,J).
            //
            for (i = 0; i <= ntab - 1; i++)
            {
                diftab[i] = 0.0;
            }

            diftab[ival] = 1.0;
            //
            //  Compute the IVAL-th Lagrange basis polynomial.
            //
            Data.data_to_dif(ntab, xtab, diftab, ref diftab);
        }

        public static void dif_deriv_table(int nd, double[] xd, double[] yd, ref double[] xdp,
                double[] ydp, int ydpIndex = 0)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DIF_DERIV_TABLE computes the divided difference table for a derivative.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    01 June 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Carl deBoor,
            //    A Practical Guide to Splines,
            //    Springer, 2001,
            //    ISBN: 0387953663,
            //    LC: QA1.A647.v27.
            //
            //  Parameters:
            //
            //    Input, int ND, the size of the input table.
            //
            //    Input, double XD[ND], the abscissas for the divided
            //    difference table.
            //
            //    Input, double YD[ND], the divided difference table.
            //
            //    Output, double XDP[ND-1], the abscissas for the divided
            //    difference table for the derivative.
            //
            //    Output, double YDP[ND-1], the divided difference
            //    table for the derivative.
            //
        {
            int i;
            double[] xd_temp;
            double[] yd_temp;
            //
            //  Using a temporary copy of the difference table, shift the
            //  abscissas to zero.
            //
            xd_temp = new double[nd];
            yd_temp = new double[nd];

            for (i = 0; i < nd; i++)
            {
                xd_temp[i] = xd[i];
            }

            for (i = 0; i < nd; i++)
            {
                yd_temp[i] = yd[i];
            }

            dif_shift_zero(nd, ref xd_temp, ref yd_temp);
            //
            //  Construct the derivative.
            //
            for (i = 0; i < nd - 1; i++)
            {
                xdp[i] = 0.0;
            }

            for (i = 0; i < nd - 1; i++)
            {
                ydp[ydpIndex + i] = (double) (i + 1) * yd_temp[i + 1];
            }
        }

        public static void dif_derivk_table(int nd, double[] xd, double[] dd, int k,
                double[] xdk, ref double[] ddk, int ddkIndex = 0)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DIF_DERIVK_TABLE computes the divided difference table for K-th derivative.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    01 June 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Carl deBoor,
            //    A Practical Guide to Splines,
            //    Springer, 2001,
            //    ISBN: 0387953663,
            //    LC: QA1.A647.v27.
            //
            //  Parameters:
            //
            //    Input, int ND, the size of the input table.
            //
            //    Input, double XD[ND], the abscissas for the divided
            //    difference table.
            //
            //    Input, double DD[ND], the divided difference table.
            //
            //    Input, int K, the index of the derivative.  0 <= K.
            //
            //    Input, double XDK[ND-K], the abscissas for the divided
            //    difference table for the derivative.
            //
            //    Output, double DDK[NDP], the divided difference
            //    table for the derivative.
            //
        {
            double[] dd_temp;
            int i;
            int j;
            int ndk;
            double[] xd_temp;

            if (k < 0)
            {
                Console.WriteLine("");
                Console.WriteLine("DIF_DERIVK_TABLE - Fatal error!");
                Console.WriteLine("  K < 0.");
                return;
            }

            if (nd <= k)
            {
                return;
            }

            //
            //  Shift the abscissas to zero.
            //
            ndk = nd;

            xd_temp = new double[ndk];
            dd_temp = new double[ndk];

            for (i = 0; i < ndk; i++)
            {
                xd_temp[i] = xd[i];
            }

            for (i = 0; i < ndk; i++)
            {
                dd_temp[i] = dd[i];
            }

            dif_shift_zero(ndk, ref xd_temp, ref dd_temp);
            //
            //  Repeatedly differentiate.
            //
            for (j = 1; j <= k; j++)
            {
                ndk = ndk - 1;

                for (i = 0; i < ndk; i++)
                {
                    dd_temp[i] = (double) (i + 1) * dd_temp[i + 1];
                }
            }

            for (i = 0; i < ndk; i++)
            {
                ddk[i + ddkIndex] = dd_temp[i];
                xdk[i] = 0.0;
            }
        }

        public static void dif_print(int ntab, double[] xtab, double[] diftab, string title)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DIF_PRINT prints the polynomial represented by a divided difference table.
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
            //    Input, int NTAB, the dimension of the arrays DIFTAB and XTAB.
            //
            //    Input, double XTAB[NTAB], the X values for the polynomial.
            //
            //    Input, double DIFTAB[NTAB], the divided difference table
            //    for the polynomial.
            //
            //    Input, string TITLE, a title.
            //
        {
            int i;

            Console.WriteLine("");
            Console.WriteLine(title + "");
            Console.WriteLine("");
            Console.WriteLine("  p(x) =                       "
                              + diftab[0].ToString().PadLeft(14) + "");

            for (i = 1; i < ntab; i++)
            {
                Console.WriteLine("       + ( x - "
                                  + xtab[i - 1].ToString().PadLeft(10) + ") * ( "
                                  + diftab[i].ToString().PadLeft(14) + "");
            }

            string cout = "        ";
            for (i = 1; i < ntab; i++)
            {
                cout += ")";
            }

            Console.WriteLine(cout);
        }

        public static void dif_shift_x(int nd, ref double[] xd, ref double[] yd, double xv)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DIF_SHIFT_X replaces one abscissa of a divided difference table.
            //
            //  Discussion:
            //
            //    This routine shifts the representation of a divided difference polynomial
            //    by dropping the last X value in XD, and adding a new X value to the
            //    beginning of the Xd array, suitably modifying the coefficients stored
            //    in YD.
            //
            //    The representation of the polynomial is changed, but the polynomial itself
            //    should be identical.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 June 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Carl deBoor,
            //    A Practical Guide to Splines,
            //    Springer, 2001,
            //    ISBN: 0387953663,
            //    LC: QA1.A647.v27.
            //
            //  Parameters:
            //
            //    Input, int ND, the number of divided difference coefficients, and
            //    the number of entries in XD.
            //
            //    Input/output, double XD[ND], the X values used in the representation of
            //    the divided difference polynomial.  After a call to this routine, the 
            //    last entry of XD has been dropped, the other
            //    entries have shifted up one index, and XV has been inserted at the
            //    beginning of the array.
            //
            //    Input/output, double YD[ND], the divided difference coefficients
            //    corresponding to the XD array.  On output, this array has been
            //    adjusted.
            //
            //    Input, double XV, a new X value which is to be used in the representation
            //    of the polynomial.  On output, XD[0] equals XV and the representation
            //    of the polynomial has been suitably changed.
            //    Note that XV does not have to be distinct from any of the original XD
            //    values.
            //
        {
            int i;
            //
            //  Recompute the divided difference coefficients.
            //
            for (i = nd - 2; 0 <= i; i--)
            {
                yd[i] = yd[i] + (xv - xd[i]) * yd[i + 1];
            }

            //
            //  Shift the X values up one position and insert XV.
            //
            for (i = nd - 1; 0 < i; i--)
            {
                xd[i] = xd[i - 1];
            }

            xd[0] = xv;
        }

        public static void dif_shift_zero(int nd, ref double[] xd, ref double[] yd)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DIF_SHIFT_ZERO shifts a divided difference table so all abscissas are zero.
            //
            //  Discussion:
            //
            //    When the abscissas are changed, the coefficients naturally
            //    must also be changed.
            //
            //    The resulting pair (XD, YD) still represents the
            //    same polynomial, but the entries in YD are now the
            //    standard polynomial coefficients.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    01 November 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Carl deBoor,
            //    A Practical Guide to Splines,
            //    Springer, 2001,
            //    ISBN: 0387953663,
            //    LC: QA1.A647.v27.
            //
            //  Parameters:
            //
            //    Input, int ND, the length of the XD and YD arrays.
            //
            //    Input/output, double XD[ND], the X values that correspond to the
            //    divided difference table.  On output, XD contains only zeroes.
            //
            //    Input/output, double YD[ND], the divided difference table
            //    for the polynomial.  On output, YD is also
            //    the coefficient array for the standard representation
            //    of the polynomial.
            //
        {
            int i;
            int j;

            for (j = 1; j <= nd; j++)
            {
                //
                //  Recompute the divided difference coefficients.
                //
                for (i = nd - 2; 0 <= i; i--)
                {
                    yd[i] = yd[i] - xd[i] * yd[i + 1];
                }

                //
                //  Shift the XD values up one position and insert XV.
                //
                for (i = nd - 1; 0 < i; i--)
                {
                    xd[i] = xd[i - 1];
                }

                xd[0] = 0.0;
            }

            return;
        }

        public static void dif_to_r8poly(int ntab, double[] xtab, double[] diftab, ref double[] c, int diftabIndex = 0)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DIF_TO_R8POLY converts a divided difference table to a standard polynomial.
            //
            //  Discussion:
            //
            //    The vector DIFTAB, containing the divided difference polynomial
            //    coefficients is overwritten with the standard form polynomial
            //    coefficients, but the abscissa vector XTAB is unchanged.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 September 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Carl deBoor,
            //    A Practical Guide to Splines,
            //    Springer, 2001,
            //    ISBN: 0387953663,
            //    LC: QA1.A647.v27.
            //
            //  Parameters:
            //
            //    Input, int NTAB, the number of coefficients, and abscissas.
            //
            //    Input, double XTAB[NTAB], the X values used in the divided difference
            //    representation of the polynomial.
            //
            //    Input, double DIFTAB[NTAB], the divided difference table.
            //
            //    Output, double C[NTAB], the standard form polyomial coefficients.
            //    C[0] is the constant term, and C[NTAB-1] is the coefficient
            //    of X**(NTAB-1).
            //
        {
            int i;
            int j;

            for (i = 0; i < ntab; i++)
            {
                c[i] = diftab[diftabIndex + i];
            }

            //
            //  Recompute the divided difference coefficients.
            //
            for (j = 1; j <= ntab - 1; j++)
            {
                for (i = 1; i <= ntab - j; i++)
                {
                    c[ntab - i - 1] = c[ntab - i - 1] - xtab[ntab - i - j] * c[ntab - i];
                }
            }
        }

        public static double dif_val(int ntab, double[] xtab, double[] diftab, double xv, int diftabIndex = 0)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DIF_VAL evaluates a divided difference polynomial at a point.
            //
            //  Discussion:
            //
            //    DATA_TO_DIF must be called first to set up the divided difference table.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    05 September 2004
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Carl deBoor,
            //    A Practical Guide to Splines,
            //    Springer, 2001,
            //    ISBN: 0387953663,
            //    LC: QA1.A647.v27.
            //
            //  Parameters:
            //
            //    Input, integer NTAB, the number of divided difference
            //    coefficients in DIFTAB, and the number of points XTAB.
            //
            //    Input, double XTAB[NTAB], the X values upon which the
            //    divided difference polynomial is based.
            //
            //    Input, double DIFTAB[NTAB], the divided difference table.
            //
            //    Input, double XV, a value of X at which the polynomial
            //    is to be evaluated.
            //
            //    Output, double DIF_VAL, the value of the polynomial at XV.
            //
        {
            int i;
            double yv;

            yv = diftab[diftabIndex + (ntab - 1)];
            for (i = 2; i <= ntab; i++)
            {
                yv = diftab[diftabIndex + (ntab - i)] + (xv - xtab[ntab - i]) * yv;
            }

            return yv;
        }

        public static double[] dif_vals(int nd, double[] xd, double[] yd, int nv, double[] xv)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    DIF_VALS evaluates a divided difference polynomial at a set of points.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    23 May 2011
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Reference:
            //
            //    Carl deBoor,
            //    A Practical Guide to Splines,
            //    Springer, 2001,
            //    ISBN: 0387953663,
            //    LC: QA1.A647.v27.
            //
            //  Parameters:
            //
            //    Input, int ND, the order of the difference table.
            //
            //    Input, double XD[ND], the X values of the difference table.
            //
            //    Input, double YD[ND], the divided differences.
            //
            //    Input, int NV, the number of evaluation points.
            //
            //    Input, double XV[NV], the evaluation points.
            //
            //    Output, double DIF_VALS[NV], the value of the divided difference
            //    polynomial at the evaluation points.
            //
        {
            int i;
            int j;
            double[] yv;

            yv = new double[nv];

            for (j = 0; j < nv; j++)
            {
                yv[j] = yd[nd - 1];
                for (i = 2; i <= nd; i++)
                {
                    yv[j] = yd[nd - i] + (xv[j] - xd[nd - i]) * yv[j];
                }
            }

            return yv;
        }
    }
}