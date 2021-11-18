using System;
using Burkardt.Composition;
using Burkardt.MonomialNS;
using Burkardt.Types;
using Burkardt.PolynomialNS;

namespace Burkardt.Interpolation;

using Polynomial = Polynomial;
public static class MultivariateLagrange
{
        
    public static double[] interpolant_value(int d, int r, int pn, int[] po, double[] pc,
            int[] pe, double[] pd, int ni, double[] xi )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INTERPOLANT_VALUE evaluates a Lagrange interpolant.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    29 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int D, the spatial dimension.
        //
        //    Input, int R, the maximum number of terms in a polynomial.
        //
        //    Input, int PN, the number of polynomials.
        //
        //    Input, int PO[PN], the "order" of the polynomials.
        //
        //    Input, double PC[PN*R], the coefficients of the polynomial.
        //
        //    Input, int PE[PN*R], the indices of the exponents of the polynomial.
        //
        //    Input, double PD[PN], the coefficient of each polynomial.  
        //    For a Lagrange interpolant, this is the data value at each Lagrange point.
        //
        //    Input, int NI, the number of interpolant evaluation points.
        //
        //    Input, double XI[D*NI], the coordinates of the interpolation evaluation points.
        //
        //    Output, double INTERPOLANT_VALUE[NI], the value of the interpolant at XI.
        //
    {
        double[] c;
        int[] e;
        int i;
        int j;
        int oj;
        double[] value;
        double[] yi;

        yi = new double[ni];

        for (i = 0; i < ni; i++)
        {
            yi[i] = 0.0;
        }

        c = new double[r];
        e = new int[r];

        for (j = 0; j < pn; j++)
        {
            oj = po[j];
            for (i = 0; i < oj; i++)
            {
                c[i] = pc[j + i * pn];
                e[i] = pe[j + i * pn];
            }

            value = Polynomial.polynomial_value(d, oj, c, e, ni, xi);
            for (i = 0; i < ni; i++)
            {
                yi[i] += pd[j] * value[i];
            }

        }

        return yi;
    }

    public static void lagrange_complete(int d, int n, int r, int nd, double[] xd, ref int[] po,
            ref double[] pc, ref int[] pe )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGRANGE_COMPLETE: Complete Lagrange polynomial basis from data.
        //
        //  Discussion:
        //
        //    This function represents algorithm 4.1 in the reference.
        //
        //    This function is given XD, a set of ND distinct data points in a 
        //    D dimensional space, and returns information defining a set of 
        //    ND Lagrange polynomials L(i)(X) with the property that:
        //
        //      L(i)(XD(j)) = delta(i,j)
        //
        //    In order for this computation to be carried out, it is necessary that
        //    ND, the number of data points, is equal to R, the dimension of the 
        //    space of polynomials in D dimensions and total degree N or less, that is:
        //
        //      ND = R = Choose ( N + D, N )
        //
        //    There will be ND polynomials returned.  Each polynomial can have
        //    as many as R coefficients.
        //
        //    Each polynomial is given as a vector, with each entry corresponding
        //    to a nonzero coefficient.  In particular, for polynomial L(i)(X):
        //
        //      PO(i) is the order, that is, the number of nonzero coefficients;
        //      PC(i,j), for 1 <= j <= PO(i), is the coefficient of the J-th term.
        //      PE(i,j), for 1 <= j <= PO(i), encodes the exponents of the J-th term.
        //
        //    The exponent codes are a compact way of recording the exponent vector
        //    associated with each monomial.  If PE(i,j) = k, then the corresponding
        //    vector of D exponents can be determined by:
        //
        //      E = mono_unrank_grlex ( D, k );
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 February 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Tomas Sauer, Yuan Xu,
        //    On multivariate Lagrange interpolation,
        //    Mathematics of Computation,
        //    Volume 64, Number 211, July 1995, pages 1147-1170.
        //
        //  Parameters:
        //
        //    Input, int D, the spatial dimension.
        //
        //    Input, int N, the maximum total degree.
        //
        //    Input, int R, the number of monomials in D dimensions 
        //    of total degree N or less.
        //
        //    Input, int ND, the number of data points.
        //    This function requires that the ND is equal to R.
        //
        //    Input, double XD[D*ND], the data points, which must be distinct.
        //
        //    Output, int PO[ND], the order (number of nonzero coefficients),
        //    for the Lagrange basis polynomials.
        //
        //    Output, double PC[ND*R], the coefficients for the 
        //    Lagrange basis polynomials.
        //
        //    Output, int PE[ND*R], the  exponent indices for the 
        //    Lagrange basis polynomials.
        //
    {
        double[] c;
        double[] cj;
        double[] ck;
        double d_max = 0;
        double d_min = 0;
        double d_tol;
        int[] e;
        int[] ej;
        int[] ek;
        int i;
        int j;
        int k;
        int l;
        int o;
        int oj;
        int ok = 0;
        double[] qc;
        int[] qe;
        int[] qo;
        double[] value = new double[1];
        //
        //  Verify that R is correct.
        //
        if (r != Monomial.mono_upto_enum(d, n))
        {
            Console.WriteLine("");
            Console.WriteLine("LAGRANGE_COMPLETE - Fatal error!");
            Console.WriteLine("  The value R is not correct.");
            return;
        }

        if (r != nd)
        {
            Console.WriteLine("");
            Console.WriteLine("LAGRANGE_COMPLETE - Fatal error!");
            Console.WriteLine("  The value R = " + r + "");
            Console.WriteLine("  does not equal ND = " + nd + "");
            return;
        }

        //
        //  Verify that the points are sufficiently distinct.
        //
        typeMethods.r8col_separation(d, nd, xd, ref d_min, ref d_max);
        d_tol = Math.Sqrt(typeMethods.r8_epsilon());

        if (d_min < d_tol)
        {
            Console.WriteLine("");
            Console.WriteLine("LAGRANGE_COMPLETE - Fatal error!");
            Console.WriteLine("  Some points are too close!");
            Console.WriteLine("  Minimum data point separation is = " + d_min + "");
            return;
        }

        //
        //  Make some work space.
        //
        c = new double[r];
        cj = new double[r];
        ck = new double[r];
        e = new int[r];
        ej = new int[r];
        ek = new int[r];
        //
        //  Initialize the polynomials Q, which span the space of
        //  N-th degree polynomials.
        //
        //  Default option: 
        //  * all ND-dimensional monomials of degree N or less.
        //    in 2D, this might be 1, x, y, x^2, xy, y^2, ...
        //
        qo = new int[r];
        qc = new double[r * r];
        qe = new int[r * r];

        for (k = 0; k < r; k++)
        {
            qo[k] = 1;
            qc[k + 0 * r] = 1.0;
            qe[k + 0 * r] = k + 1;
            for (j = 1; j < r; j++)
            {
                qc[k + j * r] = 0.0;
                qe[k + j * r] = 0;
            }
        }

        //
        //  Now set up the P polynomials.
        //
        for (k = 0; k < r; k++)
        {
            po[k] = 0;
            for (j = 0; j < r; j++)
            {
                pc[k + j * r] = 0.0;
                pe[k + j * r] = 0;
            }
        }

        for (k = 0; k < nd; k++)
        {
            //
            //  Find the first polynomial Q(K:R)(X) which is nonzero at X(K).
            //
            i = r + 1;

            for (j = k; j < r; j++)
            {
                o = qo[j];
                for (l = 0; l < o; l++)
                {
                    c[l] = qc[j + l * r];
                    e[l] = qe[j + l * r];
                }

                value = Polynomial.polynomial_value(d, o, c, e, 1, xd, xIndex: + k * d);

                if (value[0] != 0.0)
                {
                    i = j;
                    break;
                }
            }

            if (i == r + 1)
            {
                Console.WriteLine("");
                Console.WriteLine("LAGRANGE_COMPLETE - Fatal error!");
                Console.WriteLine("  I = R+1.");
                return;
            }

            //
            //  Define P(K)(X) = Q(I)(X) / Q(I)(X(k)
            //
            o = qo[i];
            po[k] = qo[i];
            for (l = 0; l < o; l++)
            {
                pc[k + l * r] = qc[i + l * r] / value[0];
                pe[k + l * r] = qe[i + l * r];
            }

            //
            //  Modify P(1:k-1)(X).
            //
            for (j = 0; j < k; j++)
            {
                oj = po[j];
                for (l = 0; l < oj; l++)
                {
                    cj[l] = pc[j + l * r];
                    ej[l] = pe[j + l * r];
                }

                value = Polynomial.polynomial_value(d, oj, cj, ej, 1, xd, xIndex: + k * d);

                ok = po[k];
                for (l = 0; l < ok; l++)
                {
                    ck[l] = pc[k + l * r];
                    ek[l] = pe[k + l * r];
                }

                Polynomial.polynomial_axpy(-value[0], ok, ck, ek, oj, cj, ej, ref o, ref c, ref e);

                po[j] = o;
                for (l = 0; l < o; l++)
                {
                    pc[j + l * r] = c[l];
                    pe[j + l * r] = e[l];
                }

            }

            //
            //  Modify Q(I:downto:K+1)
            //
            for (j = i; k < j; j--)
            {
                oj = qo[j - 1];
                for (l = 0; l < oj; l++)
                {
                    cj[l] = qc[j - 1 + l * r];
                    ej[l] = qe[j - 1 + l * r];
                }

                value = Polynomial.polynomial_value(d, oj, cj, ej, 1, xd, xIndex: + k * d);

                ok = po[k];
                for (l = 0; l < ok; l++)
                {
                    ck[l] = pc[k + l * r];
                    ek[l] = pe[k + l * r];
                }

                Polynomial.polynomial_axpy(-value[0], ok, ck, ek, oj, cj, ej, ref o, ref c, ref e);
                    
                qo[j] = o;
                for (l = 0; l < o; l++)
                {
                    qc[j + l * r] = c[l];
                    qe[j + l * r] = e[l];
                }
            }

            //
            //  Modify Q(I+1:R)
            //
            for (j = i + 1; j < r; j++)
            {
                oj = qo[j];
                for (l = 0; l < oj; l++)
                {
                    cj[l] = qc[j + l * r];
                    ej[l] = qe[j + l * r];
                }

                value = Polynomial.polynomial_value(d, oj, cj, ej, 1, xd, xIndex: + k * d);

                ok = po[k];
                for (l = 0; l < ok; l++)
                {
                    ck[l] = pc[k + l * r];
                    ek[l] = pe[k + l * r];
                }

                Polynomial.polynomial_axpy(-value[0], ok, ck, ek, oj, cj, ej, ref o, ref c, ref e);

                qo[j] = o;
                for (l = 0; l < o; l++)
                {
                    qc[j + l * r] = c[l];
                    qe[j + l * r] = e[l];
                }
            }
        }

        //
        //  Get rid of tiny coefficients.
        //
        for (i = 0; i < nd; i++)
        {
            oj = po[i];
            for (l = 0; l < oj; l++)
            {
                cj[l] = pc[i + l * r];
                ej[l] = pe[i + l * r];
            }

            Polynomial.polynomial_compress(oj, cj, ej, ref ok, ref ck, ref ek);

            po[i] = ok;
            for (l = 0; l < ok; l++)
            {
                pc[i + l * r] = ck[l];
                pe[i + l * r] = ek[l];
            }
        }
    }

    public static void lagrange_complete2(int d, int n, int r, int nd, double[] xd, ref int[] po,
            ref double[] pc, ref int[] pe )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGRANGE_COMPLETE2: Complete Lagrange polynomial basis from data.
        //
        //  Discussion:
        //
        //    This function represents algorithm 4.1 in the reference,
        //    with the further modification that a form of "pivoting" is used
        //    to select the next polynomial as the one with maximum absolute
        //    value at the current node.
        //
        //    This function is given XD, a set of ND distinct data points in a 
        //    D dimensional space, and returns information defining a set of 
        //    ND Lagrange polynomials L(i)(X) with the property that:
        //
        //      L(i)(XD(j)) = delta(i,j)
        //
        //    In order for this computation to be carried out, it is necessary that
        //    ND, the number of data points, is equal to R, the dimension of the 
        //    space of polynomials in D dimensions and total degree N or less, that is:
        //
        //      ND = R = Choose ( N + D, N )
        //
        //    There will be ND polynomials returned.  Each polynomial can have
        //    as many as R coefficients.
        //
        //    Each polynomial is given as a vector, with each entry corresponding
        //    to a nonzero coefficient.  In particular, for polynomial L(i)(X):
        //
        //      PO(i) is the order, that is, the number of nonzero coefficients;
        //      PC(i,j), for 1 <= j <= PO(i), is the coefficient of the J-th term.
        //      PE(i,j), for 1 <= j <= PO(i), encodes the exponents of the J-th term.
        //
        //    The exponent codes are a compact way of recording the exponent vector
        //    associated with each monomial.  If PE(i,j) = k, then the corresponding
        //    vector of D exponents can be determined by:
        //
        //      E = mono_unrank_grlex ( D, k );
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 February 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Tomas Sauer, Yuan Xu,
        //    On multivariate Lagrange interpolation,
        //    Mathematics of Computation,
        //    Volume 64, Number 211, July 1995, pages 1147-1170.
        //
        //  Parameters:
        //
        //    Input, int D, the spatial dimension.
        //
        //    Input, int N, the maximum total degree.
        //
        //    Input, int R, the number of monomials in D dimensions 
        //    of total degree N or less.
        //
        //    Input, int ND, the number of data points.
        //    This function requires that the ND is equal to R.
        //
        //    Input, double XD[D*ND], the data points, which must be distinct.
        //
        //    Output, int PO[ND], the order (number of nonzero coefficients),
        //    for the Lagrange basis polynomials.
        //
        //    Output, double PC[ND*R], the coefficients for the 
        //    Lagrange basis polynomials.
        //
        //    Output, int PE[ND*R], the  exponent indices for the 
        //    Lagrange basis polynomials.
        //
    {
        double[] c;
        double[] cj;
        double[] ck;
        double d_max = 0;
        double d_min = 0;
        double d_tol;
        int[] e;
        int[] ej;
        int[] ek;
        int i;
        int j;
        int k;
        int l;
        int o;
        int oj;
        int ok = 0;
        double[] qc;
        int[] qe;
        int[] qo;
        double[] value = new double[1];
        double value_max;
        //
        //  Verify that R is correct.
        //
        if (r != Monomial.mono_upto_enum(d, n))
        {
            Console.WriteLine("");
            Console.WriteLine("LAGRANGE_COMPLETE2 - Fatal error!");
            Console.WriteLine("  The value R is not correct.");
            return;
        }

        if (r != nd)
        {
            Console.WriteLine("");
            Console.WriteLine("LAGRANGE_COMPLETE2 - Fatal error!");
            Console.WriteLine("  The value R = " + r + "");
            Console.WriteLine("  does not equal ND = " + nd + "");
            return;
        }

        //
        //  Verify that the points are sufficiently distinct.
        //
        typeMethods.r8col_separation(d, nd, xd, ref d_min, ref d_max);
        d_tol = Math.Sqrt(typeMethods.r8_epsilon());

        if (d_min < d_tol)
        {
            Console.WriteLine("");
            Console.WriteLine("LAGRANGE_COMPLETE2 - Fatal error!");
            Console.WriteLine("  Some points are too close!");
            Console.WriteLine("  Minimum data point separation is = " + d_min + "");
            return;
        }

        //
        //  Make some work space.
        //
        c = new double[r];
        cj = new double[r];
        ck = new double[r];
        e = new int[r];
        ej = new int[r];
        ek = new int[r];
        //
        //  Initialize the polynomials Q, which span the space of
        //  N-th degree polynomials.
        //
        //  Default option: 
        //  * all ND-dimensional monomials of degree N or less.
        //    in 2D, this might be 1, x, y, x^2, xy, y^2, ...
        //
        qo = new int[r];
        qc = new double[r * r];
        qe = new int[r * r];

        for (k = 0; k < r; k++)
        {
            qo[k] = 1;
            qc[k + 0 * r] = 1.0;
            qe[k + 0 * r] = k + 1;
            for (j = 1; j < r; j++)
            {
                qc[k + j * r] = 0.0;
                qe[k + j * r] = 0;
            }
        }

        //
        //  Now set up the P polynomials.
        //
        for (k = 0; k < r; k++)
        {
            po[k] = 0;
            for (j = 0; j < r; j++)
            {
                pc[k + j * r] = 0.0;
                pe[k + j * r] = 0;
            }
        }

        for (k = 0; k < nd; k++)
        {
            //
            //  Find the first polynomial Q(K:R)(X) which is nonzero at X(K).
            //
            i = r + 1;
            value_max = 0.0;

            for (j = k; j < r; j++)
            {
                o = qo[j];
                for (l = 0; l < o; l++)
                {
                    c[l] = qc[j + l * r];
                    e[l] = qe[j + l * r];
                }

                value = Polynomial.polynomial_value(d, o, c, e, 1, xd, xIndex: + k * d);

                if (Math.Abs(value_max) <= Math.Abs(value[0]))
                {
                    i = j;
                    value_max = value[0];
                }

            }

            if (i == r + 1)
            {
                Console.WriteLine("");
                Console.WriteLine("LAGRANGE_COMPLETE2 - Fatal error!");
                Console.WriteLine("  I = R+1.");
                return;
            }

            value[0] = value_max;
            //
            //  Define P(K)(X) = Q(I)(X) / Q(I)(X(k)
            //
            o = qo[i];
            po[k] = qo[i];
            for (l = 0; l < o; l++)
            {
                pc[k + l * r] = qc[i + l * r] / value[0];
                pe[k + l * r] = qe[i + l * r];
            }

            //
            //  Modify P(1:k-1)(X).
            //
            for (j = 0; j < k; j++)
            {
                oj = po[j];
                for (l = 0; l < oj; l++)
                {
                    cj[l] = pc[j + l * r];
                    ej[l] = pe[j + l * r];
                }

                value = Polynomial.polynomial_value(d, oj, cj, ej, 1, xd, xIndex: + k * d);

                ok = po[k];
                for (l = 0; l < ok; l++)
                {
                    ck[l] = pc[k + l * r];
                    ek[l] = pe[k + l * r];
                }

                Polynomial.polynomial_axpy(-value[0], ok, ck, ek, oj, cj, ej, ref o, ref c, ref e);

                po[j] = o;
                for (l = 0; l < o; l++)
                {
                    pc[j + l * r] = c[l];
                    pe[j + l * r] = e[l];
                }

            }

            //
            //  Modify Q(I:downto:K+1)
            //
            for (j = i; k < j; j--)
            {
                oj = qo[j - 1];
                for (l = 0; l < oj; l++)
                {
                    cj[l] = qc[j - 1 + l * r];
                    ej[l] = qe[j - 1 + l * r];
                }

                value = Polynomial.polynomial_value(d, oj, cj, ej, 1, xd, xIndex: + k * d);

                ok = po[k];
                for (l = 0; l < ok; l++)
                {
                    ck[l] = pc[k + l * r];
                    ek[l] = pe[k + l * r];
                }

                Polynomial.polynomial_axpy(-value[0], ok, ck, ek, oj, cj, ej, ref o, ref c, ref e);

                qo[j] = o;
                for (l = 0; l < o; l++)
                {
                    qc[j + l * r] = c[l];
                    qe[j + l * r] = e[l];
                }
            }

            //
            //  Modify Q(I+1:R)
            //
            for (j = i + 1; j < r; j++)
            {
                oj = qo[j];
                for (l = 0; l < oj; l++)
                {
                    cj[l] = qc[j + l * r];
                    ej[l] = qe[j + l * r];
                }

                value = Polynomial.polynomial_value(d, oj, cj, ej, 1, xd, xIndex: + k * d);

                ok = po[k];
                for (l = 0; l < ok; l++)
                {
                    ck[l] = pc[k + l * r];
                    ek[l] = pe[k + l * r];
                }

                Polynomial.polynomial_axpy(-value[0], ok, ck, ek, oj, cj, ej, ref o, ref c, ref e);

                qo[j] = o;
                for (l = 0; l < o; l++)
                {
                    qc[j + l * r] = c[l];
                    qe[j + l * r] = e[l];
                }
            }
        }

        //
        //  Get rid of tiny coefficients.
        //
        for (i = 0; i < nd; i++)
        {
            oj = po[i];
            for (l = 0; l < oj; l++)
            {
                cj[l] = pc[i + l * r];
                ej[l] = pe[i + l * r];
            }

            Polynomial.polynomial_compress(oj, cj, ej, ref ok, ref ck, ref ek);

            po[i] = ok;
            for (l = 0; l < ok; l++)
            {
                pc[i + l * r] = ck[l];
                pe[i + l * r] = ek[l];
            }
        }
    }

    public static void lagrange_partial(int d, int n, int r, int nd, double[] xd, ref int[] po,
            ref double[] pc, ref int[] pe )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGRANGE_PARTIAL: Partial Lagrange polynomial basis from data.
        //
        //  Discussion:
        //
        //    This function represents algorithm 4.1 in the reference,
        //    modified for the case where the number of data points is less
        //    than the dimension of the desired polynomial space.
        //
        //    This function is given XD, a set of ND distinct data points in a 
        //    D dimensional space, and returns information defining a set of 
        //    ND Lagrange polynomials L(i)(X) with the property that:
        //
        //      L(i)(XD(j)) = delta(i,j)
        //
        //    This function is used in cases where ND, the number of data points, 
        //    is less than or equal to R, the dimension of the space of polynomials 
        //    in D dimensions and total degree N or less, that is:
        //
        //      ND <= R = Choose ( N + D, N )
        //
        //    There will be ND polynomials returned.  Each polynomial can have
        //    as many as R coefficients.
        //
        //    Each polynomial is given as a vector, with each entry corresponding
        //    to a nonzero coefficient.  In particular, for polynomial L(i)(X):
        //
        //      PO(i) is the order, that is, the number of nonzero coefficients;
        //      PC(i,j), for 1 <= j <= PO(i), is the coefficient of the J-th term.
        //      PE(i,j), for 1 <= j <= PO(i), encodes the exponents of the J-th term.
        //
        //    The exponent codes are a compact way of recording the exponent vector
        //    associated with each monomial.  If PE(i,j) = k, then the corresponding
        //    vector of D exponents can be determined by:
        //
        //      E = mono_unrank_grlex ( D, k );
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 February 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Tomas Sauer, Yuan Xu,
        //    On multivariate Lagrange interpolation,
        //    Mathematics of Computation,
        //    Volume 64, Number 211, July 1995, pages 1147-1170.
        //
        //  Parameters:
        //
        //    Input, int D, the spatial dimension.
        //
        //    Input, int N, the maximum total degree.
        //
        //    Input, int R, the number of monomials in D dimensions 
        //    of total degree N or less.
        //
        //    Input, int ND, the number of data points.
        //    It must be the case that ND <= R.
        //
        //    Input, double XD[D*ND], the data points, which must be distinct.
        //
        //    Output, int PO[ND], the order (number of nonzero coefficients),
        //    for the Lagrange basis polynomials.
        //
        //    Output, double PC[ND*R], the coefficients for the 
        //    Lagrange basis polynomials.
        //
        //    Output, int PE[ND*R], the  exponent indices for the 
        //    Lagrange basis polynomials.
        //
    {
        double[] c;
        double[] cj;
        double[] ck;
        double d_max = 0;
        double d_min = 0;
        double d_tol;
        int[] e;
        int[] ej;
        int[] ek;
        int i;
        int j;
        int k;
        int l;
        int o;
        int oj;
        int ok = 0;
        double[] qc;
        int[] qe;
        int[] qo;
        double[] value = new double[1];
        //
        //  Verify that R is correct.
        //
        if (r != Monomial.mono_upto_enum(d, n))
        {
            Console.WriteLine("");
            Console.WriteLine("LAGRANGE_PARTIAL - Fatal error!");
            Console.WriteLine("  The value R is not correct.");
            return;
        }

        if (r < nd)
        {
            Console.WriteLine("");
            Console.WriteLine("LAGRANGE_PARTIAL - Fatal error!");
            Console.WriteLine("  The value R = " + r + "");
            Console.WriteLine("  is less than ND = " + nd + "");
            return;
        }

        //
        //  Verify that the points are sufficiently distinct.
        //
        typeMethods.r8col_separation(d, nd, xd, ref d_min, ref d_max);
        d_tol = Math.Sqrt(typeMethods.r8_epsilon());

        if (d_min < d_tol)
        {
            Console.WriteLine("");
            Console.WriteLine("LAGRANGE_PARTIAL - Fatal error!");
            Console.WriteLine("  Some points are too close!");
            Console.WriteLine("  Minimum data point separation is = " + d_min + "");
            return;
        }

        //
        //  Make some work space.
        //
        c = new double[r];
        cj = new double[r];
        ck = new double[r];
        e = new int[r];
        ej = new int[r];
        ek = new int[r];
        //
        //  Initialize the polynomials Q, which span the space of
        //  N-th degree polynomials.
        //
        //  Default option: 
        //  * all ND-dimensional monomials of degree N or less.
        //    in 2D, this might be 1, x, y, x^2, xy, y^2, ...
        //
        qo = new int[r];
        qc = new double[r * r];
        qe = new int[r * r];

        for (k = 0; k < r; k++)
        {
            qo[k] = 1;
            qc[k + 0 * r] = 1.0;
            qe[k + 0 * r] = k + 1;
            for (j = 1; j < r; j++)
            {
                qc[k + j * r] = 0.0;
                qe[k + j * r] = 0;
            }
        }

        //
        //  Now set up the P polynomials.
        //
        for (k = 0; k < nd; k++)
        {
            po[k] = 0;
            for (j = 0; j < r; j++)
            {
                pc[k + j * nd] = 0.0;
                pe[k + j * nd] = 0;
            }
        }

        for (k = 0; k < nd; k++)
        {
            //
            //  Find the first polynomial Q(K:R)(X) which is nonzero at X(K).
            //
            i = r + 1;

            for (j = k; j < r; j++)
            {
                o = qo[j];
                for (l = 0; l < o; l++)
                {
                    c[l] = qc[j + l * r];
                    e[l] = qe[j + l * r];
                }

                value = Polynomial.polynomial_value(d, o, c, e, 1, xd, xIndex: + k * d);

                if (value[0] != 0.0)
                {
                    i = j;
                    break;
                }

            }

            if (i == r + 1)
            {
                Console.WriteLine("");
                Console.WriteLine("LAGRANGE_PARTIAL - Fatal error!");
                Console.WriteLine("  I = R+1.");
                return;
            }

            //
            //  Define P(K)(X) = Q(I)(X) / Q(I)(X(k)
            //
            o = qo[i];
            po[k] = qo[i];
            for (l = 0; l < o; l++)
            {
                pc[k + l * nd] = qc[i + l * r] / value[0];
                pe[k + l * nd] = qe[i + l * r];
            }

            //
            //  Modify P(1:k-1)(X).
            //
            for (j = 0; j < k; j++)
            {
                oj = po[j];
                for (l = 0; l < oj; l++)
                {
                    cj[l] = pc[j + l * nd];
                    ej[l] = pe[j + l * nd];
                }

                value = Polynomial.polynomial_value(d, oj, cj, ej, 1, xd, xIndex: + k * d);

                ok = po[k];
                for (l = 0; l < ok; l++)
                {
                    ck[l] = pc[k + l * nd];
                    ek[l] = pe[k + l * nd];
                }

                Polynomial.polynomial_axpy(-value[0], ok, ck, ek, oj, cj, ej, ref o, ref c, ref e);

                po[j] = o;
                for (l = 0; l < o; l++)
                {
                    pc[j + l * nd] = c[l];
                    pe[j + l * nd] = e[l];
                }

            }

            //
            //  Modify Q(I:downto:K+1)
            //
            for (j = i; k < j; j--)
            {
                oj = qo[j - 1];
                for (l = 0; l < oj; l++)
                {
                    cj[l] = qc[j - 1 + l * r];
                    ej[l] = qe[j - 1 + l * r];
                }

                value = Polynomial.polynomial_value(d, oj, cj, ej, 1, xd, xIndex: + k * d);

                ok = po[k];
                for (l = 0; l < ok; l++)
                {
                    ck[l] = pc[k + l * nd];
                    ek[l] = pe[k + l * nd];
                }

                Polynomial.polynomial_axpy(-value[0], ok, ck, ek, oj, cj, ej, ref o, ref c, ref e);
                    
                qo[j] = o;
                for (l = 0; l < o; l++)
                {
                    qc[j + l * r] = c[l];
                    qe[j + l * r] = e[l];
                }
            }

            //
            //  Modify Q(I+1:R)
            //
            for (j = i + 1; j < r; j++)
            {
                oj = qo[j];
                for (l = 0; l < oj; l++)
                {
                    cj[l] = qc[j + l * r];
                    ej[l] = qe[j + l * r];
                }

                value = Polynomial.polynomial_value(d, oj, cj, ej, 1, xd, xIndex: + k * d);

                ok = po[k];
                for (l = 0; l < ok; l++)
                {
                    ck[l] = pc[k + l * nd];
                    ek[l] = pe[k + l * nd];
                }

                Polynomial.polynomial_axpy(-value[0], ok, ck, ek, oj, cj, ej, ref o, ref c, ref e);

                qo[j] = o;
                for (l = 0; l < o; l++)
                {
                    qc[j + l * r] = c[l];
                    qe[j + l * r] = e[l];
                }
            }
        }

        //
        //  Get rid of tiny coefficients.
        //
        for (i = 0; i < nd; i++)
        {
            oj = po[i];
            for (l = 0; l < oj; l++)
            {
                cj[l] = pc[i + l * nd];
                ej[l] = pe[i + l * nd];
            }

            Polynomial.polynomial_compress(oj, cj, ej, ref ok, ref ck, ref ek);

            po[i] = ok;
            for (l = 0; l < ok; l++)
            {
                pc[i + l * nd] = ck[l];
                pe[i + l * nd] = ek[l];
            }
        }
    }

    public static void lagrange_partial2(int d, int n, int r, int nd, double[] xd, ref int[] po,
            ref double[] pc, ref int[] pe )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGRANGE_PARTIAL2: Partial Lagrange polynomial basis from data.
        //
        //  Discussion:
        //
        //    This function represents algorithm 4.1 in the reference,
        //    modified for the case where the number of data points is less
        //    than the dimension of the desired polynomial space,
        //    with the further modification that a form of "pivoting" is used
        //    to select the next polynomial as the one with maximum absolute
        //    value at the current node.
        //
        //    This function is given XD, a set of ND distinct data points in a 
        //    D dimensional space, and returns information defining a set of 
        //    ND Lagrange polynomials L(i)(X) with the property that:
        //
        //      L(i)(XD(j)) = delta(i,j)
        //
        //    This function is used in cases where ND, the number of data points, 
        //    is less than or equal to R, the dimension of the space of polynomials 
        //    in D dimensions and total degree N or less, that is:
        //
        //      ND <= R = Choose ( N + D, N )
        //
        //    There will be ND polynomials returned.  Each polynomial can have
        //    as many as R coefficients.
        //
        //    Each polynomial is given as a vector, with each entry corresponding
        //    to a nonzero coefficient.  In particular, for polynomial L(i)(X):
        //
        //      PO(i) is the order, that is, the number of nonzero coefficients;
        //      PC(i,j), for 1 <= j <= PO(i), is the coefficient of the J-th term.
        //      PE(i,j), for 1 <= j <= PO(i), encodes the exponents of the J-th term.
        //
        //    The exponent codes are a compact way of recording the exponent vector
        //    associated with each monomial.  If PE(i,j) = k, then the corresponding
        //    vector of D exponents can be determined by:
        //
        //      E = mono_unrank_grlex ( D, k );
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 February 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Tomas Sauer, Yuan Xu,
        //    On multivariate Lagrange interpolation,
        //    Mathematics of Computation,
        //    Volume 64, Number 211, July 1995, pages 1147-1170.
        //
        //  Parameters:
        //
        //    Input, int D, the spatial dimension.
        //
        //    Input, int N, the maximum total degree.
        //
        //    Input, int R, the number of monomials in D dimensions 
        //    of total degree N or less.
        //
        //    Input, int ND, the number of data points.
        //    It must be the case that ND <= R.
        //
        //    Input, double XD[D*ND], the data points, which must be distinct.
        //
        //    Output, int PO[ND], the order (number of nonzero coefficients),
        //    for the Lagrange basis polynomials.
        //
        //    Output, double PC[ND*R], the coefficients for the 
        //    Lagrange basis polynomials.
        //
        //    Output, int PE[ND*R], the  exponent indices for the 
        //    Lagrange basis polynomials.
        //
    {
        double[] c;
        double[] cj;
        double[] ck;
        double d_max = 0;
        double d_min = 0;
        double d_tol;
        int[] e;
        int[] ej;
        int[] ek;
        int i;
        int j;
        int k;
        int l;
        int o;
        int oj;
        int ok = 0;
        double[] qc;
        int[] qe;
        int[] qo;
        double[] value;
        double value_max;
        //
        //  Verify that R is correct.
        //
        if (r != Monomial.mono_upto_enum(d, n))
        {
            Console.WriteLine("");
            Console.WriteLine("LAGRANGE_PARTIAL2 - Fatal error!");
            Console.WriteLine("  The value R is not correct.");
            return;
        }

        if (r < nd)
        {
            Console.WriteLine("");
            Console.WriteLine("LAGRANGE_PARTIAL2 - Fatal error!");
            Console.WriteLine("  The value R = " + r + "");
            Console.WriteLine("  is less than ND = " + nd + "");
            return;
        }

        //
        //  Verify that the points are sufficiently distinct.
        //
        typeMethods.r8col_separation(d, nd, xd, ref d_min, ref d_max);
        d_tol = Math.Sqrt(typeMethods.r8_epsilon());

        if (d_min < d_tol)
        {
            Console.WriteLine("");
            Console.WriteLine("LAGRANGE_PARTIAL2 - Fatal error!");
            Console.WriteLine("  Some points are too close!");
            Console.WriteLine("  Minimum data point separation is = " + d_min + "");
            return;
        }

        //
        //  Make some work space.
        //
        c = new double[r];
        cj = new double[r];
        ck = new double[r];
        e = new int[r];
        ej = new int[r];
        ek = new int[r];
        //
        //  Initialize the polynomials Q, which span the space of
        //  N-th degree polynomials.
        //
        //  Default option: 
        //  * all ND-dimensional monomials of degree N or less.
        //    in 2D, this might be 1, x, y, x^2, xy, y^2, ...
        //
        qo = new int[r];
        qc = new double[r * r];
        qe = new int[r * r];

        for (k = 0; k < r; k++)
        {
            qo[k] = 1;
            qc[k + 0 * r] = 1.0;
            qe[k + 0 * r] = k + 1;
            for (j = 1; j < r; j++)
            {
                qc[k + j * r] = 0.0;
                qe[k + j * r] = 0;
            }
        }

        //
        //  Now set up the P polynomials.
        //
        for (k = 0; k < nd; k++)
        {
            po[k] = 0;
            for (j = 0; j < r; j++)
            {
                pc[k + j * nd] = 0.0;
                pe[k + j * nd] = 0;
            }
        }

        for (k = 0; k < nd; k++)
        {
            //
            //  Find the first polynomial Q(K:R)(X) which is nonzero at X(K).
            //
            i = r + 1;
            value_max = 0.0;

            for (j = k; j < r; j++)
            {
                o = qo[j];
                for (l = 0; l < o; l++)
                {
                    c[l] = qc[j + l * r];
                    e[l] = qe[j + l * r];
                }

                value = Polynomial.polynomial_value(d, o, c, e, 1, xd, xIndex: + k * d);

                if (Math.Abs(value_max) <= Math.Abs(value[0]))
                {
                    i = j;
                    value_max = value[0];
                }
            }

            if (i == r + 1)
            {
                Console.WriteLine("");
                Console.WriteLine("LAGRANGE_PARTIAL2 - Fatal error!");
                Console.WriteLine("  I = R+1.");
                return;
            }

            value = new double[1];
            value[0] = value_max;
            //
            //  Define P(K)(X) = Q(I)(X) / Q(I)(X(k)
            //
            o = qo[i];
            po[k] = qo[i];
            for (l = 0; l < o; l++)
            {
                pc[k + l * nd] = qc[i + l * r] / value[0];
                pe[k + l * nd] = qe[i + l * r];
            }

            //
            //  Modify P(1:k-1)(X).
            //
            for (j = 0; j < k; j++)
            {
                oj = po[j];
                for (l = 0; l < oj; l++)
                {
                    cj[l] = pc[j + l * nd];
                    ej[l] = pe[j + l * nd];
                }

                value = Polynomial.polynomial_value(d, oj, cj, ej, 1, xd, xIndex: + k * d);

                ok = po[k];
                for (l = 0; l < ok; l++)
                {
                    ck[l] = pc[k + l * nd];
                    ek[l] = pe[k + l * nd];
                }

                Polynomial.polynomial_axpy(-value[0], ok, ck, ek, oj, cj, ej, ref o, ref c, ref e);

                po[j] = o;
                for (l = 0; l < o; l++)
                {
                    pc[j + l * nd] = c[l];
                    pe[j + l * nd] = e[l];
                }

            }

            //
            //  Modify Q(I:downto:K+1)
            //
            for (j = i; k < j; j--)
            {
                oj = qo[j - 1];
                for (l = 0; l < oj; l++)
                {
                    cj[l] = qc[j - 1 + l * r];
                    ej[l] = qe[j - 1 + l * r];
                }

                value = Polynomial.polynomial_value(d, oj, cj, ej, 1, xd, xIndex: + k * d);

                ok = po[k];
                for (l = 0; l < ok; l++)
                {
                    ck[l] = pc[k + l * nd];
                    ek[l] = pe[k + l * nd];
                }

                Polynomial.polynomial_axpy(-value[0], ok, ck, ek, oj, cj, ej, ref o, ref c, ref e);

                qo[j] = o;
                for (l = 0; l < o; l++)
                {
                    qc[j + l * r] = c[l];
                    qe[j + l * r] = e[l];
                }
            }

            //
            //  Modify Q(I+1:R)
            //
            for (j = i + 1; j < r; j++)
            {
                oj = qo[j];
                for (l = 0; l < oj; l++)
                {
                    cj[l] = qc[j + l * r];
                    ej[l] = qe[j + l * r];
                }

                value = Polynomial.polynomial_value(d, oj, cj, ej, 1, xd, xIndex: + k * d);

                ok = po[k];
                for (l = 0; l < ok; l++)
                {
                    ck[l] = pc[k + l * nd];
                    ek[l] = pe[k + l * nd];
                }

                Polynomial.polynomial_axpy(-value[0], ok, ck, ek, oj, cj, ej, ref o, ref c, ref e);

                qo[j] = o;
                for (l = 0; l < o; l++)
                {
                    qc[j + l * r] = c[l];
                    qe[j + l * r] = e[l];
                }
            }
        }

        //
        //  Get rid of tiny coefficients.
        //
        for (i = 0; i < nd; i++)
        {
            oj = po[i];
            for (l = 0; l < oj; l++)
            {
                cj[l] = pc[i + l * nd];
                ej[l] = pe[i + l * nd];
            }

            Polynomial.polynomial_compress(oj, cj, ej, ref ok, ref ck, ref ek);

            po[i] = ok;
            for (l = 0; l < ok; l++)
            {
                pc[i + l * nd] = ck[l];
                pe[i + l * nd] = ek[l];
            }
        }
    }

    public static void lagrange_partial3(int d, int n, int nd, double[] xd, int option,
            ref int[] po, ref double[] pc, ref int[] pe, ref int n2 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGRANGE_PARTIAL3: Partial Lagrange polynomial basis from data.
        //
        //  Discussion:
        //
        //    This function, together with lagrange_partial4(), is a representation
        //    of algorithm 4.1 in the reference, modified:
        //    * for the case where the number of data points is less
        //      than the dimension of the desired polynomial space,
        //    * so that a form of "pivoting" is used
        //      to select the next polynomial as the one with maximum absolute
        //      value at the current node;
        //    * so that if the problem is not well posed, successively higher
        //      values of N are tried.
        //
        //    This function is given XD, a set of ND distinct data points in a 
        //    D dimensional space, and returns information defining a set of 
        //    ND Lagrange polynomials L(i)(X) with the property that:
        //
        //      L(i)(XD(j)) = delta(i,j)
        //
        //    This function is used in cases where ND, the number of data points, 
        //    is less than or equal to R, the dimension of the space of polynomials 
        //    in D dimensions and total degree N or less, that is:
        //
        //      ND <= R = Choose ( N + D, N )
        //
        //    There will be ND polynomials returned.  Each polynomial can have
        //    as many as R coefficients.
        //
        //    Each polynomial is given as a vector, with each entry corresponding
        //    to a nonzero coefficient.  In particular, for polynomial L(i)(X):
        //
        //      PO(i) is the order, that is, the number of nonzero coefficients;
        //      PC(i,j), for 1 <= j <= PO(i), is the coefficient of the J-th term.
        //      PE(i,j), for 1 <= j <= PO(i), encodes the exponents of the J-th term.
        //
        //    The exponent codes are a compact way of recording the exponent vector
        //    associated with each monomial.  If PE(i,j) = k, then the corresponding
        //    vector of D exponents can be determined by:
        //
        //      E = mono_unrank_grlex ( D, k );
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    24 February 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Tomas Sauer, Yuan Xu,
        //    On multivariate Lagrange interpolation,
        //    Mathematics of Computation,
        //    Volume 64, Number 211, July 1995, pages 1147-1170.
        //
        //  Parameters:
        //
        //    Input, int D, the spatial dimension.
        //
        //    Input, int N, the maximum total degree.
        //
        //    Input, int ND, the number of data points.
        //    It must be the case that ND <= R = the number of monomials 
        //    of degree N in D dimensions.
        //
        //    Input, double XD[D*ND], the data points, which must be distinct.
        //
        //    Input, int OPTION, determines the initial basis:
        //    0, monomials, 1, x, y, x^2, xy, y^2, x^3, ...
        //    1, Legendre products, 1, y, x, (3y^2-1)/2, xy, (3^x^2-1), (5y^3-3)/2, ...
        //
        //    Output, int PO[ND], the order (number of nonzero coefficients) for the 
        //    Lagrange basis polynomials.
        //
        //    Output, double **PC, the ND by R array of coefficients for the 
        //    Lagrange basis polynomials.
        //
        //    Output, int **PE, the ND by R array of exponent indices for the 
        //    Lagrange basis polynomials.
        //
        //    Output, int &N2, the adjusted value of N, which may have been
        //    increased because the interpolation problem for N was not well posed.
        //
    {
        double d_max = 0;
        double d_min = 0;
        double d_tol;
        int r;
        bool success;
        double tol;
        //
        //  Verify that the points are sufficiently distinct.
        //
        typeMethods.r8col_separation(d, nd, xd, ref d_min, ref d_max);
        d_tol = Math.Sqrt(typeMethods.r8_epsilon());

        if (d_min < d_tol)
        {
            Console.WriteLine("");
            Console.WriteLine("LAGRANGE_PARTIAL3 - Fatal error!");
            Console.WriteLine("  Some points are too close!");
            Console.WriteLine("  Minimum data point separation is = " + d_min + "");
            return;
        }

        //
        //  Search for the appropriate interpolation space.
        //
        n2 = n;
        tol = 0.0001;

        for (;;)
        {
            r = Monomial.mono_upto_enum(d, n2);

            pc = new double[nd * r];
            pe = new int[nd * r];

            success = lagrange_partial4(d, n2, r, nd, xd, option, tol, ref po, ref pc, ref pe);

            switch (success)
            {
                case true:
                    return;
                default:
                    n2 += 1;
                    Console.WriteLine("LAGRANGE_PARTIAL3 - Increase N to " + n2 + "");
                    break;
            }
        }
    }

    public static bool lagrange_partial4(int d, int n, int r, int nd, double[] xd, int option,
            double tol, ref int[] po, ref double[] pc, ref int[] pe )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LAGRANGE_PARTIAL4: Partial Lagrange polynomial basis from data.
        //
        //  Discussion:
        //
        //    This function, together with lagrange_partial3(), is a representation
        //    of algorithm 4.1 in the reference, modified:
        //    * for the case where the number of data points is less
        //      than the dimension of the desired polynomial space,
        //    * so that a form of "pivoting" is used
        //      to select the next polynomial as the one with maximum absolute
        //      value at the current node;
        //    * so that if the problem is not well posed, successively higher
        //      values of N are tried.
        //
        //    This function is given XD, a set of ND data points in a D dimensional
        //    space, and returns information defining a set of ND Lagrange polynomials
        //    L(i)(X) with the property that:
        //
        //      L(i)(XD(j)) = delta(i,j)
        //
        //    This function is used in cases where ND, the number of data points, 
        //    is less than or equal to R, the dimension of the space of polynomials 
        //    in D dimensions and total degree N or less, that is:
        //
        //      ND <= R = Choose ( N + D, N )
        //
        //    There will be ND polynomials returned.  Each polynomial can have
        //    as many as R coefficients.
        //
        //    Each polynomial is given as a vector, with each entry corresponding
        //    to a nonzero coefficient.  In particular, for polynomial L(i)(X):
        //
        //      PO(i) is the order, that is, the number of nonzero coefficients;
        //      PC(i,j), for 1 <= j <= PO(i), is the coefficient of the J-th term.
        //      PE(i,j), for 1 <= j <= PO(i), encodes the exponents of the J-th term.
        //
        //    The exponent codes are a compact way of recording the exponent vector
        //    associated with each monomial.  If PE(i,j) = k, then the corresponding
        //    vector of D exponents can be determined by:
        //
        //      E = mono_unrank_grlex ( D, k );
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 October 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Tomas Sauer, Yuan Xu,
        //    On multivariate Lagrange interpolation,
        //    Mathematics of Computation,
        //    Volume 64, Number 211, July 1995, pages 1147-1170.
        //
        //  Parameters:
        //
        //    Input, int D, the spatial dimension.
        //
        //    Input, int N, the maximum total degree.
        //
        //    Input, int R, the number of monomials in D dimensions 
        //    of total degree N or less.
        //
        //    Input, int ND, the number of data points.
        //    It must be the case that ND <= R.
        //
        //    Input, double XD[D*ND], the data points.
        //
        //    Input, int OPTION, determines the initial basis:
        //    0, monomials, 1, x, y, x^2, xy, y^2, x^3, ...
        //    1, Legendre products, 1, y, x, (3y^2-1)/2, xy, (3^x^2-1), (5y^3-3)/2, ...
        //
        //    Input, double TOL, a tolerance for the pivoting operation.
        //    If no unused polynomial can be found with a value at least TOL
        //    at the current point, the algorithm fails.
        //
        //    Output, int PO[ND], the order (number of nonzero coefficients) for the 
        //    Lagrange basis polynomials.
        //
        //    Output, double PC[ND*R], the coefficients for the 
        //    Lagrange basis polynomials.
        //
        //    Output, int PE[ND*R], the exponent indices for the 
        //    Lagrange basis polynomials.
        //
        //    Output, bool LAGRANGE_PARTIAL4, is 0 if the algorithm failed
        //    (in which case the other outputs are not useful),
        //    and 1 if it was successful.
        //
    {
        double[] c;
        double[] cj;
        double[] ck;
        int[] e;
        int[] ej;
        int[] ek;
        int i;
        int j;
        int k;
        int l;
        int[] lpp;
        int o = 0;
        int oj;
        int ok = 0;
        double[] qc;
        int[] qe;
        int[] qo;
        bool success;
        double[] value;
        double value_max;

        success = true;
        //
        //  Verify that R is acceptable.
        //
        if (r < nd)
        {
            Console.WriteLine("");
            Console.WriteLine("LAGRANGE_PARTIAL4 - Fatal error!");
            Console.WriteLine("  The value R = " + r + "");
            Console.WriteLine("  is less than ND = " + nd + "");
            return false;
        }

        //
        //  Make some work space.
        //
        c = new double[r];
        cj = new double[r];
        ck = new double[r];
        e = new int[r];
        ej = new int[r];
        ek = new int[r];
        //
        //  Initialize the polynomials Q spanning the space of N-th degree polynomials.
        //
        qo = new int[r];
        qc = new double[r * r];
        qe = new int[r * r];

        for (j = 0; j < r; j++)
        {
            qo[j] = 0;
            for (i = 0; i < r; i++)
            {
                qc[i + j * r] = 0.0;
                qe[i + j * r] = 0;
            }
        }

        //
        //  Option 0: First R D-dimensional monomials
        //  Option 1: First R D-dimensional Legendre product polynomials.
        //
        for (k = 0; k < r; k++)
        {
            switch (option)
            {
                case 0:
                    o = 1;
                    c[0] = 1.0;
                    e[0] = k + 1;
                    break;
                case 1:
                    lpp = Comp.comp_unrank_grlex(d, k + 1);
                    PolynomialNS.Legendre.lpp_to_polynomial(d, lpp, r, ref o, ref c, ref e);
                    break;
            }

            qo[k] = o;
            for (j = 0; j < o; j++)
            {
                qc[k + j * r] = c[j];
                qe[k + j * r] = e[j];
            }
        }

        //
        //  Now set up the P polynomials.
        //
        for (k = 0; k < nd; k++)
        {
            po[k] = 0;
            for (j = 0; j < r; j++)
            {
                pc[k + j * nd] = 0.0;
                pe[k + j * nd] = 0;
            }
        }

        for (k = 0; k < nd; k++)
        {
            //
            //  Find the polynomial Q(K:R)(X) which is most nonzero at X(K).
            //
            i = r + 1;
            value_max = 0.0;

            for (j = k; j < r; j++)
            {
                o = qo[j];
                for (l = 0; l < o; l++)
                {
                    c[l] = qc[j + l * r];
                    e[l] = qe[j + l * r];
                }

                value = Polynomial.polynomial_value(d, o, c, e, 1, xd, xIndex: + k * d);

                if (Math.Abs(value_max) <= Math.Abs(value[0]))
                {
                    i = j;
                    value_max = value[0];
                }
            }

            //
            //  If the best nonzero value was too small or zero, fail.
            //
            if (Math.Abs(value_max) < tol || i == r + 1)
            {
                success = false;
                Console.WriteLine("LAGRANGE_PARTIAL4 - Unacceptable VALUE_MAX = " + value_max + "");
                return success;
            }

            value = new double[1];
            value[0] = value_max;
            //
            //  Define P(K)(X) = Q(I)(X) / Q(I)(X(k)
            //
            o = qo[i];
            po[k] = qo[i];
            for (l = 0; l < o; l++)
            {
                pc[k + l * nd] = qc[i + l * r] / value[0];
                pe[k + l * nd] = qe[i + l * r];
            }

            //
            //  Modify P(1:k-1)(X).
            //
            for (j = 0; j < k; j++)
            {
                oj = po[j];
                for (l = 0; l < oj; l++)
                {
                    cj[l] = pc[j + l * nd];
                    ej[l] = pe[j + l * nd];
                }

                value = Polynomial.polynomial_value(d, oj, cj, ej, 1, xd, xIndex: + k * d);

                ok = po[k];
                for (l = 0; l < ok; l++)
                {
                    ck[l] = pc[k + l * nd];
                    ek[l] = pe[k + l * nd];
                }

                Polynomial.polynomial_axpy(-value[0], ok, ck, ek, oj, cj, ej, ref o, ref c, ref e);

                po[j] = o;
                for (l = 0; l < o; l++)
                {
                    pc[j + l * nd] = c[l];
                    pe[j + l * nd] = e[l];
                }
            }

            //
            //  Modify Q(I:downto:K+1)
            //
            for (j = i; k < j; j--)
            {
                oj = qo[j - 1];
                for (l = 0; l < oj; l++)
                {
                    cj[l] = qc[j - 1 + l * r];
                    ej[l] = qe[j - 1 + l * r];
                }

                value = Polynomial.polynomial_value(d, oj, cj, ej, 1, xd, xIndex: + k * d);

                ok = po[k];
                for (l = 0; l < ok; l++)
                {
                    ck[l] = pc[k + l * nd];
                    ek[l] = pe[k + l * nd];
                }

                Polynomial.polynomial_axpy(-value[0], ok, ck, ek, oj, cj, ej, ref o, ref c, ref e);

                qo[j] = o;
                for (l = 0; l < o; l++)
                {
                    qc[j + l * r] = c[l];
                    qe[j + l * r] = e[l];
                }
            }

            //
            //  Modify Q(I+1:R)
            //
            for (j = i + 1; j < r; j++)
            {
                oj = qo[j];
                for (l = 0; l < oj; l++)
                {
                    cj[l] = qc[j + l * r];
                    ej[l] = qe[j + l * r];
                }

                value = Polynomial.polynomial_value(d, oj, cj, ej, 1, xd, xIndex: + k * d);

                ok = po[k];
                for (l = 0; l < ok; l++)
                {
                    ck[l] = pc[k + l * nd];
                    ek[l] = pe[k + l * nd];
                }

                Polynomial.polynomial_axpy(-value[0], ok, ck, ek, oj, cj, ej, ref o, ref c, ref e);
                    
                qo[j] = o;
                for (l = 0; l < o; l++)
                {
                    qc[j + l * r] = c[l];
                    qe[j + l * r] = e[l];
                }
            }
        }

        //
        //  Get rid of tiny coefficients.
        //
        for (i = 0; i < nd; i++)
        {
            oj = po[i];
            for (l = 0; l < oj; l++)
            {
                cj[l] = pc[i + l * nd];
                ej[l] = pe[i + l * nd];
            }

            Polynomial.polynomial_compress(oj, cj, ej, ref ok, ref ck, ref ek);

            po[i] = ok;
            for (l = 0; l < ok; l++)
            {
                pc[i + l * nd] = ck[l];
                pe[i + l * nd] = ek[l];
            }
        }
        return success;
    }
}