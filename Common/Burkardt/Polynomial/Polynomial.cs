﻿using System;
using Burkardt.Types;

namespace Burkardt.PolynomialNS
{
    public static class Polynomial
    {
        public static void polynomial_axpy(double s, int o1, double[] c1, int[] e1, int o2,
                double[] c2, int[] e2, ref int o, ref double[] c, ref int[] e)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYNOMIAL_AXPY adds a multiple of one polynomial to another.
            //
            //  Discussion:
            //
            //    P(X) = S * P1(X) + P2(X)
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    21 January 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, double S, the multiplier of polynomial 1.
            //
            //    Input, int O1, the "order" of polynomial 1.
            //
            //    Input, double C1[O1], the coefficients of polynomial 1.
            //
            //    Input, int E1[O1], the indices of the exponents of 
            //    polynomial 1.
            //
            //    Input, int O2, the "order" of polynomial 2.
            //
            //    Input, double C2[O2], the coefficients of polynomial 2.
            //
            //    Input, int E2[O2], the indices of the exponents of 
            //    polynomial 2.
            //
            //    Output, int &O, the "order" of the polynomial sum.
            //
            //    Output, double C[O], the coefficients of the polynomial sum.
            //
            //    Output, int E[O], the indices of the exponents of 
            //    the polynomial sum.
            //
        {
            double[] c3;
            int[] e3;
            int i;
            int o3;
            double[] sc1;

            o3 = o1 + o2;

            c3 = new double[o3];
            e3 = new int[o3];
            sc1 = new double[o1];

            for (i = 0; i < o1; i++)
            {
                sc1[i] = s * c1[i];
            }

            typeMethods.r8vec_concatenate(o1, sc1, o2, c2, ref c3);
            typeMethods.i4vec_concatenate(o1, e1, o2, e2, ref e3);

            polynomial_sort(o3, ref c3, ref e3);
            polynomial_compress(o3, c3, e3, ref o, ref c, ref e);
        }

        public static void polynomial_compress(int o1, double[] c1, int[] e1, ref int o2, ref double[] c2,
                ref int[] e2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYNOMIAL_COMPRESS compresses a polynomial.
            //
            //  Discussion:
            //
            //    The function polynomial_sort ( ) should be called first.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    21 January 2014
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int O1, the "order" of the polynomial.
            //
            //    Input, double C1[O1], the coefficients of the polynomial.
            //
            //    Input, int E1[O1], the indices of the exponents of 
            //    the polynomial.
            //
            //    Output, int &O2, the "order" of the polynomial.
            //
            //    Output, double C2[O2], the coefficients of the polynomial.
            //
            //    Output, int E2[O2], the indices of the exponents of 
            //    the polynomial.
            //
        {
            int get;
            int put;
            const double r8_epsilon_sqrt = 0.1490116119384766E-07;

            get = 0;
            put = 0;

            while (get < o1)
            {
                get = get + 1;

                if (Math.Abs(c1[get - 1]) <= r8_epsilon_sqrt)
                {
                    continue;
                }

                if (0 == put)
                {
                    put = put + 1;
                    c2[put - 1] = c1[get - 1];
                    e2[put - 1] = e1[get - 1];
                }
                else
                {
                    if (e2[put - 1] == e1[get - 1])
                    {
                        c2[put - 1] = c2[put - 1] + c1[get - 1];
                    }
                    else
                    {
                        put = put + 1;
                        c2[put - 1] = c1[get - 1];
                        e2[put - 1] = e1[get - 1];
                    }
                }
            }

            o2 = put;

            return;
        }

        public static void polynomial_print(int d, int o, double[] c, int[] e, string title)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYNOMIAL_PRINT prints a polynomial.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    01 December 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int D, the spatial dimension.
            //
            //    Input, int O, the "order" of the polynomial, that is,
            //    simply the number of terms.
            //
            //    Input, double C[O], the coefficients.
            //
            //    Input, int E[O], the indices of the exponents.
            //
            //    Input, string TITLE, a title.
            //
        {
            int[] f;
            int i;
            int j;

            Console.WriteLine(title);

            if (o == 0)
            {
                Console.WriteLine("      0.");
            }
            else
            {
                for (j = 0; j < o; j++)
                {
                    string cout = "    ";
                    if (c[j] < 0.0)
                    {
                        cout += "- ";
                    }
                    else
                    {
                        cout += "+ ";
                    }

                    cout += Math.Abs(c[j]) + " * x^(";

                    f = Monomial.mono_unrank_grlex(d, e[j]);
                    for (i = 0; i < d; i++)
                    {
                        cout += f[i];
                        if (i < d - 1)
                        {
                            cout += ",";
                        }
                        else
                        {
                            cout += ")";
                        }
                    }

                    if (j == o - 1)
                    {
                        cout += ".";
                    }

                    Console.WriteLine(cout);
                }
            }
        }

        public static void polynomial_sort(int o, ref double[] c, ref int[] e)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYNOMIAL_SORT sorts the information in a polynomial.
            //
            //  Discussion
            //
            //    The coefficients C and exponents E are rearranged so that 
            //    the elements of E are in ascending order.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    01 December 2013
            //
            //  Author:
            //
            //    John Burkardt
            //
            //  Parameters:
            //
            //    Input, int O, the "order" of the polynomial.
            //
            //    Input/output, double C[O], the coefficients of the polynomial.
            //
            //    Input/output, int E[O], the indices of the exponents of 
            //    the polynomial.
            //
        {
            int[] indx = typeMethods.i4vec_sort_heap_index_a(o, e);

            typeMethods.i4vec_permute(o, indx, ref e);
            typeMethods.r8vec_permute(o, indx, c);
        }

        public static double[] polynomial_value(int d, int o, double[] c, int[] e, int nx,
                double[] x, int xIndex = 0)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYNOMIAL_VALUE evaluates a polynomial.
            //
            //  Discussion:
            //
            //    The polynomial is evaluated term by term, and no attempt is made to
            //    use an approach such as Horner's method to speed up the process.
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
            //  Parameters:
            //
            //    Input, int D, the spatial dimension.
            //
            //    Input, int O, the "order" of the polynomial.
            //
            //    Input, double C[O], the coefficients of the polynomial.
            //
            //    Input, int E(O), the indices of the exponents 
            //    of the polynomial.
            //
            //    Input, int NX, the number of evaluation points.
            //
            //    Input, double X[D*NX], the coordinates of the evaluation points.
            //
            //    Output, double POLYNOMIAL_VALUE[NX], the value of the polynomial at X.
            //
        {
            int[] f;
            int j;
            int k;
            double[] p;
            double[] v;

            p = new double[nx];

            for (k = 0; k < nx; k++)
            {
                p[k] = 0.0;
            }

            for (j = 0; j < o; j++)
            {
                f = Monomial.mono_unrank_grlex(d, e[j]);
                v = Monomial.mono_value(d, nx, f, x, xIndex);
                for (k = 0; k < nx; k++)
                {
                    p[k] = p[k] + c[j] * v[k];
                }
            }

            return p;
        }

        public static void poly(string code, ref int[] rexp, ref int[] sexp )

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    POLY returns the polynomial terms associated with any available element.
        //
        //  Formula:
        //
        //    Given coefficients A(I), the polynomial interpolant at (R,S) is
        //
        //      P(R,S) = sum ( 1 <= I <= N ) A(I) * R**REXP(I) * S**SEXP(I) 
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    02 April 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, string CODE, identifies the element desired.
        //    Legal values include 'Q4', 'Q8', 'Q9', 'Q12', 'Q16', 'QL', 'T3', 
        //    'T4', 'T6' and 'T10'.
        //
        //    Output, int REXP(N), SEXP(N), the powers of R and S associated
        //    with each polynomial.
        //
        {
            if (code == "Q4")
            {
                poly_q4(ref rexp, ref sexp);
            }
            else if (code == "Q8")
            {
                poly_q8(ref rexp, ref sexp);
            }
            else if (code == "Q9")
            {
                poly_q9(ref rexp, ref sexp);
            }
            else if (code == "Q12")
            {
                poly_q12(ref rexp, ref sexp);
            }
            else if (code == "Q16")
            {
                poly_q16(ref rexp, ref sexp);
            }
            else if (code == "QL")
            {
                poly_ql(ref rexp, ref sexp);
            }
            else if (code == "T3")
            {
                poly_t3(ref rexp, ref sexp);
            }
            else if (code == "T4")
            {
                Console.WriteLine("");
                Console.WriteLine("POLY - Fatal error!");
                Console.WriteLine("  The T4 element does not follow the pattern!");
            }
            else if (code == "T6")
            {
                poly_t6(ref rexp, ref sexp);
            }
            else if (code == "T10")
            {
                poly_t10(ref rexp, ref sexp);
            }
            else
            {
                Console.WriteLine("");
                Console.WriteLine("POLY - Fatal error!");
                Console.WriteLine("  Illegal value of CODE = " + code + "");
            }
        }

        public static void poly_q4(ref int[] rexp, ref int[] sexp )

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    POLY_Q4 returns the monomials associated with a 4 node quadrilateral.
        //
        //  Element Q4:
        //
        //    |
        //    1  4-----3
        //    |  |     |
        //    |  |     |
        //    S  |     |
        //    |  |     |
        //    |  |     |
        //    0  1-----2
        //    |
        //    +--0--R--1-->
        //
        //  Formula:
        //
        //    Given coefficients A(I), the polynomial interpolant at (R,S) is
        //
        //      P(R,S) = sum ( 1 <= I <= N ) A(I) * R**REXP(I) * S**SEXP(I) 
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 March 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, int REXP[4], SEXP[4], the powers of R and S associated
        //    with each polynomial.
        //
        {
            rexp[0] = 0;
            rexp[1] = 0;
            rexp[2] = 1;
            rexp[3] = 1;

            sexp[0] = 0;
            sexp[1] = 1;
            sexp[2] = 0;
            sexp[3] = 1;
        }

        public static void poly_q8(ref int[] rexp, ref int[] sexp )

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    POLY_Q8 returns the monomials associated with an 8 node quadrilateral.
        //
        //  Element Q8:
        //
        //    |
        //    1  4--7--3
        //    |  |     |
        //    |  |     |
        //    S  8     6
        //    |  |     |
        //    |  |     |
        //    0  1--5--2
        //    |
        //    +--0--R--1-->
        //
        //  Formula:
        //
        //    Given coefficients A(I), the polynomial interpolant at (R,S) is
        //
        //      P(R,S) = sum ( 1 <= I <= N ) A(I) * R**REXP(I) * S**SEXP(I) 
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 March 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, int REXP[8], SEXP[8], the powers of R and S associated
        //    with each monomial.
        //
        {
            rexp[0] = 0;
            rexp[1] = 0;
            rexp[2] = 1;
            rexp[3] = 0;
            rexp[4] = 1;
            rexp[5] = 2;
            rexp[6] = 1;
            rexp[7] = 2;

            sexp[0] = 0;
            sexp[1] = 1;
            sexp[2] = 0;
            sexp[3] = 2;
            sexp[4] = 1;
            sexp[5] = 0;
            sexp[6] = 2;
            sexp[7] = 1;
        }

        public static void poly_q9(ref int[] rexp, ref int[] sexp )

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    POLY_Q9 returns the monomials associated with a 9 node quadrilateral.
        //
        //  Element Q9:
        //
        //    |
        //    1  4--7--3
        //    |  |     |
        //    |  |     |
        //    S  8  9  6
        //    |  |     |
        //    |  |     |
        //    0  1--5--2
        //    |
        //    +--0--R--1-->
        //
        //  Formula:
        //
        //    Given coefficients A(I), the polynomial interpolant at (R,S) is
        //
        //      P(R,S) = sum ( 1 <= I <= N ) A(I) * R**REXP(I) * S**SEXP(I) 
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 March 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, int REXP[9], SEXP[9], the powers of R and S associated
        //    with each monomial.
        //
        {
            rexp[0] = 0;
            rexp[1] = 0;
            rexp[2] = 1;
            rexp[3] = 0;
            rexp[4] = 1;
            rexp[5] = 2;
            rexp[6] = 1;
            rexp[7] = 2;
            rexp[8] = 2;

            sexp[0] = 0;
            sexp[1] = 1;
            sexp[2] = 0;
            sexp[3] = 2;
            sexp[4] = 1;
            sexp[5] = 0;
            sexp[6] = 2;
            sexp[7] = 1;
            sexp[8] = 2;
        }

        public static void poly_q12(ref int[] rexp, ref int[] sexp )

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    POLY_Q12 returns the monomials associated with a 12 node quadrilateral.
        //
        //  Element Q12:
        //
        //    |
        //    1  9-10-11-12
        //    |  |        |
        //    |  7        8
        //    S  |        |
        //    |  5        6
        //    |  |        |
        //    0  1--2--3--4
        //    |
        //    +--0---R---1-->
        //
        //  Formula:
        //
        //    Given coefficients A(I), the polynomial interpolant at (R,S) is
        //
        //      P(R,S) = sum ( 1 <= I <= N ) A(I) * R**REXP(I) * S**SEXP(I) 
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 March 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, int REXP[12], SEXP[12], the powers of R and S associated
        //    with each monomial.
        //
        {
            rexp[0] = 0;
            rexp[1] = 0;
            rexp[2] = 1;
            rexp[3] = 0;
            rexp[4] = 1;
            rexp[5] = 2;
            rexp[6] = 0;
            rexp[7] = 1;
            rexp[8] = 2;
            rexp[9] = 3;
            rexp[10] = 1;
            rexp[11] = 3;

            sexp[0] = 0;
            sexp[1] = 1;
            sexp[2] = 0;
            sexp[3] = 2;
            sexp[4] = 1;
            sexp[5] = 0;
            sexp[6] = 3;
            sexp[7] = 2;
            sexp[8] = 1;
            sexp[9] = 0;
            sexp[10] = 3;
            sexp[11] = 1;
        }

        public static void poly_q16(ref int[] rexp, ref int[] sexp )

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    POLY_Q16 returns the monomials associated with a 16 node quadrilateral.
        //
        //  Element Q16:
        //
        //    |
        //    1 13--14--15--16
        //    |  |   :   :   |
        //    |  |   :   :   |
        //    |  9..10..11..12
        //    S  |   :   :   |
        //    |  |   :   :   |
        //    |  5...6...7...8
        //    |  |   :   :   |
        //    |  |   :   :   |  
        //    0  1---2---3---4
        //    |
        //    +--0-----R-----1-->
        //
        //  Formula:
        //
        //    Given coefficients A(I), the polynomial interpolant at (R,S) is
        //
        //      P(R,S) = sum ( 1 <= I <= N ) A(I) * R**REXP(I) * S**SEXP(I) 
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 March 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, int REXP[16], SEXP[16], the powers of R and S associated
        //    with each monomial.
        //
        {
            rexp[0] = 0;
            rexp[1] = 0;
            rexp[2] = 1;
            rexp[3] = 0;
            rexp[4] = 1;
            rexp[5] = 2;
            rexp[6] = 0;
            rexp[7] = 1;
            rexp[8] = 2;
            rexp[9] = 3;
            rexp[10] = 1;
            rexp[11] = 2;
            rexp[12] = 3;
            rexp[13] = 2;
            rexp[14] = 3;
            rexp[15] = 3;

            sexp[0] = 0;
            sexp[1] = 1;
            sexp[2] = 0;
            sexp[3] = 2;
            sexp[4] = 1;
            sexp[5] = 0;
            sexp[6] = 3;
            sexp[7] = 2;
            sexp[8] = 1;
            sexp[9] = 0;
            sexp[10] = 3;
            sexp[11] = 2;
            sexp[12] = 1;
            sexp[13] = 3;
            sexp[14] = 2;
            sexp[15] = 3;
        }

        public static void poly_ql(ref int[] rexp, ref int[] sexp )

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    POLY_QL returns the monomials for a quadratic/linear quadrilateral.
        //
        //  Element QL:
        //
        //    |
        //    1  4---5---6
        //    |  |       |
        //    |  |       |
        //    S  |       |
        //    |  |       |
        //    |  |       |
        //    0  1---2---3
        //    |
        //    +--0---R---1-->
        //
        //  Formula:
        //
        //    Given coefficients A(I), the polynomial interpolant at (R,S) is
        //
        //      P(R,S) = sum ( 1 <= I <= N ) A(I) * R**REXP(I) * S**SEXP(I) 
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 March 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, int REXP[6], SEXP[6], the powers of R and S associated
        //    with each monomial.
        //
        {
            rexp[0] = 0;
            rexp[1] = 0;
            rexp[2] = 1;
            rexp[3] = 1;
            rexp[4] = 2;
            rexp[5] = 2;

            sexp[0] = 0;
            sexp[1] = 1;
            sexp[2] = 0;
            sexp[3] = 1;
            sexp[4] = 0;
            sexp[5] = 1;
        }

        public static void poly_t3(ref int[] rexp, ref int[] sexp )

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    POLY_T3 returns the monomials associated with a 3 node triangle.
        //
        //  Element T3:
        //
        //    |
        //    1  3
        //    |  ..
        //    |  . .
        //    S  .  .
        //    |  .   .
        //    |  .    .
        //    0  1-----2
        //    |
        //    +--0--R--1-->
        //
        //  Formula:
        //
        //    Given coefficients A(I), the polynomial interpolant at (R,S) is
        //
        //      P(R,S) = sum ( 1 <= I <= N ) A(I) * R**REXP(I) * S**SEXP(I) 
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 March 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, int REXP[3], SEXP[3], the powers of R and S associated
        //    with each monomial.
        //
        {
            rexp[0] = 0;
            rexp[1] = 0;
            rexp[2] = 1;

            sexp[0] = 0;
            sexp[1] = 1;
            sexp[2] = 0;
        }

        public static void poly_t6(ref int[] rexp, ref int[] sexp )

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    POLY_T6 returns the monomials associated with a 6 node triangle.
        //
        //  Element T6:
        //
        //    |
        //    1  3
        //    |  ..
        //    |  . .
        //    S  6  5
        //    |  .   .
        //    |  .    .
        //    0  1--4--2
        //    |
        //    +--0--R--1-->
        //
        //  Formula:
        //
        //    Given coefficients A(I), the polynomial interpolant at (R,S) is
        //
        //      P(R,S) = sum ( 1 <= I <= N ) A(I) * R**REXP(I) * S**SEXP(I) 
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 March 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, int REXP[6], SEXP[6], the powers of R and S associated
        //    with each monomial.
        //
        {
            rexp[0] = 0;
            rexp[1] = 0;
            rexp[2] = 1;
            rexp[3] = 0;
            rexp[4] = 1;
            rexp[5] = 2;

            sexp[0] = 0;
            sexp[1] = 1;
            sexp[2] = 0;
            sexp[3] = 2;
            sexp[4] = 1;
            sexp[5] = 0;
        }

        public static void poly_t10(ref int[] rexp, ref int[] sexp )

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    POLY_T10 returns the monomials associated with a 10 node triangle.
        //
        //  Element T10:
        //
        //    |
        //    1  10
        //    |  ..
        //    |  . .
        //    |  8  9
        //    |  .   .
        //    S  .    .
        //    |  5  6  7
        //    |  .      .
        //    |  .       .
        //    0  1--2--3--4
        //    |
        //    +--0----R---1-->
        //
        //  Formula:
        //
        //    Given coefficients A(I), the polynomial interpolant at (R,S) is
        //
        //      P(R,S) = sum ( 1 <= I <= N ) A(I) * R**REXP(I) * S**SEXP(I) 
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    31 March 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Output, int REXP[10], SEXP[10], the powers of R and S associated
        //    with each monomial.
        //
        {
            rexp[0] = 0;
            rexp[1] = 0;
            rexp[2] = 1;
            rexp[3] = 0;
            rexp[4] = 1;
            rexp[5] = 2;
            rexp[6] = 0;
            rexp[7] = 1;
            rexp[8] = 2;
            rexp[9] = 3;

            sexp[0] = 0;
            sexp[1] = 1;
            sexp[2] = 0;
            sexp[3] = 2;
            sexp[4] = 1;
            sexp[5] = 0;
            sexp[6] = 3;
            sexp[7] = 2;
            sexp[8] = 1;
            sexp[9] = 0;
        }

    }
}