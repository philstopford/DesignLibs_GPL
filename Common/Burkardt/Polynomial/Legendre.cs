using System;
using Burkardt.PolynomialNS;

namespace Burkardt.PolynomialNS
{
    public static class Legendre
    {
        public static void lp_coefficients(int n, ref int o, ref double[] c, ref int[] f )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LP_COEFFICIENTS: coefficients of Legendre polynomials P(n,x).
        //
        //  First terms:
        //
        //     1
        //     0     1
        //    -1/2   0      3/2
        //     0    -3/2    0     5/2
        //     3/8   0    -30/8   0     35/8
        //     0    15/8    0   -70/8    0     63/8
        //    -5/16  0    105/16  0   -315/16   0    231/16
        //     0   -35/16   0   315/16   0   -693/16   0    429/16
        //
        //     1.00000
        //     0.00000  1.00000
        //    -0.50000  0.00000  1.50000
        //     0.00000 -1.50000  0.00000  2.5000
        //     0.37500  0.00000 -3.75000  0.00000  4.37500
        //     0.00000  1.87500  0.00000 -8.75000  0.00000  7.87500
        //    -0.31250  0.00000  6.56250  0.00000 -19.6875  0.00000  14.4375
        //     0.00000 -2.1875   0.00000  19.6875  0.00000 -43.3215  0.00000  26.8125
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 September 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Reference:
        //
        //    Milton Abramowitz, Irene Stegun,
        //    Handbook of Mathematical Functions,
        //    National Bureau of Standards, 1964,
        //    ISBN: 0-486-61272-4,
        //    LC: QA47.A34.
        //
        //    Daniel Zwillinger, editor,
        //    CRC Standard Mathematical Tables and Formulae,
        //    30th Edition,
        //    CRC Press, 1996.
        //
        //  Parameters:
        //
        //    Input, int N, the highest order polynomial to evaluate.
        //    Note that polynomials 0 through N will be evaluated.
        //
        //    Output, int &O, the number of coefficients.
        //
        //    Output, double C[(N+2)/2], the coefficients of the Legendre
        //    polynomial of degree N.
        //
        //    Output, int F[(N+2)/2], the exponents.
        //
        {
            double[] ctable;
            int i;
            int j;
            int k;

            ctable = new double[(n + 1) * (n + 1)];

            for (i = 0; i <= n; i++)
            {
                for (j = 0; j <= n; j++)
                {
                    ctable[i + j * (n + 1)] = 0.0;
                }
            }

            ctable[0 + 0 * (n + 1)] = 1.0;

            if (0 < n)
            {
                ctable[1 + 1 * (n + 1)] = 1.0;

                for (i = 2; i <= n; i++)
                {
                    for (j = 0; j <= i - 2; j++)
                    {
                        ctable[i + j * (n + 1)] =
                            (double) (-i + 1) * ctable[i - 2 + j * (n + 1)] / (double) i;
                    }

                    for (j = 1; j <= i; j++)
                    {
                        ctable[i + j * (n + 1)] = ctable[i + j * (n + 1)]
                                                  + (double) (i + i - 1) * ctable[i - 1 + (j - 1) * (n + 1)] /
                                                  (double) i;
                    }
                }
            }

            //
            //  Extract the nonzero data from the alternating columns of the last row.
            //
            o = (n + 2) / 2;

            k = o;
            for (j = n; 0 <= j; j = j - 2)
            {
                k = k - 1;
                c[k] = ctable[n + j * (n + 1)];
                f[k] = j;
            }
        }

        public static void lpp_to_polynomial(int m, int[] l, int o_max, ref int o, ref double[] c, ref int[] e )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LPP_TO_POLYNOMIAL writes a Legendre Product Polynomial as a polynomial.
        //
        //  Discussion:
        //
        //    For example, if 
        //      M = 3,
        //      L = ( 1, 0, 2 ),
        //    then
        //      L(1,0,2)(X,Y,Z) 
        //      = L(1)(X) * L(0)(Y) * L(2)(Z)
        //      = X * 1 * ( 3Z^2-1)/2
        //      = - 1/2 X + (3/2) X Z^2
        //    so
        //      O = 2 (2 nonzero terms)
        //      C = -0.5
        //           1.5
        //      E = 4    <-- index in 3-space of exponent (1,0,0)
        //          15   <-- index in 3-space of exponent (1,0,2)
        //
        //    The output value of O is no greater than
        //      O_MAX = product ( 1 <= I <= M ) (L(I)+2)/2
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 September 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the spatial dimension.
        //
        //    Input, int L[M], the index of each Legendre product 
        //    polynomial factor.  0 <= L(*).
        //
        //    Input, int O_MAX, an upper limit on the size of the 
        //    output arrays.
        //      O_MAX = product ( 1 <= I <= M ) (L(I)+2)/2.
        //
        //    Output, int &O, the "order" of the polynomial product.
        //
        //    Output, double C[O], the coefficients of the polynomial product.
        //
        //    Output, int E[O], the indices of the exponents of the 
        //    polynomial product.
        //
        {
            double[] c1;
            double[] c2;
            int[] e1;
            int[] e2;
            int[] f2;
            int i;
            int i1;
            int i2;
            int j1;
            int j2;
            int o1;
            int o2 = 0;
            int[] p = new int[1];
            int[] pp;

            c1 = new double[o_max];
            c2 = new double[o_max];
            e1 = new int[o_max];
            e2 = new int[o_max];
            f2 = new int[o_max];
            pp = new int[m];

            o1 = 1;
            c1[0] = 1.0;
            e1[0] = 1;
            //
            //  Implicate one factor at a time.
            //
            for (i = 0; i < m; i++)
            {
                lp_coefficients(l[i], ref o2, ref c2, ref f2);

                o = 0;

                for (j2 = 0; j2 < o2; j2++)
                {
                    for (j1 = 0; j1 < o1; j1++)
                    {
                        c[o] = c1[j1] * c2[j2];
                        if (0 < i)
                        {
                            p = Monomial.mono_unrank_grlex(i, e1[j1]);
                        }

                        for (i2 = 0; i2 < i; i2++)
                        {
                            pp[i2] = p[i2];
                        }

                        pp[i] = f2[j2];
                        e[o] = Monomial.mono_rank_grlex(i + 1, pp);
                        o = o + 1;
                        if (0 < i)
                        {
                            p = null;
                        }
                    }
                }

                Polynomial.polynomial_sort(o, c, e);
                Polynomial.polynomial_compress(o, c, e, ref o, ref c, ref e);

                o1 = o;
                for (i1 = 0; i1 < o; i1++)
                {
                    c1[i1] = c[i1];
                    e1[i1] = e[i1];
                }
            }
        }
    }
}