using System;
using Burkardt.MonomialNS;
using Burkardt.Types;

namespace Burkardt.PolynomialNS
{
    public static class Polynomial
    {
        public static double[] kjacoypols3(double x, double y, double a, double b, int n)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    KJACOYPOLS3 evaluates modified Jacobi polynomials.
            //
            //  Discussion:
            //
            //    This procedure evaluates Jacobi polynomials multiplied by
            //    specific polynomials given by the formula
            //      P_n^{(a,b)} (x/y) y^n
            //    at the user-provided point x/y.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU GPL license.
            //
            //  Modified:
            //
            //    11 July 2014
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
            //    This C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Hong Xiao, Zydrunas Gimbutas,
            //    A numerical algorithm for the construction of efficient quadrature
            //    rules in two and higher dimensions,
            //    Computers and Mathematics with Applications,
            //    Volume 59, 2010, pages 663-676.
            //
            //  Parameters:
            //
            //    Input, double X, Y, define the evaluation point X/Y.
            //
            //    Input, double A, B, the parameters.
            //
            //    Input, int N, the highest degree to compute.
            //
            //    Output, double KJACOYPOLS3[N+1], the polynomial values.
            //
        {
            double alpha;
            double beta;
            int k;
            double k_r8;
            double pk;
            double pkm1;
            double pkp1;
            double[] pols;

            pols = new double[n + 1];

            pkp1 = 1.0;
            pols[0] = pkp1;

            if (n == 0)
            {
                return pols;
            }

            pk = pkp1;
            pkp1 = 0.5 * ((a - b) * y + (2.0 + a + b) * x);
            pols[1] = pkp1;

            if (n == 1)
            {
                return pols;
            }

            for (k = 2; k <= n; k++)
            {
                k_r8 = (double)(k);

                alpha = (2.0 * k_r8 + a + b - 1.0)
                        * (a - b) * (a + b) * y
                        + (2.0 * k_r8 + a + b - 1.0)
                        * (2.0 * k_r8 + a + b - 2.0)
                        * (2.0 * k_r8 + a + b) * x;

                beta = 2.0 * (k_r8 + a - 1.0)
                           * (k_r8 + b - 1.0)
                           * (2.0 * k_r8 + a + b) * y * y;

                pkm1 = pk;
                pk = pkp1;
                pkp1 = (alpha * pk - beta * pkm1)
                       / (2.0 * k_r8 * (k_r8 + a + b)
                          * (2.0 * k_r8 + a + b - 2.0));

                pols[k] = pkp1;
            }

            return pols;
        }

        public static double[] ortho3eva(int degree, double[] xyz)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    ORTHO3EVA evaluates polynomials orthogonal in the reference triangle.
            //
            //  Discussion:
            //
            //    This procedure evaluates the Koornwinder's orthogonal polynomial
            //    up to order DEGREE on the reference tetrahedron with vertices
            //      (-1, -1/Sqrt(3), -1/Sqrt(6)),
            //      ( 0,  2/Sqrt(3), -1/Sqrt(6)),
            //      ( 1, -1/Sqrt(3), -1/Sqrt(6)),
            //      ( 0,  0,      3/Sqrt(6))
            //
            //    The polynomials are ordered by their order, and in each order,
            //    they are ordered lexicographically in terms of their indices
            //    (m,n,k).
            //
            //    The total number of polynomials should be
            //      NVALS = ( ( DEGREE + 1 ) * ( DEGREE + 2 ) * ( DEGREE + 3 ) ) / 6.
            //
            //    The calculations are based on Koornwinder's representation
            //    of the orthogonal polynomials on the right triangle
            //      (-1,-1,-1), (-1,1,-1), (1,-1,-1),(-1,-1,1)
            //    given by:
            //      K_mnk(x,y,z) =
            //        P_m ((2x+2+y+z)/(-y-z)) * ((-y-z)/2)^m *
            //        P_n^{2m+1,0}((2y+1+z)/(1-z)) * ((1-z)/2)^{n}
            //        P_k^{2m+2n+2,0} (z)
            //    with the input (x,y,z) transformed as
            //      x = -1/2 + xold -yold/s3 - zold/s6
            //      y = -1/2 + 2/s3 * yold - zold/s6
            //      z = -1/2 + s6/2 * zold
            //    where
            //      s3=sqrt(3)
            //      s6=sqrt(6)
            //    and
            //      P_m is the Legendre polynomial of degree m,
            //      P_n^{2m+1,0} the Jacobi polynomial of degree n, parameters 2*m+1 and 0.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license.
            //
            //  Modified:
            //
            //    11 July 2014
            //
            //  Author:
            //
            //    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
            //    This C++ version by John Burkardt.
            //
            //  Reference:
            //
            //    Hong Xiao, Zydrunas Gimbutas,
            //    A numerical algorithm for the construction of efficient quadrature
            //    rules in two and higher dimensions,
            //    Computers and Mathematics with Applications,
            //    Volume 59, 2010, pages 663-676.
            //
            //  Parameters:
            //
            //    Input, int DEGREE, the maximum degree.
            //
            //    Input, double XYZ[3], the evaluation point.
            //
            //    Output, double ORTHO3EVA[NVALS], the polynomial values.
            //
        {
            double[] f;
            double[] f1;
            double[] f2s;
            double[] f3s;
            double[] fvals;
            int i;
            int j;
            int k;
            int m;
            int mmax;
            int n;
            int ncount;
            int nvals;
            double p1;
            double p2;
            double scale;
            double[] uvw;
            double x1;
            double y1;
            //
            //  Convert coordinates from reference to Koornwinder tetrahedron.
            //
            uvw = typeMethods.ref_to_koorn(xyz);
            //
            //  Compute F1.
            //
            p1 = 0.5 * (2.0 * uvw[0] + 2.0 + uvw[1] + uvw[2]);
            p2 = -0.5 * (uvw[1] + uvw[2]);

            f1 = LegendreScaled.klegeypols(p1, p2, degree);
            //
            //  Compute F2S.
            //
            f2s = new double[(degree + 1) * (degree + 1)];

            for (j = 1; j <= degree + 1; j++)
            {
                for (i = 1; i <= degree + 1; i++)
                {
                    f2s[i - 1 + (j - 1) * (degree + 1)] = 0.0;
                }
            }

            for (m = 0; m <= degree; m++)
            {
                x1 = 0.5 * (2.0 * uvw[1] + 1.0 + uvw[2]);
                y1 = 0.5 * (1.0 - uvw[2]);
                p1 = (double)(2 * m + 1);
                p2 = 0.0;

                f = kjacoypols3(x1, y1, p1, p2, degree - m);

                for (i = 1; i <= degree + 1 - m; i++)
                {
                    f2s[i - 1 + m * (degree + 1)] = f[i - 1];
                }
            }

            //
            //  Compute F3S.
            //
            f3s = new double[(degree + 1) * (degree + 1)];

            for (j = 1; j <= degree + 1; j++)
            {
                for (i = 1; i <= degree + 1; i++)
                {
                    f3s[i - 1 + (j - 1) * (degree + 1)] = 0.0;
                }
            }

            x1 = uvw[2];
            y1 = 1.0;
            p2 = 0.0;

            for (j = 1; j <= degree + 1; j++)
            {
                p1 = (double)(2 * j);

                f = kjacoypols3(x1, y1, p1, p2, degree + 1 - j);

                for (i = 1; i <= degree + 2 - j; i++)
                {
                    f3s[i - 1 + (j - 1) * (degree + 1)] = f[i - 1];
                }
            }

            //
            //  Construct FVALS.
            //
            nvals = ((degree + 1) * (degree + 2) * (degree + 3)) / 6;
            fvals = new double[nvals];

            ncount = 0;

            for (mmax = 0; mmax <= degree; mmax++)
            {
                for (m = 0; m <= mmax; m++)
                {
                    for (n = 0; n <= mmax - m; n++)
                    {
                        k = mmax - m - n;

                        scale = Math.Sqrt
                        (
                            4.0
                            / (double)(2 * mmax + 3)
                            / (double)(2 * m + 1)
                            / (double)(n + m + 1)
                            / Math.Sqrt(2.0)
                        );

                        fvals[ncount] =
                            f1[m] *
                            f2s[n + m * (degree + 1)] *
                            f3s[k + (m + n) * (degree + 1)] / scale;

                        ncount = ncount + 1;
                    }
                }
            }

            return fvals;
        }

        public static double ts_mult ( double[] u, double h, int n )

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    TS_MULT evaluates a polynomial.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    17 May 2013
            //
            //  Author:
            //
            //    Original C++ version by Nick Hale.
            //    This C++ version by John Burkardt.
            //
            //  Parameters:
            //
            //    Input, double U[N+1], the polynomial coefficients.
            //    U[0] is ignored.
            //
            //    Input, double H, the polynomial argument.
            //
            //    Input, int N, the number of terms to compute.
            //
            //    Output, double TS_MULT, the value of the polynomial.
            //
        {
            double hk;
            int k;
            double ts;
  
            ts = 0.0;
            hk = 1.0;
            for ( k = 1; k<= n; k++ )
            {
                ts = ts + u[k] * hk;
                hk = hk * h;
            }
            return ts;
        }
        
        public static void polynomial_add ( int o1, double[] c1, int[] e1, int o2, double[] c2, 
        int[] e2, ref int o, ref double[] c, ref int[] e )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYNOMIAL_ADD adds two polynomials.
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
            o = o1 + o2;
            typeMethods.r8vec_concatenate ( o1, c1, o2, c2, ref c );
            typeMethods.i4vec_concatenate ( o1, e1, o2, e2, ref e );

            polynomial_sort ( o, ref c, ref e );
            polynomial_compress ( o, c, e, ref o, ref c, ref e );
        }
        
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
        }

        public static void polynomial_dif(int m, int o1, double[] c1, int[] e1, int[] dif,
                ref int o2, ref double[] c2, ref int[] e2)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYNOMIAL_DIF differentiates a polynomial.
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
            //    Input, int M, the spatial dimension.
            //
            //    Input, int O1, the "order" of polynomial 1.
            //
            //    Input, double C1[O1], the coefficients of polynomial 1.
            //
            //    Input, int E1[O1], the indices of the exponents of 
            //    polynomial 1.
            //
            //    Input, int DIF[M], indicates the number of 
            //    differentiations in each component.
            //
            //    Output, int &O2, the "order" of the polynomial derivative.
            //
            //    Output, double C2[O2], the coefficients of the polynomial 
            //    derivative.
            //
            //    Output, int E2[O2], the indices of the exponents of the
            //    polynomial derivative.
            //
        {
            int[] f1;
            int i;
            int j;

            o2 = o1;
            for (j = 0; j < o1; j++)
            {
                c2[j] = c1[j];
            }

            for (j = 0; j < o1; j++)
            {
                f1 = Monomial.mono_unrank_grlex(m, e1[j]);
                for (i = 0; i < m; i++)
                {
                    c2[j] = c2[j] * typeMethods.i4_fall(f1[i], dif[i]);
                    f1[i] = Math.Max(f1[i] - dif[i], 0);
                }

                e2[j] = Monomial.mono_rank_grlex(m, f1);
            }

            polynomial_sort(o2, ref c2, ref e2);

            polynomial_compress(o2, c2, e2, ref o2, ref c2, ref e2);

        }

        public static void polynomial_mul(int m, int o1, double[] c1, int[] e1, int o2, double[] c2,
                int[] e2, ref int o, ref double[] c, ref int[] e)

            //****************************************************************************80
            //
            //  Purpose:
            //
            //    POLYNOMIAL_MUL multiplies two polynomials.
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
            //    Input, int M, the spatial dimension.
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
            //    Output, int &O, the "order" of the polynomial product.
            //
            //    Output, double C[O], the coefficients of the polynomial product.
            //
            //    Output, int E[O], the indices of the exponents of the 
            //    polynomial product.
            //
        {
            int[] f;
            int[] f1;
            int[] f2;
            int i;
            int j;
            int k;

            f = new int[m];

            o = 0;
            for (j = 0; j < o2; j++)
            {
                for (i = 0; i < o1; i++)
                {
                    c[o] = c1[i] * c2[j];
                    f1 = Monomial.mono_unrank_grlex(m, e1[i]);
                    f2 = Monomial.mono_unrank_grlex(m, e2[j]);
                    for (k = 0; k < m; k++)
                    {
                        f[k] = f1[k] + f2[k];
                    }

                    e[o] = Monomial.mono_rank_grlex(m, f);
                    o = o + 1;
                }
            }

            polynomial_sort(o, ref c, ref e);
            polynomial_compress(o, c, e, ref o, ref c, ref e);
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
        
        public static void polynomial_scale ( double s, int m, int o, ref double[] c, int[] e )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLYNOMIAL_SCALE scales a polynomial.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 January 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double S, the scale factor.
        //
        //    Input, int M, the spatial dimension.
        //
        //    Input, int O, the "order" of the polynomial.
        //
        //    Input/output, double C[O], the coefficients of the polynomial.
        //
        //    Input, int E[O], the indices of the exponents of the
        //    polynomial.
        //
        {
            int i;

            for ( i = 0; i < o; i++ )
            {
                c[i] = c[i] * s;
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
            typeMethods.r8vec_permute(o, indx, ref c);
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

        public static double[] r8poly_values ( int m, double[] c, int n, double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    r8poly_values evaluates a polynomial using a naive method.
        //
        //  Discussion:
        //
        //    The polynomial 
        //
        //      p(x) = c0 + c1 * x + c2 * x^2 + ... + cm * x^m
        //
        //    is to be evaluated at the values X.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 April 2020
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Input:
        //
        //    int M, the degree of the polynomial.
        //
        //    double C[M+1], the coefficients of the polynomial.
        //    A[0] is the constant term.
        //
        //    double X[N], the points at which the polynomial is to be evaluated.
        //
        //  Output:
        //
        //    double R8POLY_VALUE[N], the values of the polynomial at X.
        //
        {
            int i;
            int j;
            double[] value;
            double xi;

            value = new double[n];

            for ( j = 0; j < n; j++ )
            {
                value[j] = c[0];
                xi = 1.0;
                for ( i = 1; i <= m; i++ )
                {
                    xi = xi * x[j];
                    value[j] = value[j] + c[i] * xi;
                }
            }

            return value;
        }
        
        public static double[] r8poly_values_2d ( int m, double[] c, int n, double[] x, double[] y )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    R8POLY_VALUES_2D evaluates a polynomial in 2 variables, X and Y.
        //
        //  Discussion:
        //
        //    We assume the polynomial is of total degree M, and has the form:
        //
        //      p(x,y) = c00 
        //             + c10 * x                + c01 * y
        //             + c20 * x^2   + c11 * xy + c02 * y^2
        //             + ...
        //             + cm0 * x^(m) + ...      + c0m * y^m.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    23 September 2012
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int M, the degree of the polynomial.
        //
        //    Input, double C[T(M+1)], the polynomial coefficients.  
        //    C[0] is the constant term.  T(M+1) is the M+1-th triangular number.
        //    The coefficients are stored consistent with the following ordering
        //    of monomials: 1, X, Y, X^2, XY, Y^2, X^3, X^2Y, XY^2, Y^3, X^4, ...
        //
        //    Input, int N, the number of evaluation points.
        //
        //    Input, double X[N], Y[N], the evaluation points.
        //
        //    Output, double R8POLY_VALUE_2D[N], the value of the polynomial at the 
        //    evaluation points.
        //
        {
            int ex;
            int ey;
            int i;
            int j;
            double[] p;
            int s;

            p = new double[n];

            for ( i = 0; i < n; i++ )
            {
                p[i] = 0.0;
            }

            j = 0;
            for ( s = 0; s <= m; s++ )
            {
                for ( ex = s; 0 <= ex; ex-- )
                {
                    ey = s - ex;
                    for ( i = 0; i < n; i++ )
                    {
                        p[i] = p[i] + c[j] * Math.Pow ( x[i], ex ) * Math.Pow ( y[i], ey );
                    }
                    j = j + 1;
                }
            }
            return p;
        }
    }
}