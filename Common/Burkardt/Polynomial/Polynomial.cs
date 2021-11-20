using System;
using Burkardt.MonomialNS;
using Burkardt.Types;

namespace Burkardt.PolynomialNS;

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
        //    C++ version by John Burkardt.
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
        int k;

        double[] pols = new double[n + 1];

        double pkp1 = 1.0;
        pols[0] = pkp1;

        switch (n)
        {
            case 0:
                return pols;
        }

        double pk = pkp1;
        pkp1 = 0.5 * ((a - b) * y + (2.0 + a + b) * x);
        pols[1] = pkp1;

        switch (n)
        {
            case 1:
                return pols;
        }

        for (k = 2; k <= n; k++)
        {
            double k_r8 = k;

            double alpha = (2.0 * k_r8 + a + b - 1.0)
                           * (a - b) * (a + b) * y
                           + (2.0 * k_r8 + a + b - 1.0)
                           * (2.0 * k_r8 + a + b - 2.0)
                           * (2.0 * k_r8 + a + b) * x;

            double beta = 2.0 * (k_r8 + a - 1.0)
                              * (k_r8 + b - 1.0)
                              * (2.0 * k_r8 + a + b) * y * y;

            double pkm1 = pk;
            pk = pkp1;
            pkp1 = (alpha * pk - beta * pkm1)
                   / (2.0 * k_r8 * (k_r8 + a + b)
                      * (2.0 * k_r8 + a + b - 2.0));

            pols[k] = pkp1;
        }

        return pols;
    }

    public static double[] lege2eva(int degree, double[] z)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LEGE2EVA evaluates orthogonal polynomials on the symmetric square.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 July 2014
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
        //    C++ version by John Burkardt.
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
        //    Input, int DEGREE, the maximum degree of the polynomials.
        //
        //    Input, double Z[2], the evaluation point.
        //
        //    Output, double LEGE2EVA[(DEGREE+1)*(DEGREE+2)/2], the orthogonal
        //    polynomials evaluated at Z.
        //
    {
        int m;

        int npols = (degree + 1) * (degree + 2) / 2;
        double[] pols = new double[npols];

        double[] f1 = llegepols1(degree, z[0]);
        double[] f2 = llegepols1(degree, z[1]);

        int kk = 0;
        for (m = 0; m <= degree; m++)
        {
            int n;
            for (n = 0; n <= m; n++)
            {
                pols[kk] = f1[m - n] * f2[n];
                double scale = (1 + 2 * n) * (1 + 2 * (m - n));
                scale = 0.5 * Math.Sqrt(scale);
                pols[kk] *= scale;
                kk += 1;
            }
        }

        return pols;
    }

    public static double[] llegepols1(int degree, double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LLEGEPOLS1 evaluates orthogonal polynomials on the symmetric interval.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 July 2014
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
        //    C++ version by John Burkardt.
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
        //    Input, double X, the evaluation point.
        //
        //    Output, double LLEGEPOLS1[DEGREE+1], the orthogonal
        //    polynomials evaluated at X.
        //
    {
        int k;

        double[] pols = new double[degree + 1];
        double pkp1 = 1.0;
        pols[0] = pkp1;

        switch (degree)
        {
            case 0:
                return pols;
        }

        double pk = pkp1;
        pkp1 = x;
        pols[1] = pkp1;

        switch (degree)
        {
            case 1:
                return pols;
        }

        for (k = 1; k <= degree - 1; k++)
        {
            double pkm1 = pk;
            pk = pkp1;
            pkp1 = ((2 * k + 1) * x * pk
                    - k * pkm1)
                   / (k + 1);

            pols[k + 1] = pkp1;
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
        //    C++ version by John Burkardt.
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
        int i;
        int j;
        int m;
        int mmax;
        double x1;
        double y1;
        //
        //  Convert coordinates from reference to Koornwinder tetrahedron.
        //
        double[] uvw = typeMethods.ref_to_koorn(xyz);
        //
        //  Compute F1.
        //
        double p1 = 0.5 * (2.0 * uvw[0] + 2.0 + uvw[1] + uvw[2]);
        double p2 = -0.5 * (uvw[1] + uvw[2]);

        double[] f1 = LegendreScaled.klegeypols(p1, p2, degree);
        //
        //  Compute F2S.
        //
        double[] f2s = new double[(degree + 1) * (degree + 1)];

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
            p1 = 2 * m + 1;
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
        double[] f3s = new double[(degree + 1) * (degree + 1)];

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
            p1 = 2 * j;

            f = kjacoypols3(x1, y1, p1, p2, degree + 1 - j);

            for (i = 1; i <= degree + 2 - j; i++)
            {
                f3s[i - 1 + (j - 1) * (degree + 1)] = f[i - 1];
            }
        }

        //
        //  Construct FVALS.
        //
        int nvals = (degree + 1) * (degree + 2) * (degree + 3) / 6;
        double[] fvals = new double[nvals];

        int ncount = 0;

        for (mmax = 0; mmax <= degree; mmax++)
        {
            for (m = 0; m <= mmax; m++)
            {
                int n;
                for (n = 0; n <= mmax - m; n++)
                {
                    int k = mmax - m - n;

                    double scale = Math.Sqrt
                    (
                        4.0
                        / (2 * mmax + 3)
                        / (2 * m + 1)
                        / (n + m + 1)
                        / Math.Sqrt(2.0)
                    );

                    fvals[ncount] =
                        f1[m] *
                        f2s[n + m * (degree + 1)] *
                        f3s[k + (m + n) * (degree + 1)] / scale;

                    ncount += 1;
                }
            }
        }

        return fvals;
    }

    public static double ts_mult(double[] u, double h, int n)

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
        //    C++ version by John Burkardt.
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
        int k;

        double ts = 0.0;
        double hk = 1.0;
        for (k = 1; k <= n; k++)
        {
            ts += u[k] * hk;
            hk *= h;
        }

        return ts;
    }

    public static void polynomial_add(int o1, double[] c1, int[] e1, int o2, double[] c2,
            int[] e2, ref int o, ref double[] c, ref int[] e)

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
        typeMethods.r8vec_concatenate(o1, c1, o2, c2, ref c);
        typeMethods.i4vec_concatenate(o1, e1, o2, e2, ref e);

        polynomial_sort(o, ref c, ref e);
        polynomial_compress(o, c, e, ref o, ref c, ref e);
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
        int i;

        int o3 = o1 + o2;

        double[] c3 = new double[o3];
        int[] e3 = new int[o3];
        double[] sc1 = new double[o1];

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
        const double r8_epsilon_sqrt = 0.1490116119384766E-07;

        int get = 0;
        int put = 0;

        while (get < o1)
        {
            get += 1;

            switch (Math.Abs(c1[get - 1]))
            {
                case <= r8_epsilon_sqrt:
                    continue;
            }

            switch (put)
            {
                case 0:
                    put += 1;
                    c2[put - 1] = c1[get - 1];
                    e2[put - 1] = e1[get - 1];
                    break;
                default:
                {
                    if (e2[put - 1] == e1[get - 1])
                    {
                        c2[put - 1] += c1[get - 1];
                    }
                    else
                    {
                        put += 1;
                        c2[put - 1] = c1[get - 1];
                        e2[put - 1] = e1[get - 1];
                    }

                    break;
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
        int j;

        o2 = o1;
        for (j = 0; j < o1; j++)
        {
            c2[j] = c1[j];
        }

        for (j = 0; j < o1; j++)
        {
            int[] f1 = Monomial.mono_unrank_grlex(m, e1[j]);
            int i;
            for (i = 0; i < m; i++)
            {
                c2[j] *= typeMethods.i4_fall(f1[i], dif[i]);
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
        int j;

        int[] f = new int[m];

        o = 0;
        for (j = 0; j < o2; j++)
        {
            int i;
            for (i = 0; i < o1; i++)
            {
                c[o] = c1[i] * c2[j];
                int[] f1 = Monomial.mono_unrank_grlex(m, e1[i]);
                int[] f2 = Monomial.mono_unrank_grlex(m, e2[j]);
                int k;
                for (k = 0; k < m; k++)
                {
                    f[k] = f1[k] + f2[k];
                }

                e[o] = Monomial.mono_rank_grlex(m, f);
                o += 1;
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
        Console.WriteLine(title);

        switch (o)
        {
            case 0:
                Console.WriteLine("      0.");
                break;
            default:
            {
                int j;
                for (j = 0; j < o; j++)
                {
                    string cout = "    ";
                    cout += c[j] switch
                    {
                        < 0.0 => "- ",
                        _ => "+ "
                    };

                    cout += Math.Abs(c[j]) + " * x^(";

                    int[] f = Monomial.mono_unrank_grlex(d, e[j]);
                    int i;
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

                break;
            }
        }
    }

    public static void polynomial_scale(double s, int m, int o, ref double[] c, int[] e)

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

        for (i = 0; i < o; i++)
        {
            c[i] *= s;
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
        int j;
        int k;

        double[] p = new double[nx];

        for (k = 0; k < nx; k++)
        {
            p[k] = 0.0;
        }

        for (j = 0; j < o; j++)
        {
            int[] f = Monomial.mono_unrank_grlex(d, e[j]);
            double[] v = Monomial.mono_value(d, nx, f, x, xIndex);
            for (k = 0; k < nx; k++)
            {
                p[k] += c[j] * v[k];
            }
        }

        return p;
    }

    public static void poly(string code, ref int[] rexp, ref int[] sexp)

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
        switch (code)
        {
            case "Q4":
                poly_q4(ref rexp, ref sexp);
                break;
            case "Q8":
                poly_q8(ref rexp, ref sexp);
                break;
            case "Q9":
                poly_q9(ref rexp, ref sexp);
                break;
            case "Q12":
                poly_q12(ref rexp, ref sexp);
                break;
            case "Q16":
                poly_q16(ref rexp, ref sexp);
                break;
            case "QL":
                poly_ql(ref rexp, ref sexp);
                break;
            case "T3":
                poly_t3(ref rexp, ref sexp);
                break;
            case "T4":
                Console.WriteLine("");
                Console.WriteLine("POLY - Fatal error!");
                Console.WriteLine("  The T4 element does not follow the pattern!");
                break;
            case "T6":
                poly_t6(ref rexp, ref sexp);
                break;
            case "T10":
                poly_t10(ref rexp, ref sexp);
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("POLY - Fatal error!");
                Console.WriteLine("  Illegal value of CODE = " + code + "");
                break;
        }
    }

    public static void poly_q4(ref int[] rexp, ref int[] sexp)

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

    public static void poly_q8(ref int[] rexp, ref int[] sexp)

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

    public static void poly_q9(ref int[] rexp, ref int[] sexp)

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

    public static int poly_coef_count(int dim, int degree)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLY_COEF_COUNT: polynomial coefficient count given dimension and degree.
        //
        //  Discussion:
        //
        //    To count all monomials of degree 5 or less in dimension 3,
        //    we can count all monomials of degree 5 in dimension 4.
        //
        //    To count all monomials of degree 5 in dimension 4, we imagine
        //    that each of the variables X, Y, Z and W is a "box" and that
        //    we need to drop 5 pebbles into these boxes.  Every distinct
        //    way of doing this represents a degree 5 monomial in dimension 4.
        //    Ignoring W gives us monomials up to degree five in dimension 3.
        //
        //    To count them, we draw 3 lines as separators to indicate the
        //    4 boxes, and then imagine all disctinct sequences involving
        //    the three lines and the 5 pebbles.  Indicate the lines by 1's
        //    and the pebbles by 0's and we're asking for the number of
        //    permutations of 3 1's and 5 0's, which is 8! / (3! 5!)
        //
        //    In other words, 56 = 8! / (3! 5!) is:
        //    * the number of monomials of degree exactly 5 in dimension 4, 
        //    * the number of monomials of degree 5 or less in dimension 3, 
        //    * the number of polynomial coefficients of a polynomial of 
        //      degree 5 in (X,Y,Z).
        //
        //    In general, the formula for the number of monomials of degree DEG
        //    or less in dimension DIM is
        //
        //      (DEG+DIM)! / (DEG! * DIM!)
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    22 June 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM, the dimension of the polynomial.
        //    0 <= DIM.
        //
        //    Input, int DEGREE, the degree of the polynomnial
        //    0 <= DEGREE
        //
        //    Output, int POLY_COEF_COUNT, the number of coefficients 
        //    in the general polynomial of dimension DIM and degree DEGREE.
        //
    {
        int value;

        switch (dim)
        {
            case < 0:
                value = -1;
                break;
            default:
            {
                value = degree switch
                {
                    < 0 => -1,
                    _ => typeMethods.i4_choose(degree + dim, degree)
                };

                break;
            }
        }

        return value;
    }

    public static void poly_q12(ref int[] rexp, ref int[] sexp)

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

    public static void poly_q16(ref int[] rexp, ref int[] sexp)

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

    public static void poly_ql(ref int[] rexp, ref int[] sexp)

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

    public static void poly_t3(ref int[] rexp, ref int[] sexp)

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

    public static void poly_t6(ref int[] rexp, ref int[] sexp)

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

    public static void poly_t10(ref int[] rexp, ref int[] sexp)

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

    public static double[] r8poly_values(int m, double[] c, int n, double[] x)

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
        int j;

        double[] value = new double[n];

        for (j = 0; j < n; j++)
        {
            value[j] = c[0];
            double xi = 1.0;
            int i;
            for (i = 1; i <= m; i++)
            {
                xi *= x[j];
                value[j] += c[i] * xi;
            }
        }

        return value;
    }

    public static double[] r8poly_values_2d(int m, double[] c, int n, double[] x, double[] y)

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
        int i;
        int s;

        double[] p = new double[n];

        for (i = 0; i < n; i++)
        {
            p[i] = 0.0;
        }

        int j = 0;
        for (s = 0; s <= m; s++)
        {
            int ex;
            for (ex = s; 0 <= ex; ex--)
            {
                int ey = s - ex;
                for (i = 0; i < n; i++)
                {
                    p[i] += c[j] * Math.Pow(x[i], ex) * Math.Pow(y[i], ey);
                }

                j += 1;
            }
        }

        return p;
    }

    public static int[] itop ( ref PLY ply, int in_ )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ITOP converts an integer to a polynomial in the field of order Q.
        //
        //  Discussion:
        //
        //    A nonnegative integer IN can be decomposed into a polynomial in
        //    powers of Q, with coefficients between 0 and P-1, by setting:
        //
        //      J = 0
        //      do while ( 0 < IN )
        //        POLY(J) = mod ( IN, Q )
        //        J = J + 1
        //        IN = IN / Q
        //      end do
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 September 2007
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
        //    Programs to Generate Niederreiter's Low-Discrepancy Sequences,
        //    ACM Transactions on Mathematical Software,
        //    Volume 20, Number 4, 1994, pages 494-495.
        //
        //  Parameters:
        //
        //    Input, int IN, the (nonnegative) integer containing the 
        //    polynomial information.
        //
        //    Output, int ITOP[DEG_MAX+2], the polynomial information.
        //    ITOP[0] contains the degree of the polynomial.  ITOP[I+1] contains
        //    the coefficient of degree I.  Each coefficient is an element of
        //    the field of order Q; in other words, each coefficient is
        //    between 0 and Q-1.
        //
    {
        int j;

        int[] poly = new int[PLY.DEG_MAX+2];

        for ( j = 0; j < PLY.DEG_MAX + 2; j++ )
        {
            poly[j] = 0;
        }

        int i = in_;
        j = -1;

        while ( 0 < i )
        {
            j += 1;

            if ( PLY.DEG_MAX < j )
            {
                Console.WriteLine("");
                Console.WriteLine("ITOP - Fatal error!");
                Console.WriteLine("  The polynomial degree exceeds DEG_MAX.");
                return null;
            }
            poly[j+1] = i % ply.Q;
            i /= ply.Q;
        }

        poly[0] = j;

        return poly;
    }

    public static int[] itop(int in_, int p)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ITOP converts an integer to a polynomial in the field of order P.
        //
        //  Discussion:
        //
        //    A nonnegative integer IN can be decomposed into a polynomial in
        //    powers of P, with coefficients between 0 and P-1, by setting:
        //
        //      J = 0
        //      do while ( 0 < IN )
        //        POLY(J) = mod ( IN, P )
        //        J = J + 1
        //        IN = IN / P
        //      end do
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
        //    Paul Bratley, Bennet Fox, Harald Niederreiter.
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
        //    Input, int IN, the (nonnegative) integer containing the 
        //    polynomial information.
        //
        //    Input, int P, the order of the field.
        //
        //    Output, int ITOP[DEG_MAX+2], the polynomial information.
        //    ITOP[0] contains the degree of the polynomial.  ITOP[I+1] contains
        //    the coefficient of degree I.  Each coefficient is an element of
        //    the field of order P; in other words, each coefficient is
        //    between 0 and P-1.
        //
    {
        int j;

        int[] poly = new int[PLY.DEG_MAX + 2];

        for (j = 0; j < PLY.DEG_MAX + 2; j++)
        {
            poly[j] = 0;
        }

        int i = in_;
        j = -1;

        while (0 < i)
        {
            j += 1;

            if (PLY.DEG_MAX < j)
            {
                Console.WriteLine("");
                Console.WriteLine("ITOP - Fatal error!");
                Console.WriteLine("  The polynomial degree exceeds DEG_MAX.");
                return null;
            }

            poly[j + 1] = i % p;
            i /= p;
        }

        poly[0] = j;

        return poly;
    }
        
    public class PLY
    {
        public const int DEG_MAX = 50;
        public int P;
        public int Q;
        public const int Q_MAX = 50;

        public int[,] add = new int[Q_MAX, Q_MAX];
        public int[,] mul = new int[Q_MAX, Q_MAX];
        public int[,] sub = new int[Q_MAX, Q_MAX];
    }

    public static int[] plyadd(ref PLY ply, int[] pa, int[] pb)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PLYADD adds two polynomials.
        //
        //  Discussion:
        //
        //    POLY[0] contains the degree of the polynomial.  POLY[I+1] contains
        //    the coefficient of degree I.  Each coefficient is an element of
        //    the field of order Q; in other words, each coefficient is
        //    between 0 and Q-1.
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
        //    Paul Bratley, Bennet Fox, Harald Niederreiter.
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
        //    Output, int PLYADD[DEG_MAX+2], the sum polynomial.
        //
    {
        int i;

        int[] pc = new int[PLY.DEG_MAX + 2];

        int maxab = Math.Max(pa[0], pb[0]);

        int degc = -1;

        for (i = 0; i <= maxab; i++)
        {
            pc[i + 1] = ply.add[pa[i + 1], pb[i + 1]];

            if (pc[i + 1] != 0)
            {
                degc = i;
            }
        }

        pc[0] = degc;

        for (i = maxab + 1; i <= PLY.DEG_MAX; i++)
        {
            pc[i + 1] = 0;
        }

        return pc;
    }

    public static void plydiv(ref PLY ply, int[] pa, int[] pb, ref int[] pq, ref int[] pr)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PLYDIV divides one polynomial by another.
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
        //    Paul Bratley, Bennet Fox, Harald Niederreiter.
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
        //    Output, int PQ[DEG_MAX+2], the quotient polynomial.
        //
        //    Output, int PR[DEG_MAX+2], the remainder polynomial.
        //
    {
        int binv = 0;
        int d;
        int degq = 0;
        int i;

        switch (pb[0])
        {
            case -1:
                Console.WriteLine("");
                Console.WriteLine("PLYDIV -  Fatal error!");
                Console.WriteLine("  Division by zero polynomial.");
                return;
        }

        for (i = -1; i <= PLY.DEG_MAX; i++)
        {
            pq[i + 1] = 0;
            pr[i + 1] = pa[i + 1];
        }

        int degr = pa[0];
        int degb = pb[0];

        degq = degq switch
        {
            < 0 => -1,
            _ => degr - degb
        };

        //
        //  Find the inverse of the leading coefficient of PB.
        //
        int j = pb[degb + 1];

        for (i = 1; i <= ply.P - 1; i++)
        {
            binv = ply.mul[i, j] switch
            {
                1 => i,
                _ => binv
            };
        }

        for (d = degq; 0 <= d; d--)
        {
            int m = ply.mul[pr[degr + 1], binv];
            for (i = degb; 0 <= i; i--)
            {
                pr[degr + i - degb + 1] = ply.sub[pr[degr + i - degb + 1], ply.mul[m, pb[i + 1]]];
            }

            degr -= 1;
            pq[d + 1] = m;
        }

        pq[0] = degq;

        while (pr[degr + 1] == 0 && 0 <= degr)
        {
            degr -= 1;
        }

        pr[0] = degr;

    }

    public static int[] plymul(ref PLY ply, int[] pa, int[] pb)

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
        //    Paul Bratley, Bennet Fox, Harald Niederreiter.
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
        int degc;
        int i;

        int[] pc = new int[PLY.DEG_MAX + 2];

        int dega = pa[0];
        int degb = pb[0];

        if (dega == -1 || degb == -1)
        {
            degc = -1;
        }
        else
        {
            degc = dega + degb;
        }

        if (PLY.DEG_MAX < degc)
        {
            Console.WriteLine("");
            Console.WriteLine("PLYMUL - Fatal error!");
            Console.WriteLine("  The degree of the product exceeds DEG_MAX.");
            return null;
        }

        for (i = 0; i <= degc; i++)
        {
            int term = 0;
            int j;
            for (j = Math.Max(0, i - dega); j <= Math.Min(degb, i); j++)
            {
                term = ply.add[term, ply.mul[pa[i - j + 1], pb[j + 1]]];
            }

            pc[i + 1] = term;
        }

        pc[0] = degc;

        for (i = degc + 1; i <= PLY.DEG_MAX; i++)
        {
            pc[i + 1] = 0;
        }

        return pc;
    }

    public static int ptoi(int[] poly, int q)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PTOI converts a polynomial in the field of order Q to an integer.
        //
        //  Discussion:
        //
        //    A polynomial with coefficients A(*) in the field of order Q
        //    can also be stored in an integer I, with
        //
        //      I = AN*Q**N + ... + A0.
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
        //    Paul Bratley, Bennet Fox, Harald Niederreiter.
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
        //    Input, int POLY[DEG_MAX+2], the polynomial information.
        //    POLY[0] contains the degree of the polynomial.  POLY[I] contains
        //    the coefficient of degree I-1.  Each coefficient is an element of
        //    the field of order Q; in other words, each coefficient is
        //    between 0 and Q-1.
        //
        //    Input, int Q, the order of the field.
        //
        //    Output, int PTOI, the (nonnegative) integer containing the 
        //    polynomial information.
        //
    {
        int j;

        int degree = poly[0];

        int i = 0;
        for (j = degree; 0 <= j; j--)
        {
            i = i * q + poly[j + 1];
        }

        return i;
    }
        
    public static int ptoi ( ref PLY ply, int[] poly )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PTOI converts a polynomial in the field of order Q to an integer.
        //
        //  Discussion:
        //
        //    A polynomial with coefficients A(*) in the field of order Q
        //    can also be stored in an integer I, with
        //
        //      I = AN*Q**N + ... + A0.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    09 September 2007
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
        //    Programs to Generate Niederreiter's Low-Discrepancy Sequences,
        //    ACM Transactions on Mathematical Software,
        //    Volume 20, Number 4, 1994, pages 494-495.
        //
        //  Parameters:
        //
        //    Input, int POLY[DEG_MAX+2], the polynomial information.
        //    POLY[0] contains the degree of the polynomial.  POLY[I] contains
        //    the coefficient of degree I-1.  Each coefficient is an element of
        //    the field of order Q; in other words, each coefficient is
        //    between 0 and Q-1.
        //
        //    Output, int PTOI, the (nonnegative) integer containing the 
        //    polynomial information.
        //
    {
        int j;

        int degree = poly[0];
  
        int i = 0;
        for ( j = degree; 0 <= j; j-- )
        {
            i = i * ply.Q + poly[j+1];
        }

        return i;
    }

}