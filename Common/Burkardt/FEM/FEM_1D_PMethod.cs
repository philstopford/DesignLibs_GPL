using System;
using System.Globalization;

namespace Burkardt.FEM;

public static class FEM_1D_PMethod
{
    public static void alpbet(ref double[] a, ref double[] alpha, ref double[] beta, int np, int problem,
        int quad_num, double[] quad_w, double[] quad_x)
//****************************************************************************80
//
//  Purpose:
//
//    ALPBET calculates the coefficients in the recurrence relationship.
//
//  Discussion:
//
//    ALPHA and BETA are the coefficients in the three
//    term recurrence relation for the orthogonal basis functions
//    on [-1,1].
//
//    The routine also calculates A, the square of the norm of each basis
//    function.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 November 2006
//
//  Author:
//
//    Original FORTRAN77 version by Max Gunzburger, Teresa Hodge.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double A(0:NP), the squares of the norms of the 
//    basis functions.
//
//    Output, double ALPHA(NP), BETA(NP), the basis function 
//    recurrence coefficients.
//
//    Input, int NP.
//    The highest degree polynomial to use.
//
//    Input, int PROBLEM, indicates the problem being solved.
//    1, U=1-x^4, P=1, Q=1, F=1.0+12.0*x^2-x^4.
//    2, U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x).
//
//    Input, int QUAD_NUM, the order of the quadrature rule.
//
//    Input, double QUAD_W(QUAD_NUM), the quadrature weights.
//
//    Input, double QUAD_X(QUAD_NUM), the quadrature abscissas.
//
    {
        double s;
        double u;
        double x;

        double ss = 0.0;
        double su = 0.0;

        for (int iq = 0; iq < quad_num; iq++)
        {
            x = quad_x[iq];

            s = 4.0 * pp(x, problem) * x * x
                + qq(x, problem) * (1.0 - x * x) * (1.0 - x * x);

            u = 2.0 * pp(x, problem) * x * (3.0 * x * x - 1.0)
                + x * qq(x, problem) * (1.0 - x * x) * (1.0 - x * x);

            ss += s * quad_w[iq];
            su += u * quad_w[iq];
        }

        a[0] = ss;
        alpha[0] = su / ss;
        beta[0] = 0.0;

        for (int i = 1; i <= np; i++)
        {
            ss = 0.0;
            su = 0.0;
            double sv = 0.0;

            for (int iq = 0; iq < quad_num; iq++)
            {
                x = quad_x[iq];
//
//  Three term recurrence for Q.
//
                double qm1 = 0.0;
                double q = 1.0;
                for (int k = 0; k <= i - 1; k++)
                {
                    double qm2 = qm1;
                    qm1 = q;
                    q = (x - alpha[k]) * qm1 - beta[k] * qm2;
                }

//
//  Three term recurrence for Q'.
//
                double qm1x = 0.0;
                double qx = 0.0;
                for (int k = 0; k <= i - 1; k++)
                {
                    double qm2x = qm1x;
                    qm1x = qx;
                    qx = qm1 + (x - alpha[k]) * qm1x - beta[k] * qm2x;
                }

                double t = 1.0 - x * x;
//
//  The basis function PHI = ( 1 - x^2 ) * q.
//
//     s = pp * ( phi(i) )' * ( phi(i) )' + qq * phi(i) * phi(i)
//
                s = pp(x, problem) * Math.Pow(t * qx - 2.0 * x * q, 2)
                    + qq(x, problem) * t * t * q * q;
//
//     u = pp * ( x * phi(i) )' * phi(i)' + qq * x * phi(i) * phi(i)
//
                u = pp(x, problem)
                    * (x * t * qx + (1.0 - 3.0 * x * x) * q)
                    * (t * qx - 2.0 * x * q) + x * qq(x, problem)
                                                 * t * t * q * q;
//
//     v = pp * ( x * phi(i) )' * phi(i-1) + qq * x * phi(i) * phi(i-1)
//
                double v = pp(x, problem)
                           * (x * t * qx + (1.0 - 3.0 * x * x) * q)
                           * (t * qm1x - 2.0 * x * qm1)
                           + x * qq(x, problem) * t * t * q * qm1;
//
//  SS(i) = <   phi(i), phi(i)   > = integral ( S )
//  SU(i) = < x phi(i), phi(i)   > = integral ( U )
//  SV(i) = < x phi(i), phi(i-1) > = integral ( V )
//
                ss += s * quad_w[iq];
                su += u * quad_w[iq];
                sv += v * quad_w[iq];
            }

            a[i] = ss;
//
//  ALPHA(i) = SU(i) / SS(i)
//  BETA(i)  = SV(i) / SS(i-1)
//
            if (i >= np)
            {
                continue;
            }

            alpha[i] = su / ss;
            beta[i] = sv / a[i - 1];
        }
    }

    public static void exact(double[] alpha, double[] beta, double[] f, int np, int nprint,
        int problem, int quad_num, double[] quad_w, double[] quad_x)
//****************************************************************************80
//
//  Purpose:
//
//    EXACT compares the computed and exact solutions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 November 2006
//
//  Author:
//
//    Original FORTRAN77 version by Max Gunzburger, Teresa Hodge.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, double ALPHA(NP), BETA(NP), the basis function 
//    recurrence coefficients.
//
//    Input, double F(0:NP).
//    F contains the basis function coefficients that form the
//    representation of the solution U.  That is,
//      U(X)  =  SUM (I=0 to NP) F(I) * BASIS(I)(X)
//    where "BASIS(I)(X)" means the I-th basis function
//    evaluated at the point X.
//
//    Input, int NP.
//    The highest degree polynomial to use.
//
//    Input, int NPRINT.
//    The number of points at which the computed solution
//    should be printed out at the end of the computation.
//
//    Input, int PROBLEM, indicates the problem being solved.
//    1, U=1-x^4, P=1, Q=1, F=1.0+12.0*x^2-x^4.
//    2, U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x).
//
//    Input, int QUAD_NUM, the order of the quadrature rule.
//
//    Input, double QUAD_W(QUAD_NUM), the quadrature weights.
//
//    Input, double QUAD_X(QUAD_NUM), the quadrature abscissas.
//
    {
        const int nsub = 10;
        double phii = 0;
        double phiix = 0;
        double x;

        Console.WriteLine("");
        Console.WriteLine("  Comparison of computed and exact solutions:");
        Console.WriteLine("");
        Console.WriteLine("    X        U computed    U exact     Difference");
        Console.WriteLine("");

        for (int i = 0; i <= nprint; i++)
        {
            x = (2 * i - nprint) / (double) nprint;
            double ue = uex(x, problem);
            double up = 0.0;
            for (int j = 0; j <= np; j++)
            {
                phi(alpha, beta, j, np, ref phii, ref phiix, x);
                up += phii * f[j];
            }

            Console.WriteLine("  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + up.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + ue.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + (ue - up).ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

//
//  Compute the big L2 error.
//
        double big_l2 = 0.0;

        for (int i = 1; i <= nsub; i++)
        {
            double xl = (2 * i - nsub - 1) / (double) nsub;
            double xr = (2 * i - nsub) / (double) nsub;

            for (int j = 0; j < quad_num; j++)
            {
                x = (xl * (1.0 - quad_x[j])
                     + xr * (1.0 + quad_x[j])) / 2.0;

                double up = 0.0;
                for (int k = 0; k <= np; k++)
                {
                    phi(alpha, beta, k, np, ref phii, ref phiix, x);
                    up += phii * f[k];
                }

                big_l2 += Math.Pow(up - uex(x, problem), 2) * quad_w[j]
                                                            * (xr - xl) / 2.0;
            }
        }

        big_l2 = Math.Sqrt(big_l2);

        Console.WriteLine("");
        Console.WriteLine("  Big L2 error = " + big_l2 + "");

    }

    public static double ff(double x, int problem)
//****************************************************************************80
//
//  Purpose:
//
//    FF evaluates the right hand side function F(X) at any point X.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 November 2006
//
//  Author:
//
//    Original FORTRAN77 version by Max Gunzburger, Teresa Hodge.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//
//    Input, int PROBLEM, indicates the problem being solved.
//    1, U=1-x^4, P=1, Q=1, F=1.0+12.0*x^2-x^4.
//    2, U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x).
//
//    Output, double FF, the value of F(X).
//
    {
        double value = problem switch
        {
//
//  Test problem 1
//
            1 => 1.0 + 12.0 * x * x - x * x * x * x,
//
//  Test problem 2
//
            2 => 0.25 * Math.PI * Math.PI * Math.Cos(0.5 * Math.PI * x),
            _ => 0
        };

        return value;
    }

    public static void ortho(double[] a, double[] alpha, double[] beta, int np, int problem,
        int quad_num, double[] quad_w, double[] quad_x)
//****************************************************************************80
//
//  Purpose:
//
//    ORTHO tests the basis functions for orthogonality.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 November 2006
//
//  Author:
//
//    Original FORTRAN77 version by Max Gunzburger, Teresa Hodge.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, double A(0:NP), the squares of the norms of the 
//    basis functions.
//
//    Input, double ALPHA(NP), BETA(NP), the basis function 
//    recurrence coefficients.
//
//    Input, int NP.
//    The highest degree polynomial to use.
//
//    Input, int PROBLEM, indicates the problem being solved.
//    1, U=1-x^4, P=1, Q=1, F=1.0+12.0*x^2-x^4.
//    2, U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x).
//
//    Input, int QUAD_NUM, the order of the quadrature rule.
//
//    Input, double QUAD_W(QUAD_NUM), the quadrature weights.
//
//    Input, double QUAD_X(QUAD_NUM), the quadrature abscissas.
//
    {
        double phii = 0;
        double phiix = 0;
        double phij = 0;
        double phijx = 0;
//
//  Zero out the B array, so we can start summing up the dot products.
//
        double[] b = new double[(np + 1) * (np + 1)];

        for (int j = 0; j <= np; j++)
        {
            for (int i = 0; i <= np; i++)
            {
                b[i + j * (np + 1)] = 0.0;
            }
        }

//
//  Approximate the integral of the product of basis function
//  I and basis function J over the interval [-1,1].
//
//  We expect to get zero, except when I and J are equal,
//  when we should get A(I).
//
        for (int iq = 0; iq < quad_num; iq++)
        {
            double x = quad_x[iq];
            for (int i = 0; i <= np; i++)
            {
                phi(alpha, beta, i, np, ref phii, ref phiix, x);
                for (int j = 0; j <= np; j++)
                {
                    phi(alpha, beta, j, np, ref phij, ref phijx, x);

                    double bij = pp(x, problem) * phiix * phijx
                                 + qq(x, problem) * phii * phij;

                    b[i + j * (np + 1)] += bij * quad_w[iq];
                }
            }
        }

//
//  Print out the results of the test.
//
        Console.WriteLine("");
        Console.WriteLine("  Basis function orthogonality test:");
        Console.WriteLine("");
        Console.WriteLine("   i   j     b(i,j)/a(i)");
        Console.WriteLine("");
        for (int i = 0; i <= np; i++)
        {
            Console.WriteLine("");
            for (int j = 0; j <= np; j++)
            {
                Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                       + "  " + j.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                       + "  " + (b[i + j * (np + 1)] / a[i]).ToString(CultureInfo.InvariantCulture)
                                       .PadLeft(12) + "");
            }
        }
    }

    public static void out_(double[] alpha, double[] beta, double[] f, int np, int nprint)
//****************************************************************************80
//
//  Purpose:
//
//    OUT prints the computed solution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 November 2006
//
//  Author:
//
//    Original FORTRAN77 version by Max Gunzburger, Teresa Hodge.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, double ALPHA(NP), BETA(NP), the basis function 
//    recurrence coefficients.
//
//    Input, double F(0:NP).
//    F contains the basis function coefficients that form the
//    representation of the solution U.  That is,
//      U(X)  =  SUM (I=0 to NP) F(I) * BASIS(I)(X)
//    where "BASIS(I)(X)" means the I-th basis function
//    evaluated at the point X.
//
//    Input, int NP.
//    The highest degree polynomial to use.
//
//    Input, int NPRINT.
//    The number of points at which the computed solution
//    should be printed out at the end of the computation.
//
    {
        double phii = 0;
        double phiix = 0;

        Console.WriteLine("");
        Console.WriteLine("  Representation of solution:");
        Console.WriteLine("");
        Console.WriteLine("  Basis function coefficients:");
        Console.WriteLine("");
        for (int i = 0; i <= np; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + f[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

        Console.WriteLine("");
        Console.WriteLine("");
        Console.WriteLine("       X     Approximate Solution");
        Console.WriteLine("");
        for (int ip = 0; ip <= nprint; ip++)
        {
            double x = (2 * ip - nprint) / (double) nprint;
            double up = 0.0;
            for (int i = 0; i <= np; i++)
            {
                phi(alpha, beta, i, np, ref phii, ref phiix, x);
                up += phii * f[i];
            }

            Console.WriteLine("  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + up.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }

        Console.WriteLine("");
    }

    public static void phi(double[] alpha, double[] beta, int i, int np, ref double phii,
        ref double phiix, double x)
//****************************************************************************80
//
//  Purpose:
//
//    PHI evaluates the I-th basis function at the point X.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 November 2006
//
//  Author:
//
//    Original FORTRAN77 version by Max Gunzburger, Teresa Hodge.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, double ALPHA(NP), BETA(NP), the basis function 
//    recurrence coefficients.
//
//    Input, int I, the index of the basis function.
//
//    Input, int NP.
//    The highest degree polynomial to use.
//
//    Output, double PHII, PHIIX, the value of the basis
//    function and its derivative.
//
//    Input, double X, the evaluation point.
//
    {
        double t;

        double qm1 = 0.0;
        double q = 1.0;
        double qm1x = 0.0;
        double qx = 0.0;

        for (int j = 1; j <= i; j++)
        {
            double qm2 = qm1;
            qm1 = q;
            double qm2x = qm1x;
            qm1x = qx;
            t = x - alpha[j - 1];
            q = t * qm1 - beta[j - 1] * qm2;
            qx = qm1 + t * qm1x - beta[j - 1] * qm2x;
        }

        t = 1.0 - x * x;
        phii = t * q;
        phiix = t * qx - 2.0 * x * q;
    }

    public static double pp(double x, int problem)
//****************************************************************************80
//
//  Purpose:
//
//    PP returns the value of the coefficient function P(X).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 November 2006
//
//  Author:
//
//    Original FORTRAN77 version by Max Gunzburger, Teresa Hodge.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//
//    Input, int PROBLEM, indicates the problem being solved.
//    1, U=1-x^4, P=1, Q=1, F=1.0+12.0*x^2-x^4.
//    2, U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x).
//
//    Output, double PP, the value of P(X).
//
    {
        double value = 0;
        switch (problem)
        {

//
//  Test problem 1
//
            case 1:

//
//  Test problem 2
//
            case 2:
                value = 1.0;
                break;
        }

        return value;
    }

    public static void quad(int quad_num, ref double[] quad_w, ref double[] quad_x)
//****************************************************************************80
//
//  Purpose:
//
//    QUAD returns the abscissas and weights for gaussian quadrature on [-1,1].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 November 2006
//
//  Author:
//
//    Original FORTRAN77 version by Max Gunzburger, Teresa Hodge.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, int QUAD_NUM, the order of the quadrature rule.
//
//    Output, double QUAD_W(QUAD_NUM), the quadrature weights.
//
//    Output, double QUAD_X(QUAD_NUM), the quadrature abscissas.
//
    {
//
//  Quadrature points on [-1,1]
//
        quad_x[0] = -0.973906528517172;
        quad_x[1] = -0.865063366688985;
        quad_x[2] = -0.679409568299024;
        quad_x[3] = -0.433395394129247;
        quad_x[4] = -0.148874338981631;
        quad_x[5] = 0.148874338981631;
        quad_x[6] = 0.433395394129247;
        quad_x[7] = 0.679409568299024;
        quad_x[8] = 0.865063366688985;
        quad_x[9] = 0.973906528517172;
//
//  Weight factors
//
        quad_w[0] = 0.066671344308688;
        quad_w[1] = 0.149451349150581;
        quad_w[2] = 0.219086362515982;
        quad_w[3] = 0.269266719309996;
        quad_w[4] = 0.295524224714753;
        quad_w[5] = 0.295524224714753;
        quad_w[6] = 0.269266719309996;
        quad_w[7] = 0.219086362515982;
        quad_w[8] = 0.149451349150581;
        quad_w[9] = 0.066671344308688;

    }

    public static double qq(double x, int problem)
//****************************************************************************80
//
//  Purpose:
//
//    QQ returns the value of the coefficient function Q(X).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 November 2006
//
//  Author:
//
//    Original FORTRAN77 version by Max Gunzburger, Teresa Hodge.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//
//    Input, int PROBLEM, indicates the problem being solved.
//    1, U=1-x^4, P=1, Q=1, F=1.0+12.0*x^2-x^4.
//    2, U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x).
//
//    Output, double QQ, the value of Q(X).
//
    {
        double value = problem switch
        {
//
//  Test problem 1
//
            1 => 1.0,
//
//  Test problem 2
//
            2 => 0.0,
            _ => 0
        };

        return value;
    }

    public static void sol(double[] a, double[] alpha, double[] beta, ref double[] f, int np,
        int problem, int quad_num, double[] quad_w, double[] quad_x)
//****************************************************************************80
//
//  Purpose:
//
//    SOL solves a linear system for the finite element coefficients.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 November 2006
//
//  Author:
//
//    Original FORTRAN77 version by Max Gunzburger, Teresa Hodge.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, double A(0:NP), the squares of the norms of the 
//    basis functions.
//
//    Input, double ALPHA(NP), BETA(NP), the basis function 
//    recurrence coefficients.
//
//    Output, double F(0:NP).
//    F contains the basis function coefficients that form the
//    representation of the solution U.  That is,
//      U(X)  =  SUM (I=0 to NP) F(I) * BASIS(I)(X)
//    where "BASIS(I)(X)" means the I-th basis function
//    evaluated at the point X.
//
//    Input, int NP.
//    The highest degree polynomial to use.
//
//    Input, int PROBLEM, indicates the problem being solved.
//    1, U=1-x^4, P=1, Q=1, F=1.0+12.0*x^2-x^4.
//    2, U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x).
//
//    Input, int QUAD_NUM, the order of the quadrature rule.
//
//    Input, double QUAD_W(QUAD_NUM), the quadrature weights.
//
//    Input, double QUAD_X(QUAD_NUM), the quadrature abscissas.
//
    {
        double phii = 0;
        double phiix = 0;

        for (int i = 0; i <= np; i++)
        {
            f[i] = 0.0;
        }

        for (int iq = 0; iq < quad_num; iq++)
        {
            double x = quad_x[iq];
            double t = ff(x, problem) * quad_w[iq];
            for (int i = 0; i <= np; i++)
            {
                phi(alpha, beta, i, np, ref phii, ref phiix, x);
                f[i] += phii * t;
            }
        }

        for (int i = 0; i <= np; i++)
        {
            f[i] /= a[i];
        }
    }


    public static double uex(double x, int problem)
//****************************************************************************80
//
//  Purpose:
//
//    UEX returns the value of the exact solution at a point X.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 November 2006
//
//  Author:
//
//    Original FORTRAN77 version by Max Gunzburger, Teresa Hodge.
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//
//    Input, int PROBLEM, indicates the problem being solved.
//    1, U=1-x^4, P=1, Q=1, F=1.0+12.0*x^2-x^4.
//    2, U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x).
//
//    Output, double UEX, the exact value of U(X).
//
    {
        double value = problem switch
        {
//
//  Test problem 1
//
            1 => 1.0 - Math.Pow(x, 4),
//
//  Test problem 2
//
            2 => Math.Cos(0.5 * Math.PI * x),
            _ => 0
        };

        return value;
    }
}