using System;
using Burkardt.FEM;

namespace Burkardt.FEM1DPMethodTest
{
    class Program
    {
        static void Main(string[] args)
//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for FEM1D_PMETHOD.
//
//  Discussion:
//
//    FEM1D_PMETHOD implements the P-version of the finite element method.
//
//    Program to solve the one dimensional problem:
//
//      - d/dX (P dU/dX) + Q U  =  F
//
//    by the finite-element method using a sequence of polynomials
//    which satisfy the boundary conditions and are orthogonal
//    with respect to the inner product:
//
//      (U,V)  =  Integral (-1 to 1) P U' V' + Q U V dx
//
//    Here U is an unknown scalar function of X defined on the
//    interval [-1,1], and P, Q and F are given functions of X.
//
//    The boundary values are U(-1) = U(1)=0.
//
//    Sample problem #1:
//
//      U=1-x^4,        P=1, Q=1, F=1.0+12.0*x^2-x^4
//
//    Sample problem #2:
//
//      U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x)
//
//    The program should be able to get the exact solution for
//    the first problem, using NP = 2.
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
//  Local Parameters:
//
//    Local, double A[NP+1], the squares of the norms of the 
//    basis functions.
//
//    Local, double ALPHA[NP], BETA[NP], the basis function 
//    recurrence coefficients.
//
//    Local, double F[NP+1].
//    F contains the basis function coefficients that form the
//    representation of the solution U.  That is,
//      U(X)  =  SUM (I=0 to NP) F(I) * BASIS(I)(X)
//    where "BASIS(I)(X)" means the I-th basis function
//    evaluated at the point X.
//
//    Local, int NP.
//    The highest degree polynomial to use.
//
//    Local, int NPRINT.
//    The number of points at which the computed solution
//    should be printed out at the end of the computation.
//
//    Local, int PROBLEM, indicates the problem being solved.
//    1, U=1-x^4, P=1, Q=1, F=1.0+12.0*x^2-x^4.
//    2, U=cos(0.5*pi*x), P=1, Q=0, F=0.25*pi*pi*cos(0.5*pi*x).
//
//    Local, int QUAD_NUM, the order of the quadrature rule.
//
//    Local, double QUAD_W[QUAD_NUM], the quadrature weights.
//
//    Local, double QUAD_X[QUAD_NUM], the quadrature abscissas.
//
        {
            int NP = 2;
            int QUAD_NUM = 10;

            double[] a = new double[NP + 1];
            double[] alpha = new double[NP];
            double[] beta = new double[NP];
            double[] f = new double[NP + 1];
            int nprint = 10;
            int problem = 2;
            double[] quad_w = new double[QUAD_NUM];
            double[] quad_x = new double[QUAD_NUM];

            Console.WriteLine("");
            Console.WriteLine("FEM1D_PMETHOD");
            Console.WriteLine("  C++ version");
            Console.WriteLine("");
            Console.WriteLine("  Solve the two-point boundary value problem");
            Console.WriteLine("");
            Console.WriteLine("  - d/dX (P dU/dX) + Q U  =  F");
            Console.WriteLine("");
            Console.WriteLine("  on the interval [-1,1], with");
            Console.WriteLine("  U(-1) = U(1) = 0.");
            Console.WriteLine("");
            Console.WriteLine("  The P method is used, which represents U as");
            Console.WriteLine("  a weighted sum of orthogonal polynomials.");
            Console.WriteLine("");
            Console.WriteLine("");
            Console.WriteLine("  Highest degree polynomial to use is " + NP + "");
            Console.WriteLine("  Number of points to be used for output = " + nprint + "");

            if (problem == 1)
            {
                Console.WriteLine("");
                Console.WriteLine("  Problem #1:");
                Console.WriteLine("  U=1-x^4,");
                Console.WriteLine("  P=1,");
                Console.WriteLine("  Q=1,");
                Console.WriteLine("  F=1 + 12 * x^2 - x^4");
            }
            else if (problem == 2)
            {
                Console.WriteLine("");
                Console.WriteLine("  Problem #2:");
                Console.WriteLine("  U=cos(0.5*pi*x),");
                Console.WriteLine("  P=1,");
                Console.WriteLine("  Q=0,");
                Console.WriteLine("  F=0.25*pi*pi*cos(0.5*pi*x)");
            }

//
//  Get quadrature abscissas and weights for interval [-1,1].
//
            FEM_1D_PMethod.quad(QUAD_NUM, ref quad_w, ref quad_x);
//
//  Compute the constants for the recurrence relationship
//  that defines the basis functions.
//
            FEM_1D_PMethod.alpbet(ref a, ref alpha, ref beta, NP, problem, QUAD_NUM, quad_w, quad_x);
//
//  Test the orthogonality of the basis functions.
//
            FEM_1D_PMethod.ortho(a, alpha, beta, NP, problem, QUAD_NUM, quad_w, quad_x);
//
//  Solve for the solution of the problem, in terms of coefficients
//  of the basis functions.
//
            FEM_1D_PMethod.sol(a, alpha, beta, ref f, NP, problem, QUAD_NUM, quad_w, quad_x);
//
//  Print out the solution, evaluated at each of the NPRINT points.
//
            FEM_1D_PMethod.out_(alpha, beta, f, NP, nprint);
//
//  Compare the computed and exact solutions.
//
            FEM_1D_PMethod.exact(alpha, beta, f, NP, nprint, problem, QUAD_NUM, quad_w, quad_x);
//
//  Terminate.
//
            Console.WriteLine("");
            Console.WriteLine("FEM1D_PMETHOD");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }
    }
}