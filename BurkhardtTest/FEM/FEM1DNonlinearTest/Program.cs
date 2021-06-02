using System;
using Burkardt.FEM;

namespace Burkardt.FEM1DNonlinearTest
{
    class Program
    {
        static void Main(string[] args)
//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for FEM1D_NONLINEAR.
//
//  Discussion:
//
//    FEM1D_NONLINLEAR solves a nonlinear one dimensional boundary value problem.
//
//    The differential equation has the form:
//
//      -d/dx (p(x) du/dx) + q(x)*u + u*u' =  f(x)
//
//    The finite-element method uses piecewise linear basis functions.
//
//    Here U is an unknown scalar function of X defined on the
//    interval [XL,XR], and P, Q and F are given functions of X.
//
//    The values of U or U' at XL and XR are also specified.
//
//    Sample problem #1:
//
//    u(x)  = x,
//    p(x)  = 1,
//    q(x)  = 0,
//    f(x)  = x,
//    u(0)  = 0,
//    u'(1) = 1.
//    The code should solve this problem exactly.
//
//    Sample problem #2:
//
//    u(x)  = 2*(1-cos(0.5*pi*x))/pi,
//    p(x)  = 1,
//    q(x)  = 0,
//    f(x)  = -0.5*pi*cos(0.5*pi*x) + 2*sin(0.5*pi*x)*(1-cos(0.5*pi*x)/pi
//    u(0)  = 0,
//    u'(1) = 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    29 April 2007
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    double ADIAG(NU).
//    ADIAG(I) is the "diagonal" coefficient of the I-th
//    equation in the linear system.  That is, ADIAG(I) is
//    the coefficient of the I-th unknown in the I-th equation.
//
//    double ALEFT(NU).
//    ALEFT(I) is the "left hand" coefficient of the I-th
//    equation in the linear system.  That is, ALEFT(I) is the
//    coefficient of the (I-1)-th unknown in the I-th equation.
//    There is no value in ALEFT(1), since the first equation
//    does not refer to a "0-th" unknown.
//
//    double ARITE(NU).
//    ARITE(I) is the "right hand" coefficient of the I-th
//    equation in the linear system.  ARITE(I) is the coefficient
//    of the (I+1)-th unknown in the I-th equation.  There is
//    no value in ARITE(NU) because the NU-th equation does not
//    refer to an "NU+1"-th unknown.
//
//    double F(NU).
//    ASSEMBLE stores into F the right hand side of the linear
//    equations.
//    SOLVE replaces those values of F by the solution of the
//    linear equations.
//
//    double FOLD(NU).
//    FOLD contains the value of F from the previous iteration,
//    and is used in ASSEMBLE to add correction terms to the
//    matrix and right hand side.
//
//    double H(N), the length of the subintervals.  
//
//    int IBC.
//    IBC declares what the boundary conditions are.
//    1, at the left endpoint, U has the value UL,
//       at the right endpoint, U' has the value UR.
//    2, at the left endpoint, U' has the value UL,
//       at the right endpoint, U has the value UR.
//    3, at the left endpoint, U has the value UL,
//       and at the right endpoint, U has the value UR.
//    4, at the left endpoint, U' has the value UL,
//       at the right endpoint U' has the value UR.
//
//    int IMAX.
//    The number of Newton iterations to carry out.
//
//    int INDX(0:N).
//    For a node I, INDX(I) is the index of the unknown
//    associated with node I.
//    If INDX(I) is equal to -1, then no unknown is associated
//    with the node, because a boundary condition fixing the
//    value of U has been applied at the node instead.
//    Unknowns are numbered beginning with 1.
//    If IBC is 2 or 4, then there is an unknown value of U
//    at node 0, which will be unknown number 1.  Otherwise,
//    unknown number 1 will be associated with node 1.
//    If IBC is 1 or 4, then there is an unknown value of U
//    at node N, which will be unknown N or N+1,
//    depending on whether there was an unknown at node 0.
//
//    int NL.
//    The number of basis functions used in a single
//    subinterval.  (NL-1) is the degree of the polynomials
//    used.  For this code, NL is fixed at 2, meaning that
//    piecewise linear functions are used as the basis.
//
//    int NODE(NL,N).
//    For each subinterval I:
//    NODE(1,I) is the number of the left node, and
//    NODE(2,I) is the number of the right node.
//
//    int NPRINT.
//    The number of points at which the computed solution
//    should be printed out when compared to the exact solution.
//
//    int NQUAD.
//    The number of quadrature points used in a subinterval.
//    This code uses NQUAD = 1.
//
//    int N.
//    The number of subintervals into which the interval
//    [XL,XR] is broken.
//
//    int NU.
//    NU is the number of unknowns in the linear system.
//    Depending on the value of IBC, there will be N-1,
//    N, or N+1 unknown values, which are the coefficients
//    of basis functions.
//
//    int PROBLEM, indicates which problem to be solved.
//    * 1, u = x, p=1, q=0, f=x, u(0)=0, u'(1)=1.
//    * 2, u = 2*(1-cos(0.5*pi*x))/pi, p=1, q=0,
//    f = -0.5*pi*cos(0.5*pi*x) + 2*sin(0.5*pi*x)*(1-cos(0.5*pi*x)/pi
//    u(0) = 0, u'(1)=1.
//
//    double UL.
//    If IBC is 1 or 3, UL is the value that U is required
//    to have at X = XL.
//    If IBC is 2 or 4, UL is the value that U' is required
//    to have at X = XL.
//
//    double UR.
//    If IBC is 2 or 3, UR is the value that U is required
//    to have at X = XR.
//    If IBC is 1 or 4, UR is the value that U' is required
//    to have at X = XR.
//
//    double XL.
//    XL is the left endpoint of the interval over which the
//    differential equation is being solved.
//
//    double XN(0:N).
//    XN(I) is the location of the I-th node.  XN(0) is XL,
//    and XN(N) is XR.
//
//    double XQUAD(N)
//    XQUAD(I) is the location of the single quadrature point
//    in interval I.
//
//    double XR.
//    XR is the right endpoint of the interval over which the
//    differential equation is being solved.
//
        {
            int N = 10;
            int NL = 2;

            double[] adiag = new double[N + 1];
            double[] aleft = new double[N + 1];
            double[] arite = new double[N + 1];
            double[] f = new double[N + 1];
            double[] fold = new double[N + 1];
            double[] h = new double[N];
            int ibc = 0;
            int imax = 0;
            int[] indx = new int[N + 1];
            int[] node = new int[NL * N];
            int nprint = 0;
            int nquad = 0;
            int nu = 0;
            int problem = 0;
            double ul = 0;
            double ur = 0;
            double xl = 0;
            double[] xn = new double[N + 1];
            double[] xquad = new double[N];
            double xr = 0;

            Console.WriteLine("");
            Console.WriteLine("FEM1D_NONLINEAR");
            Console.WriteLine("");
            Console.WriteLine("  Solve a nonlinear boundary value problem:");
            Console.WriteLine("");
            Console.WriteLine("    -d/dx (p(x) du/dx) + q(x)*u + u*u' = f(x)");
            Console.WriteLine("");
            Console.WriteLine("  on an interval [xl,xr], with the values of");
            Console.WriteLine("  u or u' specified at xl and xr.");
            Console.WriteLine("");
            Console.WriteLine("  The interval [XL,XR] is broken into N = "
                 + N + " subintervals");
            Console.WriteLine("  Number of basis functions per element is NL = "
                 + NL + "");
//
//  Initialize variables that define the problem.
//
            FEM_1D_Nonlinear.init(ref ibc, ref imax, ref nprint, ref nquad, ref problem, ref ul, ref ur, ref xl, ref xr);
//
//  Compute the quantities that describe the geometry of the problem.
//
            FEM_1D_Nonlinear.geometry(ref h, ibc, ref indx, NL, ref node, N, ref nu, xl, ref xn, ref xquad, xr);
//
//  Initialize the "previous" solution to 0.
//
            for (int i = 0; i < nu; i++)
            {
                fold[i] = 0.0;
            }

//
//  Begin the iteration.
//
            for (int i = 1; i <= imax; i++)
            {
//
//  Is it time for full nonlinear Newton iteration?
//
                if (i <= 3)
                {
                    FEM_1D_Nonlinear.assemble_picard(ref adiag, ref aleft, ref arite, ref f, fold, h, indx, N, NL, node,
                        nquad, nu, problem, ul, ur, xn, xquad);
                }
                else
                {
                    FEM_1D_Nonlinear.assemble_newton(ref adiag, ref aleft, ref arite, ref f, fold, h, indx, N, NL, node,
                        nquad, nu, problem, ul, ur, xn, xquad);
                }

//
//  Print out the linear system, just once.
//
                if (i == 1)
                {
                    FEM_1D_Nonlinear.prsys(adiag, aleft, arite, f, nu);
                }

//
//  Solve the linear system.
//
                FEM_1D_Nonlinear.solve(ref adiag, ref aleft, ref arite, ref f, nu);
//
//  Print the current solution.
//
                FEM_1D_Nonlinear.output(f, ibc, indx, N, nu, ul, ur, xn);
//
//  Save a copy of the current solution in FOLD.
//
                for (int j = 0; j < nu; j++)
                {
                    fold[j] = f[j];
                }
            }

//
//  Compare the solution to the exact solution.
//
            FEM_1D_Nonlinear.compare(f, indx, N, NL, node, nprint, nu, problem, ul, ur, xl, xn, xr);
//
//  Terminate.
//
            Console.WriteLine("");
            Console.WriteLine("FEM1D_NONLINEAR:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }
    }
}