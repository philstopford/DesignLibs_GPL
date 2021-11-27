using System;
using Burkardt.FEM;

namespace FEM1DTest;

internal static class Program
{
    private static void Main()
//****************************************************************************80
//
//  Purpose:
//
//    MAIN is the main program for FEM1D.
//
//  Discussion:
//
//    FEM1D solves a one dimensional ODE using the finite element method.
//
//    The differential equation solved is
//
//      - d/dX (P dU/dX) + Q U  =  F
//
//    The finite-element method uses piecewise linear basis functions.
//
//    Here U is an unknown scalar function of X defined on the
//    interval [XL,XR], and P, Q and F are given functions of X.
//
//    The values of U or U' at XL and XR are also specified.
//
//
//    The interval [XL,XR] is "meshed" with NSUB+1 points,
//
//    XN(0) = XL, XN(1)=XL+H, XN(2)=XL+2*H, ..., XN(NSUB)=XR.
//
//    This creates NSUB subintervals, with interval number 1
//    having endpoints XN(0) and XN(1), and so on up to interval
//    NSUB, which has endpoints XN(NSUB-1) and XN(NSUB).
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    31 October 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    double ADIAG(NU), the "diagonal" coefficients.  That is, ADIAG(I) is
//    the coefficient of the I-th unknown in the I-th equation.
//
//    double ALEFT(NU), the "left hand" coefficients.  That is, ALEFT(I) 
//    is the coefficient of the (I-1)-th unknown in the I-th equation.
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
//    double H(NSUB)
//    H(I) is the length of subinterval I.  This code uses
//    equal spacing for all the subintervals.
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
//    int INDX[NSUB+1].
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
//    at node NSUB, which will be unknown NSUB or NSUB+1,
//    depending on whether there was an unknown at node 0.
//
//    int NL.
//    The number of basis functions used in a single
//    subinterval.  (NL-1) is the degree of the polynomials
//    used.  For this code, NL is fixed at 2, meaning that
//    piecewise linear functions are used as the basis.
//
//    int NODE[NL*NSUB].
//    For each subinterval I:
//    NODE[0+I*2] is the number of the left node, and
//    NODE[1+I*2] is the number of the right node.
//
//    int NQUAD.
//    The number of quadrature points used in a subinterval.
//    This code uses NQUAD = 1.
//
//    int NSUB.
//    The number of subintervals into which the interval [XL,XR] is broken.
//
//    int NU.
//    NU is the number of unknowns in the linear system.
//    Depending on the value of IBC, there will be NSUB-1,
//    NSUB, or NSUB+1 unknown values, which are the coefficients
//    of basis functions.
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
//    double XN(0:NSUB).
//    XN(I) is the location of the I-th node.  XN(0) is XL,
//    and XN(NSUB) is XR.
//
//    double XQUAD(NSUB)
//    XQUAD(I) is the location of the single quadrature point
//    in interval I.
//
//    double XR.
//    XR is the right endpoint of the interval over which the
//    differential equation is being solved.
//
    {
        const int NSUB = 5;
        const int NL = 2;

        double[] adiag = new double[NSUB + 1];
        double[] aleft = new double[NSUB + 1];
        double[] arite = new double[NSUB + 1];
        double[] f = new double[NSUB + 1];
        double[] h = new double[NSUB];
        int ibc = 0;
        int[] indx = new int[NSUB + 1];
        int[] node = new int[NL * NSUB];
        int nquad = 0;
        int nu = 0;
        double ul = 0;
        double ur = 0;
        double xl = 0;
        double[] xn = new double[NSUB + 1];
        double[] xquad = new double[NSUB];
        double xr = 0;
            
        Console.WriteLine("");
        Console.WriteLine("FEM1D");
        Console.WriteLine("");
        Console.WriteLine("  Solve the two-point boundary value problem");
        Console.WriteLine("");
        Console.WriteLine("  - d/dX (P dU/dX) + Q U  =  F");
        Console.WriteLine("");
        Console.WriteLine("  on the interval [XL,XR], specifying");
        Console.WriteLine("  the value of U or U' at each end.");
        Console.WriteLine("");
        Console.WriteLine("  The interval [XL,XR] is broken into NSUB = "
                          + NSUB + " subintervals");
        Console.WriteLine("  Number of basis functions per element is NL = "
                          + NL + "");
//
//  Initialize the data that defines the problem.
//
        FEM_1D.init(ref ibc, ref nquad, ref ul, ref ur, ref xl, ref xr);
//
//  Compute the quantities which define the geometry of the
//  problem.
//
        FEM_1D.geometry(h, ibc, ref indx, NL, ref node, NSUB, ref nu, xl, ref xn, ref xquad, xr);
//
//  Assemble the linear system.
//
        FEM_1D.assemble(ref adiag, ref aleft, ref arite, ref f, h, indx, NL, node, nu, nquad,
            NSUB, ul, ur, xn, xquad);
//
//  Print out the linear system.
//
        FEM_1D.prsys(adiag, aleft, arite, f, nu);
//
//  Solve the linear system.
//
        FEM_1D.solve(ref adiag, ref aleft, ref arite, ref f, nu);
//
//  Print out the solution.
//
        FEM_1D.output(f, ibc, indx, NSUB, nu, ul, ur, xn);
//
//  Terminate.
//
        Console.WriteLine("");
        Console.WriteLine("FEM1D:");
        Console.WriteLine("  Normal end of execution.");

        Console.WriteLine("");
    }
}