using System;
using Burkardt.FEM;

namespace FEM1DAdaptiveTest
{
    class Program
    {
        static void Main(string[] args)
            //****************************************************************************80
            //
            //  Purpose:
            //
            //    MAIN is the main program for FEM1D_ADAPTIVE.
            //
            //  Discussion:
            //
            //    FEM1D_ADAPTIVE solves a 1D problem using an adaptive finite element method.
            //
            //    The equation to be treated is:
            //
            //      -d/dx ( P(x) * dU(x)/dx ) + Q(x) * U(x)  =  F(x)
            //
            //    by the finite-element method using piecewise linear basis
            //    functions.
            //
            //    An adaptive method is used to try to reduce the maximum
            //    error by refining the mesh in certain places.
            //
            //    Here U is an unknown scalar function of X defined on the
            //    interval [XL,XR], and P, Q and F are given functions of X.
            //
            //    The values of U at XL and XR are also specified.
            //
            //    The interval [XL,XR] is "meshed" with N+1 points,
            //
            //      XN(0) = XL, XN(1) = XL+H, XN(2) = XL+2*H, ..., XN(N) = XR.
            //
            //    This creates N subintervals, with interval I having endpoints 
            //    XN(I-1) and XN(I).
            //
            //
            //    The algorithm tries to guarantee a certain amount
            //    of accuracy by examining the current solution, estimating the error
            //    in each subinterval, and, if necessary, subdividing one or more
            //    subintervals and repeating the calculation.
            //
            //    We can think of the adaptive part of the algorithm as a refined
            //    problem.  The program re-solves the problem on the pair of
            //    intervals J and J+1, which extend from node J-1 to node J+1.
            //    The values of U that were just computed at nodes J-1 and J+1
            //    will be used as the boundary values for this refined problem.
            //    The intervals J and J+1 will each be evenly divided into NY
            //    smaller subintervals.  This boundary value problem is solved,
            //    and the derivatives of the original and refined solutions are
            //    then compared to get an estimate of the error.
            //
            //  Licensing:
            //
            //    This code is distributed under the GNU LGPL license. 
            //
            //  Modified:
            //
            //    04 November 2006
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
            //    double ETA(N).
            //    ETA(I) is the error estimate for interval I.  It is computed
            //    as the sum of two quantities, one associated with the left
            //    and one with the right node of the interval.
            //
            //    double F(NU).
            //    ASSEMBLE stores into F the right hand side of the linear
            //    equations.
            //    SOLVE replaces those values of F by the solution of the
            //    linear equations.
            //
            //    double FY(M).
            //    FY is the right hand side of the linear system of the refined
            //    problem.
            //
            //    double H(N)
            //    H(I) is the length of subinterval I.  This code uses
            //    equal spacing for all the subintervals.
            //
            //    double HY(M).
            //    HY(I) is the length of subinterval I in the refined problem.
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
            //    int IBCY.
            //    IBCY declares the boundary conditions for the refined problem
            //    which should always be that the value of U is specified at
            //    both the left and right endpoints.  This corresponds to a
            //    value of IBCY = 3.
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
            //    int INDY(0:M).
            //    INDY(I) records the index of the unknown associated with
            //    node I for the refined problem.
            //
            //    int JADD(N).
            //    JADD(I) is 1 if the error estimates show that interval I
            //    should be subdivided.
            //
            //    int KOUNT, the number of adaptive steps that have been taken.
            //
            //    int M.
            //    M is the number of subintervals used in the refined problem.
            //    M is equal to NY for computations centered at node 0 or node N,
            //    and otherwise, M is equal to 2*NY.
            //
            //    int N
            //    The number of subintervals into which the interval
            //    [XL,XR] is broken.
            //
            //    int NL.
            //    The number of basis functions used in a single
            //    subinterval.  (NL-1) is the degree of the polynomials
            //    used.  For this code, NL is fixed at 2, meaning that
            //    piecewise linear functions are used as the basis.
            //
            //    int NMAX, the maximum number of unknowns that can be handled.
            //
            //    int NODE(NL,N).
            //    For each subinterval I:
            //    NODE(1,I) is the number of the left node, and
            //    NODE(2,I) is the number of the right node.
            //
            //    int NODEY(NL,M).
            //    NODEY performs the same function for the refined problem that
            //    NODE performs for the full problem, recording the node numbers
            //    associated with a particular subinterval.
            //
            //    int NQUAD
            //    The number of quadrature points used in a subinterval.
            //
            //    int NU.
            //    NU is the number of unknowns in the linear system.
            //    Depending on the value of IBC, there will be N-1,
            //    N, or N+1 unknown values, which are the coefficients
            //    of basis functions.
            //
            //    int NUY.
            //    The number of unknowns in the refined problem.
            //
            //    int NY.
            //    NY is the number of subintervals into which a given interval
            //    will be subdivided, before solving the refined probelm.
            //
            //    int PROBLEM, chooses the problem to be solved.
            //    The user must choose this value by setting it in routine GETPRB.
            //    * 1, u = x, p = 1, q = 0, f = 0, ibc = 3, ul = 0, ur = 1.
            //    The program should find the solution exactly, and the
            //    adaptive code should find that there is no reason to
            //    subdivide any interval.
            //    * 2, u = x*x, p = 1, q = 0, f = -2, ibc = 3, ul = 0, ur = 1.
            //    This problem should find the solution exactly, and
            //    the adaptive code should again find there is nothing
            //    to do.
            //    *3, u = sin(pi*x/2), p = 1, q = 0, ibc = 3, f = 0.25*pi*pi*sin(pi*x/2), 
            //    ul = 0, ur = 1.
            //    *4, u = cos(pi*x/2), p = 1, q = 0, ibc = 3, f = 0.25*pi*pi*cos(pi*x/2), 
            //    ul = 1, ur = 0.
            //    *5: u = x**(beta+2)/((beta+2)*(beta+1)), p = 1, q = 1, ibc = 3, 
            //    f = -x**beta + (x**(beta+2))/((beta+2)*(beta+1)),
            //    ul = 0, ur = 1/((beta+2)*(beta+1))
            //    (beta must be greater than -2, and not equal to -1)
            //    *6: u = atan((x-0.5)/alpha), p = 1, q = 0, ibc = 3, 
            //    f =  2*alpha*(x-0.5) / (alpha**2 + (x-0.5)**2) **2,
            //    ul = u(0), ur = u(1)
            //
            //    int STATUS, reports status of subdivision.
            //    0, a new subdivision was carried out.
            //    1, no more subdivisions are needed.
            //    -1, no more subdivisions can be carried out.
            //
            //    double TOL.
            //    A tolerance that is used to determine whether the estimated
            //    error in an interval is so large that it should be subdivided
            //    and the problem solved again.
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
            //    double WQUAD(NQUAD).
            //    WQUAD(I) is the weight associated with the I-th point
            //    of an NQUAD point Gaussian quadrature rule.
            //
            //    double XL.
            //    XL is the left endpoint of the interval over which the
            //    differential equation is being solved.
            //
            //    double XN(0:N).
            //    XN(I) is the location of the I-th node.  XN(0) is XL,
            //    and XN(N) is XR.
            //
            //    double XQUAD(NQUAD,NMAX), the I-th quadrature point
            //    in interval J.
            //
            //    double XQUADY(NQUAD,NMAY ), the I-th quadrature point
            //    in subinterval J of the refined problem.
            //
            //    double XR.
            //    XR is the right endpoint of the interval over which the
            //    differential equation is being solved.
            //
            //    Workspace, double precision XT(0:NMAX), used to compute a new
            //    set of nodes.
            //
            //    double YN(0:M).
            //    YN(I) is the location of the I-th node in the refined
            //    problem.
            //
        {
            int NL = 2;
            int NMAX = 30;
            int NQUAD = 2;

            double[] adiag = new double[NMAX + 1];
            double[] aleft = new double[NMAX + 1];
            double[] arite = new double[NMAX + 1];
            double[] eta = new double[NMAX];
            double[] f = new double[NMAX + 1];
            double[] h = new double[NMAX];
            int ibc = 0;
            int[] indx = new int[NMAX + 1];
            int[] node = new int[NL * NMAX];
            double tol = 0;
            double ul = 0;
            double ur = 0;
            double[] wquad = new double[NQUAD];
            double xl = 0;
            double[] xn = new double[NMAX + 1];
            double[] xquad = new double[NQUAD * NMAX];
            double xr = 0;

            Console.WriteLine("");
            Console.WriteLine("FEM1D_ADAPTIVE");
            Console.WriteLine("");
            Console.WriteLine("Solve the two-point boundary value problem:");
            Console.WriteLine("");
            Console.WriteLine("  -d/dx ( P(x) * dU(x)/dx ) + Q(x) * U(x)  =  F(x)");
            Console.WriteLine("");
            Console.WriteLine("on the interval [0,1], specifying the value");
            Console.WriteLine("of U at each endpoint.");
            Console.WriteLine("");
            Console.WriteLine("  The number of basis functions per element is " + NL + "");
            Console.WriteLine("");
            Console.WriteLine("  The number of quadrature points per element is " + NQUAD + "");

            int problem = FEM_1D_Adaptive.get_problem();

            Console.WriteLine("");
            Console.WriteLine("  Problem index = " + problem + "");
            Console.WriteLine("");

            if (problem == 1)
            {
                Console.WriteLine("");
                Console.WriteLine("  \"Linear\" problem:");
                Console.WriteLine("  (No refinement needed)");
                Console.WriteLine("");
                Console.WriteLine("  U(X) =  X");
                Console.WriteLine("  P(X) =  1.0");
                Console.WriteLine("  Q(X) =  0.0");
                Console.WriteLine("  F(X) =  0.0");
                Console.WriteLine("  IBC  =  3");
                Console.WriteLine("  UL   =  0.0");
                Console.WriteLine("  UR   =  1.0");
            }
            else if (problem == 2)
            {
                Console.WriteLine("");
                Console.WriteLine("  \"Quadratic\" problem:");
                Console.WriteLine("  (No refinement needed)");
                Console.WriteLine("");
                Console.WriteLine("  U(X) =  X*X");
                Console.WriteLine("  P(X) =  1.0");
                Console.WriteLine("  Q(X) =  0.0");
                Console.WriteLine("  F(X) = -2.0");
                Console.WriteLine("  IBC  =  3");
                Console.WriteLine("  UL   =  0.0");
                Console.WriteLine("  UR   =  1.");
            }
            else if (problem == 3)
            {
                Console.WriteLine("");
                Console.WriteLine("  \"SINE\" problem:");
                Console.WriteLine("");
                Console.WriteLine("  U(X) =  SIN(PI*X/2)");
                Console.WriteLine("  P(X) =  1.0");
                Console.WriteLine("  Q(X) =  0.0");
                Console.WriteLine("  F(X) =  PI*PI*SIN(PI*X/2)/4");
                Console.WriteLine("  IBC  =  3");
                Console.WriteLine("  UL   =  0.0");
                Console.WriteLine("  UR   =  1.0");
            }
            else if (problem == 4)
            {
                Console.WriteLine("");
                Console.WriteLine("  \"COSINE\" problem:");
                Console.WriteLine("");
                Console.WriteLine("  U(X) =  COS(PI*X/2)");
                Console.WriteLine("  P(X) =  1.0");
                Console.WriteLine("  Q(X) =  0.0");
                Console.WriteLine("  F(X) =  PI*PI*COS(PI*X/2)/4");
                Console.WriteLine("  IBC  =  3");
                Console.WriteLine("  UL   =  0.0");
                Console.WriteLine("  UR   =  1.0");
            }
            else if (problem == 5)
            {
                double beta = FEM_1D_Adaptive.get_beta();

                Console.WriteLine("");
                Console.WriteLine("  \"RHEINBOLDT\" problem:");
                Console.WriteLine("");
                Console.WriteLine("  U(X) =  X**(B+2)/((B+2)*(B+1))");
                Console.WriteLine("  P(X) =  1.0");
                Console.WriteLine("  Q(X) =  1.0");
                Console.WriteLine("  F(X) =  -X**B+(X**B+2))/((B+2)*(B+1))");
                Console.WriteLine("  IBC  =  3");
                Console.WriteLine("  UL   =  0.0");
                Console.WriteLine("  UR   =  1/((B+2)*(B+1))");
                Console.WriteLine("  B    = " + beta + "");
            }
            else if (problem == 6)
            {
                double alpha = FEM_1D_Adaptive.get_alpha();

                Console.WriteLine("");
                Console.WriteLine("  \"ARCTAN\" problem:");
                Console.WriteLine("");
                Console.WriteLine("  U(X) =  ATAN((X-0.5)/A)");
                Console.WriteLine("  P(X) =  1.0");
                Console.WriteLine("  Q(X) =  0.0");
                Console.WriteLine("  F(X) =  2*A*(X-0.5)/(A**2+(X-0.5)**2)**2");
                Console.WriteLine("  IBC  =  3");
                Console.WriteLine("  UL   =  ATAN(-0.5/A)");
                Console.WriteLine("  UR   =  ATAN( 0.5/A)");
                Console.WriteLine("  A    = " + alpha + "");
            }

            //
            //  Start out with just 4 subintervals.
            //
            int n = 4;
            //
            //  Initialize values that define the problem.
            //
            FEM_1D_Adaptive.init(ref ibc, n, ref tol, ref ul, ref ur, ref xl, ref xn, ref xr);
            //
            //  Start the iteration counter off at 0.
            //
            int kount = 0;
            //
            //  Begin the next iteration.
            //
            for (;;)
            {
                kount = kount + 1;

                Console.WriteLine("");
                Console.WriteLine("  Begin new iteration with " + n + " nodes.");
                Console.WriteLine("");
                //
                //  Solve the regular problem.
                //
                int nu = 0;
                FEM_1D_Adaptive.solvex(adiag, aleft, arite, ref f, ref h, ibc, indx, kount, n, NL, NMAX,
                    node, NQUAD, ref nu, ul, ur, wquad, xn, xquad);
                //
                //  Solve N subproblems to get the error estimators.
                //
                FEM_1D_Adaptive.solvey(ref eta, f, h, n, nu, ul, ur, xn);
                //
                //  Examine the error estimators, and see how many intervals should
                //  be subdivided.
                //
                int status = FEM_1D_Adaptive.subdiv(eta, kount, ref n, NMAX, tol, ref xn);

                if (status != 0)
                {
                    break;
                }
            }

            Console.WriteLine("");
            Console.WriteLine("FEM1D_ADAPTIVE:");
            Console.WriteLine("  Normal end of execution.");

            Console.WriteLine("");
        }
    }
}