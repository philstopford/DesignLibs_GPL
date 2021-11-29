using System;
using System.Globalization;

namespace Burkardt.FEM;

public static class FEM_1D_Adaptive
{
    public static void assemble(ref double[] adiag, ref double[] aleft, ref double[] arite, ref double[] f,
            double[] h, int n, int[] indx, int[] node, int nu, int nl, int nquad,
            int nmax, double ul, double ur, double[] wquad, double[] xn, double[] xquad)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ASSEMBLE assembles the global matrix.
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
        //    Output, double ADIAG(NU).
        //    ADIAG(I) is the "diagonal" coefficient of the I-th
        //    equation in the linear system.  That is, ADIAG(I) is
        //    the coefficient of the I-th unknown in the I-th equation.
        //
        //    Output, double ALEFT(NU).
        //    ALEFT(I) is the "left hand" coefficient of the I-th
        //    equation in the linear system.  That is, ALEFT(I) is the
        //    coefficient of the (I-1)-th unknown in the I-th equation.
        //    There is no value in ALEFT(1), since the first equation
        //    does not refer to a "0-th" unknown.
        //
        //    Output, double ARITE(NU).
        //    ARITE(I) is the "right hand" coefficient of the I-th
        //    equation in the linear system.  ARITE(I) is the coefficient
        //    of the (I+1)-th unknown in the I-th equation.  There is
        //    no value in ARITE(NU) because the NU-th equation does not
        //    refer to an "NU+1"-th unknown.
        //
        //    Output, double F(NU).
        //    ASSEMBLE stores into F the right hand side of the linear
        //    equations.
        //    SOLVE replaces those values of F by the solution of the
        //    linear equations.
        //
        //    Input, double H(N)
        //    H(I) is the length of subinterval I.  This code uses
        //    equal spacing for all the subintervals.
        //
        //    Input, int N
        //    The number of subintervals into which the interval
        //    [XL,XR] is broken.
        //
        //    Input, int INDX(0:N).
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
        //    Input, int NODE(NL,N).
        //    For each subinterval I:
        //    NODE(1,I) is the number of the left node, and
        //    NODE(2,I) is the number of the right node.
        //
        //    Input, int NU.
        //    NU is the number of unknowns in the linear system.
        //    Depending on the value of IBC, there will be N-1,
        //    N, or N+1 unknown values, which are the coefficients
        //    of basis functions.
        //
        //    Input, int NL.
        //    The number of basis functions used in a single
        //    subinterval.  (NL-1) is the degree of the polynomials
        //    used.  For this code, NL is fixed at 2, meaning that
        //    piecewise linear functions are used as the basis.
        //
        //    Input, int NQUAD
        //    The number of quadrature points used in a subinterval.
        //
        //    Input, int NMAX, the maximum number of unknowns that can be handled.
        //
        //    Input, double UL.
        //    If IBC is 1 or 3, UL is the value that U is required
        //    to have at X = XL.
        //    If IBC is 2 or 4, UL is the value that U' is required
        //    to have at X = XL.
        //
        //    Input, double UR.
        //    If IBC is 2 or 3, UR is the value that U is required
        //    to have at X = XR.
        //    If IBC is 1 or 4, UR is the value that U' is required
        //    to have at X = XR.
        //
        //    Input, double WQUAD(NQUAD).
        //    WQUAD(I) is the weight associated with the I-th point
        //    of an NQUAD point Gaussian quadrature rule.
        //
        //    Input, double XN(0:N).
        //    XN(I) is the location of the I-th node.  XN(0) is XL,
        //    and XN(N) is XR.
        //
        //    Input, double XQUAD(NQUAD,NMAX), the I-th quadrature point
        //    in interval J.
        //
    {
        double phii = 0;
        double phiix = 0;
        double phij = 0;
        double phijx = 0;
        //
        //  Zero out the entries.
        //
        for (int i = 0; i < nu; i++)
        {
            f[i] = 0.0;
            aleft[i] = 0.0;
            arite[i] = 0.0;
            adiag[i] = 0.0;
        }

        //
        //  For each interval,
        //
        for (int ie = 0; ie < n; ie++)
        {
            double he = h[ie];
            double xleft = xn[node[0 + ie * 2]];
            double xrite = xn[node[1 + ie * 2]];
            //
            //  For each quadrature point in the interval,
            //
            for (int iq = 0; iq < nquad; iq++)
            {
                double xquade = xquad[iq + ie * nquad];
                double wquade = wquad[iq];
                //
                //  Pick a basis function which defines the equation,
                //
                for (int il = 1; il <= nl; il++)
                {
                    int ig = node[il - 1 + ie * nl];
                    int iu = indx[ig];

                    switch (iu)
                    {
                        case > 0:
                        {
                            phi(il, xquade, ref phii, ref phiix, xleft, xrite);
                            f[iu - 1] += he * wquade * ff(xquade) * phii;
                            //
                            //  Take care of boundary conditions specifying the value of U'.
                            //
                            double x;
                            switch (ig)
                            {
                                case 0:
                                    x = 0.0;
                                    f[iu - 1] -= pp(x) * ul;
                                    break;
                                default:
                                {
                                    if (ig == n)
                                    {
                                        x = 1.0;
                                        f[iu - 1] += pp(x) * ur;
                                    }

                                    break;
                                }
                            }

                            //
                            //  Pick a basis function which defines the coefficient
                            //  being computed.
                            //
                            for (int jl = 1; jl <= nl; jl++)
                            {
                                int jg = node[jl - 1 + ie * nl];
                                int ju = indx[jg];
                                phi(jl, xquade, ref phij, ref phijx, xleft, xrite);

                                double aij = he * wquade *
                                             (pp(xquade) * phiix * phijx
                                              + qq(xquade) * phii * phij);
                                switch (ju)
                                {
                                    //
                                    //  Decide where the coefficient is to be added.
                                    //
                                    case <= 0 when jg == 0:
                                        f[iu - 1] -= aij * ul;
                                        break;
                                    case <= 0:
                                    {
                                        if (jg == n)
                                        {
                                            f[iu - 1] -= aij * ur;
                                        }

                                        break;
                                    }
                                    default:
                                    {
                                        if (iu == ju)
                                        {
                                            adiag[iu - 1] += aij;
                                        }
                                        else if (ju < iu)
                                        {
                                            aleft[iu - 1] += aij;
                                        }
                                        else
                                        {
                                            arite[iu - 1] += aij;
                                        }

                                        break;
                                    }
                                }
                            }

                            break;
                        }
                    }
                }
            }
        }
    }

    public static double ff(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FF evaluates the function F in the differential equation.
        //
        //  Discussion:
        //
        //    This is the function F(X) that appears on the right hand
        //    side of the equation:
        //
        //      -d/dx ( P(x) * dU(x)/dx ) + Q(x) * U(x)  =  F(x)
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
        //    Input, double X, the evaluation point.
        //
        //    Output, double FF, the value of F(X).
        //
    {
        double value = 0;
        //
        //  Find out which problem we're working on.
        //
        int problem = get_problem();

        switch (problem)
        {
            case 1:
                value = 0.0;
                break;
            case 2:
                value = -2.0 * x;
                break;
            case 3:
                value = 0.25 * Math.PI * Math.PI * Math.Sin(0.5 * Math.PI * x);
                break;
            case 4:
                value = 0.25 * Math.PI * Math.PI * Math.Cos(0.5 * Math.PI * x);
                break;
            case 5:
            {
                double beta = get_beta();

                value = -Math.Pow(x, beta) + Math.Pow(x, beta + 2.0)
                    / ((beta + 2.0) * (beta + 1.0));
                break;
            }
            case 6:
            {
                double alpha = get_alpha();

                value = 2.0 * alpha * (x - 0.5)
                        / Math.Pow(alpha * alpha + (x - 0.5) * (x - 0.5), 2);
                break;
            }
        }

        return value;
    }

    public static void geometry(ref double[] h, int ibc, ref int[] indx, int n, int nl, int nmax,
            ref int[] node, int nquad, ref int nu, ref double[] wquad, double[] xn, ref double[] xquad)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GEOMETRY sets up some of the geometric information for the problem.  
        //
        //  Discussion:
        //
        //    Note, however, that the location of the nodes
        //    is done outside of this routine, and, in fact, before this
        //    routine is called.
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
        //    Output, double H(N)
        //    H(I) is the length of subinterval I.  This code uses
        //    equal spacing for all the subintervals.
        //
        //    Input, int IBC.
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
        //    Output, int INDX(0:N).
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
        //    Input, int N
        //    The number of subintervals into which the interval
        //    [XL,XR] is broken.
        //
        //    Input, int NL.
        //    The number of basis functions used in a single
        //    subinterval.  (NL-1) is the degree of the polynomials
        //    used.  For this code, NL is fixed at 2, meaning that
        //    piecewise linear functions are used as the basis.
        //
        //    Input, int NMAX, the maximum number of unknowns that can be handled.
        //
        //    Output, int NODE(NL,N).
        //    For each subinterval I:
        //    NODE(1,I) is the number of the left node, and
        //    NODE(2,I) is the number of the right node.
        //
        //    Input, int NQUAD
        //    The number of quadrature points used in a subinterval.
        //
        //    Output, int NU.
        //    NU is the number of unknowns in the linear system.
        //    Depending on the value of IBC, there will be N-1,
        //    N, or N+1 unknown values, which are the coefficients
        //    of basis functions.
        //
        //    Output, double WQUAD(NQUAD).
        //    WQUAD(I) is the weight associated with the I-th point
        //    of an NQUAD point Gaussian quadrature rule.
        //
        //    Input, double XN(0:N).
        //    XN(I) is the location of the I-th node.  XN(0) is XL,
        //    and XN(N) is XR.
        //
        //    Output, double XQUAD(NQUAD,NMAX), the I-th quadrature point
        //    in interval J.
        //
    {
        //
        //  Store in NODE the fact that interval I has node I-1
        //  as its left endpoint, and node I as its right endpoint.
        //
        for (int i = 0; i < n; i++)
        {
            node[0 + i * 2] = i;
            node[1 + i * 2] = i + 1;
        }

        //
        //  For every node that is associated with an unknown, we
        //  record the number of the unknown in INDX.
        //
        nu = 0;
        for (int i = 0; i <= n; i++)
        {
            switch (i)
            {
                case 0 when ibc is 1 or 3:
                    indx[i] = -1;
                    break;
                default:
                {
                    if (i == n && ibc is 2 or 3)
                    {
                        indx[i] = -1;
                    }
                    else
                    {
                        nu += 1;
                        indx[i] = nu;
                    }

                    break;
                }
            }
        }

        //
        //  We compute the width of each interval.
        //
        for (int i = 0; i < n; i++)
        {
            int igl = node[0 + i * 2];
            int igr = node[1 + i * 2];
            h[i] = xn[igr] - xn[igl];
        }

        //
        //  We compute the location of the quadrature points in each
        //  interval.
        //
        for (int i = 0; i < n; i++)
        {
            double xl = xn[node[0 + i * 2]];
            double xr = xn[node[1 + i * 2]];

            switch (nquad)
            {
                case 1:
                    xquad[0 + i * nquad] = 0.5 * (xl + xr);
                    break;
                default:
                {
                    double alfa;
                    switch (nquad)
                    {
                        case 2:
                            alfa = -0.577350;
                            xquad[0 + i * nquad] = ((1.0 - alfa) * xl
                                                    + (1.0 + alfa) * xr)
                                                   / 2.0;
                            alfa = +0.577350;
                            xquad[1 + i * nquad] = ((1.0 - alfa) * xl
                                                    + (1.0 + alfa) * xr)
                                                   / 2.0;
                            break;
                        case 3:
                            alfa = -0.774597;
                            xquad[0 + i * nquad] = ((1.0 - alfa) * xl
                                                    + (1.0 + alfa) * xr)
                                                   / 2.0;
                            xquad[1 + i * nquad] = 0.5 * (xl + xr);
                            alfa = +0.774597;
                            xquad[2 + i * nquad] = ((1.0 - alfa) * xl
                                                    + (1.0 + alfa) * xr)
                                                   / 2.0;
                            break;
                    }

                    break;
                }
            }
        }

        switch (nquad)
        {
            //
            //  Store the weights for the quadrature rule.
            //
            case 1:
                wquad[0] = 1.0;
                break;
            case 2:
                wquad[0] = 0.5;
                wquad[1] = 0.5;
                break;
            case 3:
                wquad[0] = 4.0 / 9.0;
                wquad[1] = 5.0 / 18.0;
                wquad[2] = 4.0 / 9.0;
                break;
        }
    }

    public static double get_alpha()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GET_ALPHA returns the value of ALPHA, for use by problem 6.
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
        //    Output, double GET_ALPHA, the value of ALPHA.
        //
    {
        const double value = 0.01;

        return value;
    }

    public static double get_beta()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GET_BETA returns the value of BETA, for use by problem 5.
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
        //    Output, double VALUE, the value of BETA.
        //
    {
        const double value = -0.9;

        return value;
    }

    public static int get_problem()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    GETPRB returns the value of the current problem number.
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
        //    Output, int GET_PROBLEM, the index of the problem.
        //
    {
        const int value = 6;

        return value;
    }

    public static void init(ref int ibc, int n, ref double tol, ref double ul, ref double ur, ref double xl,
            ref double[] xn, ref double xr)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INIT initializes some parameters that define the problem.
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
        //    Output, int *IBC.
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
        //    Input, int N
        //    The number of subintervals into which the interval
        //    [XL,XR] is broken.
        //
        //    Output, double *TOL.
        //    A tolerance that is used to determine whether the estimated
        //    error in an interval is so large that it should be subdivided
        //    and the problem solved again.
        //
        //    Output, double *UL.
        //    If IBC is 1 or 3, UL is the value that U is required
        //    to have at X = XL.
        //    If IBC is 2 or 4, UL is the value that U' is required
        //    to have at X = XL.
        //
        //    Output, double *UR.
        //    If IBC is 2 or 3, UR is the value that U is required
        //    to have at X = XR.
        //    If IBC is 1 or 4, UR is the value that U' is required
        //    to have at X = XR.
        //
        //    Output, double *XL.
        //    XL is the left endpoint of the interval over which the
        //    differential equation is being solved.
        //
        //    Output, double XN(0:N).
        //    XN(I) is the location of the I-th node.  XN(0) is XL,
        //    and XN(N) is XR.
        //
        //    Output, double *XR.
        //    XR is the right endpoint of the interval over which the
        //    differential equation is being solved.
        //
    {
        tol = 0.01;
        //
        //  Find out which problem we're working on.
        //
        int problem = get_problem();
        switch (problem)
        {
            //
            //  Set the boundary conditions for the problem, and
            //  print out its title.
            //
            case 1:
                ibc = 3;
                ul = 0.0;
                ur = 1.0;
                xl = 0.0;
                xr = 1.0;
                Console.WriteLine("");
                Console.WriteLine("Exact solution is U = X");
                break;
            case 2:
                ibc = 3;
                ul = 0.0;
                ur = 1.0;
                xl = 0.0;
                xr = 1.0;
                Console.WriteLine("");
                Console.WriteLine("Exact solution is U = X*X");
                break;
            case 3:
                ibc = 3;
                ul = 0.0;
                ur = 1.0;
                xl = 0.0;
                xr = 1.0;
                Console.WriteLine("");
                Console.WriteLine("Exact solution is U = SIN(PI*X/2)");
                break;
            case 4:
                ibc = 3;
                ul = 1.0;
                ur = 0.0;
                xl = 0.0;
                xr = 1.0;
                Console.WriteLine("");
                Console.WriteLine("Exact solution is U = COS(PI*X/2)");
                break;
            case 5:
            {
                ibc = 3;
                double beta = get_beta();
                ul = 0.0;
                ur = 1.0 / ((beta + 2.0) * (beta + 1.0));
                xl = 0.0;
                xr = 1.0;
                Console.WriteLine("");
                Console.WriteLine("Rheinboldt problem");
                break;
            }
            case 6:
                ibc = 3;
                xl = 0.0;
                xr = 1.0;
                ul = uexact(xl);
                ur = uexact(xr);
                Console.WriteLine("");
                Console.WriteLine("Arctangent problem");
                break;
        }

        //
        //  The nodes are defined here, and not in the geometry routine.
        //  This is because each new iteration chooses the location
        //  of the new nodes in a special way.
        //
        for (int i = 0; i <= n; i++)
        {
            xn[i] = ((n - i) * xl
                     + i * xr)
                    / n;
        }

        Console.WriteLine("The equation is to be solved for ");
        Console.WriteLine("X greater than " + xl + "");
        Console.WriteLine(" and less than " + xr + "");
        Console.WriteLine("");
        Console.WriteLine("The boundary conditions are:");
        Console.WriteLine("");

        switch (ibc)
        {
            case 1:
            case 3:
                Console.WriteLine("  At X = XL, U = " + ul + "");
                break;
            default:
                Console.WriteLine("  At X = XL, U' = " + ul + "");
                break;
        }

        switch (ibc)
        {
            case 2:
            case 3:
                Console.WriteLine("  At X = XR, U= " + ur + "");
                break;
            default:
                Console.WriteLine("  At X = XR, U' = " + ur + "");
                break;
        }
    }

    public static void output(double[] f, int ibc, int[] indx, int n, int nu, double ul,
            double ur, double[] xn)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    OUTPUT prints out the computed solution.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 November 2006
        //
        //  Author:
        //
        //    C++ version by John Burkardt
        //
        //  Parameters:
        //
        //    Input, double F(NU).
        //    ASSEMBLE stores into F the right hand side of the linear
        //    equations.
        //    SOLVE replaces those values of F by the solution of the
        //    linear equations.
        //
        //    Input, int IBC.
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
        //    Input, int INDX(0:N).
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
        //    Input, int N
        //    The number of subintervals into which the interval
        //    [XL,XR] is broken.
        //
        //    Input, int NU.
        //    NU is the number of unknowns in the linear system.
        //    Depending on the value of IBC, there will be N-1,
        //    N, or N+1 unknown values, which are the coefficients
        //    of basis functions.
        //
        //    Input, double UL.
        //    If IBC is 1 or 3, UL is the value that U is required
        //    to have at X = XL.
        //    If IBC is 2 or 4, UL is the value that U' is required
        //    to have at X = XL.
        //
        //    Input, double UR.
        //    If IBC is 2 or 3, UR is the value that U is required
        //    to have at X = XR.
        //    If IBC is 1 or 4, UR is the value that U' is required
        //    to have at X = XR.
        //
        //    Input, double XN(0:N).
        //    XN(I) is the location of the I-th node.  XN(0) is XL,
        //    and XN(N) is XR.
        //
    {
        Console.WriteLine("");
        Console.WriteLine("Node    X(I)        U(X(I))        Uexact        Error");
        Console.WriteLine("");

        for (int i = 0; i <= n; i++)
        {
            double u;
            switch (i)
            {
                case 0 when ibc is 1 or 3:
                    u = ul;
                    break;
                case 0:
                    u = f[indx[i] - 1];
                    break;
                default:
                {
                    if (i == n)
                    {
                        switch (ibc)
                        {
                            case 2:
                            case 3:
                                u = ur;
                                break;
                            default:
                                u = f[indx[i] - 1];
                                break;
                        }
                    }
                    else
                    {
                        u = f[indx[i] - 1];
                    }

                    break;
                }
            }

            double uex = uexact(xn[i]);
            double error = u - uex;

            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(4)
                                   + "  " + xn[i].ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + u.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + uex.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + error.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    public static void phi(int il, double x, ref double phii, ref double phiix, double xleft,
            double xrite)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PHI evaluates a linear basis function and its derivative.
        //
        //  Discussion:
        //
        //    The functions are evaluated at a point X in an interval.  In any
        //    interval, there are just two basis functions.  The first
        //    basis function is a line which is 1 at the left endpoint
        //    and 0 at the right.  The second basis function is 0 at
        //    the left endpoint and 1 at the right.
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
        //    Input, int IL, the local index of the basis function.
        //
        //    Input, double X, the evaluation point.
        //
        //    Output, double *PHII, *PHIIX, the value of the basis function
        //    and its derivative.
        //
        //    Input, double XLEFT, XRITE, the endpoints of the interval.
        //
    {
        if (xleft <= x && x <= xrite)
        {
            switch (il)
            {
                case 1:
                    phii = (xrite - x) / (xrite - xleft);
                    phiix = -1.0 / (xrite - xleft);
                    break;
                default:
                    phii = (x - xleft) / (xrite - xleft);
                    phiix = 1.0 / (xrite - xleft);
                    break;
            }
        }
        //
        //  If X is outside of the interval, then the basis function
        //  is always zero.
        //
        else
        {
            phii = 0.0;
            phiix = 0.0;
        }
    }

    public static double pp(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PP evaluates the function P in the differential equation.
        //
        //  Discussion:
        //
        //    The function P(X) occurs in the differential equation:
        //
        //      -d/dx ( P(x) * dU(x)/dx ) + Q(x) * U(x)  =  F(x)
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
        //    Input, double X, the evaluation point.
        //
        //    Output, double PP, the value of P(X).
        //
    {
        double value = 0;
        //
        //  Find out which problem we're working on.
        //
        int problem = get_problem();

        switch (problem)
        {
            case 1:
            case 2:
            case 3:
            case 4:
            case 5:
            case 6:
                value = 1.0;
                break;
        }

        return value;
    }

    public static void prsys(double[] adiag, double[] aleft, double[] arite, double[] f,
            int nu)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PRSYS prints out the tridiagonal linear system.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    06 November 2006
        //
        //  Author:
        //
        //    C++ version by John Burkardt
        //
        //  Parameters:
        //
        //    Input, double ADIAG(NU).
        //    ADIAG(I) is the "diagonal" coefficient of the I-th
        //    equation in the linear system.  That is, ADIAG(I) is
        //    the coefficient of the I-th unknown in the I-th equation.
        //
        //    Input, double ALEFT(NU).
        //    ALEFT(I) is the "left hand" coefficient of the I-th
        //    equation in the linear system.  That is, ALEFT(I) is the
        //    coefficient of the (I-1)-th unknown in the I-th equation.
        //    There is no value in ALEFT(1), since the first equation
        //    does not refer to a "0-th" unknown.
        //
        //    Input, double ARITE(NU).
        //    ARITE(I) is the "right hand" coefficient of the I-th
        //    equation in the linear system.  ARITE(I) is the coefficient
        //    of the (I+1)-th unknown in the I-th equation.  There is
        //    no value in ARITE(NU) because the NU-th equation does not
        //    refer to an "NU+1"-th unknown.
        //
        //    Input, double F(NU).
        //    ASSEMBLE stores into F the right hand side of the linear
        //    equations.
        //    SOLVE replaces those values of F by the solution of the
        //    linear equations.
        //
        //    int NU.
        //    NU is the number of unknowns in the linear system.
        //    Depending on the value of IBC, there will be N-1,
        //    N, or N+1 unknown values, which are the coefficients
        //    of basis functions.
        //
    {
        Console.WriteLine("");
        Console.WriteLine("Printout of tridiagonal linear system:");
        Console.WriteLine("");
        Console.WriteLine("Equation  A-Left  A-Diag  A-Rite  RHS");
        Console.WriteLine("");

        for (int i = 0; i < nu; i++)
        {
            string cout = (i + 1).ToString(CultureInfo.InvariantCulture).PadLeft(4);
            cout += i switch
            {
                0 => "              ",
                _ => "  " + aleft[i].ToString(CultureInfo.InvariantCulture).PadLeft(12)
            };

            cout += "  " + adiag[i].ToString(CultureInfo.InvariantCulture).PadLeft(12);
            if (i < nu - 1)
            {
                cout += "  " + arite[i].ToString(CultureInfo.InvariantCulture).PadLeft(12);
            }
            else
            {
                cout += "              ";
            }

            cout += "  " + f[i].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "";

            Console.WriteLine(cout);
        }
    }

    public static double qq(double x)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    QQ evaluates the function Q in the differential equation.
        //
        //  Discussion:
        //
        //    The function Q(X) occurs in the differential equation:
        //
        //      -d/dx ( P(x) * dU(x)/dx ) + Q(x) * U(x)  =  F(x)
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
        //    Input, double X, the evaluation point.
        //
        //    Output, double QQ, the value of Q(X).
        //
    {
        double value = 0;
        //
        //  Find out which problem we're working on.
        //
        int problem = get_problem();

        switch (problem)
        {
            case 1:
            case 2:
            case 3:
            case 4:
                value = 0.0;
                break;
            case 5:
                value = 1.0;
                break;
            case 6:
                value = 0.0;
                break;
        }

        return value;
    }

    public static void solve(ref double[] adiag, ref double[] aleft, ref double[] arite, ref double[] f,
            int nu)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SOLVE solves a tridiagonal matrix system of the form A*x = b.
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
        //    Input/output, double ADIAG(NU), ALEFT(NU), ARITE(NU).
        //    On input, ADIAG, ALEFT, and ARITE contain the diagonal,
        //    left and right entries of the equations.
        //    On output, ADIAG and ARITE have been changed in order
        //    to compute the solution.
        //    Note that for the first equation, there is no ALEFT
        //    coefficient, and for the last, there is no ARITE.
        //    So there is no need to store a value in ALEFT(1), nor
        //    in ARITE(NU).
        //
        //    Input/output, double F(NU).
        //    On input, F contains the right hand side of the linear
        //    system to be solve
        //    On output, F contains the solution of the linear system.
        //
        //    Input, int NU, the number of equations to be solved.
        //
    {
        switch (nu)
        {
            //
            //  Handle the special case of a single equation.
            //
            case 1:
                f[0] /= adiag[0];
                break;
            //
            default:
            {
                arite[0] /= adiag[0];
                for (int i = 2; i <= nu - 1; i++)
                {
                    adiag[i - 1] -= aleft[i - 1] * arite[i - 2];
                    arite[i - 1] /= adiag[i - 1];
                }

                adiag[nu - 1] -= aleft[nu - 1] * arite[nu - 2];

                f[0] /= adiag[0];
                for (int i = 2; i <= nu; i++)
                {
                    f[i - 1] = (f[i - 1] - aleft[i - 1] * f[i - 2]) / adiag[i - 1];
                }

                for (int i = nu - 1; 1 <= i; i--)
                {
                    f[i - 1] -= arite[i - 1] * f[i];
                }

                break;
            }
        }
    }

    public static void solvex(double[] adiag, double[] aleft, double[] arite, ref double[] f,
            ref double[] h, int ibc, int[] indx, int kount, int n, int nl, int nmax,
            int[] node, int nquad, ref int nu, double ul, double ur, double[] wquad,
            double[] xn, double[] xquad)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SOLVEX discretizes and solves a differential equation given the nodes.
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
        //    Workspace, double ADIAG(NU).
        //    ADIAG(I) is the "diagonal" coefficient of the I-th
        //    equation in the linear system.  That is, ADIAG(I) is
        //    the coefficient of the I-th unknown in the I-th equation.
        //
        //    Workspace, double ALEFT(NU).
        //    ALEFT(I) is the "left hand" coefficient of the I-th
        //    equation in the linear system.  That is, ALEFT(I) is the
        //    coefficient of the (I-1)-th unknown in the I-th equation.
        //    There is no value in ALEFT(1), since the first equation
        //    does not refer to a "0-th" unknown.
        //
        //    Workspace, double ARITE(NU).
        //    ARITE(I) is the "right hand" coefficient of the I-th
        //    equation in the linear system.  ARITE(I) is the coefficient
        //    of the (I+1)-th unknown in the I-th equation.  There is
        //    no value in ARITE(NU) because the NU-th equation does not
        //    refer to an "NU+1"-th unknown.
        //
        //    Output, double F(NU).
        //    ASSEMBLE stores into F the right hand side of the linear
        //    equations.
        //    SOLVE replaces those values of F by the solution of the
        //    linear equations.
        //
        //    Output, double H(N)
        //    H(I) is the length of subinterval I.  This code uses
        //    equal spacing for all the subintervals.
        //
        //    Input, int IBC.
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
        //    Workspace, int INDX(0:N).
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
        //    Input, int KOUNT, the number of adaptive steps that have been taken.
        //
        //    Input, int N
        //    The number of subintervals into which the interval
        //    [XL,XR] is broken.
        //
        //    Input, int NL.
        //    The number of basis functions used in a single
        //    subinterval.  (NL-1) is the degree of the polynomials
        //    used.  For this code, NL is fixed at 2, meaning that
        //    piecewise linear functions are used as the basis.
        //
        //    Input, int NMAX, the maximum number of unknowns that can be handled.
        //
        //    Workspace, int NODE(NL,N).
        //    For each subinterval I:
        //    NODE(1,I) is the number of the left node, and
        //    NODE(2,I) is the number of the right node.
        //
        //    Workspace, int NQUAD
        //    The number of quadrature points used in a subinterval.
        //
        //    Output, int *NU.
        //    NU is the number of unknowns in the linear system.
        //    Depending on the value of IBC, there will be N-1,
        //    N, or N+1 unknown values, which are the coefficients
        //    of basis functions.
        //
        //    Input, double UL.
        //    If IBC is 1 or 3, UL is the value that U is required
        //    to have at X = XL.
        //    If IBC is 2 or 4, UL is the value that U' is required
        //    to have at X = XL.
        //
        //    Input, double UR.
        //    If IBC is 2 or 3, UR is the value that U is required
        //    to have at X = XR.
        //    If IBC is 1 or 4, UR is the value that U' is required
        //    to have at X = XR.
        //
        //    Workspace, double WQUAD(NQUAD).
        //    WQUAD(I) is the weight associated with the I-th point
        //    of an NQUAD point Gaussian quadrature rule.
        //
        //    Input, double XN(0:N).
        //    XN(I) is the location of the I-th node.  XN(0) is XL,
        //    and XN(N) is XR.
        //
        //    Workspace, double XQUAD(NQUAD,NMAX), the I-th quadrature point
        //    in interval J.
        //
    {
        //
        //  Given a set of N nodes (where N increases on each iteration),
        //  compute the other geometric information.
        //
        geometry(ref h, ibc, ref indx, n, nl, nmax, ref node, nquad, ref nu, ref wquad, xn, ref xquad);
        //
        //  Assemble the linear system.
        //
        assemble(ref adiag, ref aleft, ref arite, ref f, h, n, indx, node, nu, nl,
            nquad, nmax, ul, ur, wquad, xn, xquad);
        switch (kount)
        {
            //
            //  Print out the linear system, just once.
            //
            case 1:
                prsys(adiag, aleft, arite, f, nu);
                break;
        }

        //
        //  Solve the linear system.
        //
        solve(ref adiag, ref aleft, ref arite, ref f, nu);
        //
        //  Print out the solution.
        //
        Console.WriteLine("");
        Console.WriteLine("Basic solution");

        output(f, ibc, indx, n, nu, ul, ur, xn);
    }

    public static void solvey(ref double[] eta, double[] f, double[] h, int n, int nu, double ul,
            double ur, double[] xn)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SOLVEY computes error estimators for a finite element solution.
        //
        //  Discussion:
        //
        //    SOLVEY accepts information about the solution of a finite element
        //    problem on a grid of nodes with coordinates XN.  It then starts
        //    at node 0, and for each node, computes two "error estimators",
        //    one for the left, and one for the right interval associated with the
        //    node.  These estimators are found by solving a finite element problem
        //    over the two intervals, using the known values of the original
        //    solution as boundary data, and using a mesh that is "slightly"
        //    refined over the original one.
        //
        //    Note that the computations at the 0-th and N-th nodes only involve
        //    a single interval.
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
        //    Output, double ETA(N).
        //    ETA(I) is the error estimate for interval I.  It is computed
        //    as the sum of two quantities, one associated with the left
        //    and one with the right node of the interval.
        //
        //    Input, double F(NU).
        //    ASSEMBLE stores into F the right hand side of the linear
        //    equations.
        //    SOLVE replaces those values of F by the solution of the
        //    linear equations.
        //
        //    Input, double H(N)
        //    H(I) is the length of subinterval I.  This code uses
        //    equal spacing for all the subintervals.
        //
        //    Input, int N
        //    The number of subintervals into which the interval
        //    [XL,XR] is broken.
        //
        //    Input, int NU.
        //    NU is the number of unknowns in the linear system.
        //    Depending on the value of IBC, there will be N-1,
        //    N, or N+1 unknown values, which are the coefficients
        //    of basis functions.
        //
        //    Input, double UL.
        //    If IBC is 1 or 3, UL is the value that U is required
        //    to have at X = XL.
        //    If IBC is 2 or 4, UL is the value that U' is required
        //    to have at X = XL.
        //
        //    Input, double UR.
        //    If IBC is 2 or 3, UR is the value that U is required
        //    to have at X = XR.
        //    If IBC is 1 or 4, UR is the value that U' is required
        //    to have at X = XR.
        //
        //    Input, double XN(0:N).
        //    XN(I) is the location of the I-th node.  XN(0) is XL,
        //    and XN(N) is XR.
        //
    {
        const int NL = 2;
        const int NY = 2;
        const int NQUAD = 2;
        const int NMAY = 2 * NY;

        double[] adiag = new double[NMAY];
        double[] aleft = new double[NMAY];
        double[] arite = new double[NMAY];
        double[] fy = new double[NMAY];
        double[] hy = new double[NMAY];
        int[] indy = new int[NMAY + 1];
        int[] nodey = new int[NL * NMAY];
        int nuy = 0;
        double[] wquad = new double[NQUAD];
        double[] xquady = new double[NQUAD * NMAY];
        double[] yn = new double[NMAY + 1];
        //
        //  Initialize the error estimators to zero.
        //
        for (int j = 0; j < n; j++)
        {
            eta[j] = 0.0;
        }

        //
        //  Set the boundary conditions for each subproblem to be
        //  known values of U at the left and right.
        //
        //
        //  For each node, subdivide its left and right hand intervals
        //  into NY subintervals.
        //
        //  Set up and solve the differential equation again on this
        //  smaller region.
        //
        //  The 0-th and N-th nodes are special cases.
        //
        const int ibcy = 3;

        for (int j = 0; j <= n; j++)
        {
            int jhi;
            int jlo;
            int jmid;
            int m;
            switch (j)
            {
                case 0:
                    m = NY;
                    jlo = j;
                    jmid = j + 1;
                    jhi = j + 1;
                    break;
                default:
                {
                    if (j == n)
                    {
                        m = NY;
                        jlo = j - 1;
                        jmid = j;
                        jhi = j;
                    }
                    else
                    {
                        m = 2 * NY;
                        jlo = j - 1;
                        jmid = j;
                        jhi = j + 1;
                    }

                    break;
                }
            }

            //
            //  Set the location of the nodes in the subintervals.
            //
            double yl = xn[jlo];
            double ym = xn[jmid];
            double yr = xn[jhi];

            for (int i = 0; i <= NY; i++)
            {
                yn[i] = ((NY - i) * yl
                         + i * ym)
                        / NY;
            }

            for (int i = NY + 1; i <= m; i++)
            {
                yn[i] = ((m - i) * ym
                         + (i - NY) * yr)
                        / (m - NY);
            }

            //
            //  Set up the geometry of the sub-problem.
            //
            geometry(ref hy, ibcy, ref indy, m, NL, NMAY, ref nodey, NQUAD, ref nuy,
                ref wquad, yn, ref xquady);
            //
            //  Set the boundary values for the sub-problem.
            //
            double uly = j switch
            {
                <= 1 => ul,
                _ => f[j - 2]
            };

            double ury = n - 1 <= j ? ur : f[j];

            //
            //  Assemble the matrix for the sub-problem.
            //
            assemble(ref adiag, ref aleft, ref arite, ref fy, hy, m, indy, nodey, nuy, NL,
                NQUAD, NMAY, uly, ury, wquad, yn, xquady);
            //
            //  Solve the system.
            //
            solve(ref adiag, ref aleft, ref arite, ref fy, nuy);
            //
            //  Compute the weighted sum of the squares of the differences
            //  of the original computed slope and the refined computed slopes.
            //
            //  Calculation for left interval.
            //
            double total;
            double uleft;
            double ulval;
            double uprime;
            double vrval;
            double uval;
            double vprime;
            double y;
            double vval;
            double urite;
            double urval;
            double vlval;
            switch (j)
            {
                case >= 1:
                {
                    switch (j)
                    {
                        case <= 1:
                            uleft = ul;
                            urite = f[0];
                            break;
                        default:
                        {
                            if (j == n)
                            {
                                uleft = f[j - 2];
                                urite = ur;
                            }
                            else
                            {
                                uleft = f[j - 2];
                                urite = f[j - 1];
                            }

                            break;
                        }
                    }

                    uprime = (urite - uleft) / h[j - 1];

                    total = 0.0;

                    for (int i = 1; i <= NY; i++)
                    {
                        yl = yn[i - 1];
                        yr = yn[i];

                        switch (i)
                        {
                            case 1:
                                vlval = uly;
                                vrval = fy[i - 1];
                                break;
                            default:
                            {
                                if (i == m)
                                {
                                    vlval = fy[i - 2];
                                    vrval = ury;
                                }
                                else
                                {
                                    vlval = fy[i - 2];
                                    vrval = fy[i - 1];
                                }

                                break;
                            }
                        }

                        vprime = (vrval - vlval) / hy[i - 1];

                        ulval = ((NY - i + 1) * uleft
                                 + (i - 1) * urite)
                                / NY;

                        urval = ((NY - i) * uleft
                                 + i * urite)
                                / NY;
                        //
                        //  Compute the integral of
                        //
                        //    p(x)*(u'(x)-v'(x))**2 + q(x)*(u(x)-v(x))**2
                        //
                        for (int k = 0; k < NQUAD; k++)
                        {
                            y = xquady[k + (i - 1) * NQUAD];

                            uval = ((yl - y) * urval
                                    + (y - yr) * ulval)
                                   / (yl - yr);

                            vval = ((yl - y) * vrval
                                    + (y - yr) * vlval)
                                   / (yl - yr);

                            total += 0.5 * wquad[k] * hy[i - 1] *
                                     (pp(y) * Math.Pow(uprime - vprime, 2)
                                      + qq(y) * Math.Pow(uval - vval, 2));
                        }
                    }

                    eta[j - 1] += 0.5 * Math.Sqrt(total);
                    break;
                }
            }

            //
            //  Calculation for right interval.
            //
            if (j > n - 1)
            {
                continue;
            }

            {
                switch (j)
                {
                    case 0:
                        uleft = ul;
                        urite = f[j];
                        break;
                    default:
                    {
                        if (n - 1 <= j)
                        {
                            uleft = f[j - 1];
                            urite = ur;
                        }
                        else
                        {
                            uleft = f[j - 1];
                            urite = f[j];
                        }

                        break;
                    }
                }

                uprime = (urite - uleft) / h[j];

                total = 0.0;

                for (int i = m + 1 - NY; i <= m; i++)
                {
                    yl = yn[i - 1];
                    yr = yn[i];

                    switch (i)
                    {
                        case 1:
                            vlval = uly;
                            vrval = fy[i - 1];
                            break;
                        default:
                        {
                            if (i == m)
                            {
                                vlval = fy[i - 2];
                                vrval = ury;
                            }
                            else
                            {
                                vlval = fy[i - 2];
                                vrval = fy[i - 1];
                            }

                            break;
                        }
                    }

                    vprime = (vrval - vlval) / hy[i - 1];

                    ulval = ((m - i + 1) * uleft
                             + (NY - m + i - 1) * urite)
                            / NY;

                    urval = ((m - i) * uleft
                             + (NY - m + i) * urite)
                            / NY;
                    //
                    //  Compute the integral of
                    //
                    //    p(x)*(u'(x)-v'(x))**2 + q(x)*(u(x)-v(x))**2
                    //
                    for (int k = 0; k < NQUAD; k++)
                    {
                        y = xquady[k + (i - 1) * NQUAD];

                        uval = ((yl - y) * urval
                                + (y - yr) * ulval)
                               / (yl - yr);

                        vval = ((yl - y) * vrval
                                + (y - yr) * vlval)
                               / (yl - yr);

                        total += 0.5 * wquad[k] * hy[i - 1] *
                                 (pp(y) * Math.Pow(uprime - vprime, 2)
                                  + qq(y) * Math.Pow(uval - vval, 2));
                    }
                }

                eta[j] += 0.5 * Math.Sqrt(total);
            }
        }

        //
        //  Print out the error estimators.
        //
        Console.WriteLine("");
        Console.WriteLine("ETA");
        Console.WriteLine("");
        for (int j = 0; j < n; j++)
        {
            Console.WriteLine(eta[j].ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    public static int subdiv(double[] eta, int kount, ref int n, int nmax, double tol,
            ref double[] xn)
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    SUBDIV decides which intervals should be subdivided.
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
        //    Input, double ETA(N).
        //    ETA(I) is the error estimate for interval I.  It is computed
        //    as the sum of two quantities, one associated with the left
        //    and one with the right node of the interval.
        //
        //    Input, int KOUNT, the number of adaptive steps that have been taken.
        //
        //    Input/output, int N
        //    The number of subintervals into which the interval
        //    [XL,XR] is broken.  
        //
        //    Input, int NMAX, the maximum number of unknowns that can be handled.
        //
        //    Input, double TOL.
        //    A tolerance that is used to determine whether the estimated
        //    error in an interval is so large that it should be subdivided
        //    and the problem solved again.
        //
        //    Input/output, double XN(0:N).
        //    XN(I) is the location of the I-th node.  XN(0) is XL,
        //    and XN(N) is XR.
        //
        //    Output, int SUBDIV, reports status of subdivision.
        //    0, a new subdivision was carried out.
        //    1, no more subdivisions are needed.
        //    -1, no more subdivisions can be carried out.
        //
        //  Local Parameters:
        //
        //    Local, int JADD(N).
        //    JADD(I) is 1 if the error estimates show that interval I
        //    should be subdivided.
        //
    {
        int status = 0;
        //
        //  Add up the ETA's, and get their average.
        //
        double ave = 0.0;
        for (int j = 0; j < n; j++)
        {
            ave += eta[j];
        }

        ave /= n;
        //
        //  Look for intervals whose ETA value is relatively large,
        //  and note in JADD that these intervals should be subdivided.
        //
        int k = 0;
        double temp = Math.Max(1.2 * ave + 0.00001, tol * tol / n);

        Console.WriteLine("");
        Console.WriteLine("Tolerance = " + temp + "");
        Console.WriteLine("");

        int[] jadd = new int[n];

        for (int j = 0; j < n; j++)
        {
            if (temp < eta[j])
            {
                k += 1;
                jadd[j] = 1;
                Console.WriteLine("Subdivide interval " + (j + 1) + "");
            }
            else
            {
                jadd[j] = 0;
            }
        }

        switch (k)
        {
            //
            //  If no subdivisions needed, we're done.
            //
            case <= 0:
                Console.WriteLine("Success on step " + kount + "");
                status = 1;
                return status;
        }

        //
        //  See if we're about to go over our limit.
        //
        if (nmax < n + k)
        {
            Console.WriteLine("");
            Console.WriteLine("The iterations did not reach their goal.");
            Console.WriteLine("The next value of N is " + (n + k) + "");
            Console.WriteLine("which exceeds NMAX = " + nmax + "");
            status = -1;
            return status;
        }

        //
        //  Insert new nodes where needed.
        //
        double[] xt = new double[nmax + 1];

        k = 0;
        xt[0] = xn[0];
        for (int j = 0; j < n; j++)
        {
            switch (jadd[j])
            {
                case > 0:
                    xt[j + 1 + k] = 0.5 * (xn[j + 1] + xn[j]);
                    k += 1;
                    break;
            }

            xt[j + 1 + k] = xn[j + 1];
        }

        //
        //  Update the value of N, and copy the new nodes into XN.
        //
        n += k;

        for (int j = 0; j <= n; j++)
        {
            xn[j] = xt[j];
        }

        return status;
    }

    public static double uexact(double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    UEXACT returns the value of the exact solution at any point X.
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
        //    Input, double X, the evaluation point.
        //
        //    Output, double UEXACT, the value of the exact solution at X.
        //
    {
        double value = 0;
        //
        //  Find out which problem we're working on.
        //
        int problem = get_problem();

        switch (problem)
        {
            case 1:
                value = x;
                break;
            case 2:
                value = x * x;
                break;
            case 3:
                value = Math.Sin(Math.PI * x / 2.0);
                break;
            case 4:
                value = Math.Cos(Math.PI * x / 2.0);
                break;
            case 5:
            {
                double beta = get_beta();

                value = Math.Pow(x, beta + 2.0) / ((beta + 2.0) * (beta + 1.0));
                break;
            }
            case 6:
            {
                double alpha = get_alpha();
                value = Math.Atan((x - 0.5) / alpha);
                break;
            }
            default:
                value = 0.0;
                break;
        }

        return value;
    }
}