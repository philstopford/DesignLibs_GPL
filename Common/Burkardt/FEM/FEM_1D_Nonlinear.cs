using System;
using System.Globalization;

namespace Burkardt.FEM;

public static class FEM_1D_Nonlinear
{
    public static void assemble_newton(ref double[] adiag, ref double[] aleft, ref double[] arite,
        ref double[] f, double[] fold, double[] h, int[] indx, int n, int nl,
        int[] node, int nquad, int nu, int problem, double ul, double ur,
        double[] xn, double[] xquad )
//****************************************************************************80
//
//  Purpose:
//
//    ASSEMBLE_NEWTON assembles the Newton linear system.
//
//  Discussion:
//
//    The linear system being solved here is for the Newton correction
//    to an approximate solution of a nonlinear system.
//
//    Thus, we suppose that we have a nonlinear function F(X),
//    and an approximate solution X0.  If we can suppose there is an
//    exact solution X* that is "nearby", and in fact close enough
//    that Taylor's theorem gives us a useful estimate, then we
//    may write:
//
//      F(X*) = F(X0) + F'(X0) * ( X* - X0 ) + Order ( X* - X0 )^2 
//
//    and by rearranging, we get the Newton system (which is only
//    approximately correct):
//
//      F'(X0) * ( X* - X0 ) = - F(X0)
//
//    We solve this system and add the solution to X0 to get a
//    new approximate solution that, we hope, is much improved.
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
//    Output, double F(NU), the right hand side of the linear
//    equations.
//
//    Input, double FOLD(NU), the solution value 
//    from the previous iteration,
//
//    Input, double H(N), the length of the subintervals.  
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
//    Input, int N.
//    The number of subintervals into which the interval
//    [XL,XR] is broken.
//
//    Input, int NL.
//    The number of basis functions used in a single
//    subinterval.  (NL-1) is the degree of the polynomials
//    used.  For this code, NL is fixed at 2, meaning that
//    piecewise linear functions are used as the basis.
//
//    Input, int NODE(NL,N).
//    For each subinterval I:
//    NODE(1,I) is the number of the left node, and
//    NODE(2,I) is the number of the right node.
//
//    Input, int NQUAD.
//    The number of quadrature points used in a subinterval.
//    This code uses NQUAD = 1.
//
//    Input, int NU.
//    NU is the number of unknowns in the linear system.
//    Depending on the value of IBC, there will be N-1,
//    N, or N+1 unknown values, which are the coefficients
//    of basis functions.
//
//    Input, int PROBLEM, indicates which problem to be solved.
//    * 1, u = x, p=1, q=0, f=x, u(0)=0, u'(1)=1.
//    * 2, u = 2*(1-cos(0.5*pi*x))/pi, p=1, q=0,
//    f = -0.5*pi*cos(0.5*pi*x) + 2*sin(0.5*pi*x)*(1-cos(0.5*pi*x)/pi
//    u(0) = 0, u'(1)=1.
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
//    Input, double XQUAD(N)
//    XQUAD(I) is the location of the single quadrature point
//    in interval I.
//
    {
        double phii = 0;
        double phiix = 0;
        double phij = 0;
        double phijx = 0;

        for (int i = 0; i < nu; i++)
        {
            f[i] = 0.0;
            adiag[i] = 0.0;
            aleft[i] = 0.0;
            arite[i] = 0.0;
        }

//
//  For element IE...
//
        for (int ie = 0; ie < n; ie++)
        {
            double he = h[ie];
            double xleft = xn[node[0 + ie * 2]];
            double xrite = xn[node[1 + ie * 2]];
//
//  For quadrature point IQ...
//
            for (int iq = 0; iq < nquad; iq++)
            {
                double xqe = xquad[ie];
//
//  Compute value of U for previous solution.
//
                double total = 0.0;

                int ig;
                int iu;
                for (int il = 1; il <= nl; il++)
                {
                    ig = node[il - 1 + ie * 2];
                    iu = indx[ig] - 1;

                    total += iu switch
                    {
                        < 0 when il == 1 => ul,
                        < 0 => ur,
                        _ => fold[iu]
                    };
                }

                double uold = total / nl;
//
//  Compute value of U' for previous solution.
//
                int jl = node[0 + ie * 2];
                int jr = node[1 + ie * 2];
                int iul = indx[jl] - 1;
                int iur = indx[jr] - 1;

                double uoldx;
                switch (iul)
                {
                    case < 0:
                        uoldx = (fold[iur] - ul) / he;
                        break;
                    default:
                    {
                        uoldx = iur switch
                        {
                            < 0 => (ur - fold[iul]) / he,
                            _ => (fold[iur] - fold[iul]) / he
                        };

                        break;
                    }
                }

//
//  For basis function IL...
//
                for (int il = 1; il <= nl; il++)
                {
                    ig = node[il - 1 + ie * 2];
                    iu = indx[ig] - 1;

                    switch (iu)
                    {
                        case >= 0:
                        {
                            phi(il, xqe, ref phii, ref phiix, xleft, xrite);

                            f[iu] += he * phii * (ff(xqe, problem) + uold * uoldx);
//
//  Handle boundary conditions that prescribe the value of U'.
//
                            double x;
                            switch (ig)
                            {
                                case 0:
                                    x = 0.0;
                                    f[iu] -= pp(x, problem) * ul;
                                    break;
                                default:
                                {
                                    if (ig == n)
                                    {
                                        x = 1.0;
                                        f[iu] += pp(x, problem) * ur;
                                    }

                                    break;
                                }
                            }

//
//  For basis function JL...
//
                            for (jl = 1; jl <= nl; jl++)
                            {
                                int jg = node[jl - 1 + ie * 2];
                                int ju = indx[jg] - 1;

                                phi(jl, xqe, ref phij, ref phijx, xleft, xrite);

                                double aij = he * (pp(xqe, problem) * phiix * phijx
                                                   + qq(xqe, problem) * phii * phij
                                                   + uold * phii * phijx
                                                   + uoldx * phij * phii);

                                switch (ju)
                                {
                                    case < 0 when jg == 0:
                                        f[iu] -= aij * ul;
                                        break;
                                    case < 0:
                                    {
                                        if (jg == n)
                                        {
                                            f[iu] -= aij * ur;
                                        }

                                        break;
                                    }
                                    default:
                                    {
                                        if (iu == ju)
                                        {
                                            adiag[iu] += aij;
                                        }
                                        else if (ju < iu)
                                        {
                                            aleft[iu] += aij;
                                        }
                                        else
                                        {
                                            arite[iu] += aij;
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

    public static void assemble_picard(ref double[] adiag, ref double[] aleft, ref double[] arite,
        ref double[] f, double[] fold, double[] h, int[] indx, int n, int nl,
        int[] node, int nquad, int nu, int problem, double ul, double ur,
        double[] xn, double[] xquad )
//****************************************************************************80
//
//  Purpose:
//
//    ASSEMBLE_PICARD assembles the Picard linear system.
//
//  Discussion:
//
//    The equation we are trying to solve has the form:
//
//      -d/dx ( p(x) du/dx ) + q(x) * u + u * u' = f(x)
//
//    For the Picard iteration, we need to modify the nonlinear term u * u'
//    so that it is linear in the unknown u, and any other factors of u are
//    lagged.  One way to do this gives us the following equation:
//
//      -d/dx ( p(x) du/dx ) + q(x) * u + u * uold' = f(x)
//
//    where uold is the previous iterate.
//
//    Now we can formulate this system as a (linear) finite element problem
//
//      A * u = rhs
//
//    to be solved for the new approximate solution u.
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
//    Output, double F(NU), the right hand side of the linear
//    equations.
//
//    Input, double FOLD(NU), the solution value 
//    from the previous iteration,
//
//    Input, double H(N), the length of the subintervals.  
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
//    Input, int N.
//    The number of subintervals into which the interval
//    [XL,XR] is broken.
//
//    Input, int NL.
//    The number of basis functions used in a single
//    subinterval.  (NL-1) is the degree of the polynomials
//    used.  For this code, NL is fixed at 2, meaning that
//    piecewise linear functions are used as the basis.
//
//    Input, int NODE(NL,N).
//    For each subinterval I:
//    NODE(1,I) is the number of the left node, and
//    NODE(2,I) is the number of the right node.
//
//    Input, int NQUAD.
//    The number of quadrature points used in a subinterval.
//    This code uses NQUAD = 1.
//
//    Input, int NU.
//    NU is the number of unknowns in the linear system.
//    Depending on the value of IBC, there will be N-1,
//    N, or N+1 unknown values, which are the coefficients
//    of basis functions.
//
//    Input, int PROBLEM, indicates which problem to be solved.
//    * 1, u = x, p=1, q=0, f=x, u(0)=0, u'(1)=1.
//    * 2, u = 2*(1-cos(0.5*pi*x))/pi, p=1, q=0,
//    f = -0.5*pi*cos(0.5*pi*x) + 2*sin(0.5*pi*x)*(1-cos(0.5*pi*x)/pi
//    u(0) = 0, u'(1)=1.
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
//    Input, double XQUAD(N)
//    XQUAD(I) is the location of the single quadrature point
//    in interval I.
//
    {
        double phii = 0;
        double phiix = 0;
        double phij = 0;
        double phijx = 0;

        for (int i = 0; i < nu; i++)
        {
            f[i] = 0.0;
            adiag[i] = 0.0;
            aleft[i] = 0.0;
            arite[i] = 0.0;
        }

//
//  For element IE...
//
        for (int ie = 0; ie < n; ie++)
        {
            double he = h[ie];
            double xleft = xn[node[0 + ie * 2]];
            double xrite = xn[node[1 + ie * 2]];
//
//  For quadrature point IQ...
//
            for (int iq = 0; iq < nquad; iq++)
            {
                double xqe = xquad[ie];
//
//  Compute value of U for previous solution.
//
                double total = 0.0;

                int ig;
                int iu;
                for (int il = 1; il <= nl; il++)
                {
                    ig = node[il - 1 + ie * 2];
                    iu = indx[ig] - 1;

                    total += iu switch
                    {
                        < 0 when il == 1 => ul,
                        < 0 => ur,
                        _ => fold[iu]
                    };
                }

                double uold = total / nl;
//
//  Compute value of U' for previous solution.
//
                int jl = node[0 + ie * 2];
//
//  For basis function IL...
//
                for (int il = 1; il <= nl; il++)
                {
                    ig = node[il - 1 + ie * 2];
                    iu = indx[ig] - 1;

                    switch (iu)
                    {
                        case >= 0:
                        {
                            phi(il, xqe, ref phii, ref phiix, xleft, xrite);

                            f[iu] += he * phii * ff(xqe, problem);
//
//  Handle boundary conditions that prescribe the value of U'.
//
                            double x;
                            switch (ig)
                            {
                                case 0:
                                    x = 0.0;
                                    f[iu] -= pp(x, problem) * ul;
                                    break;
                                default:
                                {
                                    if (ig == n)
                                    {
                                        x = 1.0;
                                        f[iu] += pp(x, problem) * ur;
                                    }

                                    break;
                                }
                            }

//
//  For basis function JL...
//
                            for (jl = 1; jl <= nl; jl++)
                            {
                                int jg = node[jl - 1 + ie * 2];
                                int ju = indx[jg] - 1;

                                phi(jl, xqe, ref phij, ref phijx, xleft, xrite);

                                double aij = he * (pp(xqe, problem) * phiix * phijx
                                                   + qq(xqe, problem) * phii * phij
                                                   + uold * phii * phijx);

                                switch (ju)
                                {
                                    case < 0 when jg == 0:
                                        f[iu] -= aij * ul;
                                        break;
                                    case < 0:
                                    {
                                        if (jg == n)
                                        {
                                            f[iu] -= aij * ur;
                                        }

                                        break;
                                    }
                                    default:
                                    {
                                        if (iu == ju)
                                        {
                                            adiag[iu] += aij;
                                        }
                                        else if (ju < iu)
                                        {
                                            aleft[iu] += aij;
                                        }
                                        else
                                        {
                                            arite[iu] += aij;
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

    public static void compare(double[] f, int[] indx, int n, int nl, int[] node, int nprint,
        int nu, int problem, double ul, double ur, double xl, double[] xn,
        double xr )
//****************************************************************************80
//
//  Purpose:
//
//    COMPARE compares the computed and exact solutions.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 November 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Input, double F(NU), the solution of the linear equations.
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
//    Input, int N.
//    The number of subintervals into which the interval
//    [XL,XR] is broken.
//
//    Input, int NL.
//    The number of basis functions used in a single
//    subinterval.  (NL-1) is the degree of the polynomials
//    used.  For this code, NL is fixed at 2, meaning that
//    piecewise linear functions are used as the basis.
//
//    Input, int NODE(NL,N).
//    For each subinterval I:
//    NODE(1,I) is the number of the left node, and
//    NODE(2,I) is the number of the right node.
//
//    Input, int NPRINT.
//    The number of points at which the computed solution
//    should be printed out when compared to the exact solution.
//
//    Input, int NU.
//    NU is the number of unknowns in the linear system.
//    Depending on the value of IBC, there will be N-1,
//    N, or N+1 unknown values, which are the coefficients
//    of basis functions.
//
//    Input, int PROBLEM, indicates which problem to be solved.
//    * 1, u = x, p=1, q=0, f=x, u(0)=0, u'(1)=1.
//    * 2, u = 2*(1-cos(0.5*pi*x))/pi, p=1, q=0,
//    f = -0.5*pi*cos(0.5*pi*x) + 2*sin(0.5*pi*x)*(1-cos(0.5*pi*x)/pi
//    u(0) = 0, u'(1)=1.
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
//    Input, double XL.
//    XL is the left endpoint of the interval over which the
//    differential equation is being solved.
//
//    Input, double XN(0:N).
//    XN(I) is the location of the I-th node.  XN(0) is XL,
//    and XN(N) is XR.
//
//    Input, double XR.
//    XR is the right endpoint of the interval over which the
//    differential equation is being solved.
//
    {
        double phii = 0;
        double phiix = 0;

        Console.WriteLine("");
        Console.WriteLine("Compare computed and exact solutions:");
        Console.WriteLine("");
        Console.WriteLine("      X      Computed U      Exact U");
        Console.WriteLine("");

        for (int i = 1; i <= nprint; i++)
        {
            double x = ((nprint - i) * xl
                        + (i - 1) * xr)
                       / (nprint - 1);

            double ux = FEM_Test_Methods.u_exact(x, problem);

            double u = 0;
            for (int j = 1; j <= n; j++)
            {
                double xleft = xn[j - 1];
                double xrite = xn[j];
//
//  Search for the interval that X lies in.
//
                if (xleft <= x && x <= xrite)
                {
                    u = 0.0;

                    for (int k = 1; k <= nl; k++)
                    {
                        int ig = node[k - 1 + (j - 1) * 2];
                        int iu = indx[ig];
                        phi(k, x, ref phii, ref phiix, xleft, xrite);

                        switch (iu)
                        {
                            case <= 0 when j == 1 && k == 1:
                                u += ul * phii;
                                break;
                            case <= 0:
                            {
                                if (j == n && k == nl)
                                {
                                    u += ur * phii;
                                }

                                break;
                            }
                            default:
                                u += f[iu - 1] * phii;
                                break;
                        }
                    }

                    break;
                }
            }

            Console.WriteLine("  " + x.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + u.ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + ux.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
        }
    }

    public static double ff(double x, int problem)
//****************************************************************************80
//
//  Purpose:
//
//    FF returns the right hand side of the differential equation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 November 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//
//    Input, int PROBLEM, indicates which problem to be solved.
//    * 1, u = x, p=1, q=0, f=x, u(0)=0, u'(1)=1.
//    * 2, u = 2*(1-cos(0.5*pi*x))/pi, p=1, q=0,
//    f = -0.5*pi*cos(0.5*pi*x) + 2*sin(0.5*pi*x)*(1-cos(0.5*pi*x)/pi
//    u(0) = 0, u'(1)=1.
//
//    Output, double FF, the value of F(X).
//
    {
        double value = problem switch
        {
            //
            //  Test problem 1
            //
            1 => x,
            //
            //  Test problem 2
            //
            2 => -0.5 * Math.PI * Math.Cos(0.5 * Math.PI * x) +
                 2.0 * Math.Sin(0.5 * Math.PI * x) * (1.0 - Math.Cos(0.5 * Math.PI * x)) / Math.PI,
            _ => 0
        };

        return value;
    }

    public static void geometry(ref double[] h, int ibc, ref int[] indx, int nl, ref int[] node, int nsub,
        ref int nu, double xl, ref double[] xn, ref double[] xquad, double xr )
//****************************************************************************80
//
//  Purpose:
//
//    GEOMETRY sets up the geometry for the interval [XL,XR].
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 November 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Output, double H(N), the length of the subintervals.  
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
//    Input, int NL.
//    The number of basis functions used in a single
//    subinterval.  (NL-1) is the degree of the polynomials
//    used.  For this code, NL is fixed at 2, meaning that
//    piecewise linear functions are used as the basis.
//
//    Output, int NODE(NL,N).
//    For each subinterval I:
//    NODE(1,I) is the number of the left node, and
//    NODE(2,I) is the number of the right node.
//
//    Input, int NSUB.
//    The number of subintervals into which the interval
//    [XL,XR] is broken.
//
//    Output, int NU.
//    NU is the number of unknowns in the linear system.
//    Depending on the value of IBC, there will be N-1,
//    N, or N+1 unknown values, which are the coefficients
//    of basis functions.
//
//    Input, double XL.
//    XL is the left endpoint of the interval over which the
//    differential equation is being solved.
//
//    Output, double XN(0:N).
//    XN(I) is the location of the I-th node.  XN(0) is XL,
//    and XN(N) is XR.
//
//    Output, double XQUAD(N)
//    XQUAD(I) is the location of the single quadrature point
//    in interval I.
//
//    Input, double XR.
//    XR is the right endpoint of the interval over which the
//    differential equation is being solved.
//
    {
        int i;
//
//  Set the value of XN, the locations of the nodes.
//
        Console.WriteLine("");
        Console.WriteLine("  Node      Location");
        Console.WriteLine("");
        for (i = 0; i <= nsub; i++)
        {
            xn[i] = ((nsub - i) * xl
                     + i * xr)
                    / nsub;
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + xn[i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

//
//  Set the lengths of each subinterval.
//
        Console.WriteLine("");
        Console.WriteLine("Subint    Length");
        Console.WriteLine("");
        for ( i = 0; i < nsub; i++)
        {
            h[i] = xn[i + 1] - xn[i];
            Console.WriteLine("  " + (i + 1).ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + h[i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

//
//  Set the quadrature points, each of which is the midpoint
//  of its subinterval.
//
        Console.WriteLine("");
        Console.WriteLine("Subint    Quadrature point");
        Console.WriteLine("");
        for ( i = 0; i < nsub; i++)
        {
            xquad[i] = 0.5 * (xn[i] + xn[i + 1]);
            Console.WriteLine("  " + (i + 1).ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + xquad[i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }

//
//  Set the value of NODE, which records, for each interval,
//  the node numbers at the left and right.
//
        Console.WriteLine("");
        Console.WriteLine("Subint  Left Node  Right Node");
        Console.WriteLine("");
        for (i = 0; i < nsub; i++)
        {
            node[0 + i * 2] = i;
            node[1 + i * 2] = i + 1;
            Console.WriteLine("  " + (i + 1).ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + node[0 + i * 2].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + node[1 + i * 2].ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");
        }

//
//  Starting with node 0, see if an unknown is associated with
//  the node.  If so, give it an index.
//
        nu = 0;
//
//  Handle first node.
//
        i = 0;
        switch (ibc)
        {
            case 1:
            case 3:
                indx[i] = -1;
                break;
            default:
                nu += 1;
                indx[i] = nu;
                break;
        }

//
//  Handle nodes 1 through nsub-1
//
        for (i = 1; i < nsub; i++)
        {
            nu += 1;
            indx[i] = nu;
        }

//
//  Handle the last node.
//
        i = nsub;

        switch (ibc)
        {
            case 2:
            case 3:
                indx[i] = -1;
                break;
            default:
                nu += 1;
                indx[i] = nu;
                break;
        }

        Console.WriteLine("");
        Console.WriteLine("  Number of unknowns NU = " + nu + "");
        Console.WriteLine("");
        Console.WriteLine("  Node  Unknown");
        Console.WriteLine("");
        for (i = 0; i <= nsub; i++)
        {
            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + indx[i].ToString(CultureInfo.InvariantCulture).PadLeft(8) + "");
        }
    }

    public static void init(ref int ibc, ref int imax, ref int nprint, ref int nquad, ref int problem,
        ref double ul, ref double ur, ref double xl, ref double xr)
//****************************************************************************80
//
//  Purpose:
//
//    INIT initializes variables that define the problem.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 November 2006
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
//    Output, int *IMAX.
//    The number of Newton iterations to carry out.
//
//    Output, int *NPRINT.
//    The number of points at which the computed solution
//    should be printed out when compared to the exact solution.
//
//    Output, int *NQUAD.
//    The number of quadrature points used in a subinterval.
//    This code uses NQUAD = 1.
//
//    Output, int *PROBLEM, indicates which problem to be solved.
//    * 1, u = x, p=1, q=0, f=x, u(0)=0, u'(1)=1.
//    * 2, u = 2*(1-cos(0.5*pi*x))/pi, p=1, q=0,
//    f = -0.5*pi*cos(0.5*pi*x) + 2*sin(0.5*pi*x)*(1-cos(0.5*pi*x)/pi
//    u(0) = 0, u'(1)=1.
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
//    Output, double *XR.
//    XR is the right endpoint of the interval over which the
//    differential equation is being solved.
//
    {
        ibc = 1;
        imax = 10;
        nprint = 9;
        nquad = 1;
        problem = 2;
        ul = 0.0;
        ur = 1.0;
        xl = 0.0;
        xr = 1.0;
//
//  Print out the values that have been set.
//
        Console.WriteLine("");
        Console.WriteLine("  The equation is to be solved for");
        Console.WriteLine("  X greater than XL = " + xl + "");
        Console.WriteLine("  and less than XR = " + xr + "");
        Console.WriteLine("");
        Console.WriteLine("  The boundary conditions are:");
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
                Console.WriteLine("  At X = XR, U = " + ur + "");
                break;
            default:
                Console.WriteLine("  At X = XR, U' = " + ur + "");
                break;
        }

        switch (problem)
        {
            case 1:
                Console.WriteLine("");
                Console.WriteLine("  This is test problem #1:");
                Console.WriteLine("");
                Console.WriteLine("  P(X) = 1, Q(X) = 0, F(X) = X.");
                Console.WriteLine("  Boundary conditions: U(0) = 0, U''(1) = 1.");
                Console.WriteLine("");
                Console.WriteLine("  The exact solution is U(X) = X");
                break;
            case 2:
                Console.WriteLine("");
                Console.WriteLine("  This is test problem #2:");
                Console.WriteLine("");
                Console.WriteLine("  P(X) = 1, Q(X) = 0, ");
                Console.WriteLine("  F(X) = -0.5*pi*cos(0.5*pi*X)");
                Console.WriteLine("        + 2*sin(0.5*pi*X)*(1-cos(0.5*pi*X)/pi.");
                Console.WriteLine("  Boundary conditions: U(0) = 0, U''(1) = 1.");
                Console.WriteLine("");
                Console.WriteLine("  The exact solution is U(X) = 2*(1-cos(pi*x/2))/pi");
                break;
        }

        Console.WriteLine("");
        Console.WriteLine("  Number of quadrature points per element is " + nquad + "");
        Console.WriteLine("  Number of iterations is " + imax + "");
    }

    public static void output(double[] f, int ibc, int[] indx, int nsub, int nu, double ul,
        double ur, double[] xn )
//****************************************************************************80
//
//  Purpose:
//
//    OUTPUT prints out the computed solution at the nodes.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 November 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Input, double F(NU), the solution of the linear equations.
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
//    int NSUB.
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
        double u;

        Console.WriteLine("");
        Console.WriteLine("Computed solution:");
        Console.WriteLine("");
        Console.WriteLine("Node    X(I)        U(X(I))");
        Console.WriteLine("");

        for (int i = 0; i <= nsub; i++)
        {
            switch (i)
            {
                case 0 when ibc == 1 || ibc == 3:
                    u = ul;
                    break;
                case 0:
                    u = f[indx[i] - 1];
                    break;
                default:
                {
                    if (i == nsub)
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

            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(6)
                                   + "  " + xn[i].ToString(CultureInfo.InvariantCulture).PadLeft(12)
                                   + "  " + u.ToString(CultureInfo.InvariantCulture).PadLeft(12) + "");
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
//    In any interval, there are just two basis functions.  The first
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
//    02 November 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Input, int IL, the index of the basis function.
//    1, the function which is 1 at XLEFT and 0 at XRITE.
//    2, the function which is 0 at XLEFT and 1 at XRITE.
//
//    Input, double X, the evaluation point.
//
//    Output, double PHII, PHIIX, the value of the
//    basis function and its derivative at X.
//
//    Input, double XLEFT, XRITE, the left and right
//    endpoints of the interval.
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
//  If X is outside of the interval, just set everything to 0.
//
        else
        {
            phii = 0.0;
            phiix = 0.0;
        }
    }

    public static double pp(double x, int problem)
//****************************************************************************80
//
//  Purpose:
//
//    PP evaluates the function P in the differential equation.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    02 November 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//
//    Input, int PROBLEM, indicates which problem to be solved.
//    * 1, u = x, p=1, q=0, f=x, u(0)=0, u'(1)=1.
//    * 2, u = 2*(1-cos(0.5*pi*x))/pi, p=1, q=0,
//    f = -0.5*pi*cos(0.5*pi*x) + 2*sin(0.5*pi*x)*(1-cos(0.5*pi*x)/pi
//    u(0) = 0, u'(1)=1.
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

    public static void prsys(double[] adiag, double[] aleft, double[] arite, double[] f,
        int nu )

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
//    31 October 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameter:
//
//    Input, double ADIAG(NU), the "diagonal" coefficients.  That is, 
//    ADIAG(I) is the coefficient of the I-th unknown in the I-th equation.
//
//    Input, double ALEFT(NU), the "left hand" coefficients.  That is, ALEFT(I) 
//    is the coefficient of the (I-1)-th unknown in the I-th equation.
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
//    Input, int NU.
//    NU is the number of unknowns in the linear system.
//    Depending on the value of IBC, there will be NSUB-1,
//    NSUB, or NSUB+1 unknown values, which are the coefficients
//    of basis functions.
//
    {
        int i;

        Console.WriteLine("");
        Console.WriteLine("Printout of tridiagonal linear system:");
        Console.WriteLine("");
        Console.WriteLine("Equation  ALEFT  ADIAG  ARITE  RHS");
        Console.WriteLine("");

        for (i = 0; i < nu; i++)
        {
            Console.WriteLine("  " + (i + 1).ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + aleft[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + adiag[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + arite[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + f[i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
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
//    02 November 2006
//
//  Author:
//
//    C++ version by John Burkardt
//
//  Parameters:
//
//    Input, double X, the evaluation point.
//
//    Input, int PROBLEM, indicates which problem to be solved.
//    * 1, u = x, p=1, q=0, f=x, u(0)=0, u'(1)=1.
//    * 2, u = 2*(1-cos(0.5*pi*x))/pi, p=1, q=0,
//    f = -0.5*pi*cos(0.5*pi*x) + 2*sin(0.5*pi*x)*(1-cos(0.5*pi*x)/pi
//    u(0) = 0, u'(1)=1.
//
//    Output, double QQ, the value of Q(X).
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
                value = 0.0;
                break;
        }

        return value;
    }

    public static void solve(ref double[] adiag, ref double[] aleft, ref double[] arite, ref double[] f,
        int nu )
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
//    31 October 2006
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
//    system to be solved.
//    On output, F contains the solution of the linear system.
//
//    Input, int NU, the number of equations to be solved.
//
    {
//
//  Carry out Gauss elimination on the matrix, saving information
//  needed for the backsolve.
//
        arite[0] /= adiag[0];

        for (int i = 1; i < nu - 1; i++)
        {
            adiag[i] -= aleft[i] * arite[i - 1];
            arite[i] /= adiag[i];
        }

        adiag[nu - 1] -= aleft[nu - 1] * arite[nu - 2];
//
//  Carry out the same elimination steps on F that were done to the
//  matrix.
//
        f[0] /= adiag[0];
        for (int i = 1; i < nu; i++)
        {
            f[i] = (f[i] - aleft[i] * f[i - 1]) / adiag[i];
        }

//
//  And now carry out the steps of "back substitution".
//
        for (int i = nu - 2; 0 <= i; i--)
        {
            f[i] -= arite[i] * f[i + 1];
        }
    }
}