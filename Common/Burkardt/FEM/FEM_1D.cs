using System;
using System.Globalization;

namespace Burkardt.FEM;

public static class FEM_1D
{

    public static void assemble(ref double[] adiag, ref double[] aleft, ref double[] arite, ref double[] f,
        double[] h, int[] indx, int nl, int[] node, int nu, int nquad, int nsub,
        double ul, double ur, double[] xn, double[] xquad)

//****************************************************************************80
//
//  Purpose:
//
//    ASSEMBLE assembles the matrix and right-hand-side of the linear system.
//
//  Discussion:
//
//    The linear system has the form:
//
//      K * C = F
//
//    that is to be solved for the coefficients C.
//
//    Numerical integration is used to compute the entries of K and F.
//
//    Note that a 1 point quadrature rule, which is sometimes used to
//    assemble the matrix and right hand side, is just barely accurate
//    enough for simple problems.  If you want better results, you
//    should use a quadrature rule that is more accurate.
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
//    C++ version by John Burkardt.
//
//  Parameters:
//
//    Output, double ADIAG(NU), the "diagonal" coefficients.  That is, 
//    ADIAG(I) is the coefficient of the I-th unknown in the I-th equation.
//
//    Output, double ALEFT(NU), the "left hand" coefficients.  That is, 
//    ALEFT(I) is the coefficient of the (I-1)-th unknown in the I-th equation.
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
//    Input, double H(NSUB)
//    H(I) is the length of subinterval I.  This code uses
//    equal spacing for all the subintervals.
//
//    Input, int INDX[NSUB+1].
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
//    Input, int NL.
//    The number of basis functions used in a single
//    subinterval.  (NL-1) is the degree of the polynomials
//    used.  For this code, NL is fixed at 2, meaning that
//    piecewise linear functions are used as the basis.
//
//    Input, int NODE[NL*NSUB].
//    For each subinterval I:
//    NODE[0+I*2] is the number of the left node, and
//    NODE[1+I*2] is the number of the right node.
//
//    Input, int NU.
//    NU is the number of unknowns in the linear system.
//    Depending on the value of IBC, there will be NSUB-1,
//    NSUB, or NSUB+1 unknown values, which are the coefficients
//    of basis functions.
//
//    Input, int NQUAD.
//    The number of quadrature points used in a subinterval.
//    This code uses NQUAD = 1.
//
//    Input, int NSUB.
//    The number of subintervals into which the interval [XL,XR] is broken.
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
//    Input, double XR.
//    XR is the right endpoint of the interval over which the
//    differential equation is being solved.
//
    {
        double phii = 0;
        double phiix = 0;
        double phij = 0;
        double phijx = 0;
        //
//  Zero out the arrays that hold the coefficients of the matrix
//  and the right hand side.
//
        for (int i = 0; i < nu; i++)
        {
            f[i] = 0.0;
            adiag[i] = 0.0;
            aleft[i] = 0.0;
            arite[i] = 0.0;
        }

//
//  For interval number IE,
//
        for (int ie = 0; ie < nsub; ie++)
        {
            double he = h[ie];
            double xleft = xn[node[0 + ie * 2]];
            double xrite = xn[node[1 + ie * 2]];
//
//  consider each quadrature point IQ,
//
            for (int iq = 0; iq < nquad; iq++)
            {
                double xquade = xquad[ie];
                //
//  and evaluate the integrals associated with the basis functions
//  for the left, and for the right nodes.
//
                for (int il = 1; il <= nl; il++)
                {
                    int ig = node[il - 1 + ie * 2];
                    int iu = indx[ig] - 1;

                    switch (iu)
                    {
                        case >= 0:
                        {
                            phi(il, xquade, ref phii, ref phiix, xleft, xrite);
                            f[iu] += he * ff(xquade) * phii;
//
//  Take care of boundary nodes at which U' was specified.
//
                            double x;
                            switch (ig)
                            {
                                case 0:
                                    x = 0.0;
                                    f[iu] -= pp(x) * ul;
                                    break;
                                default:
                                {
                                    if (ig == nsub)
                                    {
                                        x = 1.0;
                                        f[iu] += pp(x) * ur;
                                    }

                                    break;
                                }
                            }

//
//  Evaluate the integrals that take a product of the basis
//  function times itself, or times the other basis function
//  that is nonzero in this interval.
//
                            for (int jl = 1; jl <= nl; jl++)
                            {
                                int jg = node[jl - 1 + ie * 2];
                                int ju = indx[jg] - 1;

                                phi(jl, xquade, ref phij, ref phijx, xleft, xrite);

                                double aij = he * (pp(xquade) * phiix * phijx
                                                   + qq(xquade) * phii * phij);
                                switch (ju)
                                {
    
//
                                    //  If there is no variable associated with the node, then it's
                                    //  a specified boundary value, so we multiply the coefficient
                                    //  times the specified boundary value and subtract it from the
                                    //  right hand side.
                                    //
                                    case < 0 when jg == 0:
                                        f[iu] -= aij * ul;
                                        break;
case < 0:
{
    if (jg == nsub)
                                        {
                                            f[iu] -= aij * ur;
                                        }

    break;
}

//
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

    public static double ff ( double x )
//****************************************************************************80
//
//  Purpose:
//
//    FF evaluates the right hand side function.
//
//  Discussion:
//
//    This routine evaluates the function F(X) in the differential equation.
//
//      -d/dx (p du/dx) + q u  =  f
//
//    at the point X.
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
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the function.
//
//    Output, double FF, the value of the function.
//
    {
        double value = 0.0;

        return value;
    }
        
    public static void geometry ( double[] h, int ibc, ref int[] indx, int nl, ref int[] node, int nsub, 
        ref int nu, double xl, ref double[] xn, ref double[] xquad, double xr )
//****************************************************************************80
//
//  Purpose: 
//
//    GEOMETRY sets up the geometry for the interval [XL,XR].
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
//    Output, double H(NSUB)
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
//    Output, int INDX[NSUB+1].
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
//    Input, int NL.
//    The number of basis functions used in a single
//    subinterval.  (NL-1) is the degree of the polynomials
//    used.  For this code, NL is fixed at 2, meaning that
//    piecewise linear functions are used as the basis.
//
//    Output, int NODE[NL*NSUB].
//    For each subinterval I:
//    NODE[0+I*2] is the number of the left node, and
//    NODE[1+I*2] is the number of the right node.
//
//    Input, int NSUB.
//    The number of subintervals into which the interval [XL,XR] is broken.
//
//    Output, int *NU.
//    NU is the number of unknowns in the linear system.
//    Depending on the value of IBC, there will be NSUB-1,
//    NSUB, or NSUB+1 unknown values, which are the coefficients
//    of basis functions.
//
//    Input, double XL.
//    XL is the left endpoint of the interval over which the
//    differential equation is being solved.
//
//    Output, double XN(0:NSUB).
//    XN(I) is the location of the I-th node.  XN(0) is XL,
//    and XN(NSUB) is XR.
//
//    Output, double XQUAD(NSUB)
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
        for (i = 0; i < nsub; i++)
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
        for (i = 0; i < nsub; i++)
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


    public static void init(ref int ibc, ref int nquad, ref double ul, ref double ur, ref double xl,
        ref double xr)
//****************************************************************************80
//
//  Purpose: 
//
//    INIT assigns values to variables which define the problem.
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
//    Output, int *NQUAD.
//    The number of quadrature points used in a subinterval.
//    This code uses NQUAD = 1.
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
//
//  IBC declares what the boundary conditions are.
//
        ibc = 1;
//
//  NQUAD is the number of quadrature points per subinterval.
//  The program as currently written cannot handle any value for
//  NQUAD except 1//
//
        nquad = 1;
//
//  Set the values of U or U' at the endpoints.
//
        ul = 0.0;
        ur = 1.0;
//
//  Define the location of the endpoints of the interval.
//
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

        Console.WriteLine("");
        Console.WriteLine("  Number of quadrature points per element is " + nquad + "");
    }

    public static void output ( double[] f, int ibc, int[] indx, int nsub, int nu, double ul, 
        double ur, double[] xn )

//****************************************************************************80
//
//  Purpose:
//
//    OUTPUT prints out the computed solution.
//
//  Discussion:
//
//    We simply print out the solution vector F, except that, for
//    certain boundary conditions, we are going to have to get the
//    value of the solution at XL or XR by using the specified
//    boundary value.
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
//    Input, int INDX[NSUB+1].
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
//    Input, int NSUB.
//    The number of subintervals into which the interval [XL,XR] is broken.
//
//    Input, int NU.
//    NU is the number of unknowns in the linear system.
//    Depending on the value of IBC, there will be NSUB-1,
//    NSUB, or NSUB+1 unknown values, which are the coefficients
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
//    Input, double XN(0:NSUB).
//    XN(I) is the location of the I-th node.  XN(0) is XL,
//    and XN(NSUB) is XR.
//
    {
        double u;

        Console.WriteLine("");
        Console.WriteLine("  Computed solution coefficients:");
        Console.WriteLine("");
        Console.WriteLine("  Node    X(I)        U(X(I))");
        Console.WriteLine("");

        for (int i = 0; i <= nsub; i++)
        {
            switch (i)
            {
    
//
                //  If we're at the first node, check the boundary condition.
                //
                case 0 when ibc == 1 || ibc == 3:
                    u = ul;
                    break;
case 0:
                    u = f[indx[i] - 1];
                    break;

//
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
//
//  Any other node, we're sure the value is stored in F.
//
                    else
                    {
                        u = f[indx[i] - 1];
                    }

                    break;
                }
            }

            Console.WriteLine("  " + i.ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + xn[i].ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + u.ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    public static void phi ( int il, double x, ref double phii, ref double phiix, double xleft, 
        double xrite )
//****************************************************************************80
//
//  Purpose:
//
//    PHI evaluates a linear basis function and its derivative.
//
//  Discussion:
//
//    The evaluation is done at a point X in an interval [XLEFT,XRITE].
//
//    In this interval, there are just two nonzero basis functions.
//    The first basis function is a line which is 1 at the left
//    endpoint and 0 at the right.  The second basis function is 0 at
//    the left endpoint and 1 at the right.
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
//    Input, int IL, the index of the basis function.
//    1, the function which is 1 at XLEFT and 0 at XRITE.
//    2, the function which is 0 at XLEFT and 1 at XRITE.
//
//    Input, double X, the evaluation point.
//
//    Output, double *PHII, *PHIIX, the value of the
//    basis function and its derivative at X.
//
//    Input, double XLEFT, XRITE, the left and right
//    endpoints of the interval.
//
    {
        if ( xleft <= x && x <= xrite )
        {
            switch (il)
            {
                case 1:
                    phii = ( xrite - x ) / ( xrite - xleft );
                    phiix =         -1.0 / ( xrite - xleft );
                    break;
                default:
                    phii = ( x - xleft ) / ( xrite - xleft );
                    phiix = 1.0          / ( xrite - xleft );
                    break;
            }
        }
//
//  If X is outside of the interval, just set everything to 0.
//
        else
        {
            phii  = 0.0;
            phiix = 0.0;
        }
    }
        
        
    public static double pp ( double x )
//****************************************************************************80
//
//  Purpose:
//
//    PP evaluates the function P in the differential equation.
//
//  Discussion:
//
//    The function P appears in the differential equation as;
//
//      - d/dx (p du/dx) + q u  =  f
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
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the function.
//
//    Output, double PP, the value of the function.
//
    {
        double value = 1.0;

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

        Console.WriteLine("");
        Console.WriteLine("Printout of tridiagonal linear system:");
        Console.WriteLine("");
        Console.WriteLine("Equation  ALEFT  ADIAG  ARITE  RHS");
        Console.WriteLine("");

        for (int i = 0; i < nu; i++)
        {
            Console.WriteLine("  " + (i + 1).ToString(CultureInfo.InvariantCulture).PadLeft(8)
                                   + "  " + aleft[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + adiag[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + arite[i].ToString(CultureInfo.InvariantCulture).PadLeft(14)
                                   + "  " + f[i].ToString(CultureInfo.InvariantCulture).PadLeft(14) + "");
        }
    }

    public static double qq ( double x )
//****************************************************************************80
//
//  Purpose: 
//
//    QQ evaluates the function Q in the differential equation.
//
//  Discussion:
//
//    The function Q appears in the differential equation as:
//
//      - d/dx (p du/dx) + q u  =  f
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
//    John Burkardt
//
//  Parameters:
//
//    Input, double X, the argument of the function.
//
//    Output, double QQ, the value of the function.
//
    {
        double value = 0.0;

        return value;
    }
        
    public static void solve ( ref double[] adiag, ref double[] aleft, ref double[] arite, ref double[] f, 
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

        for (int i = 1; i < nu - 1; i++ )
        {
            adiag[i] -= aleft[i] * arite[i-1];
            arite[i] /= adiag[i];
        }
        adiag[nu-1] -= aleft[nu-1] * arite[nu-2];
//
//  Carry out the same elimination steps on F that were done to the
//  matrix.
//
        f[0] /= adiag[0];
        for (int i = 1; i < nu; i++ )
        {
            f[i] = ( f[i] - aleft[i] * f[i-1] ) / adiag[i];
        }
//
//  And now carry out the steps of "back substitution".
//
        for (int i = nu - 2; 0 <= i; i-- )
        {
            f[i] -= arite[i] * f[i+1];
        }
    }
}