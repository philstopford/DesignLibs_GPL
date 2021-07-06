using System;
using Burkardt.FEM;
using Burkardt.Types;
using Burkardt.PolynomialNS;

namespace Burkardt
{
    public static class Map
    {
        public static double[] map(string code, int element_order)

            //****************************************************************************80
            //
            //  Purpose: 
            //
            //    MAP returns the interpolation matrix for any available element.
            //
            //  Formula:
            //
            //    For an element of order N, we suppose we are given N items of data 
            //    Q associated with the nodes.
            //
            //   Let PHI(J)(R,S) be the Lagrange basis polynomial associated with 
            //   node J.  PHI(J)(R,S) is 1 at node J, and 0 at each of the other nodes.
            //
            //   Let P(R,S) be the polynomial of N terms which interpolates the
            //   data Q, that is,
            //
            //      P(R(J),S(J)) = Q(J)
            //
            //   where the coordinates of node J are (R(J),S(J)).  Then we know
            //   that we can write
            //
            //     P(R,S) = sum ( 1 <= J <= N ) Q(J) * PHI(J)(R,S)
            //
            //   But P(R,S) also has a standard representation as
            //
            //     P(R,S) = sum ( 1 <= I <= N ) A(I) * R**REXP(I) * S**SEXP(I)
            //
            //   where REXP(I) and SEXP(I) are the exponents of R and S and
            //   the A(I) are the appropriate coefficients.
            //
            //   The interpolation matrix W allows us to immediately compute
            //   the standard basis coefficients A from the data Q to be interpolated
            //   using the formula:
            //
            //      A(I) = sum ( 1 <= J <= N ) W(I,J) * Q(J)
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
            //    Input, string CODE, identifies the element.
            //    Legal values include 'Q4', 'Q8', 'Q9', 'Q12', 'Q16', 'QL', 
            //    'T3', 'T4', 'T6' and 'T10'.
            //
            //    Input, int N, the order associated with the code.
            //
            //    Output, double MAP[N*N], the interpolation matrix.
            //
        {
            double area = 0;
            int i;
            int info;
            int j;
            int[] pivot;
            double[] r;
            int[] rexp;
            double rfact;
            double[] s;
            int[] sexp;
            double sfact;
            double[] v;
            double[] w;

            pivot = new int[element_order];
            r = new double[element_order];
            rexp = new int[element_order];
            s = new double[element_order];
            sexp = new int[element_order];
            v = new double[element_order * element_order];
            //
            //  Get the (R,S) location of the nodes.
            //
            NodeReference.node_reference(code, ref r, ref s, ref area);
            //
            //  Get the associated monomials.
            //
            Polynomial.poly(code, ref rexp, ref sexp);
            //
            //  Set up the Vandermonde matrix.
            //  Factors of the form 0**0 are to be understood as 1.
            //
            for (i = 0; i < element_order; i++)
            {
                for (j = 0; j < element_order; j++)
                {
                    if (rexp[j] == 0)
                    {
                        rfact = 1.0;
                    }
                    else
                    {
                        rfact = Math.Pow(r[i], rexp[j]);
                    }

                    if (sexp[j] == 0)
                    {
                        sfact = 1.0;
                    }
                    else
                    {
                        sfact = Math.Pow(s[i], sexp[j]);
                    }

                    v[i + j * element_order] = rfact * sfact;
                }
            }

            //
            //  Factor the Vandermonde matrix.
            //
            info = typeMethods.r8ge_fa(element_order, ref v, ref pivot);

            if (info != 0)
            {
                Console.WriteLine("");
                Console.WriteLine("MAP - Fatal error!");
                Console.WriteLine("  The Vandermonde matrix is singular.");
                return null;
            }

            //
            //  Invert the Vandermonde matrix.
            //
            w = typeMethods.r8ge_inverse(element_order, v, pivot);

            return w;
        }

        public static void map_test(string code)

            //****************************************************************************80
            //
            //  Purpose: 
            //
            //    MAP_TEST tests the map routines.
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
            //    Input, string CODE, the code for the element.
            //    Legal values include 'Q4', 'Q8', 'Q9', 'Q12', 'Q16', 'QL', 
            //    'T3', 'T4', 'T6' and 'T10'.
            //
        {
            int element_order;
            double[] w;

            Console.WriteLine("");
            Console.WriteLine("  MAP_TEST: The interpolation matrix for element " + code + "");

            element_order = FEM.Order.order_code(code);

            w = map(code, element_order);

            typeMethods.r8mat_print(element_order, element_order, w,
                "  The interpolation matrix:");

        }

        public static void xy_to_rs_map ( double[] t, ref double a, ref double b, ref double c, ref double d, 
        ref double e, ref double f )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    XY_TO_RS_MAP returns the linear map from physical to reference triangle.
        //
        //  Discussion:
        //
        //    Given the vertices T of an arbitrary triangle in the (X,Y) coordinate
        //    system, this function returns the coefficients of the linear map
        //    that sends the vertices of T to (0,0), (1,0) and (0,1) respectively
        //    in the reference triangle with coordinates (R,S):
        //
        //      R = A + B * X + C * Y;
        //      S = D + E * X + F * Y.
        //
        //  Reference Element T3:
        //
        //    |
        //    1  3
        //    |  |.
        //    |  | .
        //    S  |  .
        //    |  |   .
        //    |  |    .
        //    0  1-----2
        //    |
        //    +--0--R--1-->
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    21 April 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double T[2*3], the X and Y coordinates
        //    of the vertices.  The vertices are assumed to be the images of
        //    (0,0), (1,0) and (0,1) respectively.
        //
        //    Output, double &A, &B, &C, &D, &E, &F, the mapping coefficients.
        //
        {
            double g;

            g =    ( ( t[1+2*2] - t[1+0*2] ) * ( t[0+1*2] - t[0+0*2] )   
                     - ( t[0+2*2] - t[0+0*2] ) * ( t[1+1*2] - t[1+0*2] ) );

            a = ( - ( t[1+2*2] - t[1+0*2] ) * t[0+0*2]  
                  + ( t[0+2*2] - t[0+0*2] ) * t[1+0*2] ) / g;

            b =     ( t[1+2*2] - t[1+0*2] ) / g;

            c =   - ( t[0+2*2] - t[0+0*2] ) / g;

            d = (   ( t[1+1*2] - t[1+0*2] ) * t[0+0*2] 
                    - ( t[0+1*2] - t[0+0*2] ) * t[1+0*2] ) / g;

            e =   - ( t[1+1*2] - t[1+0*2] ) / g;

            f =     ( t[0+1*2] - t[0+0*2] ) / g;

            return;
        }
        
        public static void rs_to_xy_map ( double[] t, ref double a, ref double b, ref double c, ref double d, 
                ref double e, ref double f )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RS_TO_XY_MAP returns the linear map from reference to physical triangle.
        //
        //  Discussion:
        //
        //    This function returns the coefficients of the linear map that sends
        //    the vertices of the reference triangle, (0,0), (1,0) and (0,1), to
        //    the vertices of a physical triangle T, of the form:
        //
        //      X = A + B * R + C * S;
        //      Y = D + E * R + F * S.
        //
        //  Reference Element:
        //
        //    |
        //    1  3
        //    |  |.
        //    |  | .
        //    S  |  .
        //    |  |   .
        //    |  |    .
        //    0  1-----2
        //    |
        //    +--0--R--1-->
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 April 2015
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double T[2,3], the coordinates of the vertices.  The
        //    vertices are assumed to be the images of (0,0), (1,0) and (0,1) 
        //    respectively.
        //
        //    Output, double &A, &B, &C, &D, &E, &F, the mapping coefficients.
        //
        {
            a = t[0+0*2];
            b = t[0+1*2] - t[0+0*2];
            c = t[0+2*2] - t[0+0*2];

            d = t[1+0*2];
            e = t[1+1*2] - t[1+0*2];
            f = t[1+2*2] - t[1+0*2];
        }
    }
}