namespace Burkardt
{
    public static class Map
    {
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