using System;
using System.IO;
using Burkardt.Types;

namespace Burkardt;

public static class Helpers
{
    public static double arc_sine ( double s )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ARC_SINE computes the arc sine function, with argument truncation.
        //
        //  Discussion:
        //
        //    If you call your system ASIN routine with an input argument that is
        //    outside the range [-1.0, 1.0 ], you may get an unpleasant surprise.
        //    This routine truncates arguments outside the range.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 June 2002
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double S, the argument, the sine of an angle.
        //
        //    Output, double ARC_SINE, an angle whose sine is S.
        //
    {
        double angle = s switch
        {
            <= -1.0 => -Math.PI / 2.0,
            >= 1.0 => Math.PI / 2.0,
            _ => Math.Asin(s)
        };
        return angle;
    }

    public static double atan4(double y, double x)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ATAN4 computes the inverse tangent of the ratio Y / X.
        //
        //  Discussion:
        //
        //    ATAN4 returns an angle whose tangent is ( Y / X ), a job which
        //    the built in functions ATAN and ATAN2 already do.
        //
        //    However:
        //
        //    * ATAN4 always returns a positive angle, between 0 and 2 PI,
        //      while ATAN and ATAN2 return angles in the interval [-PI/2,+PI/2]
        //      and [-PI,+PI] respectively;
        //
        //    * ATAN4 accounts for the signs of X and Y, (as does ATAN2).  The ATAN
        //     function by contrast always returns an angle in the first or fourth
        //     quadrants.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    13 June 2002
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double Y, X, two quantities which represent the tangent of
        //    an angle.  If Y is not zero, then the tangent is (Y/X).
        //
        //    Output, double ATAN4, an angle between 0 and 2 * PI, whose tangent is
        //    (Y/X), and which lies in the appropriate quadrant so that the signs
        //    of its cosine and sine match those of X and Y.
        //
    {
        switch (x)
        {
            //
            //  Special cases:
            //
            case 0.0:
                switch (y)
                {
                    case > 0.0:
                        return Math.PI / 2.0;
                    case < 0.0:
                        return 3.0 * Math.PI / 2.0;
                    case 0.0:
                        return 0.0;
                }

                break;
            default:
            {
                switch (y)
                {
                    case 0.0:
                        switch (x)
                        {
                            case > 0.0:
                                return 0.0;
                            case < 0.0:
                                return Math.PI;
                        }

                        break;
                }

                break;
            }
        }

        return x switch
        {
            //
            //  We assume that ATAN2 is reliable when both arguments are positive.
            //
            > 0.0 when 0.0 < y => Math.Atan2(y, x),
            < 0.0 when 0.0 < y => Math.PI - Math.Atan2(y, -x),
            < 0.0 when y < 0.0 => Math.PI + Math.Atan2(-y, -x),
            > 0.0 when y < 0.0 => 2.0 * Math.PI - Math.Atan2(-y, x),
            _ => 0.0
        };
    }

    public static double exact ( double x, double y )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    EXACT calculates the exact solution and its first derivatives.
        //
        //  Discussion:
        //
        //    The function specified here depends on the problem being
        //    solved.  The user must be sure to change both EXACT and RHS
        //    or the program will have inconsistent data.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 December 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double X, Y, the coordinates of a point
        //    in the region, at which the exact solution is to be evaluated.
        //
        //    Output, double EXACT, the value of the exact solution.
        //
    {
        double u = Math.Sin ( Math.PI * x ) * Math.Sin ( Math.PI * y ) + x;

        return u;
    }
        
    public static void mult_givens ( double c, double s, int k, ref double[] g )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MULT_GIVENS applies a Givens rotation to two successive entries of a vector.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    08 August 2006
        //
        //  Author:
        //
        //    Original C version by Lili Ju.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Richard Barrett, Michael Berry, Tony Chan, James Demmel,
        //    June Donato, Jack Dongarra, Victor Eijkhout, Roidan Pozo,
        //    Charles Romine, Henk van der Vorst,
        //    Templates for the Solution of Linear Systems:
        //    Building Blocks for Iterative Methods,
        //    SIAM, 1994,
        //    ISBN: 0898714710,
        //    LC: QA297.8.T45.
        //
        //    Tim Kelley,
        //    Iterative Methods for Linear and Nonlinear Equations,
        //    SIAM, 2004,
        //    ISBN: 0898713528,
        //    LC: QA297.8.K45.
        //
        //    Yousef Saad,
        //    Iterative Methods for Sparse Linear Systems,
        //    Second Edition,
        //    SIAM, 2003,
        //    ISBN: 0898715342,
        //    LC: QA188.S17.
        //
        //  Parameters:
        //
        //    Input, double C, S, the cosine and sine of a Givens
        //    rotation.
        //
        //    Input, int K, indicates the location of the first vector entry.
        //
        //    Input/output, double G[K+2], the vector to be modified.  On output,
        //    the Givens rotation has been applied to entries G(K) and G(K+1).
        //
    {
        double g1 = c * g[k] - s * g[k+1];
        double g2 = s * g[k] + c * g[k+1];

        g[k]   = g1;
        g[k+1] = g2;
    }
    public static double degrees_to_radians ( double angle )

        //****************************************************************************80
        //
        //  Purpose: 
        //
        //    DEGREES_TO_RADIANS converts an angle from degrees to radians.
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
        //    Input, double ANGLE, an angle in degrees.
        //
        //    Output, double DEGREES_TO_RADIANS, the equivalent angle
        //    in radians.
        //
    {
        return angle * Math.PI / 180.0;

    }
        
    public static double dms_to_radians ( int degrees, int minutes, int seconds )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    DMS_TO_RADIANS converts an angle from degrees/minutes/seconds to radians.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DEGREES, MINUTES, SECONDS, an angle in degrees, minutes,
        //    and seconds.
        //
        //    Output, double DMS_TO_RADIANS, the equivalent angle in radians.
        //
    {
        double angle = degrees
                       + ( minutes
                           + seconds / 60.0 ) / 60.0;

        double radians = angle / 180.0 * Math.PI;

        return radians;
    }
        
    public static double enorm ( int n, double[] x, int xIndex = 0 )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    enorm() returns the Euclidean norm of a vector.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 April 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Input:
        //
        //    int N, the number of entries in A.
        //
        //    double X[N], the vector whose norm is desired.
        //
        //  Output:
        //
        //    double ENORM, the norm of X.
        //
    {
        int i;

        double value = 0.0;
        for ( i = 0; i < n; i++ )
        {
            value += x[xIndex + i] * x[xIndex + i];
        }
        value = Math.Sqrt ( value );
        return value;
    }
        
    public static double enorm0_nd ( int dim_num, double[] x, double[] y )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ENORM0_ND computes the Euclidean norm of a (X-Y) in N space.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    18 April 1999
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the dimension of the space.
        //
        //    Input, double X[DIM_NUM], Y[DIM_NUM], the coordinates of the vectors.
        //
        //    Output, double ENORM0_ND, the Euclidean norm of the vector.
        //
    {
        int i;

        double value = 0.0;

        for ( i = 0; i < dim_num; i++ )
        {
            value += ( x[i] - y[i] ) * ( x[i] - y[i] );
        }

        value = Math.Sqrt ( value );

        return value;
    }
        
    public static void normalize ( int n, ref double[] x )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    NORMALIZE scales a vector X so its entries sum to 1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 February 2016
        //
        //  Author:
        //
        //    Original C version by Warren Smith.
        //    C++ version by John Burkardt.
        //
        //  Parameters:
        //
        //    Input, unsigned int N, indicates the size of X.
        //
        //    Input/output, double X[N+2], the vector to be normalized.
        //    Entries X[1] through X[N] will sum to 1 on output.
        //
    {
        int i;
        //
        //  Sum X.
        //
        double sum = 0.0;
        for ( i = 1; i <= n; i++ )
        {
            sum += Math.Abs ( x[i] );
        }
        //
        //  Normalize so that the new sum of X will be 1.
        //
        sum = 1.0 / sum;
        for ( i = 1; i <= n; i++ )
        {
            x[i] *= sum; 
        }
    }

    public static int points_point_near_naive_nd ( int dim_num, int nset, double[] pset,
            double[] ptest, ref double d_min )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POINTS_POINT_NEAR_NAIVE_ND finds the nearest point to a given point in ND.
        //
        //  Discussion:
        //
        //    A naive algorithm is used.  The distance to every point is calculated,
        //    in order to determine the smallest.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 June 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int DIM_NUM, the spatial dimension.
        //
        //    Input, int NSET, the number of points in the set.
        //
        //    Input, double PSET[DIM_NUM*NSET], the coordinates of the points
        //    in the set.
        //
        //    Input, double PTEST[DIM_NUM], the point whose nearest neighbor is sought.
        //
        //    Output, double *D_MIN, the distance between P and PSET(*,I_MIN).
        //
        //    Output, int POINTS_POINT_NEAR_NAIVE_ND, I_MIN, the index of the nearest
        //    point in PSET to P.
        //
    {
        int j;

        d_min = typeMethods.r8_huge ( );
        int p_min = -1;

        for ( j = 0; j < nset; j++ )
        {
            double d = 0.0;
            int i;
            for ( i = 0; i < dim_num; i++ )
            {
                d += ( ptest[i] - pset[i+j*dim_num] )
                     * ( ptest[i] - pset[i+j*dim_num] );
            }

            if (!(d < d_min))
            {
                continue;
            }

            d_min = d;
            p_min = j + 1;
        }

        d_min = Math.Sqrt ( d_min );

        return p_min;
    }

    public static int lrline(double xu, double yu, double xv1, double yv1, double xv2,
            double yv2, double dv)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    LRLINE determines where a point lies in relation to a directed line.
        //
        //  Discussion:
        //
        //    LRLINE determines whether a point is to the left of, right of,
        //    or on a directed line parallel to a line through given points.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    28 August 2003
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Barry Joe.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Barry Joe,
        //    GEOMPACK - a software package for the generation of meshes
        //    using geometric algorithms,
        //    Advances in Engineering Software,
        //    Volume 13, pages 325-331, 1991.
        //
        //  Parameters:
        //
        //    Input, double XU, YU, XV1, YV1, XV2, YV2, are vertex coordinates; the
        //    directed line is parallel to and at signed distance DV to the left of
        //    the directed line from (XV1,YV1) to (XV2,YV2); (XU,YU) is the vertex for
        //    which the position relative to the directed line is to be determined.
        //
        //    Input, double DV, the signed distance, positive for left.
        //
        //    Output, int LRLINE, is +1, 0, or -1 depending on whether (XU,YU) is
        //    to the right of, on, or left of the directed line.  LRLINE is 0 if
        //    the line degenerates to a point.
        //
    {
        const double tol = 0.0000001;
        int value = 0;

        double dx = xv2 - xv1;
        double dy = yv2 - yv1;
        double dxu = xu - xv1;
        double dyu = yu - yv1;

        double tolabs = tol * Math.Max(Math.Abs(dx),
            Math.Max(Math.Abs(dy),
                Math.Max(Math.Abs(dxu),
                    Math.Max(Math.Abs(dyu), Math.Abs(dv)))));

        double t = dy * dxu - dx * dyu + dv * Math.Sqrt(dx * dx + dy * dy);

        if (tolabs < t)
        {
            value = 1;
        }
        else if (-tolabs <= t)
        {
            value = 0;
        }
        else if (t < -tolabs)
        {
            value = -1;
        }

        return value;
    }

    public static double angle_rad_2d(double[] p1, double[] p2, double[] p3, int p1Index = 0, int p2Index = 0, int p3Index = 0)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ANGLE_RAD_2D returns the angle in radians swept out between two rays in 2D.
        //
        //  Discussion:
        //
        //      ANGLE_RAD_2D ( P1, P2, P3 ) + ANGLE_RAD_2D ( P3, P2, P1 ) = 2 * PI
        //
        //        P1
        //        /
        //       /
        //      /
        //     /
        //    P2--------->P3
        //
        //  Modified:
        //
        //    24 June 2005
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double P1[2], P2[2], P3[2], define the rays
        //    P1 - P2 and P3 - P2 which define the angle.
        //
        //    Output, double ANGLE_RAD_3D, the angle between the two rays,
        //    in radians.  This value will always be between 0 and 2*PI.  If either ray has
        //    zero length, then the angle is returned as zero.
        //
    {
        const int DIM_NUM = 2;

        double[] p = new double[DIM_NUM];
        double value = 0;

        p[0] = (p3[p3Index + 0] - p2[p2Index + 0]) * (p1[p1Index + 0] - p2[p2Index + 0])
               + (p3[p3Index + 1] - p2[p2Index + 1]) * (p1[p1Index + 1] - p2[p2Index + 1]);


        p[1] = (p3[p3Index + 0] - p2[p2Index + 0]) * (p1[p1Index + 1] - p2[p2Index + 1])
               - (p3[p3Index + 1] - p2[p2Index + 1]) * (p1[p1Index + 0] - p2[p2Index + 0]);

        switch (p[0])
        {
            case 0.0 when p[1] == 0.0:
                value = 0.0;
                return value;
        }

        value = Math.Atan2(p[1], p[0]);

        switch (value)
        {
            case < 0.0:
                value += 2.0 * Math.PI;
                break;
        }

        return value;
    }

    public static double arc_cosine(double c)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    ARC_COSINE computes the arc cosine function, with argument truncation.
        //
        //  Discussion:
        //
        //    If you call your system ACOS routine with an input argument that is
        //    outside the range [-1.0, 1.0 ], you may get an unpleasant surprise.
        //    This routine truncates arguments outside the range.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 June 2002
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double C, the argument, the cosine of an angle.
        //
        //    Output, double ARC_COSINE, an angle whose cosine is C.
        //
    {
        double value = c switch
        {
            <= -1.0 => Math.PI,
            >= 1.0 => 0.0,
            _ => Math.Acos(c)
        };

        return value;
    }

    public static double haversine ( double a )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    HAVERSINE computes the haversine of an angle.
        //
        //  Discussion:
        //
        //    haversine(A) = ( 1 - cos ( A ) ) / 2
        //
        //    The haversine is useful in spherical trigonometry.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    02 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double A, the angle.
        //
        //    Output, double HAVERSINE, the haversine of the angle.
        //
    {
        return ( 1.0 - Math.Cos ( a ) ) / 2.0;
    }

    public static double inverse_error ( int n, double[] a, double[] b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INVERSE_ERROR determines the error in an inverse matrix.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    30 October 2013
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int N, the order of the matrix.
        //
        //    Input, double A[N*N], the matrix.
        //
        //    Input, double B[N*N], the inverse.
        //
        //    Output, double ERROR_FROBENIUS, the Frobenius norm
        //    of (A*B-I) + (B*A-I).
        //
    {
        int j;

        double[] c = typeMethods.r8mat_mm_new ( n, n, n, a, b );

        for ( j = 0; j < n; j++ )
        {
            c[j+j*n] -= 1.0;
        }

        double value = typeMethods.r8mat_norm_fro ( n, n, c );

        c = typeMethods.r8mat_mm_new ( n, n, n, b, a );

        for ( j = 0; j < n; j++ )
        {
            c[j+j*n] -= 1.0;
        }

        value += typeMethods.r8mat_norm_fro ( n, n, c );

        return value;
    }
        
    public static int inverse_mod_n ( int b, int n )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    INVERSE_MOD_N computes the inverse of B mod N.
        //
        //  Discussion:
        //
        //    If 
        //
        //      Y = inverse_mod_n ( B, N )
        //
        //    then
        //
        //      mod ( B * Y, N ) = 1
        //
        //    The value Y will exist if and only if B and N are relatively prime.
        //
        //  Examples:
        //
        //    B  N  Y
        //
        //    1  2  1
        //
        //    1  3  1
        //    2  3  2
        //
        //    1  4  1
        //    2  4  0
        //    3  4  3
        //
        //    1  5  1
        //    2  5  3
        //    3  5  2
        //    4  5  4
        //
        //    1  6  1
        //    2  6  0
        //    3  6  0
        //    4  6  0
        //    5  6  5
        //
        //    1  7  1
        //    2  7  4
        //    3  7  5
        //    4  7  2
        //    5  7  3
        //    6  7  6
        //
        //    1  8  1
        //    2  8  0
        //    3  8  3
        //    4  8  0
        //    5  8  5
        //    6  8  0
        //    7  8  7
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    03 November 2009
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int B, the number whose inverse mod N is desired.
        //    B should be positive.  Normally, B < N, but this is not required.
        //
        //    Input, int N, the number with respect to which the
        //    modulus is computed.  N should be positive.
        //
        //    Output, int INVERSE_MOD_N, the inverse of B mod N, or 0 if there
        //    is not inverse for B mode N.
        //
    {
        int y;

        int n0 = n;
        int b0 = b;
        int t0 = 0;
        int t = 1;

        int q = n / b;
        int r = n - q * b;

        while ( 0 < r )
        {
            int temp = t0 - q * t;

            switch (temp)
            {
                case >= 0:
                    temp %= n;
                    break;
            }

            temp = temp switch
            {
                < 0 => n - -temp % n,
                _ => temp
            };

            t0 = t;
            t = temp;
            n0 = b0;
            b0 = r;
            q = n0 / b0;
            r = n0 - q * b0;
        }

        if ( b0 != 1 )
        {
            y = 0;
            return y;
        }

        y = t % n;

        return y;
    }

    public static void polar_to_xy ( double r, double t, ref double[] xy )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POLAR_TO_XY converts polar coordinates to XY coordinates.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 January 2007
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double R, T, the radius and angle (in radians).
        //
        //    Output, double XY[2], the Cartesian coordinates.
        //
    {
        xy[0] = r * Math.Cos ( t );
        xy[1] = r * Math.Sin ( t );

    }
    public static int power_mod(int a, int n, int m)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    POWER_MOD computes mod ( A^N, M ).
        //
        //  Discussion:
        //
        //    Some programming tricks are used to speed up the computation, and to
        //    allow computations in which A^N is much too large to store in a
        //    real word.
        //
        //    First, for efficiency, the power A^N is computed by determining
        //    the binary expansion of N, then computing A, A^2, A^4, and so on
        //    by repeated squaring, and multiplying only those factors that
        //    contribute to A^N.
        //
        //    Secondly, the intermediate products are immediately "mod'ed", which
        //    keeps them small.
        //
        //    For instance, to compute mod ( A^13, 11 ), we essentially compute
        //
        //       13 = 1 + 4 + 8
        //
        //       A**13 = A * A^4 * A^8
        //
        //       mod ( A^13, 11 ) = mod ( A, 11 ) * mod ( A^4, 11 ) * mod ( A^8, 11 ).
        //
        //    Fermat's little theorem says that if P is prime, and A is not divisible
        //    by P, then ( A^(P-1) - 1 ) is divisible by P.
        //
        //  Example:
        //
        //     A  N  M  X
        //
        //    13  0 31  1
        //    13  1 31 13
        //    13  2 31 14
        //    13  3 31 27
        //    13  4 31 10
        //    13  5 31  6
        //    13  6 31 16
        //    13  7 31 22
        //    13  8 31  7 
        //    13  9 31 29
        //    13 10 31  5
        //    13 11 31  3
        //    13 12 31  8
        //    13 13 31 11
        //    13 14 31 19
        //    13 15 31 30
        //    13 16 31 18
        //    13 17 31 17
        //    13 18 31  4
        //    13 19 31 21
        //    13 20 31 25
        //    13 21 31 15
        //    13 22 31  9
        //    13 23 31 24
        //    13 24 31  2
        //    13 25 31 26
        //    13 26 31 28
        //    13 27 31 23
        //    13 28 31 20
        //    13 29 31 12
        //    13 30 31  1
        //    13 31 31 13
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    16 November 2004
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int A, the base of the expression to be tested.
        //    A should be nonnegative.
        //
        //    Input, int N, the power to which the base is raised.
        //    N should be nonnegative.
        //
        //    Input, int M, the divisor against which the expression is tested.
        //    M should be positive.
        //
        //    Output, int POWER_MOD, the remainder when A^N is divided by M.
        //
    {
        switch (a)
        {
            case < 0:
                return -1;
        }

        switch (m)
        {
            case <= 0:
                return -1;
        }

        switch (n)
        {
            case < 0:
                return -1;
        }

        //
        //  A_SQUARE2 contains the successive squares of A.
        //
        long a_square2 = a;
        long m2 = m;
        long x2 = 1;

        while (0 < n)
        {
            int d = n % 2;

            x2 = d switch
            {
                1 => x2 * a_square2 % m2,
                _ => x2
            };

            a_square2 = a_square2 * a_square2 % m2;
            n = (n - d) / 2;
        }

        //
        //  Ensure that 0 <= X.
        //
        while (x2 < 0)
        {
            x2 += m2;
        }

        int x = (int)x2;

        return x;
    }

    public static double pythag ( double a, double b )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    PYTHAG computes SQRT ( A * A + B * B ) carefully.
        //
        //  Discussion:
        //
        //    The formula
        //
        //      PYTHAG = sqrt ( A * A + B * B )
        //
        //    is reasonably accurate, but can fail if, for example, A^2 is larger
        //    than the machine overflow.  The formula can lose most of its accuracy
        //    if the sum of the squares is very large or very small.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    08 November 2012
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Smith, Boyle, Dongarra, Garbow, Ikebe,
        //    Klema, Moler.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    James Wilkinson, Christian Reinsch,
        //    Handbook for Automatic Computation,
        //    Volume II, Linear Algebra, Part 2,
        //    Springer, 1971,
        //    ISBN: 0387054146,
        //    LC: QA251.W67.
        //
        //    Brian Smith, James Boyle, Jack Dongarra, Burton Garbow,
        //    Yasuhiko Ikebe, Virginia Klema, Cleve Moler,
        //    Matrix Eigensystem Routines, EISPACK Guide,
        //    Lecture Notes in Computer Science, Volume 6,
        //    Springer Verlag, 1976,
        //    ISBN13: 978-3540075462,
        //    LC: QA193.M37.
        //
        //  Modified:
        //
        //    08 November 2012
        //
        //  Parameters:
        //
        //    Input, double A, B, the two legs of a right triangle.
        //
        //    Output, double PYTHAG, the length of the hypotenuse.
        //
    {
        double p = Math.Max ( Math.Abs ( a ), Math.Abs ( b ) );

        if (p == 0.0)
        {
            return p;
        }

        double r = Math.Min ( Math.Abs ( a ), Math.Abs ( b ) ) / p;
        r *= r;

        while ( true )
        {
            double t = 4.0 + r;

            if ( Math.Abs(t - 4.0) <= double.Epsilon )
            {
                break;
            }

            double s = r / t;
            double u = 1.0 + 2.0 * s;
            p = u * p;
            r = s / u * ( s / u ) * r;
        }
        return p;
    }
        
    public static double radians_to_degrees ( double angle )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RADIANS_TO_DEGREES converts an angle from radians to degrees.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    27 August 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double ANGLE, an angle in radians.
        //
        //    Output, double RADIANS_TO_DEGREES, the equivalent angle in degrees.
        //
    {
        double value = angle / Math.PI * 180.0;

        return value;
    }

    public static void radians_to_dms ( double radians, ref int degrees, ref int minutes, ref int seconds )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RADIANS_TO_DMS converts an angle from radians to degrees/minutes/seconds.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    01 September 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, double RADIANS, the angle in radians.
        //
        //    Output, int *DEGREES, *MINUTES, *SECONDS, the equivalent angle in
        //    degrees, minutes, and seconds.
        //
    {
        double angle = 180.0 * Math.Abs ( radians ) / Math.PI;

        degrees = ( int ) angle;
        angle = ( angle - degrees ) * 60.0;
        minutes = ( int ) angle;
        angle = ( angle - minutes ) * 60.0;
        seconds = ( int ) angle;

        switch (radians)
        {
            case < 0.0:
                degrees = - degrees;
                minutes = - minutes;
                seconds = - seconds;
                break;
        }
    }
        
        
    public static void sincos(double input, ref double outsin, ref double outcos)
    {
        outsin = Math.Sin(input);
        outcos = Math.Cos(input);
    }
        
    public static double[] getExampleDoubleData()
    {
        string[] lines;

        int n = 2;
        int m = 100;
        double[] a = new double [m*n];
            
        try
        {
            lines = File.ReadAllLines("points_100.txt");
        }
        catch (Exception)
        {
            Console.WriteLine("Could not open points_100.txt");
            throw;
        }
            
        for (int i = 1; i <= m; i++)
        {
            r8vec values = typeMethods.s_to_r8vec(lines[i - 1], n);
            for (int j = 1; j <= n; j++)
            {
                a[i - 1 + (j - 1) * m] = values.rvec[j - 1];
            }
        }

        return a;
    }

    private const int GammaN = 10;

    private const double GammaR = 10.900511;

    private static readonly double[] GammaDk =
    {
        2.48574089138753565546e-5,
        1.05142378581721974210,
        -3.45687097222016235469,
        4.51227709466894823700,
        -2.98285225323576655721,
        1.05639711577126713077,
        -1.95428773191645869583e-1,
        1.70970543404441224307e-2,
        -5.71926117404305781283e-4,
        4.63399473359905636708e-6,
        -2.71994908488607703910e-9
    };

    public const double TwoSqrtEOverPi = 1.8603827342052657173362492472666631120594218414085755;
    public const double LogTwoSqrtEOverPi = 0.6207822376352452223455184457816472122518527279025978;
    public const double LnPi = 1.1447298858494001741434273513530587116472948129153d;
        
    public static double LogGamma(double z)
    {
        switch (z)
        {
            case < 0.5:
            {
                double s = GammaDk[0];
                for (int i = 1; i <= GammaN; i++)
                {
                    s += GammaDk[i] / (i - z);
                }

                return LnPi
                       - Math.Log(Math.Sin(Math.PI * z))
                       - Math.Log(s)
                       - LogTwoSqrtEOverPi
                       - (0.5 - z) * Math.Log((0.5 - z + GammaR) / Math.E);
            }
            default:
            {
                double s = GammaDk[0];
                for (int i = 1; i <= GammaN; i++)
                {
                    s += GammaDk[i] / (z + i - 1.0);
                }

                return Math.Log(s)
                       + LogTwoSqrtEOverPi
                       + (z - 0.5) * Math.Log((z - 0.5 + GammaR) / Math.E);
            }
        }
    }
        
    public static int fibonacci ( int k )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    FIBONACCI returns the Fibonacci number of given index.
        //
        //  Example:
        //
        //    K   Fibonacci
        //
        //    0   0
        //    1   1
        //    2   1
        //    3   2
        //    4   3
        //    5   5
        //    6   8
        //    7  13
        //    8  21
        //    9  34
        //   10  55
        //   11  89
        //   12 144
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license. 
        //
        //  Modified:
        //
        //    19 April 2003
        //
        //  Author:
        //
        //    John Burkardt
        //
        //  Parameters:
        //
        //    Input, int K, the index of the Fibonacci number to be used.
        //    K must be at least 1.
        //
        //    Output, int FIBONACCI, the value of the K-th Fibonacci number.
        //
    {
        int a;
        int kk;

        switch (k)
        {
            case < 0:
                a = - typeMethods.i4_huge ( );
                return a;
            case 0:
                a = 0;
                return a;
            case 1:
                a = 1;
                return a;
        }

        int c = 0;
        int b = 0;
        a = 1;

        for ( kk = 2; kk <= k; kk++ )
        {
            c = b;
            b = a;
            a = c + b;
        }
        return a;
    }
        
    public static double Gamma(double z)
    {
        switch (z)
        {
            case < 0.5:
            {
                double s = GammaDk[0];
                for (int i = 1; i <= GammaN; i++)
                {
                    s += GammaDk[i]/(i - z);
                }

                return Math.PI/(Math.Sin(Math.PI*z)
                                *s
                                *TwoSqrtEOverPi
                                *Math.Pow((0.5 - z + GammaR)/Math.E, 0.5 - z));
            }
            default:
            {
                double s = GammaDk[0];
                for (int i = 1; i <= GammaN; i++)
                {
                    s += GammaDk[i]/(z + i - 1.0);
                }

                return s*TwoSqrtEOverPi*Math.Pow((z - 0.5 + GammaR)/Math.E, z - 0.5);
            }
        }
    }

    public static double GammaOLD
    (
        double x    // We require x > 0
    )
    {
        switch (x)
        {
            case <= 0.0:
            {
                string msg = string.Format("Invalid input argument {0}. Argument must be positive.", x);
                throw new ArgumentOutOfRangeException(msg);
            }
            default:
            {
                // Split the function domain into three intervals:
                // (0, 0.001), [0.001, 12), and (12, infinity)

                ///////////////////////////////////////////////////////////////////////////
                // First interval: (0, 0.001)
                //
                // For small x, 1/Gamma(x) has power series x + gamma x^2  - ...
                // So in this range, 1/Gamma(x) = x + gamma x^2 with error on the order of x^3.
                // The relative error over this interval is less than 6e-7.

                const double gamma = 0.577215664901532860606512090; // Euler's gamma constant

                switch (x)
                {
                    case < 0.001:
                        return 1.0/(x*(1.0 + gamma*x));
                    ///////////////////////////////////////////////////////////////////////////
                    // Second interval: [0.001, 12)
                    case < 12.0:
                    {
                        // The algorithm directly approximates gamma over (1,2) and uses
                        // reduction identities to reduce other arguments to this interval.
		
                        double y = x;
                        int n = 0;
                        bool arg_was_less_than_one = y < 1.0;

                        switch (arg_was_less_than_one)
                        {
                            // Add or subtract integers as necessary to bring y into (1,2)
                            // Will correct for this below
                            case true:
                                y += 1.0;
                                break;
                            default:
                                n = (int) Math.Floor(y) - 1;  // will use n later
                                y -= n;
                                break;
                        }

                        // numerator coefficients for approximation over the interval (1,2)
                        double[] p =
                        {
                            -1.71618513886549492533811E+0,
                            2.47656508055759199108314E+1,
                            -3.79804256470945635097577E+2,
                            6.29331155312818442661052E+2,
                            8.66966202790413211295064E+2,
                            -3.14512729688483675254357E+4,
                            -3.61444134186911729807069E+4,
                            6.64561438202405440627855E+4
                        };

                        // denominator coefficients for approximation over the interval (1,2)
                        double[] q =
                        {
                            -3.08402300119738975254353E+1,
                            3.15350626979604161529144E+2,
                            -1.01515636749021914166146E+3,
                            -3.10777167157231109440444E+3,
                            2.25381184209801510330112E+4,
                            4.75584627752788110767815E+3,
                            -1.34659959864969306392456E+5,
                            -1.15132259675553483497211E+5
                        };

                        double num = 0.0;
                        double den = 1.0;
                        int i;

                        double z = y - 1;
                        for (i = 0; i < 8; i++)
                        {
                            num = (num + p[i])*z;
                            den = den*z + q[i];
                        }
                        double result = num/den + 1.0;

                        switch (arg_was_less_than_one)
                        {
                            // Apply correction if argument was not initially in (1,2)
                            case true:
                                // Use identity gamma(z) = gamma(z+1)/z
                                // The variable "result" now holds gamma of the original y + 1
                                // Thus we use y-1 to get back the orginal y.
                                result /= y-1.0;
                                break;
                            default:
                            {
                                // Use the identity gamma(z+n) = z*(z+1)* ... *(z+n-1)*gamma(z)
                                for (i = 0; i < n; i++)
                                    result *= y++;
                                break;
                            }
                        }

                        return result;
                    }
                    ///////////////////////////////////////////////////////////////////////////
                    // Third interval: [12, infinity)
                    case > 171.624:
                        // Correct answer too large to display. 
                        return double.PositiveInfinity;
                    default:
                        return Math.Exp(LogGammaOLD(x));
                }
            }
        }
    }

    public static double LogGammaOLD
    (
        double x    // x must be positive
    )
    {
        switch (x)
        {
            case <= 0.0:
            {
                string msg = string.Format("Invalid input argument {0}. Argument must be positive.", x);
                throw new ArgumentOutOfRangeException(msg);
            }
            case < 12.0:
                return Math.Log(Math.Abs(GammaOLD(x)));
        }

        // Abramowitz and Stegun 6.1.41
        // Asymptotic series should be good to at least 11 or 12 figures
        // For error analysis, see Whittiker and Watson
        // A Course in Modern Analysis (1927), page 252

        double[] c =
        {
            1.0/12.0,
            -1.0/360.0,
            1.0/1260.0,
            -1.0/1680.0,
            1.0/1188.0,
            -691.0/360360.0,
            1.0/156.0,
            -3617.0/122400.0
        };
        double z = 1.0/(x*x);
        double sum = c[7];
        for (int i=6; i >= 0; i--)
        {
            sum *= z;
            sum += c[i];
        }
        double series = sum/x;

        double halfLogTwoPi = 0.91893853320467274178032973640562;
        double logGamma = (x - 0.5)*Math.Log(x) - x + halfLogTwoPi + series;    
        return logGamma;
    }
        
    public static double Erf(double x)
    {
        // constants
        const double a1 = 0.254829592;
        const double a2 = -0.284496736;
        const double a3 = 1.421413741;
        const double a4 = -1.453152027;
        const double a5 = 1.061405429;
        const double p = 0.3275911;
 
        // Save the sign of x
        int sign = x switch
        {
            < 0 => -1,
            _ => 1
        };

        x = Math.Abs(x);
 
        // A&S formula 7.1.26
        double t = 1.0 / (1.0 + p*x);
        double y = 1.0 - ((((a5*t + a4)*t + a3)*t + a2)*t + a1)*t*Math.Exp(-x*x);
 
        return sign*y;
    }

    public static string[] splitStringByWhitespace(string line)
    {
        string[] tokens = line.Split(null as char[], StringSplitOptions.RemoveEmptyEntries);//line.Split(new char[0], StringSplitOptions.RemoveEmptyEntries);
        // line.Split(new char[] {' ', '\t'}, StringSplitOptions.RemoveEmptyEntries);
        return tokens;
    }
        
}