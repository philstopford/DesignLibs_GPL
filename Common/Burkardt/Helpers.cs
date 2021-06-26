using System;
using System.IO;
using Burkardt.Types;

namespace Burkardt
{
    public static partial class Helpers
    {
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
            double d;
            int i;
            int j;
            int p_min;

            d_min = typeMethods.r8_huge ( );
            p_min = -1;

            for ( j = 0; j < nset; j++ )
            {
                d = 0.0;
                for ( i = 0; i < dim_num; i++ )
                {
                    d = d + ( ptest[i] - pset[i+j*dim_num] )
                        * ( ptest[i] - pset[i+j*dim_num] );
                }
                if ( d < d_min )
                {
                    d_min = d;
                    p_min = j + 1;
                }
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
            double dx;
            double dxu;
            double dy;
            double dyu;
            double t;
            double tol = 0.0000001;
            double tolabs;
            int value = 0;

            dx = xv2 - xv1;
            dy = yv2 - yv1;
            dxu = xu - xv1;
            dyu = yu - yv1;

            tolabs = tol * Math.Max(Math.Abs(dx),
                Math.Max(Math.Abs(dy),
                    Math.Max(Math.Abs(dxu),
                        Math.Max(Math.Abs(dyu), Math.Abs(dv)))));

            t = dy * dxu - dx * dyu + dv * Math.Sqrt(dx * dx + dy * dy);

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
            int DIM_NUM = 2;

            double[] p = new double[DIM_NUM];
            double value;

            p[0] = (p3[p3Index + 0] - p2[p2Index + 0]) * (p1[p1Index + 0] - p2[p2Index + 0])
                   + (p3[p3Index + 1] - p2[p2Index + 1]) * (p1[p1Index + 1] - p2[p2Index + 1]);


            p[1] = (p3[p3Index + 0] - p2[p2Index + 0]) * (p1[p1Index + 1] - p2[p2Index + 1])
                   - (p3[p3Index + 1] - p2[p2Index + 1]) * (p1[p1Index + 0] - p2[p2Index + 0]);

            if (p[0] == 0.0 && p[1] == 0.0)
            {
                value = 0.0;
                return value;
            }

            value = Math.Atan2(p[1], p[0]);

            if (value < 0.0)
            {
                value = value + 2.0 * Math.PI;
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
            double value;

            if (c <= -1.0)
            {
                value = Math.PI;
            }
            else if (1.0 <= c)
            {
                value = 0.0;
            }
            else
            {
                value = Math.Acos(c);
            }

            return value;
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
            double[] c;
            int j;
            double value;

            c = typeMethods.r8mat_mm_new ( n, n, n, a, b );

            for ( j = 0; j < n; j++ )
            {
                c[j+j*n] = c[j+j*n] - 1.0;
            }

            value = typeMethods.r8mat_norm_fro ( n, n, c );

            c = typeMethods.r8mat_mm_new ( n, n, n, b, a );

            for ( j = 0; j < n; j++ )
            {
                c[j+j*n] = c[j+j*n] - 1.0;
            }

            value = value + typeMethods.r8mat_norm_fro ( n, n, c );

            return value;
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
            double p;
            double r;
            double s;
            double t;
            double u;

            p = Math.Max ( Math.Abs ( a ), Math.Abs ( b ) );

            if ( p != 0.0 )
            {
                r = Math.Min ( Math.Abs ( a ), Math.Abs ( b ) ) / p;
                r = r * r;

                while ( true )
                {
                    t = 4.0 + r;

                    if ( t == 4.0 )
                    {
                        break;
                    }

                    s = r / t;
                    u = 1.0 + 2.0 * s;
                    p = u * p;
                    r = ( s / u ) * ( s / u ) * r;
                }
            }
            return p;
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
            catch (Exception e)
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
        
        public static double Gamma
        (
            double x    // We require x > 0
        )
        {
            if (x <= 0.0)
            {
                string msg = string.Format("Invalid input argument {0}. Argument must be positive.", x);
                throw new ArgumentOutOfRangeException(msg);
            }

            // Split the function domain into three intervals:
            // (0, 0.001), [0.001, 12), and (12, infinity)

            ///////////////////////////////////////////////////////////////////////////
            // First interval: (0, 0.001)
	        //
	        // For small x, 1/Gamma(x) has power series x + gamma x^2  - ...
	        // So in this range, 1/Gamma(x) = x + gamma x^2 with error on the order of x^3.
	        // The relative error over this interval is less than 6e-7.

	        const double gamma = 0.577215664901532860606512090; // Euler's gamma constant

            if (x < 0.001)
                return 1.0/(x*(1.0 + gamma*x));

            ///////////////////////////////////////////////////////////////////////////
            // Second interval: [0.001, 12)
    
	        if (x < 12.0)
            {
                // The algorithm directly approximates gamma over (1,2) and uses
                // reduction identities to reduce other arguments to this interval.
		
		        double y = x;
                int n = 0;
                bool arg_was_less_than_one = (y < 1.0);

                // Add or subtract integers as necessary to bring y into (1,2)
                // Will correct for this below
                if (arg_was_less_than_one)
                {
                    y += 1.0;
                }
                else
                {
                    n = (int) (Math.Floor(y)) - 1;  // will use n later
                    y -= n;
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

                // Apply correction if argument was not initially in (1,2)
                if (arg_was_less_than_one)
                {
                    // Use identity gamma(z) = gamma(z+1)/z
                    // The variable "result" now holds gamma of the original y + 1
                    // Thus we use y-1 to get back the orginal y.
                    result /= (y-1.0);
                }
                else
                {
                    // Use the identity gamma(z+n) = z*(z+1)* ... *(z+n-1)*gamma(z)
                    for (i = 0; i < n; i++)
                        result *= y++;
                }

		        return result;
            }

            ///////////////////////////////////////////////////////////////////////////
            // Third interval: [12, infinity)

            if (x > 171.624)
            {
		        // Correct answer too large to display. 
		        return double.PositiveInfinity;
            }

            return Math.Exp(LogGamma(x));
        }

        public static double LogGamma
        (
            double x    // x must be positive
        )
        {
	        if (x <= 0.0)
	        {
		        string msg = string.Format("Invalid input argument {0}. Argument must be positive.", x);
                throw new ArgumentOutOfRangeException(msg);
	        }

            if (x < 12.0)
            {
                return Math.Log(Math.Abs(Gamma(x)));
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
            double a1 = 0.254829592;
            double a2 = -0.284496736;
            double a3 = 1.421413741;
            double a4 = -1.453152027;
            double a5 = 1.061405429;
            double p = 0.3275911;
 
            // Save the sign of x
            int sign = 1;
            if (x < 0)
                sign = -1;
            x = Math.Abs(x);
 
            // A&S formula 7.1.26
            double t = 1.0 / (1.0 + p*x);
            double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*Math.Exp(-x*x);
 
            return sign*y;
        }

        public static string[] splitStringByWhitespace(string line)
        {
            string[] tokens = line.Split(new char[0], StringSplitOptions.RemoveEmptyEntries);
            // line.Split(new char[] {' ', '\t'}, StringSplitOptions.RemoveEmptyEntries);
            return tokens;
        }
        
    }
}