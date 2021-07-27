using System;
using Burkardt.LineNS;
using Burkardt.Uniform;

namespace Burkardt.Types
{
public static partial class typeMethods
{
public static void triangle_angles_2d ( double[] t, double[] angle )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_ANGLES_2D computes the angles of a triangle in 2D.
//
//  Discussion:
//
//    The law of cosines is used:
//
//      C * C = A * A + B * B - 2 * A * B * COS ( GAMMA )
//
//    where GAMMA is the angle opposite side C.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*3], the triangle vertices.
//
//    Output, double ANGLE[3], the angles opposite
//    sides P1-P2, P2-P3 and P3-P1, in radians.
//
{
double a;
double b;
double c;
double pi = 3.141592653589793;

a = Math.Sqrt ( Math.Pow ( t[0+1*2] - t[0+0*2], 2 ) 
+ Math.Pow ( t[1+1*2] - t[1+0*2], 2 ) );

b = Math.Sqrt ( Math.Pow ( t[0+2*2] - t[0+1*2], 2 ) 
+ Math.Pow ( t[1+2*2] - t[1+1*2], 2 ) );

c = Math.Sqrt ( Math.Pow ( t[0+0*2] - t[0+2*2], 2 ) 
+ Math.Pow ( t[1+0*2] - t[1+2*2], 2 ) );
//
//  Take care of a ridiculous special case.
//
if ( a == 0.0 && b == 0.0 && c == 0.0 )
{
angle[0] = 2.0 * pi / 3.0;
angle[1] = 2.0 * pi / 3.0;
angle[2] = 2.0 * pi / 3.0;
return;
}

if ( c == 0.0 || a == 0.0 )
{
angle[0] = pi;
}
else
{
angle[0] = r8_acos ( ( c * c + a * a - b * b ) / ( 2.0 * c * a ) );
}

if ( a == 0.0 || b == 0.0 )
{
angle[1] = pi;
}
else
{
angle[1] = r8_acos ( ( a * a + b * b - c * c ) / ( 2.0 * a * b ) );
}

if ( b == 0.0 || c == 0.0 )
{
angle[2] = pi;
}
else
{
angle[2] = r8_acos ( ( b * b + c * c - a * a ) / ( 2.0 * b * c ) );
}
}
public static void triangle_reference_sample(int n, ref int seed, ref double[] p)

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_REFERENCE_SAMPLE returns random points in the reference triangle.
//
//  Diagram:
//
//       3
//    s  |.
//    i  | .
//    d  |  .
//    e  |   .  side 2
//       |    .
//    3  |     .
//       |      .
//       1-------2
//
//         side 1
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, integer N, the number of points to sample.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double P[2*N], a random point in the triangle.
//
{
int DIM_NUM = 2;

double alpha;
double beta;
int j;
double r;

for (j = 0; j < n; j++)
{
r = UniformRNG.r8_uniform_01(ref seed);
//
//  Interpret R as a percentage of the triangle's area.
//
//  Imagine a line L, parallel to side 1, so that the area between
//  vertex 1 and line L is R percent of the full triangle's area.
//
//  The line L will intersect sides 2 and 3 at a fraction
//  ALPHA = SQRT ( R ) of the distance from vertex 1 to vertices 2 and 3.
//
alpha = Math.Sqrt(r);
//
//  Now choose, uniformly at random, a point on the line L.
//
beta = UniformRNG.r8_uniform_01(ref seed);

p[0 + j * 2] = (1.0 - beta) * alpha;
p[1 + j * 2] = beta * alpha;
}
}

public static void triangle_sample(double[] t, int n, ref int seed, ref double[] p, int pIndex = 0)

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_SAMPLE returns random points in a triangle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    06 December 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*3], the triangle vertices.
//
//    Input, integer N, the number of points to sample.
//
//    Input/output, int *SEED, a seed for the random number generator.
//
//    Output, double P[2*N], a random point in the triangle.
//
{
int DIM_NUM = 2;

double alpha;
double beta;
int j;
double r;
double[] p12 = new double[DIM_NUM];
double[] p13 = new double[DIM_NUM];

for (j = 0; j < n; j++)
{
r = UniformRNG.r8_uniform_01(ref seed);
//
//  Interpret R as a percentage of the triangle's area.
//
//  Imagine a line L, parallel to side 1, so that the area between
//  vertex 1 and line L is R percent of the full triangle's area.
//
//  The line L will intersect sides 2 and 3 at a fraction
//  ALPHA = SQRT ( R ) of the distance from vertex 1 to vertices 2 and 3.
//
alpha = Math.Sqrt(r);
//
//  Determine the coordinates of the points on sides 2 and 3 intersected
//  by line L.
//
p12[0] = (1.0 - alpha) * t[0 + 0 * 2] + alpha * t[0 + 1 * 2];
p12[1] = (1.0 - alpha) * t[1 + 0 * 2] + alpha * t[1 + 1 * 2];

p13[0] = (1.0 - alpha) * t[0 + 0 * 2] + alpha * t[0 + 2 * 2];
;
p13[1] = (1.0 - alpha) * t[1 + 0 * 2] + alpha * t[1 + 2 * 2];
;
//
//  Now choose, uniformly at random, a point on the line L.
//
beta = UniformRNG.r8_uniform_01(ref seed);

p[pIndex + (0 + j * 2)] = (1.0 - beta) * p12[0] + beta * p13[0];
p[pIndex + (1 + j * 2)] = (1.0 - beta) * p12[1] + beta * p13[1];
}
}

public static double[] triangle_angles_2d_new ( double[] t )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_ANGLES_2D_NEW computes the angles of a triangle in 2D.
//
//  Discussion:
//
//    The law of cosines is used:
//
//      C * C = A * A + B * B - 2 * A * B * COS ( GAMMA )
//
//    where GAMMA is the angle opposite side C.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    30 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*3], the triangle vertices.
//
//    Output, double ANGLE[3], the angles opposite
//    sides P1-P2, P2-P3 and P3-P1, in radians.
//
{
double[] angle;
double a;
double b;
double c;
double pi = 3.141592653589793;

angle = new double[3];

a = Math.Sqrt ( Math.Pow ( t[0+1*2] - t[0+0*2], 2 )
+ Math.Pow ( t[1+1*2] - t[1+0*2], 2 ) );

b = Math.Sqrt ( Math.Pow ( t[0+2*2] - t[0+1*2], 2 )
+ Math.Pow ( t[1+2*2] - t[1+1*2], 2 ) );

c = Math.Sqrt ( Math.Pow ( t[0+0*2] - t[0+2*2], 2 )
+ Math.Pow ( t[1+0*2] - t[1+2*2], 2 ) );
//
//  Take care of a ridiculous special case.
//
if ( a == 0.0 && b == 0.0 && c == 0.0 )
{
angle[0] = 2.0 * pi / 3.0;
angle[1] = 2.0 * pi / 3.0;
angle[2] = 2.0 * pi / 3.0;
return angle;
}

if ( c == 0.0 || a == 0.0 )
{
angle[0] = pi;
}
else
{
angle[0] = Helpers.arc_cosine ( ( c * c + a * a - b * b ) / ( 2.0 * c * a ) );
}

if ( a == 0.0 || b == 0.0 )
{
angle[1] = pi;
}
else
{
angle[1] = Helpers.arc_cosine ( ( a * a + b * b - c * c ) / ( 2.0 * a * b ) );
}

if ( b == 0.0 || c == 0.0 )
{
angle[2] = pi;
}
else
{
angle[2] = Helpers.arc_cosine ( ( b * b + c * c - a * a ) / ( 2.0 * b * c ) );
}

return angle;
}

public static double triangle_area_2d ( double[] t )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_AREA_2D computes the area of a triangle in 2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    03 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*3], the vertices of the triangle.
//
//    Output, double TRIANGLE_AREA_2D, the area of the triangle.  AREA will
//    be nonnegative.
//
{
double area = Math.Abs ( 0.5 * (
t[0+0*2] * ( t[1+2*2] - t[1+1*2] ) +
t[0+1*2] * ( t[1+0*2] - t[1+2*2] ) +
t[0+2*2] * ( t[1+1*2] - t[1+0*2] ) ) );

return area;
}

public static double[] triangle_centroid_2d ( double[] t )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_CENTROID_2D computes the centroid of a triangle in 2D.
//
//  Discussion:
//
//    The centroid of a triangle can also be considered the center
//    of gravity, assuming that the triangle is made of a thin uniform
//    sheet of massy material.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Adrian Bowyer, John Woodwark,
//    A Programmer's Geometry,
//    Butterworths, 1983.
//
//  Parameters:
//
//    Input, double T[2*3], the vertices of the triangle.
//
//    Output, double TRIANGLE_CENTROID_2D[2], the coordinates of the centroid of the triangle.
//
{
double[] centroid;

centroid = new double[2];

centroid[0] = ( t[0+0*2] + t[0+1*2] + t[0+2*2] ) / 3.0;
centroid[1] = ( t[1+0*2] + t[1+1*2] + t[1+2*2] ) / 3.0;

return centroid;
}

public static double[] triangle_circumcenter_2d(double[] t)

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_CIRCUMCENTER_2D computes the circumcenter of a triangle in 2D.
//
//  Discussion:
//
//    The circumcenter of a triangle is the center of the circumcircle, the
//    circle that passes through the three vertices of the triangle.
//
//    The circumcircle contains the triangle, but it is not necessarily the
//    smallest triangle to do so.
//
//    If all angles of the triangle are no greater than 90 degrees, then
//    the center of the circumscribed circle will lie inside the triangle.
//    Otherwise, the center will lie outside the triangle.
//
//    The circumcenter is the intersection of the perpendicular bisectors
//    of the sides of the triangle.
//
//    In geometry, the circumcenter of a triangle is often symbolized by "O".
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    09 February 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*3], the triangle vertices.
//
//    Output, double *TRIANGLE_CIRCUMCENTER_2D[2], the circumcenter of
//    the triangle.
//
{
int DIM_NUM = 2;

double asq;
double bot;
double[] center;
double csq;
double top1;
double top2;

center = new double[DIM_NUM];

asq = (t[0 + 1 * 2] - t[0 + 0 * 2]) * (t[0 + 1 * 2] - t[0 + 0 * 2])
+ (t[1 + 1 * 2] - t[1 + 0 * 2]) * (t[1 + 1 * 2] - t[1 + 0 * 2]);

csq = (t[0 + 2 * 2] - t[0 + 0 * 2]) * (t[0 + 2 * 2] - t[0 + 0 * 2])
+ (t[1 + 2 * 2] - t[1 + 0 * 2]) * (t[1 + 2 * 2] - t[1 + 0 * 2]);

top1 = (t[1 + 1 * 2] - t[1 + 0 * 2]) * csq - (t[1 + 2 * 2] - t[1 + 0 * 2]) * asq;
top2 = -(t[0 + 1 * 2] - t[0 + 0 * 2]) * csq + (t[0 + 2 * 2] - t[0 + 0 * 2]) * asq;

bot = (t[1 + 1 * 2] - t[1 + 0 * 2]) * (t[0 + 2 * 2] - t[0 + 0 * 2])
- (t[1 + 2 * 2] - t[1 + 0 * 2]) * (t[0 + 1 * 2] - t[0 + 0 * 2]);

center[0] = t[0 + 0 * 2] + 0.5 * top1 / bot;
center[1] = t[1 + 0 * 2] + 0.5 * top2 / bot;

return center;
}

public static void triangle_circumcircle_2d ( double[] t, ref double r, ref double[] pc )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_CIRCUMCIRCLE_2D computes the circumcircle of a triangle in 2D.
//
//  Discussion:
//
//    The circumcenter of a triangle is the center of the circumcircle, the
//    circle that passes through the three vertices of the triangle.
//
//    The circumcircle contains the triangle, but it is not necessarily the
//    smallest triangle to do so.
//
//    If all angles of the triangle are no greater than 90 degrees, then
//    the center of the circumscribed circle will lie inside the triangle.
//    Otherwise, the center will lie outside the triangle.
//
//    The circumcenter is the intersection of the perpendicular bisectors
//    of the sides of the triangle.
//
//    In geometry, the circumcenter of a triangle is often symbolized by "O".
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*3], the triangle vertices.
//
//    Output, double *R, PC[2], the circumradius, and the coordinates of the 
//    circumcenter of the triangle.
//
{
double a;
double b;
double bot;
double c;
double top1;
double top2;
//
//  Circumradius.
//
a = Math.Sqrt ( Math.Pow ( t[0+1*2] - t[0+0*2], 2 ) + Math.Pow ( t[1+1*2] - t[1+0*2], 2 ) );
b = Math.Sqrt ( Math.Pow ( t[0+2*2] - t[0+1*2], 2 ) + Math.Pow ( t[1+2*2] - t[1+1*2], 2 ) );
c = Math.Sqrt ( Math.Pow ( t[0+0*2] - t[0+2*2], 2 ) + Math.Pow ( t[1+0*2] - t[1+2*2], 2 ) );

bot = ( a + b + c ) * ( - a + b + c ) * (   a - b + c ) * (   a + b - c );

if ( bot <= 0.0 )
{
r = -1.0;
pc[0] = 0.0;
pc[1] = 0.0;
return;
}

r = a * b * c / Math.Sqrt ( bot );
//
//  Circumcenter.
//
top1 =  ( t[1+1*2] - t[1+0*2] ) * c * c - ( t[1+2*2] - t[1+0*2] ) * a * a;
top2 =  ( t[0+1*2] - t[0+0*2] ) * c * c - ( t[0+2*2] - t[0+0*2] ) * a * a;
bot  =  ( t[1+1*2] - t[1+0*2] ) * ( t[0+2*2] - t[0+0*2] )  
- ( t[1+2*2] - t[1+0*2] ) * ( t[0+1*2] - t[0+0*2] );

pc[0] = t[0+0*2] + 0.5 * top1 / bot;
pc[1] = t[1+0*2] - 0.5 * top2 / bot;

return;
}

public static bool triangle_contains_point_2d ( double[] t, double[] p )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_CONTAINS_POINT_2D finds if a point is inside a triangle in 2D.
//
//  Discussion:
//
//    The routine assumes that the vertices are given in counter clockwise
//    order.  If the triangle vertices are actually given in clockwise
//    order, this routine will behave as though the triangle contains
//    no points whatsoever!
//
//    The routine determines if P is "to the right of" each of the lines
//    that bound the triangle.  It does this by computing the cross product
//    of vectors from a vertex to its next vertex, and to P.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license.
//
//  Modified:
//
//    07 June 2006
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*3], the triangle vertices.
//    The vertices should be given in counter clockwise order.
//
//    Input, double P[2], the point to be checked.
//
//    Output, bool TRIANGLE_CONTAINS_POINT_2D, is TRUE if P is inside
//    the triangle or on its boundary.
//
{
int j;
int k;

for ( j = 0; j < 3; j++ )
{
k = ( j + 1 ) % 3;
if ( 0.0 < ( p[0] - t[0+j*2] ) * ( t[1+k*2] - t[1+j*2] )
- ( p[1] - t[1+j*2] ) * ( t[0+k*2] - t[0+j*2] ) )
{
return false;
}
}

return true;
}

public static double[] triangle_edge_length_2d ( double[] t )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_EDGE_LENGTH_2D returns edge lengths of a triangle in 2D.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    17 August 2009
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*3], the triangle vertices.
//
//    Output, double TRIANGLE_EDGE_LENGTH[3], the length of the edges.
//
{
double[] edge_length;
int j1;
int j2;

edge_length = new double[3];

for ( j1 = 0; j1 < 3; j1++ )
{
j2 = i4_wrap ( j1 + 1, 0, 2 );
edge_length[j1] = Math.Sqrt ( Math.Pow ( t[0+j2*2] - t[0+j1*2], 2 ) 
+ Math.Pow ( t[1+j2*2] - t[1+j1*2], 2 ) );
}

return edge_length;
}

public static void triangle_incircle_2d ( double[] t, ref double[] pc, ref double r )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_INCIRCLE_2D computes the inscribed circle of a triangle in 2D.
//
//  Discussion:
//
//    The inscribed circle of a triangle is the largest circle that can
//    be drawn inside the triangle.  It is tangent to all three sides,
//    and the lines from its center to the vertices bisect the angles
//    made by each vertex.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    03 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Adrian Bowyer, John Woodwark,
//    A Programmer's Geometry,
//    Butterworths, 1983.
//
//  Parameters:
//
//    Input, double T[2*3], the triangle vertices.
//
//    Output, double PC[2], *R, the center of the inscribed circle, and its radius.
//
{
double perim;
double s12;
double s23;
double s31;

s12 = Math.Sqrt ( Math.Pow ( t[0+1*2] - t[0+0*2], 2 ) 
+ Math.Pow ( t[1+1*2] - t[1+0*2], 2 ) );
s23 = Math.Sqrt ( Math.Pow ( t[0+2*2] - t[0+1*2], 2 ) 
+ Math.Pow ( t[1+2*2] - t[1+1*2], 2 ) );
s31 = Math.Sqrt ( Math.Pow ( t[0+0*2] - t[0+2*2], 2 ) 
+ Math.Pow ( t[1+0*2] - t[1+2*2], 2 ) );

perim = s12 + s23 + s31;

if ( perim == 0.0 )
{
r = 0.0;
pc[0] = t[0+0*2];
pc[1] = t[1+0*2];
}
else
{
pc[0] = ( s23 * t[0+0*2] + s31 * t[0+1*2] + s12 * t[0+2*2] ) / perim;
pc[1] = ( s23 * t[1+0*2] + s31 * t[1+1*2] + s12 * t[1+2*2] ) / perim;

r = 0.5 * Math.Sqrt (
( - s12 + s23 + s31 )
* ( + s12 - s23 + s31 )
* ( + s12 + s23 - s31 ) / perim );
}
}

public static int triangle_orientation_2d ( double[] t )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_ORIENTATION_2D determines the orientation of a triangle in 2D.
//
//  Discussion:
//
//    Three distinct non-colinear points in the plane define a circle.
//    If the points are visited in the order (x1,y1), (x2,y2), and then
//    (x3,y3), this motion defines a clockwise or counter clockwise
//    rotation along the circle.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    04 July 2005
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, double T[2*3], the triangle vertices.
//
//    Output, int TRIANGLE_ORIENTATION_2D, reports if the three points lie
//    clockwise on the circle that passes through them.  The possible
//    return values are:
//    0, the points are distinct, noncolinear, and lie counter clockwise
//    on their circle.
//    1, the points are distinct, noncolinear, and lie clockwise
//    on their circle.
//    2, the points are distinct and colinear.
//    3, at least two of the points are identical.
//
{
double det;
int value = 0;

if ( r8vec_eq ( 2, t, t, 0*2, 1*2 ) || 
r8vec_eq ( 2, t, t,1*2, 2*2 ) || 
r8vec_eq ( 2, t, t,2*2,0*2 ) )
{
value = 3;
return value;
}

det = ( t[0+0*2] - t[0+2*2] ) * ( t[1+1*2] - t[1+2*2] ) 
- ( t[0+1*2] - t[0+2*2] ) * ( t[1+0*2] - t[1+2*2] );

if ( det == 0.0 )
{
value = 2;
}
else if ( det < 0.0 )
{
value = 1;
}
else if ( 0.0 < det )
{
value = 0;
}
return value;
}

public static void triangle_orthocenter_2d ( double[] t, double[] p, ref bool flag )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_ORTHOCENTER_2D computes the orthocenter of a triangle in 2D.
//
//  Discussion:
//
//    The orthocenter is defined as the intersection of the three altitudes
//    of a triangle.
//
//    An altitude of a triangle is the line through a vertex of the triangle
//    and perpendicular to the opposite side.
//
//    In geometry, the orthocenter of a triangle is often symbolized by "H".
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Adrian Bowyer, John Woodwark,
//    A Programmer's Geometry,
//    Butterworths, 1983.
//
//  Parameters:
//
//    Input, double T[2*3], the triangle vertices.
//
//    Output, double P[2], the coordinates of the orthocenter of the triangle.
//
//    Output, bool *FLAG, is TRUE if the point could not be computed.
//
{
int ival = 0;
double[] p23;
double[] p31;
//
//  Determine a point P23 common to the line through P2 and P3 and
//  its perpendicular through P1.
//
p23 = Line.line_exp_perp_2d ( t, t, t, ref flag, p1Index:1*2, p2Index:2*2, p3Index:0*2 );

if ( flag )
{
p[0] = r8_huge ( );
p[1] = r8_huge ( );
return;
}
//
//  Determine a point P31 common to the line through P3 and P1 and
//  its perpendicular through P2.
//
p31 = Line.line_exp_perp_2d ( t, t, t, ref flag, p1Index:2*2, p2Index:0*2, p3Index:1*2 );
if ( flag )
{
p[0] = r8_huge ( );
p[1] = r8_huge ( );
}
//
//  Determine P, the intersection of the lines through P1 and P23, and
//  through P2 and P31.
//
Line.lines_exp_int_2d ( t, p23, t, p31, ref ival, ref p, p1Index:0*2, p3Index:1*2 );

if ( ival != 1 )
{
p[0] = r8_huge ( );
p[1] = r8_huge ( );
flag = true;
}
}

public static double triangle_quality_2d ( double[] t )

//****************************************************************************80
//
//  Purpose:
//
//    TRIANGLE_QUALITY_2D: "quality" of a triangle in 2D.
//
//  Discussion:
//
//    The quality of a triangle is 2 times the ratio of the radius of the inscribed
//    circle divided by that of the circumscribed circle.  An equilateral
//    triangle achieves the maximum possible quality of 1.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    30 July 2009
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Adrian Bowyer, John Woodwark,
//    A Programmer's Geometry,
//    Butterworths, 1983.
//
//  Parameters:
//
//    Input, double T[2*3], the triangle vertices.
//
//    Output, double TRIANGLE_QUALITY_2D, the quality of the triangle.
//
{
double a;
double b;
double c;
int i;
double value;
//
//  Compute the length of each side.
//
a = 0.0;
b = 0.0;
c = 0.0;

for ( i = 0; i < 2; i++ )
{
a = a + Math.Pow ( t[i+0*2] - t[i+1*2], 2 );
b = b + Math.Pow ( t[i+1*2] - t[i+2*2], 2 );
c = c + Math.Pow ( t[i+2*2] - t[i+0*2], 2 );
}
a = Math.Sqrt ( a );
b = Math.Sqrt ( b );
c = Math.Sqrt ( c );

if ( a * b * c == 0.0 )
{
value = 0.0;
}
else
{
value = ( - a + b + c ) * ( a - b + c ) * ( a + b - c ) 
/ ( a * b * c );
}
return value;
}
}
}
