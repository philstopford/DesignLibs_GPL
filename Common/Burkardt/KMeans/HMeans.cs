using System;
using Burkardt.Types;
using Burkardt.Uniform;

namespace Burkardt.Means
{
public static class HMeans
{
public static void hmeans_01 ( int dim_num, int point_num, int cluster_num, int it_max, 
ref int it_num, double[] point, int[] cluster, double[] cluster_center, 
int[] cluster_population, double[] cluster_energy )

//****************************************************************************80
//
//  Purpose:
//
//   HMEANS_01 applies the H-Means algorithm.
//
//  Discussion:
//
//    The data for the H-Means problem is a set of N points X in
//    M-dimensions, and a desired number of clusters K.
//
//    The goal is to determine K points Z, called cluster centers, so that
//    if we associate each point X with its nearest Z value, we minimize
//    the standard deviation or cluster energy.  Writing CLUSTER(I) to
//    indicate the index of the nearest cluster center to point X(I), the
//    energy can be written as:
//
//      Energy = Sum ( 1 <= I <= N ) || X(I) - Z(CLUSTER(I)) ||^2
//
//    where
//
//      || X - Z ||^2 = Sum ( 1 <= J <= M ) ( X(J) - Z(J) )^2
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 October 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Wendy Martinez, Angel Martinez,
//    Computational Statistics Handbook with MATLAB,
//    pages 373-376,
//    Chapman and Hall / CRC, 2002.
//
//  Parameters:
//
//    Input, int DIM_NUM, the number of spatial dimensions.
//
//    Input, int POINT_NUM, the number of data points.
//
//    Input, int CLUSTER_NUM, the number of clusters.
//
//    Input, int IT_MAX, the maximum number of iterations.
//
//    Output, int &IT_NUM, the number of iterations taken.
//
//    Input, double POINT[DIM_NUM*POINT_NUM], the data points.
//
//    Input/output, int CLUSTER[POINT_NUM].  On input, the user 
//    may specify an initial cluster for each point, or leave all entrie of
//    CLUSTER set to 0.  On output, CLUSTER contains the index of the
//    cluster to which each data point belongs.
//
//    Input/output, double CLUSTER_CENTER[DIM_NUM*CLUSTER_NUM], the 
//    centers associated with the minimal energy clustering.
//
//    Output, int CLUSTER_POPULATION[CLUSTER_NUM],
//    the populuation of each cluster.
//
//    Output, double CLUSTER_ENERGY[CLUSTER_NUM], the energy
//    associated with each cluster.
//
{
double[] centroid;
int i;
int j;
int k;
int k2;
int missed;
double point_energy;
double point_energy_min;
int swap;
//
//  Data checks.
//
if ( cluster_num < 1 )
{
Console.WriteLine("");
Console.WriteLine("HMEANS_01 - Fatal error!");
Console.WriteLine("  CLUSTER_NUM < 1.");
return;
}

if ( dim_num < 1 )
{
Console.WriteLine("");
Console.WriteLine("HMEANS_01 - Fatal error!");
Console.WriteLine("  DIM_NUM < 1.");
return;
}

if ( point_num < 1 )
{
Console.WriteLine("");
Console.WriteLine("HMEANS_01 - Fatal error!");
Console.WriteLine("  POINT_NUM < 1.");
return;
}

if ( it_max < 0 )
{
Console.WriteLine("");
Console.WriteLine("HMEANS_01 - Fatal error!");
Console.WriteLine("  IT_MAX < 0.");
return;
}
//
//  On input, legal entries in CLUSTER are preserved, but
//  otherwise, each point is assigned to its nearest cluster.
//
for ( j = 0; j < point_num; j++ )
{
if ( cluster[j] < 0 || cluster_num <= cluster[j] )
{
point_energy_min = typeMethods.r8_huge ( );

for ( k = 0; k < cluster_num; k++ )
{
point_energy = 0.0;
for ( i = 0; i < dim_num; i++ )
{
point_energy = point_energy + 
Math.Pow ( point[i+j*dim_num] - cluster_center[i+k*dim_num], 2 );
}

if ( point_energy < point_energy_min )
{
point_energy_min = point_energy;
cluster[j] = k;
}
}
}
}
it_num = 0;

while ( it_num < it_max )
{
it_num = it_num + 1;
//
//  #1:
//  Assign each point to the cluster of its nearest center.
//
swap = 0;

for ( j = 0; j < point_num; j++ )
{
point_energy_min = typeMethods.r8_huge ( );
k = cluster[j];

for ( k2 = 0; k2 < cluster_num; k2++ )
{
point_energy = 0.0;
for ( i = 0; i < dim_num; i++ )
{
point_energy = point_energy + 
Math.Pow ( point[i+j*dim_num] - cluster_center[i+k2*dim_num], 2 );
}

if ( point_energy < point_energy_min )
{
point_energy_min = point_energy;
cluster[j] = k2;
}
}

if ( k != cluster[j] )
{
swap = swap + 1;
}
}
//
//  Terminate if no points were swapped.
//
if ( 1 < it_num )
{
if ( swap == 0 )
{
break;
}
}
//
//  #2:
//  Determine the total energy of the new clustering with current centroids.
//
typeMethods.r8vec_zero ( cluster_num, ref cluster_energy );

for ( j = 0; j < point_num; j++ )
{
k = cluster[j];

point_energy = 0.0;
for ( i = 0; i < dim_num; i++ )
{
point_energy = point_energy + 
Math.Pow ( point[i+j*dim_num] - cluster_center[i+k*dim_num], 2 );
}
cluster_energy[k] = cluster_energy[k] + point_energy;
}
//
//  #3:
//  Determine the centroids of the clusters.
//
centroid = typeMethods.r8vec_zero_new ( dim_num * cluster_num );
typeMethods.i4vec_zero ( cluster_num, ref cluster_population );

for ( j = 0; j < point_num; j++ )
{
k = cluster[j];
cluster_population[k] = cluster_population[k] + 1;
for ( i = 0; i < dim_num; i++ )
{
centroid[i+k*dim_num] = centroid[i+k*dim_num] + point[i+j*dim_num];
}
}
//
//  Now divide by the population to get the centroid.
//  But if a center has no population, pick a point at random.
//
missed = 0;

for ( k = 0; k < cluster_num; k++ )
{
if ( cluster_population[k] != 0 )
{
for ( i = 0; i < dim_num; i++ )
{
centroid[i+k*dim_num] = centroid[i+k*dim_num]
/ ( double ) ( cluster_population[k] );
}
}
else
{
for ( i = 0; i < dim_num; i++ )
{
centroid[i+k*dim_num] = point[i+missed*dim_num];
}
missed = missed + 1;
}
}

for ( k = 0; k < cluster_num; k++ )
{
for ( i = 0; i < dim_num; i++ )
{
cluster_center[i+k*dim_num] = centroid[i+k*dim_num];
}
}

//
//  #4:
//  Determine the total energy of the current clustering with new centroids.
//
typeMethods.r8vec_zero ( cluster_num, ref cluster_energy );

for ( j = 0; j < point_num; j++ )
{
k = cluster[j];

point_energy = 0.0;
for ( i = 0; i < dim_num; i++ )
{
point_energy = point_energy + 
Math.Pow ( point[i+j*dim_num] - cluster_center[i+k*dim_num], 2 );
}
cluster_energy[k] = cluster_energy[k] + point_energy;
}
}
}

public static void hmeans_02 ( int dim_num, int point_num, int cluster_num, int it_max, 
ref int it_num, double[] point, int[] cluster, double[] cluster_center, 
int[] cluster_population, double[] cluster_energy, ref int seed )

//****************************************************************************80
//
//  Purpose:
//
//   HMEANS_02 applies the H-Means algorithm.
//
//  Discussion:
//
//    This is a simple routine to group a set of points into K clusters,
//    each with a center point, in such a way that the total cluster 
//    energy is minimized.  The total cluster energy is the sum of the
//    squares of the distances of each point to the center of its cluster.
//
//    The algorithm begins with an initial estimate for the cluster centers:
//
//    1. The points are assigned to the nearest cluster centers.
//
//    2. The iteration exit ( 1 );s if the total energy has not changed 
//        significantly, or we have reached the maximum number of iterations.
//
//    3. Each cluster center is replaced by the centroid of the points
//       in the cluster.
//
//    4. Return to step 1.
//
//    The algorithm may fail to find the best solution.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 October 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the number of spatial dimensions.
//
//    Input, int POINT_NUM, the number of points.
//
//    Input, int CLUSTER_NUM, the number of clusters.
//
//    Input, int IT_MAX, the maximum number of iterations.
//
//    Output, int IT_NUM, the number of iterations taken.
//
//    Input, double POINT[DIM_NUM*POINT_NUM], the coordinates 
//    of the points.
//
//    Input/output, int CLUSTER[POINT_NUM].  On input, the user 
//    may specify an initial cluster for each point, or leave all entrie of
//    CLUSTER set to 0.  On output, CLUSTER contains the index of the
//    cluster to which each data point belongs.
//
//    Input/output, double CLUSTER_CENTER[DIM_NUM*CLUSTER_NUM],
//    the coordinates of the cluster centers.
//
//    Output, int CLUSTER_POPULATION[CLUSTER_NUM], the number of
//    points assigned to each cluster.
//
//    Output, double CLUSTER_ENERGY[CLUSTER_NUM], the energy of 
//    the clusters.
//
//    Input/output, int *SEED, a seed for the random
//    number generator.
//
{
bool debug = false;
int i;
int j;
int k;
int k2;
double point_energy;
double point_energy_min;
int swap;
//
//  Data checks.
//
if ( cluster_num < 1 )
{
Console.WriteLine("");
Console.WriteLine("HMEANS_02 - Fatal error!");
Console.WriteLine("  CLUSTER_NUM < 1.");
return;
}

if ( dim_num < 1 )
{
Console.WriteLine("");
Console.WriteLine("HMEANS_02 - Fatal error!");
Console.WriteLine("  DIM_NUM < 1.");
return;
}

if ( point_num < 1 )
{
Console.WriteLine("");
Console.WriteLine("HMEANS_02 - Fatal error!");
Console.WriteLine("  POINT_NUM < 1.");
return;
}

if ( it_max < 0 )
{
Console.WriteLine("");
Console.WriteLine("HMEANS_02 - Fatal error!");
Console.WriteLine("  IT_MAX < 0.");
return;
}
//
//  On input, legal entries in CLUSTER are preserved, but
//  otherwise, each point is assigned to its nearest cluster.
//
for ( j = 0; j < point_num; j++ )
{
if ( cluster[j] < 0 || cluster_num <= cluster[j] )
{
point_energy_min = typeMethods.r8_huge ( );
for ( k = 0; k < cluster_num; k++ )
{
point_energy = 0.0;
for ( i = 0; i < dim_num; i++ )
{
point_energy = point_energy + 
Math.Pow ( point[i+j*dim_num] - cluster_center[i+k*dim_num], 2 );
}
if ( point_energy < point_energy_min )
{
point_energy_min = point_energy;
cluster[j] = k;
}
}
}
}

it_num = 0;

for ( ; ; )
{
//
//  Given centers, assign points to nearest center.
//
typeMethods.i4vec_zero ( cluster_num, ref cluster_population );
typeMethods.r8vec_zero ( cluster_num, ref cluster_energy );

swap = 0;

for ( j = 0; j < point_num; j++ )
{
point_energy_min = typeMethods.r8_huge ( );
k = cluster[j];

for ( k2 = 0; k2 < cluster_num; k2++ )
{
point_energy = 0.0;
for ( i = 0; i < dim_num; i++ )
{
point_energy = point_energy + 
Math.Pow ( point[i+j*dim_num] - cluster_center[i+k2*dim_num], 2 );
}

if ( point_energy < point_energy_min )
{
point_energy_min = point_energy;
cluster[j] = k2;
}
}
if ( k != cluster[j] )
{
swap = swap + 1;
}
k = cluster[j];
cluster_energy[k] = cluster_energy[k] + point_energy_min;
cluster_population[k] = cluster_population[k] + 1;
}

if ( debug )
{
Console.WriteLine("  " + it_num.ToString().PadLeft(3)
+ "  " + typeMethods.r8vec_sum ( cluster_num, cluster_energy ).ToString().PadLeft(14) + "");
}

if ( 0 < it_num )
{
if ( swap == 0 )
{
break;
}
}

if ( it_max <= it_num )
{
break;
}

it_num = it_num + 1;
//
//  Given points in cluster, replace center by centroid.
//
typeMethods.r8vec_zero ( dim_num * cluster_num, ref cluster_center );

for ( j = 0; j < point_num; j++ )
{
k = cluster[j];
for ( i = 0; i < dim_num; i++ )
{
cluster_center[i+k*dim_num] = cluster_center[i+k*dim_num] 
+ point[i+j*dim_num];
}
}

for ( k = 0; k < cluster_num; k++ )
{
if ( cluster_population[k] != 0 )
{
for ( i = 0; i < dim_num; i++ )
{
cluster_center[i+k*dim_num] = cluster_center[i+k*dim_num] / 
( double ) ( cluster_population[k] );
}
}
else
{
j = UniformRNG.i4_uniform ( 0, point_num - 1, ref seed );
for ( i = 0; i < dim_num; i++ )
{
cluster_center[i+k*dim_num] = point[i+j*dim_num];
}
}
}
}
//
//  Compute the energy based on the final value of the cluster centers.
//
typeMethods.r8vec_zero ( cluster_num, ref cluster_energy );

for ( j = 0; j < point_num; j++ )
{
k = cluster[j];

point_energy = 0.0;
for ( i = 0; i < dim_num; i++ )
{
point_energy = point_energy + 
Math.Pow ( point[i+j*dim_num] - cluster_center[i+k*dim_num], 2 );
}
cluster_energy[k] = cluster_energy[k] + point_energy;
}
}

public static void hmeans_w_01 ( int dim_num, int point_num, int cluster_num, int it_max, 
ref int it_num, double[] point, double[] weight, int[] cluster, 
double[] cluster_center, int[] cluster_population, double[] cluster_energy )

//****************************************************************************80
//
//  Purpose:
//
//   HMEANS_W_01 applies the weighted H-Means algorithm. 
//
//  Discussion:
//
//    The input data for the weight H-Means problem includes:
//    * a set of N data points X in M dimensions, 
//    * a set of N nonnegative weights W,
//    * a desired number of clusters K.
//    * an initial set of cluster centers Z,
//    * an (optional) initial set of cluster assignments.
//
//    The goal is to determine K points Z, called cluster centers, and
//    to assign each point X(I) to some cluster Z(J), so that we minimize
//    the weighted standard deviation of the distance of each data point
//    to the center of its cluster.  Writing J = CLUSTER(I) to
//    indicate the index of the nearest cluster center Z(J) to the 
//    point X(I), the quantity we are trying to minimize is the sum
//    of the weighted cluster energies E(J), where:
//
//      E(J) = Sum ( 1 <= I <= N ) W(I) * || X(I) - Z(J) ||^2
//
//    Here, we assume that we are using the Euclidean norm, so that
//    
//      || X(I) - Z(J) ||^2 = Sum ( 1 <= K <= M )
//        ( X(I)(K) - Z(J)(K) )^2
//
//    In this notation, X(I)(K) is the K-th spatial component of the
//    I-th data point.
//
//    Note that this routine should give the same results as HMEANS_01
//    in any case in which all the entries of the WEIGHT vector are equal.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    08 October 2011
//
//  Author:
//
//    John Burkardt
//
//  Reference:
//
//    Wendy Martinez, Angel Martinez,
//    Computational Statistics Handbook with MATLAB,
//    pages 373-376,
//    Chapman and Hall / CRC, 2002.
//
//  Parameters:
//
//    Input, int DIM_NUM, the number of spatial dimensions.
//
//    Input, int POINT_NUM, the number of data points.
//
//    Input, int CLUSTER_NUM, the number of clusters.
//
//    Input, int IT_MAX, the maximum number of iterations.
//
//    Output, int &IT_NUM, the number of iterations taken.
//
//    Input, double POINT[DIM_NUM*POINT_NUM], the data points.
//
//    Input, double WEIGHT[POINT_NUM], the weights
//    assigned to the data points.  These must be nonnegative, and
//    at least one must be strictly positive.
//
//    Input/output, int CLUSTER[POINT_NUM].  On input, the user 
//    may specify an initial cluster for each point, or leave all entrie of
//    CLUSTER set to 0.  On output, CLUSTER contains the index of the
//    cluster to which each data point belongs.
//
//    Input/output, double CLUSTER_CENTER[DIM_NUM*CLUSTER_NUM], the
//    centers associated with the minimal energy clustering.
//
//    Output, double CLUSTER_ENERGY[CLUSTER_NUM], the energy
//    associated with each cluster.
//
{
double[] centroid;
double[] cluster_weight;
double energy;
int i;
int j;
int k;
int k2 = 0;
int missed;
double point_energy;
double point_energy_min;
int swap;
//
//  Data checks.
//
if ( cluster_num < 1 )
{
Console.WriteLine("");
Console.WriteLine("HMEANS_W_01 - Fatal error!");
Console.WriteLine("  CLUSTER_NUM < 1.");
return;
}

if ( dim_num < 1 )
{
Console.WriteLine("");
Console.WriteLine("HMEANS_W_01 - Fatal error!");
Console.WriteLine("  DIM_NUM < 1.");
return;
}

if ( point_num < 1 )
{
Console.WriteLine("");
Console.WriteLine("HMEANS_W_01 - Fatal error!");
Console.WriteLine("  POINT_NUM < 1.");
return;
}

if ( it_max < 0 )
{
Console.WriteLine("");
Console.WriteLine("HMEANS_W_01 - Fatal error!");
Console.WriteLine("  IT_MAX < 0.");
return;
}

if ( typeMethods.r8vec_any_negative ( point_num, weight ) )
{
Console.WriteLine("");
Console.WriteLine("HMEANS_W_01 - Fatal error!");
Console.WriteLine("  Some weight entry is negative.");
return;
}

if ( typeMethods.r8vec_all_nonpositive ( point_num, weight ) )
{
Console.WriteLine("");
Console.WriteLine("HMEANS_W_01 - Fatal error!");
Console.WriteLine("  No weight entry is positive.");
return;
}
//
//  On input, legal entries in CLUSTER are preserved, but
//  otherwise, each point is assigned to its nearest cluster.
//
for ( j = 0; j < point_num; j++ )
{
if ( cluster[j] < 0 || cluster_num <= cluster[j] )
{
point_energy_min = typeMethods.r8_huge ( );

for ( k = 0; k < cluster_num; k++ )
{
point_energy = 0.0;
for ( i = 0; i < dim_num; i++ )
{
point_energy = point_energy + 
Math.Pow ( point[i+j*dim_num] - cluster_center[i+k*dim_num], 2 );
}

if ( point_energy < point_energy_min )
{
point_energy_min = point_energy;
cluster[j] = k;
}
}
}
}
it_num = 0;

while ( it_num < it_max )
{
it_num = it_num + 1;
//
//  #1:
//  Reassign points to clusters:
//  Assign each point to the cluster whose center is nearest;
//  Count the number of points whose cluster assignment is changed.
//
swap = 0;

for ( j = 0; j < point_num; j++ )
{
point_energy_min = typeMethods.r8_huge ( );
k = cluster[j];

for ( k2 = 0; k2 < cluster_num; k2++ )
{
point_energy = 0.0;
for ( i = 0; i < dim_num; i++ )
{
point_energy = point_energy + 
Math.Pow ( point[i+j*dim_num] - cluster_center[i+k2*dim_num], 2 );
}

if ( point_energy < point_energy_min )
{
point_energy_min = point_energy;
cluster[j] = k2;
}
}

if ( k != cluster[j] )
{
swap = swap + 1;
}
}
//
//  If no point changed its cluster assignment, the algorithm can make no 
//  more improvements, so terminate.
//
if ( 1 < it_num )
{
if ( swap == 0 )
{
break;
}
}
//
//  Determine the current energy.
//
energy = 0.0;
for ( j = 0; j < point_num; j++ )
{
point_energy = 0.0;
for ( i = 0; i < dim_num; i++ )
{
point_energy = point_energy + 
Math.Pow ( point[i+j*dim_num] - cluster_center[i+k2*dim_num], 2 );
}
energy = energy + weight[j] * point_energy;
}
Console.WriteLine("  " + it_num.ToString().PadLeft(4)
+ "  " + energy.ToString().PadLeft(14) + "");
//
//  #2:
//  Determine the centroids of the clusters, and set the 
//  cluster center to the cluster centroid.
//
centroid = typeMethods.r8vec_zero_new ( dim_num * cluster_num );
cluster_weight = typeMethods.r8vec_zero_new ( cluster_num );

for ( j = 0; j < point_num; j++ )
{
k = cluster[j];
cluster_population[k] = cluster_population[k] + 1;
cluster_weight[k] = cluster_weight[k] + weight[j];
for ( i = 0; i < dim_num; i++ )
{
centroid[i+k*dim_num] = centroid[i+k*dim_num] 
+ weight[j] * point[i+j*dim_num];
}
}

missed = 0;

for ( k = 0; k < cluster_num; k++ )
{
if ( cluster_weight[k] != 0.0 )
{
for ( i = 0; i < dim_num; i++ )
{
centroid[i+k*dim_num] = centroid[i+k*dim_num] / cluster_weight[k];
}
}
else
{
for ( i = 0; i < dim_num; i++ )
{
centroid[i+k*dim_num] = point[i+missed*dim_num];
}
missed = missed + 1;
}
}

for ( k = 0; k < cluster_num; k++ )
{
for ( i = 0; i < dim_num; i++ )
{
cluster_center[i+k*dim_num] = centroid[i+k*dim_num];
}
}
}
//
//  Compute the energy based on the final value of the cluster centers.
//
typeMethods.r8vec_zero ( cluster_num, ref cluster_energy );

for ( j = 0; j < point_num; j++ )
{
k = cluster[j];

point_energy = 0.0;
for ( i = 0; i < dim_num; i++ )
{
point_energy = point_energy + 
Math.Pow ( point[i+j*dim_num] - cluster_center[i+k*dim_num], 2 );
}
cluster_energy[k] = cluster_energy[k] + weight[j] * point_energy;
}
}

public static void hmeans_w_02 ( int dim_num, int point_num, int cluster_num, int it_max, 
ref int it_num, double[] point, double[] weight, int[] cluster, 
double[] cluster_center, int[] cluster_population, double[] cluster_energy, 
ref int seed )

//****************************************************************************80
//
//  Purpose:
//
//   HMEANS_W_02 applies the weighted H-Means algorithm.
//
//  Discussion:
//
//    The input data for the weight H-Means problem includes:
//    * a set of N data points X in M dimensions, 
//    * a set of N nonnegative weights W,
//    * a desired number of clusters K.
//    * an initial set of cluster centers Z,
//    * an (optional) initial set of cluster assignments.
//
//    The goal is to determine K points Z, called cluster centers, and
//    to assign each point X(I) to some cluster Z(J), so that we minimize
//    the weighted standard deviation of the distance of each data point
//    to the center of its cluster.  Writing J = CLUSTER(I) to
//    indicate the index of the nearest cluster center Z(J) to the 
//    point X(I), the quantity we are trying to minimize is the sum
//    of the weighted cluster energies E(J), where:
//
//      E(J) = Sum ( 1 <= I <= N ) W(I) * || X(I) - Z(J) ||^2
//
//    Here, we assume that we are using the Euclidean norm, so that
//    
//      || X(I) - Z(J) ||^2 = Sum ( 1 <= K <= M )
//        ( X(I)(K) - Z(J)(K) )^2
//
//    In this notation, X(I)(K) is the K-th spatial component of the
//    I-th data point.
//
//    Note that this routine should give the same results as HMEANS_02
//    in any case in which all the entries of the WEIGHT vector are equal.
//
//  Licensing:
//
//    This code is distributed under the GNU LGPL license. 
//
//  Modified:
//
//    09 October 2011
//
//  Author:
//
//    John Burkardt
//
//  Parameters:
//
//    Input, int DIM_NUM, the number of spatial dimensions.
//
//    Input, int POINT_NUM, the number of points.
//
//    Input, int CLUSTER_NUM, the number of clusters.
//
//    Input, int IT_MAX, the maximum number of iterations.
//
//    Output, int &IT_NUM, the number of iterations taken.
//
//    Input, double POINT[DIM_NUM*POINT_NUM], the coordinates 
//    of the points.
//
//    Input, double WEIGHT[POINT_NUM], the weights
//    assigned to the data points.  These must be nonnegative, and
//    at least one must be strictly positive.
//
//    Input/output, int CLUSTER[POINT_NUM].  On input, the user 
//    may specify an initial cluster for each point, or leave all entrie of
//    CLUSTER set to 0.  On output, CLUSTER contains the index of the
//    cluster to which each data point belongs.
//
//    Input/output, double CLUSTER_CENTER[DIM_NUM*CLUSTER_NUM],
//    the coordinates of the cluster centers.
//
//    Output, int CLUSTER_POPULATION[CLUSTER_NUM], the number of
//    points assigned to each cluster.
//
//    Output, double CLUSTER_ENERGY[CLUSTER_NUM], the energy of 
//    the clusters.
//
//    Input/output, int *SEED, a seed for the random 
//    number generator.
//
{
double[] cluster_weight;
int i;
int j;
int k;
int k2;
double point_energy;
double point_energy_min;
int swap;
//
//  Data checks.
//
if ( cluster_num < 1 )
{
Console.WriteLine("");
Console.WriteLine("HMEANS_W_02 - Fatal error!");
Console.WriteLine("  CLUSTER_NUM < 1.");
return;
}

if ( dim_num < 1 )
{
Console.WriteLine("");
Console.WriteLine("HMEANS_W_02 - Fatal error!");
Console.WriteLine("  DIM_NUM < 1.");
return;
}

if ( point_num < 1 )
{
Console.WriteLine("");
Console.WriteLine("HMEANS_W_02 - Fatal error!");
Console.WriteLine("  POINT_NUM < 1.");
return;
}

if ( it_max < 0 )
{
Console.WriteLine("");
Console.WriteLine("HMEANS_W_02 - Fatal error!");
Console.WriteLine("  IT_MAX < 0.");
return;
}

if ( typeMethods.r8vec_any_negative ( point_num, weight ) )
{
Console.WriteLine("");
Console.WriteLine("HMEANS_W_02 - Fatal error!");
Console.WriteLine("  Some weight entry is negative.");
return;
}

if ( typeMethods.r8vec_all_nonpositive ( point_num, weight ) )
{
Console.WriteLine("");
Console.WriteLine("HMEANS_W_02 - Fatal error!");
Console.WriteLine("  No weight entry is positive.");
return;
}
//
//  On input, legal entries in CLUSTER are preserved, but
//  otherwise, each point is assigned to its nearest cluster.
//
for ( j = 0; j < point_num; j++ )
{
if ( cluster[j] < 0 || cluster_num <= cluster[j] )
{
point_energy_min = typeMethods.r8_huge ( );
for ( k = 0; k < cluster_num; k++ )
{
point_energy = 0.0;
for ( i = 0; i < dim_num; i++ )
{
point_energy = point_energy + 
Math.Pow ( point[i+j*dim_num] - cluster_center[i+k*dim_num], 2 );
}
if ( point_energy < point_energy_min )
{
point_energy_min = point_energy;
cluster[j] = k;
}
}
}
}
it_num = 0;

for ( ; ; )
{
//
//  Given centers, assign points to nearest center.
//
typeMethods.i4vec_zero ( cluster_num, ref cluster_population );
typeMethods.r8vec_zero ( cluster_num, ref cluster_energy );
cluster_weight = typeMethods.r8vec_zero_new ( cluster_num );

swap = 0;

for ( j = 0; j < point_num; j++ )
{
point_energy_min = typeMethods.r8_huge ( );
k = cluster[j];

for ( k2 = 0; k2 < cluster_num; k2++ )
{
point_energy = 0.0;
for ( i = 0; i < dim_num; i++ )
{
point_energy = point_energy + 
Math.Pow ( point[i+j*dim_num] - cluster_center[i+k2*dim_num], 2 );
}

if ( point_energy < point_energy_min )
{
point_energy_min = point_energy;
cluster[j] = k2;
}
}

if ( k != cluster[j] )
{
swap = swap + 1;
}
k = cluster[j];
cluster_energy[k] = cluster_energy[k] + weight[j] * point_energy_min;
cluster_population[k] = cluster_population[k] + 1;
cluster_weight[k] = cluster_weight[k] + weight[j];
}

if ( 0 < it_num )
{
if ( swap == 0 )
{
break;
}
}

if ( it_max <= it_num )
{
break;
}

it_num = it_num + 1;
//
//  Given points in cluster, replace center by weighted centroid.
//
typeMethods.r8vec_zero ( dim_num * cluster_num, ref cluster_center );

for ( j = 0; j < point_num; j++ )
{
k = cluster[j];
for ( i = 0; i < dim_num; i++ )
{
cluster_center[i+k*dim_num] = cluster_center[i+k*dim_num] 
+ weight[j] * point[i+j*dim_num];
}
}

for ( k = 0; k < cluster_num; k++ )
{
if ( cluster_weight[k] != 0.0 )
{
for ( i = 0; i < dim_num; i++ )
{
cluster_center[i+k*dim_num] = cluster_center[i+k*dim_num] / 
cluster_weight[k];
}
}
else
{
j = UniformRNG.i4_uniform ( 0, point_num - 1, ref seed );
for ( i = 0; i < dim_num; i++ )
{
cluster_center[i+k*dim_num] = point[i+j*dim_num];
}
}
}
}
//
//  Compute the energy based on the final value of the cluster centers.
//
typeMethods.r8vec_zero ( cluster_num, ref cluster_energy );

for ( j = 0; j < point_num; j++ )
{
k = cluster[j];

point_energy = 0.0;
for ( i = 0; i < dim_num; i++ )
{
point_energy = point_energy + 
Math.Pow ( point[i+j*dim_num] - cluster_center[i+k*dim_num], 2 );
}
cluster_energy[k] = cluster_energy[k] + weight[j] * point_energy;
}

return;
}
}
}