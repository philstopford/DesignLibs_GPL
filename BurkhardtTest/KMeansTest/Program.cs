using System;
using Burkardt.Means;
using Burkardt.SolveNS;
using Burkardt.Table;
using Burkardt.Types;

namespace KMeansTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for KMEANS_TEST.
        //
        //  Discussion:
        //
        //    KMEANS_TEST tests the KMEANS library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    13 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        Console.WriteLine("");
        Console.WriteLine("KMEANS_TEST");
        Console.WriteLine("  Test the KMEANS library.");

        test01();
        test02();
        test03();
        test04();
        test05();
        test06();
        test07();
        test08();
        test09();

        test10();
        test11();
        test12();
        test13();
        test14();
        test15();
        test16();

        Console.WriteLine("");
        Console.WriteLine("KMEANS_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tries out the HMEANS_01 routine.
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
    {
        const string center_filename = "test01_centers.txt";
        const string cluster_filename = "test01_clusters.txt";
        int it_num = 0;
        const string point_filename = "points_100.txt";

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Test the HMEANS_01 algorithm.");
        //
        //  Read the data points.
        //
        Console.WriteLine("");
        Console.WriteLine("  Data points will be read from \"" + point_filename + "\".");

        TableHeader h = typeMethods.r8mat_header_read(point_filename);
        int dim_num = h.m;
        int point_num = h.n;

        Console.WriteLine("");
        Console.WriteLine("  Point spatial dimension = " + dim_num + "");
        Console.WriteLine("  Number of points = " + point_num + "");

        double[] point = typeMethods.r8mat_data_read(point_filename, dim_num, point_num);
        //
        //  Clustering parameters.
        //
        const int cluster_num = 5;
        const int it_max = 20;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("  Number of iterations allowed is " + it_max + "");
        //
        //  Initialize the centers.
        //
        double[] cluster_center = Cluster.cluster_initialize_5(dim_num, point_num, cluster_num,
            point, ref seed);

        int[] cluster = typeMethods.i4vec_negone_new(point_num);
        double[] cluster_energy = new double[cluster_num];
        int[] cluster_population = new int[cluster_num];

        HMeans.hmeans_01(dim_num, point_num, cluster_num, it_max, ref it_num, point,
            cluster, cluster_center, cluster_population, cluster_energy);

        Console.WriteLine("");
        Console.WriteLine("  Number of iterations taken is " + it_num + "");

        double[] cluster_variance = Cluster.cluster_variance_compute(dim_num, point_num,
            cluster_num, point, cluster, cluster_center);

        Cluster.cluster_print_summary(point_num, cluster_num,
            cluster_population, cluster_energy, cluster_variance);

        typeMethods.r8mat_write(center_filename, dim_num, cluster_num, cluster_center);

        Console.WriteLine("");
        Console.WriteLine("  Cluster centers written to \"" + center_filename + "\".");

        typeMethods.i4mat_write(cluster_filename, 1, point_num, cluster);

        Console.WriteLine("  Cluster assignments written to \"" + cluster_filename + "\".");
    }

    private static void test02()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tries out the HMEANS_02 routine.
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
    {
        const string center_filename = "test02_centers.txt";
        const string cluster_filename = "test02_clusters.txt";
        int it_num = 0;
        const string point_filename = "points_100.txt";

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  Test the HMEANS_02 algorithm.");
        //
        //  Read the data points.
        //
        Console.WriteLine("");
        Console.WriteLine("  Data points will be read from \"" + point_filename + "\".");

        TableHeader h = typeMethods.r8mat_header_read(point_filename);
        int dim_num = h.m;
        int point_num = h.n;

        Console.WriteLine("");
        Console.WriteLine("  Point spatial dimension = " + dim_num + "");
        Console.WriteLine("  Number of points = " + point_num + "");

        double[] point = typeMethods.r8mat_data_read(point_filename, dim_num, point_num);
        //
        //  Clustering parameters.
        //
        const int cluster_num = 5;
        const int it_max = 20;
        int seed = 123456789;

        Console.WriteLine("");
        Console.WriteLine("  Number of iterations allowed is " + it_max + "");
        //
        //  Initialize the centers.
        //
        double[] cluster_center = Cluster.cluster_initialize_5(dim_num, point_num, cluster_num,
            point, ref seed);

        int[] cluster = typeMethods.i4vec_negone_new(point_num);
        double[] cluster_energy = new double[cluster_num];
        int[] cluster_population = new int[cluster_num];

        HMeans.hmeans_02(dim_num, point_num, cluster_num, it_max, ref it_num, point,
            cluster, cluster_center, cluster_population, cluster_energy, ref seed);

        Console.WriteLine("");
        Console.WriteLine("  Number of iterations taken is " + it_num + "");

        double[] cluster_variance = Cluster.cluster_variance_compute(dim_num, point_num,
            cluster_num, point, cluster, cluster_center);

        Cluster.cluster_print_summary(point_num, cluster_num,
            cluster_population, cluster_energy, cluster_variance);

        typeMethods.r8mat_write(center_filename, dim_num, cluster_num, cluster_center);

        Console.WriteLine("");
        Console.WriteLine("  Cluster centers written to \"" + center_filename + "\".");

        typeMethods.i4mat_write(cluster_filename, 1, point_num, cluster);

        Console.WriteLine("  Cluster assignments written to \"" + cluster_filename + "\".");
    }

    private static void test03()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tries out the KMEANS_01 routine.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const string center_filename = "test03_centers.txt";
        const string cluster_filename = "test03_clusters.txt";
        int it_num = 0;
        const string point_filename = "points_100.txt";

        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  Test the KMEANS_01 algorithm.");
        Console.WriteLine("  (Applied Statistics Algorithm #58)");
        //
        //  Read the data points.
        //
        Console.WriteLine("");
        Console.WriteLine("  Data points will be read from \"" + point_filename + "\".");

        TableHeader h = typeMethods.r8mat_header_read(point_filename);
        int dim_num = h.m;
        int point_num = h.n;

        Console.WriteLine("");
        Console.WriteLine("  Point spatial dimension = " + dim_num + "");
        Console.WriteLine("  Number of points = " + point_num + "");

        double[] point = typeMethods.r8mat_data_read(point_filename, dim_num, point_num);
        //
        //  Clustering parameters.
        //
        const int cluster_num = 5;
        const int it_max = 20;
        int seed = 123456789;

        int[] cluster = typeMethods.i4vec_negone_new(point_num);
        double[] cluster_energy = new double[cluster_num];
        int[] cluster_population = new int[cluster_num];

        Console.WriteLine("");
        Console.WriteLine("  Number of iterations allowed is " + it_max + "");
        //
        //  Initialize the centers.
        //
        double[] cluster_center = Cluster.cluster_initialize_5(dim_num, point_num, cluster_num, point,
            ref seed);

        KMeans.kmeans_01(dim_num, point_num, cluster_num, it_max, ref it_num, point,
            ref cluster, ref cluster_center, ref cluster_population, ref cluster_energy);

        Console.WriteLine("");
        Console.WriteLine("  Number of KMEANS_01 iterations taken is " + it_num + "");

        double[] cluster_variance = Cluster.cluster_variance_compute(dim_num, point_num, cluster_num,
            point, cluster, cluster_center);

        Cluster.cluster_print_summary(point_num, cluster_num, cluster_population,
            cluster_energy, cluster_variance);

        typeMethods.r8mat_write(center_filename, dim_num, cluster_num, cluster_center);

        Console.WriteLine("");
        Console.WriteLine("  Cluster centers written to \"" + center_filename + "\".");

        typeMethods.i4mat_write(cluster_filename, 1, point_num, cluster);

        Console.WriteLine("  Cluster assignments written to \"" + cluster_filename + "\".");
    }

    private static void test04()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 tries out the KMEANS_02 routine.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const string center_filename = "test04_centers.txt";
        const string cluster_filename = "test04_clusters.txt";
        int it_num = 0;
        const string point_filename = "points_100.txt";

        Console.WriteLine("");
        Console.WriteLine("TEST04");
        Console.WriteLine("  Test the KMEANS_02 algorithm.");
        Console.WriteLine("  (Applied Statistics Algorithm #136)");
        //
        //  Read the data points.
        //
        Console.WriteLine("");
        Console.WriteLine("  Data points will be read from \"" + point_filename + "\".");

        TableHeader h = typeMethods.r8mat_header_read(point_filename);
        int dim_num = h.m;
        int point_num = h.n;

        Console.WriteLine("");
        Console.WriteLine("  Point spatial dimension = " + dim_num + "");
        Console.WriteLine("  Number of points = " + point_num + "");

        double[] point = typeMethods.r8mat_data_read(point_filename, dim_num, point_num);
        //
        //  Clustering parameters.
        //
        const int cluster_num = 5;
        const int it_max = 20;

        int[] cluster = typeMethods.i4vec_negone_new(point_num);
        double[] cluster_energy = new double[cluster_num];
        int[] cluster_population = new int[cluster_num];

        Console.WriteLine("");
        Console.WriteLine("  Number of iterations allowed is " + it_max + "");
        //
        //  Initialize the centers.
        //
        double[] cluster_center = Cluster.cluster_initialize_1(dim_num, point_num, cluster_num, point);

        KMeans.kmeans_02(dim_num, point_num, cluster_num, it_max, ref it_num, point,
            ref cluster, ref cluster_center, ref cluster_population, ref cluster_energy);

        Console.WriteLine("");
        Console.WriteLine("  Number of iterations taken is " + it_num + "");

        double[] cluster_variance = Cluster.cluster_variance_compute(dim_num, point_num, cluster_num,
            point, cluster, cluster_center);

        Cluster.cluster_print_summary(point_num, cluster_num,
            cluster_population, cluster_energy, cluster_variance);

        typeMethods.r8mat_write(center_filename, dim_num, cluster_num, cluster_center);

        Console.WriteLine("");
        Console.WriteLine("  Cluster centers written to \"" + center_filename + "\".");

        typeMethods.i4mat_write(cluster_filename, 1, point_num, cluster);

        Console.WriteLine("  Cluster assignments written to \"" + cluster_filename + "\".");
    }

    private static void test05()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 tries out the KMEANS_03 routine.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const string center_filename = "test05_centers.txt";
        const string cluster_filename = "test05_clusters.txt";
        int it_num = 0;
        const string point_filename = "points_100.txt";

        Console.WriteLine("");
        Console.WriteLine("TEST05");
        Console.WriteLine("  Test the KMEANS_03 algorithm.");
        //
        //  Read the data points.
        //
        Console.WriteLine("");
        Console.WriteLine("  Data points will be read from \"" + point_filename + "\".");

        TableHeader h = typeMethods.r8mat_header_read(point_filename);
        int dim_num = h.m;
        int point_num = h.n;

        Console.WriteLine("");
        Console.WriteLine("  Point spatial dimension = " + dim_num + "");
        Console.WriteLine("  Number of points = " + point_num + "");

        double[] point = typeMethods.r8mat_data_read(point_filename, dim_num, point_num);
        //
        //  Clustering parameters.
        //
        const int cluster_num = 5;
        const int it_max = 20;

        int[] cluster = typeMethods.i4vec_negone_new(point_num);
        double[] cluster_energy = new double[cluster_num];
        int[] cluster_population = new int[cluster_num];

        Console.WriteLine("");
        Console.WriteLine("  Number of iterations allowed is " + it_max + "");
        //
        //  Initialize the centers.
        //
        double[] cluster_center = Cluster.cluster_initialize_1(dim_num, point_num, cluster_num,
            point);

        KMeans.kmeans_03(dim_num, point_num, cluster_num, it_max, ref it_num, point,
            ref cluster, ref cluster_center, ref cluster_population, ref cluster_energy);

        Console.WriteLine("");
        Console.WriteLine("  Number of iterations taken is " + it_num + "");

        double[] cluster_variance = Cluster.cluster_variance_compute(dim_num, point_num, cluster_num,
            point, cluster, cluster_center);

        Cluster.cluster_print_summary(point_num, cluster_num, cluster_population,
            cluster_energy, cluster_variance);

        typeMethods.r8mat_write(center_filename, dim_num, cluster_num, cluster_center);

        Console.WriteLine("");
        Console.WriteLine("  Cluster centers written to \"" + center_filename + "\".");

        typeMethods.i4mat_write(cluster_filename, 1, point_num, cluster);

        Console.WriteLine("  Cluster assignments written to \"" + cluster_filename + "\".");

    }

    private static void test06()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST06 tries out the HMEANS_01 + KMEANS_01 routine.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const string center_filename = "test06_centers.txt";
        const string cluster_filename = "test06_clusters.txt";
        int it_num = 0;
        const string point_filename = "points_100.txt";

        Console.WriteLine("");
        Console.WriteLine("TEST06");
        Console.WriteLine("  Test the HMEANS_01 + KMEANS_01 algorithm.");
        //
        //  Read the data points.
        //
        Console.WriteLine("");
        Console.WriteLine("  Data points will be read from \"" + point_filename + "\".");

        TableHeader h = typeMethods.r8mat_header_read(point_filename);
        int dim_num = h.m;
        int point_num = h.n;

        Console.WriteLine("");
        Console.WriteLine("  Point spatial dimension = " + dim_num + "");
        Console.WriteLine("  Number of points = " + point_num + "");

        double[] point = typeMethods.r8mat_data_read(point_filename, dim_num, point_num);
        //
        //  Clustering parameters.
        //
        const int cluster_num = 5;
        const int it_max = 20;
        int seed = 123456789;
        const int it_max_h = 3;
        const int it_max_k = it_max;

        int[] cluster = typeMethods.i4vec_negone_new(point_num);
        double[] cluster_energy = new double[cluster_num];
        int[] cluster_population = new int[cluster_num];

        Console.WriteLine("");
        Console.WriteLine("  Number of HMEANS_01 iterations allowed is " + it_max_h + "");
        Console.WriteLine("  Number of KMEANS_01 iterations allowed is " + it_max_k + "");
        //
        //  Initialize the centers.
        //
        double[] cluster_center = Cluster.cluster_initialize_5(dim_num, point_num, cluster_num,
            point, ref seed);

        HMeans.hmeans_01(dim_num, point_num, cluster_num, it_max_h, ref it_num, point,
            cluster, cluster_center, cluster_population, cluster_energy);

        Console.WriteLine("");
        Console.WriteLine("  Number of HMEANS_01 iterations taken is " + it_num + "");

        double[] cluster_variance = Cluster.cluster_variance_compute(dim_num, point_num, cluster_num,
            point, cluster, cluster_center);

        Cluster.cluster_print_summary(point_num, cluster_num, cluster_population,
            cluster_energy, cluster_variance);

        KMeans.kmeans_01(dim_num, point_num, cluster_num, it_max_k, ref it_num, point,
            ref cluster, ref cluster_center, ref cluster_population, ref cluster_energy);

        Console.WriteLine("");
        Console.WriteLine("  Number of KMEANS_01 iterations taken is " + it_num + "");

        cluster_variance = Cluster.cluster_variance_compute(dim_num, point_num, cluster_num,
            point, cluster, cluster_center);

        Cluster.cluster_print_summary(point_num, cluster_num, cluster_population,
            cluster_energy, cluster_variance);

        typeMethods.r8mat_write(center_filename, dim_num, cluster_num, cluster_center);

        Console.WriteLine("");
        Console.WriteLine("  Cluster centers written to \"" + center_filename + "\".");

        typeMethods.i4mat_write(cluster_filename, 1, point_num, cluster);

        Console.WriteLine("  Cluster assignments written to \"" + cluster_filename + "\".");

    }

    private static void test07()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST07 tries out the HMEANS_01 + KMEANS_02 routine.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const string center_filename = "test07_centers.txt";
        const string cluster_filename = "test07_clusters.txt";
        int it_num = 0;
        const string point_filename = "points_100.txt";

        Console.WriteLine("");
        Console.WriteLine("TEST07");
        Console.WriteLine("  Test the HMEANS_01 + KMEANS_02 algorithm.");
        //
        //  Read the data points.
        //
        Console.WriteLine("");
        Console.WriteLine("  Data points will be read from \"" + point_filename + "\".");

        TableHeader h = typeMethods.r8mat_header_read(point_filename);
        int dim_num = h.m;
        int point_num = h.n;

        Console.WriteLine("");
        Console.WriteLine("  Point spatial dimension = " + dim_num + "");
        Console.WriteLine("  Number of points = " + point_num + "");

        double[] point = typeMethods.r8mat_data_read(point_filename, dim_num, point_num);
        //
        //  Clustering parameters.
        //
        const int cluster_num = 5;
        const int it_max = 20;
        int seed = 123456789;

        const int it_max_h = 3;
        const int it_max_k = it_max;

        int[] cluster = typeMethods.i4vec_negone_new(point_num);
        double[] cluster_energy = new double[cluster_num];
        int[] cluster_population = new int[cluster_num];

        Console.WriteLine("");
        Console.WriteLine("  Number of HMEANS_01 iterations allowed is " + it_max_h + "");
        Console.WriteLine("  Number of KMEANS_02 iterations allowed is " + it_max_k + "");
        //
        //  Initialize the centers.
        //
        double[] cluster_center = Cluster.cluster_initialize_5(dim_num, point_num, cluster_num,
            point, ref seed);

        HMeans.hmeans_01(dim_num, point_num, cluster_num, it_max_h, ref it_num, point,
            cluster, cluster_center, cluster_population, cluster_energy);

        Console.WriteLine("");
        Console.WriteLine("  Number of HMEANS_01 iterations taken is " + it_num + "");

        double[] cluster_variance = Cluster.cluster_variance_compute(dim_num, point_num, cluster_num,
            point, cluster, cluster_center);

        Cluster.cluster_print_summary(point_num, cluster_num, cluster_population,
            cluster_energy, cluster_variance);

        KMeans.kmeans_02(dim_num, point_num, cluster_num, it_max_k, ref it_num, point,
            ref cluster, ref cluster_center, ref cluster_population, ref cluster_energy);

        Console.WriteLine("");
        Console.WriteLine("  Number of KMEANS_02 iterations taken is " + it_num + "");

        cluster_variance = Cluster.cluster_variance_compute(dim_num, point_num, cluster_num,
            point, cluster, cluster_center);

        Cluster.cluster_print_summary(point_num, cluster_num, cluster_population,
            cluster_energy, cluster_variance);

        typeMethods.r8mat_write(center_filename, dim_num, cluster_num, cluster_center);

        Console.WriteLine("");
        Console.WriteLine("  Cluster centers written to \"" + center_filename + "\".");

        typeMethods.i4mat_write(cluster_filename, 1, point_num, cluster);

        Console.WriteLine("  Cluster assignments written to \"" + cluster_filename + "\".");

    }

    private static void test08()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST08 tries out the HMEANS_01 + KMEANS_03 routine.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const string center_filename = "test08_centers.txt";
        const string cluster_filename = "test08_clusters.txt";
        int it_num = 0;
        const string point_filename = "points_100.txt";

        Console.WriteLine("");
        Console.WriteLine("TEST08");
        Console.WriteLine("  Test the HMEANS_01 + KMEANS_03 algorithm.");
        //
        //  Read the data points.
        //
        Console.WriteLine("");
        Console.WriteLine("  Data points will be read from \"" + point_filename + "\".");

        TableHeader h = typeMethods.r8mat_header_read(point_filename);
        int dim_num = h.m;
        int point_num = h.n;

        Console.WriteLine("");
        Console.WriteLine("  Point spatial dimension = " + dim_num + "");
        Console.WriteLine("  Number of points = " + point_num + "");

        double[] point = typeMethods.r8mat_data_read(point_filename, dim_num, point_num);
        //
        //  Clustering parameters.
        //
        const int cluster_num = 5;
        const int it_max = 20;
        int seed = 123456789;
        const int it_max_h = 3;
        const int it_max_k = it_max;

        int[] cluster = typeMethods.i4vec_negone_new(point_num);
        double[] cluster_energy = new double[cluster_num];
        int[] cluster_population = new int[cluster_num];

        Console.WriteLine("");
        Console.WriteLine("  Initialize by using a few steps of HMEANS_02:");
        Console.WriteLine("  Number of HMEANS_01 iterations allowed is " + it_max_h + "");
        Console.WriteLine("  Number of KMEANS_03 iterations allowed is " + it_max_k + "");
        //
        //  Initialize the centers.
        //
        double[] cluster_center = Cluster.cluster_initialize_5(dim_num, point_num, cluster_num,
            point, ref seed);
        //
        //  Initialize the clusters.
        //
        HMeans.hmeans_01(dim_num, point_num, cluster_num, it_max_h, ref it_num, point,
            cluster, cluster_center, cluster_population, cluster_energy);

        Console.WriteLine("");
        Console.WriteLine("  Number of HMEANS_01 iterations taken is " + it_num + "");

        double[] cluster_variance = Cluster.cluster_variance_compute(dim_num, point_num, cluster_num,
            point, cluster, cluster_center);

        Cluster.cluster_print_summary(point_num, cluster_num, cluster_population,
            cluster_energy, cluster_variance);

        KMeans.kmeans_03(dim_num, point_num, cluster_num, it_max_k, ref it_num, point,
            ref cluster, ref cluster_center, ref cluster_population, ref cluster_energy);

        Console.WriteLine("");
        Console.WriteLine("  Number of KMEANS_03 iterations taken is " + it_num + "");

        cluster_variance = Cluster.cluster_variance_compute(dim_num, point_num, cluster_num,
            point, cluster, cluster_center);

        Cluster.cluster_print_summary(point_num, cluster_num, cluster_population,
            cluster_energy, cluster_variance);

        typeMethods.r8mat_write(center_filename, dim_num, cluster_num, cluster_center);

        Console.WriteLine("");
        Console.WriteLine("  Cluster centers written to \"" + center_filename + "\".");

        typeMethods.i4mat_write(cluster_filename, 1, point_num, cluster);

        Console.WriteLine("  Cluster assignments written to \"" + cluster_filename + "\".");

    }

    private static void test09()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST09 tries out the HMEANS_W_01 routine.
        //
        //  Discussion:
        //
        //    The weights are all equal, so the results should
        //    be identical to those for HMEANS_01.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const string center_filename = "test09_centers.txt";
        const string cluster_filename = "test09_clusters.txt";
        int it_num = 0;
        const string point_filename = "points_100.txt";
        const string weight_filename = "weights_equal_100.txt";

        Console.WriteLine("");
        Console.WriteLine("TEST09");
        Console.WriteLine("  Test the HMEANS_W_01 algorithm.");
        //
        //  Read the data points.
        //
        Console.WriteLine("");
        Console.WriteLine("  Data points will be read from \"" + point_filename + "\".");

        TableHeader h = typeMethods.r8mat_header_read(point_filename);
        int dim_num = h.m;
        int point_num = h.n;

        Console.WriteLine("");
        Console.WriteLine("  Point spatial dimension = " + dim_num + "");
        Console.WriteLine("  Number of points = " + point_num + "");

        double[] point = typeMethods.r8mat_data_read(point_filename, dim_num, point_num);
        //
        //  Read the weights.
        //
        Console.WriteLine("");
        Console.WriteLine("  Weights will be read from \"" + weight_filename + "\".");

        h = typeMethods.r8mat_header_read(weight_filename);
        int weight_dim = h.m;
        int weight_num = h.n;

        if (weight_dim != 1)
        {
            Console.WriteLine("");
            Console.WriteLine("Fatal error!");
            Console.WriteLine("  Spatial dimension of weight array is not 1.");
            return;
        }

        if (weight_num != point_num)
        {
            Console.WriteLine("");
            Console.WriteLine("Fatal error!");
            Console.WriteLine("  Number of weights not equal to number of points.");
            return;
        }

        double[] weight = typeMethods.r8mat_data_read(weight_filename, weight_dim, weight_num);
        //
        //  Clustering parameters.
        //
        const int cluster_num = 5;
        const int it_max = 20;
        int seed = 123456789;

        int[] cluster = typeMethods.i4vec_negone_new(point_num);
        double[] cluster_energy = new double[cluster_num];
        int[] cluster_population = new int[cluster_num];

        Console.WriteLine("");
        Console.WriteLine("  Number of iterations allowed is " + it_max + "");
        //
        //  Initialize the cluster centers.
        //
        double[] cluster_center = Cluster.cluster_initialize_5(dim_num, point_num, cluster_num,
            point, ref seed);

        HMeans.hmeans_w_01(dim_num, point_num, cluster_num, it_max, ref it_num,
            point, weight, cluster, cluster_center, cluster_population,
            cluster_energy);

        Console.WriteLine("");
        Console.WriteLine("  Number of iterations taken is " + it_num + "");

        double[] cluster_variance = Cluster.cluster_variance_compute(dim_num, point_num, cluster_num,
            point, cluster, cluster_center);

        Cluster.cluster_print_summary(point_num, cluster_num, cluster_population,
            cluster_energy, cluster_variance);

        typeMethods.r8mat_write(center_filename, dim_num, cluster_num, cluster_center);

        Console.WriteLine("");
        Console.WriteLine("  Cluster centers written to \"" + center_filename + "\".");

        typeMethods.i4mat_write(cluster_filename, 1, point_num, cluster);

        Console.WriteLine("  Cluster assignments written to \"" + cluster_filename + "\".");

    }

    private static void test10()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST10 tries out the HMEANS_W_02 routine.
        //
        //  Discussion:
        //
        //    The weights are all equal, so the results should
        //    be identical to those for HMEANS_02.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const string center_filename = "test10_centers.txt";
        const string cluster_filename = "test10_clusters.txt";
        int it_num = 0;
        const string point_filename = "points_100.txt";
        const string weight_filename = "weights_equal_100.txt";

        Console.WriteLine("");
        Console.WriteLine("TEST10");
        Console.WriteLine("  Test the HMEANS_W_02 algorithm.");
        //
        //  Read the data points.
        //
        Console.WriteLine("");
        Console.WriteLine("  Data points will be read from \"" + point_filename + "\".");

        TableHeader h = typeMethods.r8mat_header_read(point_filename);
        int dim_num = h.m;
        int point_num = h.n;

        Console.WriteLine("");
        Console.WriteLine("  Point spatial dimension = " + dim_num + "");
        Console.WriteLine("  Number of points = " + point_num + "");

        double[] point = typeMethods.r8mat_data_read(point_filename, dim_num, point_num);
        //
        //  Read the weights.
        //
        Console.WriteLine("");
        Console.WriteLine("  Weights will be read from \"" + weight_filename + "\".");

        h = typeMethods.r8mat_header_read(weight_filename);
        int weight_dim = h.m;
        int weight_num = h.n;

        if (weight_dim != 1)
        {
            Console.WriteLine("");
            Console.WriteLine("Fatal error!");
            Console.WriteLine("  Spatial dimension of weight array is not 1.");
            return;
        }

        if (weight_num != point_num)
        {
            Console.WriteLine("");
            Console.WriteLine("Fatal error!");
            Console.WriteLine("  Number of weights not equal to number of points.");
            return;
        }

        double[] weight = typeMethods.r8mat_data_read(weight_filename, weight_dim, weight_num);
        //
        //  Clustering parameters.
        //
        const int cluster_num = 5;
        const int it_max = 20;
        int seed = 123456789;

        int[] cluster = typeMethods.i4vec_negone_new(point_num);
        double[] cluster_energy = new double[cluster_num];
        int[] cluster_population = new int[cluster_num];

        Console.WriteLine("");
        Console.WriteLine("  Number of iterations allowed is " + it_max + "");
        //
        //  Initialize the cluster centers.
        //
        double[] cluster_center = Cluster.cluster_initialize_5(dim_num, point_num, cluster_num,
            point, ref seed);

        HMeans.hmeans_w_02(dim_num, point_num, cluster_num, it_max, ref it_num,
            point, weight, cluster, cluster_center, cluster_population,
            cluster_energy, ref seed);

        Console.WriteLine("");
        Console.WriteLine("  Number of iterations taken is " + it_num + "");

        double[] cluster_variance = Cluster.cluster_variance_compute(dim_num, point_num, cluster_num,
            point, cluster, cluster_center);

        Cluster.cluster_print_summary(point_num, cluster_num, cluster_population,
            cluster_energy, cluster_variance);

        typeMethods.r8mat_write(center_filename, dim_num, cluster_num, cluster_center);

        Console.WriteLine("");
        Console.WriteLine("  Cluster centers written to \"" + center_filename + "\".");

        typeMethods.i4mat_write(cluster_filename, 1, point_num, cluster);

        Console.WriteLine("  Cluster assignments written to \"" + cluster_filename + "\".");

    }

    private static void test11()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST11 tries out the KMEANS_W_01 routine.
        //
        //  Discussion:
        //
        //   The weights are all equal, so the results should
        //    be identical to those for KMEANS_01.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    12 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const string center_filename = "test11_centers.txt";
        const string cluster_filename = "test11_clusters.txt";
        int it_num = 0;
        const string point_filename = "points_100.txt";
        const string weight_filename = "weights_equal_100.txt";

        Console.WriteLine("");
        Console.WriteLine("TEST11");
        Console.WriteLine("  Test the KMEANS_W_01 algorithm.");
        //
        //  Read the data points.
        //
        Console.WriteLine("");
        Console.WriteLine("  Data points will be read from \"" + point_filename + "\".");

        TableHeader h = typeMethods.r8mat_header_read(point_filename);
        int dim_num = h.m;
        int point_num = h.n;

        Console.WriteLine("");
        Console.WriteLine("  Point spatial dimension = " + dim_num + "");
        Console.WriteLine("  Number of points = " + point_num + "");

        double[] point = typeMethods.r8mat_data_read(point_filename, dim_num, point_num);
        //
        //  Read the weights.
        //
        Console.WriteLine("");
        Console.WriteLine("  Weights will be read from \"" + weight_filename + "\".");

        h = typeMethods.r8mat_header_read(weight_filename);
        int weight_dim = h.m;
        int weight_num = h.n;

        if (weight_dim != 1)
        {
            Console.WriteLine("");
            Console.WriteLine("Fatal error!");
            Console.WriteLine("  Spatial dimension of weight array is not 1.");
            return;
        }

        if (weight_num != point_num)
        {
            Console.WriteLine("");
            Console.WriteLine("Fatal error!");
            Console.WriteLine("  Number of weights not equal to number of points.");
            return;
        }

        double[] weight = typeMethods.r8mat_data_read(weight_filename, weight_dim, weight_num);
        //
        //  Clustering parameters.
        //
        const int cluster_num = 5;
        const int it_max = 20;
        int seed = 123456789;

        int[] cluster = typeMethods.i4vec_negone_new(point_num);
        double[] cluster_energy = new double[cluster_num];
        int[] cluster_population = new int[cluster_num];

        Console.WriteLine("");
        Console.WriteLine("  Number of iterations allowed is " + it_max + "");
        //
        //  Initialize the cluster centers.
        //
        double[] cluster_center = Cluster.cluster_initialize_5(dim_num, point_num, cluster_num,
            point, ref seed);

        KMeans.kmeans_w_01(dim_num, point_num, cluster_num, it_max, ref it_num,
            point, weight, ref cluster, ref cluster_center, ref cluster_population,
            ref cluster_energy);

        Console.WriteLine("");
        Console.WriteLine("  Number of iterations taken is " + it_num + "");

        double[] cluster_variance = Cluster.cluster_variance_compute(dim_num, point_num, cluster_num,
            point, cluster, cluster_center);

        Cluster.cluster_print_summary(point_num, cluster_num,
            cluster_population, cluster_energy, cluster_variance);

        typeMethods.r8mat_write(center_filename, dim_num, cluster_num, cluster_center);

        Console.WriteLine("");
        Console.WriteLine("  Cluster centers written to \"" + center_filename + "\".");

        typeMethods.i4mat_write(cluster_filename, 1, point_num, cluster);

        Console.WriteLine("  Cluster assignments written to \"" + cluster_filename + "\".");

    }

    private static void test12()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST12 tries out the KMEANS_W_03 routine.
        //
        //  Discussion:
        //
        //    The weights are all equal, so the results should
        //    be identical to those for KMEANS_03.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const string center_filename = "test12_centers.txt";
        const string cluster_filename = "test12_clusters.txt";
        int it_num = 0;
        const string point_filename = "points_100.txt";
        const string weight_filename = "weights_equal_100.txt";

        Console.WriteLine("");
        Console.WriteLine("TEST12");
        Console.WriteLine("  Test the KMEANS_W_03 algorithm.");
        //
        //  Read the data points.
        //
        Console.WriteLine("");
        Console.WriteLine("  Data points will be read from \"" + point_filename + "\".");

        TableHeader h = typeMethods.r8mat_header_read(point_filename);
        int dim_num = h.m;
        int point_num = h.n;

        Console.WriteLine("");
        Console.WriteLine("  Point spatial dimension = " + dim_num + "");
        Console.WriteLine("  Number of points = " + point_num + "");

        double[] point = typeMethods.r8mat_data_read(point_filename, dim_num, point_num);
        //
        //  Read the weights.
        //
        Console.WriteLine("");
        Console.WriteLine("  Weights will be read from \"" + weight_filename + "\".");

        h = typeMethods.r8mat_header_read(weight_filename);
        int weight_dim = h.m;
        int weight_num = h.n;

        if (weight_dim != 1)
        {
            Console.WriteLine("");
            Console.WriteLine("Fatal error!");
            Console.WriteLine("  Spatial dimension of weight array is not 1.");
            return;
        }

        if (weight_num != point_num)
        {
            Console.WriteLine("");
            Console.WriteLine("Fatal error!");
            Console.WriteLine("  Number of weights not equal to number of points.");
            return;
        }

        double[] weight = typeMethods.r8mat_data_read(weight_filename, weight_dim, weight_num);
        //
        //  Clustering parameters.
        //
        const int cluster_num = 5;
        const int it_max = 20;

        int[] cluster = typeMethods.i4vec_negone_new(point_num);
        double[] cluster_energy = new double[cluster_num];
        int[] cluster_population = new int[cluster_num];

        Console.WriteLine("");
        Console.WriteLine("  Number of iterations allowed is " + it_max + "");
        //
        //  Initialize the cluster centers.
        //
        double[] cluster_center = Cluster.cluster_initialize_1(dim_num, point_num, cluster_num,
            point);

        KMeans.kmeans_w_03(dim_num, point_num, cluster_num, it_max, ref it_num,
            point, weight, ref cluster, ref cluster_center, ref cluster_population,
            ref cluster_energy);

        Console.WriteLine("");
        Console.WriteLine("  Number of iterations taken is " + it_num + "");

        double[] cluster_variance = Cluster.cluster_variance_compute(dim_num, point_num, cluster_num,
            point, cluster, cluster_center);

        Cluster.cluster_print_summary(point_num, cluster_num,
            cluster_population, cluster_energy, cluster_variance);

        typeMethods.r8mat_write(center_filename, dim_num, cluster_num, cluster_center);

        Console.WriteLine("");
        Console.WriteLine("  Cluster centers written to \"" + center_filename + "\".");

        typeMethods.i4mat_write(cluster_filename, 1, point_num, cluster);

        Console.WriteLine("  Cluster assignments written to \"" + cluster_filename + "\".");

    }

    private static void test13()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST13 tries out the HMEANS_W_01 routine.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const string center_filename = "test13_centers.txt";
        const string cluster_filename = "test13_clusters.txt";
        int it_num = 0;
        const string point_filename = "points_100.txt";
        const string weight_filename = "weights_unequal_100.txt";

        Console.WriteLine("");
        Console.WriteLine("TEST13");
        Console.WriteLine("  Test the HMEANS_W_01 algorithm.");
        //
        //  Read the data points.
        //
        Console.WriteLine("");
        Console.WriteLine("  Data points will be read from \"" + point_filename + "\".");

        TableHeader h = typeMethods.r8mat_header_read(point_filename);
        int dim_num = h.m;
        int point_num = h.n;

        Console.WriteLine("");
        Console.WriteLine("  Point spatial dimension = " + dim_num + "");
        Console.WriteLine("  Number of points = " + point_num + "");

        double[] point = typeMethods.r8mat_data_read(point_filename, dim_num, point_num);
        //
        //  Read the weights.
        //
        Console.WriteLine("");
        Console.WriteLine("  Weights will be read from \"" + weight_filename + "\".");

        h = typeMethods.r8mat_header_read(weight_filename);
        int weight_dim = h.m;
        int weight_num = h.n;

        if (weight_dim != 1)
        {
            Console.WriteLine("");
            Console.WriteLine("Fatal error!");
            Console.WriteLine("  Spatial dimension of weight array is not 1.");
            return;
        }

        if (weight_num != point_num)
        {
            Console.WriteLine("");
            Console.WriteLine("Fatal error!");
            Console.WriteLine("  Number of weights not equal to number of points.");
            return;
        }

        double[] weight = typeMethods.r8mat_data_read(weight_filename, weight_dim, weight_num);
        //
        //  Clustering parameters.
        //
        const int cluster_num = 5;
        const int it_max = 20;
        int seed = 123456789;

        int[] cluster = typeMethods.i4vec_negone_new(point_num);
        double[] cluster_energy = new double[cluster_num];
        int[] cluster_population = new int[cluster_num];

        Console.WriteLine("");
        Console.WriteLine("  Number of iterations allowed is " + it_max + "");
        //
        //  Initialize the cluster centers.
        //
        double[] cluster_center = Cluster.cluster_initialize_5(dim_num, point_num, cluster_num,
            point, ref seed);

        HMeans.hmeans_w_01(dim_num, point_num, cluster_num, it_max, ref it_num,
            point, weight, cluster, cluster_center, cluster_population,
            cluster_energy);

        Console.WriteLine("");
        Console.WriteLine("  Number of iterations taken is " + it_num + "");

        double[] cluster_variance = Cluster.cluster_variance_compute(dim_num, point_num, cluster_num,
            point, cluster, cluster_center);

        Cluster.cluster_print_summary(point_num, cluster_num,
            cluster_population, cluster_energy, cluster_variance);

        typeMethods.r8mat_write(center_filename, dim_num, cluster_num, cluster_center);

        Console.WriteLine("");
        Console.WriteLine("  Cluster centers written to \"" + center_filename + "\".");

        typeMethods.i4mat_write(cluster_filename, 1, point_num, cluster);

        Console.WriteLine("  Cluster assignments written to \"" + cluster_filename + "\".");

    }

    private static void test14()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST14 tries out the HMEANS_W_02 routine.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const string center_filename = "test14_centers.txt";
        const string cluster_filename = "test14_clusters.txt";
        int it_num = 0;
        const string point_filename = "points_100.txt";
        const string weight_filename = "weights_unequal_100.txt";

        Console.WriteLine("");
        Console.WriteLine("TEST14");
        Console.WriteLine("  Test the HMEANS_W_02 algorithm.");
        //
        //  Read the data points.
        //
        Console.WriteLine("");
        Console.WriteLine("  Data points will be read from \"" + point_filename + "\".");

        TableHeader h = typeMethods.r8mat_header_read(point_filename);
        int dim_num = h.m;
        int point_num = h.n;

        Console.WriteLine("");
        Console.WriteLine("  Point spatial dimension = " + dim_num + "");
        Console.WriteLine("  Number of points = " + point_num + "");

        double[] point = typeMethods.r8mat_data_read(point_filename, dim_num, point_num);
        //
        //  Read the weights.
        //
        Console.WriteLine("");
        Console.WriteLine("  Weights will be read from \"" + weight_filename + "\".");

        h = typeMethods.r8mat_header_read(weight_filename);
        int weight_dim = h.m;
        int weight_num = h.n;

        if (weight_dim != 1)
        {
            Console.WriteLine("");
            Console.WriteLine("Fatal error!");
            Console.WriteLine("  Spatial dimension of weight array is not 1.");
            return;
        }

        if (weight_num != point_num)
        {
            Console.WriteLine("");
            Console.WriteLine("Fatal error!");
            Console.WriteLine("  Number of weights not equal to number of points.");
            return;
        }

        double[] weight = typeMethods.r8mat_data_read(weight_filename, weight_dim, weight_num);
        //
        //  Clustering parameters.
        //
        const int cluster_num = 5;
        const int it_max = 20;
        int seed = 123456789;

        int[] cluster = typeMethods.i4vec_negone_new(point_num);
        double[] cluster_energy = new double[cluster_num];
        int[] cluster_population = new int[cluster_num];

        Console.WriteLine("");
        Console.WriteLine("  Number of iterations allowed is " + it_max + "");
        //
        //  Initialize the cluster centers.
        //
        double[] cluster_center = Cluster.cluster_initialize_5(dim_num, point_num, cluster_num,
            point, ref seed);

        HMeans.hmeans_w_02(dim_num, point_num, cluster_num, it_max, ref it_num,
            point, weight, cluster, cluster_center, cluster_population,
            cluster_energy, ref seed);

        Console.WriteLine("");
        Console.WriteLine("  Number of iterations taken is " + it_num + "");

        double[] cluster_variance = Cluster.cluster_variance_compute(dim_num, point_num, cluster_num,
            point, cluster, cluster_center);

        Cluster.cluster_print_summary(point_num, cluster_num,
            cluster_population, cluster_energy, cluster_variance);

        typeMethods.r8mat_write(center_filename, dim_num, cluster_num, cluster_center);

        Console.WriteLine("");
        Console.WriteLine("  Cluster centers written to \"" + center_filename + "\".");

        typeMethods.i4mat_write(cluster_filename, 1, point_num, cluster);

        Console.WriteLine("  Cluster assignments written to \"" + cluster_filename + "\".");

    }

    private static void test15()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST15 tries out the KMEANS_W_01 routine.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const string center_filename = "test15_centers.txt";
        const string cluster_filename = "test15_clusters.txt";
        int it_num = 0;
        const string point_filename = "points_100.txt";
        const string weight_filename = "weights_unequal_100.txt";

        Console.WriteLine("");
        Console.WriteLine("TEST15");
        Console.WriteLine("  Test the KMEANS_W_01 algorithm.");
        //
        //  Read the data points.
        //
        Console.WriteLine("");
        Console.WriteLine("  Data points will be read from \"" + point_filename + "\".");

        TableHeader h = typeMethods.r8mat_header_read(point_filename);
        int dim_num = h.m;
        int point_num = h.n;

        Console.WriteLine("");
        Console.WriteLine("  Point spatial dimension = " + dim_num + "");
        Console.WriteLine("  Number of points = " + point_num + "");

        double[] point = typeMethods.r8mat_data_read(point_filename, dim_num, point_num);
        //
        //  Read the weights.
        //
        Console.WriteLine("");
        Console.WriteLine("  Weights will be read from \"" + weight_filename + "\".");

        h = typeMethods.r8mat_header_read(weight_filename);
        int weight_dim = h.m;
        int weight_num = h.n;

        if (weight_dim != 1)
        {
            Console.WriteLine("");
            Console.WriteLine("Fatal error!");
            Console.WriteLine("  Spatial dimension of weight array is not 1.");
            return;
        }

        if (weight_num != point_num)
        {
            Console.WriteLine("");
            Console.WriteLine("Fatal error!");
            Console.WriteLine("  Number of weights not equal to number of points.");
            return;
        }

        double[] weight = typeMethods.r8mat_data_read(weight_filename, weight_dim, weight_num);
        //
        //  Clustering parameters.
        //
        const int cluster_num = 5;
        const int it_max = 20;
        int seed = 123456789;

        int[] cluster = typeMethods.i4vec_negone_new(point_num);
        double[] cluster_energy = new double[cluster_num];
        int[] cluster_population = new int[cluster_num];

        Console.WriteLine("");
        Console.WriteLine("  Number of iterations allowed is " + it_max + "");
        //
        //  Initialize the cluster centers.
        //
        double[] cluster_center = Cluster.cluster_initialize_5(dim_num, point_num, cluster_num,
            point, ref seed);

        KMeans.kmeans_w_01(dim_num, point_num, cluster_num, it_max, ref it_num,
            point, weight, ref cluster, ref cluster_center, ref cluster_population,
            ref cluster_energy);

        Console.WriteLine("");
        Console.WriteLine("  Number of iterations taken is " + it_num + "");

        double[] cluster_variance = Cluster.cluster_variance_compute(dim_num, point_num, cluster_num,
            point, cluster, cluster_center);

        Cluster.cluster_print_summary(point_num, cluster_num,
            cluster_population, cluster_energy, cluster_variance);

        typeMethods.r8mat_write(center_filename, dim_num, cluster_num, cluster_center);

        Console.WriteLine("");
        Console.WriteLine("  Cluster centers written to \"" + center_filename + "\".");

        typeMethods.i4mat_write(cluster_filename, 1, point_num, cluster);

        Console.WriteLine("  Cluster assignments written to \"" + cluster_filename + "\".");
    }

    private static void test16()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST16 tries out the KMEANS_W_03 routine.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    11 October 2011
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        const string center_filename = "test16_centers.txt";
        const string cluster_filename = "test16_clusters.txt";
        int it_num = 0;
        const string point_filename = "points_100.txt";
        const string weight_filename = "weights_unequal_100.txt";

        Console.WriteLine("");
        Console.WriteLine("TEST16");
        Console.WriteLine("  Test the KMEANS_W_03 algorithm.");
        //
        //  Read the data points.
        //
        Console.WriteLine("");
        Console.WriteLine("  Data points will be read from \"" + point_filename + "\".");

        TableHeader h = typeMethods.r8mat_header_read(point_filename);
        int dim_num = h.m;
        int point_num = h.n;

        Console.WriteLine("");
        Console.WriteLine("  Point spatial dimension = " + dim_num + "");
        Console.WriteLine("  Number of points = " + point_num + "");

        double[] point = typeMethods.r8mat_data_read(point_filename, dim_num, point_num);
        //
        //  Read the weights.
        //
        Console.WriteLine("");
        Console.WriteLine("  Weights will be read from \"" + weight_filename + "\".");

        h = typeMethods.r8mat_header_read(weight_filename);
        int weight_dim = h.m;
        int weight_num = h.n;

        if (weight_dim != 1)
        {
            Console.WriteLine("");
            Console.WriteLine("Fatal error!");
            Console.WriteLine("  Spatial dimension of weight array is not 1.");
            return;
        }

        if (weight_num != point_num)
        {
            Console.WriteLine("");
            Console.WriteLine("Fatal error!");
            Console.WriteLine("  Number of weights not equal to number of points.");
            return;
        }

        double[] weight = typeMethods.r8mat_data_read(weight_filename, weight_dim, weight_num);
        //
        //  Clustering parameters.
        //
        const int cluster_num = 5;
        const int it_max = 20;

        int[] cluster = typeMethods.i4vec_negone_new(point_num);
        double[] cluster_energy = new double[cluster_num];
        int[] cluster_population = new int[cluster_num];

        Console.WriteLine("");
        Console.WriteLine("  Number of iterations allowed is " + it_max + "");
        //
        //  Initialize the cluster centers.
        //
        double[] cluster_center = Cluster.cluster_initialize_1(dim_num, point_num, cluster_num, point);

        KMeans.kmeans_w_03(dim_num, point_num, cluster_num, it_max, ref it_num,
            point, weight, ref cluster, ref cluster_center, ref cluster_population,
            ref cluster_energy);

        Console.WriteLine("");
        Console.WriteLine("  Number of iterations taken is " + it_num + "");

        double[] cluster_variance = Cluster.cluster_variance_compute(dim_num, point_num, cluster_num,
            point, cluster, cluster_center);

        Cluster.cluster_print_summary(point_num, cluster_num,
            cluster_population, cluster_energy, cluster_variance);

        typeMethods.r8mat_write(center_filename, dim_num, cluster_num, cluster_center);

        Console.WriteLine("");
        Console.WriteLine("  Cluster centers written to \"" + center_filename + "\".");

        typeMethods.i4mat_write(cluster_filename, 1, point_num, cluster);

        Console.WriteLine("  Cluster assignments written to \"" + cluster_filename + "\".");
    }
}