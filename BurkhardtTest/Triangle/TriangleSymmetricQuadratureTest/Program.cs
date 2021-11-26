﻿using System;
using System.Collections.Generic;
using System.IO;
using Burkardt.PolynomialNS;
using Burkardt.TriangleNS;
using Burkardt.Types;
using Burkardt.Uniform;

namespace TriangleSymmetricQuadratureTest;

internal static class Program
{
    private static void Main()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    MAIN is the main program for TRIANGLE_SYMQ_RULE_TEST.
        //
        //  Discussion:
        //
        //    TRIANGLE_SYMQ_RULE_TEST tests the TRIANGLE_SYMQ_RULE library.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU GPL license.
        //
        //  Modified:
        //
        //    30 June 2014
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Hong Xiao, Zydrunas Gimbutas,
        //    A numerical algorithm for the construction of efficient quadrature
        //    rules in two and higher dimensions,
        //    Computers and Mathematics with Applications,
        //    Volume 59, 2010, pages 663-676.
        //
    {
        int degree = 0;
        string header = "";
        int itype;
        double[] vert1 = new double[2];
        double[] vert2 = new double[2];
        double[] vert3 = new double[2];

        Console.WriteLine("");
        Console.WriteLine("TRIANGLE_SYMQ_RULE_TEST");
        Console.WriteLine("  Test the TRIANGLE_SYMQ_RULE library.");

        test01();

        for (itype = 0; itype <= 2; itype++)
        {
            switch (itype)
            {
                case 0:
                    Console.WriteLine("");
                    Console.WriteLine("  Region is user-defined triangle.");
                    vert1[0] = 1.0;
                    vert1[1] = 0.0;
                    vert2[0] = 4.0;
                    vert2[1] = 4.0;
                    vert3[0] = 0.0;
                    vert3[1] = 3.0;
                    header = "user08";
                    degree = 8;
                    break;
                case 1:
                    Console.WriteLine("");
                    Console.WriteLine("  Region is standard equilateral triangle.");
                    vert1[0] = -1.0;
                    vert1[1] = -1.0 / Math.Sqrt(3.0);
                    vert2[0] = +1.0;
                    vert2[1] = -1.0 / Math.Sqrt(3.0);
                    vert3[0] = 0.0;
                    vert3[1] = 2.0 / Math.Sqrt(3.0);
                    header = "equi08";
                    degree = 8;
                    break;
                case 2:
                    Console.WriteLine("");
                    Console.WriteLine("  Region is the simplex (0,0),(1,0),(0,1).");
                    vert1[0] = 0.0;
                    vert1[1] = 0.0;
                    vert2[0] = 1.0;
                    vert2[1] = 0.0;
                    vert3[0] = 0.0;
                    vert3[1] = 1.0;
                    header = "simp08";
                    degree = 8;
                    break;
            }

            Console.WriteLine("");
            Console.WriteLine("  Triangle:");
            Console.WriteLine("");
            Console.WriteLine(vert1[0] + "  " + vert1[1] + "");
            Console.WriteLine(vert2[0] + "  " + vert2[1] + "");
            Console.WriteLine(vert3[0] + "  " + vert3[1] + "");
            //
            //  Determine the size of the rule.
            //
            int numnodes = rule_full_size(degree);
            //
            //  Retrieve a rule and print it.
            //
            test02(degree, numnodes, vert1, vert2, vert3);
            //
            //  Get a rule, and write data files that gnuplot can use to plot the points.
            //
            test03(degree, numnodes, vert1, vert2, vert3, header);

            test04(degree, numnodes, vert1, vert2, vert3, header);

            test05(degree, numnodes, vert1, vert2, vert3);
        }

        Console.WriteLine("");
        Console.WriteLine("TRIANGLE_SYMQ_RULE_TEST");
        Console.WriteLine("  Normal end of execution.");
        Console.WriteLine("");
    }

    private static int rule_full_size ( int mmax )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    RULE_FULL_SIZE returns the full size of the requested quadrature rule.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU GPL license.
        //
        //  Modified:
        //
        //    30 June 2014
        //
        //  Author:
        //
        //    John Burkardt.
        //
        //  Reference:
        //
        //    Hong Xiao, Zydrunas Gimbutas,
        //    A numerical algorithm for the construction of efficient quadrature
        //    rules in two and higher dimensions,
        //    Computers and Mathematics with Applications,
        //    Volume 59, 2010, pages 663-676.
        //
        //  Parameters:
        //
        //    Input, int MMAX, the degree of the quadrature (the maximum degree of
        //    the polynomials of two variables that are integrated
        //    exactly.  1 <= MMAX <= 50.
        //
        //    Output, int RULE_FULL_SIZE, the number of nodes in the full rule.
        //
    {
        int[] nnodes = {
            1,   3,   6,   6,   7,  12,  15,  16,  19,  25, 
            28,  33,  37,  42,  49,  55,  60,  67,  73,  79, 
            87,  96, 103, 112, 120, 130, 141, 150, 159, 171, 
            181, 193, 204, 214, 228, 243, 252, 267, 282, 295, 
            309, 324, 339, 354, 370, 385, 399, 423, 435, 453 };
        int npts = 0;

        switch (mmax)
        {
            case >= 1 and <= 50:
                npts = nnodes[mmax-1];
                break;
            default:
                Console.WriteLine("");
                Console.WriteLine("RULE_FULL_SIZE - Fatal error!");
                Console.WriteLine("  Degree MMAX must be between 1 and 50.");
                break;
        }

        return npts;
    }

    public static void test01()

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests TRIANGLE_TO_SIMPLEX, TRIANGLE_TO_REF, REF_TO_TRIANGLE, SIMPLEX_TO_TRIANGLE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU GPL license.
        //
        //  Modified:
        //
        //    30 June 2014
        //
        //  Author:
        //
        //    John Burkardt
        //
    {
        int i;
        double[] tv1 = new double[2];
        double[] tv2 = new double[2];
        double[] tv3 = new double[2];

        Console.WriteLine("");
        Console.WriteLine("TEST01");
        Console.WriteLine("  Map points from one triangle to another.");
        Console.WriteLine("");
        Console.WriteLine("  R = reference triangle");
        Console.WriteLine("  S = simplex");
        Console.WriteLine("  T = user-defined triangle.");
        Console.WriteLine("  REF_TO_TRIANGLE:     R => T");
        Console.WriteLine("  SIMPLEX_TO_TRIANGLE: S => T");
        Console.WriteLine("  TRIANGLE_TO_REF:     T => R");
        Console.WriteLine("  TRIANGLE_TO_SIMPLEX: T => S");
        //
        //  Reference triangle
        //
        //  rv1[0] = -1.0;
        //  rv1[1] = -1.0 / Math.Sqrt ( 3.0 );
        //  rv2[0] = +1.0;
        //  rv2[1] = -1.0 / Math.Sqrt ( 3.0 );
        //  rv3[0] =  0.0;
        //  rv3[1] =  2.0 / Math.Sqrt ( 3.0 );
        //
        //  Simplex
        //
        //  sv1[0] = 0.0;
        //  sv1[1] = 0.0;
        //  sv2[0] = 1.0;
        //  sv2[1] = 0.0;
        //  sv3[0] = 0.0;
        //  sv3[1] = 1.0;
        //
        //  User triangle.
        //
        tv1[0] = 1.0;
        tv1[1] = 0.0;
        tv2[0] = 4.0;
        tv2[1] = 4.0;
        tv3[0] = 0.0;
        tv3[1] = 3.0;

        int seed = 123456789;

        for (i = 1; i <= 5; i++)
        {
            double[] sp1 = UniformRNG.r8vec_uniform_01_new(2, ref seed);

            switch (sp1[0] + sp1[1])
            {
                case > 1.0:
                    sp1[0] = 1.0 - sp1[0];
                    sp1[1] = 1.0 - sp1[1];
                    break;
            }

            double[] tp1 = Burkardt.SimplexNS.Simplex.simplex_to_triangle(tv1, tv2, tv3, sp1);
            double[] rp1 = typeMethods.triangle_to_ref(tv1, tv2, tv3, tp1);
            double[] tp2 = typeMethods.ref_to_triangle(tv1, tv2, tv3, rp1);
            double[] sp2 = typeMethods.triangle_to_simplex(tv1, tv2, tv3, tp2);

            Console.WriteLine("");
            Console.WriteLine("  SP1: " + sp1[0] + "  " + sp1[1] + "");
            Console.WriteLine("  TP1: " + tp1[0] + "  " + tp1[1] + "");
            Console.WriteLine("  RP1: " + rp1[0] + "  " + rp1[1] + "");
            Console.WriteLine("  TP2: " + tp2[0] + "  " + tp2[1] + "");
            Console.WriteLine("  SP2: " + sp2[0] + "  " + sp2[1] + "");
        }
    }

    public static void test02(int degree, int numnodes, double[] vert1, double[] vert2,
            double[] vert3)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 calls TRIASYMQ for a quadrature rule of given order and region.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU GPL license.
        //
        //  Modified:
        //
        //    28 June 2014
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Hong Xiao, Zydrunas Gimbutas,
        //    A numerical algorithm for the construction of efficient quadrature
        //    rules in two and higher dimensions,
        //    Computers and Mathematics with Applications,
        //    Volume 59, 2010, pages 663-676.
        //
        //  Parameters:
        //
        //    Input, int DEGREE, the desired total polynomial degree exactness
        //    of the quadrature rule.  0 <= DEGREE <= 50.
        //
        //    Input, int NUMNODES, the number of nodes to be used by the rule.
        //
        //    Input, double VERT1[2], VERT2[2], VERT3[2], the
        //    vertices of the triangle.
        //
    {
        int j;

        Console.WriteLine("");
        Console.WriteLine("TEST02");
        Console.WriteLine("  Symmetric quadrature rule for a triangle.");
        Console.WriteLine("  Polynomial exactness degree DEGREE = " + degree + "");

        double area = typeMethods.triangle_area(vert1, vert2, vert3);
        //
        //  Retrieve and print a symmetric quadrature rule.
        //
        double[] rnodes = new double[2 * numnodes];
        double[] weights = new double[numnodes];

        SymmetricQuadrature.triasymq(degree, vert1, vert2, vert3, ref rnodes, ref weights, numnodes);

        Console.WriteLine("");
        Console.WriteLine("  NUMNODES = " + numnodes + "");

        Console.WriteLine("");
        Console.WriteLine("     J      W               X               Y");
        Console.WriteLine("");
        for (j = 0; j < numnodes; j++)
        {
            Console.WriteLine(j + "  "
                                + weights[j] + "  "
                                + rnodes[0 + j * 2] + "  "
                                + rnodes[1 + j * 2] + "");
        }

        double d = typeMethods.r8vec_sum(numnodes, weights);

        Console.WriteLine("   Sum  " + d + "");
        Console.WriteLine("  Area  " + area + "");
    }

    public static void test03(int degree, int numnodes, double[] vert1, double[] vert2,
            double[] vert3, string header)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 calls TRIASYMQ_GNUPLOT to generate graphics files.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU GPL license.
        //
        //  Modified:
        //
        //    30 June 2014
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Hong Xiao, Zydrunas Gimbutas,
        //    A numerical algorithm for the construction of efficient quadrature
        //    rules in two and higher dimensions,
        //    Computers and Mathematics with Applications,
        //    Volume 59, 2010, pages 663-676.
        //
        //  Parameters:
        //
        //    Input, int DEGREE, the desired total polynomial degree exactness
        //    of the quadrature rule.  0 <= DEGREE <= 50.
        //
        //    Input, int NUMNODES, the number of nodes to be used by the rule.
        //
        //    Input, double VERT1[2], VERT2[2], VERT3[2], the
        //    vertices of the triangle.
        //
        //    Input, string HEADER, an identifier for the graphics filenames.
        //
    {
        Console.WriteLine("");
        Console.WriteLine("TEST03");
        Console.WriteLine("  TRIASYMQ_GNUPLOT creates gnuplot graphics files.");
        Console.WriteLine("  Polynomial exactness degree DEGREE = " + degree + "");

        double[] rnodes = new double[2 * numnodes];
        double[] weights = new double[numnodes];

        SymmetricQuadrature.triasymq(degree, vert1, vert2, vert3, ref rnodes, ref weights, numnodes);

        Console.WriteLine("  Number of nodes = " + numnodes + "");

        SymmetricQuadrature.triasymq_gnuplot(vert1, vert2, vert3, numnodes, rnodes, header);

    }

    public static void test04(int degree, int numnodes, double[] vert1, double[] vert2,
            double[] vert3, string header)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 gets a rule and writes it to a file.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU GPL license.
        //
        //  Modified:
        //
        //    30 June 2014
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Hong Xiao, Zydrunas Gimbutas,
        //    A numerical algorithm for the construction of efficient quadrature
        //    rules in two and higher dimensions,
        //    Computers and Mathematics with Applications,
        //    Volume 59, 2010, pages 663-676.
        //
        //  Parameters:
        //
        //    Input, int DEGREE, the desired total polynomial degree exactness
        //    of the quadrature rule.  0 <= DEGREE <= 50.
        //
        //    Input, int NUMNODES, the number of nodes to be used by the rule.
        //
        //    Input, double VERT1[2], VERT2[2], VERT3[2], the
        //    vertices of the triangle.
        //
        //    Input, string HEADER, an identifier for the filenames.
        //
    {
        int j;
        List<string> rule_unit = new();

        Console.WriteLine("");
        Console.WriteLine("TEST04");
        Console.WriteLine("  Get a quadrature rule for a triangle.");
        Console.WriteLine("  Then write it to a file.");
        Console.WriteLine("  Polynomial exactness degree DEGREE = " + degree + "");
        //
        //  Retrieve a symmetric quadrature rule.
        //
        double[] rnodes = new double[2 * numnodes];
        double[] weights = new double[numnodes];

        SymmetricQuadrature.triasymq(degree, vert1, vert2, vert3, ref rnodes, ref weights, numnodes);
        //
        //  Write the points and weights to a file.
        //
        string rule_filename = header + ".txt";

        for (j = 0; j < numnodes; j++)
        {
            rule_unit.Add(rnodes[0 + j * 2] + "  "
                                            + rnodes[1 + j * 2] + "  "
                                            + weights[j] + "");
        }

        File.WriteAllLines(rule_filename, rule_unit);

        Console.WriteLine("");
        Console.WriteLine("  Quadrature rule written to file '" + rule_filename + "'");

    }

    public static void test05(int degree, int numnodes, double[] vert1, double[] vert2,
            double[] vert3)

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 calls TRIASYMQ for a quadrature rule of given order and region.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU GPL license.
        //
        //  Modified:
        //
        //    28 June 2014
        //
        //  Author:
        //
        //    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
        //    C++ version by John Burkardt.
        //
        //  Reference:
        //
        //    Hong Xiao, Zydrunas Gimbutas,
        //    A numerical algorithm for the construction of efficient quadrature
        //    rules in two and higher dimensions,
        //    Computers and Mathematics with Applications,
        //    Volume 59, 2010, pages 663-676.
        //
        //  Parameters:
        //
        //    Input, int DEGREE, the desired total polynomial degree 
        //    exactness of the quadrature rule.  0 <= DEGREE <= 50.
        //
        //    Input, int NUMNODES, the number of nodes to be used by the rule.
        //
        //    Input, double VERT1[2], VERT2[2], VERT3[2], the
        //    vertices of the triangle.
        //
    {
        int i;
        int j;
        double[] z = new double[2];

        Console.WriteLine("");
        Console.WriteLine("TEST05");
        Console.WriteLine("  Compute a quadrature rule for a triangle.");
        Console.WriteLine("  Check it by integrating orthonormal polynomials.");
        Console.WriteLine("  Polynomial exactness degree DEGREE = " + degree + "");

        double area = typeMethods.triangle_area(vert1, vert2, vert3);
        //
        //  Retrieve a symmetric quadrature rule.
        //
        double[] rnodes = new double[2 * numnodes];
        double[] weights = new double[numnodes];

        SymmetricQuadrature.triasymq(degree, vert1, vert2, vert3, ref rnodes, ref weights, numnodes);
        //
        //  Construct the matrix of values of the orthogonal polynomials
        //  at the user-provided nodes
        //
        int npols = (degree + 1) * (degree + 2) / 2;
        double[] rints = new double[npols];

        for (j = 0; j < npols; j++)
        {
            rints[j] = 0.0;
        }

        for (i = 0; i < numnodes; i++)
        {
            z[0] = rnodes[0 + i * 2];
            z[1] = rnodes[1 + i * 2];
            double[] r = typeMethods.triangle_to_ref(vert1, vert2, vert3, z);
            double[] pols = Orthogonal.ortho2eva(degree, r);
            for (j = 0; j < npols; j++)
            {
                rints[j] += weights[i] * pols[j];
            }
        }

        double scale = Math.Sqrt(Math.Sqrt(3.0)) / Math.Sqrt(area);

        for (j = 0; j < npols; j++)
        {
            rints[j] *= scale;
        }

        double d = Math.Pow(rints[0] - Math.Sqrt(area), 2);
        for (j = 1; j < npols; j++)
        {
            d += rints[j] * rints[j];
        }

        d = Math.Sqrt(d) / npols;

        Console.WriteLine("");
        Console.WriteLine("  RMS integration error = " + d + "");
    }
}