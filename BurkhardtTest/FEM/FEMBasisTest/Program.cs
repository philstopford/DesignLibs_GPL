using System;
using Burkardt.FEM;
using Burkardt.Types;

namespace Burkardt.FEMBasisTest
{
    class Program
    {
        static void Main(string[] args)
        {
            Console.WriteLine("");
            Console.WriteLine("FEM_BASIS_TEST:");
            Console.WriteLine("  Test the FEM_BASIS library.");

            test01();
            test02();
            test03();
            test04();
            test05();
            test06();
            test07();
            test08();

            Console.WriteLine("");
            Console.WriteLine("FEM_BASIS_TEST:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine("");
        }

        static void test01()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests FEM_BASIS_1D
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 October 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            int d;
            int i1;
            int i2;
            int j1;
            int j2;
            double lij;
            double x1;
            double x2;

            Console.WriteLine("");
            Console.WriteLine("TEST01");
            Console.WriteLine("  FEM_BASIS_1D evaluates an arbitrary");
            Console.WriteLine("  basis function over an interval.");

            i1 = 2;
            j1 = 1;
            d = i1 + j1;
            x1 = typeMethods.r8_fraction(i1, d);
            Console.WriteLine("");
            Console.WriteLine("   I   J        X      L(I,J)(X)");
            Console.WriteLine("");
            Console.WriteLine("  " + i1.ToString().PadLeft(2)
                + "  " + j1.ToString().PadLeft(2)
                + "  " + x1.ToString().PadLeft(10)
                + "  " + 1.0.ToString().PadLeft(14) + "");
            Console.WriteLine("");
            for (i2 = 0; i2 <= d; i2++)
            {
                j2 = d - i2;
                x2 = typeMethods.r8_fraction(i2, d);
                lij = FEM_basis.fem_basis_1d(i1, j1, x2);
                Console.WriteLine("  " + i2.ToString().PadLeft(2)
                    + "  " + j2.ToString().PadLeft(2)
                    + "  " + x2.ToString().PadLeft(10)
                    + "  " + lij.ToString().PadLeft(14) + "");
            }
        }

        static void test02()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests FEM_BASIS_2D
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 October 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            int d;
            int i1;
            int i2;
            int j1;
            int j2;
            int k1;
            int k2;
            double lijk;
            double x1;
            double x2;
            double y1;
            double y2;

            Console.WriteLine("");
            Console.WriteLine("TEST02");
            Console.WriteLine("  FEM_BASIS_2D evaluates an arbitrary triangular");
            Console.WriteLine("  basis function.");

            i1 = 1;
            j1 = 0;
            k1 = 2;
            d = i1 + j1 + k1;
            x1 = typeMethods.r8_fraction(i1, d);
            y1 = typeMethods.r8_fraction(j1, d);
            Console.WriteLine("");
            Console.WriteLine("   I   J   K        X           Y      L(I,J,K)(X,Y)");
            Console.WriteLine("");
            Console.WriteLine("  " + i1.ToString().PadLeft(2)
                + "  " + j1.ToString().PadLeft(2)
                + "  " + k1.ToString().PadLeft(2)
                + "  " + x1.ToString().PadLeft(10)
                + "  " + y1.ToString().PadLeft(10)
                + "  " + 1.0.ToString().PadLeft(14) + "");
            Console.WriteLine("");
            for (j2 = 0; j2 <= d; j2++)
            {
                for (i2 = 0; i2 <= d - j2; i2++)
                {
                    k2 = d - i2 - j2;
                    x2 = typeMethods.r8_fraction(i2, d);
                    y2 =  typeMethods.r8_fraction(j2, d);
                    lijk = FEM_basis.fem_basis_2d(i1, j1, k1, x2, y2);
                    Console.WriteLine("  " + i2.ToString().PadLeft(2)
                        + "  " + j2.ToString().PadLeft(2)
                        + "  " + k2.ToString().PadLeft(2)
                        + "  " + x2.ToString().PadLeft(10)
                        + "  " + y2.ToString().PadLeft(10)
                        + "  " + lijk.ToString().PadLeft(14) + "");
                }
            }
        }

        static void test03()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests FEM_BASIS_3D
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    05 October 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            int d;
            int i1;
            int i2;
            int j1;
            int j2;
            int k1;
            int k2;
            int l1;
            int l2;
            double lijkl;
            double x1;
            double x2;
            double y1;
            double y2;
            double z1;
            double z2;

            Console.WriteLine("");
            Console.WriteLine("TEST03");
            Console.WriteLine("  FEM_BASIS_3D evaluates an arbitrary tetrahedral");
            Console.WriteLine("  basis function.");

            i1 = 1;
            j1 = 0;
            k1 = 2;
            l1 = 1;
            d = i1 + j1 + k1 + l1;
            x1 = typeMethods.r8_fraction(i1, d);
            y1 = typeMethods.r8_fraction(j1, d);
            z1 = typeMethods.r8_fraction(k1, d);
            Console.WriteLine("");
            Console.WriteLine("   I   J   K   L        X           Y           Z      L(I,J,K,L)(X,Y,Z)");
            Console.WriteLine("");
            Console.WriteLine("  " + i1.ToString().PadLeft(2)
                + "  " + j1.ToString().PadLeft(2)
                + "  " + k1.ToString().PadLeft(2)
                + "  " + l1.ToString().PadLeft(2)
                + "  " + x1.ToString().PadLeft(10)
                + "  " + y1.ToString().PadLeft(10)
                + "  " + z1.ToString().PadLeft(10)
                + "  " + 1.0.ToString().PadLeft(14) + "");
            Console.WriteLine("");
            for (k2 = 0; k2 <= d; k2++)
            {
                for (j2 = 0; j2 <= d - k2; j2++)
                {
                    for (i2 = 0; i2 <= d - j2 - k2; i2++)
                    {
                        l2 = d - i2 - j2 - k2;
                        x2 = typeMethods.r8_fraction(i2, d);
                        y2 = typeMethods.r8_fraction(j2, d);
                        z2 = typeMethods.r8_fraction(k2, d);
                        lijkl = FEM_basis.fem_basis_3d(i1, j1, k1, l1, x2, y2, z2);
                        Console.WriteLine("  " + i2.ToString().PadLeft(2)
                            + "  " + j2.ToString().PadLeft(2)
                            + "  " + k2.ToString().PadLeft(2)
                            + "  " + l2.ToString().PadLeft(2)
                            + "  " + x2.ToString().PadLeft(10)
                            + "  " + y2.ToString().PadLeft(10)
                            + "  " + z2.ToString().PadLeft(10)
                            + "  " + lijkl.ToString().PadLeft(14) + "");
                    }
                }
            }
        }

        static void test04()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 tests FEM_BASIS_MD, repeating TEST01.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 October 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            int d;
            int i;
            int[] i1 = new int[2];
            int[] i2 = new int[2];
            double l;
            int m = 1;
            int p1;
            double[] x1 = new double[1];
            double[] x2 = new double[1];

            Console.WriteLine("");
            Console.WriteLine("TEST04");
            Console.WriteLine("  FEM_BASIS_MD evaluates an arbitrary");
            Console.WriteLine("  basis function over an M-dimensional simplex.");

            i1[0] = 2;
            i1[1] = 1;
            d = typeMethods.i4vec_sum(m + 1, i1);
            for (i = 0; i < m; i++)
            {
                x1[i] = typeMethods.r8_fraction(i1[i], d);
            }

            Console.WriteLine("");
            Console.WriteLine("   I   J        X      L(I,J)(X)");
            Console.WriteLine("");
            Console.WriteLine("  " + i1[0].ToString().PadLeft(2)
                + "  " + i1[1].ToString().PadLeft(2)
                + "  " + x1[0].ToString().PadLeft(2)
                + "  " + 1.0.ToString().PadLeft(2) + "");
            Console.WriteLine("");
            for (p1 = 0; p1 <= d; p1++)
            {
                i2[0] = p1;
                i2[1] = d - i2[0];
                for (i = 0; i < m; i++)
                {
                    x2[i] = typeMethods.r8_fraction(i2[i], d);
                }

                l = FEM_basis.fem_basis_md(m, i1, x2);
                Console.WriteLine("  " + i2[0].ToString().PadLeft(2)
                    + "  " + i2[1].ToString().PadLeft(2)
                    + "  " + x2[0].ToString().PadLeft(2)
                    + "  " + l.ToString().PadLeft(2) + "");
            }
        }

        static void test05()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 tests FEM_BASIS_MD, repeating TEST02.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 October 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            int d;
            int i;
            int[] i1 = new int[3];
            int[] i2 = new int[3];
            double l;
            int m = 2;
            int p1;
            int p2;
            double[] x1 = new double[2];
            double[] x2 = new double[2];

            Console.WriteLine("");
            Console.WriteLine("TEST05");
            Console.WriteLine("  FEM_BASIS_MD evaluates an arbitrary");
            Console.WriteLine("  basis function over an M-dimensional simplex.");

            i1[0] = 1;
            i1[1] = 0;
            i1[2] = 2;
            d = typeMethods.i4vec_sum(m + 1, i1);
            for (i = 0; i < m; i++)
            {
                x1[i] = typeMethods.r8_fraction(i1[i], d);
            }

            Console.WriteLine("");
            Console.WriteLine("   I   J   K        X           Y      L(I,J,K)(X,Y)");
            Console.WriteLine("");
            Console.WriteLine("  " + i1[0].ToString().PadLeft(2)
                + "  " + i1[1].ToString().PadLeft(2)
                + "  " + i1[2].ToString().PadLeft(2)
                + "  " + x1[0].ToString().PadLeft(10)
                + "  " + x1[1].ToString().PadLeft(10)
                + "  " + 1.0.ToString().PadLeft(14) + "");
            Console.WriteLine("");
            for (p2 = 0; p2 <= d; p2++)
            {
                i2[1] = p2;
                for (p1 = 0; p1 <= d - p2; p1++)
                {
                    i2[0] = p1;
                    i2[2] = d - i2[0] - i2[1];
                    for (i = 0; i < m; i++)
                    {
                        x2[i] = typeMethods.r8_fraction(i2[i], d);
                    }

                    l = FEM_basis.fem_basis_md(m, i1, x2);
                    Console.WriteLine("  " + i2[0].ToString().PadLeft(2)
                        + "  " + i2[1].ToString().PadLeft(2)
                        + "  " + i2[2].ToString().PadLeft(2)
                        + "  " + x2[0].ToString().PadLeft(10)
                        + "  " + x2[1].ToString().PadLeft(10)
                        + "  " + l.ToString().PadLeft(14) + "");
                }
            }
        }

        static void test06()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST06 tests FEM_BASIS_MD, repeating TEST03.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    10 October 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            int d;
            int i;
            int[] i1 = new int[4];
            int[] i2 = new int[4];
            double l;
            int m = 3;
            int p1;
            int p2;
            int p3;
            double[] x1 = new double[3];
            double[] x2 = new double[3];

            Console.WriteLine("");
            Console.WriteLine("TEST06");
            Console.WriteLine("  FEM_BASIS_MD evaluates an arbitrary");
            Console.WriteLine("  basis function over an M-dimensional simplex.");

            i1[0] = 1;
            i1[1] = 0;
            i1[2] = 2;
            i1[3] = 1;
            d = typeMethods.i4vec_sum(m + 1, i1);
            for (i = 0; i < m; i++)
            {
                x1[i] = typeMethods.r8_fraction(i1[i], d);
            }

            Console.WriteLine("");
            Console.WriteLine("   I   J   K   L        X           Y           Z      L(I,J,K,L)(X,Y,Z)");
            Console.WriteLine("");
            Console.WriteLine("  " + i1[0].ToString().PadLeft(2)
                + "  " + i1[1].ToString().PadLeft(2)
                + "  " + i1[2].ToString().PadLeft(2)
                + "  " + i1[3].ToString().PadLeft(2)
                + "  " + x1[0].ToString().PadLeft(10)
                + "  " + x1[1].ToString().PadLeft(10)
                + "  " + x1[2].ToString().PadLeft(10)
                + "  " + 1.0.ToString().PadLeft(14) + "");
            Console.WriteLine("");
            for (p3 = 0; p3 <= d; p3++)
            {
                i2[2] = p3;
                for (p2 = 0; p2 <= d - p3; p2++)
                {
                    i2[1] = p2;
                    for (p1 = 0; p1 <= d - p3 - p2; p1++)
                    {
                        i2[0] = p1;
                        i2[3] = d - i2[0] - i2[1] - i2[2];
                        for (i = 0; i < m; i++)
                        {
                            x2[i] = typeMethods.r8_fraction(i2[i], d);
                        }

                        l = FEM_basis.fem_basis_md(m, i1, x2);
                        Console.WriteLine("  " + i2[0].ToString().PadLeft(2)
                            + "  " + i2[1].ToString().PadLeft(2)
                            + "  " + i2[2].ToString().PadLeft(2)
                            + "  " + i2[3].ToString().PadLeft(2)
                            + "  " + x2[0].ToString().PadLeft(10)
                            + "  " + x2[1].ToString().PadLeft(10)
                            + "  " + x2[2].ToString().PadLeft(10)
                            + "  " + l.ToString().PadLeft(14) + "");
                    }
                }
            }
        }

        static void test07()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST07 tests FEM_BASIS_PRISM_TRIANGLE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 October 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            double b;
            int di;
            int dj;
            int[] i1 =  {
                2, 0, 0
            }
            ;
            int[] i2 = new int[3];
            int i_0;
            int i_1;
            int[] j1 =  {
                1, 1
            }
            ;
            int[] j2 = new int[2];
            int j_0;
            double[] xyz1 = new double[3];
            double[] xyz2 = new double[3];

            Console.WriteLine("");
            Console.WriteLine("TEST07");
            Console.WriteLine("  FEM_BASIS_PRISM_TRIANGLE evaluates an arbitrary");
            Console.WriteLine("  basis function over a right triangular prism.");
            Console.WriteLine("");
            Console.WriteLine("  Here, we generate basis functions which can be");
            Console.WriteLine("  up to degree 2 in X and Y, and up to degree 2 in Z.");
            Console.WriteLine("");
            Console.WriteLine("  Choose a node N1, define the basis function associated");
            Console.WriteLine("  with that node, and then evaluate it at all other nodes.");

            di = typeMethods.i4vec_sum(3, i1);
            xyz1[0] = typeMethods.r8_fraction(i1[0], di);
            xyz1[1] = typeMethods.r8_fraction(i1[1], di);

            dj = typeMethods.i4vec_sum(2, j1);
            xyz1[2] = typeMethods.r8_fraction(j1[0], dj);

            Console.WriteLine("");
            Console.WriteLine("  I1  I2  I3  J1  J2        X           Y           Z          B(X,Y,Z)");
            Console.WriteLine("");
            Console.WriteLine("  " + i1[0].ToString().PadLeft(2)
                + "  " + i1[1].ToString().PadLeft(2)
                + "  " + i1[2].ToString().PadLeft(2)
                + "  " + j1[0].ToString().PadLeft(2)
                + "  " + j1[1].ToString().PadLeft(2)
                + "  " + xyz1[0].ToString().PadLeft(10)
                + "  " + xyz1[1].ToString().PadLeft(10)
                + "  " + xyz1[2].ToString().PadLeft(10)
                + "  " + 1.0.ToString().PadLeft(10) + "");

            Console.WriteLine("");

            for (i_0 = 0; i_0 <= di; i_0++)
            {
                i2[0] = i_0;
                xyz2[0] = typeMethods.r8_fraction(i2[0], di);
                for (i_1 = 0; i_1 <= di - i2[0]; i_1++)
                {
                    i2[1] = i_1;
                    xyz2[1] = typeMethods.r8_fraction(i2[1], di);
                    i2[2] = di - i2[0] - i2[1];
                    for (j_0 = 0; j_0 <= dj; j_0++)
                    {
                        j2[0] = j_0;
                        j2[1] = dj - j2[0];
                        xyz2[2] = typeMethods.r8_fraction(j2[0], dj);

                        b = FEM_basis.fem_basis_prism_triangle(i1, j1, xyz2);

                        Console.WriteLine("  " + i2[0].ToString().PadLeft(2)
                            + "  " + i2[1].ToString().PadLeft(2)
                            + "  " + i2[2].ToString().PadLeft(2)
                            + "  " + j2[0].ToString().PadLeft(2)
                            + "  " + j2[1].ToString().PadLeft(2)
                            + "  " + xyz2[0].ToString().PadLeft(10)
                            + "  " + xyz2[1].ToString().PadLeft(10)
                            + "  " + xyz2[2].ToString().PadLeft(10)
                            + "  " + b.ToString().PadLeft(10) + "");
                    }
                }
            }
        }

        static void test08()
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST08 tests FEM_BASIS_PRISM_TRIANGLE.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    20 October 2010
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            double b;
            int di;
            int dj;
            int[] i1 =  {
                2, 0, 1
            }
            ;
            int[] i2 = new int[3];
            int i_0;
            int i_1;
            int[] j1 =  {
                1, 0
            }
            ;
            int[] j2 = new int[2];
            int j_0;
            double[] xyz1 = new double[3];
            double[] xyz2 = new double[3];

            Console.WriteLine("");
            Console.WriteLine("TEST08");
            Console.WriteLine("  FEM_BASIS_PRISM_TRIANGLE evaluates an arbitrary");
            Console.WriteLine("  basis function over a right triangular prism.");
            Console.WriteLine("");
            Console.WriteLine("  Here, we generate basis functions which can be");
            Console.WriteLine("  up to degree 3 in X and Y, and up to degree 1 in Z.");
            Console.WriteLine("");
            Console.WriteLine("  Choose a node N1, define the basis function associated");
            Console.WriteLine("  with that node, and then evaluate it at all other nodes.");

            di = typeMethods.i4vec_sum(3, i1);
            xyz1[0] = typeMethods.r8_fraction(i1[0], di);
            xyz1[1] = typeMethods.r8_fraction(i1[1], di);

            dj = typeMethods.i4vec_sum(2, j1);
            xyz1[2] = typeMethods.r8_fraction(j1[0], dj);

            Console.WriteLine("");
            Console.WriteLine("  I1  I2  I3  J1  J2        X           Y           Z          B(X,Y,Z)");
            Console.WriteLine("");
            Console.WriteLine("  " + i1[0].ToString().PadLeft(2)
                + "  " + i1[1].ToString().PadLeft(2)
                + "  " + i1[2].ToString().PadLeft(2)
                + "  " + j1[0].ToString().PadLeft(2)
                + "  " + j1[1].ToString().PadLeft(2)
                + "  " + xyz1[0].ToString().PadLeft(10)
                + "  " + xyz1[1].ToString().PadLeft(10)
                + "  " + xyz1[2].ToString().PadLeft(10)
                + "  " + 1.0.ToString().PadLeft(10) + "");
            Console.WriteLine("");

            for (i_0 = 0; i_0 <= di; i_0++)
            {
                i2[0] = i_0;
                xyz2[0] = typeMethods.r8_fraction(i2[0], di);
                for (i_1 = 0; i_1 <= di - i2[0]; i_1++)
                {
                    i2[1] = i_1;
                    xyz2[1] = typeMethods.r8_fraction(i2[1], di);
                    i2[2] = di - i2[0] - i2[1];
                    for (j_0 = 0; j_0 <= dj; j_0++)
                    {
                        j2[0] = j_0;
                        j2[1] = dj - j2[0];
                        xyz2[2] = typeMethods.r8_fraction(j2[0], dj);

                        b = FEM_basis.fem_basis_prism_triangle(i1, j1, xyz2);

                        Console.WriteLine("  " + i2[0].ToString().PadLeft(2)
                            + "  " + i2[1].ToString().PadLeft(2)
                            + "  " + i2[2].ToString().PadLeft(2)
                            + "  " + j2[0].ToString().PadLeft(2)
                            + "  " + j2[1].ToString().PadLeft(2)
                            + "  " + xyz2[0].ToString().PadLeft(10)
                            + "  " + xyz2[1].ToString().PadLeft(10)
                            + "  " + xyz2[2].ToString().PadLeft(10)
                            + "  " + b.ToString().PadLeft(10) + "");
                    }
                }
            }
        }
    }
}