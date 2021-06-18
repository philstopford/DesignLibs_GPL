using System;
using Burkardt.Blend;

namespace BlendTest
{
    partial class Program
    {
        static void Main(string[] args)
        {
            Console.WriteLine();
            Console.WriteLine("BLEND_TEST");
            Console.WriteLine("  Test the BLEND library.");

            test01 ( );
            test02 ( );
            test03 ( );
            test04 ( );
            test05 ( );
            test06 ( );
            test07 ( );
            test08 ( );
            test09 ( );

            test10 ( );
            test11 ( );
            test12 ( );
            blend_i_0d1_test ( );
            blend_ij_0d1_test ( );
            blend_ij_1d1_test ( );
            blend_ijk_0d1_test ( );
            blend_ijk_1d1_test ( );
            blend_ijk_2d1_test ( );

            Console.WriteLine();
            Console.WriteLine("BLEND_TEST:");
            Console.WriteLine("  Normal end of execution.");
            Console.WriteLine();
        }
        
        
        static void test01 ( )

        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST01 tests BLEND_R_0DN on the identity map.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 October 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            double[] x = new double[1];

            Console.WriteLine();
            Console.WriteLine("TEST01");
            Console.WriteLine("  Identity test on BLEND_R_0DN.");

            int n = 1;

            double r = 0.0;
            x = Blend.blend_r_0dn ( r, x, n, identity_r );
            string cout = "  ";
            string t = r.ToString().PadLeft(8) + "  ";
            cout += t;
            cout += x[0].ToString().PadLeft(8);
            Console.WriteLine(cout);

            r = 1.0;
            x = Blend.blend_r_0dn ( r, x, n, identity_r );
            cout = "  ";
            t = r.ToString().PadLeft(8) + "  ";
            cout += t;
            cout += x[0].ToString().PadLeft(8);
            Console.WriteLine(cout);

            r = 0.5;
            x = Blend.blend_r_0dn ( r, x, n, identity_r );
            cout = "  ";
            t = r.ToString().PadLeft(8) + "  ";
            cout += t;
            cout += x[0].ToString().PadLeft(8);
            Console.WriteLine(cout);
        }

        
        static void test02 ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST02 tests BLEND_RS_0DN on the identity map.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 October 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            double[] x = new double[2];

            Console.WriteLine();
            Console.WriteLine("TEST02");
            Console.WriteLine("  Identity test on BLEND_RS_0DN.");

            int n = 2;

            double r = 0.0;
            double s = 0.0;
            x = Blend.blend_rs_0dn ( r, s, x, n, identity_rs );
            string cout =  "  ";
            string t = r.ToString().PadLeft(8) + "  ";
            cout += t;
            t = s.ToString().PadLeft(8) + "  ";
            cout += t;
            t = x[0].ToString().PadLeft(8) + "  ";
            cout += t;
            t = x[1].ToString().PadLeft(8);
            cout += t.ToString().PadLeft(8);
            Console.WriteLine(cout);

            r = 1.0;
            s = 0.0;
            x = Blend.blend_rs_0dn ( r, s, x, n, identity_rs );
            cout =  "  ";
            t = r.ToString().PadLeft(8) + "  ";
            cout += t;
            t = s.ToString().PadLeft(8) + "  ";
            cout += t;
            t = x[0].ToString().PadLeft(8) + "  ";
            cout += t;
            t = x[1].ToString().PadLeft(8);
            cout += t;
            Console.WriteLine(cout);

            r = 0.0;
            s = 1.0;
            x = Blend.blend_rs_0dn ( r, s, x, n, identity_rs );
            cout =  "  ";
            t = r.ToString().PadLeft(8) + "  ";
            cout += t;
            t = s.ToString().PadLeft(8) + "  ";
            cout += t;
            t = x[0].ToString().PadLeft(8) + "  ";
            cout += t;
            t = x[1].ToString().PadLeft(8);
            cout += t;
            Console.WriteLine(cout);

            r = 1.0;
            s = 1.0;
            x = Blend.blend_rs_0dn ( r, s, x, n, identity_rs );
            cout =  "  ";
            t = r.ToString().PadLeft(8) + "  ";
            cout += t;
            t = s.ToString().PadLeft(8) + "  ";
            cout += t;
            t = x[0].ToString().PadLeft(8) + "  ";
            cout += t;
            t = x[1].ToString().PadLeft(8);
            cout += t;
            Console.WriteLine(cout);

            r = 0.5;
            s = 0.5;
            x = Blend.blend_rs_0dn ( r, s, x, n, identity_rs );
            cout =  "  ";
            t = r.ToString().PadLeft(8) + "  ";
            cout += t;
            t = s.ToString().PadLeft(8) + "  ";
            cout += t;
            t = x[0].ToString().PadLeft(8) + "  ";
            cout += t;
            t = x[1].ToString().PadLeft(8);
            cout += t;
            Console.WriteLine(cout);
        }


        static void test03 ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST03 tests BLEND_RS_1DN on the identity map.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 October 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            double[] x = new double[2];

            Console.WriteLine();
            Console.WriteLine("TEST03");
            Console.WriteLine("  Identity test on BLEND_RS_1DN.");

            int n = 2;

            double r = 0.0;
            double s = 0.0;
            x = Blend.blend_rs_1dn ( r, s, x, n, identity_rs );
            string cout =  "  ";
            string t = r.ToString().PadLeft(8) + "  ";
            cout += t;
            t = s.ToString().PadLeft(8) + "  ";
            cout += t;
            t = x[0].ToString().PadLeft(8) + "  ";
            cout += t;
            t = x[1].ToString().PadLeft(8);
            cout += t;
            Console.WriteLine(cout);

            r = 1.0;
            s = 0.0;
            x = Blend.blend_rs_1dn ( r, s, x, n, identity_rs );
            cout =  "  ";
            t = r.ToString().PadLeft(8) + "  ";
            cout += t;
            t = s.ToString().PadLeft(8) + "  ";
            cout += t;
            t = x[0].ToString().PadLeft(8) + "  ";
            cout += t;
            t = x[1].ToString().PadLeft(8);
            cout += t;
            Console.WriteLine(cout);

            r = 0.0;
            s = 1.0;
            x = Blend.blend_rs_1dn ( r, s, x, n, identity_rs );
            cout =  "  ";
            t = r.ToString().PadLeft(8) + "  ";
            cout += t;
            t = s.ToString().PadLeft(8) + "  ";
            cout += t;
            t = x[0].ToString().PadLeft(8) + "  ";
            cout += t;
            t = x[1].ToString().PadLeft(8);
            cout += t;
            Console.WriteLine(cout);

            r = 1.0;
            s = 1.0;
            x = Blend.blend_rs_1dn ( r, s, x, n, identity_rs );
            cout =  "  ";
            t = r.ToString().PadLeft(8) + "  ";
            cout += t;
            t = s.ToString().PadLeft(8) + "  ";
            cout += t;
            t = x[0].ToString().PadLeft(8) + "  ";
            cout += t;
            t = x[1].ToString().PadLeft(8);
            cout += t;
            Console.WriteLine(cout);

            r = 0.5;
            s = 0.5;
            x = Blend.blend_rs_1dn ( r, s, x, n, identity_rs );
            cout =  "  ";
            t = r.ToString().PadLeft(8) + "  ";
            cout += t;
            t = s.ToString().PadLeft(8) + "  ";
            cout += t;
            t = x[0].ToString().PadLeft(8) + "  ";
            cout += t;
            t = x[1].ToString().PadLeft(8);
            cout += t;
            Console.WriteLine(cout);
        }

        
        static void test04 ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST04 tests BLEND_RST_0DN on the identity map.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 October 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            double[] x = new double[3];

            Console.WriteLine();
            Console.WriteLine("TEST04");
            Console.WriteLine("  Identity test on BLEND_RST_0DN.");

            int n = 3;

            double r = 0.0;
            double s = 0.0;
            double t = 0.0;
            x = Blend.blend_rst_0dn ( r, s, t, x, n, identity_rst );
            string cout =  "  ";
            string ts = r.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = s.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = t.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[0].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[1].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[2].ToString().PadLeft(8);
            cout += ts;
            Console.WriteLine(cout);

            r = 1.0;
            s = 0.0;
            t = 0.0;
            x = Blend.blend_rst_0dn ( r, s, t, x, n, identity_rst );
            cout =  "  ";
            ts = r.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = s.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = t.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[0].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[1].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[2].ToString().PadLeft(8);
            cout += ts;
            Console.WriteLine(cout);

            r = 0.0;
            s = 1.0;
            t = 0.0;
            x = Blend.blend_rst_0dn ( r, s, t, x, n, identity_rst );
            cout =  "  ";
            ts = r.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = s.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = t.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[0].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[1].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[2].ToString().PadLeft(8);
            cout += ts;
            Console.WriteLine(cout);

            r = 0.0;
            s = 0.0;
            t = 1.0;
            x = Blend.blend_rst_0dn ( r, s, t, x, n, identity_rst );
            cout =  "  ";
            ts = r.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = s.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = t.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[0].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[1].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[2].ToString().PadLeft(8);
            cout += ts;
            Console.WriteLine(cout);

            r = 1.0;
            s = 1.0;
            t = 1.0;
            x = Blend.blend_rst_0dn ( r, s, t, x, n, identity_rst );
            cout =  "  ";
            ts = r.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = s.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = t.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[0].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[1].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[2].ToString().PadLeft(8);
            cout += ts;
            Console.WriteLine(cout);

            r = 0.5;
            s = 0.5;
            t = 0.5;
            x = Blend.blend_rst_0dn ( r, s, t, x, n, identity_rst );
            cout =  "  ";
            ts = r.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = s.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = t.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[0].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[1].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[2].ToString().PadLeft(8);
            cout += ts;
            Console.WriteLine(cout);
        }

        
        static void test05 ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST05 tests BLEND_RST_1DN on the identity map.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 October 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            double[] x = new double[3];

            Console.WriteLine();
            Console.WriteLine("TEST05");
            Console.WriteLine("  Identity test on BLEND_RST_1DN.");

            int n = 3;

            double r = 0.0;
            double s = 0.0;
            double t = 0.0;
            x = Blend.blend_rst_1dn ( r, s, t, x, n, identity_rst );
            string cout =  "  ";
            string ts = r.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = s.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = t.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[0].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[1].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[2].ToString().PadLeft(8);
            cout += ts;
            Console.WriteLine(cout);

            r = 1.0;
            s = 0.0;
            t = 0.0;
            x = Blend.blend_rst_1dn ( r, s, t, x, n, identity_rst );
            cout =  "  ";
            ts = r.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = s.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = t.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[0].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[1].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[2].ToString().PadLeft(8);
            cout += ts;
            Console.WriteLine(cout);

            r = 0.0;
            s = 1.0;
            t = 0.0;
            x = Blend.blend_rst_1dn ( r, s, t, x, n, identity_rst );
            cout =  "  ";
            ts = r.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = s.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = t.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[0].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[1].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[2].ToString().PadLeft(8);
            cout += ts;
            Console.WriteLine(cout);

            r = 0.0;
            s = 0.0;
            t = 1.0;
            x = Blend.blend_rst_1dn ( r, s, t, x, n, identity_rst );
            cout =  "  ";
            ts = r.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = s.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = t.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[0].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[1].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[2].ToString().PadLeft(8);
            cout += ts;
            Console.WriteLine(cout);

            r = 1.0;
            s = 1.0;
            t = 1.0;
            x = Blend.blend_rst_1dn ( r, s, t, x, n, identity_rst );
            cout =  "  ";
            ts = r.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = s.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = t.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[0].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[1].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[2].ToString().PadLeft(8);
            cout += ts;
            Console.WriteLine(cout);

            r = 0.5;
            s = 0.5;
            t = 0.5;
            x = Blend.blend_rst_1dn ( r, s, t, x, n, identity_rst );
            cout =  "  ";
            ts = r.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = s.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = t.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[0].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[1].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[2].ToString().PadLeft(8);
            cout += ts;
            Console.WriteLine(cout);
        }

        
        static void test06 ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST06 tests BLEND_RST_2DN on the identity map.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 October 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            
            double[] x = new double[3];

            Console.WriteLine();
            Console.WriteLine("TEST06");
            Console.WriteLine("  Identity test on BLEND_RST_2DN.");

            int n = 3;

            double r = 0.0;
            double s = 0.0;
            double t = 0.0;
            x = Blend.blend_rst_2dn ( r, s, t, x, n, identity_rst );
            string cout =  "  ";
            string ts = r.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = s.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = t.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[0].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[1].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[2].ToString().PadLeft(8);
            cout += ts.PadLeft(8);
            Console.WriteLine(cout);
            
            r = 1.0;
            s = 0.0;
            t = 0.0;
            x = Blend.blend_rst_2dn ( r, s, t, x, n, identity_rst );
            cout =  "  ";
            ts = r.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = s.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = t.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[0].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[1].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[2].ToString().PadLeft(8);
            cout += ts;
            Console.WriteLine(cout);

            r = 0.0;
            s = 1.0;
            t = 0.0;
            x = Blend.blend_rst_2dn ( r, s, t, x, n, identity_rst );
            cout =  "  ";
            ts = r.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = s.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = t.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[0].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[1].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[2].ToString().PadLeft(8);
            cout += ts;
            Console.WriteLine(cout);

            r = 0.0;
            s = 0.0;
            t = 1.0;
            x = Blend.blend_rst_2dn ( r, s, t, x, n, identity_rst );
            cout =  "  ";
            ts = r.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = s.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = t.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[0].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[1].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[2].ToString().PadLeft(8);
            cout += ts;
            Console.WriteLine(cout);

            r = 1.0;
            s = 1.0;
            t = 1.0;
            x = Blend.blend_rst_2dn ( r, s, t, x, n, identity_rst );
            cout =  "  ";
            ts = r.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = s.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = t.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[0].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[1].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[2].ToString().PadLeft(8);
            cout += ts;
            Console.WriteLine(cout);

            r = 0.5;
            s = 0.5;
            t = 0.5;
            x = Blend.blend_rst_2dn ( r, s, t, x, n, identity_rst );
            cout =  "  ";
            ts = r.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = s.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = t.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[0].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[1].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[2].ToString().PadLeft(8);
            cout += ts;
            Console.WriteLine(cout);
        }
        
        
        static void test07 ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST07 tests BLEND_R_0DN on the stretch map.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 October 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            double[] x = new double[1];

            Console.WriteLine();
            Console.WriteLine("TEST07");
            Console.WriteLine("  Stretch test on BLEND_R_0DN.");

            int n = 1;

            double r = 0.0;
            x = Blend.blend_r_0dn ( r, x, n, stretch_r );
            string cout =  "  ";
            string ts = r.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[0].ToString().PadLeft(8) + "  ";
            cout += ts;
            Console.WriteLine(cout);

            r = 1.0;
            x = Blend.blend_r_0dn ( r, x, n, stretch_r );
            cout =  "  ";
            ts = r.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[0].ToString().PadLeft(8) + "  ";
            cout += ts;
            Console.WriteLine(cout);

            r = 0.5;
            x = Blend.blend_r_0dn ( r, x, n, stretch_r );
            cout =  "  ";
            ts = r.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[0].ToString().PadLeft(8) + "  ";
            cout += ts;
            Console.WriteLine(cout);
        }

        
        static void test08 ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST08 tests BLEND_RS_0DN on the stretch map.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 October 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            double[] x = new double[2];

            Console.WriteLine();
            Console.WriteLine("TEST08");
            Console.WriteLine("  Stretch test on BLEND_RS_0DN.");

            int n = 2;

            double r = 0.0;
            double s = 0.0;
            x = Blend.blend_rs_0dn ( r, s, x, n, stretch_rs );
            string cout =  "  ";
            string ts = r.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = s.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[0].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[1].ToString().PadLeft(8);
            cout += ts;
            Console.WriteLine(cout);

            r = 1.0;
            s = 0.0;
            x = Blend.blend_rs_0dn ( r, s, x, n, stretch_rs );
            cout =  "  ";
            ts = r.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = s.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[0].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[1].ToString().PadLeft(8);
            cout += ts;
            Console.WriteLine(cout);

            r = 0.0;
            s = 1.0;
            x = Blend.blend_rs_0dn ( r, s, x, n, stretch_rs );
            cout =  "  ";
            ts = r.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = s.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[0].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[1].ToString().PadLeft(8);
            cout += ts;
            Console.WriteLine(cout);

            r = 1.0;
            s = 1.0;
            x = Blend.blend_rs_0dn ( r, s, x, n, stretch_rs );
            cout =  "  ";
            ts = r.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = s.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[0].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[1].ToString().PadLeft(8);
            cout += ts;
            Console.WriteLine(cout);

            r = 0.5;
            s = 0.5;
            x = Blend.blend_rs_0dn ( r, s, x, n, stretch_rs );
            cout =  "  ";
            ts = r.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = s.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[0].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[1].ToString().PadLeft(8);
            cout += ts;
            Console.WriteLine(cout);
        }


        static void test09 ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST09 tests BLEND_RS_1DN on the stretch map.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 October 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            double[] x = new double[2];

            Console.WriteLine();
            Console.WriteLine("TEST09");
            Console.WriteLine("  Stretch test on BLEND_RS_1DN.");
            
            int n = 2;

            double r = 0.0;
            double s = 0.0;
            x = Blend.blend_rs_1dn ( r, s, x, n, stretch_rs );
            string cout =  "  ";
            string ts = r.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = s.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[0].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[1].ToString().PadLeft(8);
            cout += ts;
            Console.WriteLine(cout);

            r = 1.0;
            s = 0.0;
            x = Blend.blend_rs_1dn ( r, s, x, n, stretch_rs );
            cout =  "  ";
            ts = r.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = s.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[0].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[1].ToString().PadLeft(8);
            cout += ts;
            Console.WriteLine(cout);

            r = 0.0;
            s = 1.0;
            x = Blend.blend_rs_1dn ( r, s, x, n, stretch_rs );
            cout =  "  ";
            ts = r.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = s.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[0].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[1].ToString().PadLeft(8);
            cout += ts;
            Console.WriteLine(cout);

            r = 1.0;
            s = 1.0;
            x = Blend.blend_rs_1dn ( r, s, x, n, stretch_rs );
            cout =  "  ";
            ts = r.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = s.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[0].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[1].ToString().PadLeft(8);
            cout += ts;
            Console.WriteLine(cout);

            r = 0.5;
            s = 0.5;
            x = Blend.blend_rs_1dn ( r, s, x, n, stretch_rs );
            cout =  "  ";
            ts = r.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = s.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[0].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[1].ToString().PadLeft(8);
            cout += ts;
            Console.WriteLine(cout);
        }
        
        
        static void test10 ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST10 tests BLEND_RST_0DN on the stretch map.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 October 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            double[] x = new double[3];

            Console.WriteLine();
            Console.WriteLine("TEST10");
            Console.WriteLine("  Stretch test on BLEND_RST_0DN.");

            int n = 3;

            double r = 0.0;
            double s = 0.0;
            double t = 0.0;
            x = Blend.blend_rst_0dn ( r, s, t, x, n, stretch_rst );
            string cout =  "  ";
            string ts = r.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = s.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = t.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[0].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[1].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[2].ToString().PadLeft(8);
            cout += ts;
            Console.WriteLine(cout);
            
            r = 1.0;
            s = 0.0;
            t = 0.0;
            x = Blend.blend_rst_0dn ( r, s, t, x, n, stretch_rst );
            cout =  "  ";
            ts = r.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = s.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = t.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[0].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[1].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[2].ToString().PadLeft(8);
            cout += ts;
            Console.WriteLine(cout);

            r = 0.0;
            s = 1.0;
            t = 0.0;
            x = Blend.blend_rst_0dn ( r, s, t, x, n, stretch_rst );
            cout =  "  ";
            ts = r.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = s.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = t.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[0].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[1].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[2].ToString().PadLeft(8);
            cout += ts;
            Console.WriteLine(cout);

            r = 0.0;
            s = 0.0;
            t = 1.0;
            x = Blend.blend_rst_0dn ( r, s, t, x, n, stretch_rst );
            cout =  "  ";
            ts = r.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = s.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = t.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[0].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[1].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[2].ToString().PadLeft(8);
            cout += ts;
            Console.WriteLine(cout);

            r = 1.0;
            s = 1.0;
            t = 1.0;
            x = Blend.blend_rst_0dn ( r, s, t, x, n, stretch_rst );
            cout =  "  ";
            ts = r.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = s.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = t.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[0].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[1].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[2].ToString().PadLeft(8);
            cout += ts;
            Console.WriteLine(cout);

            r = 0.5;
            s = 0.5;
            t = 0.5;
            x = Blend.blend_rst_0dn ( r, s, t, x, n, stretch_rst );
            cout =  "  ";
            ts = r.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = s.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = t.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[0].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[1].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[2].ToString().PadLeft(8);
            cout += ts;
            Console.WriteLine(cout);
        }

        
        static void test11 ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST11 tests BLEND_RST_1DN on the stretch map.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 October 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            double[] x = new double[3];

            Console.WriteLine();
            Console.WriteLine("TEST11");
            Console.WriteLine("  Stretch test on BLEND_RST_1DN.");

            int n = 3;

            double r = 0.0;
            double s = 0.0;
            double t = 0.0;
            x = Blend.blend_rst_1dn ( r, s, t, x, n, stretch_rst );
            string cout =  "  ";
            string ts = r.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = s.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = t.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[0].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[1].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[2].ToString().PadLeft(8);
            cout += ts;
            Console.WriteLine(cout);

            r = 1.0;
            s = 0.0;
            t = 0.0;
            x = Blend.blend_rst_1dn ( r, s, t, x, n, stretch_rst );
            cout =  "  ";
            ts = r.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = s.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = t.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[0].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[1].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[2].ToString().PadLeft(8);
            cout += ts;
            Console.WriteLine(cout);

            r = 0.0;
            s = 1.0;
            t = 0.0;
            x = Blend.blend_rst_1dn ( r, s, t, x, n, stretch_rst );
            cout =  "  ";
            ts = r.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = s.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = t.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[0].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[1].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[2].ToString().PadLeft(8);
            cout += ts;
            Console.WriteLine(cout);

            r = 0.0;
            s = 0.0;
            t = 1.0;
            x = Blend.blend_rst_1dn ( r, s, t, x, n, stretch_rst );
            cout =  "  ";
            ts = r.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = s.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = t.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[0].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[1].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[2].ToString().PadLeft(8);
            cout += ts;
            Console.WriteLine(cout);

            r = 1.0;
            s = 1.0;
            t = 1.0;
            x = Blend.blend_rst_1dn ( r, s, t, x, n, stretch_rst );
            cout =  "  ";
            ts = r.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = s.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = t.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[0].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[1].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[2].ToString().PadLeft(8);
            cout += ts;
            Console.WriteLine(cout);

            r = 0.5;
            s = 0.5;
            t = 0.5;
            x = Blend.blend_rst_1dn ( r, s, t, x, n, stretch_rst );
            cout =  "  ";
            ts = r.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = s.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = t.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[0].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[1].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[2].ToString().PadLeft(8);
            cout += ts;
            Console.WriteLine(cout);
        }
        
        
        static void test12 ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    TEST12 tests BLEND_RST_2DN on the stretch map.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 October 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            double[] x = new double[3];

            Console.WriteLine();
            Console.WriteLine("TEST12");
            Console.WriteLine("  Stretch test on BLEND_RST_2DN.");

            int n = 3;

            double r = 0.0;
            double s = 0.0;
            double t = 0.0;
            x = Blend.blend_rst_2dn ( r, s, t, x, n, stretch_rst );
            string cout =  "  ";
            string ts = r.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = s.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = t.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[0].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[1].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[2].ToString().PadLeft(8);
            cout += ts;
            Console.WriteLine(cout);

            r = 1.0;
            s = 0.0;
            t = 0.0;
            x = Blend.blend_rst_2dn ( r, s, t, x, n, stretch_rst );
            cout =  "  ";
            ts = r.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = s.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = t.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[0].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[1].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[2].ToString().PadLeft(8);
            cout += ts;
            Console.WriteLine(cout);

            r = 0.0;
            s = 1.0;
            t = 0.0;
            x = Blend.blend_rst_2dn ( r, s, t, x, n, stretch_rst );
            cout =  "  ";
            ts = r.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = s.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = t.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[0].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[1].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[2].ToString().PadLeft(8);
            cout += ts;
            Console.WriteLine(cout);

            r = 0.0;
            s = 0.0;
            t = 1.0;
            x = Blend.blend_rst_2dn ( r, s, t, x, n, stretch_rst );
            cout =  "  ";
            ts = r.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = s.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = t.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[0].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[1].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[2].ToString().PadLeft(8);
            cout += ts;
            Console.WriteLine(cout);

            r = 1.0;
            s = 1.0;
            t = 1.0;
            x = Blend.blend_rst_2dn ( r, s, t, x, n, stretch_rst );
            cout =  "  ";
            ts = r.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = s.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = t.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[0].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[1].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[2].ToString().PadLeft(8);
            cout += ts;
            Console.WriteLine(cout);

            r = 0.5;
            s = 0.5;
            t = 0.5;
            x = Blend.blend_rst_2dn ( r, s, t, x, n, stretch_rst );
            cout =  "  ";
            ts = r.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = s.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = t.ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[0].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[1].ToString().PadLeft(8) + "  ";
            cout += ts;
            ts = x[2].ToString().PadLeft(8);
            cout += ts;
            Console.WriteLine(cout);
        }


        static void blend_i_0d1_test ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BLEND_I_0D1_TEST tests BLEND_I_0D1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 October 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            double[] x = new double[5];

            int m = 5;
            x[0] = 100.0;
            x[m-1] = 100.0 + ( m - 1 ) * 5;

            Console.WriteLine();
            Console.WriteLine("BLEND_I_0D1_TEST");
            Console.WriteLine("  BLEND_I_0D1 interpolates data in a vector.");
            Console.WriteLine();
            Console.WriteLine("  X[0] = " + x[0]);
            Console.WriteLine("  X(" + (m-1) + ")= " + x[m-1]);
            Console.WriteLine();
            Console.WriteLine("  Interpolated values:");
            Console.WriteLine();

            double[] x_ = Blend.blend_i_0d1 ( ref x, m );

            for (int i = 0; i < m; i++ )
            {
                string cout = "  ";
                string t = i.ToString().PadLeft(6) + "  ";
                cout += t;
                t = x_[i].ToString().PadLeft(8);
                cout += t;
                Console.WriteLine(cout);
            }
        }


        static void blend_ij_0d1_test ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BLEND_IJ_0D1_TEST tests BLEND_IJ_0D1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 October 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            int m1 = 5;
            int m2 = 4;
            double[] x = new double[20];

            Console.WriteLine();
            Console.WriteLine("BLEND_IJ_0D1_TEST");
            Console.WriteLine("  BLEND_IJ_0D1 interpolates data in a table,");
            Console.WriteLine("  from corner data.");
            Console.WriteLine();
            Console.WriteLine("  The table is " + m1 + " rows by " + m2 + " columns.");
            //
            //  Load data in the corners only.
            //
            int i = 0;
            int j = 0;
            double r = i / ( double ) ( m1 - 1 );
            double s = j / ( double ) ( m2 - 1 );
            double temp = cubic_rs ( r, s, 1 );
            x[i*m2+j] = temp;

            i = m1 - 1;
            j = 0;
            r = i / ( double ) ( m1 - 1 );
            s = j / ( double ) ( m2 - 1 );
            temp = cubic_rs ( r, s, 1 );
            x[i*m2+j] = temp;

            i = 0;
            j = m2 - 1;
            r = i / ( double ) ( m1 - 1 );
            s = j / ( double ) ( m2 - 1 );
            temp = cubic_rs ( r, s, 1 );
            x[i*m2+j] = temp;

            i = m1 - 1;
            j = m2 - 1;
            r = i / ( double ) ( m1 - 1 );
            s = j / ( double ) ( m2 - 1 );
            temp = cubic_rs ( r, s, 1 );
            x[i*m2+j] = temp;

            double[] x_ = Blend.blend_ij_0d1 ( x, m1, m2 );

            Console.WriteLine();
            Console.WriteLine("  Values interpolated by BLEND_IJ_0D1:");
            Console.WriteLine();

            for ( i = 0; i < m1; i++ )
            {
                string cout = "  ";
                string t = x_[i*m2].ToString().PadLeft(8) + "  ";
                cout += t;
                t = x_[i*m2+1].ToString().PadLeft(8) + "  ";
                cout += t;
                t = x_[i*m2+2].ToString().PadLeft(8) + "  ";
                cout += t;
                t = x_[i*m2+3].ToString().PadLeft(8);
                cout += t;
                Console.WriteLine(cout);
            }
        }

        
        static void blend_ij_1d1_test ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BLEND_IJ_1D1_TEST tests BLEND_IJ_1D1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 October 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            int m1 = 5;
            int m2 = 4;
            double[] x = new double[20];

            Console.WriteLine();
            Console.WriteLine("BLEND_IJ_1D1_TEST");
            Console.WriteLine("  BLEND_IJ_1D1 interpolates data in a 2D table,");
            Console.WriteLine("  from edge data.");
            Console.WriteLine();
            Console.WriteLine("  The table is " + m1 + " rows by " + m2 + " columns.");
            //
            //  Load data in the edges only.
            //
            for (int i = 0; i < m1; i++ )
            {
                double r = i / ( double ) ( m1 - 1 );

                int j = 0;
                double s = j / ( double ) ( m2 - 1 );
                double temp = cubic_rs ( r, s, 1 );
                x[i*m2+j] = temp;

                j = m2 - 1;
                s = j / ( double ) ( m2 - 1 );
                temp = cubic_rs ( r, s, 1 );
                x[i*m2+j] = temp;
            }

            for (int j = 0; j < m2; j++ )
            {
                double s = j / ( double ) ( m2 - 1 );

                int i = 0;
                double r = i / ( double ) ( m1 - 1 );
                double temp = cubic_rs ( r, s, 1 );
                x[i*m2+j] = temp;

                i = m1 - 1;
                r = i / ( double ) ( m1 - 1 );
                temp = cubic_rs ( r, s, 1 );
                x[i*m2+j] = temp;
            }

            double[] x_ = Blend.blend_ij_1d1 ( x, m1, m2 );

            Console.WriteLine();
            Console.WriteLine("  Values interpolated by BLEND_IJ_1D1:");
            Console.WriteLine();

            for (int i = 0; i < m1; i++ )
            {
                string cout = "  ";
                string t = x_[i*m2].ToString().PadLeft(8) + "  ";
                cout += t;
                t = x_[i*m2+1].ToString().PadLeft(8) + "  ";
                cout += t;
                t = x_[i*m2+2].ToString().PadLeft(8) + "  ";
                cout += t;
                t = x_[i*m2+3].ToString().PadLeft(8);
                cout += t;
                Console.WriteLine(cout);
            }
        }
        
        
        static void blend_ijk_0d1_test ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BLEND_IJK_0D1_TEST tests BLEND_IJK_0D1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 October 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            int m1 = 4;
            int m2 = 3;
            int m3 = 3;
            double temp;
            double[] x = new double[36];

            Console.WriteLine();
            Console.WriteLine("BLEND_IJK_0D1_TEST");
            Console.WriteLine("  BLEND_IJK_0D1_TEST interpolates data in a 3D table,");
            Console.WriteLine("  from corner data.");
            Console.WriteLine();
            Console.WriteLine("  The table is " + m1 + " rows by " + m2 + " columns by " + m3 + " layers.");

            //
            //  Load data in the faces only.
            //
            for (int i = 0; i < m1; i++ )
            {
                double r = i / ( double ) ( m1 - 1 );
                for (int j = 0; j < m2; j++ )
                {
                    double s = j / ( double ) ( m2 - 1 );
                    for (int k = 0; k < m3; k++ )
                    {
                        double t = k / ( double ) ( m3 - 1 );
                        int num_extreme = 0;
                        if ( i == 0 || i == m1 - 1 )
                        {
                            num_extreme = num_extreme + 1;
                        }
                        if ( j == 0 || j == m2 - 1 )
                        {
                            num_extreme = num_extreme + 1;
                        }
                        if ( k == 0 || k == m3 - 1 )
                        {
                            num_extreme = num_extreme + 1;
                        }
                        if ( num_extreme >= 3 )
                        {
                            temp = quad_rst ( r, s, t, 1 );
                        }
                        else
                        {
                            temp = 0.0;
                        }
                        x[(i*m3+j)*m2+k] = temp;
                    }
                }
            }

            Console.WriteLine();
            Console.WriteLine("  Data given to BLEND_IJK_0D1:");
            Console.WriteLine();
            
            for (int k = 0; k < m3; k++ )
            {
                Console.WriteLine();
                Console.WriteLine("  Layer K = " + k);
                Console.WriteLine();

                for (int i = 0; i < m1; i++)
                {
                    string cout = "  ";
                    string t = x[(i * m3 + 0) * m2 + k].ToString().PadLeft(8) + "  ";
                    cout += t;
                    t = x[(i * m3 + 1) * m2 + k].ToString().PadLeft(8) + "  ";
                    cout += t;
                    t = x[(i * m3 + 2) * m2 + k].ToString().PadLeft(8);
                    cout += t;
                    Console.WriteLine(cout);
                }
            }
            x = Blend.blend_ijk_0d1 ( x, m1, m2, m3 );

            Console.WriteLine();
            Console.WriteLine("  Values interpolated by BLEND_IJK_0D1:");
            Console.WriteLine();

            for (int k = 0; k < m3; k++ )
            {
                Console.WriteLine();
                Console.WriteLine("  Layer K = " + k);
                Console.WriteLine();

                for (int i = 0; i < m1; i++ )
                {
                    string cout = "  ";
                    string t = x[(i * m3 + 0) * m2 + k].ToString().PadLeft(8) + "  ";
                    cout += t;
                    t = x[(i * m3 + 1) * m2 + k].ToString().PadLeft(8) + "  ";
                    cout += t;
                    t = x[(i * m3 + 2) * m2 + k].ToString().PadLeft(8);
                    cout += t;
                    Console.WriteLine(cout);
                }
            }

            for (int i = 0; i < m1; i++ )
            {
                double r = i / ( double ) ( m1 - 1 );
                for (int j = 0; j < m2; j++ )
                {
                    double s = j / ( double ) ( m2 - 1 );
                    for (int k = 0; k < m3; k++ )
                    {
                        double t = k / ( double ) ( m3 - 1 );
                        temp = quad_rst ( r, s, t, 1 );
                        x[(i*m3+j)*m2+k] = temp;
                    }
                }
            }

            Console.WriteLine();
            Console.WriteLine("  Exact data:");
            Console.WriteLine();

            for (int k = 0; k < m3; k++ )
            {
                Console.WriteLine();
                Console.WriteLine("  Layer K = " + k);
                Console.WriteLine();

                for (int i = 0; i < m1; i++ )
                {
                    string cout = "  ";
                    string t = x[(i * m3 + 0) * m2 + k].ToString().PadLeft(8) + "  ";
                    cout += t;
                    t = x[(i * m3 + 1) * m2 + k].ToString().PadLeft(8) + "  ";
                    cout += t;
                    t = x[(i * m3 + 2) * m2 + k].ToString().PadLeft(8);
                    cout += t;
                    Console.WriteLine(cout);
                }
            }
        }


        static void blend_ijk_1d1_test ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BLEND_IJK_1D1_TEST tests BLEND_IJK_1D1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 October 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            int m1 = 4;
            int m2 = 3;
            int m3 = 3;
            double temp;
            double[] x = new double[36];

            Console.WriteLine();
            Console.WriteLine("BLEND_IJK_1D1_TEST");
            Console.WriteLine("  BLEND_IJK_1D1_TEST interpolates data in a 3D table,");
            Console.WriteLine("  from edge data.");
            Console.WriteLine();
            Console.WriteLine("  The table is " + m1 + " rows by " + m2 + " columns by " + m3 + " layers.");
            //
            //  Load data in the faces only.
            //
            for (int i = 0; i < m1; i++ )
            {
                double r = i / ( double ) ( m1 - 1 );
                for (int j = 0; j < m2; j++ )
                {
                    double s = j / ( double ) ( m2 - 1 );
                    for (int k = 0; k < m3; k++ )
                    {
                        double t = k / ( double ) ( m3 - 1 );
                        int num_extreme = 0;
                        if ( i == 0 || i == m1 - 1 )
                        {
                            num_extreme = num_extreme + 1;
                        }
                        if ( j == 0 || j == m2 - 1 )
                        {
                            num_extreme = num_extreme + 1;
                        }
                        if ( k == 0 || k == m3 - 1 )
                        {
                            num_extreme = num_extreme + 1;
                        }
                        if ( num_extreme >= 2 )
                        {
                            temp = quad_rst ( r, s, t, 1 );
                        }
                        else
                        {
                            temp = 0.0;
                        }
                        x[(i*m3+j)*m2+k] = temp;
                    }
                }
            }

            Console.WriteLine();
            Console.WriteLine("  Data given to BLEND_IJK_1D1:");
            Console.WriteLine();

            for (int k = 0; k < m3; k++ )
            {
                Console.WriteLine();
                Console.WriteLine("  Layer K = " + k);
                Console.WriteLine();

                for (int i = 0; i < m1; i++ )
                {
                    string cout = "  ";
                    string t = x[(i * m3 + 0) * m2 + k].ToString().PadLeft(8) + "  ";
                    cout += t;
                    t = x[(i * m3 + 1) * m2 + k].ToString().PadLeft(8) + "  ";
                    cout += t;
                    t = x[(i * m3 + 2) * m2 + k].ToString().PadLeft(8);
                    cout += t;
                    Console.WriteLine(cout);
                }
            }
            x = Blend.blend_ijk_1d1 ( x, m1, m2, m3 );

            Console.WriteLine();
            Console.WriteLine("  Values interpolated by BLEND_IJK_1D1:");
            Console.WriteLine();

            for (int k = 0; k < m3; k++ )
            {
                Console.WriteLine();
                Console.WriteLine("  Layer K = " + k);
                Console.WriteLine();

                for (int i = 0; i < m1; i++ )
                {
                    string cout = "  ";
                    string t = x[(i * m3 + 0) * m2 + k].ToString().PadLeft(8) + "  ";
                    cout += t;
                    t = x[(i * m3 + 1) * m2 + k].ToString().PadLeft(8) + "  ";
                    cout += t;
                    t = x[(i * m3 + 2) * m2 + k].ToString().PadLeft(8);
                    cout += t;
                    Console.WriteLine(cout);
                }
            }

            for (int i = 0; i < m1; i++ )
            {
                double r = i / ( double ) ( m1 - 1 );
                for (int j = 0; j < m2; j++ )
                {
                    double s = j / ( double ) ( m2 - 1 );
                    for (int k = 0; k < m3; k++ )
                    {
                        double t = k / ( double ) ( m3 - 1 );
                        temp = quad_rst ( r, s, t, 1 );
                        x[(i*m3+j)*m2+k] = temp;
                    }
                }
            }

            Console.WriteLine();
            Console.WriteLine("  Exact data:");
            Console.WriteLine();

            for (int k = 0; k < m3; k++ )
            {
                Console.WriteLine();
                Console.WriteLine("  Layer K = " + k);
                Console.WriteLine();

                for (int i = 0; i < m1; i++ )
                {
                    string cout = "  ";
                    string t = x[(i * m3 + 0) * m2 + k].ToString().PadLeft(8) + "  ";
                    cout += t;
                    t = x[(i * m3 + 1) * m2 + k].ToString().PadLeft(8) + "  ";
                    cout += t;
                    t = x[(i * m3 + 2) * m2 + k].ToString().PadLeft(8);
                    cout += t;
                    Console.WriteLine(cout);
                }
            }
        }

        
        static void blend_ijk_2d1_test ( )
        //****************************************************************************80
        //
        //  Purpose:
        //
        //    BLEND_IJK_2D1_TEST tests BLEND_IJK_2D1.
        //
        //  Licensing:
        //
        //    This code is distributed under the GNU LGPL license.
        //
        //  Modified:
        //
        //    22 October 2008
        //
        //  Author:
        //
        //    John Burkardt
        //
        {
            int m1 = 4;
            int m2 = 3;
            int m3 = 3;
            double temp;
            double[] x = new double[36];

            Console.WriteLine();
            Console.WriteLine("BLEND_IJK_2D1_TEST");
            Console.WriteLine("  BLEND_IJK_2D1_TEST interpolates data in a 3D table,");
            Console.WriteLine("  from face data.");
            Console.WriteLine();
            Console.WriteLine("  The table is " + m1 + " rows by " + m2 + " columns by " + m3 + " layers.");

            //
            //  Load data in the faces only.
            //
            for ( int i = 0; i < m1; i++ )
            {
                double r = i / ( double ) ( m1 - 1 );
                for ( int j = 0; j < m2; j++ )
                {
                    double s = j / ( double ) ( m2 - 1 );
                    for ( int k = 0; k < m3; k++ )
                    {
                        double t = k / ( double ) ( m3 - 1 );
                        if ( i == 0 || i == m1 - 1 ||
                        j == 0 || j == m2 - 1 ||
                        k == 0 || k == m3 - 1 )
                        {
                            temp = quad_rst ( r, s, t, 1 );
                        }
                        else
                        {
                            temp = 0.0;
                        }
                        x[(i*m3+j)*m2+k] = temp;
                    }
                }
            }

            Console.WriteLine();
            Console.WriteLine("  Data given to BLEND_IJK_2D1:");
            Console.WriteLine();

            for ( int k = 0; k < m3; k++ )
            {
                Console.WriteLine();
                Console.WriteLine("  Layer K = " + k);
                Console.WriteLine();

                for ( int i = 0; i < m1; i++ )
                {
                    string cout = "  ";
                    string t = x[(i * m3 + 0) * m2 + k].ToString().PadLeft(8) + "  ";
                    cout += t;
                    t = x[(i * m3 + 1) * m2 + k].ToString().PadLeft(8) + "  ";
                    cout += t;
                    t = x[(i * m3 + 2) * m2 + k].ToString().PadLeft(8);
                    cout += t;
                    Console.WriteLine(cout);
                }
            }
            x = Blend.blend_ijk_2d1 ( x, m1, m2, m3 );

            Console.WriteLine();
            Console.WriteLine("  Values interpolated by BLEND_IJK_2D1:");
            Console.WriteLine();

            for ( int k = 0; k < m3; k++ )
            {
                Console.WriteLine();
                Console.WriteLine("  Layer K = " + k);
                Console.WriteLine();

                for ( int i = 0; i < m1; i++ )
                {
                    string cout = "  ";
                    string t = x[(i * m3 + 0) * m2 + k].ToString().PadLeft(8) + "  ";
                    cout += t;
                    t = x[(i * m3 + 1) * m2 + k].ToString().PadLeft(8) + "  ";
                    cout += t;
                    t = x[(i * m3 + 2) * m2 + k].ToString().PadLeft(8);
                    cout += t;
                    Console.WriteLine(cout);
                }
            }

            for ( int i = 0; i < m1; i++ )
            {
                double r = i / ( double ) ( m1 - 1 );
                for ( int j = 0; j < m2; j++ )
                {
                    double s = j / ( double ) ( m2 - 1 );
                    for ( int k = 0; k < m3; k++ )
                    {
                        double t = k / ( double ) ( m3 - 1 );
                        temp = quad_rst ( r, s, t, 1 );
                        x[(i*m3+j)*m2+k] = temp;
                    }
                }
            }

            Console.WriteLine();
            Console.WriteLine("  Exact data:");
            Console.WriteLine();

            for ( int k = 0; k < m3; k++ )
            {
                Console.WriteLine();
                Console.WriteLine("  Layer K = " + k);
                Console.WriteLine();

                for ( int i = 0; i < m1; i++ )
                {
                    string cout = "  ";
                    string t = x[(i * m3 + 0) * m2 + k].ToString().PadLeft(8) + "  ";
                    cout += t;
                    t = x[(i * m3 + 1) * m2 + k].ToString().PadLeft(8) + "  ";
                    cout += t;
                    t = x[(i * m3 + 2) * m2 + k].ToString().PadLeft(8);
                    cout += t;
                    Console.WriteLine(cout);
                }
            }
        }

        
        
    }
}