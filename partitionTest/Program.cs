using System;
using System.Collections.Generic;
using System.Linq;
using System.Net.Security;
using ClipperLib;
using geoLib;
using geoWrangler;
using geoCoreLib;

namespace partitionTest
{
    using Path = List<IntPoint>;
    using Paths = List<List<IntPoint>>;
    class Program
    {
        /*
        // Unusued
        static void ZFillCallback(IntPoint bot1, IntPoint top1, IntPoint bot2, IntPoint top2, ref IntPoint pt)
        {
            pt.Z = -1; // Tag our intersection points.
        }
        */

        static int scaling = 10000;

        static void Main(string[] args)
        {
            // L
            GeoLibPoint[] L = new GeoLibPoint[] {

            new GeoLibPoint( 0, 0),
            new GeoLibPoint( 0, 50),
            new GeoLibPoint( 10, 50),
            new GeoLibPoint( 10, 20),
            new GeoLibPoint( 60, 20),
            new GeoLibPoint( 60, 0),
            new GeoLibPoint( 0, 0)

            };

            GeoLibPoint[] rL = new GeoLibPoint[] {

            new GeoLibPoint( 0, 0),
            new GeoLibPoint( 0, 20),
            new GeoLibPoint( 10, 20),
            new GeoLibPoint( 10, 50),
            new GeoLibPoint( 60, 50),
            new GeoLibPoint( 60, 0),
            new GeoLibPoint( 0, 0)

            };

            // Reversed orientation.
            GeoLibPoint[] L_ccw = L.ToArray();
            Array.Reverse(L_ccw);


            // U
            GeoLibPoint[] U = new GeoLibPoint[] {

            new GeoLibPoint( 0, 0),
            new GeoLibPoint( 0, 50),
            new GeoLibPoint( 10, 50),
            new GeoLibPoint( 10, 20),
            new GeoLibPoint( 60, 20),
            new GeoLibPoint( 60, 80),
            new GeoLibPoint( 120, 80),
            new GeoLibPoint( 120, 0),
            new GeoLibPoint( 0, 0)

            };

            // T
            GeoLibPoint[] T = new GeoLibPoint[] {

            new GeoLibPoint( 0, 50),
            new GeoLibPoint( 0, 80),
            new GeoLibPoint( 80, 80),
            new GeoLibPoint( 80, 50),
            new GeoLibPoint( 60, 50),
            new GeoLibPoint( 60, 0),
            new GeoLibPoint( 40, 0),
            new GeoLibPoint( 40, 50),
            new GeoLibPoint( 0, 50)

            };


            // X
            GeoLibPoint[] X = new GeoLibPoint[] {

            new GeoLibPoint( 0, 50),
            new GeoLibPoint( 0, 80),
            new GeoLibPoint( 60, 80),
            new GeoLibPoint( 60, 100),
            new GeoLibPoint( 80, 100),
            new GeoLibPoint( 80, 80),
            new GeoLibPoint( 100, 80),
            new GeoLibPoint( 100, 50),
            new GeoLibPoint( 80, 50),
            new GeoLibPoint( 80, 20),
            new GeoLibPoint( 60, 20),
            new GeoLibPoint( 60, 50),
            new GeoLibPoint( 0, 50)

            };


            // S
            GeoLibPoint[] S = new GeoLibPoint[] {

            new GeoLibPoint( 0, 0),
            new GeoLibPoint( 0, 20),
            new GeoLibPoint( 20, 20),
            new GeoLibPoint( 20, 50),
            new GeoLibPoint( 0, 50),
            new GeoLibPoint( 0, 110),
            new GeoLibPoint( 100, 110),
            new GeoLibPoint( 100, 80),
            new GeoLibPoint( 80, 80),
            new GeoLibPoint( 80, 60),
            new GeoLibPoint( 100, 60),
            new GeoLibPoint( 100, 0),
            new GeoLibPoint( 0, 0)

            };

            // Negative S
            GeoLibPoint[] nS = new GeoLibPoint[] {

            new GeoLibPoint( 0, 0),
            new GeoLibPoint( 0, 20),
            new GeoLibPoint( 20, 20),
            new GeoLibPoint( 20, 50),
            new GeoLibPoint( 0, 50),
            new GeoLibPoint( 0, 110),
            new GeoLibPoint( 100, 110),
            new GeoLibPoint( 100, 80),
            new GeoLibPoint( 80, 80),
            new GeoLibPoint( 80, 60),
            new GeoLibPoint( 100, 60),
            new GeoLibPoint( 100, 0),
            new GeoLibPoint( 0, 0)

            };

            for (int i = 0; i < nS.Length; i++)
            {
                nS[i] = new GeoLibPoint(nS[i].X, nS[i].Y - 200);
            }

            // Complex 1
            GeoLibPoint[] C1 = new GeoLibPoint[] {
            new GeoLibPoint(10,  -61),
            new GeoLibPoint(10,  -50),
            new GeoLibPoint(0,   -50),
            new GeoLibPoint(0,   -60),
            new GeoLibPoint(-10, -60),
            new GeoLibPoint(-10, -30),
            new GeoLibPoint(0,   -30),
            new GeoLibPoint(0,   -40),
            new GeoLibPoint(10,  -40),
            new GeoLibPoint(10,  0),
            new GeoLibPoint(0,   0),
            new GeoLibPoint(0,   10),
            new GeoLibPoint(-60, 10),
            new GeoLibPoint(-60, 20),
            new GeoLibPoint(-50, 20),
            new GeoLibPoint(-50, 30),
            new GeoLibPoint(-60, 30),
            new GeoLibPoint(-60, 40),
            new GeoLibPoint(-30, 40),
            new GeoLibPoint(-30, 30),
            new GeoLibPoint(-40, 30),
            new GeoLibPoint(-40, 20),
            new GeoLibPoint(0,   20),
            new GeoLibPoint(0,   30),
            new GeoLibPoint(10,  30),
            new GeoLibPoint(10,  91),
            new GeoLibPoint(20,  91),
            new GeoLibPoint(20,  80),
            new GeoLibPoint(30,  80),
            new GeoLibPoint(30,  90),
            new GeoLibPoint(40,  90),
            new GeoLibPoint(40,  60),
            new GeoLibPoint(30,  60),
            new GeoLibPoint(30,  70),
            new GeoLibPoint(20,  70),
            new GeoLibPoint(20,  30),
            new GeoLibPoint(30,  30),
            new GeoLibPoint(30,  20),
            new GeoLibPoint(90,  20),
            new GeoLibPoint(90,  10),
            new GeoLibPoint(80,  10),
            new GeoLibPoint(80,  0),
            new GeoLibPoint(90,  0),
            new GeoLibPoint(90,  -10),
            new GeoLibPoint(60,  -10),
            new GeoLibPoint(60,  0),
            new GeoLibPoint(70,  0),
            new GeoLibPoint(70,  10),
            new GeoLibPoint(30,  10),
            new GeoLibPoint(30,  0),
            new GeoLibPoint(20,  0),
            new GeoLibPoint(20,  -61),
            new GeoLibPoint(10,  -61)
            };
            C1 = GeoWrangler.clockwiseAndReorder(C1);


            // Complex 2
            GeoLibPoint[] C2 = new GeoLibPoint[]
            {
            new GeoLibPoint(10, -213),
            new GeoLibPoint(10, -202),
            new GeoLibPoint(0,  -202),
            new GeoLibPoint(0,  -212),
            new GeoLibPoint(-10,-212),
            new GeoLibPoint(-10,-182),
            new GeoLibPoint(0,  -182),
            new GeoLibPoint(0,  -192),
            new GeoLibPoint(10, -192),
            new GeoLibPoint(10, -152),
            new GeoLibPoint(0,  -152),
            new GeoLibPoint(0,  -142),
            new GeoLibPoint(-60,-142),
            new GeoLibPoint(-60,-132),
            new GeoLibPoint(-50,-132),
            new GeoLibPoint(-50,-122),
            new GeoLibPoint(-60,-122),
            new GeoLibPoint(-60,-112),
            new GeoLibPoint(-30,-112),
            new GeoLibPoint(-30,-122),
            new GeoLibPoint(-40,-122),
            new GeoLibPoint(-40,-132),
            new GeoLibPoint(0,  -132),
            new GeoLibPoint(0,  -122),
            new GeoLibPoint(10, -122),
            new GeoLibPoint(10, -50),
            new GeoLibPoint(0,  -50),
            new GeoLibPoint(0,  -60),
            new GeoLibPoint(-10,-60),
            new GeoLibPoint(-10,-30),
            new GeoLibPoint(0,  -30),
            new GeoLibPoint(0,  -40),
            new GeoLibPoint(10, -40),
            new GeoLibPoint(10, 0),
            new GeoLibPoint(0,  0),
            new GeoLibPoint(0,  10),
            new GeoLibPoint(-70,10),
            new GeoLibPoint(-70,0),
            new GeoLibPoint(-60,0),
            new GeoLibPoint(-60,-10),
            new GeoLibPoint(-90,-10),
            new GeoLibPoint(-90,0),
            new GeoLibPoint(-80,0),
            new GeoLibPoint(-80,10),
            new GeoLibPoint(-120,  10),
            new GeoLibPoint(-120,  0),
            new GeoLibPoint(-130,  0),
            new GeoLibPoint(-130,  -61),
            new GeoLibPoint(-140,  -61),
            new GeoLibPoint(-140,  -50),
            new GeoLibPoint(-150,  -50),
            new GeoLibPoint(-150,  -60),
            new GeoLibPoint(-160,  -60),
            new GeoLibPoint(-160,  -30),
            new GeoLibPoint(-150,  -30),
            new GeoLibPoint(-150,  -40),
            new GeoLibPoint(-140,  -40),
            new GeoLibPoint(-140,  0),
            new GeoLibPoint(-150,  0),
            new GeoLibPoint(-150,  10),
            new GeoLibPoint(-210,  10),
            new GeoLibPoint(-210,  20),
            new GeoLibPoint(-200,  20),
            new GeoLibPoint(-200,  30),
            new GeoLibPoint(-210,  30),
            new GeoLibPoint(-210,  40),
            new GeoLibPoint(-180,  40),
            new GeoLibPoint(-180,  30),
            new GeoLibPoint(-190,  30),
            new GeoLibPoint(-190,  20),
            new GeoLibPoint(-150,  20),
            new GeoLibPoint(-150,  30),
            new GeoLibPoint(-140,  30),
            new GeoLibPoint(-140,  91),
            new GeoLibPoint(-130,  91),
            new GeoLibPoint(-130,  80),
            new GeoLibPoint(-120,  80),
            new GeoLibPoint(-120,  90),
            new GeoLibPoint(-110,  90),
            new GeoLibPoint(-110,  60),
            new GeoLibPoint(-120,  60),
            new GeoLibPoint(-120,  70),
            new GeoLibPoint(-130,  70),
            new GeoLibPoint(-130,  30),
            new GeoLibPoint(-120,  30),
            new GeoLibPoint(-120,  20),
            new GeoLibPoint(-50,20),
            new GeoLibPoint(-50,30),
            new GeoLibPoint(-60,30),
            new GeoLibPoint(-60,40),
            new GeoLibPoint(-30,40),
            new GeoLibPoint(-30,30),
            new GeoLibPoint(-40,30),
            new GeoLibPoint(-40,20),
            new GeoLibPoint(0,  20),
            new GeoLibPoint(0,  30),
            new GeoLibPoint(10, 30),
            new GeoLibPoint(10, 102),
            new GeoLibPoint(0,  102),
            new GeoLibPoint(0,  92),
            new GeoLibPoint(-10,92),
            new GeoLibPoint(-10,122),
            new GeoLibPoint(0,  122),
            new GeoLibPoint(0,  112),
            new GeoLibPoint(10, 112),
            new GeoLibPoint(10, 152),
            new GeoLibPoint(0,  152),
            new GeoLibPoint(0,  162),
            new GeoLibPoint(-60,162),
            new GeoLibPoint(-60,172),
            new GeoLibPoint(-50,172),
            new GeoLibPoint(-50,182),
            new GeoLibPoint(-60,182),
            new GeoLibPoint(-60,192),
            new GeoLibPoint(-30,192),
            new GeoLibPoint(-30,182),
            new GeoLibPoint(-40,182),
            new GeoLibPoint(-40,172),
            new GeoLibPoint(0,  172),
            new GeoLibPoint(0,  182),
            new GeoLibPoint(10, 182),
            new GeoLibPoint(10, 243),
            new GeoLibPoint(20, 243),
            new GeoLibPoint(20, 232),
            new GeoLibPoint(30, 232),
            new GeoLibPoint(30, 242),
            new GeoLibPoint(40, 242),
            new GeoLibPoint(40, 212),
            new GeoLibPoint(30, 212),
            new GeoLibPoint(30, 222),
            new GeoLibPoint(20, 222),
            new GeoLibPoint(20, 182),
            new GeoLibPoint(30, 182),
            new GeoLibPoint(30, 172),
            new GeoLibPoint(90, 172),
            new GeoLibPoint(90, 162),
            new GeoLibPoint(80, 162),
            new GeoLibPoint(80, 152),
            new GeoLibPoint(90, 152),
            new GeoLibPoint(90, 142),
            new GeoLibPoint(60, 142),
            new GeoLibPoint(60, 152),
            new GeoLibPoint(70, 152),
            new GeoLibPoint(70, 162),
            new GeoLibPoint(30, 162),
            new GeoLibPoint(30, 152),
            new GeoLibPoint(20, 152),
            new GeoLibPoint(20, 80),
            new GeoLibPoint(30, 80),
            new GeoLibPoint(30, 90),
            new GeoLibPoint(40, 90),
            new GeoLibPoint(40, 60),
            new GeoLibPoint(30, 60),
            new GeoLibPoint(30, 70),
            new GeoLibPoint(20, 70),
            new GeoLibPoint(20, 30),
            new GeoLibPoint(30, 30),
            new GeoLibPoint(30, 20),
            new GeoLibPoint(100,20),
            new GeoLibPoint(100,30),
            new GeoLibPoint(90, 30),
            new GeoLibPoint(90, 40),
            new GeoLibPoint(120,40),
            new GeoLibPoint(120,30),
            new GeoLibPoint(110,30),
            new GeoLibPoint(110,20),
            new GeoLibPoint(150,20),
            new GeoLibPoint(150,30),
            new GeoLibPoint(160,30),
            new GeoLibPoint(160,91),
            new GeoLibPoint(170,91),
            new GeoLibPoint(170,80),
            new GeoLibPoint(180,80),
            new GeoLibPoint(180,90),
            new GeoLibPoint(190,90),
            new GeoLibPoint(190,60),
            new GeoLibPoint(180,60),
            new GeoLibPoint(180,70),
            new GeoLibPoint(170,70),
            new GeoLibPoint(170,30),
            new GeoLibPoint(180,30),
            new GeoLibPoint(180,20),
            new GeoLibPoint(240,20),
            new GeoLibPoint(240,10),
            new GeoLibPoint(230,10),
            new GeoLibPoint(230,0),
            new GeoLibPoint(240,0),
            new GeoLibPoint(240,-10),
            new GeoLibPoint(210,-10),
            new GeoLibPoint(210,0),
            new GeoLibPoint(220,0),
            new GeoLibPoint(220,10),
            new GeoLibPoint(180,10),
            new GeoLibPoint(180,0),
            new GeoLibPoint(170,0),
            new GeoLibPoint(170,-61),
            new GeoLibPoint(160,-61),
            new GeoLibPoint(160,-50),
            new GeoLibPoint(150,-50),
            new GeoLibPoint(150,-60),
            new GeoLibPoint(140,-60),
            new GeoLibPoint(140,-30),
            new GeoLibPoint(150,-30),
            new GeoLibPoint(150,-40),
            new GeoLibPoint(160,-40),
            new GeoLibPoint(160,0),
            new GeoLibPoint(150,0),
            new GeoLibPoint(150,10),
            new GeoLibPoint(80, 10),
            new GeoLibPoint(80, 0),
            new GeoLibPoint(90, 0),
            new GeoLibPoint(90, -10),
            new GeoLibPoint(60, -10),
            new GeoLibPoint(60, 0),
            new GeoLibPoint(70, 0),
            new GeoLibPoint(70, 10),
            new GeoLibPoint(30, 10),
            new GeoLibPoint(30, 0),
            new GeoLibPoint(20, 0),
            new GeoLibPoint(20, -72),
            new GeoLibPoint(30, -72),
            new GeoLibPoint(30, -62),
            new GeoLibPoint(40, -62),
            new GeoLibPoint(40, -92),
            new GeoLibPoint(30, -92),
            new GeoLibPoint(30, -82),
            new GeoLibPoint(20, -82),
            new GeoLibPoint(20, -122),
            new GeoLibPoint(30, -122),
            new GeoLibPoint(30, -132),
            new GeoLibPoint(90, -132),
            new GeoLibPoint(90, -142),
            new GeoLibPoint(80, -142),
            new GeoLibPoint(80, -152),
            new GeoLibPoint(90, -152),
            new GeoLibPoint(90, -162),
            new GeoLibPoint(60, -162),
            new GeoLibPoint(60, -152),
            new GeoLibPoint(70, -152),
            new GeoLibPoint(70, -142),
            new GeoLibPoint(30, -142),
            new GeoLibPoint(30, -152),
            new GeoLibPoint(20, -152),
            new GeoLibPoint(20, -213),
            new GeoLibPoint(10, -213),
            };

            C2 = GeoWrangler.clockwiseAndReorder(C2);

            // Complex 3
            GeoLibPoint[] C3 = new GeoLibPoint[]
            {
            new GeoLibPoint(-63, -63),
            new GeoLibPoint(-63, 63),
            new GeoLibPoint(2, 63),
            new GeoLibPoint(2, 57),
            new GeoLibPoint(-57, 57),
            new GeoLibPoint(-57, -57),
            new GeoLibPoint(57, -57),
            new GeoLibPoint(57, 37),
            new GeoLibPoint(-37, 37),
            new GeoLibPoint(-37, -37),
            new GeoLibPoint(37, -37),
            new GeoLibPoint(37, 17),
            new GeoLibPoint(-17, 17),
            new GeoLibPoint(-17, -17),
            new GeoLibPoint(17, -17),
            new GeoLibPoint(17, -3),
            new GeoLibPoint(-2, -3),
            new GeoLibPoint(-2, 3),
            new GeoLibPoint(23, 3),
            new GeoLibPoint(23, -23),
            new GeoLibPoint(-23, -23),
            new GeoLibPoint(-23, 23),
            new GeoLibPoint(43, 23),
            new GeoLibPoint(43, -43),
            new GeoLibPoint(-43, -43),
            new GeoLibPoint(-43, 43),
            new GeoLibPoint(63, 43),
            new GeoLibPoint(63, -63),
            new GeoLibPoint(-63, -63),
            };

            // C3 = GeoWrangler.clockwiseAndReorder(C3);

            bool orth = GeoWrangler.orthogonal(C2);
            bool orth2 = GeoWrangler.orthogonal(C3);

            // Staircase
            GeoLibPoint[] S1 = new GeoLibPoint[] {

            new GeoLibPoint( -50, -50),
            new GeoLibPoint( 200,  -50),
            new GeoLibPoint( 200,   300),
            new GeoLibPoint( 150,   300),
            new GeoLibPoint( 150,   200),
            new GeoLibPoint( 100,   200),
            new GeoLibPoint( 100,   120),
            new GeoLibPoint( 0,   120),
            new GeoLibPoint( 0,   0),
            new GeoLibPoint( -50,   0),
            new GeoLibPoint( -50,   -50),
        };

            S1 = GeoWrangler.clockwiseAndReorder(S1);

            // Rectangular decomposition will return non-orthogonal polygons in the output when encountered.

            List<GeoLibPoint[]> l = GeoWrangler.rectangular_decomposition(L, maxRayLength: 1000);

            writeToLayout("l", L, l);

            List<GeoLibPoint[]> lccw = GeoWrangler.rectangular_decomposition(L_ccw);

            writeToLayout("lccw", L_ccw, lccw);

            List<GeoLibPoint[]> rl = GeoWrangler.rectangular_decomposition(rL);

            writeToLayout("rl", rL, rl);

            List<GeoLibPoint[]> u = GeoWrangler.rectangular_decomposition(U);

            writeToLayout("u", U, u);

            List<GeoLibPoint[]> t = GeoWrangler.rectangular_decomposition(T);

            writeToLayout("t", T, t);

            List<GeoLibPoint[]> x = GeoWrangler.rectangular_decomposition(X);

            writeToLayout("x", X, x);

            List<GeoLibPoint[]> s = GeoWrangler.rectangular_decomposition(S);

            writeToLayout("s", S, s);

            List<GeoLibPoint[]> ns = GeoWrangler.rectangular_decomposition(nS);

            writeToLayout("ns", nS, ns);

            List<GeoLibPoint[]> c1 = GeoWrangler.rectangular_decomposition(C1);

            writeToLayout("c1", C1, c1);

            List<GeoLibPoint[]> c2 = GeoWrangler.rectangular_decomposition(C2);

            writeToLayout("c2", C2, c2);

            List<GeoLibPoint[]> c3 = GeoWrangler.rectangular_decomposition(C3);

            writeToLayout("c3", C3, c3);

            List<GeoLibPoint[]> s1 = GeoWrangler.rectangular_decomposition(S1);

            writeToLayout("s1", S1, s1);

            Console.WriteLine("Hello World!");
        }


        static void writeToLayout(string filename, GeoLibPoint[] orig, List<GeoLibPoint[]> decomped)
        {
            // Can the system define geometry and write it correctly to Oasis and GDS files.
            GeoCore g = new GeoCore();
            g.reset();
            GCDrawingfield drawing_ = new GCDrawingfield("");
            drawing_.accyear = 2018;
            drawing_.accmonth = 12;
            drawing_.accday = 5;
            drawing_.acchour = 2;
            drawing_.accmin = 10;
            drawing_.accsec = 10;
            drawing_.modyear = 2018;
            drawing_.modmonth = 12;
            drawing_.modday = 5;
            drawing_.modhour = 2;
            drawing_.modmin = 10;
            drawing_.modsec = 10;
            drawing_.databaseunits = 1000;
            drawing_.userunits = 1E-3;// 0.001 / 1E-6;
            drawing_.libname = "noname";

            GCCell gcell = drawing_.addCell();
            gcell.accyear = 2018;
            gcell.accmonth = 12;
            gcell.accday = 5;
            gcell.acchour = 2;
            gcell.accmin = 10;
            gcell.accsec = 10;
            gcell.modyear = 2018;
            gcell.modmonth = 12;
            gcell.modday = 5;
            gcell.modhour = 2;
            gcell.modmin = 10;
            gcell.modsec = 10;

            gcell.cellName = "test";

            gcell.addPolygon(orig, 1, 0);

            for (int i = 0; i < decomped.Count; i++)
            {
                gcell.addBox(decomped[i], 1, 1);
            }

            g.setDrawing(drawing_);
            g.setValid(true);

            gds.gdsWriter gw = new gds.gdsWriter(g, "../../../../decomp_out/" + filename + "_partitiontest.gds");
            gw.save();

            oasis.oasWriter ow = new oasis.oasWriter(g, "../../../../decomp_out/" + filename + "_partitiontest.oas");
            ow.save();
        }
    }
}
