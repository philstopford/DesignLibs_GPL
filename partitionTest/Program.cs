using System;
using System.Collections.Generic;
using System.Linq;
using System.Net.Security;
using ClipperLib;
using geoLib;
using geoWrangler;

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

            List<GeoLibPoint[]> l = GeoWrangler.rectangular_decomposition(L);

            List<GeoLibPoint[]> rl = GeoWrangler.rectangular_decomposition(rL);

            List<GeoLibPoint[]> u = GeoWrangler.rectangular_decomposition(U);

            List<GeoLibPoint[]> t = GeoWrangler.rectangular_decomposition(T);

            List<GeoLibPoint[]> x = GeoWrangler.rectangular_decomposition(X);

            List<GeoLibPoint[]> s = GeoWrangler.rectangular_decomposition(S);

            Console.WriteLine("Hello World!");
        }
    }
}
