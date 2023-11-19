using NUnit.Framework;

namespace ClipperLibTest;

using Paths64 = Clipper2Lib.Paths64;
using Paths = List<List<ClipperLib1.IntPoint>>;

public static class ComparisonTest
{
    private static readonly Paths64 collisionGeometry = new () {
        new () {
            new (0,-250000),
            new (0,-230000),
            new (0,-210000),
            new (0,-190000),
            new (0,-170000),
            new (0,-150000),
            new (0,-130000),
            new (0,-110000),
            new (0,-90000),
            new (0,-70000),
            new (0,-50000),
            new (4323,-29663),
            new (16543,-12843),
            new (34549,-2447),
            new (50000,0),
            new (70337,-4323),
            new (87157,-16543),
            new (97553,-34549),
            new (100000,-50000),
            new (100000,-70000),
            new (100000,-90000),
            new (100000,-110000),
            new (100000,-130000),
            new (100000,-150000),
            new (100000,-170000),
            new (100000,-190000),
            new (100000,-210000),
            new (100000,-230000),
            new (100000,-250000),
            new (95677,-270337),
            new (83457,-287157),
            new (65451,-297553),
            new (50000,-300000),
            new (29663,-295677),
            new (12843,-283457),
            new (0,-250000),
        },
        new () {
            new (200000,-250000),
            new (200000,-230000),
            new (200000,-210000),
            new (200000,-190000),
            new (200000,-170000),
            new (200000,-150000),
            new (200000,-130000),
            new (200000,-110000),
            new (200000,-90000),
            new (200000,-70000),
            new (200000,-50000),
            new (204323,-29663),
            new (216543,-12843),
            new (234549,-2447),
            new (250000,0),
            new (270337,-4323),
            new (287157,-16543),
            new (297553,-34549),
            new (300000,-50000),
            new (300000,-70000),
            new (300000,-90000),
            new (300000,-110000),
            new (300000,-130000),
            new (300000,-150000),
            new (300000,-170000),
            new (300000,-190000),
            new (300000,-210000),
            new (300000,-230000),
            new (300000,-250000),
            new (295677,-270337),
            new (283457,-287157),
            new (265451,-297553),
            new (250000,-300000),
            new (229663,-295677),
            new (212843,-283457),
            new (200000,-250000),
        }
    };
    public static void lineClipTest1()
    {
        Paths64 castLines = new() {
        new() {
        new (0,-250000),
        new (-291700,-320076),
        },
        new() {
        new (0,-230000),
        new (-300000,-230000),
        },
        new() {
        new (0,-210000),
        new (-300000,-210000),
        },
        new() {
        new (0,-190000),
        new (-300000,-190000),
        },
        new() {
        new (0,-170000),
        new (-300000,-170000),
        },
        new() {
        new (0,-150000),
        new (-300000,-150000),
        },
        new() {
        new (0,-130000),
        new (-300000,-130000),
        },
        new() {
        new (0,-110000),
        new (-300000,-110000),
        },
        new() {
        new (0,-90000),
        new (-300000,-90000),
        },
        new() {
        new (0,-70000),
        new (-300000,-70000),
        },
        new() {
        new (0,-50000),
        new (-298292,-18037),
        },
        new() {
        new (4323,-29663),
        new (-269743,92353),
        },
        new() {
        new (16543,-12843),
        new (-184197,210100),
        },
        new() {
        new (34549,-2447),
        new (-72957,277629),
        },
        new() {
        new (50000,0),
        new (65704,299588),
        },
        new() {
        new (70337,-4323),
        new (192352,269743),
        },
        new() {
        new (87157,-16543),
        new (310099,184197),
        },
        new() {
        new (97553,-34549),
        new (377628,72957),
        },
        new() {
        new (100000,-50000),
        new (399288,-29350),
        },
        new() {
        new (100000,-70000),
        new (400000,-70000),
        },
        new() {
        new (100000,-90000),
        new (400000,-90000),
        },
        new() {
        new (100000,-110000),
        new (400000,-110000),
        },
        new() {
        new (100000,-130000),
        new (400000,-130000),
        },
        new() {
        new (100000,-150000),
        new (400000,-150000),
        },
        new() {
        new (100000,-170000),
        new (400000,-170000),
        },
        new() {
        new (100000,-190000),
        new (400000,-190000),
        },
        new() {
        new (100000,-210000),
        new (400000,-210000),
        },
        new() {
        new (100000,-230000),
        new (400000,-230000),
        },
        new() {
        new (100000,-250000),
        new (398292,-281963),
        },
        new() {
        new (95677,-270337),
        new (369743,-392353),
        },
        new() {
        new (83457,-287157),
        new (284197,-510100),
        },
        new() {
        new (65451,-297553),
        new (172957,-577629),
        },
        new() {
        new (50000,-300000),
        new (34296,-599588),
        },
        new() {
        new (29663,-295677),
        new (-92352,-569743),
        },
        new() {
        new (12843,-283457),
        new (-238759,-446847),
        },
        new() {
        new (0,-250000),
        new (-291700,-320076),}
        };

        // Create our Clipper1 collision geometry from the Clipper2 definition
        Paths collsionGeometry1 = new();
        foreach (List<Clipper2Lib.Point64> t in collisionGeometry)
        {
            List<ClipperLib1.IntPoint> a = new List<ClipperLib1.IntPoint>();
            for (int j = 0; j < t.Count; j++)
            {
                a.Add(new ClipperLib1.IntPoint(t[j].X, t[j].Y));
            }
            collsionGeometry1.Add(a);
        }

        foreach (Clipper2Lib.Path64 t in castLines)
        {
            Clipper2Lib.Clipper64 c2 = new();
            c2.AddOpenSubject(t);
            c2.AddClip(collisionGeometry);
            Paths64 unused2 = new();
            Paths64 clipped2 = new();
            c2.Execute(Clipper2Lib.ClipType.Difference, Clipper2Lib.FillRule.EvenOdd, unused2, clipped2);

            ClipperLib1.Clipper c1 = new();
            List<ClipperLib1.IntPoint> c1Ray = new()
                {
                    new (t[0].X, t[0].Y),
                    new (t[1].X, t[1].Y),
                }
                ;
            c1.AddPath(c1Ray, ClipperLib1.PolyType.ptSubject, false);
            c1.AddPaths(collsionGeometry1, ClipperLib1.PolyType.ptClip, true);
            ClipperLib1.PolyTree pt = new();
            c1.Execute(ClipperLib1.ClipType.ctDifference, pt);
            Paths clipped1 = ClipperLib1.Clipper.OpenPathsFromPolyTree(pt);
            
            // Results ought to be ~consistent. ClipperLib2 reverses the direction compared to ClipperLib1
            Assert.True(((clipped1[0][0].X == clipped2[0][0].X) && (clipped1[0][0].Y == clipped2[0][0].Y)) ||
                        ((clipped1[0][0].X == clipped2[0][1].X) && (clipped1[0][0].Y == clipped2[0][1].Y)));
        }
    }

    public static void lineClipTest2()
    {
        Paths64 castLines = new() {
            new() {
                new (200000,-250000),
                new (-91700,-320076),
            },
            new() {
                new (200000,-230000),
                new (-100000,-230000),
            },
            new() {
                new (200000,-210000),
                new (-100000,-210000),
            },
            new() {
                new (200000,-190000),
                new (-100000,-190000),
            },
            new() {
                new (200000,-170000),
                new (-100000,-170000),
            },
            new() {
                new (200000,-150000),
                new (-100000,-150000),
            },
            new() {
                new (200000,-130000),
                new (-100000,-130000),
            },
            new() {
                new (200000,-110000),
                new (-100000,-110000),
            },
            new() {
                new (200000,-90000),
                new (-100000,-90000),
            },
            new() {
                new (200000,-70000),
                new (-100000,-70000),
            },
            new() {
                new (200000,-50000),
                new (-98292,-18037),
            },
            new() {
                new (204323,-29663),
                new (-69743,92353),
            },
            new() {
                new (216543,-12843),
                new (15803,210100),
            },
            new() {
                new (234549,-2447),
                new (127043,277629),
            },
            new() {
                new (250000,0),
                new (265704,299588),
            },
            new() {
                new (270337,-4323),
                new (392352,269743),
            },
            new() {
                new (287157,-16543),
                new (510099,184197),
            },
            new() {
                new (297553,-34549),
                new (577628,72957),
            },
            new() {
                new (300000,-50000),
                new (599288,-29350),
            },
            new() {
                new (300000,-70000),
                new (600000,-70000),
            },
            new() {
                new (300000,-90000),
                new (600000,-90000),
            },
            new() {
                new (300000,-110000),
                new (600000,-110000),
            },
            new() {
                new (300000,-130000),
                new (600000,-130000),
            },
            new() {
                new (300000,-150000),
                new (600000,-150000),
            },
            new() {
                new (300000,-170000),
                new (600000,-170000),
            },
            new() {
                new (300000,-190000),
                new (600000,-190000),
            },
            new() {
                new (300000,-210000),
                new (600000,-210000),
            },
            new() {
                new (300000,-230000),
                new (600000,-230000),
            },
            new() {
                new (300000,-250000),
                new (598292,-281963),
            },
            new() {
                new (295677,-270337),
                new (569743,-392353),
            },
            new() {
                new (283457,-287157),
                new (484197,-510100),
            },
            new() {
                new (265451,-297553),
                new (372957,-577629),
            },
            new() {
                new (250000,-300000),
                new (234296,-599588),
            },
            new() {
                new (229663,-295677),
                new (107648,-569743),
            },
            new() {
                new (212843,-283457),
                new (-38759,-446847),
            },
            new() {
                new (200000,-250000),
                new (-91700,-320076)
            }
        };        
        
        // Create our Clipper1 collision geometry from the Clipper2 definition
        Paths collsionGeometry1 = new();
        foreach (List<Clipper2Lib.Point64> t in collisionGeometry)
        {
            List<ClipperLib1.IntPoint> a = new List<ClipperLib1.IntPoint>();
            for (int j = 0; j < t.Count; j++)
            {
                a.Add(new ClipperLib1.IntPoint(t[j].X, t[j].Y));
            }
            collsionGeometry1.Add(a);
        }

        foreach (Clipper2Lib.Path64 t in castLines)
        {
            Clipper2Lib.Clipper64 c2 = new();
            c2.AddOpenSubject(t);
            c2.AddClip(collisionGeometry);
            Paths64 unused2 = new();
            Paths64 clipped2 = new();
            c2.Execute(Clipper2Lib.ClipType.Difference, Clipper2Lib.FillRule.EvenOdd, unused2, clipped2);

            ClipperLib1.Clipper c1 = new();
            List<ClipperLib1.IntPoint> c1Ray = new()
                {
                    new (t[0].X, t[0].Y),
                    new (t[1].X, t[1].Y),
                }
                ;
            c1.AddPath(c1Ray, ClipperLib1.PolyType.ptSubject, false);
            c1.AddPaths(collsionGeometry1, ClipperLib1.PolyType.ptClip, true);
            ClipperLib1.PolyTree pt = new();
            c1.Execute(ClipperLib1.ClipType.ctDifference, pt);
            Paths clipped1 = ClipperLib1.Clipper.OpenPathsFromPolyTree(pt);
            
            // Results ought to be ~consistent. ClipperLib2 reverses the direction compared to ClipperLib1
            Assert.True(((clipped1[0][0].X == clipped2[0][0].X) && (clipped1[0][0].Y == clipped2[0][0].Y)) ||
                        ((clipped1[0][0].X == clipped2[0][1].X) && (clipped1[0][0].Y == clipped2[0][1].Y)));
        }
    }
}