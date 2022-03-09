
using ClipperLib1;

namespace ClipperLib1Test;

using Path = List<IntPoint>;
using Paths = List<List<IntPoint>>;

public static class Clipper1Test
{
    const double keyhole_sizing = 500;

    public static void notTest()
    {
        Path firstPath = new() {
        new (-250000,-250000),
        new (-250000,250000),
        new (250000,250000),
        new (250000,-250000),
        new (-250000,-250000),
        };
        
        Path secondPath = new() {
        new (-150000,-150000),
        new (-150000,150000),
        new (150000,150000),
        new (150000,-150000),
        new (-150000,-150000),
        new (-2147483647,-2147483647),
        new (-2147483647,2147483647),
        new (2147483647,2147483647),
        new (2147483647,-2147483647),
        new (-2147483647,-2147483647),
        new (-150000,-150000),
        };
        
        
        Clipper c = new();
        
        c.AddPath(firstPath, PolyType.ptSubject, true);
        c.AddPath(secondPath,PolyType.ptClip, true);
        
        Paths outputPoints = new();
        
        c.Execute(ClipType.ctIntersection, outputPoints);
    }
    public static void edgeOffsetTest()
    {
        Path edge = new()
        {
            new(-100000, 99500),
            new(-100000, 200500)
        };

        ClipperOffset co = new();
        co.AddPath(edge, JoinType.jtMiter, EndType.etOpenSquare);
        Paths p = new();
        co.Execute(ref p, 500);

    }
    private static void zFillTest(IntPoint bot1, IntPoint top1, IntPoint bot2, IntPoint top2, ref IntPoint pt)
    {
        pt.Z = -1;
    }
    public static void zFillCallbackTest()
    {
        Path outer = new()
        {
            new(-1000, -1000),
            new(-1000, 1000),
            new(1000, 1000),
            new(1000, -1000)
        };
        
        Path cutter = new()
        {
            new(-100, -1100),
            new(-100, -900),
            new(100, -900),
            new(100, -1100),
        };

        Clipper c = new();

        c.ZFillFunction = zFillTest;

        c.AddPath(outer, PolyType.ptSubject, true);
        c.AddPath(cutter, PolyType.ptClip, true);

        Paths solution = new();
        c.Execute(ClipType.ctIntersection, solution);
    }
    
    public static void coincident_openPathTest()
    {
        Path lPoly = new()
        {
            new(-1000, -1000),
            new(-1000, 1000),
            new(1000, 1000),
            new(1000, -1000)
        };

        Path t = new()
        {
            new (-1000, -1100),
            new (-1000, 500),
        };

        Clipper c = new()
        {
            PreserveCollinear = true,
            StrictlySimple = false
        };
        c.AddPath(lPoly, PolyType.ptClip, true);
        c.AddPath(t, PolyType.ptSubject, false);

        PolyTree pt = new();
        c.Execute(ClipType.ctIntersection, pt);
        Paths solution = Clipper.OpenPathsFromPolyTree(pt);
        
        Path t2 = new()
        {
            new (-1000, -1100),
            new (-1000, 500),
            new (-900, 500),
        };

        Clipper c2 = new()
        {
            PreserveCollinear = true,
            StrictlySimple = false
        };
        
        c2.AddPath(lPoly, PolyType.ptClip, true);
        c2.AddPath(t2, PolyType.ptSubject, false);

        PolyTree pt2 = new();
        c2.Execute(ClipType.ctIntersection, pt2);
        Paths solution2 = Clipper.OpenPathsFromPolyTree(pt2);

        Path t2b = new()
        {
            new (-900, 500),
            new (-1000, 500),
            new (-1000, -1100),
        };

        Clipper c2b = new()
        {
            PreserveCollinear = true,
            StrictlySimple = false
        };
        
        c2b.AddPath(lPoly, PolyType.ptClip, true);
        c2b.AddPath(t2b, PolyType.ptSubject, false);

        PolyTree pt2b = new();
        c2b.Execute(ClipType.ctIntersection, pt2b);
        Paths solution2b = Clipper.OpenPathsFromPolyTree(pt2b);


        Path t3 = new();
        int x = 0;
        int y = -1100;
        while (y < 1200)
        {
            t3.Add(new()
            {
                X = x,
                Y = y
            });
            y += 100;
        }
        Clipper c3 = new()
        {
            PreserveCollinear = true,
            StrictlySimple = false
        };
        c3.ZFillFunction = zFillTest;
        c3.AddPath(lPoly, PolyType.ptClip, true);
        c3.AddPath(t3, PolyType.ptSubject, false);

        PolyTree pt3 = new();
        c3.Execute(ClipType.ctIntersection, pt3);
        Paths solution3 = Clipper.OpenPathsFromPolyTree(pt3);

    }
    public static void keyHole_test2()
    {
        Path lPoly = new()
        {
            new IntPoint(200000, 0),
            new IntPoint(200000, 1100000),
            new IntPoint(1000000, 1100000),
            new IntPoint(1000000, 800000),
            new IntPoint(800000, 800000),
            new IntPoint(800000, 0),
            new IntPoint(200000, 0)
        };

        Paths t = new()
        {
            new Path()
            {
                new IntPoint(800000, 800000),
                new IntPoint(800000, 1100000)
            }
        };
        
        // Turn the new edges into cutters and slice. Not terribly elegant and we're relying on rounding to squash notches later.
        ClipperOffset co = new();
        co.AddPaths(t, JoinType.jtMiter, EndType.etOpenSquare);
        PolyTree tp = new();
        co.Execute(ref tp, 1.0);

        Paths cutters = Clipper.ClosedPathsFromPolyTree(tp);

        Clipper c = new();

        c.AddPath(lPoly, PolyType.ptSubject, true);

        // Take first cutter only - we only cut once, no matter how many potential cutters we have.
        c.AddPath(cutters[0], PolyType.ptClip, true);
        Paths f = new();
        c.Execute(ClipType.ctDifference, f, PolyFillType.pftEvenOdd, PolyFillType.pftEvenOdd);
    }
    
    public static void openPath_clipTest1()
    {
        Path lPoly = new() {
        new IntPoint(0,0),
        new IntPoint(0,200000),
        new IntPoint(200000,200000),
        new IntPoint(200000,500000),
        new IntPoint(0,500000),
        new IntPoint(0,1100000),
        new IntPoint(1000000,1100000),
        new IntPoint(1000000,800000),
        new IntPoint(800000,800000),
        new IntPoint(800000,600000),
        new IntPoint(1000000,600000),
        new IntPoint(1000000,0),
        new IntPoint(0,0)
        };

        Path t = new () {
        new IntPoint(0,200000),
        new IntPoint(0,-9800000)
        };

        Clipper c = new();

        c.AddPath(t, PolyType.ptSubject, false);
        c.AddPath(lPoly, PolyType.ptClip, true);

        PolyTree pt = new();
        Paths p = new();

        c.Execute(ClipType.ctIntersection, pt);
        Paths solution = Clipper.OpenPathsFromPolyTree(pt);
 
    }

    public static void keyHole_test1()
    {
        Console.WriteLine("Clipper1 Test1");
        Path outer = new()
        {
            new IntPoint(-200000, -200000),
            new IntPoint(200000, -200000),
            new IntPoint(200000, 200000),
            new IntPoint(-200000, 200000),
            new IntPoint(-200000, -200000)
        };

        Path inner1 = new()
        {
            new IntPoint(-100000, -100000),
            new IntPoint(-100000, 100000),
            new IntPoint(100000, 100000),
            new IntPoint(100000, -100000),
            new IntPoint(-100000, -100000)
        };
        Paths kHSource = new()
        {
            outer,
            inner1
        };
        
        ClipperOffset co = new();
        co.AddPaths(kHSource, JoinType.jtMiter, EndType.etClosedPolygon);
        Paths out_ = new ();
        co.Execute(ref out_, keyhole_sizing);

        Console.WriteLine("Out count: " + out_.Count);
        
    }

    public static void openPath_clipTest2()
    {
        Paths rays = new()
        {
            new Path() {new IntPoint(100000, 200000), new IntPoint(100000, -9800000)}
        };

        Paths collisionPaths = new()
        {
            new Path()
            {
                new IntPoint(0, 0),
                new IntPoint(0, 500000),
                new IntPoint(100000, 500000),
                new IntPoint(100000, 200000),
                new IntPoint(600000, 200000),
                new IntPoint(600000, 800000),
                new IntPoint(1200000, 800000),
                new IntPoint(1200000, 0),
                new IntPoint(0, 0)

            }
        };

        Clipper c = new Clipper();
        c.AddPaths(rays, PolyType.ptSubject, false);
        c.AddPaths(collisionPaths, PolyType.ptClip, true);
        PolyTree pt = new PolyTree();
        c.Execute(ClipType.ctIntersection, pt);
        Paths solution = Clipper.OpenPathsFromPolyTree(pt);
    }
    public static void offsetTest()
    {
        Path lPoly = new ()
        {
            new IntPoint(0, 0),
            new IntPoint(0, 500000),
            new IntPoint(100000, 500000),
            new IntPoint(100000, 200000),
            new IntPoint(600000, 200000),
            new IntPoint(600000, 0),
            new IntPoint(0, 0)
        };

        Path newEdge = new()
        {
            new IntPoint(100000, 200000),
            new IntPoint(100000, 0)
        };

        Paths newEdges = new()
        {
            newEdge
        };

        ClipperOffset co = new ClipperOffset();
        co.AddPaths(newEdges, JoinType.jtMiter, EndType.etOpenSquare);
        PolyTree tp = new();
        co.Execute(ref tp, 1.0);

        Paths cutters = Clipper.ClosedPathsFromPolyTree(tp);

        Clipper c = new Clipper();
        c.AddPath(lPoly, PolyType.ptSubject, true);
        c.AddPaths(cutters, PolyType.ptClip, true);
        Paths solution = new();
        c.Execute(ClipType.ctDifference, solution);

    }
    
    public static void leftChordTest()
    {
        Path testPath = new List<IntPoint>() {
         new(-200000,0),
         new(-300000,0),
         new(-310453,-548),
         new(-320791,-2185),
         new(-330902,-4894),
         new(-340674,-8645),
         new(-350000,-13397),
         new(-358779,-19098),
         new(-366913,-25686),
         new(-374314,-33087),
         new(-380902,-41221),
         new(-386603,-50000),
         new(-391355,-59326),
         new(-395106,-69098),
         new(-397815,-79209),
         new(-399452,-89547),
         new(-400000,-100000),
         new(-400000,-700000),
         new(-400000,-700000),
         new(-398907,-710396),
         new(-397815,-720791),
         new(-395106,-730902),
         new(-391355,-740674),
         new(-386603,-750000),
         new(-380902,-758779),
         new(-374314,-766913),
         new(-366913,-774314),
         new(-358779,-780902),
         new(-350000,-786603),
         new(-340674,-791355),
         new(-330902,-795106),
         new(-320791,-797815),
         new(-310453,-799452),
         new(-300000,-800000),
         };

        Path a = new List<IntPoint>() {
          new(-460000,-390000),
          new(-459955,-395235),
          new(-459594,-405679),
          new(-458873,-416047),
          new(-457796,-426288),
          new(-456369,-436353),
          new(-454102,-448610),
          new(-451315,-460421),
          new(-448030,-471696),
          new(-444272,-482349),
          new(-440068,-492300),
          new(-435451,-501472),
          new(-429416,-511353),
          new(-422903,-519904),
          new(-414796,-528076),
          new(-406257,-534189),
          new(-396132,-538540),
          new(-385806,-540000),
          new(-381067,-540036),
          new(-369255,-540441),
          new(-357570,-541292),
          new(-346100,-542583),
          new(-334932,-544305),
          new(-324151,-546443),
          new(-313840,-548983),
          new(-304076,-551904),
          new(-293187,-555882),
          new(-283312,-560333),
          new(-273218,-566059),
          new(-264802,-572279),
          new(-257399,-579871),
          new(-252063,-588852),
          new(-250000,-599117),
          new(-242658,-614841),
          new(-234819,-621423),
          new(-225801,-626830),
          new(-216572,-631139),
          new(-206066,-635097),
          new(-196418,-638096),
          new(-186036,-640798),
          new(-175000,-643183),
          new(-163393,-645233),
          new(-151303,-646931),
          new(-141346,-648029),
          new(-131187,-648888),
          new(-120876,-649505),
          new(-110463,-649876),
          new(-100000,-650000),
          new(-94765,-649967),
          new(-84321,-649701),
          new(-73953,-649169),
          new(-63712,-648376),
          new(-53647,-647324),
          new(-41390,-645654),
          new(-29579,-643601),
          new(-18304,-641180),
          new(-7651,-638411),
          new(2300,-635314),
          new(13206,-631197),
          new(22873,-626688),
          new(32442,-620997),
          new(40954,-614029),
          new(47244,-605763),
          new(50000,-595332),
          new(50438,-591472),
          new(55347,-581946),
          new(63107,-574604),
          new(72568,-568506),
          new(82553,-563595),
          new(92112,-559765),
          new(102721,-556206),
          new(114298,-552945),
          new(124199,-550567),
          new(134615,-548408),
          new(145495,-546477),
          new(156787,-544784),
          new(168436,-543337),
          new(180385,-542143),
          new(192576,-541209),
          new(204949,-540538),
          new(217444,-540135),
          new(230000,-540000),
          new(237542,-540000),
          new(241117,-539909),
          new(251801,-538540),
          new(262329,-535544),
          new(272585,-530954),
          new(282456,-524819),
          new(290312,-518575),
          new(297765,-511353),
          new(304760,-503206),
          new(311244,-494199),
          new(317167,-484398),
          new(322483,-473879),
          new(327154,-462721),
          new(331142,-451010),
          new(333821,-441303),
          new(336031,-431346),
          new(337761,-421187),
          new(339003,-410876),
          new(339750,-400463),
          new(340000,-390000),
          new(339909,-384765),
          new(339178,-374321),
          new(337721,-363953),
          new(335544,-353712),
          new(332658,-343647),
          new(329078,-333809),
          new(324819,-324244),
          new(319904,-315000),
          new(314356,-306121),
          new(308202,-297651),
          new(301472,-289630),
          new(294199,-282099),
          new(286418,-275093),
          new(278168,-268647),
          new(269488,-262793),
          new(260421,-257558),
          new(251010,-252968),
          new(241303,-249046),
          new(231346,-245811),
          new(221187,-243278),
          new(210876,-241460),
          new(200463,-240365),
          new(190000,-240000),
          new(189927,-240000),
          new(180166,-239562),
          new(168038,-237784),
          new(156076,-234653),
          new(146687,-231190),
          new(137509,-226893),
          new(128587,-221783),
          new(119964,-215885),
          new(111681,-209227),
          new(103779,-201842),
          new(96298,-193766),
          new(89272,-185039),
          new(82737,-175702),
          new(76724,-165801),
          new(71262,-155385),
          new(66379,-144505),
          new(62097,-133213),
          new(58439,-121564),
          new(55421,-109615),
          new(53058,-97424),
          new(51362,-85051),
          new(50341,-72556),
          new(50000,-60000),
          new(50000,-30000),
          new(49909,-24765),
          new(49178,-14321),
          new(47721,-3953),
          new(45544,6288),
          new(42658,16353),
          new(39078,26191),
          new(34819,35756),
          new(29904,45000),
          new(24356,53879),
          new(18202,62349),
          new(11472,70370),
          new(4199,77901),
          new(-3582,84907),
          new(-11832,91353),
          new(-20512,97207),
          new(-29579,102442),
          new(-38990,107032),
          new(-48697,110954),
          new(-58654,114189),
          new(-68813,116722),
          new(-79124,118540),
          new(-89537,119635),
          new(-100000,120000),
          new(-105235,119909),
          new(-115679,119178),
          new(-126047,117721),
          new(-136288,115544),
          new(-146353,112658),
          new(-156191,109078),
          new(-165756,104819),
          new(-175000,99904),
          new(-183879,94356),
          new(-192349,88202),
          new(-200370,81472),
          new(-207901,74199),
          new(-214907,66418),
          new(-221353,58168),
          new(-227207,49488),
          new(-232442,40421),
          new(-237032,31010),
          new(-240954,21303),
          new(-244189,11346),
          new(-246722,1187),
          new(-248540,-9124),
          new(-249635,-19537),
          new(-250000,-30000),
          new(-250000,-60000),
          new(-250062,-66282),
          new(-250555,-78815),
          new(-251539,-91257),
          new(-253010,-103546),
          new(-254959,-115623),
          new(-257378,-127429),
          new(-260255,-138907),
          new(-263575,-150000),
          new(-267323,-160655),
          new(-271480,-170819),
          new(-276026,-180444),
          new(-280939,-189481),
          new(-287560,-199886),
          new(-294665,-209227),
          new(-302202,-217432),
          new(-310113,-224438),
          new(-320015,-231190),
          new(-330260,-236067),
          new(-340735,-239014),
          new(-351326,-240000),
          new(-357014,-240206),
          new(-368327,-241847),
          new(-379453,-245111),
          new(-390272,-249963),
          new(-398966,-255181),
          new(-407297,-261425),
          new(-415203,-268647),
          new(-422623,-276794),
          new(-429500,-285801),
          new(-435782,-295602),
          new(-441421,-306121),
          new(-446374,-317279),
          new(-450605,-328990),
          new(-453446,-338697),
          new(-455790,-348654),
          new(-457625,-358813),
          new(-458942,-369124),
          new(-459735,-379537),
          new(-460000,-390000),
          };

        // Let's see if we can track the origin of the chords.
        for (int p = 0; p < a.Count; p++)
        {
            a[p] = new()
            {
                X = a[p].X,
                Y = a[p].Y,
                Z = 1
            };
        }

        for (int p = 0; p < testPath.Count; p++)
        {
            testPath[p] = new()
            {
                X = testPath[p].X,
                Y = testPath[p].Y,
                Z = 2
            };
        }
        
        Clipper c = new();
        c.ZFillFunction = zFillTest;
        c.AddPath(testPath, PolyType.ptSubject, false);
        c.AddPath(a, PolyType.ptClip, true);

        PolyTree pt = new();
        c.Execute(ClipType.ctIntersection, pt);
        Paths open = Clipper.OpenPathsFromPolyTree(pt);
    }
    
    public static void rightChordTest()
    {
        Path testPath = new ()
        {
            new(-400000, -700000),
            new(-398907, -710396),
            new(-397815, -720791),
            new(-395106, -730902),
            new(-391355, -740674),
            new(-386603, -750000),
            new(-380902, -758779),
            new(-374314, -766913),
            new(-366913, -774314),
            new(-358779, -780902),
            new(-350000, -786603),
            new(-340674, -791355),
            new(-330902, -795106),
            new(-320791, -797815),
            new(-310453, -799452),
            new(-300000, -800000),
            new(-200000, -800000),
            new(-189547, -799452),
            new(-179209, -797815),
            new(-169098, -795106),
            new(-159326, -791355),
            new(-150000, -786603),
            new(-141221, -780902),
            new(-133087, -774314),
            new(-125686, -766913),
            new(-119098, -758779),
            new(-113397, -750000),
            new(-108645, -740674),
            new(-104894, -730902),
            new(-102185, -720791),
            new(-100548, -710453),
            new(-100000, -700000),
            new(-100000, -100000),
            new(-100548, -89547),
            new(-102185, -79209),
            new(-104894, -69098),
            new(-108645, -59326),
            new(-113397, -50000),
            new(-119098, -41221),
            new(-125686, -33087),
            new(-133087, -25686),
            new(-141221, -19098),
            new(-150000, -13397),
            new(-159326, -8645),
            new(-169098, -4894),
            new(-179209, -2185),
            new(-189547, -548),
            new(-200000, 0)
        };

        Path a = new ()
        {
            new(-460000, -390000),
            new(-459955, -395235),
            new(-459594, -405679),
            new(-458873, -416047),
            new(-457796, -426288),
            new(-456369, -436353),
            new(-454102, -448610),
            new(-451315, -460421),
            new(-448030, -471696),
            new(-444272, -482349),
            new(-440068, -492300),
            new(-435451, -501472),
            new(-429416, -511353),
            new(-422903, -519904),
            new(-414796, -528076),
            new(-406257, -534189),
            new(-396132, -538540),
            new(-385806, -540000),
            new(-381067, -540036),
            new(-369255, -540441),
            new(-357570, -541292),
            new(-346100, -542583),
            new(-334932, -544305),
            new(-324151, -546443),
            new(-313840, -548983),
            new(-304076, -551904),
            new(-293187, -555882),
            new(-283312, -560333),
            new(-273218, -566059),
            new(-264802, -572279),
            new(-257399, -579871),
            new(-252063, -588852),
            new(-250000, -599117),
            new(-242658, -614841),
            new(-234819, -621423),
            new(-225801, -626830),
            new(-216572, -631139),
            new(-206066, -635097),
            new(-196418, -638096),
            new(-186036, -640798),
            new(-175000, -643183),
            new(-163393, -645233),
            new(-151303, -646931),
            new(-141346, -648029),
            new(-131187, -648888),
            new(-120876, -649505),
            new(-110463, -649876),
            new(-100000, -650000),
            new(-94765, -649967),
            new(-84321, -649701),
            new(-73953, -649169),
            new(-63712, -648376),
            new(-53647, -647324),
            new(-41390, -645654),
            new(-29579, -643601),
            new(-18304, -641180),
            new(-7651, -638411),
            new(2300, -635314),
            new(13206, -631197),
            new(22873, -626688),
            new(32442, -620997),
            new(40954, -614029),
            new(47244, -605763),
            new(50000, -595332),
            new(50438, -591472),
            new(55347, -581946),
            new(63107, -574604),
            new(72568, -568506),
            new(82553, -563595),
            new(92112, -559765),
            new(102721, -556206),
            new(114298, -552945),
            new(124199, -550567),
            new(134615, -548408),
            new(145495, -546477),
            new(156787, -544784),
            new(168436, -543337),
            new(180385, -542143),
            new(192576, -541209),
            new(204949, -540538),
            new(217444, -540135),
            new(230000, -540000),
            new(237542, -540000),
            new(241117, -539909),
            new(251801, -538540),
            new(262329, -535544),
            new(272585, -530954),
            new(282456, -524819),
            new(290312, -518575),
            new(297765, -511353),
            new(304760, -503206),
            new(311244, -494199),
            new(317167, -484398),
            new(322483, -473879),
            new(327154, -462721),
            new(331142, -451010),
            new(333821, -441303),
            new(336031, -431346),
            new(337761, -421187),
            new(339003, -410876),
            new(339750, -400463),
            new(340000, -390000),
            new(339909, -384765),
            new(339178, -374321),
            new(337721, -363953),
            new(335544, -353712),
            new(332658, -343647),
            new(329078, -333809),
            new(324819, -324244),
            new(319904, -315000),
            new(314356, -306121),
            new(308202, -297651),
            new(301472, -289630),
            new(294199, -282099),
            new(286418, -275093),
            new(278168, -268647),
            new(269488, -262793),
            new(260421, -257558),
            new(251010, -252968),
            new(241303, -249046),
            new(231346, -245811),
            new(221187, -243278),
            new(210876, -241460),
            new(200463, -240365),
            new(190000, -240000),
            new(189927, -240000),
            new(180166, -239562),
            new(168038, -237784),
            new(156076, -234653),
            new(146687, -231190),
            new(137509, -226893),
            new(128587, -221783),
            new(119964, -215885),
            new(111681, -209227),
            new(103779, -201842),
            new(96298, -193766),
            new(89272, -185039),
            new(82737, -175702),
            new(76724, -165801),
            new(71262, -155385),
            new(66379, -144505),
            new(62097, -133213),
            new(58439, -121564),
            new(55421, -109615),
            new(53058, -97424),
            new(51362, -85051),
            new(50341, -72556),
            new(50000, -60000),
            new(50000, -30000),
            new(49909, -24765),
            new(49178, -14321),
            new(47721, -3953),
            new(45544, 6288),
            new(42658, 16353),
            new(39078, 26191),
            new(34819, 35756),
            new(29904, 45000),
            new(24356, 53879),
            new(18202, 62349),
            new(11472, 70370),
            new(4199, 77901),
            new(-3582, 84907),
            new(-11832, 91353),
            new(-20512, 97207),
            new(-29579, 102442),
            new(-38990, 107032),
            new(-48697, 110954),
            new(-58654, 114189),
            new(-68813, 116722),
            new(-79124, 118540),
            new(-89537, 119635),
            new(-100000, 120000),
            new(-105235, 119909),
            new(-115679, 119178),
            new(-126047, 117721),
            new(-136288, 115544),
            new(-146353, 112658),
            new(-156191, 109078),
            new(-165756, 104819),
            new(-175000, 99904),
            new(-183879, 94356),
            new(-192349, 88202),
            new(-200370, 81472),
            new(-207901, 74199),
            new(-214907, 66418),
            new(-221353, 58168),
            new(-227207, 49488),
            new(-232442, 40421),
            new(-237032, 31010),
            new(-240954, 21303),
            new(-244189, 11346),
            new(-246722, 1187),
            new(-248540, -9124),
            new(-249635, -19537),
            new(-250000, -30000),
            new(-250000, -60000),
            new(-250062, -66282),
            new(-250555, -78815),
            new(-251539, -91257),
            new(-253010, -103546),
            new(-254959, -115623),
            new(-257378, -127429),
            new(-260255, -138907),
            new(-263575, -150000),
            new(-267323, -160655),
            new(-271480, -170819),
            new(-276026, -180444),
            new(-280939, -189481),
            new(-287560, -199886),
            new(-294665, -209227),
            new(-302202, -217432),
            new(-310113, -224438),
            new(-320015, -231190),
            new(-330260, -236067),
            new(-340735, -239014),
            new(-351326, -240000),
            new(-357014, -240206),
            new(-368327, -241847),
            new(-379453, -245111),
            new(-390272, -249963),
            new(-398966, -255181),
            new(-407297, -261425),
            new(-415203, -268647),
            new(-422623, -276794),
            new(-429500, -285801),
            new(-435782, -295602),
            new(-441421, -306121),
            new(-446374, -317279),
            new(-450605, -328990),
            new(-453446, -338697),
            new(-455790, -348654),
            new(-457625, -358813),
            new(-458942, -369124),
            new(-459735, -379537),
            new(-460000, -390000),
        };

        // Let's see if we can track the origin of the chords.
        for (int p = 0; p < a.Count; p++)
        {
            a[p] = new()
            {
                X = a[p].X,
                Y = a[p].Y,
                Z = 1
            };
        }

        for (int p = 0; p < testPath.Count; p++)
        {
            testPath[p] = new()
            {
                X = testPath[p].X,
                Y = testPath[p].Y,
                Z = 2
            };
        }

        Clipper c = new();
        c.ZFillFunction = zFillTest;
        c.AddPath(testPath, PolyType.ptSubject, false);
        c.AddPath(a, PolyType.ptClip, true);

        PolyTree pt = new();
        c.Execute(ClipType.ctIntersection, pt);
        Paths open = Clipper.OpenPathsFromPolyTree(pt);
    }

}