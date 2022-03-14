
using ClipperLib2;

namespace ClipperLib2Test;

using Path64 = List<Point64>;
using Paths64 = List<List<Point64>>;

public static class Clipper2Test
{
    const double keyhole_sizing = 500;

    public static void colinearTest()
    {
        Path64 colinear = new()
        {
            new(-10, -10),
            new(-10, 0),
            new(-10, 10),
            new(0, 10),
            new(10, 10),
            new(10, 0),
            new(10, -10),
            new(0, -10),
            new(-10, -10),
        };

        Clipper c = new() {PreserveCollinear = true};
        c.AddSubject(colinear);
        c.AddClip(colinear);
        Paths64 output = new();
        c.Execute(ClipType.Union, FillRule.EvenOdd, output);
    }

    public static void colinearOffsetTest()
    {
        Path64 colinear = new()
        {
            new(-10, -10),
            new(-10, 0),
            new(-10, 10),
            new(0, 10),
            new(10, 10),
            new(10, 0),
            new(10, -10),
            new(0, -10),
            new(-10, -10),
        };

        ClipperOffset co = new();
        co.AddPath(colinear, JoinType.Miter, EndType.Polygon);
        Paths64 temp = ClipperFunc.Paths64(co.Execute(1.0));

    }

    public static void unionTest()
    {
        Path64 simpleFirstPath = new()
        {
            new(100000, -300000),
            new(100000, 0),
            new(200000, 0),
            new(200000, -300000),
            new(100000, -300000),
        };
        
        Path64 simpleSecondPath = new()
        {
            new (0,-200000),
            new (0,-100000),
            new (300000,-100000),
            new (300000,-200000),
            new (0,-200000),
        };
        
        Clipper cs = new();
        
        cs.AddSubject(simpleFirstPath);
        cs.AddClip(simpleSecondPath);
        
        Paths64 simpleOutputPoints = new();
        cs.Execute(ClipType.Union, FillRule.EvenOdd, simpleOutputPoints);

        Path64 firstPath = new() {
        new (100000,-300000),
        new (100000,-290000),
        new (100000,-280000),
        new (100000,-270000),
        new (100000,-260000),
        new (100000,-250000),
        new (100000,-240000),
        new (100000,-230000),
        new (100000,-220000),
        new (100000,-210000),
        new (100000,-200000),
        new (100000,-190000),
        new (100000,-180000),
        new (100000,-170000),
        new (100000,-160000),
        new (100000,-150000),
        new (100000,-140000),
        new (100000,-130000),
        new (100000,-120000),
        new (100000,-110000),
        new (100000,-100000),
        new (100000,-90000),
        new (100000,-80000),
        new (100000,-70000),
        new (100000,-60000),
        new (100000,-50000),
        new (100000,-40000),
        new (100000,-30000),
        new (100000,-20000),
        new (100000,-10000),
        new (100000,0),
        new (110000,0),
        new (120000,0),
        new (130000,0),
        new (140000,0),
        new (150000,0),
        new (160000,0),
        new (170000,0),
        new (180000,0),
        new (190000,0),
        new (200000,0),
        new (200000,-10000),
        new (200000,-20000),
        new (200000,-30000),
        new (200000,-40000),
        new (200000,-50000),
        new (200000,-60000),
        new (200000,-70000),
        new (200000,-80000),
        new (200000,-90000),
        new (200000,-100000),
        new (200000,-110000),
        new (200000,-120000),
        new (200000,-130000),
        new (200000,-140000),
        new (200000,-150000),
        new (200000,-160000),
        new (200000,-170000),
        new (200000,-180000),
        new (200000,-190000),
        new (200000,-200000),
        new (200000,-210000),
        new (200000,-220000),
        new (200000,-230000),
        new (200000,-240000),
        new (200000,-250000),
        new (200000,-260000),
        new (200000,-270000),
        new (200000,-280000),
        new (200000,-290000),
        new (200000,-300000),
        new (190000,-300000),
        new (180000,-300000),
        new (170000,-300000),
        new (160000,-300000),
        new (150000,-300000),
        new (140000,-300000),
        new (130000,-300000),
        new (120000,-300000),
        new (110000,-300000),
        new (100000,-300000),
        };

        Path64 secondPath = new() {
        new (0,-200000),
        new (0,-190000),
        new (0,-180000),
        new (0,-170000),
        new (0,-160000),
        new (0,-150000),
        new (0,-140000),
        new (0,-130000),
        new (0,-120000),
        new (0,-110000),
        new (0,-100000),
        new (10000,-100000),
        new (20000,-100000),
        new (30000,-100000),
        new (40000,-100000),
        new (50000,-100000),
        new (60000,-100000),
        new (70000,-100000),
        new (80000,-100000),
        new (90000,-100000),
        new (100000,-100000),
        new (110000,-100000),
        new (120000,-100000),
        new (130000,-100000),
        new (140000,-100000),
        new (150000,-100000),
        new (160000,-100000),
        new (170000,-100000),
        new (180000,-100000),
        new (190000,-100000),
        new (200000,-100000),
        new (210000,-100000),
        new (220000,-100000),
        new (230000,-100000),
        new (240000,-100000),
        new (250000,-100000),
        new (260000,-100000),
        new (270000,-100000),
        new (280000,-100000),
        new (290000,-100000),
        new (300000,-100000),
        new (300000,-110000),
        new (300000,-120000),
        new (300000,-130000),
        new (300000,-140000),
        new (300000,-150000),
        new (300000,-160000),
        new (300000,-170000),
        new (300000,-180000),
        new (300000,-190000),
        new (300000,-200000),
        new (290000,-200000),
        new (280000,-200000),
        new (270000,-200000),
        new (260000,-200000),
        new (250000,-200000),
        new (240000,-200000),
        new (230000,-200000),
        new (220000,-200000),
        new (210000,-200000),
        new (200000,-200000),
        new (190000,-200000),
        new (180000,-200000),
        new (170000,-200000),
        new (160000,-200000),
        new (150000,-200000),
        new (140000,-200000),
        new (130000,-200000),
        new (120000,-200000),
        new (110000,-200000),
        new (100000,-200000),
        new (90000,-200000),
        new (80000,-200000),
        new (70000,-200000),
        new (60000,-200000),
        new (50000,-200000),
        new (40000,-200000),
        new (30000,-200000),
        new (20000,-200000),
        new (10000,-200000),
        new (0,-200000),
        };

        Clipper c = new();

        c.AddSubject(firstPath);
        c.AddClip(secondPath);

        Paths64 outputPoints = new();
        c.Execute(ClipType.Union, FillRule.EvenOdd, outputPoints);        
        /*
        
        yields this in ClipperLib2:
        
        outputPoints = {List<List<Point64>>} Count = 1
        new Count = 14
        new (100000,-200000),
        new (100000,-300000),
        new (200000,-300000),
        new (200000,-200000),
        new (300000,-200000),
        new (300000,-100000),
        new (200000,-100000),
        new (200000,0),
        new (100000,0),
        new (100000,-100000),
        new (0,-100000),
        new (0,-200000),
        new (10000,-200000),
        new (-9223372036854765808,-300000),
        */
    }

    public static void notTest()
    {
        Path64 firstPath = new() {
            new (-250000,-250000),
            new (-250000,250000),
            new (250000,250000),
            new (250000,-250000),
            new (-250000,-250000),
        };
        
        Path64 secondPath = new() {
        new (-150000,-150000),
        new (-150000,150000),
        new (150000,150000),
        new (150000,-150000),
        new (-150000,-150000),
        };
        
        
        Clipper c = new();
        
        c.AddSubject(firstPath);
        c.AddClip(secondPath);
        
        Paths64 outputPoints = new();
        
        c.Execute(ClipType.Difference, FillRule.EvenOdd, outputPoints);
    }

    public static void edgeOffsetTest()
    {
        Path64 edge = new()
        {
            new(-100000, 99500),
            new(-100000, 200500)
        };

        ClipperOffset co = new();
        co.AddPath(edge, JoinType.Miter, EndType.Square);
        Paths64 p = ClipperFunc.Paths64(co.Execute(500));

    }
    private static void zFillTest(Point64 bot1, Point64 top1, Point64 bot2, Point64 top2, ref Point64 pt)
    {
        pt.Z = -1;
    }
    public static void zFillCallbackTest()
    {
        Path64 outer = new()
        {
            new(-1000, -1000),
            new(-1000, 1000),
            new(1000, 1000),
            new(1000, -1000),
        };
        
        Path64 cutter = new()
        {
            new(-100, -1100),
            new(-100, -900),
            new(100, -900),
            new(100, -1100),
        };

        Clipper c = new();

        c.ZFillFunc = zFillTest;

        c.AddSubject(outer);
        c.AddClip(cutter);

        Paths64 solution = new();
        c.Execute(ClipType.Intersection, FillRule.EvenOdd, solution);
    }
    public static void coincident_openPathTest()
    {
        Path64 lPoly = new()
        {
            new(-1000, -1000),
            new(-1000, 1000),
            new(1000, 1000),
            new(1000, -1000)
        };

        Path64 t = new()
        {
            new (-1000, -1100),
            new (-1000, 500),
        };

        Clipper c = new();
        c.AddClip(lPoly);
        c.AddOpenSubject(t);

        Paths64 open = new();
        Paths64 solution = new();
        c.Execute(ClipType.Intersection, FillRule.EvenOdd, solution, open);
        
        Path64 t2 = new()
        {
            new (-1000, -1100),
            new (-1000, 500),
            new (-900, 500),
        };

        Clipper c2 = new();
        c2.AddClip(lPoly);
        c2.AddOpenSubject(t2);

        Paths64 open2 = new();
        Paths64 solution2 = new();
        c2.Execute(ClipType.Intersection, FillRule.EvenOdd, solution2, open2);
        
        Path64 t2b = new()
        {
            new (-900, 500),
            new (-1000, 500),
            new (-1000, -1100),
        };

        Clipper c2b = new();
        c2b.AddClip(lPoly);
        c2b.AddOpenSubject(t2b);

        Paths64 open2b = new();
        Paths64 solution2b = new();
        c2b.Execute(ClipType.Intersection, FillRule.EvenOdd, solution2b, open2b);
        
        Path64 t3 = new();
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
        Clipper c3 = new();
        c3.ZFillFunc = zFillTest;
        
        c3.AddClip(lPoly);
        c3.AddOpenSubject(t3);

        Paths64 open3 = new();
        Paths64 solution3 = new();
        c3.Execute(ClipType.Intersection, FillRule.EvenOdd, solution3, open3 );
        
    }
    
    public static void keyHole_test2()
    {
        Path64 lPoly = new()
        {
            new Point64(200000, 0),
            new Point64(200000, 1100000),
            new Point64(1000000, 1100000),
            new Point64(1000000, 800000),
            new Point64(800000, 800000),
            new Point64(800000, 0),
            new Point64(200000, 0)
        };

        Paths64 t = new()
        {
            new Path64()
            {
                new Point64(800000, 800000),
                new Point64(800000, 1100000)
            }
        };
        
        // Turn the new edges into cutters and slice. Not terribly elegant and we're relying on rounding to squash notches later.
        ClipperOffset co = new();
        co.AddPaths(t, JoinType.Miter, EndType.Square);

        Paths64 cutters = ClipperFunc.Paths64(co.Execute(2.0));

        Clipper c = new();

        c.AddSubject(lPoly);

        // Take first cutter only - we only cut once, no matter how many potential cutters we have.
        c.AddClip(cutters[0]);
        Paths64 f = new();
        c.Execute(ClipType.Difference, FillRule.EvenOdd, f);
    }
    public static void openPath_clipTest1()
    {
        Path64 lPoly = new Path64() {
        new Point64(0,0),
        new Point64(0,200000),
        new Point64(200000,200000),
        new Point64(200000,500000),
        new Point64(0,500000),
        new Point64(0,1100000),
        new Point64(1000000,1100000),
        new Point64(1000000,800000),
        new Point64(800000,800000),
        new Point64(800000,600000),
        new Point64(1000000,600000),
        new Point64(1000000,0),
        new Point64(0,0)
        };

        Path64 t = new () {
        new Point64(0,200000),
        new Point64(0,-9800000)
        };

        Clipper c = new();

        c.AddOpenSubject(t);
        c.AddClip(lPoly);

        PolyTree pt = new();
        Paths64 p = new();

        c.Execute(ClipType.Intersection, FillRule.EvenOdd, pt, p);
 
    }
    
    public static void openPath_clipTest2()
    {
        Paths64 rays2 = new()
        {
            new Path64() {new Point64(100000, 200000), new Point64(100000, -9800000)}
        };

        Paths64 rays = new()
        {
            new Path64() {new Point64(100000, 500000), new Point64(10100000, 500000)}
        };

        Paths64 collisionPaths = new()
        {
            new Path64()
            {
                new Point64(0, 0),
                new Point64(0, 500000),
                new Point64(100000, 500000),
                new Point64(100000, 200000),
                new Point64(600000, 200000),
                new Point64(600000, 800000),
                new Point64(1200000, 800000),
                new Point64(1200000, 0),
                new Point64(0, 0)

            }
        };

        Clipper c = new Clipper();
        c.AddOpenSubject(rays);
        c.AddClip(collisionPaths);
        PolyTree pt = new PolyTree();
        Paths64 solution = new();
        c.Execute(ClipType.Intersection, FillRule.EvenOdd, pt, solution);
    }
    
    public static void keyHole_test1()
    {
        Console.WriteLine("Clipper2 Test1");
        Path64 outer = new()
        {
            new Point64(-200000, -200000),
            new Point64(200000, -200000),
            new Point64(200000, 200000),
            new Point64(-200000, 200000),
            new Point64(-200000, -200000)
        };

        Path64 inner1 = new()
        {
            new Point64(-100000, -100000),
            new Point64(-100000, 100000),
            new Point64(100000, 100000),
            new Point64(100000, -100000),
            new Point64(-100000, -100000)
        };

        Paths64 kHSource = new()
        {
            outer,
            inner1
        };
        
        ClipperOffset co = new();
        co.AddPaths(kHSource, JoinType.Miter, EndType.Polygon);
        // ClipperLib2 specifies full width in offset for open path, unlike version 1
        Paths64 out_ = ClipperFunc.Paths64(co.Execute(2*keyhole_sizing));
        
        Console.WriteLine("Out count: " + out_.Count);

    }
    
    public static void offsetTest()
    {
        Path64 lPoly = new ()
        {
            new Point64(0, 0),
            new Point64(0, 500000),
            new Point64(100000, 500000),
            new Point64(100000, 200000),
            new Point64(600000, 200000),
            new Point64(600000, 0),
            new Point64(0, 0)
        };

        Path64 newEdge = new()
        {
            new Point64(100000, 200000),
            new Point64(100000, 0)
        };

        Paths64 newEdges = new()
        {
            newEdge
        };

        ClipperOffset co = new ClipperOffset();
        co.AddPaths(newEdges, JoinType.Miter, EndType.Square);
        PolyTree tp = new();
        // ClipperLib2 specifies full width in offset for open path, unlike version 1
        Paths64 cutters = ClipperFunc.Paths64(co.Execute( 2.0));
        
        Clipper c = new Clipper();
        c.AddSubject(lPoly);
        c.AddClip(cutters);
        c.Execute(ClipType.Difference, FillRule.EvenOdd, tp);
        Paths64 solution = ClipperFunc.PolyTreeToPaths(tp);

    }
    
    public static void leftChordTest()
    {
        Path64 testPath = new List<Point64>() {
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

        Path64 a = new List<Point64>() {
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
        c.ZFillFunc = zFillTest;
        c.AddOpenSubject(testPath);
        c.AddClip(a);

        Paths64 unused = new();
        Paths64 open = new();
        c.Execute(ClipType.Intersection, FillRule.EvenOdd, unused, open);

    }

    public static void rightChordTest()
    {
        Path64 testPath = new ()
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

        Path64 a = new ()
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
        c.ZFillFunc = zFillTest;
        c.AddOpenSubject(testPath);
        c.AddClip(a);

        Paths64 unused = new();
        Paths64 open = new();
        c.Execute(ClipType.Intersection, FillRule.EvenOdd, unused, open);

    }


}