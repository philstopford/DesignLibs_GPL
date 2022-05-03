namespace ClipperLibTest;

using Paths = List<List<ClipperLib1.IntPoint>>;

public static class SubtractionTest
{
    static Clipper2Lib.Paths64 layerAPaths = new()
    {
        new()
        {
            new(0, 0),
            new(0, 900000),
            new(500000, 900000),
            new(500000, 0),
            new(0, 0),
        }
    };


    static Clipper2Lib.Paths64 layerBPaths = new()
    {
        new()
        {
            new(0, 0),
            new(0, 11900),
            new(0, 23800),
            new(0, 35700),
            new(0, 47600),
            new(0, 59500),
            new(12500, 59500),
            new(25000, 59500),
            new(37500, 59500),
            new(50000, 59500),
            new(50000, 50000),
            new(350000, 50000),
            new(350000, 350000),
            new(50000, 350000),
            new(50000, 339661),
            new(50000, 329321),
            new(50000, 318982),
            new(50000, 308643),
            new(50000, 298304),
            new(50000, 287964),
            new(50000, 277625),
            new(50000, 267286),
            new(50000, 256946),
            new(50000, 246607),
            new(50000, 236268),
            new(50000, 225929),
            new(50000, 215589),
            new(50000, 205250),
            new(50000, 194911),
            new(50000, 184571),
            new(50000, 174232),
            new(50000, 163893),
            new(50000, 153554),
            new(50000, 143214),
            new(50000, 132875),
            new(50000, 122536),
            new(50000, 112196),
            new(50000, 101857),
            new(50000, 91518),
            new(50000, 81179),
            new(50000, 70839),
            new(50000, 60500),
            new(37500, 60500),
            new(25000, 60500),
            new(12500, 60500),
            new(0, 60500),
            new(0, 70731),
            new(0, 80962),
            new(0, 91192),
            new(0, 101423),
            new(0, 111654),
            new(0, 121885),
            new(0, 132115),
            new(0, 142346),
            new(0, 152577),
            new(0, 162808),
            new(0, 173038),
            new(0, 183269),
            new(0, 193500),
            new(0, 203731),
            new(0, 213962),
            new(0, 224192),
            new(0, 234423),
            new(0, 244654),
            new(0, 254885),
            new(0, 265115),
            new(0, 275346),
            new(0, 285577),
            new(0, 295808),
            new(0, 306038),
            new(0, 316269),
            new(0, 326500),
            new(0, 336731),
            new(0, 346962),
            new(0, 357192),
            new(0, 367423),
            new(0, 377654),
            new(0, 387885),
            new(0, 398115),
            new(0, 408346),
            new(0, 418577),
            new(0, 428808),
            new(0, 439038),
            new(0, 449269),
            new(0, 459500),
            new(12500, 459500),
            new(25000, 459500),
            new(37500, 459500),
            new(50000, 459500),
            new(50000, 450000),
            new(350000, 450000),
            new(350000, 750000),
            new(50000, 750000),
            new(50000, 739661),
            new(50000, 729321),
            new(50000, 718982),
            new(50000, 708643),
            new(50000, 698304),
            new(50000, 687964),
            new(50000, 677625),
            new(50000, 667286),
            new(50000, 656946),
            new(50000, 646607),
            new(50000, 636268),
            new(50000, 625929),
            new(50000, 615589),
            new(50000, 605250),
            new(50000, 594911),
            new(50000, 584571),
            new(50000, 574232),
            new(50000, 563893),
            new(50000, 553554),
            new(50000, 543214),
            new(50000, 532875),
            new(50000, 522536),
            new(50000, 512196),
            new(50000, 501857),
            new(50000, 491518),
            new(50000, 481179),
            new(50000, 470839),
            new(50000, 460500),
            new(37500, 460500),
            new(25000, 460500),
            new(12500, 460500),
            new(0, 460500),
            new(0, 470721),
            new(0, 480942),
            new(0, 491163),
            new(0, 501384),
            new(0, 511605),
            new(0, 521826),
            new(0, 532047),
            new(0, 542267),
            new(0, 552488),
            new(0, 562709),
            new(0, 572930),
            new(0, 583151),
            new(0, 593372),
            new(0, 603593),
            new(0, 613814),
            new(0, 624035),
            new(0, 634256),
            new(0, 644477),
            new(0, 654698),
            new(0, 664919),
            new(0, 675140),
            new(0, 685360),
            new(0, 695581),
            new(0, 705802),
            new(0, 716023),
            new(0, 726244),
            new(0, 736465),
            new(0, 746686),
            new(0, 756907),
            new(0, 767128),
            new(0, 777349),
            new(0, 787570),
            new(0, 797791),
            new(0, 808012),
            new(0, 818233),
            new(0, 828453),
            new(0, 838674),
            new(0, 848895),
            new(0, 859116),
            new(0, 869337),
            new(0, 879558),
            new(0, 889779),
            new(0, 900000),
            new(500000, 900000),
            new(500000, 0),
            new(0, 0),
            new(-2147483647, -2147483647),
            new(-2147483647, 2147483647),
            new(2147483647, 2147483647),
            new(2147483647, -2147483647),
            new(-2147483647, -2147483647),
            new(0, 0),
        }
    };

    public static void compare()
    {
        Clipper2Lib.Clipper c = new()
        {
            PreserveCollinear = true
        };

        c.AddSubject(layerBPaths);
        c.AddClip(layerAPaths);

        Clipper2Lib.Paths64 solution_cl = new();

        c.Execute(Clipper2Lib.ClipType.Difference, Clipper2Lib.FillRule.EvenOdd, solution_cl);

        c.PreserveCollinear = false;
        
        Clipper2Lib.Paths64 solution_ncl = new();

        c.Execute(Clipper2Lib.ClipType.Difference, Clipper2Lib.FillRule.EvenOdd, solution_ncl);

        ClipperLib1.Clipper c1 = new()
        {
            PreserveCollinear = true
        };

        Paths layerBPaths1 = new();

        foreach (List<Clipper2Lib.Point64> t in layerBPaths)
        {
            List<ClipperLib1.IntPoint> r = new();
            for (int pt = 0; pt < t.Count; pt++)
            {
                r.Add(new(t[pt].X, t[pt].Y));
            }
            
            layerBPaths1.Add(r);
        }

        Paths layerAPaths1 = new();
        foreach (List<Clipper2Lib.Point64> t in layerAPaths)
        {
            List<ClipperLib1.IntPoint> r = new();
            for (int pt = 0; pt < t.Count; pt++)
            {
                r.Add(new(t[pt].X, t[pt].Y));
            }
            
            layerAPaths1.Add(r);
        }

        c1.AddPaths(layerBPaths1, ClipperLib1.PolyType.ptSubject, true);
        c1.AddPaths(layerAPaths1, ClipperLib1.PolyType.ptClip, true);

        Paths solution_cl1 = new();

        c1.Execute(ClipperLib1.ClipType.ctDifference, solution_cl1);

        c1.PreserveCollinear = false;
        
        Paths solution_ncl1 = new();

        c1.Execute(ClipperLib1.ClipType.ctDifference, solution_ncl1);

    }
}