using Clipper2Lib;

namespace PartitionTestGeometrySource;

public static partial class TestGeometry
{
    public static PathD getL()
    {
        // L
        return Clipper.MakePath(new double[] {
            0, 0,
            0, 50,
            10, 50,
            10, 20,
            60, 20,
            60, 0,
            0, 0
        });
    }

    public static PathD getRL()
    {
        return Clipper.MakePath(new double[] {
            0, 0,
            0, 20,
            10, 20,
            10, 50,
            60, 50,
            60, 0,
            0, 0
        });
    }

    public static PathD getU()
    {
        return Clipper.MakePath(new double[] {
            0, 0,
            0, 50,
            10, 50,
            10, 20,
            60, 20,
            60, 80,
            120, 80,
            120, 0,
            0, 0
        });
    }

    public static PathD getT()
    {
        return Clipper.MakePath(new double[] {
            0, 50,
            0, 80,
            80, 80,
            80, 50,
            60, 50,
            60, 0,
            40, 0,
            40, 50,
            0, 50
        });
    }

    public static PathD getX()
    {
        return Clipper.MakePath(new double[] {
            0, 50,
            0, 80,
            60, 80,
            60, 100,
            80, 100,
            80, 80,
            100, 80,
            100, 50,
            80, 50,
            80, 20,
            60, 20,
            60, 50,
            0, 50
        });
    }

    public static PathD getS()
    {
        return Clipper.MakePath(new double[] {
            0, 0,
            0, 20,
            20, 20,
            20, 50,
            0, 50,
            0, 110,
            100, 110,
            100, 80,
            80, 80,
            80, 60,
            100, 60,
            100, 0,
            0, 0
        });
    }

    public static PathD getnegS()
    {
        PathD S = Clipper.MakePath(new double[] {
            0, 0,
            0, 20,
            20, 20,
            20, 50,
            0, 50,
            0, 110,
            100, 110,
            100, 80,
            80, 80,
            80, 60,
            100, 60,
            100, 0,
            0, 0
        });

        for (int i = 0; i < S.Count; i++)
        {
            S[i] = new(S[i].x, S[i].y - 200);
        }

        return S;
    }

}