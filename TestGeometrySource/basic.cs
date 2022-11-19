using System.IO;
using Clipper2Lib;
using geoLib;

namespace PartitionTestGeometrySource;

public static partial class TestGeometry
{
    public static Path64 getL()
    {
        // L
        Path64 L = new () {

            new( 0, 0),
            new( 0, 50),
            new( 10, 50),
            new( 10, 20),
            new( 60, 20),
            new( 60, 0),
            new( 0, 0)

        };

        return L;
    }

    public static Path64 getRL()
    {
        Path64 rL = new () {

            new( 0, 0),
            new( 0, 20),
            new( 10, 20),
            new( 10, 50),
            new( 60, 50),
            new( 60, 0),
            new( 0, 0)

        };

        return rL;
    }

    public static Path64 getU()
    {
        Path64 U = new () {

            new( 0, 0),
            new( 0, 50),
            new( 10, 50),
            new( 10, 20),
            new( 60, 20),
            new( 60, 80),
            new( 120, 80),
            new( 120, 0),
            new( 0, 0)

        };

        return U;
    }

    public static Path64 getT()
    {
        Path64 T = new () {

            new( 0, 50),
            new( 0, 80),
            new( 80, 80),
            new( 80, 50),
            new( 60, 50),
            new( 60, 0),
            new( 40, 0),
            new( 40, 50),
            new( 0, 50)

        };

        return T;
    }

    public static Path64 getX()
    {
        Path64 X = new () {

            new( 0, 50),
            new( 0, 80),
            new( 60, 80),
            new( 60, 100),
            new( 80, 100),
            new( 80, 80),
            new( 100, 80),
            new( 100, 50),
            new( 80, 50),
            new( 80, 20),
            new( 60, 20),
            new( 60, 50),
            new( 0, 50)

        };

        return X;
    }

    public static Path64 getS()
    {
        Path64 S = new () {

            new( 0, 0),
            new( 0, 20),
            new( 20, 20),
            new( 20, 50),
            new( 0, 50),
            new( 0, 110),
            new( 100, 110),
            new( 100, 80),
            new( 80, 80),
            new( 80, 60),
            new( 100, 60),
            new( 100, 0),
            new( 0, 0)

        };

        return S;
    }

    public static Path64 getnegS()
    {
        Path64 S = new () {

            new( 0, 0),
            new( 0, 20),
            new( 20, 20),
            new( 20, 50),
            new( 0, 50),
            new( 0, 110),
            new( 100, 110),
            new( 100, 80),
            new( 80, 80),
            new( 80, 60),
            new( 100, 60),
            new( 100, 0),
            new( 0, 0)

        };

        for (int i = 0; i < S.Count; i++)
        {
            S[i] = new (S[i].X, S[i].Y - 200);
        }

        return S;
    }

}