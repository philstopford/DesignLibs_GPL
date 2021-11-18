using geoLib;

namespace PartitionTestGeometrySource;

public static partial class TestGeometry
{
    public static GeoLibPoint[] getL()
    {
        // L
        GeoLibPoint[] L = new GeoLibPoint[] {

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

    public static GeoLibPoint[] getRL()
    {
        GeoLibPoint[] rL = new GeoLibPoint[] {

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

    public static GeoLibPoint[] getU()
    {
        GeoLibPoint[] U = new GeoLibPoint[] {

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

    public static GeoLibPoint[] getT()
    {
        GeoLibPoint[] T = new GeoLibPoint[] {

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

    public static GeoLibPoint[] getX()
    {
        GeoLibPoint[] X = new GeoLibPoint[] {

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

    public static GeoLibPoint[] getS()
    {
        GeoLibPoint[] S = new GeoLibPoint[] {

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

    public static GeoLibPoint[] getnegS()
    {
        GeoLibPoint[] S = new GeoLibPoint[] {

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

        for (int i = 0; i < S.Length; i++)
        {
            S[i] = new GeoLibPoint(S[i].X, S[i].Y - 200);
        }

        return S;
    }

}