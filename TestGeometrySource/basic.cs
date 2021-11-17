using System;
using geoLib;

namespace PartitionTestGeometrySource;

public static partial class TestGeometry
{
    public static GeoLibPoint[] getL()
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

        return L;
    }

    public static GeoLibPoint[] getRL()
    {
        GeoLibPoint[] rL = new GeoLibPoint[] {

            new GeoLibPoint( 0, 0),
            new GeoLibPoint( 0, 20),
            new GeoLibPoint( 10, 20),
            new GeoLibPoint( 10, 50),
            new GeoLibPoint( 60, 50),
            new GeoLibPoint( 60, 0),
            new GeoLibPoint( 0, 0)

        };

        return rL;
    }

    public static GeoLibPoint[] getU()
    {
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

        return U;
    }

    public static GeoLibPoint[] getT()
    {
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

        return T;
    }

    public static GeoLibPoint[] getX()
    {
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

        return X;
    }

    public static GeoLibPoint[] getS()
    {
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

        return S;
    }

    public static GeoLibPoint[] getnegS()
    {
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

        for (int i = 0; i < S.Length; i++)
        {
            S[i] = new GeoLibPoint(S[i].X, S[i].Y - 200);
        }

        return S;
    }

}