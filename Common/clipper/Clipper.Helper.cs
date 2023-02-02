namespace Clipper2Lib;

public static class Helper
{
    public static Path64 initedPath64(int count)
    {
        Path64 ret = new(count);
        for (int i = 0; i < count; i++)
        {
            ret.Add(new());
        }

        return ret;
    }
    
    public static Paths64 initedPaths64(int count)
    {
        Paths64 ret = new(count);
        for (int i = 0; i < count; i++)
        {
            ret.Add(new());
        }

        return ret;
    }
    
    public static PathD initedPathD(int count)
    {
        PathD ret = new(count);
        for (int i = 0; i < count; i++)
        {
            ret.Add(new());
        }

        return ret;
    }
    
    public static PathsD initedPathsD(int count)
    {
        PathsD ret = new(count);
        for (int i = 0; i < count; i++)
        {
            ret.Add(new());
        }

        return ret;
    }
}