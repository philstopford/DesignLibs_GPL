using System;
using System.Globalization;
using System.IO;
using System.IO.Compression;
using System.Runtime.Serialization;
using System.Security.Cryptography;
using System.Threading;

namespace utility;

public static class Utils
{
    // Performance optimization: Cache for commonly used hash algorithm instances
    private static readonly ThreadLocal<MD5> _threadLocalMD5 = new(() => MD5.Create());
    private static readonly ThreadLocal<SHA1> _threadLocalSHA1 = new(() => SHA1.Create());
    private static readonly ThreadLocal<SHA256> _threadLocalSHA256 = new(() => SHA256.Create());
    
    // Performance optimization: Lookup table for small integer powers
    private static readonly double[] PowersOfTwo = new double[32];
    private static readonly double[] PowersOfTen = new double[16];
    
    static Utils()
    {
        // Pre-calculate common powers for fast lookup
        for (int i = 0; i < PowersOfTwo.Length; i++)
        {
            PowersOfTwo[i] = Math.Pow(2.0, i);
        }
        for (int i = 0; i < PowersOfTen.Length; i++)
        {
            PowersOfTen[i] = Math.Pow(10.0, i);
        }
    }
    public static string friendlyNumber(int number)
    {
        return pFriendlyNumber(number);
    }

    public static string friendlyNumber(long number)
    {
        return pFriendlyNumber(number);
    }

    public static string friendlyNumber(double number)
    {
        return pFriendlyNumber(number);
    }

    private static string pFriendlyNumber(double number)
    {
        string ret = number switch
        {
            >= 1E12 => (number / 1E12).ToString("0.##") + " trillion",
            >= 1E9 => (number / 1E9).ToString("0.##") + " billion",
            >= 1E6 => (number / 1E6).ToString("0.##") + " million",
            >= 1E3 => (number / 1E3).ToString("0.##") + " thousand",
            >= 1E2 => (number / 1E2).ToString("0.##") + " hundred",
            _ => number.ToString(CultureInfo.InvariantCulture)
        };

        return ret;
    }

    public static double myPow(double num, int exp)
    {
        // Performance optimization: Handle common cases with lookup tables
        if (num == 2.0 && exp >= 0 && exp < PowersOfTwo.Length)
        {
            return PowersOfTwo[exp];
        }
        
        if (num == 10.0 && exp >= 0 && exp < PowersOfTen.Length)
        {
            return PowersOfTen[exp];
        }

        switch (exp)
        {
            // Performance optimization: Handle more common cases
            case 0:
                return 1.0;
            case 1:
                return num;
            case 2:
                return num * num;
            case 3:
                return num * num * num;
            case 4:
                var numSquared = num * num;
                return numSquared * numSquared;
        }

        double result = 1.0;
        bool invertResult = false;

        switch (exp)
        {
            case < 0:
                invertResult = true;
                exp = Math.Abs(exp);
                break;
        }

        // Performance optimization: Use binary exponentiation (exponentiation by squaring)
        while (exp > 0)
        {
            if ((exp & 1) == 1)
            {
                result *= num;
            }
            exp >>= 1;
            num *= num;
        }

        result = invertResult switch
        {
            true => 1.0 / result,
            _ => result
        };

        return result;
    }

    public static double toRadians(double deg)
    {
        return deg * (Math.PI / 180.0f);
    }

    public static double toDegrees(double rad)
    {
        return rad / (Math.PI / 180.0f);
    }

    public static byte[] compress(byte[] data)
    {
        MemoryStream compressedStream = new();
        ZLibStream zipStream = new(compressedStream, CompressionMode.Compress);
        zipStream.Write(data, 0, data.Length);
        zipStream.Close();
        return compressedStream.ToArray();
    }

    public static byte[] decompress(byte[] data)
    {
        MemoryStream compressedStream = new(data);
        DeflateStream zipStream = new(compressedStream, CompressionMode.Decompress);
        MemoryStream resultStream = new();
        zipStream.CopyTo(resultStream);
        return resultStream.ToArray();
    }

    private static string GetHash(HashAlgorithm hashAlgorithm, object input)
    {
        using MemoryStream memoryStream = new();
        DataContractSerializer serializer = new(input.GetType());
        serializer.WriteObject(memoryStream, input);
        byte[] inputBytes = memoryStream.ToArray();
        byte[] hashBytes = hashAlgorithm.ComputeHash(inputBytes);
        return Convert.ToBase64String(hashBytes);
    }

    public static string GetMD5Hash(object input)
    {
        var hasher = _threadLocalMD5.Value!;
        return GetHash(hasher, input);
    }
        
    public static string GetSHA1Hash(object input)
    {
        var hasher = _threadLocalSHA1.Value!;
        return GetHash(hasher, input);
    }
        
    public static string GetSHA256Hash(object input)
    {
        var hasher = _threadLocalSHA256.Value!;
        return GetHash(hasher, input);
    }
}