using System;
using System.Globalization;
using System.IO;
using System.IO.Compression;
using System.Runtime.CompilerServices;
using System.Runtime.Serialization;
using System.Security.Cryptography;

namespace utility;

public static class Utils
{
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

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double myPow(double num, int exp)
    {
        return exp switch
        {
            // Micro-optimizations for common cases
            0 => 1.0,
            1 => num,
            2 => num * num,
            3 => num * num * num,
            4 => num * num * num * num,
            -1 => 1.0 / num,
            -2 => 1.0 / (num * num),
            _ => FastPowGeneric(num, exp)
        };
    }

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static double FastPowGeneric(double num, int exp)
    {
        double result = 1.0;
        bool invertResult = false;

        if (exp < 0)
        {
            invertResult = true;
            exp = -exp; // More efficient than Math.Abs for int
        }

        // Fast exponentiation by squaring
        while (exp > 0)
        {
            if ((exp & 1) == 1) // Check if odd using bitwise AND
            {
                result *= num;
            }
            exp >>= 1;
            num *= num;
        }

        return invertResult ? 1.0 / result : result;
    }

    /// <summary>
    /// Optimized degrees to radians conversion with precomputed constant
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double toRadians(double deg)
    {
        return deg * 0.017453292519943295; // Math.PI / 180.0 precomputed
    }

    /// <summary>
    /// Optimized radians to degrees conversion with precomputed constant
    /// </summary>
    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double toDegrees(double rad)
    {
        return rad * 57.29577951308232; // 180.0 / Math.PI precomputed
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
        hashAlgorithm.ComputeHash(memoryStream.ToArray());
        return hashAlgorithm.Hash != null ? Convert.ToBase64String(hashAlgorithm.Hash) : "";
    }

    public static string GetMD5Hash(object input)
    {
        MD5 hasher = MD5.Create();
        return GetHash(hasher, input);
    }
        
    public static string GetSHA1Hash(object input)
    {
        SHA1 hasher = SHA1.Create();
        return GetHash(hasher, input);
    }
        
    public static string GetSHA256Hash(object input)
    {
        SHA256 hasher = SHA256.Create();
        return GetHash(hasher, input);
    }
}