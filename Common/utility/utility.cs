using System;
using System.Collections.Concurrent;
using System.Globalization;
using System.IO;
using System.IO.Compression;
using System.Runtime.CompilerServices;
using System.Runtime.Serialization;
using System.Security.Cryptography;
using System.Text;
using System.Threading;

namespace utility;

public static class Utils
{
    // Performance optimization: Cache for commonly used hash algorithm instances
    private static readonly ThreadLocal<MD5> _threadLocalMD5 = new(() => MD5.Create());
    private static readonly ThreadLocal<SHA1> _threadLocalSHA1 = new(() => SHA1.Create());
    private static readonly ThreadLocal<SHA256> _threadLocalSHA256 = new(() => SHA256.Create());
    
    // Performance optimization: Cache for serializers to avoid repeated creation
    private static readonly ConcurrentDictionary<Type, DataContractSerializer> _serializerCache = new();
    
    // Performance optimization: Lookup table for small integer powers
    private static readonly double[] PowersOfTwo = new double[32];
    private static readonly double[] PowersOfTen = new double[16];
    private static readonly double[] PowersOfE = new double[16]; // e^x for small x
    private static readonly double[] PowersOfHalf = new double[16]; // 0.5^x for small x
    
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
        for (int i = 0; i < PowersOfE.Length; i++)
        {
            PowersOfE[i] = Math.Pow(Math.E, i);
        }
        for (int i = 0; i < PowersOfHalf.Length; i++)
        {
            PowersOfHalf[i] = Math.Pow(0.5, i);
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

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    public static double myPow(double num, int exp)
    {
        // Performance optimization: Handle common cases with lookup tables
        if (exp >= 0)
        {
            if (num == 2.0 && exp < PowersOfTwo.Length)
            {
                return PowersOfTwo[exp];
            }
            
            if (num == 10.0 && exp < PowersOfTen.Length)
            {
                return PowersOfTen[exp];
            }
            
            if (Math.Abs(num - Math.E) < 1e-15 && exp < PowersOfE.Length)
            {
                return PowersOfE[exp];
            }
            
            if (num == 0.5 && exp < PowersOfHalf.Length)
            {
                return PowersOfHalf[exp];
            }
        }

        // Performance optimization: Handle more common cases first
        switch (exp)
        {
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
            case 5:
                var num2 = num * num;
                return num2 * num2 * num;
            case 6:
                var num3 = num * num * num;
                return num3 * num3;
            case 8:
                var num4 = num * num;
                num4 *= num4;
                return num4 * num4;
        }

        double result = 1.0;
        bool invertResult = false;

        if (exp < 0)
        {
            invertResult = true;
            exp = -exp; // Use subtraction instead of Math.Abs for better performance
        }

        // Performance optimization: Use binary exponentiation (exponentiation by squaring)
        double baseValue = num;
        while (exp > 0)
        {
            if ((exp & 1) == 1)
            {
                result *= baseValue;
            }
            exp >>= 1;
            baseValue *= baseValue;
        }

        return invertResult ? 1.0 / result : result;
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

    [MethodImpl(MethodImplOptions.AggressiveInlining)]
    private static string GetHash(HashAlgorithm hashAlgorithm, object input)
    {
        byte[] inputBytes;
        
        // Performance optimization: Fast paths for common types to avoid serialization
        switch (input)
        {
            case string str:
                inputBytes = Encoding.UTF8.GetBytes(str);
                break;
            case int intVal:
                inputBytes = BitConverter.GetBytes(intVal);
                break;
            case long longVal:
                inputBytes = BitConverter.GetBytes(longVal);
                break;
            case double doubleVal:
                inputBytes = BitConverter.GetBytes(doubleVal);
                break;
            case float floatVal:
                inputBytes = BitConverter.GetBytes(floatVal);
                break;
            case byte[] byteArray:
                inputBytes = byteArray;
                break;
            default:
                {
                    // Fall back to serialization for complex types with cached serializer
                    var serializer = _serializerCache.GetOrAdd(input.GetType(), type => new DataContractSerializer(type));
                    using var memoryStream = new MemoryStream();
                    serializer.WriteObject(memoryStream, input);
                    inputBytes = memoryStream.ToArray();
                    break;
                }
        }
        
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