using System;
using System.Globalization;
using System.IO;
using System.IO.Compression;
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

    public static double myPow(double num, int exp)
    {
        switch (exp)
        {
            // Micro-optimization
            case 2:
                return num * num;
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

        while (exp > 0)
        {
            switch (exp % 2)
            {
                case 1:
                    result *= num;
                    break;
            }
            exp >>= 1;
            num *= num;
        }

        result = invertResult switch
        {
            true => 1.0f / result,
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
        GZipStream zipStream = new(compressedStream, CompressionMode.Compress);
        zipStream.Write(data, 0, data.Length);
        zipStream.Close();
        return compressedStream.ToArray();
    }

    public static byte[] decompress(byte[] data)
    {
        MemoryStream compressedStream = new(data);
        GZipStream zipStream = new(compressedStream, CompressionMode.Decompress);
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
        return Convert.ToBase64String(hashAlgorithm.Hash);
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