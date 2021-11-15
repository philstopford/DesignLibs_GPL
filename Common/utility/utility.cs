using System;
using System.Globalization;
using System.IO;
using System.IO.Compression;
using System.Runtime.Serialization;
using System.Security.Cryptography;
using System.Text;

namespace utility
{
    public static partial class Utils
    {
        public static string friendlyNumber(Int32 number)
        {
            return pFriendlyNumber(number);
        }

        public static string friendlyNumber(Int64 number)
        {
            return pFriendlyNumber(number);
        }

        public static string friendlyNumber(double number)
        {
            return pFriendlyNumber(number);
        }

        static string pFriendlyNumber(double number)
        {
            string ret;
            if (number >= 1E12)
            {
                ret = (number / 1E12).ToString("0.##") + " trillion";
            }
            else if (number >= 1E9)
            {
                ret = (number / 1E9).ToString("0.##") + " billion";
            }
            else if (number >= 1E6)
            {
                ret = (number / 1E6).ToString("0.##") + " million";
            }
            else if (number >= 1E3)
            {
                ret = (number / 1E3).ToString("0.##") + " thousand";
            }
            else if (number >= 1E2)
            {
                ret = (number / 1E2).ToString("0.##") + " hundred";
            }
            else
            {
                ret = number.ToString(CultureInfo.InvariantCulture);
            }

            return ret;
        }

        public static double myPow(double num, int exp)
        {
            // Micro-optimization
            if (exp == 2)
            {
                return num * num;
            }

            double result = 1.0;
            bool invertResult = false;

            if (exp < 0)
            {
                invertResult = true;
                exp = Math.Abs(exp);
            }

            while (exp > 0)
            {
                if (exp % 2 == 1)
                    result *= num;
                exp >>= 1;
                num *= num;
            }

            if (invertResult)
            {
                result = 1.0f / result;
            }

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
            MemoryStream compressedStream = new MemoryStream();
            GZipStream zipStream = new GZipStream(compressedStream, CompressionMode.Compress);
            zipStream.Write(data, 0, data.Length);
            zipStream.Close();
            return compressedStream.ToArray();
        }

        public static byte[] decompress(byte[] data)
        {
            MemoryStream compressedStream = new MemoryStream(data);
            GZipStream zipStream = new GZipStream(compressedStream, CompressionMode.Decompress);
            MemoryStream resultStream = new MemoryStream();
            zipStream.CopyTo(resultStream);
            return resultStream.ToArray();
        }

        private static string GetHash(HashAlgorithm hashAlgorithm, object input)
        {
            using (MemoryStream memoryStream = new MemoryStream())
            {
                DataContractSerializer serializer = new DataContractSerializer(input.GetType());
                serializer.WriteObject(memoryStream, input);
                hashAlgorithm.ComputeHash(memoryStream.ToArray());
                return Convert.ToBase64String(hashAlgorithm.Hash);
            }
        }

        public static string GetMD5Hash(object input)
        {
            MD5 hasher = MD5.Create();
            string hash = GetHash(hasher, input);
            return hash;
        }
        
        public static string GetMD5Hash(object input)
        {
            MD5 hasher = MD5.Create();
            string hash = GetHash(hasher, input);
            return hash;
        }
    }
    
}
