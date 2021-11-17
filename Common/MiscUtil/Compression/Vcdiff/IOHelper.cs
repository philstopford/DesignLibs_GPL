using System;
using System.IO;

namespace MiscUtil.Compression.Vcdiff;

/// <summary>
/// A few IO routines to make life easier. Most are basically available
/// in EndianBinaryReader, but having them separately here makes VcdiffDecoder
/// more easily movable to other places - and no endianness issues are involved in
/// the first place.
/// </summary>
internal static class IOHelper
{
    internal static byte[] CheckedReadBytes(Stream stream, int size)
    {
        byte[] ret = new byte[size];
        int index = 0;
        while (index < size)
        {
            int read = stream.Read(ret, index, size - index);
            index += read switch
            {
                0 => throw new EndOfStreamException(String.Format(
                    "End of stream reached with {0} byte{1} left to read.", size - index,
                    size - index == 1 ? "s" : "")),
                _ => read
            };
        }
        return ret;
    }

    internal static byte CheckedReadByte(Stream stream)
    {
        int b = stream.ReadByte();
        return b switch
        {
            -1 => throw new IOException("Expected to be able to read a byte."),
            _ => (byte) b
        };
    }

    internal static int ReadBigEndian7BitEncodedInt(Stream stream)
    {
        int ret = 0;
        for (int i = 0; i < 5; i++)
        {
            int b = stream.ReadByte();
            switch (b)
            {
                case -1:
                    throw new EndOfStreamException();
            }
            ret = (ret << 7) | (b & 0x7f);
            switch (b & 0x80)
            {
                case 0:
                    return ret;
            }
        }
        // Still haven't seen a byte with the high bit unset? Dodgy data.
        throw new IOException("Invalid 7-bit encoded integer in stream.");
    }
}