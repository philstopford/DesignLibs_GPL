using MiscUtil.Conversion;
using System;
using System.IO;
using System.Text;

namespace MiscUtil.IO;

/// <summary>
/// Equivalent of System.IO.BinaryWriter, but with either endianness, depending on
/// the EndianBitConverter it is constructed with.
/// </summary>
public class EndianBinaryWriter : IDisposable
{
    #region Fields not directly related to properties
    /// <summary>
    /// Whether or not this writer has been disposed yet.
    /// </summary>
    private bool disposed;
    /// <summary>
    /// Buffer used for temporary storage during conversion from primitives
    /// </summary>
    private readonly byte[] buffer = new byte[16];
    /// <summary>
    /// Buffer used for Write(char)
    /// </summary>
    private readonly char[] charBuffer = new char[1];
    #endregion

    #region Constructors
    /// <summary>
    /// Constructs a new binary writer with the given bit converter, writing
    /// to the given stream, using UTF-8 encoding.
    /// </summary>
    /// <param name="bitConverter">Converter to use when writing data</param>
    /// <param name="stream">Stream to write data to</param>
    public EndianBinaryWriter(EndianBitConverter bitConverter,
        Stream stream) : this(bitConverter, stream, Encoding.UTF8)
    {
    }

    /// <summary>
    /// Constructs a new binary writer with the given bit converter, writing
    /// to the given stream, using the given encoding.
    /// </summary>
    /// <param name="bitConverter">Converter to use when writing data</param>
    /// <param name="stream">Stream to write data to</param>
    /// <param name="encoding">Encoding to use when writing character data</param>
    public EndianBinaryWriter(EndianBitConverter bitConverter, Stream stream, Encoding encoding)
    {
        switch (bitConverter)
        {
            case null:
                throw new ArgumentNullException(nameof(bitConverter));
        }
        switch (stream)
        {
            case null:
                throw new ArgumentNullException(nameof(stream));
        }
        switch (encoding)
        {
            case null:
                throw new ArgumentNullException(nameof(encoding));
        }
        switch (stream.CanWrite)
        {
            case false:
                throw new ArgumentException("Stream isn't writable", nameof(stream));
        }
        BaseStream = stream;
        BitConverter = bitConverter;
        Encoding = encoding;
    }
    #endregion

    #region Properties

    /// <summary>
    /// The bit converter used to write values to the stream
    /// </summary>
    public EndianBitConverter BitConverter { get; }

    /// <summary>
    /// The encoding used to write strings
    /// </summary>
    public Encoding Encoding { get; }

    /// <summary>
    /// Gets the underlying stream of the EndianBinaryWriter.
    /// </summary>
    public Stream BaseStream { get; }

    #endregion

    #region Public methods
    /// <summary>
    /// Closes the writer, including the underlying stream.
    /// </summary>
    public void Close()
    {
        Dispose();
    }

    /// <summary>
    /// Flushes the underlying stream.
    /// </summary>
    public void Flush()
    {
        CheckDisposed();
        BaseStream.Flush();
    }

    /// <summary>
    /// Seeks within the stream.
    /// </summary>
    /// <param name="offset">Offset to seek to.</param>
    /// <param name="origin">Origin of seek operation.</param>
    public void Seek(int offset, SeekOrigin origin)
    {
        CheckDisposed();
        BaseStream.Seek(offset, origin);
    }

    /// <summary>
    /// Writes a boolean value to the stream. 1 byte is written.
    /// </summary>
    /// <param name="value">The value to write</param>
    public void Write(bool value)
    {
        BitConverter.CopyBytes(value, buffer, 0);
        WriteInternal(buffer, 1);
    }

    /// <summary>
    /// Writes a 16-bit signed integer to the stream, using the bit converter
    /// for this writer. 2 bytes are written.
    /// </summary>
    /// <param name="value">The value to write</param>
    public void Write(short value)
    {
        BitConverter.CopyBytes(value, buffer, 0);
        WriteInternal(buffer, 2);
    }

    /// <summary>
    /// Writes a 32-bit signed integer to the stream, using the bit converter
    /// for this writer. 4 bytes are written.
    /// </summary>
    /// <param name="value">The value to write</param>
    public void Write(int value)
    {
        BitConverter.CopyBytes(value, buffer, 0);
        WriteInternal(buffer, 4);
    }

    /// <summary>
    /// Writes a 64-bit signed integer to the stream, using the bit converter
    /// for this writer. 8 bytes are written.
    /// </summary>
    /// <param name="value">The value to write</param>
    public void Write(long value)
    {
        BitConverter.CopyBytes(value, buffer, 0);
        WriteInternal(buffer, 8);
    }

    /// <summary>
    /// Writes a 16-bit unsigned integer to the stream, using the bit converter
    /// for this writer. 2 bytes are written.
    /// </summary>
    /// <param name="value">The value to write</param>
    public void Write(ushort value)
    {
        BitConverter.CopyBytes(value, buffer, 0);
        WriteInternal(buffer, 2);
    }

    /// <summary>
    /// Writes a 32-bit unsigned integer to the stream, using the bit converter
    /// for this writer. 4 bytes are written.
    /// </summary>
    /// <param name="value">The value to write</param>
    public void Write(uint value)
    {
        BitConverter.CopyBytes(value, buffer, 0);
        WriteInternal(buffer, 4);
    }

    /// <summary>
    /// Writes a 64-bit unsigned integer to the stream, using the bit converter
    /// for this writer. 8 bytes are written.
    /// </summary>
    /// <param name="value">The value to write</param>
    public void Write(ulong value)
    {
        BitConverter.CopyBytes(value, buffer, 0);
        WriteInternal(buffer, 8);
    }

    /// <summary>
    /// Writes a single-precision floating-point value to the stream, using the bit converter
    /// for this writer. 4 bytes are written.
    /// </summary>
    /// <param name="value">The value to write</param>
    public void Write(float value)
    {
        BitConverter.CopyBytes(value, buffer, 0);
        WriteInternal(buffer, 4);
    }

    /// <summary>
    /// Writes a double-precision floating-point value to the stream, using the bit converter
    /// for this writer. 8 bytes are written.
    /// </summary>
    /// <param name="value">The value to write</param>
    public void Write(double value)
    {
        BitConverter.CopyBytes(value, buffer, 0);
        WriteInternal(buffer, 8);
    }

    /// <summary>
    /// Writes a decimal value to the stream, using the bit converter for this writer.
    /// 16 bytes are written.
    /// </summary>
    /// <param name="value">The value to write</param>
    public void Write(decimal value)
    {
        BitConverter.CopyBytes(value, buffer, 0);
        WriteInternal(buffer, 16);
    }

    /// <summary>
    /// Writes a signed byte to the stream.
    /// </summary>
    /// <param name="value">The value to write</param>
    public void Write(byte value)
    {
        buffer[0] = value;
        WriteInternal(buffer, 1);
    }

    /// <summary>
    /// Writes an unsigned byte to the stream.
    /// </summary>
    /// <param name="value">The value to write</param>
    public void Write(sbyte value)
    {
        buffer[0] = unchecked((byte)value);
        WriteInternal(buffer, 1);
    }

    /// <summary>
    /// Writes an array of bytes to the stream.
    /// </summary>
    /// <param name="value">The values to write</param>
    public void Write(byte[] value)
    {
        switch (value)
        {
            case null:
                throw new ArgumentNullException(nameof(value));
            default:
                WriteInternal(value, value.Length);
                break;
        }
    }

    /// <summary>
    /// Writes a portion of an array of bytes to the stream.
    /// </summary>
    /// <param name="value">An array containing the bytes to write</param>
    /// <param name="offset">The index of the first byte to write within the array</param>
    /// <param name="count">The number of bytes to write</param>
    public void Write(byte[] value, int offset, int count)
    {
        CheckDisposed();
        BaseStream.Write(value, offset, count);
    }

    /// <summary>
    /// Writes a single character to the stream, using the encoding for this writer.
    /// </summary>
    /// <param name="value">The value to write</param>
    public void Write(char value)
    {
        charBuffer[0] = value;
        Write(charBuffer);
    }

    /// <summary>
    /// Writes an array of characters to the stream, using the encoding for this writer.
    /// </summary>
    /// <param name="value">An array containing the characters to write</param>
    public void Write(char[] value)
    {
        switch (value)
        {
            case null:
                throw new ArgumentNullException(nameof(value));
        }
        CheckDisposed();
        byte[] data = Encoding.GetBytes(value, 0, value.Length);
        WriteInternal(data, data.Length);
    }

    /// <summary>
    /// Writes a string to the stream, using the encoding for this writer.
    /// </summary>
    /// <param name="value">The value to write. Must not be null.</param>
    /// <exception cref="ArgumentNullException">value is null</exception>
    public void Write(string value)
    {
        switch (value)
        {
            case null:
                throw new ArgumentNullException(nameof(value));
        }
        CheckDisposed();
        byte[] data = Encoding.GetBytes(value);
        Write7BitEncodedInt(data.Length);
        WriteInternal(data, data.Length);
    }

    /// <summary>
    /// Writes a 7-bit encoded integer from the stream. This is stored with the least significant
    /// information first, with 7 bits of information per byte of value, and the top
    /// bit as a continuation flag.
    /// </summary>
    /// <param name="value">The 7-bit encoded integer to write to the stream</param>
    public void Write7BitEncodedInt(int value)
    {
        CheckDisposed();
        switch (value)
        {
            case < 0:
                throw new ArgumentOutOfRangeException(nameof(value), "Value must be greater than or equal to 0.");
        }
        int index = 0;
        while (value >= 128)
        {
            buffer[index++] = (byte)((value & 0x7f) | 0x80);
            value >>= 7;
            index++;
        }
        buffer[index++] = (byte)value;
        BaseStream.Write(buffer, 0, index);
    }

    #endregion

    #region Private methods
    /// <summary>
    /// Checks whether or not the writer has been disposed, throwing an exception if so.
    /// </summary>
    private void CheckDisposed()
    {
        switch (disposed)
        {
            case true:
                throw new ObjectDisposedException("EndianBinaryWriter");
        }
    }

    /// <summary>
    /// Writes the specified number of bytes from the start of the given byte array,
    /// after checking whether or not the writer has been disposed.
    /// </summary>
    /// <param name="bytes">The array of bytes to write from</param>
    /// <param name="length">The number of bytes to write</param>
    private void WriteInternal(byte[] bytes, int length)
    {
        CheckDisposed();
        BaseStream.Write(bytes, 0, length);
    }
    #endregion

    #region IDisposable Members
    /// <summary>
    /// Disposes of the underlying stream.
    /// </summary>
    public void Dispose()
    {
        switch (disposed)
        {
            case false:
                Flush();
                disposed = true;
                ((IDisposable)BaseStream).Dispose();
                break;
        }
    }
    #endregion
}