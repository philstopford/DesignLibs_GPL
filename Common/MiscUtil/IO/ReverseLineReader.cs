﻿using System;
using System.Collections;
using System.Collections.Generic;
using System.IO;
using System.Text;

namespace MiscUtil.IO;

/// <summary>
/// Takes an encoding (defaulting to UTF-8) and a function which produces a seekable stream
/// (or a filename for convenience) and yields lines from the end of the stream backwards.
/// Only single byte encodings, and UTF-8 and Unicode, are supported. The stream
/// returned by the function must be seekable.
/// </summary>
public sealed class ReverseLineReader : IEnumerable<string>
{
    /// <summary>
    /// Buffer size to use by default. Classes with internal access can specify
    /// a different buffer size - this is useful for testing.
    /// </summary>
    private const int DefaultBufferSize = 4096;

    /// <summary>
    /// Means of creating a Stream to read from.
    /// </summary>
    private readonly Func<Stream> streamSource;

    /// <summary>
    /// Encoding to use when converting bytes to text
    /// </summary>
    private readonly Encoding encoding;

    /// <summary>
    /// Size of buffer (in bytes) to read each time we read from the
    /// stream. This must be at least as big as the maximum number of
    /// bytes for a single character.
    /// </summary>
    private readonly int bufferSize;

    /// <summary>
    /// Function which, when given a position within a file and a byte, states whether
    /// or not the byte represents the start of a character.
    /// </summary>
    private readonly Func<long, byte, bool> characterStartDetector;

    /// <summary>
    /// Creates a LineReader from a stream source. The delegate is only
    /// called when the enumerator is fetched. UTF-8 is used to decode
    /// the stream into text.
    /// </summary>
    /// <param name="streamSource">Data source</param>
    public ReverseLineReader(Func<Stream> streamSource)
        : this(streamSource, Encoding.UTF8)
    {
    }

    /// <summary>
    /// Creates a LineReader from a filename. The file is only opened
    /// (or even checked for existence) when the enumerator is fetched.
    /// UTF8 is used to decode the file into text.
    /// </summary>
    /// <param name="filename">File to read from</param>
    public ReverseLineReader(string filename)
        : this(filename, Encoding.UTF8)
    {
    }

    /// <summary>
    /// Creates a LineReader from a filename. The file is only opened
    /// (or even checked for existence) when the enumerator is fetched.
    /// </summary>
    /// <param name="filename">File to read from</param>
    /// <param name="encoding">Encoding to use to decode the file into text</param>
    public ReverseLineReader(string filename, Encoding encoding)
        : this(() => File.OpenRead(filename), encoding)
    {
    }

    /// <summary>
    /// Creates a LineReader from a stream source. The delegate is only
    /// called when the enumerator is fetched.
    /// </summary>
    /// <param name="streamSource">Data source</param>
    /// <param name="encoding">Encoding to use to decode the stream into text</param>
    public ReverseLineReader(Func<Stream> streamSource, Encoding encoding)
        : this(streamSource, encoding, DefaultBufferSize)
    {
    }

    internal ReverseLineReader(Func<Stream> streamSource, Encoding encoding, int bufferSize)
    {
        this.streamSource = streamSource;
        this.encoding = encoding;
        this.bufferSize = bufferSize;
        characterStartDetector = encoding.IsSingleByte switch
        {
            true =>
                // For a single byte encoding, every byte is the start (and end) of a character
                (pos, data) => true,
            _ => encoding switch
            {
                UnicodeEncoding =>
                    // For UTF-16, even-numbered positions are the start of a character
                    (pos, data) => (pos & 1) == 0,
                UTF8Encoding =>
                    // For UTF-8, bytes with the top bit clear or the second bit set are the start of a character
                    // See http://www.cl.cam.ac.uk/~mgk25/unicode.html
                    (pos, data) => (data & 0x80) == 0 || (data & 0x40) != 0,
                _ => throw new ArgumentException("Only single byte, UTF-8 and Unicode encodings are permitted")
            }
        };
    }

    /// <summary>
    /// Returns the enumerator reading strings backwards. If this method discovers that
    /// the returned stream is either unreadable or unseekable, a NotSupportedException is thrown.
    /// </summary>
    public IEnumerator<string> GetEnumerator()
    {
        Stream stream = streamSource();
        switch (stream.CanSeek)
        {
            case false:
                stream.Dispose();
                throw new NotSupportedException("Unable to seek within stream");
        }
        switch (stream.CanRead)
        {
            case false:
                stream.Dispose();
                throw new NotSupportedException("Unable to read within stream");
            default:
                return GetEnumeratorImpl(stream);
        }
    }

    private IEnumerator<string> GetEnumeratorImpl(Stream stream)
    {
        using (stream)
        {
            long position = stream.Length;

            switch (encoding)
            {
                case UnicodeEncoding when (position & 1) != 0:
                    throw new InvalidDataException("UTF-16 encoding provided, but stream has odd length.");
            }

            // Allow up to two bytes for data from the start of the previous
            // read which didn't quite make it as full characters
            byte[] buffer = new byte[bufferSize + 2];
            char[] charBuffer = new char[encoding.GetMaxCharCount(buffer.Length)];
            int leftOverData = 0;
            string previousEnd = null;
            // TextReader doesn't return an empty string if there's line break at the end
            // of the data. Therefore we don't return an empty string if it's our *first*
            // return.
            bool firstYield = true;

            // A line-feed at the start of the previous buffer means we need to swallow
            // the carriage-return at the end of this buffer - hence this needs declaring
            // way up here!
            bool swallowCarriageReturn = false;

            while (position > 0)
            {
                int bytesToRead = Math.Min(position > int.MaxValue ? bufferSize : (int)position, bufferSize);

                position -= bytesToRead;
                stream.Position = position;
                StreamUtil.ReadExactly(stream, buffer, bytesToRead);
                switch (leftOverData)
                {
                    // If we haven't read a full buffer, but we had bytes left
                    // over from before, copy them to the end of the buffer
                    case > 0 when bytesToRead != bufferSize:
                        // Buffer.BlockCopy doesn't document its behaviour with respect
                        // to overlapping data: we *might* just have read 7 bytes instead of
                        // 8, and have two bytes to copy...
                        Array.Copy(buffer, bufferSize, buffer, bytesToRead, leftOverData);
                        break;
                }
                // We've now *effectively* read this much data.
                bytesToRead += leftOverData;

                int firstCharPosition = 0;
                while (!characterStartDetector(position + firstCharPosition, buffer[firstCharPosition]))
                {
                    firstCharPosition++;
                    // Bad UTF-8 sequences could trigger this. For UTF-8 we should always
                    // see a valid character start in every 3 bytes, and if this is the start of the file
                    // so we've done a short read, we should have the character start
                    // somewhere in the usable buffer.
                    if (firstCharPosition == 3 || firstCharPosition == bytesToRead)
                    {
                        throw new InvalidDataException("Invalid UTF-8 data");
                    }
                }
                leftOverData = firstCharPosition;

                int charsRead = encoding.GetChars(buffer, firstCharPosition, bytesToRead - firstCharPosition, charBuffer, 0);
                int endExclusive = charsRead;

                for (int i = charsRead - 1; i >= 0; i--)
                {
                    char lookingAt = charBuffer[i];
                    switch (swallowCarriageReturn)
                    {
                        case true:
                        {
                            swallowCarriageReturn = false;
                            switch (lookingAt)
                            {
                                case '\r':
                                    endExclusive--;
                                    continue;
                            }

                            break;
                        }
                    }
                    // Anything non-line-breaking, just keep looking backwards
                    if (lookingAt != '\n' && lookingAt != '\r')
                    {
                        continue;
                    }

                    swallowCarriageReturn = lookingAt switch
                    {
                        // End of CRLF? Swallow the preceding CR
                        '\n' => true,
                        _ => swallowCarriageReturn
                    };
                    int start = i + 1;
                    string bufferContents = new(charBuffer, start, endExclusive - start);
                    endExclusive = i;
                    string stringToYield = previousEnd == null ? bufferContents : bufferContents + previousEnd; ;
                    if (!firstYield || stringToYield.Length != 0)
                    {
                        yield return stringToYield;
                    }
                    firstYield = false;
                    previousEnd = null;
                }

                previousEnd = endExclusive == 0 ? null : new string(charBuffer, 0, endExclusive) + previousEnd;

                // If we didn't decode the start of the array, put it at the end for next time
                if (leftOverData != 0)
                {
                    Buffer.BlockCopy(buffer, 0, buffer, bufferSize, leftOverData);
                }
            }
            if (leftOverData != 0)
            {
                // At the start of the final buffer, we had the end of another character.
                throw new InvalidDataException("Invalid UTF-8 data at start of stream");
            }
            switch (firstYield)
            {
                case true when string.IsNullOrEmpty(previousEnd):
                    yield break;
                default:
                    yield return previousEnd ?? "";
                    break;
            }
        }
    }

    IEnumerator IEnumerable.GetEnumerator()
    {
        return GetEnumerator();
    }
}