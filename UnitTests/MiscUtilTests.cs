using MiscUtil;
using MiscUtil.Collections;
using NUnit.Framework;
using System;
using System.Collections.Generic;
using System.Linq;

namespace UnitTests;

/// <summary>
/// Tests for the MiscUtil library functionality, covering utility classes,
/// collections, and helper functions.
/// </summary>
[TestFixture]
public class MiscUtilTests
{
    #region StaticRandom Tests

    [Test]
    public void StaticRandom_Next_ShouldGenerateNonNegativeValues()
    {
        for (int i = 0; i < 10; i++)
        {
            int value = StaticRandom.Next();
            Assert.That(value, Is.GreaterThanOrEqualTo(0));
        }
    }

    [Test]
    public void StaticRandom_NextWithMax_ShouldRespectBounds()
    {
        int max = 100;
        for (int i = 0; i < 50; i++)
        {
            int value = StaticRandom.Next(max);
            Assert.That(value, Is.GreaterThanOrEqualTo(0).And.LessThan(max));
        }
    }

    [Test]
    public void StaticRandom_NextWithMinMax_ShouldRespectBounds()
    {
        int min = 10;
        int max = 20;
        for (int i = 0; i < 50; i++)
        {
            int value = StaticRandom.Next(min, max);
            Assert.That(value, Is.GreaterThanOrEqualTo(min).And.LessThan(max));
        }
    }

    [Test]
    public void StaticRandom_NextDouble_ShouldGenerateValuesInRange()
    {
        for (int i = 0; i < 50; i++)
        {
            double value = StaticRandom.NextDouble();
            Assert.That(value, Is.GreaterThanOrEqualTo(0.0).And.LessThan(1.0));
        }
    }

    #endregion

    #region Range Tests

    [Test]
    public void Range_Constructor_ShouldCreateInclusiveRange()
    {
        var range = new Range<int>(1, 10);

        Assert.That(range.Start, Is.EqualTo(1));
        Assert.That(range.End, Is.EqualTo(10));
        Assert.That(range.IncludesStart, Is.True);
        Assert.That(range.IncludesEnd, Is.True);
    }

    [Test]
    public void Range_Contains_ShouldCheckValueInRange()
    {
        var range = new Range<int>(1, 10);

        Assert.That(range.Contains(1), Is.True);  // Start value
        Assert.That(range.Contains(5), Is.True);  // Middle value
        Assert.That(range.Contains(10), Is.True); // End value
        Assert.That(range.Contains(0), Is.False); // Below range
        Assert.That(range.Contains(11), Is.False); // Above range
    }

    [Test]
    public void Range_ExcludeStart_ShouldCreateNewRangeWithoutStart()
    {
        var range = new Range<int>(1, 10);
        var newRange = range.ExcludeStart();

        Assert.That(newRange.IncludesStart, Is.False);
        Assert.That(newRange.IncludesEnd, Is.True);
        Assert.That(newRange.Contains(1), Is.False);
        Assert.That(newRange.Contains(2), Is.True);
    }

    [Test]
    public void Range_ExcludeEnd_ShouldCreateNewRangeWithoutEnd()
    {
        var range = new Range<int>(1, 10);
        var newRange = range.ExcludeEnd();

        Assert.That(newRange.IncludesStart, Is.True);
        Assert.That(newRange.IncludesEnd, Is.False);
        Assert.That(newRange.Contains(10), Is.False);
        Assert.That(newRange.Contains(9), Is.True);
    }

    #endregion

    #region SmartEnumerable Tests

    [Test]
    public void SmartEnumerable_Create_ShouldEnumerateWithIndexInformation()
    {
        var source = new[] { "A", "B", "C" };
        var smart = SmartEnumerable.Create(source);

        var entries = smart.ToList();

        Assert.That(entries.Count, Is.EqualTo(3));

        // Check first entry
        Assert.That(entries[0].Value, Is.EqualTo("A"));
        Assert.That(entries[0].IsFirst, Is.True);
        Assert.That(entries[0].IsLast, Is.False);
        Assert.That(entries[0].Index, Is.EqualTo(0));

        // Check middle entry
        Assert.That(entries[1].Value, Is.EqualTo("B"));
        Assert.That(entries[1].IsFirst, Is.False);
        Assert.That(entries[1].IsLast, Is.False);
        Assert.That(entries[1].Index, Is.EqualTo(1));

        // Check last entry
        Assert.That(entries[2].Value, Is.EqualTo("C"));
        Assert.That(entries[2].IsFirst, Is.False);
        Assert.That(entries[2].IsLast, Is.True);
        Assert.That(entries[2].Index, Is.EqualTo(2));
    }

    [Test]
    public void SmartEnumerable_WithSingleItem_ShouldMarkAsFirstAndLast()
    {
        var source = new[] { "OnlyItem" };
        var smart = SmartEnumerable.Create(source);

        var entry = smart.First();

        Assert.That(entry.Value, Is.EqualTo("OnlyItem"));
        Assert.That(entry.IsFirst, Is.True);
        Assert.That(entry.IsLast, Is.True);
        Assert.That(entry.Index, Is.EqualTo(0));
    }

    [Test]
    public void SmartEnumerable_WithEmptyCollection_ShouldProduceNoEntries()
    {
        var source = new string[0];
        var smart = SmartEnumerable.Create(source);

        Assert.That(smart.Count(), Is.EqualTo(0));
    }

    [Test]
    public void SmartEnumerable_Constructor_WithNullSource_ShouldThrowArgumentNullException()
    {
        Assert.Throws<ArgumentNullException>(() => new SmartEnumerable<string>(null));
    }

    #endregion

    #region RandomAccessQueue Tests

    [Test]
    public void RandomAccessQueue_Constructor_ShouldCreateEmptyQueue()
    {
        var queue = new RandomAccessQueue<int>();

        Assert.That(queue.Count, Is.EqualTo(0));
    }

    [Test]
    public void RandomAccessQueue_Enqueue_ShouldIncreaseCount()
    {
        var queue = new RandomAccessQueue<int>();

        queue.Enqueue(1);
        queue.Enqueue(2);

        Assert.That(queue.Count, Is.EqualTo(2));
    }

    [Test]
    public void RandomAccessQueue_Dequeue_ShouldDecreaseCount()
    {
        var queue = new RandomAccessQueue<int>();
        queue.Enqueue(1);
        queue.Enqueue(2);

        int removed = queue.Dequeue();

        Assert.That(removed, Is.EqualTo(1)); // FIFO behavior
        Assert.That(queue.Count, Is.EqualTo(1));
        Assert.That(queue[0], Is.EqualTo(2));
    }

    [Test]
    public void RandomAccessQueue_IndexAccess_ShouldGetCorrectValues()
    {
        var queue = new RandomAccessQueue<string>();
        queue.Enqueue("First");
        queue.Enqueue("Second");
        queue.Enqueue("Third");

        Assert.That(queue[0], Is.EqualTo("First"));
        Assert.That(queue[1], Is.EqualTo("Second"));
        Assert.That(queue[2], Is.EqualTo("Third"));
    }

    [Test]
    public void RandomAccessQueue_RemoveAt_ShouldRemoveSpecificIndex()
    {
        var queue = new RandomAccessQueue<int>();
        queue.Enqueue(10);
        queue.Enqueue(20);
        queue.Enqueue(30);

        int removed = queue.RemoveAt(1); // Remove middle element

        Assert.That(removed, Is.EqualTo(20));
        Assert.That(queue.Count, Is.EqualTo(2));
        Assert.That(queue[0], Is.EqualTo(10));
        Assert.That(queue[1], Is.EqualTo(30));
    }

    #endregion
}