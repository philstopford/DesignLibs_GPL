// See https://aka.ms/new-console-template for more information

using ClipperLib1Test;
using ClipperLib2Test;

Console.WriteLine("Comparison tests for Clipper 1 and Clipper 2");

Clipper1Test.offsetTest();
Clipper2Test.offsetTest();

Clipper1Test.test1();
Clipper2Test.test1();
