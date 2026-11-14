// Build against pinned upstream with USINGZ and -ffp-contract=off, then pass Tests/Polygons.txt
// or Tests/Lines.txt. Emits complete input/output fixtures, without rounding area comparisons.
#include "clipper2/clipper.h"
#include "ClipFileLoad.h"
#include <iostream>
using namespace Clipper2Lib;
static void paths(const Paths64& pp) {
  std::cout << '[';
  bool first = true;
  for (const auto& p : pp) {
    if (!first) std::cout << ','; first = false;
    std::cout << '[';
    bool firstPt = true;
    for (const auto& pt : p) {
      if (!firstPt) std::cout << ','; firstPt = false;
      std::cout << '[' << pt.x << ',' << pt.y << ']';
    }
    std::cout << ']';
  }
  std::cout << ']';
}
int main(int argc, char** argv) {
  if (argc != 2) return 1;
  std::ifstream input(argv[1]);
  if (!input) return 2;
  std::cout << '[';
  for (int n = 1;; ++n) {
    Paths64 subject, open, clip, closedResult, openResult;
    int64_t area, count;
    ClipType ct; FillRule fr;
    if (!LoadTestNum(input, n, subject, open, clip, area, count, ct, fr)) break;
    Clipper64 c;
    c.AddSubject(subject); c.AddOpenSubject(open); c.AddClip(clip);
    if (!c.Execute(ct, fr, closedResult, openResult)) return 3;
    if (n != 1) std::cout << ',';
    std::cout << "{\"number\":" << n << ",\"clipType\":" << int(ct) << ",\"fillRule\":" << int(fr);
    std::cout << ",\"subject\":"; paths(subject);
    std::cout << ",\"open\":"; paths(open);
    std::cout << ",\"clip\":"; paths(clip);
    std::cout << ",\"closedResult\":"; paths(closedResult);
    std::cout << ",\"openResult\":"; paths(openResult);
    std::cout << '}';
  }
  std::cout << "]\n";
}
