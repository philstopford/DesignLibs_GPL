import numpy
import gdstk

def testf():
    lib = gdstk.Library("Library")
    cell = lib.new_cell("Base")
    pol = gdstk.Polygon(((0, 0), (0, 1), (0.2, 1), (0.2, 0.8), (0.1, 0.8), (0.1, 0.6), (0.2, 0.6), (0.2, 0.4), (0.1, 0.4), (0.1, 0)))
    pol.repetition = gdstk.Repetition(x_offsets=(1, 3, -2))
    cell.add(pol)
    lib.write_gds("/d/Downloads/f_rep.gds")
    lib.write_oas("/d/Downloads/f_rep.oas")

testf()

def testf2():
    lib = gdstk.Library("Library")
    cell = lib.new_cell("Base")
    pol = gdstk.Polygon(((0, 0), (0, 1), (0.2, 1), (0.2, 0.8), (0.1, 0.8), (0.1, 0.6), (0.2, 0.6), (0.2, 0.4), (0.1, 0.4), (0.1, 0)))
    pol.repetition = gdstk.Repetition(columns=2, rows=1, v1=(3,3), v2=(4,4))
    cell.add(pol)
    lib.write_gds("/d/Downloads/f_rep2.gds")
    lib.write_oas("/d/Downloads/f_rep2.oas")

testf2()

def testf3():
    lib = gdstk.Library("Library")
    cell = lib.new_cell("Base")
    pol = gdstk.Polygon(((0, 0), (0, 1), (0.2, 1), (0.2, 0.8), (0.1, 0.8), (0.1, 0.6), (0.2, 0.6), (0.2, 0.4), (0.1, 0.4), (0.1, 0)))
    pol.repetition = gdstk.Repetition(columns=2, rows=3, spacing=(5,5))
    cell.add(pol)
    lib.write_gds("/d/Downloads/f_rep3.gds")
    lib.write_oas("/d/Downloads/f_rep3.oas")

testf3()

def testf4():
    lib = gdstk.Library("Library")
    cell = lib.new_cell("Base")
    pol = gdstk.Polygon(((0, 0), (0, 1), (0.2, 1), (0.2, 0.8), (0.1, 0.8), (0.1, 0.6), (0.2, 0.6), (0.2, 0.4), (0.1, 0.4), (0.1, 0)))
    pol.repetition = gdstk.Repetition(offsets=[(0.5, 1), (2, 0), (1.5, 0.5)])
    cell.add(pol)
    lib.write_gds("/d/Downloads/f_rep4.gds")
    lib.write_oas("/d/Downloads/f_rep4.oas")

testf4()

def testf5():
    lib = gdstk.Library("Library")
    cell = lib.new_cell("Base")
    pol = gdstk.Polygon(((0, 0), (0, 1), (0.2, 1), (0.2, 0.8), (0.1, 0.8), (0.1, 0.6), (0.2, 0.6), (0.2, 0.4), (0.1, 0.4), (0.1, 0)))
    pol.repetition = gdstk.Repetition(columns=1, rows=2, v1=(3,3), v2=(4,4))
    cell.add(pol)
    lib.write_gds("/d/Downloads/f_rep5.gds")
    lib.write_oas("/d/Downloads/f_rep5.oas")

testf5()

def test():
    lib = gdstk.Library("Library")
    cell = lib.new_cell("Base")
    pol = gdstk.rectangle((0, 0), (1, 1))
    pol.repetition = gdstk.Repetition(x_offsets=(1, 3, -2))
    cell.add(pol)
    lib.write_gds("/d/Downloads/rectangle_rep.gds")
    lib.write_oas("/d/Downloads/rectangle_rep.oas")

test()

def test2():
    lib = gdstk.Library("Library")
    cell = lib.new_cell("Base")
    pol = gdstk.rectangle((0, 0), (1, 1))
    pol.repetition = gdstk.Repetition(columns=2, rows=1, v1=(3,3), v2=(4,4))
    cell.add(pol)
    lib.write_gds("/d/Downloads/rectangle_rep2.gds")
    lib.write_oas("/d/Downloads/rectangle_rep2.oas")

test2()

def test3():
    lib = gdstk.Library("Library")
    cell = lib.new_cell("Base")
    pol = gdstk.rectangle((0, 0), (1, 1))
    pol.repetition = gdstk.Repetition(columns=2, rows=3, spacing=(5,5))
    cell.add(pol)
    lib.write_gds("/d/Downloads/rectangle_rep3.gds")
    lib.write_oas("/d/Downloads/rectangle_rep3.oas")

test3()

def test4():
    lib = gdstk.Library("Library")
    cell = lib.new_cell("Base")
    pol = gdstk.rectangle((0, 0), (1, 1))
    pol.repetition = gdstk.Repetition(offsets=[(0.5, 1), (2, 0), (1.5, 0.5)])
    cell.add(pol)
    lib.write_gds("/d/Downloads/rectangle_rep4.gds")
    lib.write_oas("/d/Downloads/rectangle_rep4.oas")

test4()

def test5():
    lib = gdstk.Library("Library")
    cell = lib.new_cell("Base")
    pol = gdstk.rectangle((0, 0), (1, 1))
    pol.repetition = gdstk.Repetition(columns=1, rows=2, v1=(3,3), v2=(4,4))
    cell.add(pol)
    lib.write_gds("/d/Downloads/rectangle_rep5.gds")
    lib.write_oas("/d/Downloads/rectangle_rep5.oas")

test5()

def ref_testf():
    lib = gdstk.Library("Library")
    ref_cell = lib.new_cell("Base")
    pol = gdstk.Polygon(((0, 0), (0, 1), (0.2, 1), (0.2, 0.8), (0.1, 0.8), (0.1, 0.6), (0.2, 0.6), (0.2, 0.4), (0.1, 0.4), (0.1, 0)))
    ref_cell.add(pol)
    cell = lib.new_cell("Ref")
    ref = gdstk.Reference(ref_cell,(0,0))
    ref.repetition = gdstk.Repetition(x_offsets=(1, 3, -2))
    cell.add(ref)
    lib.write_gds("/d/Downloads/ref_f_rep.gds")
    lib.write_oas("/d/Downloads/ref_f_rep.oas")

ref_testf()

def ref_testf2():
    lib = gdstk.Library("Library")
    ref_cell = lib.new_cell("Base")
    pol = gdstk.Polygon(((0, 0), (0, 1), (0.2, 1), (0.2, 0.8), (0.1, 0.8), (0.1, 0.6), (0.2, 0.6), (0.2, 0.4), (0.1, 0.4), (0.1, 0)))
    ref_cell.add(pol)
    cell = lib.new_cell("Ref")
    ref = gdstk.Reference(ref_cell,(0,0))
    ref.repetition = gdstk.Repetition(columns=2, rows=1, v1=(3,3), v2=(4,4))
    cell.add(ref)
    lib.write_gds("/d/Downloads/ref_f_rep2.gds")
    lib.write_oas("/d/Downloads/ref_f_rep2.oas")

ref_testf2()

def ref_testf3():
    lib = gdstk.Library("Library")
    ref_cell = lib.new_cell("Base")
    pol = gdstk.Polygon(((0, 0), (0, 1), (0.2, 1), (0.2, 0.8), (0.1, 0.8), (0.1, 0.6), (0.2, 0.6), (0.2, 0.4), (0.1, 0.4), (0.1, 0)))
    ref_cell.add(pol)
    cell = lib.new_cell("Ref")
    ref = gdstk.Reference(ref_cell,(0,0))
    ref.repetition = gdstk.Repetition(columns=2, rows=3, spacing=(5,5))
    cell.add(ref)
    lib.write_gds("/d/Downloads/ref_f_rep3.gds")
    lib.write_oas("/d/Downloads/ref_f_rep3.oas")

ref_testf3()

def ref_testf4():
    lib = gdstk.Library("Library")
    ref_cell = lib.new_cell("Base")
    pol = gdstk.Polygon(((0, 0), (0, 1), (0.2, 1), (0.2, 0.8), (0.1, 0.8), (0.1, 0.6), (0.2, 0.6), (0.2, 0.4), (0.1, 0.4), (0.1, 0)))
    ref_cell.add(pol)
    cell = lib.new_cell("Ref")
    ref = gdstk.Reference(ref_cell,(0,0))
    ref.repetition = gdstk.Repetition(offsets=[(0.5, 1), (2, 0), (1.5, 0.5)])
    cell.add(ref)
    lib.write_gds("/d/Downloads/ref_f_rep4.gds")
    lib.write_oas("/d/Downloads/ref_f_rep4.oas")

ref_testf4()

def ref_testf5():
    lib = gdstk.Library("Library")
    ref_cell = lib.new_cell("Base")
    pol = gdstk.Polygon(((0, 0), (0, 1), (0.2, 1), (0.2, 0.8), (0.1, 0.8), (0.1, 0.6), (0.2, 0.6), (0.2, 0.4), (0.1, 0.4), (0.1, 0)))
    ref_cell.add(pol)
    cell = lib.new_cell("Ref")
    ref = gdstk.Reference(ref_cell,(0,0))
    ref.repetition = gdstk.Repetition(columns=1, rows=2, v1=(3,3), v2=(4,4))
    cell.add(ref)
    lib.write_gds("/d/Downloads/ref_f_rep5.gds")
    lib.write_oas("/d/Downloads/ref_f_rep5.oas")

ref_testf5()

def ref_test():
    lib = gdstk.Library("Library")
    ref_cell = lib.new_cell("Base")
    pol = gdstk.rectangle((0, 0), (1, 1))
    ref_cell.add(pol)
    cell = lib.new_cell("Ref")
    ref = gdstk.Reference(ref_cell,(0,0))
    ref.repetition = gdstk.Repetition(x_offsets=(1, 3, -2))
    cell.add(ref)
    lib.write_gds("/d/Downloads/ref_rectangle_rep.gds")
    lib.write_oas("/d/Downloads/ref_rectangle_rep.oas")

ref_test()

def ref_test2():
    lib = gdstk.Library("Library")
    ref_cell = lib.new_cell("Base")
    pol = gdstk.rectangle((0, 0), (1, 1))
    ref_cell.add(pol)
    cell = lib.new_cell("Ref")
    ref = gdstk.Reference(ref_cell,(0,0))
    ref.repetition = gdstk.Repetition(columns=2, rows=1, v1=(3,3), v2=(4,4))
    cell.add(ref)
    lib.write_gds("/d/Downloads/ref_rectangle_rep2.gds")
    lib.write_oas("/d/Downloads/ref_rectangle_rep2.oas")

ref_test2()

def ref_test3():
    lib = gdstk.Library("Library")
    ref_cell = lib.new_cell("Base")
    pol = gdstk.rectangle((0, 0), (1, 1))
    ref_cell.add(pol)
    cell = lib.new_cell("Ref")
    ref = gdstk.Reference(ref_cell,(0,0))
    ref.repetition = gdstk.Repetition(columns=2, rows=3, spacing=(5,5))
    cell.add(ref)
    lib.write_gds("/d/Downloads/ref_rectangle_rep3.gds")
    lib.write_oas("/d/Downloads/ref_rectangle_rep3.oas")

ref_test3()

def ref_test4():
    lib = gdstk.Library("Library")
    ref_cell = lib.new_cell("Base")
    pol = gdstk.rectangle((0, 0), (1, 1))
    ref_cell.add(pol)
    cell = lib.new_cell("Ref")
    ref = gdstk.Reference(ref_cell,(0,0))
    ref.repetition = gdstk.Repetition(offsets=[(0.5, 1), (2, 0), (1.5, 0.5)])
    cell.add(ref)
    lib.write_gds("/d/Downloads/ref_rectangle_rep4.gds")
    lib.write_oas("/d/Downloads/ref_rectangle_rep4.oas")

ref_test4()

def ref_test5():
    lib = gdstk.Library("Library")
    ref_cell = lib.new_cell("Base")
    pol = gdstk.rectangle((0, 0), (1, 1))
    ref_cell.add(pol)
    cell = lib.new_cell("Ref")
    ref = gdstk.Reference(ref_cell,(0,0))
    ref.repetition = gdstk.Repetition(columns=1, rows=2, v1=(3,3), v2=(4,4))
    cell.add(ref)
    lib.write_gds("/d/Downloads/ref_rectangle_rep5.gds")
    lib.write_oas("/d/Downloads/ref_rectangle_rep5.oas")

ref_test5()

