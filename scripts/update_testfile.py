import os
from pathlib import Path

'''
generates the file svg_testcases.hpp.
'''

TESTCASES_SVG_DIR = "testcases/geometry"
TESTCASE_NAME = "test/svg_testcases.cpp"

def construct_tst_case(name):

    return f'''TEST_CASE("{name}"){{\n    do_test("{name}");\n}}\n\n'''

def construct_tst_file():
    subfolders = [Path(f.path).name for f in os.scandir(TESTCASES_SVG_DIR) if f.is_dir() ]
    subfolders.sort(reverse=True)

    cases_str= ''.join([construct_tst_case(i) for i in subfolders])

    includes = '''#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN\n#include "test_utilities.hpp"\n#include "doctest.h"\n\n'''

    tst_start = 'TEST_SUITE_BEGIN("svg geometry testcases");\n\n'
    comment = "//This file was auto-generated from `scripts/update_testfile.py`.\n\n"

    tst_end = "TEST_SUITE_END;\n"

    return comment + includes + tst_start+ cases_str + tst_end


write_p = Path(TESTCASE_NAME)
with write_p.open("w", encoding="utf-8") as f:
    f.write(construct_tst_file())

