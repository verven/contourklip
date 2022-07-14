import subprocess
import textwrap

'''
This script creates the single-header file. 
It is quite manual since it does not automatically collect and sort the dependencies.
it uses the utility "unifdef" (https://dotat.at/prog/unifdef/) which is essentially a
preprocessor which can conditionally define / undefine symbols. Unlike eg. gcc -E, it will
only expand macros and not modify the rest, so it is suitable to generate a readable header.
'''

PREAMBLE = '''
/* 
Contourklip, a contour clipping library which supports cubic beziers.

Copyright (C) 2022 verven [ vervencode@protonmail.com ]

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

'''

ROOT_HEADER = "include/polyclip.hpp"
OUTPUT = "single_include/contourklip.hpp"
HEADER_NAME = "CONTOURKLIP_CONTOURKLIP_HPP"

# combine files
# unifdef -o single_include/contourklip.hpp -UDEBUG -DCONTOURKLIP_SVG_IO_HPP single_include/contourklip.hpp
UDEBUG_COMMAND = ['unifdef', '-o', OUTPUT, '-UDEBUG', '-DCONTOURKLIP_SVG_IO_HPP', OUTPUT]

ORDERED_FILES = [ f"include/{i}" for i in [
    "direct_solvers.hpp",
    "polynomial_solver.hpp",
    "geometry_base.hpp",
    "bezier_utils.hpp",
    "sweep_point.hpp",
    "contour_postprocessing.hpp",
    "polyclip.hpp",
    ]
]


with open(OUTPUT, 'w+') as outfile:
    for f in ORDERED_FILES:
        outfile.write("\n")
        with open(f) as infile:
            for line in infile:
                outfile.write(line)
    #
    outfile.write("\n")

subprocess.run(UDEBUG_COMMAND)

std_includes = set()
data = []

with open(OUTPUT, "r+") as outfile:
    for line in outfile:
        if not "#include" in line:
            data.append(line)
        if "#include <" in line:
            incl = line[:line.find(">")+1].strip()
            std_includes.add(incl)

    outfile.seek(0)
    outfile.truncate()
    outfile.write(PREAMBLE + "#ifndef " + HEADER_NAME + "\n#define " + HEADER_NAME
            + "\n" )
    for i in sorted(std_includes):
            outfile.write(f"{i}\n")
    for i in data:
        outfile.write(i)
    outfile.write("\n#endif //" + HEADER_NAME)