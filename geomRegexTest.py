__author__ = 'Thomas Heavey'

import re

filename = "testg.out"

def findgeoms(filename):
    """A function that takes a file name and returns a list of
    geometries."""
    relevantelem = [1,3,4,5]
    xyzformat = '{:>2}      {: f}   {: f}   {: f}'
    geomregex = re.compile(
        r'(?:Standard orientation)' # non-capturing (nc) start of geometry
        r'(?:.+?)'                  # nc geometry header
        r'((?:(?:\s+\d+\s+)'        # nc atom number
        r'(\d+\s+)'                 # (capturing) atomic number
        r'(?:\d+\s+)'               # nc atomic type
        r'(-?\d+\.\d+\s*){3,3}'     # 3 cartesian coordinates (x,y,z)
        r')+)'                      # repeat for at least one atom
        r'(?:-)'                    # nc end at line of dashes
        , re.DOTALL)

    with open(filename, 'r') as file:
        geoms = geomregex.search(file.read())
        print(geoms.group(1))
        mlgeoms = geoms.group(1)
        for line in mlgeoms.split('\n'):
            # Ignore blank lines:
            if len(line) < 2:
                continue
            xyzelemstring = [line.split()[i] for i in relevantelem]
            xyzelemnum = [float(i) for i in xyzelemstring]
            xyzelemnum[0] = int(xyzelemstring[0])
            print(xyzformat.format(*xyzelemnum))

findgeoms(filename)