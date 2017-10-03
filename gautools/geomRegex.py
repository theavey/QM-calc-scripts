__author__ = 'Thomas Heavey'

import re

filename = "testg.out"

def findallgeoms(filename):
    """A function that takes a file name and returns a list of
    geometries. Works with Gaussian output, haven't checked with
    Q-Chem."""
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
        allxyz = []
        geoms = geomregex.finditer(file.read())
        for geom in geoms:
            thisxyz = []
            mlgeom = geom.group(1)
            for line in mlgeom.split('\n'):
                # Ignore blank lines:
                if len(line) < 2:
                    continue
                xyzelemstring = [line.split()[i] for i in relevantelem]
                xyzelemnum = [float(i) for i in xyzelemstring]
                xyzelemnum[0] = int(xyzelemstring[0])
                thisxyz.append(xyzformat.format(*xyzelemnum))
            allxyz.append(thisxyz)

    return(allxyz)
# I don't know if I like this format. It would be reasonable for
# Mathematica, but somewhat odd for Python. I guess for outputting
# it though it won't be terrible because I can just double
# iterate over the nested list, writing lines from the strings.
# I'll need to pick a separator for between geometries maybe but that's
# not a problem. Also with this format, should be easy to count number
# of atoms.

# Still need to have way to just find stationary points

if __name__ == "__main__":
    print(findallgeoms(filename))
# Ugly because returned as list of list of strings
