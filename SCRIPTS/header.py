# Local Variables:
# mode: python
# coding: utf-8
# End:
import subprocess
from pprint import pprint

FOR_HDR = """!! This file is part of SEM
!!
!! Copyright CEA, ECP, IPGP
!!
"""
"""
!! The fact that you are presently reading this means that you have had
!! knowledge of the CeCILL-C license and that you accept its terms.
!! See the file LICENSE.txt included with this distribution for details.
!!
"""

FOR_FOOT = """
!! Local Variables:
!! mode: f90
!! show-trailing-whitespace: t
!! coding: utf-8
!! f90-do-indent: 4
!! f90-if-indent: 4
!! f90-type-indent: 4
!! f90-program-indent: 4
!! f90-continuation-indent: 4
!! End:
!! vim: set sw=4 ts=8 et tw=80 smartindent :
"""

C_FOOT = """
/* Local Variables:                                                        */
/* mode: c++                                                               */
/* show-trailing-whitespace: t                                             */
/* coding: utf-8                                                           */
/* c-file-style: "stroustrup"                                              */
/* End:                                                                    */
/* vim: set sw=4 ts=8 et tw=80 smartindent :                               */
"""

PY_FOOT = ""

def for_to_c(s):
    return "\n".join([ l.replace("!!","/*").ljust(75)+"*/" for l in s.splitlines() ])

def for_to_py(s):
    return "\n".join([ l.replace("!!","#") for l in s.splitlines() ])
C_HDR = for_to_c(FOR_HDR)
PY_HDR = for_to_py(FOR_HDR)

EXCLUDES = [
    "COMMON/splib.F90",
    "COMMON/blas.F90",
    "COMMON/file_scan.c",
    "COMMON/file_scan.h",
    "COMMON/mpif.h",
    ]

fnames = subprocess.check_output(["git","ls-files"])

fortranfiles = []
cfiles = []
pyfiles = []
external = []
for fn in fnames.split():
    if fn.startswith("MESH/metis/"):
        external.append(fn)
        continue
    if fn.startswith("COMMON/netlib/"):
        external.append(fn)
        continue
    if fn.startswith("SEM3D/SRC/BLAS"):
        external.append(fn)
        continue
    if fn in EXCLUDES:
        external.append(fn)
        continue
    lwfn = fn.lower()
    if lwfn.endswith(".f90") or lwfn.endswith(".f"):
        fortranfiles.append(fn)
    if lwfn.endswith(".c") or lwfn.endswith(".h") or lwfn.endswith(".cpp") or lwfn.endswith(".hpp"):
        cfiles.append(fn)
    if lwfn.endswith(".py"):
        pyfiles.append(fn)

print "FORTRAN:", len(fortranfiles)
pprint(fortranfiles)
print
print "C/C++:", len(cfiles)
pprint(cfiles)
print
print "PYTHON:", len(pyfiles)
pprint(pyfiles)

def replace_head(lines, hdrlines, comment):
    lines[0:0] = hdrlines

def is_for_comment(line):
    return line.strip().startswith("!")

def is_c_comment(line):
    return line.strip().startswith("//") or line.strip().startswith("/*")

def is_py_comment(line):
    return line.strip().startswith("#")

def replace_foot(lines, footlines, is_comment):
    # find the last line that is not blank or a comment
    k = len(lines)-1
    while k>=0:
        if is_comment(lines[k]):
            k=k-1
            continue
        break
    found = False
    # Now move forward to find "Local Variables:" emacs tag
    while k<len(lines):
        if "Local Variables:" in lines[k]:
            found = True
            break
        k = k+1
    lines[k:] = footlines

def replace_head_foot(fn, hdr, foot, comment):
    lines = [ l[:-1] for l in file(fn).readlines() ]
    footlines = foot.splitlines()
    nf = len(footlines)
    hdrlines = hdr.splitlines()
    nh = len(hdrlines)
    changed = False
    if lines[:nh] != hdrlines:
        replace_head(lines, hdrlines, comment)
        changed = True
    if lines[-nf:]!= footlines:
        replace_foot(lines, footlines, comment)
        changed = True
    if changed:
        file(fn,"w").write("\n".join(lines))

def add_newline_at_eof(fn):
    data = file(fn).read()
    file(fn,"w").write(data+"\n")

#for fn in fortranfiles:
#    replace_head_foot(fn, FOR_HDR, FOR_FOOT, is_for_comment)

for fn in fortranfiles+cfiles:
    add_newline_at_eof(fn)

#for fn in cfiles:
#    replace_head_foot(fn, C_HDR, C_FOOT, is_c_comment)

print C_HDR
print C_FOOT
