#!/usr/bin/python
# encoding: utf-8
#
# Ceres Solver - A fast non-linear least squares minimizer
# Copyright 2015 Google Inc. All rights reserved.
# http://ceres-solver.org/
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice,
#   this list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
# * Neither the name of Google Inc. nor the names of its contributors may be
#   used to endorse or promote products derived from this software without
#   specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.
#
# Author: sameeragarwal@google.com (Sameer Agarwal)
#
# Note: You will need Sphinx and Pygments installed for this to work.

from __future__ import print_function
import glob
import io
import os
import sys

# Number of arguments
N = len(sys.argv)

if N < 3:
  print('make_docs.py src_root destination_root')
  sys.exit(1)

src_dir    = sys.argv[1] + '/docs/source'
build_root = sys.argv[2]
cache_dir  = build_root + '/doctrees'
html_dir   = build_root + '/html'

# Called from Command Line
if N == 3:
  sphinx_exe = 'sphinx-build'

# Called from CMake (using the SPHINX_EXECUTABLE found)
elif N == 4:
  sphinx_exe = sys.argv[3]

# Run Sphinx to build the documentation.
os.system('%s -b html -d %s %s %s' %(sphinx_exe, cache_dir, src_dir, html_dir))

replacements = [
  # By default MathJax uses does not use TeX fonts. This simple search
  # and replace fixes that.
  ('''config=TeX-AMS-MML_HTMLorMML"></script>''',
   '''config=TeX-AMS_HTML">
      MathJax.Hub.Config({
          "HTML-CSS": {
            availableFonts: ["TeX"]
          }
        });
      </script>'''),

  # The title for the homepage is not ideal, so change it.
  ('<title>Ceres Solver &mdash; Ceres Solver</title>',
   '<title>Ceres Solver &mdash; A Large Scale Non-linear Optimization Library</title>')
]

# This is a nasty hack to strip the breadcrumb navigation. A better strategy is
# to fork the upstream template, but that is no fun either. Whitespace matters!
# This doesn't use regular expressions since the escaping makes it untenable.
breadcrumb_start_other = \
'''<div role="navigation" aria-label="breadcrumbs navigation">
  <ul class="wy-breadcrumbs">
    <li><a href="index.html">Docs</a> &raquo;</li>

    <li>'''

# The index page has a slightly different breadcrumb.
breadcrumb_start_index = breadcrumb_start_other.replace('index.html', '#')

breadcrumb_end = \
'''</li>
      <li class="wy-breadcrumbs-aside">

      </li>
  </ul>
  <hr/>
</div>'''

for name in glob.glob('%s/*.html' % html_dir):
  print('Postprocessing: ', name)
  with io.open(name, encoding="utf-8") as fptr:
    out = fptr.read()

  for input_pattern, output_pattern in replacements:
    out = out.replace(input_pattern, output_pattern)

  try:
    breadcrumb_start = breadcrumb_start_index \
                       if name.endswith('index.html') \
                       else breadcrumb_start_other
    pre_breadcrumb_start, post_breadcrumb_start = out.split(breadcrumb_start)
    title, post_breadcrumb_end = post_breadcrumb_start.split(breadcrumb_end)
    print('Stripping breadcrumb for -', title)
    out = pre_breadcrumb_start + post_breadcrumb_end
  except ValueError:
    print('Skipping breadcrumb strip for', name)

  with io.open(name, 'w', encoding="utf-8") as fptr:
    fptr.write(out)
