#!/usr/bin/env python

import spiraltraj

print('basic test: spiral out trajectory calculation with default settings:')
res = spiraltraj.calc_traj(spiraltype=1)
print('shape of returned list: [', len(res), ', ', len(res[0]), ']')
