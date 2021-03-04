from distutils.core import setup, Extension

module1 = Extension('spiraltraj',
                    sources = ['gen_spiral.cpp', 'vdspiral.cpp', 'nonCartesianTraj.cpp'])

setup (name = 'SpiralTraj',
       version = '1.0',
       description = 'This package calculate a spiral trajectory',
       ext_modules = [module1])
