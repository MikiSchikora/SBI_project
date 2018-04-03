from distutils.core import setup
import os

setup(name='BYCProject', version='1.0', description='Build a biological macrocomplex from pairwise interactions',
	author='Miquel Angel Schikora Tamarit and Marina Lleal Custey',
	author_email='miquelangel.schikora01@estudiant.upf.edu, marina.lleal01@estudiant.upf.edu',
	url='https://github.com/MikiSchikora/SBI_project',
	packages=['byc'], scripts=['byc/run_byc.py','byc/refine_model.py'])