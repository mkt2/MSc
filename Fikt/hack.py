import collections, sys
sys.path.insert(0, '/home/markus/Project-hack')

def foo(x):
	for i in x:
		print i,i

if __name__ == "__main__":
	foo(sys.argv)
	k = 3

#Command line with ipython:
"""
ipython
%run hack.py a b c d
a
b
c
d
breyta linu 5: print i,i. save-a
reload(hack)
from hack import foo
foo(["a","b","c","d"])
a a
b b
c c
d d
"""