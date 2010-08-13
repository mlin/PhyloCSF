from sys import argv, exit
from os import system

escape = False
argv2 = []
for arg in argv[1:len(argv)]:
	if arg == "-pp":
		escape = True
	elif escape:
		escape = False
	else:
		argv2.append(arg)
cmd = ' '.join(['ocamllex'] + argv2)
exit(system(cmd))
