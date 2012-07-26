import os, sys


in1 = open(('water' if len(sys.argv)<3 else sys.argv[1]) + '.xyz')
in2 = open(('test' if len(sys.argv)<3 else sys.argv[2]) + '.xyz')
out = open('compare.xyz', 'w')

for i,line1 in enumerate(in1):
	line2 = in2.readline()
	if len(line1) <= 6:
		if i > 10000: break
		if not line1 or not line2: continue
		if line1[0].isdigit(): out.write(str(int(line1[0])+int(line2[0]))+'\n')
		else: out.write('Atoms\n')
	else:
		out.write( line1 )
		out.write( line2 )

out.close()
os.system("/fs/cbe/g_pc/vmd-1.9 compare.xyz")
