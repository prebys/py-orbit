#!/usr/bin/env python

#---------------------------------------------------------
#This will replace 0.0000000 by 0.0 
#---------------------------------------------------------

import sys

if( len(sys.argv) != 3 ):
	print "Usage: >python ",sys.argv[0]," transition.txt  new_transition.txt"
	sys.exit(1)
	
fl_in = open(sys.argv[1],"r")
fl_out = open(sys.argv[2],"w")

count = 0
for line in fl_in:
	if(count%1000 == 0): print "line #=",count
	trans_part = line[0:25]
	data_arr = line[25:].split()
	#print "trans_part=",trans_part
	st = trans_part[:]
	for ind in range(len(data_arr)):
		dt = data_arr[ind]
		#print "ind=",ind," str=",dt
		if(dt == "0.000000000000000000000000000000" ):
			st = st + " 0.0 "
		else:
			st = st + " " + dt + "  "
	fl_out.write(st+"\n")
	#print "res s=",st
	count = count + 1

fl_out.close()
fl_in.close()
