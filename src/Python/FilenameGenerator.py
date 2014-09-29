import sys
import os

if len(sys.argv) < 4:
    print "Not enough arguments"
    sys.exit(1)



# gather the args
num_entries     = int(sys.argv[1])
template        = sys.argv[2]
output_filename = sys.argv[3]


# check the template
if template.find("%d") < 0:
    print "Template not in the correct format, must have a %d in it"
    sys.exit(1)


# open the file
f = open(output_filename, 'w')

for n in range(num_entries):
    fname = template % n
    f.write(fname + '\n')

f.close()


