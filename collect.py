#!/usr/bin/python

import os, sys
import math

def rstda(fname):

    sys.stdout.write("#%9s%10s%10s\n"%("e (eV)", "f", "S"))

    excs = []
    with open(fname) as pfname:

        for line in pfname:

            if line.find("excitation energies") >= 0:
                line = pfname.next()
                line = pfname.next()
                words = line.split()

                while len(words) > 0:
                    #          0  6        18 24        36 42
                    # Ci line[47:53], line[65:71], line[83:89]
                    #          7 11
                    # Oc line[54:58], line[72:76], line[90:94]
                    #         13 17
                    # Vi line[60:64], line[78:82], line[96:100]

                    S = 0.
                    wmax = max([float(line[47:53])**2., float(line[65:71])**2., float(line[83:89])**2.])
                    for i in range(3):
                        wi = float(line[47+i*18:53+i*18])**2. / wmax
                        S += wi

                    if float(words[3]) > 0.:
                        sys.stdout.write("%10.4f%10.4f%10.4f\n"%(float(line[5:14]), float(line[22:33]), S/3.))

                    line = pfname.next()
                    words = line.split()
                    
    return 0

def main():

    stda = None
    
    i = 1
    while i < len(sys.argv):
    
        arg = sys.argv[i]
    
        if arg == '-stda':
            i += 1
            arg = sys.argv[i]
            if os.path.isfile(arg):
                stda = arg
        else:
            sys.stderr.write("%s: %s: keyword does not exist\n"%(sys.argv[0], arg))
            sys.exit(1)
    
        i += 1
    
    if stda == None:
        sys.stderr.write("%s: stda input file needed\n"%(sys.argv[0]))
        sys.exit(1)

    rstda(stda)

    return 0

if __name__ == "__main__":

    main()

