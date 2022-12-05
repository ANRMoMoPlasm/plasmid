#!/usr/bin/python

import os, sys
import math

def readtyp(line):

    words = line.split()
    otyp = words[0]
    nfnc = int(words[1])

    return otyp, nfnc

def selectorbital(fbasis):

    ngauss = {"s": 1, "p": 3, "d": 6, "f": 10}

    nfncs = 0
    select = []
    for fb in fbasis:

        for i in range(len(fb)-1):

            nfncs += ngauss[fb[i]]

            if fb[i] == "s": 
                if fb[i+1] == "p":
                    select.append(nfncs)

        nfncs += ngauss[fb[i+1]]

    return nfncs, select

def rmolden(fname):

    with open(fname) as pfname:

        for line in pfname:

            if line.find("[Atoms]") >= 0:

                line = pfname.next()

                ## count the number of atoms
                natom = 0
                while line.find("[GTO]") < 0:
                    natom += 1
                    line = pfname.next()

                fbasis = []
                ## count the number of basis functions
                for iatom in range(natom):
                    line = pfname.next()
                    line = pfname.next()
                    words = line.split()
   
                    fb = [] 
                    while len(words) > 0:
                        otyp, nfnc = readtyp(line)
                        fb.append(otyp)
                        for i in range(nfnc):
                            line = pfname.next()
    
                        line = pfname.next()
                        words = line.split()

                    fbasis.append(fb)
   
                nfncs, supats = selectorbital(fbasis)

                ## start to read the nfncs MO containing nfncs basis
                line = pfname.next()
                line = pfname.next()
 
                Smo = []; ihomo = 0
                for imo in range(nfncs): 
                    #print "MO=%i"%(imo+1)
                    for i in range(3):
                        line = pfname.next()

                    if float(line.split()[1]) == 2.:
                        ihomo += 1
                    else:
                        break

                    ci = 0.; ctot = 0.
                    for ibf in range(nfncs):
                        line = pfname.next()
                        words = line.split()

                        if int(words[0]) in supats:
                            ci += float(words[1])**2.

                        ctot += float(words[1])**2.

                    #print "MO=%i -> %20.10f"%(imo+1, math.sqrt(ci/ctot))
                    Smo.append(ci/ctot)

                    if imo < nfncs-1:
                        line = pfname.next()

    return ihomo, Smo

def rstda(fname, ihomo, Smo):

    sys.stdout.write("#%9s%10s%10s\n"%("e (eV)", "f", "S"))

    excs = []
    with open(fname) as pfname:

        for line in pfname:

            if line.find("oMOs in TDA:") >= 0:
                words = line.split()
                nocc = int(words[-1])

            elif line.find("vMOs in TDA:") >= 0:
                words = line.split()
                nvir = int(words[-1])

            elif line.find("excitation energies") >= 0:
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
                    wtot = float(line[47:53])**2. + float(line[65:71])**2. + float(line[83:89])**2.
                    for i in range(3):
                        iocc = int(line[54+i*18:58+i*18]) - nocc + ihomo
                        wi = float(line[47+i*18:53+i*18])**2. / wtot
                        S += wi * Smo[iocc-1]

                    if float(words[3]) > 0.:
                        if S >= 0.1:
                            sys.stdout.write("%10.4f%10.4f%10.4f\n"%(float(line[5:14]), float(line[22:33]), S))

                    line = pfname.next()
                    words = line.split()
                    
    return 0

def main():

    stda = None
    mldn = None
    
    i = 1
    while i < len(sys.argv):
    
        arg = sys.argv[i]
    
        if arg == '-stda':
            i += 1
            arg = sys.argv[i]
            if os.path.isfile(arg):
                stda = arg
        elif arg == '-molden':
            i += 1
            arg = sys.argv[i]
            if os.path.isfile(arg):
                mldn = arg
        else:
            sys.stderr.write("%s: %s: keyword does not exist\n"%(sys.argv[0], arg))
            sys.exit(1)
    
        i += 1
    
    if stda == None:
        sys.stderr.write("%s: stda input file needed\n"%(sys.argv[0]))
        sys.exit(1)

    if mldn == None:
        sys.stderr.write("%s: molden input file needed\n"%(sys.argv[0]))
        sys.exit(1)

    ihomo, Smo = rmolden(mldn)
    rstda(stda, ihomo, Smo)

    return 0

if __name__ == "__main__":

    main()

