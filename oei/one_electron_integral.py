#!---------------------------------------------------------------------!
#! Written by Madu Manathunga on 07/01/2021                            !
#!                                                                     !
#! Copyright (C) 2020-2021 Merz lab                                    !
#! Copyright (C) 2020-2021 Götz lab                                    !
#!                                                                     !
#! This Source Code Form is subject to the terms of the Mozilla Public !
#! License, v. 2.0. If a copy of the MPL was not distributed with this !
#! file, You can obtain one at http://mozilla.org/MPL/2.0/.            !
#!_____________________________________________________________________!

#!---------------------------------------------------------------------!
#! This source file contains classes necessary for generating one      !
#! electron integrals.                                                 !
#!---------------------------------------------------------------------!

import params

# parent class for all one electron integrals
class OEint:
    # set file handlers
    fhc = 0 # file handler for class declarations
    fhd = 0 # file handler for function implementations

    # max_m ranges from 0 to a+b; where a and b are the angular momentum of i and j
    # of the integral being considered <i|j>. If max_m=0, the integral is a true integral
    # and m>0 results in auxliary integrals. 
    def __init__(self,max_m=0):
        self.max_m=max_m
        if self.fhc == 0 or self.fhd == 0:
            print("Warning: file handlers in OEint class are not set. \n")

    # generate source code for integrals
    def gen_int(self):
        pass 


# <s|s> class, a subclass of OEint
class SSint(OEint):
    pass

# <p|s> class, a subclass of OEint
class PSint(OEint):
    def gen_int(self):
            # write code paths for integrals. Note that we use C++ classes here. 
            for m in range(0,self.max_m+1):
                if m == 0:
                    self.fhc.write("\n/* PS true integral, m=%d */ \n" % (m)) 
                    self.fhd.write("\n/* PS true integral, m=%d */ \n" % (m)) 
                else:
                    self.fhc.write("\n/* PS auxilary integral, m=%d */ \n" % (m))
                    self.fhd.write("\n/* PS auxilary integral, m=%d */ \n" % (m))  

                self.fhc.write("class PSint_%d { \n" % (m))
                self.fhc.write("public: \n")

                # set labels for improving readability 
                lbl=["Px", "Py", "Pz"]

                # write class variables
                for i in range(0,3):
                    self.fhc.write("  QUICKDouble x_%d_%d_%d; // %s, %s \n" % (i, 0, m, lbl[i], "S"))

                # write class functions
                self.fhc.write("  __device__ __inline__ PSint_%d(QUICKDouble* AA, QUICKDouble* CC, QUICKDouble* PP); \n" % (m))
                self.fhc.write("}; \n")

                self.fhd.write("__device__ __inline__ PSint_%d::PSint_%d(QUICKDouble* AA, QUICKDouble* CC, QUICKDouble* PP){ \n\n" % (m, m))
                for i in range(0,3):
                    self.fhd.write("  x_%d_%d_%d=(PP[%d]-AA[%d]) * LOCVY(0, 0, %d) - (PP[%d]-CC[%d]) * LOCVY(0, 0, %d)\n" % (i, 0, m, i, i, m, i, i, m))

                self.fhd.write("} \n\n")                


# <s|p> class, a subclass of OEint
class SPint(OEint):
    def gen_int(self):
            # write code paths for integrals. Note that we use C++ classes here. 
            for m in range(0,self.max_m+1):
                if m==0:
                    self.fhc.write("\n/* SP true integral, m=%d */ \n" % (m)) 
                    self.fhd.write("\n/* SP true integral, m=%d */ \n" % (m)) 
                else:
                    self.fhc.write("\n/* SP auxilary integral, m=%d */ \n" % (m))
                    self.fhd.write("\n/* SP auxilary integral, m=%d */ \n" % (m))                

                self.fhc.write("class SPint_%d { \n" % (m))
                self.fhc.write("public: \n")

                # set labels for improving readability 
                lbl=["Px", "Py", "Pz"]

                # write class variables
                for i in range(0,3):
                    self.fhc.write("  QUICKDouble x_%d_%d_%d; // %s, %s \n" % (0, i, m, "S", lbl[i]))

                # write class functions
                self.fhc.write("  __device__ __inline__ SPint_%d(QUICKDouble* BB, QUICKDouble* CC, QUICKDouble* PP); \n" % (m))
                self.fhc.write("}; \n")

                self.fhd.write("__device__ __inline__ SPint_%d::SPint_%d(QUICKDouble* BB, QUICKDouble* CC, QUICKDouble* PP){ \n\n" % (m, m))

                for i in range(0,3):
                    self.fhd.write("  x_%d_%d_%d=(PP[%d]-BB[%d]) * LOCVY(0, 0, %d) - (PP[%d]-CC[%d]) * LOCVY(0, 0, %d)\n" % (0, i, m, i, i, m, i, i, m))

                self.fhd.write("} \n\n") 


# <p|p> class, subclass of OEint
class PPint(OEint):
    def gen_int(self):
        # write code paths for integrals. Note that we use C++ classes here.
        for m in range(0,self.max_m+1):
            if m == 0:
                self.fhc.write("\n/* PP true integral, m=%d */ \n" % (m)) 
                self.fhd.write("\n/* PP true integral, m=%d */ \n" % (m)) 
            else:
                self.fhc.write("\n/* PP auxilary integral, m=%d */ \n" % (m))
                self.fhd.write("\n/* PP auxilary integral, m=%d */ \n" % (m))          
        
            self.fhc.write("class PPint_%d { \n" % (m))
            self.fhc.write("public: \n")

            # set labels for improving readability 
            lbl=["Px", "Py", "Pz"]

            # write class variables
            for i in range(0,3):
                for j in range(0,3):
                    self.fhc.write("  QUICKDouble x_%d_%d_%d; // %s, %s \n" % (i, j, m, lbl[i], lbl[j]))

            # write class functions
            self.fhc.write("  __device__ __inline__ PPint_%d(QUICKDouble* BB, QUICKDouble* CC, QUICKDouble* PP, QUICKDouble ABCD); \n" % (m))          
            self.fhc.write("}; \n")

            # write function definitions
            self.fhd.write("__device__ __inline__ PPint_%d::PPint_%d(QUICKDouble* BB, QUICKDouble* CC, QUICKDouble* PP, QUICKDouble ABCD){ \n\n" % (m, m))
            self.fhd.write("  PSint ps_%d(BB, CC, PP); // construct [p|s] for m=%d \n" % (m, m))
            self.fhd.write("  PSint ps_%d(BB, CC, PP); // construct [p|s] for m=%d \n" % (m+1, m+1))

            idx=1
            for i in range(0,3):
                for j in range(0,3):
                    self.fhd.write("  x_%d_%d_%d = (PP[%d]-BB[%d]) * ps_%d.x_%d_%d_%d - (PP[%d]-CC[%d]) * ps_%d.x_%d_%d_%d \n" % (i, j, m, j, j, m, i, 0, m,\
                    j, j, m+1, i, 0, m+1))

                    if i == j:
                        self.fhd.write("  x_%d_%d_%d += 0.5/ABCD * (LOCVY(0, 0, %d) - LOCVY(0, 0, %d)) \n" % (i, j, m, m, m+1))

                    idx += 1
            self.fhd.write("\n } \n")


# <d|s> class, subclass of OEint
class DSint(OEint):
    def gen_int(self):
        # write code paths for integrals. Note that we use C++ classes here.
        for m in range(0,self.max_m+1):
            if m == 0:
                self.fhc.write("\n/* DS true integral, m=%d */ \n" % (m)) 
                self.fhd.write("\n/* DS true integral, m=%d */ \n" % (m)) 
            else:
                self.fhc.write("\n/* DS auxilary integral, m=%d */ \n" % (m))
                self.fhd.write("\n/* DS auxilary integral, m=%d */ \n" % (m))              

            self.fhc.write("class DSint_%d { \n" % (m))
            self.fhc.write("public: \n")

            # set labels for improving readability 
            lbl=["Dxy", "Dyz", "Dxz", "Dxx", "Dyy", "Dzz"]

            # write class variables
            for i in range(0,6):
                self.fhc.write("  QUICKDouble x_%d_%d_%d; // %s, %s \n" % (i, 0, m, lbl[i], "S"))

            # write class functions
            self.fhc.write("  __device__ __inline__ DSint_%d(QUICKDouble* AA, QUICKDouble* CC, QUICKDouble* PP, QUICKDouble ABCD); \n" % (m))          
            self.fhc.write("}; \n")

            # write function definitions
            self.fhd.write("__device__ __inline__ DSint_%d::DSint_%d(QUICKDouble* AA, QUICKDouble* CC, QUICKDouble* PP, QUICKDouble ABCD){ \n\n" % (m, m))
            self.fhd.write("  PSint ps_%d(AA, CC, PP); // construct [p|s] for m=%d \n" % (m, m))
            self.fhd.write("  PSint ps_%d(AA, CC, PP); // construct [p|s] for m=%d \n" % (m+1, m+1))
            self.fhd.write("  PSint ps_%d(AA, CC, PP); // construct [p|s] for m=%d \n" % (m+2, m+2))

            for i in range(0,6):
                tmp_mcal=[params.Mcal[i+4][0], params.Mcal[i+4][1], params.Mcal[i+4][2]]
                for j in range(0,3):
                    if params.Mcal[i+4][j] != 0:
                        tmp_mcal[j] -= 1
                        tmp_i=params.trans[tmp_mcal[0]][tmp_mcal[1]][tmp_mcal[2]]
                        self.fhd.write("  x_%d_%d_%d = (PP[%d]-AA[%d]) * ps_%d.x_%d_%d_%d - (PP[%d]-CC[%d]) * ps_%d.x_%d_%d_%d \n" % (i, 0, m, j, j, m, tmp_i-1, 0, m,\
                        j, j, m+1, tmp_i-1, 0, m+1))

                        if params.Mcal[i+4][j] == 2:
                            self.fhd.write("  x_%d_%d_%d += 0.5/ABCD * (LOCVY(0, 0, %d) - LOCVY(0, 0, %d)) \n" % (i, 0, m, m, m+1))

                        break
            self.fhd.write("\n } \n")

def write_oei():

    # set files
    OEint.fhc = open('gpu_oei_classes.h','w')
    OEint.fhd = open('gpu_oei_definitions.h','w')

    # generate integrals for systems containing only s and p functions.
    # <s|s> is trivial and we wont auto generate it.

    # generate <p|s>
    ps=PSint(2)
    ps.gen_int() 

    # generate <s|p>
    sp=SPint(0)
    sp.gen_int()

    # generate <p|p>
    pp=PPint(0)
    pp.gen_int()

    # generate <d|s>
    ds=DSint(0)
    ds.gen_int()

    OEint.fhc.close()
    OEint.fhd.close()

write_oei()
