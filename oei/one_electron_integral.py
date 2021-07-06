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
                self.fhc.write("  QUICKDouble x_%d_1; // Px, S \n" % (m))
                self.fhc.write("  QUICKDouble x_%d_2; // Py, S \n" % (m))
                self.fhc.write("  QUICKDouble x_%d_3; // Pz, S \n" % (m))
                self.fhc.write("  __device__ __inline__ PSint_%d(QUICKDouble* BB, QUICKDouble* CC, QUICKDouble* PP); \n" % (m))
                self.fhc.write("}; \n")

                self.fhd.write("__device__ __inline__ PSint_%d::PSint_%d(QUICKDouble* BB, QUICKDouble* CC, QUICKDouble* PP){ \n\n" % (m, m))
                self.fhd.write("  QUICKDouble x_%d_1=(PP[0]-BB[0]) * LOCVY(0, 0, %d)\
 - (PP[0]-CC[0]) * LOCVY(0, 0, %d)\n" % (m, m, m))
                self.fhd.write("  QUICKDouble x_%d_2=(PP[1]-BB[1]) * LOCVY(0, 0, %d)\
 - (PP[1]-CC[1]) * LOCVY(0, 0, %d)\n" % (m, m, m))
                self.fhd.write("  QUICKDouble x_%d_3=(PP[2]-BB[2]) * LOCVY(0, 0, %d)\
 - (PP[2]-CC[2]) * LOCVY(0, 0, %d)\n" % (m, m, m))
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
                self.fhc.write("  QUICKDouble x_%d_1; // S, Px \n" % (m))
                self.fhc.write("  QUICKDouble x_%d_2; // S, Py \n" % (m))
                self.fhc.write("  QUICKDouble x_%d_3; // S, Pz \n" % (m))
                self.fhc.write("  __device__ __inline__ SPint_%d(QUICKDouble* BB, QUICKDouble* CC, QUICKDouble* PP); \n" % (m))
                self.fhc.write("}; \n")

                self.fhd.write("__device__ __inline__ SPint_%d::SPint_%d(QUICKDouble* BB, QUICKDouble* CC, QUICKDouble* PP){ \n\n" % (m, m))
                self.fhd.write("  QUICKDouble x_%d_1=(PP[0]-BB[0]) * LOCVY(0, 0, %d)\
 - (PP[0]-CC[0]) * LOCVY(0, 0, %d)\n" % (m, m, m))
                self.fhd.write("  QUICKDouble x_%d_2=(PP[1]-BB[1]) * LOCVY(0, 0, %d)\
 - (PP[1]-CC[1]) * LOCVY(0, 0, %d)\n" % (m, m, m))
                self.fhd.write("  QUICKDouble x_%d_3=(PP[2]-BB[2]) * LOCVY(0, 0, %d)\
 - (PP[2]-CC[2]) * LOCVY(0, 0, %d)\n" % (m, m, m))
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
            self.fhc.write("  QUICKDouble x_%d_1; // Px, Px \n" % (m))
            self.fhc.write("  QUICKDouble x_%d_2; // Px, Py \n" % (m))
            self.fhc.write("  QUICKDouble x_%d_3; // Px, Pz \n" % (m))
            self.fhc.write("  QUICKDouble x_%d_4; // Py, Px \n" % (m))
            self.fhc.write("  QUICKDouble x_%d_5; // Py, Py \n" % (m))
            self.fhc.write("  QUICKDouble x_%d_6; // Py, Pz \n" % (m))
            self.fhc.write("  QUICKDouble x_%d_7; // Pz, Px \n" % (m))
            self.fhc.write("  QUICKDouble x_%d_8; // Pz, Py \n" % (m))
            self.fhc.write("  QUICKDouble x_%d_9; // Pz, Pz \n" % (m))
            self.fhc.write("  __device__ __inline__ PPint_%d(QUICKDouble* BB, QUICKDouble* CC, QUICKDouble* PP, QUICKDouble ABCD); \n" % (m))          
            self.fhc.write("}; \n")

            self.fhd.write("__device__ __inline__ PPauxint_%d::PPauxint_%d(QUICKDouble* BB, QUICKDouble* CC, QUICKDouble* PP, QUICKDouble ABCD){ \n\n" % (m, m))
            self.fhd.write("  PSint ps_%d(BB, CC, PP); // construct [p|s] for m=%d \n" % (m, m))
            self.fhd.write("  PSint ps_%d(BB, CC, PP); // construct [p|s] for m=%d \n" % (m+1, m+1))

            idx=1
            for i in range(0,3):
                for j in range(0,3):
                    self.fhd.write("  x_%d_%d = (PP[%d]-BB[%d]) * ps_%d.x_%d_%d - (PP[%d]-CC[%d]) * ps_%d.x_%d_%d \n" % (m, idx, j, j, m, m, i+1,\
                    j, j, m+1, m+1, i+1))

                    if i == j:
                        self.fhd.write("  x_%d_%d += 0.5/ABCD * (LOCVY(0, 0, %d) - LOCVY(0, 0, %d)) \n" % (m, idx, m, m+1))

                    idx += 1
            self.fhd.write("\n } \n")

def write_oei():

    # set files
    OEint.fhc = open('gpu_oei_classes.h','w')
    OEint.fhd = open('gpu_oei_definitions.h','w')

    # generate integrals for systems containing only s and p functions.
    # <s|s> is trivial and we wont auto generate it.

    # generate <p|s>
    ps=PSint(1)
    ps.gen_int() 

    # generate <s|p>
    sp=SPint(0)

    # generate <p|p>
    pp=PPint(0)
    pp.gen_int()

    OEint.fhc.close()
    OEint.fhd.close()

write_oei()
