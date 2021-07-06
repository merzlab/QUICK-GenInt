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

    # m ranges from 0 to a+b; where a and b are the angular momentum of i and j
    # of the integral being considered <i|j>. If m=0, the integral is a true integral
    # and m>0 results in auxliary integrals. 
    def __init__(self,m=0):
        self.m=m
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
        # write code for true integral (i.e. m=0). We will store the computed integral value directly inside store array.
        #if self.m == 0:
        #    self.fhd.write("__device__ void PSint(QUICKDouble* BB, QUICKDouble* CC, QUICKDouble* PP){\n")
        #    for i in range(0,3):
        #        self.fhd.write("  LOCSTORE(store, 0, %d, STOREDIM, STOREDIM)=(PP[%d]-BB[%d]) * LOCVY(0, 0, %d)\
        #- (PP[%d]-CC[%d]) * LOCVY(0, 0, %d)\n" % (i, i, i, 0, i, i, 0))
        #    self.fhd.write("}\n")

        #else:
            # write code paths for auxilary integrals. Note that we use C++ classes here. 
            for i in range(0,self.m+1):
                self.fhc.write("/* PS auxilary integral m=%d */\n" % (i))
                self.fhc.write("class PSint_%d { \n" % (i))
                self.fhc.write("public: \n")
                self.fhc.write("  QUICKDouble x_%d_1; // Px, S \n" % (i))
                self.fhc.write("  QUICKDouble x_%d_2; // Py, S \n" % (i))
                self.fhc.write("  QUICKDouble x_%d_3; // Pz, S \n" % (i))
                self.fhc.write("  __device__ __inline__ PSint_%d(QUICKDouble* BB, QUICKDouble* CC, QUICKDouble* PP); \n" % (i))
                self.fhc.write("}; \n")

                self.fhd.write("/* PS auxilary integral m=%d */ \n" % (i))
                self.fhd.write("  __device__ __inline__ PSint_%d::PSint_%d(QUICKDouble* BB, QUICKDouble* CC, QUICKDouble* PP){ \n" % (i, i))
                self.fhd.write("    QUICKDouble x_%d_1=(PP[0]-BB[0]) * LOCVY(0, 0, %d)\
 - (PP[0]-CC[0]) * LOCVY(0, 0, %d)\n" % (i, i, i))
                self.fhd.write("    QUICKDouble x_%d_2=(PP[1]-BB[1]) * LOCVY(0, 0, %d)\
 - (PP[1]-CC[1]) * LOCVY(0, 0, %d)\n" % (i, i, i))
                self.fhd.write("    QUICKDouble x_%d_3=(PP[2]-BB[2]) * LOCVY(0, 0, %d)\
 - (PP[2]-CC[2]) * LOCVY(0, 0, %d)\n" % (i, i, i))
                self.fhd.write("  } \n")                


# <s|p> class, a subclass of OEint
class SPint(OEint):
    def gen_int(self):
        # write code for true integral (i.e. m=0). We will store the computed integral value directly inside store array.
        #if self.m == 0:
        #    self.fhd.write("__device__ void SPint(QUICKDouble* BB, QUICKDouble* CC, QUICKDouble* PP){\n")
        #    for i in range(0,3):
        #        self.fhd.write("  LOCSTORE(store, %d, 0, STOREDIM, STOREDIM)=(PP[%d]-BB[%d]) * LOCVY(0, 0, %d)\
        #- (PP[%d]-CC[%d]) * LOCVY(0, 0, %d)\n" % (i, i, i, 0, i, i, 0))
        #    self.fhd.write("}\n")

        #else:
            # write code paths for auxilary integrals. Note that we use C++ classes here. 
            for i in range(0,self.m+1):
                self.fhc.write("/* SP auxilary integral m=%d */ \n" % (i))
                self.fhc.write("class SPint_%d { \n" % (i))
                self.fhc.write("public: \n")
                self.fhc.write("  QUICKDouble x_%d_1; // S, Px \n" % (i))
                self.fhc.write("  QUICKDouble x_%d_2; // S, Py \n" % (i))
                self.fhc.write("  QUICKDouble x_%d_3; // S, Pz \n" % (i))
                self.fhc.write("  __device__ __inline__ SPint_%d(QUICKDouble* BB, QUICKDouble* CC, QUICKDouble* PP); \n" % (i))
                self.fhc.write("}; \n")

                self.fhd.write("/* PS auxilary integral m=%d */ \n" % (i))
                self.fhd.write("  __device__ __inline__ SPint_%d::SPint_%d(QUICKDouble* BB, QUICKDouble* CC, QUICKDouble* PP){ \n" % (i, i))
                self.fhd.write("    QUICKDouble x_%d_1=(PP[0]-BB[0]) * LOCVY(0, 0, %d)\
 - (PP[0]-CC[0]) * LOCVY(0, 0, %d)\n" % (i, i, i))
                self.fhd.write("    QUICKDouble x_%d_2=(PP[1]-BB[1]) * LOCVY(0, 0, %d)\
 - (PP[1]-CC[1]) * LOCVY(0, 0, %d)\n" % (i, i, i))
                self.fhd.write("    QUICKDouble x_%d_3=(PP[2]-BB[2]) * LOCVY(0, 0, %d)\
 - (PP[2]-CC[2]) * LOCVY(0, 0, %d)\n" % (i, i, i))
                self.fhd.write("  } \n") 


# <p|p> class, subclass of OEint
class PPint(OEint):
    def gen_int(self):
        for i in range(0,self.m+1):
            if i == 0:
                self.fhc.write("/* PP true integral m=%d */ \n" % (i)) 
                self.fhd.write("/* PP true integral m=%d */ \n" % (i)) 
            else:
                self.fhc.write("/* PP auxilary integral m=%d */ \n" % (i))
                self.fhd.write("/* PP auxilary integral m=%d */ \n" % (i))          
        
            self.fhc.write("class PPint_%d { \n" % (i))
            self.fhc.write("public: \n")
            self.fhc.write("  QUICKDouble x_%d_1; // Px, Px \n" % (i))
            self.fhc.write("  QUICKDouble x_%d_2; // Px, Py \n" % (i))
            self.fhc.write("  QUICKDouble x_%d_3; // Px, Pz \n" % (i))
            self.fhc.write("  QUICKDouble x_%d_4; // Py, Px \n" % (i))
            self.fhc.write("  QUICKDouble x_%d_5; // Py, Py \n" % (i))
            self.fhc.write("  QUICKDouble x_%d_6; // Py, Pz \n" % (i))
            self.fhc.write("  QUICKDouble x_%d_7; // Pz, Px \n" % (i))
            self.fhc.write("  QUICKDouble x_%d_8; // Pz, Py \n" % (i))
            self.fhc.write("  QUICKDouble x_%d_9; // Pz, Pz \n" % (i))
            self.fhc.write("  __device__ __inline__ PPint_%d(QUICKDouble* BB, QUICKDouble* CC, QUICKDouble* PP, QUICKDouble ABCD); \n" % (i))          
            self.fhc.write("}; \n")

            self.fhd.write("  __device__ __inline__ PPauxint_%d::PPauxint_%d(QUICKDouble* BB, QUICKDouble* CC, QUICKDouble* PP, QUICKDouble ABCD){ \n" % (i, i))
            self.fhd.write("    PSint ps_%d(BB, CC, PP); // construct [p|s] for m=%d \n" % (i, i))
            self.fhd.write("    PSint ps_%d(BB, CC, PP); // construct [p|s] for m=%d \n" % (i+1, i+1))

            idx=1
            for idim in range(0,3):
                for jdim in range(0,3):
                    self.fhd.write("    x_%d_%d = (PP[%d]-BB[%d]) * ps_%d.x_%d_%d - (PP[%d]-CC[%d]) * ps_%d.x_%d_%d \n" % (i, idx, jdim, jdim, i, i, idim+1,\
                    jdim, jdim, i+1, i+1, idim+1))

                    if idim == jdim:
                        self.fhd.write("    x_%d_%d = x_%d_%d + 0.5/ABCD * (LOCVY(0, 0, %d) - LOCVY(0, 0, %d)) \n" % (i, idx, i, idx, i, i+1))

                    idx += 1
            self.fhd.write("  } \n")

def write_oei():

    # set files
    OEint.fhc = open('gpu_oei_classes.h','w')
    OEint.fhd = open('gpu_oei_definitions.h','w')

    pp=PPint(0)
    pp.gen_int()

    OEint.fhc.close()
    OEint.fhd.close()

write_oei()
