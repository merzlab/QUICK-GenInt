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
        if self.m == 0:
            self.fhd.write("__device__ void PSint(QUICKDouble* BB, QUICKDouble* CC, QUICKDouble* PP){\n")
            for i in range(0,3):
                self.fhd.write("  LOCSTORE(store, 0, %d, STOREDIM, STOREDIM)=(PP[%d]-BB[%d]) * LOCVY(YVerticalTemp, 0, 0, %d, VDIM1, VDIM2, VDIM3)\
 - (PP[%d]-CC[%d]) * LOCVY(YVerticalTemp, 0, 0, %d, VDIM1, VDIM2, VDIM3)\n" % (i, i, i, 0, i, i, 0))
            self.fhd.write("}\n")

        else:
            # write code paths for auxilary integrals. Note that we use C++ classes here. 
            for i in range(0,self.m):
                self.fhc.write("PS auxilary integral m=%d \n" % (i))
                self.fhc.write("class PSauxint_%d { \n" % (i))
                self.fhc.write("public:")
                self.fhc.write("  QUICKDouble x_%d_1; \n" % (i))
                self.fhc.write("  QUICKDouble x_%d_2; \n" % (i))
                self.fhc.write("  QUICKDouble x_%d_3; \n" % (i))
                self.fhc.write("  __device__ __inline__ PSauxint_%d(QUICKDouble* BB, QUICKDouble* CC, QUICKDouble* PP); \n" % (i))
                self.fhc.write("}; \n")

                self.fhd.write("PS auxilary integral m=%d \n" % (i))
                self.fhd.write("  __device__ __inline__ PSauxint_%d::PSauxint_%d(QUICKDouble* BB, QUICKDouble* CC, QUICKDouble* PP){ \n" % (i, i))
                self.fhd.write("    QUICKDouble x_%d_1=(PP[0]-BB[0]) * LOCVY(YVerticalTemp, 0, 0, %d, VDIM1, VDIM2, VDIM3)\
 - (PP[0]-CC[0]) * LOCVY(YVerticalTemp, 0, 0, %d, VDIM1, VDIM2, VDIM3)\n" % (i, i, i))
                self.fhd.write("    QUICKDouble x_%d_2=(PP[1]-BB[1]) * LOCVY(YVerticalTemp, 0, 0, %d, VDIM1, VDIM2, VDIM3)\
 - (PP[1]-CC[1]) * LOCVY(YVerticalTemp, 0, 0, %d, VDIM1, VDIM2, VDIM3)\n" % (i, i, i))
                self.fhd.write("    QUICKDouble x_%d_3=(PP[2]-BB[2]) * LOCVY(YVerticalTemp, 0, 0, %d, VDIM1, VDIM2, VDIM3)\
 - (PP[2]-CC[2]) * LOCVY(YVerticalTemp, 0, 0, %d, VDIM1, VDIM2, VDIM3)\n" % (i, i, i))
                self.fhd.write("  } \n")                


def write_oei():

    # set files
    OEint.fhc = open('gpu_oei_classes.h','w')
    OEint.fhd = open('gpu_oei_definitions.h','w')

    ps=PSint(6)
    ps.gen_int()

    OEint.fhc.close()
    OEint.fhd.close()

#write_oei()
