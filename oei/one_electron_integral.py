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

                self.fhc.write("class PSint_%d{ \n" % (m))
                self.fhc.write("public: \n")

                # set labels for improving readability 
                lbl=["Px", "Py", "Pz"]

                # write class variables; convention being used is s=0, p=1-3, d=4-9, f=10-19, g=20-34
                for i in range(0,3):
                    self.fhc.write("  QUICKDouble x_%d_%d; // %s, %s \n" % (i+1, 0, lbl[i], "S"))

                # write class functions
                self.fhc.write("  __device__ __inline__ PSint_%d(QUICKDouble* AA, QUICKDouble* CC, QUICKDouble* PP); \n" % (m))
                self.fhc.write("}; \n")

                self.fhd.write("__device__ __inline__ PSint_%d::PSint_%d(QUICKDouble* AA, QUICKDouble* CC, QUICKDouble* PP){ \n\n" % (m, m))
                for i in range(0,3):
                    self.fhd.write("  x_%d_%d=(PP[%d]-AA[%d]) * LOCVY(0, 0, %d) - (PP[%d]-CC[%d]) * LOCVY(0, 0, %d)\n" % (i+1, 0, i, i, m, i, i, m))

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

                self.fhc.write("class SPint_%d{ \n" % (m))
                self.fhc.write("public: \n")

                # set labels for improving readability 
                lbl=["Px", "Py", "Pz"]

                # write class variables; convention being used is s=0, p=1-3, d=4-9, f=10-19, g=20-34
                for i in range(0,3):
                    self.fhc.write("  QUICKDouble x_%d_%d; // %s, %s \n" % (0, i+1, "S", lbl[i]))

                # write class functions
                self.fhc.write("  __device__ __inline__ SPint_%d(QUICKDouble* BB, QUICKDouble* CC, QUICKDouble* PP); \n" % (m))
                self.fhc.write("}; \n")

                self.fhd.write("__device__ __inline__ SPint_%d::SPint_%d(QUICKDouble* BB, QUICKDouble* CC, QUICKDouble* PP){ \n\n" % (m, m))

                for i in range(0,3):
                    self.fhd.write("  x_%d_%d=(PP[%d]-BB[%d]) * LOCVY(0, 0, %d) - (PP[%d]-CC[%d]) * LOCVY(0, 0, %d)\n" % (0, i+1, i, i, m, i, i, m))

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
        
            self.fhc.write("class PPint_%d{ \n" % (m))
            self.fhc.write("public: \n")

            # set labels for improving readability 
            lbl=["Px", "Py", "Pz"]

            # write class variables; convention being used is s=0, p=1-3, d=4-9, f=10-19, g=20-34
            for i in range(0,3):
                for j in range(0,3):
                    self.fhc.write("  QUICKDouble x_%d_%d; // %s, %s \n" % (i+1, j+1, lbl[i], lbl[j]))

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
                    self.fhd.write("  x_%d_%d = (PP[%d]-BB[%d]) * ps_%d.x_%d_%d - (PP[%d]-CC[%d]) * ps_%d.x_%d_%d \n" % (i+1, j+1, j, j, m, i+1, 0,\
                    j, j, m+1, i+1, 0))

                    if i == j:
                        self.fhd.write("  x_%d_%d += 0.5/ABCD * (LOCVY(0, 0, %d) - LOCVY(0, 0, %d)) \n" % (i+1, j+1, m, m+1))

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

            self.fhc.write("class DSint_%d{ \n" % (m))
            self.fhc.write("public: \n")

            # set labels for improving readability 
            lbl=["Dxy", "Dyz", "Dxz", "Dxx", "Dyy", "Dzz"]

            # write class variables; convention being used is s=0, p=1-3, d=4-9, f=10-19, g=20-34
            for i in range(0,6):
                self.fhc.write("  QUICKDouble x_%d_%d; // %s, %s \n" % (i+4, 0, lbl[i], "S"))

            # write class functions
            self.fhc.write("  __device__ __inline__ DSint_%d(QUICKDouble* AA, QUICKDouble* CC, QUICKDouble* PP, QUICKDouble ABCD); \n" % (m))          
            self.fhc.write("}; \n")

            # write function definitions
            self.fhd.write("__device__ __inline__ DSint_%d::DSint_%d(QUICKDouble* AA, QUICKDouble* CC, QUICKDouble* PP, QUICKDouble ABCD){ \n\n" % (m, m))
            self.fhd.write("  PSint ps_%d(AA, CC, PP); // construct [p|s] for m=%d \n" % (m, m))
            self.fhd.write("  PSint ps_%d(AA, CC, PP); // construct [p|s] for m=%d \n" % (m+1, m+1))

            for i in range(0,6):
                tmp_mcal=[params.Mcal[i+4][0], params.Mcal[i+4][1], params.Mcal[i+4][2]]
                for j in range(0,3):
                    if params.Mcal[i+4][j] != 0:
                        tmp_mcal[j] -= 1
                        tmp_i=params.trans[tmp_mcal[0]][tmp_mcal[1]][tmp_mcal[2]]
                        self.fhd.write("  x_%d_%d = (PP[%d]-AA[%d]) * ps_%d.x_%d_%d - (PP[%d]-CC[%d]) * ps_%d.x_%d_%d \n" % (i+4, 0, j, j, m, tmp_i-1, 0,\
                        j, j, m+1, tmp_i-1, 0))

                        if params.Mcal[i+4][j] == 2:
                            self.fhd.write("  x_%d_%d += 0.5/ABCD * (LOCVY(0, 0, %d) - LOCVY(0, 0, %d)) \n" % (i+4, 0, m, m+1))

                        break
            self.fhd.write("\n } \n")


# <s|d> class, subclass of OEint
class SDint(OEint):
    def gen_int(self):
        # write code paths for integrals. Note that we use C++ classes here.
        for m in range(0,self.max_m+1):
            if m == 0:
                self.fhc.write("\n/* SD true integral, m=%d */ \n" % (m)) 
                self.fhd.write("\n/* SD true integral, m=%d */ \n" % (m)) 
            else:
                self.fhc.write("\n/* SD auxilary integral, m=%d */ \n" % (m))
                self.fhd.write("\n/* SD auxilary integral, m=%d */ \n" % (m))              

            self.fhc.write("class SDint_%d{ \n" % (m))
            self.fhc.write("public: \n")

            # set labels for improving readability 
            lbl=["Dxy", "Dyz", "Dxz", "Dxx", "Dyy", "Dzz"]

            # write class variables; convention being used is s=0, p=1-3, d=4-9, f=10-19, g=20-34
            for i in range(0,6):
                self.fhc.write("  QUICKDouble x_%d_%d; // %s, %s \n" % (0, i+4, "S", lbl[i]))

            # write class functions
            self.fhc.write("  __device__ __inline__ SDint_%d(QUICKDouble* AA, QUICKDouble* CC, QUICKDouble* PP, QUICKDouble ABCD); \n" % (m))          
            self.fhc.write("}; \n")

            # write function definitions
            self.fhd.write("__device__ __inline__ SDint_%d::SDint_%d(QUICKDouble* AA, QUICKDouble* CC, QUICKDouble* PP, QUICKDouble ABCD){ \n\n" % (m, m))
            self.fhd.write("  SPint sp_%d(BB, CC, PP); // construct [s|p] for m=%d \n" % (m, m))
            self.fhd.write("  SPint sp_%d(BB, CC, PP); // construct [s|p] for m=%d \n" % (m+1, m+1))

            for i in range(0,6):
                tmp_mcal=[params.Mcal[i+4][0], params.Mcal[i+4][1], params.Mcal[i+4][2]]
                for j in range(0,3):
                    if params.Mcal[i+4][j] != 0:
                        tmp_mcal[j] -= 1
                        tmp_i=params.trans[tmp_mcal[0]][tmp_mcal[1]][tmp_mcal[2]]
                        self.fhd.write("  x_%d_%d = (PP[%d]-BB[%d]) * ps_%d.x_%d_%d - (PP[%d]-CC[%d]) * ps_%d.x_%d_%d \n" % (0, i+4, j, j, m, 0, tmp_i-1,\
                        j, j, m+1, 0, tmp_i-1))

                        if params.Mcal[i+4][j] == 2:
                            self.fhd.write("  x_%d_%d += 0.5/ABCD * (LOCVY(0, 0, %d) - LOCVY(0, 0, %d)) \n" % (0, i+4, m, m+1))

                        break
            self.fhd.write("\n } \n")


# <d|p> class, subclass of OEint
class DPint(OEint):
    def gen_int(self):
        # write code paths for integrals. Note that we use C++ classes here.
        for m in range(0,self.max_m+1):
            if m == 0:
                self.fhc.write("\n/* DP true integral, m=%d */ \n" % (m)) 
                self.fhd.write("\n/* DP true integral, m=%d */ \n" % (m)) 
            else:
                self.fhc.write("\n/* DP auxilary integral, m=%d */ \n" % (m))
                self.fhd.write("\n/* DP auxilary integral, m=%d */ \n" % (m))              

            self.fhc.write("class DPint_%d{ \n" % (m))
            self.fhc.write("public: \n")

            # set labels for improving readability 
            lbl=["Dxy", "Dyz", "Dxz", "Dxx", "Dyy", "Dzz"]
            lbl2=["Px", "Py", "Pz"]

            # write class variables; convention being used is s=0, p=1-3, d=4-9, f=10-19, g=20-34
            for i in range(0,6):
                for j in range(0,3):
                    self.fhc.write("  QUICKDouble x_%d_%d; // %s, %s \n" % (i+4, j+1, lbl[i], lbl2[j]))

            # write class functions
            self.fhc.write("  __device__ __inline__ DPint_%d(QUICKDouble* BB, QUICKDouble* CC, QUICKDouble* PP, QUICKDouble ABCD); \n" % (m))          
            self.fhc.write("}; \n")

            # write function definitions
            self.fhd.write("__device__ __inline__ DPint_%d::DPint_%d(QUICKDouble* BB, QUICKDouble* CC, QUICKDouble* PP, QUICKDouble ABCD){ \n\n" % (m, m))
            self.fhd.write("  DSint ds_%d(BB, CC, PP); // construct [d|s] for m=%d \n" % (m, m))
            self.fhd.write("  DSint ds_%d(BB, CC, PP); // construct [d|s] for m=%d \n" % (m+1, m+1))

            for i in range(0,6):
                for j in range(0,3):
                    tmp_mcal=[params.Mcal[i+4][0], params.Mcal[i+4][1], params.Mcal[i+4][2]]
                    for k in range(0,3):
                        if params.Mcal[j+1][k] != 0:
                            self.fhd.write("  x_%d_%d = (PP[%d]-BB[%d]) * ds_%d.x_%d_%d - (PP[%d]-CC[%d]) * ds_%d.x_%d_%d \n" % (i+4, j+1, k, k, m, i+4, 0,\
                            k, k, m+1, i+4, 0))


                            if tmp_mcal[k] != 0:
                                tmp_mcal[k] -= 1
                                tmp_i=params.trans[tmp_mcal[0]][tmp_mcal[1]][tmp_mcal[2]]
                                self.fhd.write("  x_%d_%d += 0.5/ABCD * %f * (ds_%d.x_%d_%d - ds_%d.x_%d_%d) \n" % (i+4, j+1, params.Mcal[i+4][k], m, tmp_i-1, 0, m+1, tmp_i-1, 0))

                            break
            self.fhd.write("\n } \n")



# <p|d> class, subclass of OEint
class PDint(OEint):
    def gen_int(self):
        # write code paths for integrals. Note that we use C++ classes here.
        for m in range(0,self.max_m+1):
            if m == 0:
                self.fhc.write("\n/* PD true integral, m=%d */ \n" % (m)) 
                self.fhd.write("\n/* PD true integral, m=%d */ \n" % (m)) 
            else:
                self.fhc.write("\n/* PD auxilary integral, m=%d */ \n" % (m))
                self.fhd.write("\n/* PD auxilary integral, m=%d */ \n" % (m))              

            self.fhc.write("class PDint_%d{ \n" % (m))
            self.fhc.write("public: \n")

            # set labels for improving readability 
            lbl=["Dxy", "Dyz", "Dxz", "Dxx", "Dyy", "Dzz"]
            lbl2=["Px", "Py", "Pz"]

            # write class variables; convention being used is s=0, p=1-3, d=4-9, f=10-19, g=20-34
            for i in range(0,6):
                for j in range(0,3):
                    self.fhc.write("  QUICKDouble x_%d_%d; // %s, %s \n" % (j+1, i+4, lbl2[j], lbl[i]))

            # write class functions
            self.fhc.write("  __device__ __inline__ PDint_%d(QUICKDouble* AA, QUICKDouble* CC, QUICKDouble* PP, QUICKDouble ABCD); \n" % (m))          
            self.fhc.write("}; \n")

            # write function definitions
            self.fhd.write("__device__ __inline__ PDint_%d::PDint_%d(QUICKDouble* AA, QUICKDouble* CC, QUICKDouble* PP, QUICKDouble ABCD){ \n\n" % (m, m))
            self.fhd.write("  SDint sd_%d(AA, CC, PP); // construct [s|d] for m=%d \n" % (m, m))
            self.fhd.write("  SDint sd_%d(AA, CC, PP); // construct [s|d] for m=%d \n" % (m+1, m+1))

            for i in range(0,6):
                for j in range(0,3):
                    tmp_mcal=[params.Mcal[i+4][0], params.Mcal[i+4][1], params.Mcal[i+4][2]]
                    for k in range(0,3):
                        if params.Mcal[j+1][k] != 0:
                            self.fhd.write("  x_%d_%d = (PP[%d]-AA[%d]) * sd_%d.x_%d_%d - (PP[%d]-CC[%d]) * sd_%d.x_%d_%d \n" % (j+1, i+4, k, k, m, 0, i+4,\
                            k, k, m+1, 0, i+4))


                            if tmp_mcal[k] != 0:
                                tmp_mcal[k] -= 1
                                tmp_i=params.trans[tmp_mcal[0]][tmp_mcal[1]][tmp_mcal[2]]
                                self.fhd.write("  x_%d_%d += 0.5/ABCD * %f * (sd_%d.x_%d_%d - sd_%d.x_%d_%d) \n" % (j+1, i+4, params.Mcal[i+4][k], m, 0, tmp_i-1, m+1, 0, tmp_i-1))
                            break
            self.fhd.write("\n } \n")


# <d|d> class, subclass of OEint
class DDint(OEint):
    def gen_int(self):
        # write code paths for integrals. Note that we use C++ classes here.
        for m in range(0,self.max_m+1):
            if m == 0:
                self.fhc.write("\n/* DD true integral, m=%d */ \n" % (m)) 
                self.fhd.write("\n/* DD true integral, m=%d */ \n" % (m)) 
            else:
                self.fhc.write("\n/* DD auxilary integral, m=%d */ \n" % (m))
                self.fhd.write("\n/* DD auxilary integral, m=%d */ \n" % (m))              

            self.fhc.write("class DDint_%d{ \n" % (m))
            self.fhc.write("public: \n")

            # set labels for improving readability 
            lbl=["Dxy", "Dyz", "Dxz", "Dxx", "Dyy", "Dzz"]

            # write class variables; convention being used is s=0, p=1-3, d=4-9, f=10-19, g=20-34
            for i in range(0,6):
                for j in range(0,6):
                    self.fhc.write("  QUICKDouble x_%d_%d; // %s, %s \n" % (i+4, j+4, lbl[i], lbl[j]))

            # write class functions
            self.fhc.write("  __device__ __inline__ DDint_%d(QUICKDouble* BB, QUICKDouble* CC, QUICKDouble* PP, QUICKDouble ABCD); \n" % (m))          
            self.fhc.write("}; \n")

            # write function definitions
            self.fhd.write("__device__ __inline__ DDint_%d::DDint_%d(QUICKDouble* BB, QUICKDouble* CC, QUICKDouble* PP, QUICKDouble ABCD){ \n\n" % (m, m))
            self.fhd.write("  PPint pp_%d(BB, CC, PP); // construct [p|p] for m=%d \n" % (m, m))
            self.fhd.write("  DSint ds_%d(BB, CC, PP); // construct [d|s] for m=%d \n" % (m, m))
            self.fhd.write("  DPint dp_%d(BB, CC, PP); // construct [d|p] for m=%d \n" % (m, m))            
            self.fhd.write("  PPint pp_%d(BB, CC, PP); // construct [p|p] for m=%d \n" % (m+1, m+1))
            self.fhd.write("  DSint ds_%d(BB, CC, PP); // construct [d|s] for m=%d \n" % (m+1, m+1))
            self.fhd.write("  DPint dp_%d(BB, CC, PP); // construct [d|p] for m=%d \n" % (m+1, m+1))             

            for i in range(0,6):
                for j in range(0,6):
                    tmp_mcal1=[params.Mcal[i+4][0], params.Mcal[i+4][1], params.Mcal[i+4][2]]
                    tmp_mcal2=[params.Mcal[j+4][0], params.Mcal[j+4][1], params.Mcal[j+4][2]]
                    for k in range(0,3):
                        if params.Mcal[j+4][k] != 0:
                            tmp_mcal2[k] -= 1
                            tmp_j=params.trans[tmp_mcal2[0]][tmp_mcal2[1]][tmp_mcal2[2]]

                            iclass_obj="dp"
                            self.fhd.write("  x_%d_%d = (PP[%d]-BB[%d]) * %s_%d.x_%d_%d - (PP[%d]-CC[%d]) * %s_%d.x_%d_%d \n" % (j+4, i+4, k, k, iclass_obj, m, i+4, tmp_j-1,\
                            k, k, iclass_obj, m+1, i+4, tmp_j-1))

                            if params.Mcal[j+4][k] > 0:
                                iclass_obj="ds"
                                self.fhd.write("  x_%d_%d += 0.5/ABCD * %f * (%s_%d.x_%d_%d - %s_%d.x_%d_%d) \n" % (j+4, i+4, params.Mcal[j+4][k], iclass_obj, m, i+4, 0, \
                                iclass_obj, m+1, i+4, 0))

                            if tmp_mcal1[k] > 0:
                                tmp_mcal1[k] -= 1
                                tmp_i=params.trans[tmp_mcal1[0]][tmp_mcal1[1]][tmp_mcal1[2]]
                                iclass_obj="pp"
                                self.fhd.write("  x_%d_%d += 0.5/ABCD * %f * (%s_%d.x_%d_%d - %s_%d.x_%d_%d) \n" % (j+4, i+4, params.Mcal[i+4][k], iclass_obj, m, tmp_i-1, tmp_j-1,\
                                iclass_obj, m+1, tmp_i-1, tmp_j-1))
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
    sp=SPint(2)
    sp.gen_int()

    # generate <p|p>
    pp=PPint(1)
    pp.gen_int()

    # generate <d|s>
    ds=DSint(1)
    ds.gen_int()

    # generate <s|d>
    sd=SDint(1)
    sd.gen_int() 

    # generate <d|p>
    dp=DPint(1)
    dp.gen_int()

    # generate <p|d>
    pd=PDint(0)
    pd.gen_int()

    # generate <d|d>
    dd=DDint(0)
    dd.gen_int()

    OEint.fhc.close()
    OEint.fhd.close()

write_oei()
