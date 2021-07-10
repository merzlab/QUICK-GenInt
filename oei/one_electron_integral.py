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
#! electron integrals. Note that we use vertical recurrence relations  !
#! algorithm developed by Obara and Saika. See J. Chem. Phys. 1986, 84,!
#! 3963−3974 paper for theoretical details.                            !
#!                                                                     !
#!---------------------------------------------------------------------!

import params

# parent class for all one electron integrals
class OEint:
    # set file handlers
    fhc = 0 # file handler for class declarations
    fhd = 0 # file handler for function implementations
    fhdr= 0 # file handler for driver implementations

    # max_m ranges from 0 to a+b; where a and b are the angular momentum of i and j
    # of the integral being considered [i|j]. If max_m=0, the integral is a true integral
    # and m>0 results in auxliary integrals. 
    def __init__(self,max_m=0):
        self.max_m=max_m
        if self.fhc == 0 or self.fhd == 0:
            print("Warning: file handlers in OEint class are not set. \n")
        
        # set orbital labels for improving readability 
        self.p_lbl=("Px", "Py", "Pz")
        self.d_lbl=("Dxy", "Dyz", "Dxz", "Dxx", "Dyy", "Dzz")

        # set labels for useful variables
        self.AA=("Ax", "Ay", "Az") # position of center A
        self.BB=("Bx", "By", "Bz") # position of center B
        self.CC=("Cx", "Cy", "Cz") # position of center C
        self.PP=("Px", "Py", "Pz") # quantity P in OS eqn 15
        self.PA=("PAx", "PAy", "PAz") # Px-Ax, Py-Ay, Pz-Az
        self.PB=("PBx", "PBy", "PBz") # Px-Bx, Py-By, Pz-Bz
        self.PC=("PCx", "PCy", "PCz") # Px-Cx, Py-Cy, Pz-Cz

    # generate source code for integrals
    def gen_int(self):
        pass 

    # save true integrals into store array
    def save_int(self):
        pass

# [s|s] class, a subclass of OEint
class SSint(OEint):

    # generate code to save computed [s|s] integral
    def save_int(self):
        self.fhdr.write("\n  /* SS integral, m=%d */ \n" % (0))
        self.fhdr.write("  if(I == 0 && J == 0){ \n")
        self.fhdr.write("    LOC2(store, 0, 0, STOREDIM, STOREDIM) = VY(0, 0, 0);\n")
        self.fhdr.write("  } \n")

# [p|s] class, a subclass of OEint
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
            self.fhc.write("  __device__ __inline__ PSint_%d(QUICKDouble PAx, QUICKDouble PAy, QUICKDouble PAz,\
                QUICKDouble PCx, QUICKDouble PCy, QUICKDouble PCz, QUICKDouble* YVerticalTemp); \n" % (m))
            self.fhc.write("}; \n")

            self.fhd.write("__device__ __inline__ PSint_%d::PSint_%d(QUICKDouble PAx, QUICKDouble PAy, QUICKDouble PAz,\
                QUICKDouble PCx, QUICKDouble PCy, QUICKDouble PCz, QUICKDouble* YVerticalTemp){ \n\n" % (m, m))
            for i in range(0,3):
                self.fhd.write("  x_%d_%d = %s * VY(0, 0, %d) - %s * VY(0, 0, %d);\n" % (i+1, 0, self.PA[i], m, self.PC[i], m))

            self.fhd.write("} \n\n")

    # generate code to save computed [p|s] integral
    def save_int(self):
        self.fhdr.write("\n  /* PS integral, m=%d */ \n" % (0))
        self.fhdr.write("  if(I == 1 && J == 0){ \n")
        self.fhdr.write("    PSint_0 ps(PAx, PAy, PAz, PCx, PCy, PCz, YVerticalTemp); \n")
        for i in range(0,3):                
            self.fhdr.write("    LOC2(store, %d, %d, STOREDIM, STOREDIM) = ps.x_%d_%d;\n" % (i+1, 0, i+1, 0))
        self.fhdr.write("  } \n")

# [s|p] class, a subclass of OEint
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
            self.fhc.write("  __device__ __inline__ SPint_%d(QUICKDouble* BB, QUICKDouble* CC, QUICKDouble* PP, QUICKDouble* YVerticalTemp); \n" % (m))
            self.fhc.write("}; \n")

            self.fhd.write("__device__ __inline__ SPint_%d::SPint_%d(QUICKDouble* BB, QUICKDouble* CC, QUICKDouble* PP, QUICKDouble* YVerticalTemp){ \n\n" % (m, m))

            for i in range(0,3):
                self.fhd.write("  x_%d_%d=(PP[%d]-BB[%d]) * VY(0, 0, %d) - (PP[%d]-CC[%d]) * VY(0, 0, %d);\n" % (0, i+1, i, i, m, i, i, m))

            self.fhd.write("} \n\n") 

    # generate code to save computed [s|p] integral
    def save_int(self):
        self.fhdr.write("\n  /* SP integral, m=%d */ \n" % (0))
        self.fhdr.write("  if(I == 0 && J == 1){ \n")
        self.fhdr.write("    SPint_0 sp(BB, CC, PP, YVerticalTemp); \n")
        for i in range(0,3):                
            self.fhdr.write("    LOC2(store, %d, %d, STOREDIM, STOREDIM) = sp.x_%d_%d;\n" % (0, i+1, 0, i+1))
        self.fhdr.write("  } \n")

# [p|p] class, subclass of OEint
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
            self.fhc.write("  __device__ __inline__ PPint_%d(QUICKDouble* BB, QUICKDouble* CC, QUICKDouble* PP, QUICKDouble ABCD, QUICKDouble* YVerticalTemp); \n" % (m))          
            self.fhc.write("}; \n")

            # write function definitions
            self.fhd.write("__device__ __inline__ PPint_%d::PPint_%d(QUICKDouble* BB, QUICKDouble* CC, QUICKDouble* PP, QUICKDouble ABCD, QUICKDouble* YVerticalTemp){ \n\n" % (m, m))
            self.fhd.write("  PSint_%d ps_%d(PAx, PAy, PAz, PCx, PCy, PCz, YVerticalTemp); // construct [p|s] for m=%d \n" % (m, m, m))
            self.fhd.write("  PSint_%d ps_%d(PAx, PAy, PAz, PCx, PCy, PCz, YVerticalTemp); // construct [p|s] for m=%d \n\n" % (m+1, m+1, m+1))

            idx=1
            for i in range(0,3):
                for j in range(0,3):
                    self.fhd.write("  x_%d_%d = (PP[%d]-BB[%d]) * ps_%d.x_%d_%d - (PP[%d]-CC[%d]) * ps_%d.x_%d_%d; \n" % (i+1, j+1, j, j, m, i+1, 0,\
                    j, j, m+1, i+1, 0))

                    if i == j:
                        self.fhd.write("  x_%d_%d += 0.5/ABCD * (VY(0, 0, %d) - VY(0, 0, %d)); \n" % (i+1, j+1, m, m+1))

                    idx += 1
            self.fhd.write("\n } \n")

    # generate code to save computed [p|p] integral
    def save_int(self):
        self.fhdr.write("\n  /* PP integral, m=%d */ \n" % (0))
        self.fhdr.write("  if(I == 1 && J == 1){ \n")
        self.fhdr.write("    PPint_0 pp(BB, CC, PP, ABCD, YVerticalTemp); \n")
        for i in range(0,3):
            for j in range(0,3):
                self.fhdr.write("    LOC2(store, %d, %d, STOREDIM, STOREDIM) = pp.x_%d_%d;\n" % (i+1, j+1, i+1, j+1))
        self.fhdr.write("  } \n")

# [d|s] class, subclass of OEint
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
            self.fhc.write("  __device__ __inline__ DSint_%d(QUICKDouble* AA, QUICKDouble* CC, QUICKDouble* PP, QUICKDouble ABCD, QUICKDouble* YVerticalTemp); \n" % (m))          
            self.fhc.write("}; \n")

            # write function definitions
            self.fhd.write("__device__ __inline__ DSint_%d::DSint_%d(QUICKDouble* AA, QUICKDouble* CC, QUICKDouble* PP, QUICKDouble ABCD, QUICKDouble* YVerticalTemp){ \n\n" % (m, m))
            self.fhd.write("  PSint_%d ps_%d(PAx, PAy, PAz, PCx, PCy, PCz, YVerticalTemp); // construct [p|s] for m=%d \n" % (m, m, m))
            self.fhd.write("  PSint_%d ps_%d(PAx, PAy, PAz, PCx, PCy, PCz, YVerticalTemp); // construct [p|s] for m=%d \n\n" % (m+1, m+1, m+1))

            for i in range(0,6):
                tmp_mcal=[params.Mcal[i+4][0], params.Mcal[i+4][1], params.Mcal[i+4][2]]
                for j in range(0,3):
                    if params.Mcal[i+4][j] != 0:
                        tmp_mcal[j] -= 1
                        tmp_i=params.trans[tmp_mcal[0]][tmp_mcal[1]][tmp_mcal[2]]
                        self.fhd.write("  x_%d_%d = (PP[%d]-AA[%d]) * ps_%d.x_%d_%d - (PP[%d]-CC[%d]) * ps_%d.x_%d_%d; \n" % (i+4, 0, j, j, m, tmp_i-1, 0,\
                        j, j, m+1, tmp_i-1, 0))

                        if params.Mcal[i+4][j] == 2:
                            self.fhd.write("  x_%d_%d += 0.5/ABCD * (VY(0, 0, %d) - VY(0, 0, %d)); \n" % (i+4, 0, m, m+1))

                        break
            self.fhd.write("\n } \n")

    # generate code to save computed [d|s] integral
    def save_int(self):
        self.fhdr.write("\n  /* DS integral, m=%d */ \n" % (0))
        self.fhdr.write("  if(I == 2 && J == 0){ \n")
        self.fhdr.write("    DSint_0 ds(AA, CC, PP, ABCD, YVerticalTemp); \n")
        for i in range(0,6):
            self.fhdr.write("    LOC2(store, %d, %d, STOREDIM, STOREDIM) = ds.x_%d_%d;\n" % (i+4, 0, i+4, 0))
        self.fhdr.write("  } \n")

# [s|d] class, subclass of OEint
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
            self.fhc.write("  __device__ __inline__ SDint_%d(QUICKDouble* BB, QUICKDouble* CC, QUICKDouble* PP, QUICKDouble ABCD, QUICKDouble* YVerticalTemp); \n" % (m))          
            self.fhc.write("}; \n")

            # write function definitions
            self.fhd.write("__device__ __inline__ SDint_%d::SDint_%d(QUICKDouble* BB, QUICKDouble* CC, QUICKDouble* PP, QUICKDouble ABCD, QUICKDouble* YVerticalTemp){ \n\n" % (m, m))
            self.fhd.write("  SPint_%d sp_%d(BB, CC, PP, YVerticalTemp); // construct [s|p] for m=%d \n" % (m, m, m))
            self.fhd.write("  SPint_%d sp_%d(BB, CC, PP, YVerticalTemp); // construct [s|p] for m=%d \n\n" % (m+1, m+1, m+1))

            for i in range(0,6):
                tmp_mcal=[params.Mcal[i+4][0], params.Mcal[i+4][1], params.Mcal[i+4][2]]
                for j in range(0,3):
                    if params.Mcal[i+4][j] != 0:
                        tmp_mcal[j] -= 1
                        tmp_i=params.trans[tmp_mcal[0]][tmp_mcal[1]][tmp_mcal[2]]
                        self.fhd.write("  x_%d_%d = (PP[%d]-BB[%d]) * sp_%d.x_%d_%d - (PP[%d]-CC[%d]) * sp_%d.x_%d_%d; \n" % (0, i+4, j, j, m, 0, tmp_i-1,\
                        j, j, m+1, 0, tmp_i-1))

                        if params.Mcal[i+4][j] == 2:
                            self.fhd.write("  x_%d_%d += 0.5/ABCD * (VY(0, 0, %d) - VY(0, 0, %d)); \n" % (0, i+4, m, m+1))

                        break
            self.fhd.write("\n } \n")

    # generate code to save computed [s|d] integral
    def save_int(self):
        self.fhdr.write("\n  /* SD integral, m=%d */ \n" % (0))
        self.fhdr.write("  if(I == 0 && J == 2){ \n")
        self.fhdr.write("    SDint_0 sd(AA, CC, PP, ABCD, YVerticalTemp); \n")
        for i in range(0,6):
            self.fhdr.write("    LOC2(store, %d, %d, STOREDIM, STOREDIM) = sd.x_%d_%d;\n" % (0, i+4, 0, i+4))
        self.fhdr.write("  } \n")

# [d|p] class, subclass of OEint
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
            self.fhc.write("  __device__ __inline__ DPint_%d(QUICKDouble* BB, QUICKDouble* CC, QUICKDouble* PP, QUICKDouble ABCD, QUICKDouble* YVerticalTemp); \n" % (m))          
            self.fhc.write("}; \n")

            # write function definitions
            self.fhd.write("__device__ __inline__ DPint_%d::DPint_%d(QUICKDouble* BB, QUICKDouble* CC, QUICKDouble* PP, QUICKDouble ABCD, QUICKDouble* YVerticalTemp){ \n\n" % (m, m))
            self.fhd.write("  PSint_%d ps_%d(PAx, PAy, PAz, PCx, PCy, PCz, YVerticalTemp); // construct [p|s] for m=%d \n" % (m, m, m))
            self.fhd.write("  PSint_%d ps_%d(PAx, PAy, PAz, PCx, PCy, PCz, YVerticalTemp); // construct [p|s] for m=%d \n" % (m+1, m+1, m+1))
            self.fhd.write("  DSint_%d ds_%d(BB, CC, PP, ABCD, YVerticalTemp); // construct [d|s] for m=%d \n" % (m, m, m))            
            self.fhd.write("  DSint_%d ds_%d(BB, CC, PP, ABCD, YVerticalTemp); // construct [d|s] for m=%d \n\n" % (m+1, m+1, m+1))

            for i in range(0,6):
                for j in range(0,3):
                    tmp_mcal=[params.Mcal[i+4][0], params.Mcal[i+4][1], params.Mcal[i+4][2]]
                    for k in range(0,3):
                        if params.Mcal[j+1][k] != 0:
                            self.fhd.write("  x_%d_%d = (PP[%d]-BB[%d]) * ds_%d.x_%d_%d - (PP[%d]-CC[%d]) * ds_%d.x_%d_%d; \n" % (i+4, j+1, k, k, m, i+4, 0,\
                            k, k, m+1, i+4, 0))


                            if tmp_mcal[k] != 0:
                                tmp_mcal[k] -= 1
                                tmp_i=params.trans[tmp_mcal[0]][tmp_mcal[1]][tmp_mcal[2]]
                                self.fhd.write("  x_%d_%d += 0.5/ABCD * %f * (ps_%d.x_%d_%d - ps_%d.x_%d_%d); \n" % (i+4, j+1, params.Mcal[i+4][k], m, tmp_i-1, 0, m+1, tmp_i-1, 0))

                            break
            self.fhd.write("\n } \n")

    # generate code to save computed [d|p] integral
    def save_int(self):
        self.fhdr.write("\n  /* DP integral, m=%d */ \n" % (0))
        self.fhdr.write("  if(I == 2 && J == 1){ \n")
        self.fhdr.write("    DPint_0 dp(BB, CC, PP, ABCD, YVerticalTemp); \n")
        for i in range(0,6):
            for j in range(0,3):
                self.fhdr.write("    LOC2(store, %d, %d, STOREDIM, STOREDIM) = dp.x_%d_%d;\n" % (i+4, j+1, i+4, j+1))
        self.fhdr.write("  } \n")


# [p|d] class, subclass of OEint
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
            self.fhc.write("  __device__ __inline__ PDint_%d(QUICKDouble* AA, QUICKDouble* CC, QUICKDouble* PP, QUICKDouble ABCD, QUICKDouble* YVerticalTemp); \n" % (m))          
            self.fhc.write("}; \n")

            # write function definitions
            self.fhd.write("__device__ __inline__ PDint_%d::PDint_%d(QUICKDouble* AA, QUICKDouble* CC, QUICKDouble* PP, QUICKDouble ABCD, QUICKDouble* YVerticalTemp){ \n\n" % (m, m))
            self.fhd.write("  SPint_%d sp_%d(AA, CC, PP, YVerticalTemp); // construct [s|p] for m=%d \n" % (m, m, m))
            self.fhd.write("  SPint_%d sp_%d(AA, CC, PP, YVerticalTemp); // construct [s|p] for m=%d \n" % (m+1, m+1, m+1))
            self.fhd.write("  SDint_%d sd_%d(AA, CC, PP, ABCD, YVerticalTemp); // construct [s|d] for m=%d \n" % (m, m, m))
            self.fhd.write("  SDint_%d sd_%d(AA, CC, PP, ABCD, YVerticalTemp); // construct [s|d] for m=%d \n\n" % (m+1, m+1, m+1))

            for i in range(0,6):
                for j in range(0,3):
                    tmp_mcal=[params.Mcal[i+4][0], params.Mcal[i+4][1], params.Mcal[i+4][2]]
                    for k in range(0,3):
                        if params.Mcal[j+1][k] != 0:
                            self.fhd.write("  x_%d_%d = (PP[%d]-AA[%d]) * sd_%d.x_%d_%d - (PP[%d]-CC[%d]) * sd_%d.x_%d_%d; \n" % (j+1, i+4, k, k, m, 0, i+4,\
                            k, k, m+1, 0, i+4))


                            if tmp_mcal[k] != 0:
                                tmp_mcal[k] -= 1
                                tmp_i=params.trans[tmp_mcal[0]][tmp_mcal[1]][tmp_mcal[2]]
                                self.fhd.write("  x_%d_%d += 0.5/ABCD * %f * (sp_%d.x_%d_%d - sp_%d.x_%d_%d); \n" % (j+1, i+4, params.Mcal[i+4][k], m, 0, tmp_i-1, m+1, 0, tmp_i-1))
                            break
            self.fhd.write("\n } \n")

    # generate code to save computed [p|d] integral
    def save_int(self):
        self.fhdr.write("\n  /* PD integral, m=%d */ \n" % (0))
        self.fhdr.write("  if(I == 1 && J == 2){ \n")
        self.fhdr.write("    PDint_0 pd(AA, CC, PP, ABCD, YVerticalTemp); \n")
        for i in range(0,6):
            for j in range(0,3):
                self.fhdr.write("    LOC2(store, %d, %d, STOREDIM, STOREDIM) = pd.x_%d_%d;\n" % (j+1, i+4, j+1, i+4))
        self.fhdr.write("  } \n")


# [d|d] class, subclass of OEint
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
            self.fhc.write("  __device__ __inline__ DDint_%d(QUICKDouble* BB, QUICKDouble* CC, QUICKDouble* PP, QUICKDouble ABCD, QUICKDouble* YVerticalTemp); \n" % (m))          
            self.fhc.write("}; \n")

            # write function definitions
            self.fhd.write("__device__ __inline__ DDint_%d::DDint_%d(QUICKDouble* BB, QUICKDouble* CC, QUICKDouble* PP, QUICKDouble ABCD, QUICKDouble* YVerticalTemp){ \n\n" % (m, m))
            self.fhd.write("  PPint_%d pp_%d(BB, CC, PP, ABCD, YVerticalTemp); // construct [p|p] for m=%d \n" % (m, m, m))
            self.fhd.write("  DSint_%d ds_%d(BB, CC, PP, ABCD, YVerticalTemp); // construct [d|s] for m=%d \n" % (m, m, m))
            self.fhd.write("  DPint_%d dp_%d(BB, CC, PP, ABCD, YVerticalTemp); // construct [d|p] for m=%d \n" % (m, m, m))            
            self.fhd.write("  PPint_%d pp_%d(BB, CC, PP, ABCD, YVerticalTemp); // construct [p|p] for m=%d \n" % (m+1, m+1, m+1))
            self.fhd.write("  DSint_%d ds_%d(BB, CC, PP, ABCD, YVerticalTemp); // construct [d|s] for m=%d \n" % (m+1, m+1, m+1))
            self.fhd.write("  DPint_%d dp_%d(BB, CC, PP, ABCD, YVerticalTemp); // construct [d|p] for m=%d \n\n" % (m+1, m+1, m+1))             

            for i in range(0,6):
                for j in range(0,6):
                    tmp_mcal1=[params.Mcal[i+4][0], params.Mcal[i+4][1], params.Mcal[i+4][2]]
                    tmp_mcal2=[params.Mcal[j+4][0], params.Mcal[j+4][1], params.Mcal[j+4][2]]
                    for k in range(0,3):
                        if params.Mcal[j+4][k] != 0:
                            tmp_mcal2[k] -= 1
                            tmp_j=params.trans[tmp_mcal2[0]][tmp_mcal2[1]][tmp_mcal2[2]]

                            iclass_obj="dp"
                            self.fhd.write("  x_%d_%d = (PP[%d]-BB[%d]) * %s_%d.x_%d_%d - (PP[%d]-CC[%d]) * %s_%d.x_%d_%d; \n" % (j+4, i+4, k, k, iclass_obj, m, i+4, tmp_j-1,\
                            k, k, iclass_obj, m+1, i+4, tmp_j-1))

                            if params.Mcal[j+4][k] > 0:
                                iclass_obj="ds"
                                self.fhd.write("  x_%d_%d += 0.5/ABCD * %f * (%s_%d.x_%d_%d - %s_%d.x_%d_%d); \n" % (j+4, i+4, params.Mcal[j+4][k], iclass_obj, m, i+4, 0, \
                                iclass_obj, m+1, i+4, 0))

                            if tmp_mcal1[k] > 0:
                                tmp_mcal1[k] -= 1
                                tmp_i=params.trans[tmp_mcal1[0]][tmp_mcal1[1]][tmp_mcal1[2]]
                                iclass_obj="pp"
                                self.fhd.write("  x_%d_%d += 0.5/ABCD * %f * (%s_%d.x_%d_%d - %s_%d.x_%d_%d); \n" % (j+4, i+4, params.Mcal[i+4][k], iclass_obj, m, tmp_i-1, tmp_j-1,\
                                iclass_obj, m+1, tmp_i-1, tmp_j-1))
                            break
            self.fhd.write("\n } \n")

    # generate code to save computed [d|d] integral
    def save_int(self):
        self.fhdr.write("\n  /* DD integral, m=%d */ \n" % (0))
        self.fhdr.write("  if(I == 2 && J == 2){ \n")
        self.fhdr.write("    DDint_0 dd(BB, CC, PP, ABCD, YVerticalTemp); \n")
        for i in range(0,6):
            for j in range(0,6):
                self.fhdr.write("    LOC2(store, %d, %d, STOREDIM, STOREDIM) = dd.x_%d_%d;\n" % (j+4, i+4, j+4, i+4))
        self.fhdr.write("  } \n")

def write_oei():

    # set files
    OEint.fhc = open('gpu_oei_classes.h','w')
    OEint.fhd = open('gpu_oei_definitions.h','w')
    OEint.fhdr= open('gpu_oei_drivers.h','w')

    # generate integral classes for systems containing only s, p and d functions.
    # generate [s|s], this is trivial and we will directly save the integral value from the driver. 
    ss=SSint()

    # generate [p|s]
    ps=PSint(3)
    ps.gen_int() 

    # generate [s|p]
    sp=SPint(3)
    sp.gen_int()

    # generate [p|p]
    pp=PPint(1)
    pp.gen_int()

    # generate [d|s]
    ds=DSint(2)
    ds.gen_int()

    # generate [s|d]
    sd=SDint(2)
    sd.gen_int() 

    # generate [d|p]
    dp=DPint(1)
    dp.gen_int()

    # generate [p|d]
    pd=PDint(0)
    pd.gen_int()

    # generate [d|d]
    dd=DDint(0)
    dd.gen_int()

    # write driver to use classes and save computed primitive integrals
    OEint.fhdr.write("__device__ __inline__ void OEint_vertical(int I, int J, QUICKDouble* AA, QUICKDouble* BB, QUICKDouble* CC, QUICKDouble* PP, QUICKDouble ABCD, \
        QUICKDouble* store, QUICKDouble* YVerticalTemp){ \n")
    ss.save_int()
    ps.save_int()
    sp.save_int()
    pp.save_int()
    ds.save_int()
    sd.save_int()
    dp.save_int()
    pd.save_int()
    dd.save_int()
    OEint.fhdr.write("\n } \n")

    OEint.fhc.close()
    OEint.fhd.close()
    OEint.fhdr.close()

write_oei()
