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
import file_handler

# parent class for all one electron integrals
class OEint:
    # set file handlers
    fhc = 0 # file handler for class declarations
    fhd = 0 # file handler for function implementations
    fha= 0  # file handler for integral assembler
    debug=1 # include debug info in generated code, 0=no, 1=yes 

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
        self.fha.write("\n  /* SS integral, m=%d */ \n" % (0))
        self.fha.write("  if(I == 0 && J == 0){ \n")
        self.fha.write("    LOC2(store, 0, 0, STOREDIM, STOREDIM) = VY(0, 0, 0);\n")

        # include print statements if debug option is on 
        if OEint.debug == 1:
            self.fha.write("\n#ifdef DEBUG_OEI \n")
            self.fha.write("    printf(\"II %%d JJ %%d %s store[%d,%d] = %%f \\n\", II, JJ, LOC2(store, 0, 0, STOREDIM, STOREDIM)); \n" % ( "SS", 0, 0))
            self.fha.write("#endif \n\n")

        self.fha.write("  } \n")

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

            # write class variables; convention being used is s=0, p=1-3, d=4-9, f=10-19, g=20-34
            for i in range(0,3):
                self.fhc.write("  QUICKDouble x_%d_%d; // %s, %s \n" % (i+1, 0, self.p_lbl[i], "S"))

            # write class functions
            self.fhc.write("  __device__ __inline__ PSint_%d(QUICKDouble PAx, QUICKDouble PAy, QUICKDouble PAz,\n\
                QUICKDouble PCx, QUICKDouble PCy, QUICKDouble PCz, QUICKDouble* YVerticalTemp); \n" % (m))
            self.fhc.write("}; \n")

            self.fhd.write("__device__ __inline__ PSint_%d::PSint_%d(QUICKDouble PAx, QUICKDouble PAy, QUICKDouble PAz,\n\
                QUICKDouble PCx, QUICKDouble PCy, QUICKDouble PCz, QUICKDouble* YVerticalTemp){ \n\n" % (m, m))
            for i in range(0,3):
                self.fhd.write("  x_%d_%d = %s * VY(0, 0, %d) - %s * VY(0, 0, %d);\n" % (i+1, 0, self.PA[i], m, self.PC[i], m+1))

            self.fhd.write("} \n\n")

    # generate code to save computed [p|s] integral
    def save_int(self):
        self.fha.write("\n  /* PS integral, m=%d */ \n" % (0))
        self.fha.write("  if(I == 1 && J == 0){ \n")
        self.fha.write("    PSint_0 ps(PAx, PAy, PAz, PCx, PCy, PCz, YVerticalTemp); \n")
        for i in range(0,3):                
            self.fha.write("    LOC2(store, %d, %d, STOREDIM, STOREDIM) = ps.x_%d_%d;\n" % (i+1, 0, i+1, 0))

        # include print statements if debug option is on 
        if OEint.debug == 1:
            self.fha.write("\n#ifdef DEBUG_OEI \n")
            for i in range(0,3):
                self.fha.write("    printf(\"II %%d JJ %%d %s store[%d,%d] = %%f \\n\", II, JJ, LOC2(store, %d, %d, STOREDIM, STOREDIM)); \n" % ( "PS", i+1, 0, i+1, 0))
            self.fha.write("#endif \n\n")

        self.fha.write("  } \n")

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

            # write class variables; convention being used is s=0, p=1-3, d=4-9, f=10-19, g=20-34
            for i in range(0,3):
                self.fhc.write("  QUICKDouble x_%d_%d; // %s, %s \n" % (0, i+1, "S", self.p_lbl[i]))

            # write class functions
            self.fhc.write("  __device__ __inline__ SPint_%d(QUICKDouble PBx, QUICKDouble PBy, QUICKDouble PBz,\n\
                QUICKDouble PCx, QUICKDouble PCy, QUICKDouble PCz, QUICKDouble* YVerticalTemp); \n" % (m))
            self.fhc.write("}; \n")

            self.fhd.write("__device__ __inline__ SPint_%d::SPint_%d(QUICKDouble PBx, QUICKDouble PBy, QUICKDouble PBz,\n\
                QUICKDouble PCx, QUICKDouble PCy, QUICKDouble PCz, QUICKDouble* YVerticalTemp){ \n\n" % (m, m))

            for i in range(0,3):
                self.fhd.write("  x_%d_%d = %s * VY(0, 0, %d) - %s * VY(0, 0, %d);\n" % (0, i+1, self.PB[i], m, self.PC[i], m+1))

            self.fhd.write("} \n\n") 

    # generate code to save computed [s|p] integral
    def save_int(self):
        self.fha.write("\n  /* SP integral, m=%d */ \n" % (0))
        self.fha.write("  if(I == 0 && J == 1){ \n")
        self.fha.write("    SPint_0 sp(PBx, PBy, PBz, PCx, PCy, PCz, YVerticalTemp); \n")
        for i in range(0,3):                
            self.fha.write("    LOC2(store, %d, %d, STOREDIM, STOREDIM) = sp.x_%d_%d;\n" % (0, i+1, 0, i+1))

        # include print statements if debug option is on    
        if OEint.debug == 1:
            self.fha.write("\n#ifdef DEBUG_OEI \n")
            for i in range(0,3):
                self.fha.write("    printf(\"II %%d JJ %%d %s store[%d,%d] = %%f \\n\", II, JJ, LOC2(store, %d, %d, STOREDIM, STOREDIM)); \n" % ( "SP", 0, i+1, 0, i+1))
            self.fha.write("#endif \n\n")

        self.fha.write("  } \n")

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

            # write class variables; convention being used is s=0, p=1-3, d=4-9, f=10-19, g=20-34
            for i in range(0,3):
                for j in range(0,3):
                    self.fhc.write("  QUICKDouble x_%d_%d; // %s, %s \n" % (i+1, j+1, self.p_lbl[i], self.p_lbl[j]))

            # write class functions
            self.fhc.write("  __device__ __inline__ PPint_%d(QUICKDouble PAx, QUICKDouble PAy, QUICKDouble PAz,\n\
                QUICKDouble PBx, QUICKDouble PBy, QUICKDouble PBz, QUICKDouble PCx, QUICKDouble PCy, QUICKDouble PCz,\n\
                QUICKDouble Zeta, QUICKDouble* YVerticalTemp); \n" % (m))          
            self.fhc.write("}; \n")

            # write function definitions
            self.fhd.write("__device__ __inline__ PPint_%d::PPint_%d(QUICKDouble PAx, QUICKDouble PAy, QUICKDouble PAz,\n\
                QUICKDouble PBx, QUICKDouble PBy, QUICKDouble PBz, QUICKDouble PCx, QUICKDouble PCy, QUICKDouble PCz,\n\
                QUICKDouble Zeta, QUICKDouble* YVerticalTemp){ \n\n" % (m, m))
            self.fhd.write("  PSint_%d ps_%d(PAx, PAy, PAz, PCx, PCy, PCz, YVerticalTemp); // construct [p|s] for m=%d \n" % (m, m, m))
            self.fhd.write("  PSint_%d ps_%d(PAx, PAy, PAz, PCx, PCy, PCz, YVerticalTemp); // construct [p|s] for m=%d \n\n" % (m+1, m+1, m+1))

            idx=1
            for i in range(0,3):
                for j in range(0,3):
                    self.fhd.write("  x_%d_%d = %s * ps_%d.x_%d_%d - %s * ps_%d.x_%d_%d; \n" % (i+1, j+1, self.PB[j], m, i+1, 0,\
                    self.PC[j], m+1, i+1, 0))

                    if i == j:
                        self.fhd.write("  x_%d_%d += 0.5/Zeta * (VY(0, 0, %d) - VY(0, 0, %d)); \n" % (i+1, j+1, m, m+1))

                    idx += 1
            self.fhd.write("\n } \n")

    # generate code to save computed [p|p] integral
    def save_int(self):
        self.fha.write("\n  /* PP integral, m=%d */ \n" % (0))
        self.fha.write("  if(I == 1 && J == 1){ \n")
        self.fha.write("    PPint_0 pp(PAx, PAy, PAz, PBx, PBy, PBz, PCx, PCy, PCz, Zeta, YVerticalTemp); \n")
        for i in range(0,3):
            for j in range(0,3):
                self.fha.write("    LOC2(store, %d, %d, STOREDIM, STOREDIM) = pp.x_%d_%d;\n" % (i+1, j+1, i+1, j+1))

        # include print statements if debug option is on    
        if OEint.debug == 1:
            self.fha.write("\n#ifdef DEBUG_OEI \n")
            for i in range(0,3):
                for j in range(0,3):
                    self.fha.write("    printf(\"II %%d JJ %%d %s store[%d,%d] = %%f \\n\", II, JJ, LOC2(store, %d, %d, STOREDIM, STOREDIM)); \n" % ( "PP", i+1, j+1, i+1, j+1))
            self.fha.write("#endif \n\n")

        self.fha.write("  } \n")

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

            # write class variables; convention being used is s=0, p=1-3, d=4-9, f=10-19, g=20-34
            for i in range(0,6):
                self.fhc.write("  QUICKDouble x_%d_%d; // %s, %s \n" % (i+4, 0, self.d_lbl[i], "S"))

            # write class functions
            self.fhc.write("  __device__ __inline__ DSint_%d(QUICKDouble PAx, QUICKDouble PAy, QUICKDouble PAz,\n\
                QUICKDouble PCx, QUICKDouble PCy, QUICKDouble PCz, QUICKDouble Zeta, QUICKDouble* YVerticalTemp); \n" % (m))          
            self.fhc.write("}; \n")

            # write function definitions
            self.fhd.write("__device__ __inline__ DSint_%d::DSint_%d(QUICKDouble PAx, QUICKDouble PAy, QUICKDouble PAz,\n\
                QUICKDouble PCx, QUICKDouble PCy, QUICKDouble PCz, QUICKDouble Zeta, QUICKDouble* YVerticalTemp){ \n\n" % (m, m))
            self.fhd.write("  PSint_%d ps_%d(PAx, PAy, PAz, PCx, PCy, PCz, YVerticalTemp); // construct [p|s] for m=%d \n" % (m, m, m))
            self.fhd.write("  PSint_%d ps_%d(PAx, PAy, PAz, PCx, PCy, PCz, YVerticalTemp); // construct [p|s] for m=%d \n\n" % (m+1, m+1, m+1))

            for i in range(0,6):
                tmp_mcal=[params.Mcal[i+4][0], params.Mcal[i+4][1], params.Mcal[i+4][2]]
                for j in range(0,3):
                    if params.Mcal[i+4][j] != 0:
                        tmp_mcal[j] -= 1
                        tmp_i=params.trans[tmp_mcal[0]][tmp_mcal[1]][tmp_mcal[2]]
                        self.fhd.write("  x_%d_%d = %s * ps_%d.x_%d_%d - %s * ps_%d.x_%d_%d; \n" % (i+4, 0, self.PA[j], m, tmp_i-1, 0,\
                        self.PC[j], m+1, tmp_i-1, 0))

                        if params.Mcal[i+4][j] == 2:
                            self.fhd.write("  x_%d_%d += 0.5/Zeta * (VY(0, 0, %d) - VY(0, 0, %d)); \n" % (i+4, 0, m, m+1))

                        break
            self.fhd.write("\n } \n")

    # generate code to save computed [d|s] integral
    def save_int(self):
        self.fha.write("\n  /* DS integral, m=%d */ \n" % (0))
        self.fha.write("  if(I == 2 && J == 0){ \n")
        self.fha.write("    DSint_0 ds(PAx, PAy, PAz, PCx, PCy, PCz, Zeta, YVerticalTemp); \n")
        for i in range(0,6):
            self.fha.write("    LOC2(store, %d, %d, STOREDIM, STOREDIM) = ds.x_%d_%d;\n" % (i+4, 0, i+4, 0))

        # include print statements if debug option is on    
        if OEint.debug == 1:
            self.fha.write("\n#ifdef DEBUG_OEI \n")
            for i in range(0,6):
                self.fha.write("    printf(\"II %%d JJ %%d %s store[%d,%d] = %%f \\n\", II, JJ, LOC2(store, %d, %d, STOREDIM, STOREDIM)); \n" % ( "DS", i+4, 0, i+4, 0))
            self.fha.write("#endif \n\n")

        self.fha.write("  } \n")

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

            # write class variables; convention being used is s=0, p=1-3, d=4-9, f=10-19, g=20-34
            for i in range(0,6):
                self.fhc.write("  QUICKDouble x_%d_%d; // %s, %s \n" % (0, i+4, "S", self.d_lbl[i]))

            # write class functions
            self.fhc.write("  __device__ __inline__ SDint_%d(QUICKDouble PBx, QUICKDouble PBy, QUICKDouble PBz,\n\
                QUICKDouble PCx, QUICKDouble PCy, QUICKDouble PCz, QUICKDouble Zeta, QUICKDouble* YVerticalTemp); \n" % (m))          
            self.fhc.write("}; \n")

            # write function definitions
            self.fhd.write("__device__ __inline__ SDint_%d::SDint_%d(QUICKDouble PBx, QUICKDouble PBy, QUICKDouble PBz,\n\
                QUICKDouble PCx, QUICKDouble PCy, QUICKDouble PCz, QUICKDouble Zeta, QUICKDouble* YVerticalTemp){ \n\n" % (m, m))
            self.fhd.write("  SPint_%d sp_%d(PBx, PBy, PBz, PCx, PCy, PCz, YVerticalTemp); // construct [s|p] for m=%d \n" % (m, m, m))
            self.fhd.write("  SPint_%d sp_%d(PBx, PBy, PBz, PCx, PCy, PCz, YVerticalTemp); // construct [s|p] for m=%d \n\n" % (m+1, m+1, m+1))

            for i in range(0,6):
                tmp_mcal=[params.Mcal[i+4][0], params.Mcal[i+4][1], params.Mcal[i+4][2]]
                for j in range(0,3):
                    if params.Mcal[i+4][j] != 0:
                        tmp_mcal[j] -= 1
                        tmp_i=params.trans[tmp_mcal[0]][tmp_mcal[1]][tmp_mcal[2]]
                        self.fhd.write("  x_%d_%d = %s * sp_%d.x_%d_%d - %s * sp_%d.x_%d_%d; \n" % (0, i+4, self.PB[j], m, 0, tmp_i-1,\
                        self.PC[j], m+1, 0, tmp_i-1))

                        if params.Mcal[i+4][j] == 2:
                            self.fhd.write("  x_%d_%d += 0.5/Zeta * (VY(0, 0, %d) - VY(0, 0, %d)); \n" % (0, i+4, m, m+1))

                        break
            self.fhd.write("\n } \n")

    # generate code to save computed [s|d] integral
    def save_int(self):
        self.fha.write("\n  /* SD integral, m=%d */ \n" % (0))
        self.fha.write("  if(I == 0 && J == 2){ \n")
        self.fha.write("    SDint_0 sd(PAx, PAy, PAz, PCx, PCy, PCz, Zeta, YVerticalTemp); \n")
        for i in range(0,6):
            self.fha.write("    LOC2(store, %d, %d, STOREDIM, STOREDIM) = sd.x_%d_%d;\n" % (0, i+4, 0, i+4))

        # include print statements if debug option is on    
        if OEint.debug == 1:
            self.fha.write("\n#ifdef DEBUG_OEI \n")
            for i in range(0,6):
                self.fha.write("    printf(\"II %%d JJ %%d %s store[%d,%d] = %%f \\n\", II, JJ, LOC2(store, %d, %d, STOREDIM, STOREDIM)); \n" % ( "SD", 0, i+4, 0, i+4))
            self.fha.write("#endif \n\n")

        self.fha.write("  } \n")

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

            # write class variables; convention being used is s=0, p=1-3, d=4-9, f=10-19, g=20-34
            for i in range(0,6):
                for j in range(0,3):
                    self.fhc.write("  QUICKDouble x_%d_%d; // %s, %s \n" % (i+4, j+1, self.d_lbl[i], self.p_lbl[j]))

            # write class functions
            self.fhc.write("  __device__ __inline__ DPint_%d(QUICKDouble PAx, QUICKDouble PAy, QUICKDouble PAz,\n\
                QUICKDouble PBx, QUICKDouble PBy, QUICKDouble PBz, QUICKDouble PCx, QUICKDouble PCy, QUICKDouble PCz,\n\
                QUICKDouble Zeta, QUICKDouble* YVerticalTemp); \n" % (m))          
            self.fhc.write("}; \n")

            # write function definitions
            self.fhd.write("__device__ __inline__ DPint_%d::DPint_%d(QUICKDouble PAx, QUICKDouble PAy, QUICKDouble PAz,\n\
                QUICKDouble PBx, QUICKDouble PBy, QUICKDouble PBz, QUICKDouble PCx, QUICKDouble PCy, QUICKDouble PCz,\n\
                QUICKDouble Zeta, QUICKDouble* YVerticalTemp){ \n\n" % (m, m))
            self.fhd.write("  PSint_%d ps_%d(PAx, PAy, PAz, PCx, PCy, PCz, YVerticalTemp); // construct [p|s] for m=%d \n" % (m, m, m))
            self.fhd.write("  PSint_%d ps_%d(PAx, PAy, PAz, PCx, PCy, PCz, YVerticalTemp); // construct [p|s] for m=%d \n" % (m+1, m+1, m+1))
            self.fhd.write("  DSint_%d ds_%d(PAx, PAy, PAz, PCx, PCy, PCz, Zeta, YVerticalTemp); // construct [d|s] for m=%d \n" % (m, m, m))            
            self.fhd.write("  DSint_%d ds_%d(PAx, PAy, PAz, PCx, PCy, PCz, Zeta, YVerticalTemp); // construct [d|s] for m=%d \n\n" % (m+1, m+1, m+1))

            for i in range(0,6):
                for j in range(0,3):
                    tmp_mcal=[params.Mcal[i+4][0], params.Mcal[i+4][1], params.Mcal[i+4][2]]
                    for k in range(0,3):
                        if params.Mcal[j+1][k] != 0:
                            self.fhd.write("  x_%d_%d = %s * ds_%d.x_%d_%d - %s * ds_%d.x_%d_%d; \n" % (i+4, j+1, self.PB[k], m, i+4, 0,\
                            self.PC[k], m+1, i+4, 0))


                            if tmp_mcal[k] != 0:
                                tmp_mcal[k] -= 1
                                tmp_i=params.trans[tmp_mcal[0]][tmp_mcal[1]][tmp_mcal[2]]
                                self.fhd.write("  x_%d_%d += 0.5/Zeta * %f * (ps_%d.x_%d_%d - ps_%d.x_%d_%d); \n" % (i+4, j+1, params.Mcal[i+4][k], m, tmp_i-1, 0, m+1, tmp_i-1, 0))

                            break
            self.fhd.write("\n } \n")

    # generate code to save computed [d|p] integral
    def save_int(self):
        self.fha.write("\n  /* DP integral, m=%d */ \n" % (0))
        self.fha.write("  if(I == 2 && J == 1){ \n")
        self.fha.write("    DPint_0 dp(PAx, PAy, PAz, PBx, PBy, PBz, PCx, PCy, PCz, Zeta, YVerticalTemp); \n")
        for i in range(0,6):
            for j in range(0,3):
                self.fha.write("    LOC2(store, %d, %d, STOREDIM, STOREDIM) = dp.x_%d_%d;\n" % (i+4, j+1, i+4, j+1))

        # include print statements if debug option is on    
        if OEint.debug == 1:
            self.fha.write("\n#ifdef DEBUG_OEI \n")
            for i in range(0,6):
                for j in range(0,3):
                    self.fha.write("    printf(\"II %%d JJ %%d %s store[%d,%d] = %%f \\n\", II, JJ, LOC2(store, %d, %d, STOREDIM, STOREDIM)); \n" % ( "DP", i+4, j+1, i+4, j+1))
            self.fha.write("#endif \n\n")

        self.fha.write("  } \n")


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

            # write class variables; convention being used is s=0, p=1-3, d=4-9, f=10-19, g=20-34
            for i in range(0,6):
                for j in range(0,3):
                    self.fhc.write("  QUICKDouble x_%d_%d; // %s, %s \n" % (j+1, i+4, self.p_lbl[j], self.d_lbl[i]))

            # write class functions
            self.fhc.write("  __device__ __inline__ PDint_%d(QUICKDouble PAx, QUICKDouble PAy, QUICKDouble PAz,\n\
                QUICKDouble PBx, QUICKDouble PBy, QUICKDouble PBz, QUICKDouble PCx, QUICKDouble PCy, QUICKDouble PCz,\n\
                QUICKDouble Zeta, QUICKDouble* YVerticalTemp); \n" % (m))          
            self.fhc.write("}; \n")

            # write function definitions
            self.fhd.write("__device__ __inline__ PDint_%d::PDint_%d(QUICKDouble PAx, QUICKDouble PAy, QUICKDouble PAz,\n\
                QUICKDouble PBx, QUICKDouble PBy, QUICKDouble PBz, QUICKDouble PCx, QUICKDouble PCy, QUICKDouble PCz,\n\
                QUICKDouble Zeta, QUICKDouble* YVerticalTemp){ \n\n" % (m, m))
            self.fhd.write("  SPint_%d sp_%d(PBx, PBy, PBz, PCx, PCy, PCz, YVerticalTemp); // construct [s|p] for m=%d \n" % (m, m, m))
            self.fhd.write("  SPint_%d sp_%d(PBx, PBy, PBz, PCx, PCy, PCz, YVerticalTemp); // construct [s|p] for m=%d \n" % (m+1, m+1, m+1))
            self.fhd.write("  SDint_%d sd_%d(PBx, PBy, PBz, PCx, PCy, PCz, Zeta, YVerticalTemp); // construct [s|d] for m=%d \n" % (m, m, m))
            self.fhd.write("  SDint_%d sd_%d(PBx, PBy, PBz, PCx, PCy, PCz, Zeta, YVerticalTemp); // construct [s|d] for m=%d \n\n" % (m+1, m+1, m+1))

            for i in range(0,6):
                for j in range(0,3):
                    tmp_mcal=[params.Mcal[i+4][0], params.Mcal[i+4][1], params.Mcal[i+4][2]]
                    
                    for k in range(0,3):
                        #self.fhd.write("a(i,j) %d %d %d %d %d\n" % (tmp_mcal[0], tmp_mcal[1], tmp_mcal[2], params.Mcal[j+1][k], tmp_mcal[k]))
                        if params.Mcal[j+1][k] != 0:
                            self.fhd.write("  x_%d_%d = %s * sd_%d.x_%d_%d - %s * sd_%d.x_%d_%d; \n" % (j+1, i+4, self.PA[k], m, 0, i+4,\
                            self.PC[k], m+1, 0, i+4))

                            #self.fhd.write("  printf(\" x_%d_%d: %s= %%f, sd_%d.x_%d_%d= %%f, %s= %%f, sd_%d.x_%d_%d= %%f \\n\", %s, sd_%d.x_%d_%d, %s, sd_%d.x_%d_%d); \n" \
                            #% ( j+1, i+4, self.PA[k], m, 0, i+4, self.PC[k], m+1, 0, i+4, self.PA[k], m, 0, i+4, self.PC[k], m+1, 0, i+4))

                            if tmp_mcal[k] != 0:                                
                                tmp_mcal[k] -= 1
                                tmp_i=params.trans[tmp_mcal[0]][tmp_mcal[1]][tmp_mcal[2]]
                                self.fhd.write("  x_%d_%d += 0.5/Zeta * %f * (sp_%d.x_%d_%d - sp_%d.x_%d_%d); \n" % (j+1, i+4, params.Mcal[i+4][k], m, 0, tmp_i-1, m+1, 0, tmp_i-1))

                            # added only for debugging
                            if tmp_mcal[k]+1 != 0:
                                self.fhd.write("  printf(\" x_%d_%d: %s= %%f, sd_%d.x_%d_%d= %%f, %s= %%f, sd_%d.x_%d_%d= %%f, Zeta= %%f, sp_%d.x_%d_%d= %%f, sp_%d.x_%d_%d= %%f\\n\", %s, sd_%d.x_%d_%d, %s, sd_%d.x_%d_%d, Zeta, sp_%d.x_%d_%d, sp_%d.x_%d_%d); \n" \
                                % ( j+1, i+4, self.PA[k], m, 0, i+4, self.PC[k], m+1, 0, i+4, m, 0, tmp_i-1, m+1, 0, tmp_i-1, self.PA[k], m, 0, i+4, self.PC[k], m+1, 0, i+4, m, 0, tmp_i-1, m+1, 0, tmp_i-1))
                            else:
                                self.fhd.write("  printf(\" x_%d_%d: %s= %%f, sd_%d.x_%d_%d= %%f, %s= %%f, sd_%d.x_%d_%d= %%f \\n\", %s, sd_%d.x_%d_%d, %s, sd_%d.x_%d_%d); \n" \
                                % ( j+1, i+4, self.PA[k], m, 0, i+4, self.PC[k], m+1, 0, i+4, self.PA[k], m, 0, i+4, self.PC[k], m+1, 0, i+4))                                

                            break
            self.fhd.write("\n } \n")

    # generate code to save computed [p|d] integral
    def save_int(self):
        self.fha.write("\n  /* PD integral, m=%d */ \n" % (0))
        self.fha.write("  if(I == 1 && J == 2){ \n")
        self.fha.write("    PDint_0 pd(PAx, PAy, PAz, PBx, PBy, PBz, PCx, PCy, PCz, Zeta, YVerticalTemp); \n")
        for i in range(0,6):
            for j in range(0,3):
                self.fha.write("    LOC2(store, %d, %d, STOREDIM, STOREDIM) = pd.x_%d_%d;\n" % (j+1, i+4, j+1, i+4))

        # include print statements if debug option is on    
        if OEint.debug == 1:
            self.fha.write("\n#ifdef DEBUG_OEI \n")
            for i in range(0,6):
                for j in range(0,3):
                    self.fha.write("    printf(\"II %%d JJ %%d %s store[%d,%d] = %%f \\n\", II, JJ, LOC2(store, %d, %d, STOREDIM, STOREDIM)); \n" % ( "PD", j+1, i+4, j+1, i+4))
            self.fha.write("#endif \n\n")

        self.fha.write("  } \n")


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

            # write class variables; convention being used is s=0, p=1-3, d=4-9, f=10-19, g=20-34
            for i in range(0,6):
                for j in range(0,6):
                    self.fhc.write("  QUICKDouble x_%d_%d; // %s, %s \n" % (i+4, j+4, self.d_lbl[i], self.d_lbl[j]))

            # write class functions
            self.fhc.write("  __device__ __inline__ DDint_%d(QUICKDouble PAx, QUICKDouble PAy, QUICKDouble PAz,\n\
                QUICKDouble PBx, QUICKDouble PBy, QUICKDouble PBz, QUICKDouble PCx, QUICKDouble PCy, QUICKDouble PCz,\n\
                QUICKDouble Zeta, QUICKDouble* YVerticalTemp); \n" % (m))          
            self.fhc.write("}; \n")

            # write function definitions
            self.fhd.write("__device__ __inline__ DDint_%d::DDint_%d(QUICKDouble PAx, QUICKDouble PAy, QUICKDouble PAz,\n\
                QUICKDouble PBx, QUICKDouble PBy, QUICKDouble PBz, QUICKDouble PCx, QUICKDouble PCy, QUICKDouble PCz,\n\
                QUICKDouble Zeta, QUICKDouble* YVerticalTemp){ \n\n" % (m, m))
            self.fhd.write("  PPint_%d pp_%d(PAx, PAy, PAz, PBx, PBy, PBz, PCx, PCy, PCz, Zeta, YVerticalTemp); // construct [p|p] for m=%d \n" % (m, m, m))
            self.fhd.write("  DSint_%d ds_%d(PAx, PAy, PAz, PCx, PCy, PCz, Zeta, YVerticalTemp); // construct [d|s] for m=%d \n" % (m, m, m))
            self.fhd.write("  DPint_%d dp_%d(PAx, PAy, PAz, PBx, PBy, PBz, PCx, PCy, PCz, Zeta, YVerticalTemp); // construct [d|p] for m=%d \n" % (m, m, m))            
            self.fhd.write("  PPint_%d pp_%d(PAx, PAy, PAz, PBx, PBy, PBz, PCx, PCy, PCz, Zeta, YVerticalTemp); // construct [p|p] for m=%d \n" % (m+1, m+1, m+1))
            self.fhd.write("  DSint_%d ds_%d(PAx, PAy, PAz, PCx, PCy, PCz, Zeta, YVerticalTemp); // construct [d|s] for m=%d \n" % (m+1, m+1, m+1))
            self.fhd.write("  DPint_%d dp_%d(PAx, PAy, PAz, PBx, PBy, PBz, PCx, PCy, PCz, Zeta, YVerticalTemp); // construct [d|p] for m=%d \n\n" % (m+1, m+1, m+1))             

            for i in range(0,6):
                for j in range(0,6):
                    tmp_mcal1=[params.Mcal[i+4][0], params.Mcal[i+4][1], params.Mcal[i+4][2]]
                    tmp_mcal2=[params.Mcal[j+4][0], params.Mcal[j+4][1], params.Mcal[j+4][2]]
                    for k in range(0,3):
                        if params.Mcal[j+4][k] != 0:
                            tmp_mcal2[k] -= 1
                            tmp_j=params.trans[tmp_mcal2[0]][tmp_mcal2[1]][tmp_mcal2[2]]

                            iclass_obj="dp"
                            self.fhd.write("  x_%d_%d = %s * %s_%d.x_%d_%d - %s * %s_%d.x_%d_%d; \n" % (j+4, i+4, self.PB[k], iclass_obj, m, i+4, tmp_j-1,\
                            self.PC[k], iclass_obj, m+1, i+4, tmp_j-1))

                            if params.Mcal[j+4][k] > 1:
                                iclass_obj="ds"
                                self.fhd.write("  x_%d_%d += 0.5/Zeta * %f * (%s_%d.x_%d_%d - %s_%d.x_%d_%d); \n" % (j+4, i+4, params.Mcal[j+4][k]-1, iclass_obj, m, i+4, 0, \
                                iclass_obj, m+1, i+4, 0))

                            if tmp_mcal1[k] > 0:
                                tmp_mcal1[k] -= 1
                                tmp_i=params.trans[tmp_mcal1[0]][tmp_mcal1[1]][tmp_mcal1[2]]
                                iclass_obj="pp"
                                self.fhd.write("  x_%d_%d += 0.5/Zeta * %f * (%s_%d.x_%d_%d - %s_%d.x_%d_%d); \n" % (j+4, i+4, params.Mcal[i+4][k], iclass_obj, m, tmp_i-1, tmp_j-1,\
                                iclass_obj, m+1, tmp_i-1, tmp_j-1))
                            break
            self.fhd.write("\n } \n")

    # generate code to save computed [d|d] integral
    def save_int(self):
        self.fha.write("\n  /* DD integral, m=%d */ \n" % (0))
        self.fha.write("  if(I == 2 && J == 2){ \n")
        self.fha.write("    DDint_0 dd(PAx, PAy, PAz, PBx, PBy, PBz, PCx, PCy, PCz, Zeta, YVerticalTemp); \n")
        for i in range(0,6):
            for j in range(0,6):
                self.fha.write("    LOC2(store, %d, %d, STOREDIM, STOREDIM) = dd.x_%d_%d;\n" % (j+4, i+4, j+4, i+4))

        # include print statements if debug option is on    
        if OEint.debug == 1:
            self.fha.write("\n#ifdef DEBUG_OEI \n")
            for i in range(0,6):
                for j in range(0,6):
                    self.fha.write("    printf(\"II %%d JJ %%d %s store[%d,%d] = %%f \\n\", II, JJ, LOC2(store, %d, %d, STOREDIM, STOREDIM)); \n" % ( "DD", j+4, i+4, j+4, i+4))
            self.fha.write("#endif \n\n")

        self.fha.write("  } \n")

def write_oei():

    # set files
    OEint.fhc = open('gpu_oei_classes.h','w')
    OEint.fhd = open('gpu_oei_definitions.h','w')
    OEint.fha= open('gpu_oei_assembler.h','w')

    # write license info
    file_handler.write_license(OEint.fhc)
    file_handler.write_license(OEint.fhd)
    file_handler.write_license(OEint.fha)

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
    OEint.fha.write("__device__ __inline__ void OEint_vertical(int I, int J, QUICKDouble PAx, QUICKDouble PAy, QUICKDouble PAz,\n\
        QUICKDouble PBx, QUICKDouble PBy, QUICKDouble PBz, QUICKDouble PCx, QUICKDouble PCy, QUICKDouble PCz, QUICKDouble Zeta,\n\
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
    OEint.fha.write("\n } \n")

    OEint.fhc.close()
    OEint.fhd.close()
    OEint.fha.close()

write_oei()
