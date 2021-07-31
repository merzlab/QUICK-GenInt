#!---------------------------------------------------------------------!
#! Written by Madu Manathunga on 07/01/2021                            !
#!                                                                     !
#! Copyright (C) 2020-2021 Merz lab                                    !
#! Copyright (C) 2020-2021 Götz lab                                    !
#!                                                                     !
#! This source file is a part of QUICK-GenInt code generator and       !
#! is subjected to the terms of the Mozilla Public License, v. 2.0.    !
#! If a copy of the MPL was not distributed with this file, you can    !
#! obtain one at http://mozilla.org/MPL/2.0/.                          !
#!_____________________________________________________________________!

#!---------------------------------------------------------------------!
#! This source file contains classes necessary for generating one      !
#! electron integrals. Note that we use vertical recurrence relations  !
#! algorithm developed by Obara and Saika. See J. Chem. Phys. 1986, 84,!
#! 3963−3974 paper for theoretical details.                            !
#!                                                                     !
#!---------------------------------------------------------------------!

import src.common.params as params
import src.common.file_handler as file_handler
from src.oei.iclass.OEint import OEint

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
            self.fha.write("    LOC2(store, %d, %d, STOREDIM, STOREDIM) += sp.x_%d_%d;\n" % (0, i+1, 0, i+1))

        # include print statements if debug option is on    
        if OEint.debug == 1:
            self.fha.write("\n#ifdef DEBUG_OEI \n")
            for i in range(0,3):
                self.fha.write("    printf(\"II %%d JJ %%d %s store[%d,%d] = %%f \\n\", II, JJ, LOC2(store, %d, %d, STOREDIM, STOREDIM)); \n" % ( "SP", 0, i+1, 0, i+1))
            self.fha.write("#endif \n\n")

        self.fha.write("  } \n")

    # generate code to save [s|p] integral gradients
    def save_int_grad(self):
        self.fhga.write("\n  /* SP integral gradient, m=%d */ \n" % (0))
        self.fhga.write("  if(I == 0 && J == 1){ \n")
        self.fhga.write("    PPint_0 pp(PAx, PAy, PAz, PBx, PBy, PBz, PCx, PCy, PCz, Zeta, YVerticalTemp); \n")
        self.fhga.write("    SDint_0 sd(PBx, PBy, PBz, PCx, PCy, PCz, Zeta, YVerticalTemp); \n\n")

        self.fhga.write("    LOC2(store, 0, 0, STOREDIM, STOREDIM) += VY(0, 0, 0);\n")

        for i in range(0,3):
            for j in range(0,3):
                self.fhga.write("    LOC2(store, %d, %d, STOREDIM, STOREDIM) += pp.x_%d_%d;\n" % (i+1, j+1, i+1, j+1))  

        for i in range(0,6):
            self.fhga.write("    LOC2(store, %d, %d, STOREDIM, STOREDIM) += sd.x_%d_%d;\n" % (0, i+4, 0, i+4))      

        if OEint.debug == 1:
            self.fhga.write("\n#ifdef DEBUG_OEI \n")
            self.fhga.write("    printf(\"II %%d JJ %%d %s store[%d,%d] = %%f \\n\", II, JJ, LOC2(store, 0, 0, STOREDIM, STOREDIM)); \n" % ( "SS", 0, 0))
            for i in range(0,3):
                for j in range(0,3):
                    self.fhga.write("    printf(\"II %%d JJ %%d %s store[%d,%d] = %%f \\n\", II, JJ, LOC2(store, %d, %d, STOREDIM, STOREDIM)); \n" % ( "PP", i+1, j+1, i+1, j+1))

            for i in range(0,6):
                self.fhga.write("    printf(\"II %%d JJ %%d %s store[%d,%d] = %%f \\n\", II, JJ, LOC2(store, %d, %d, STOREDIM, STOREDIM)); \n" % ( "SD", 0, i+4, 0, i+4))
            self.fhga.write("#endif \n\n")

        self.fhga.write("  } \n")
