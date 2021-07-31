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
                self.fha.write("    LOC2(store, %d, %d, STOREDIM, STOREDIM) += pp.x_%d_%d;\n" % (i+1, j+1, i+1, j+1))

        # include print statements if debug option is on    
        if OEint.debug == 1:
            self.fha.write("\n#ifdef DEBUG_OEI \n")
            for i in range(0,3):
                for j in range(0,3):
                    self.fha.write("    printf(\"II %%d JJ %%d %s store[%d,%d] = %%f \\n\", II, JJ, LOC2(store, %d, %d, STOREDIM, STOREDIM)); \n" % ( "PP", i+1, j+1, i+1, j+1))
            self.fha.write("#endif \n\n")

        self.fha.write("  } \n")

    # generate code to save [p|p] integral gradients
    def save_int_grad(self):
        self.fhga.write("\n  /* PP integral gradient, m=%d */ \n" % (0))
        self.fhga.write("  if(I == 1 && J == 1){ \n")
        self.fhga.write("    SPint_0 sp(PBx, PBy, PBz, PCx, PCy, PCz, YVerticalTemp); \n")
        self.fhga.write("    PSint_0 ps(PAx, PAy, PAz, PCx, PCy, PCz, YVerticalTemp); \n")
        self.fhga.write("    DPint_0 dp(PAx, PAy, PAz, PBx, PBy, PBz, PCx, PCy, PCz, Zeta, YVerticalTemp); \n")
        self.fhga.write("    PDint_0 pd(PAx, PAy, PAz, PBx, PBy, PBz, PCx, PCy, PCz, Zeta, YVerticalTemp); \n\n")

        for i in range(0,3):                
            self.fhga.write("    LOC2(store, %d, %d, STOREDIM, STOREDIM) += sp.x_%d_%d;\n" % (0, i+1, 0, i+1))

        for i in range(0,3):                
            self.fhga.write("    LOC2(store, %d, %d, STOREDIM, STOREDIM) += ps.x_%d_%d;\n" % (i+1, 0, i+1, 0))

        for i in range(0,6):
            for j in range(0,3):
                self.fhga.write("    LOC2(store, %d, %d, STOREDIM, STOREDIM) += dp.x_%d_%d;\n" % (i+4, j+1, i+4, j+1)) 

        for i in range(0,6):
            for j in range(0,3):
                self.fhga.write("    LOC2(store, %d, %d, STOREDIM, STOREDIM) += pd.x_%d_%d;\n" % (j+1, i+4, j+1, i+4))
     
        if OEint.debug == 1:
            self.fhga.write("\n#ifdef DEBUG_OEI \n")

            for i in range(0,3):
                self.fhga.write("    printf(\"II %%d JJ %%d %s store[%d,%d] = %%f \\n\", II, JJ, LOC2(store, %d, %d, STOREDIM, STOREDIM)); \n" % ( "SP", 0, i+1, 0, i+1))

            for i in range(0,3):
                self.fhga.write("    printf(\"II %%d JJ %%d %s store[%d,%d] = %%f \\n\", II, JJ, LOC2(store, %d, %d, STOREDIM, STOREDIM)); \n" % ( "PS", i+1, 0, i+1, 0))

            for i in range(0,6):
                for j in range(0,3):
                    self.fhga.write("    printf(\"II %%d JJ %%d %s store[%d,%d] = %%f \\n\", II, JJ, LOC2(store, %d, %d, STOREDIM, STOREDIM)); \n" % ( "DP", i+4, j+1, i+4, j+1))

            for i in range(0,6):
                for j in range(0,3):
                    self.fhga.write("    printf(\"II %%d JJ %%d %s store[%d,%d] = %%f \\n\", II, JJ, LOC2(store, %d, %d, STOREDIM, STOREDIM)); \n" % ( "PD", j+1, i+4, j+1, i+4))
                    
            self.fhga.write("#endif \n\n")

        self.fhga.write("  } \n")

