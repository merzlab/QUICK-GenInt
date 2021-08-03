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
            self.fha.write("    LOC2(store, %d, %d, STOREDIM, STOREDIM) += ds.x_%d_%d;\n" % (i+4, 0, i+4, 0))

        # include print statements if debug option is on    
        if OEint.debug == 1:
            self.fha.write("\n#ifdef DEBUG_OEI \n")
            for i in range(0,6):
                self.fha.write("    printf(\"II %%d JJ %%d %s store[%d,%d] = %%f \\n\", II, JJ, LOC2(store, %d, %d, STOREDIM, STOREDIM)); \n" % ( "DS", i+4, 0, i+4, 0))
            self.fha.write("#endif \n\n")

        self.fha.write("  } \n")

    # generate code to save [d|s] integral gradients
    def save_int_grad(self):
        self.fhga.write("\n  /* DS integral gradient, m=%d */ \n" % (0))
        self.fhga.write("  if(I == 2 && J == 0){ \n")
        self.fhga.write("    PSint_0 ps(PAx, PAy, PAz, PCx, PCy, PCz, YVerticalTemp); \n")
        self.fhga.write("    DPint_0 dp(PAx, PAy, PAz, PBx, PBy, PBz, PCx, PCy, PCz, Zeta, YVerticalTemp); \n")
        self.fhga.write("    FSint_0 fs(PAx, PAy, PAz, PCx, PCy, PCz, Zeta, YVerticalTemp); \n\n")

        for i in range(0,3):                
            self.fhga.write("    LOC2(store2, %d, %d, STOREDIM, STOREDIM) = ps.x_%d_%d;\n" % (i+1, 0, i+1, 0))

        for i in range(0,6):
            for j in range(0,3):
                self.fhga.write("    LOC2(store2, %d, %d, STOREDIM, STOREDIM) = dp.x_%d_%d;\n" % (i+4, j+1, i+4, j+1))

        for i in range(0,10):
            self.fhga.write("    LOC2(store2, %d, %d, STOREDIM, STOREDIM) = fs.x_%d_%d;\n" % (i+10, 0, i+10, 0))

        if OEint.debug == 1:
            self.fhga.write("\n#ifdef DEBUG_OEI \n")

            for i in range(0,3):
                self.fhga.write("    printf(\"II %%d JJ %%d %s store2[%d,%d] = %%f \\n\", II, JJ, LOC2(store2, %d, %d, STOREDIM, STOREDIM)); \n" % ( "PS", i+1, 0, i+1, 0))

            for i in range(0,6):
                for j in range(0,3):
                    self.fhga.write("    printf(\"II %%d JJ %%d %s store2[%d,%d] = %%f \\n\", II, JJ, LOC2(store2, %d, %d, STOREDIM, STOREDIM)); \n" % ( "DP", i+4, j+1, i+4, j+1))

            for i in range(0,10):
                self.fhga.write("    printf(\"II %%d JJ %%d %s store2[%d,%d] = %%f \\n\", II, JJ, LOC2(store2, %d, %d, STOREDIM, STOREDIM)); \n" % ( "FS", i+10, 0, i+10, 0))
            self.fhga.write("#endif \n\n")

        self.fhga.write("  } \n")
