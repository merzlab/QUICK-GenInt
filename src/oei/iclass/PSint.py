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
            self.fha.write("    LOC2(store, %d, %d, STOREDIM, STOREDIM) += ps.x_%d_%d;\n" % (i+1, 0, i+1, 0))

        # include print statements if debug option is on 
        if OEint.debug == 1:
            self.fha.write("\n#ifdef DEBUG_OEI \n")
            for i in range(0,3):
                self.fha.write("    printf(\"II %%d JJ %%d %s store[%d,%d] = %%f \\n\", II, JJ, LOC2(store, %d, %d, STOREDIM, STOREDIM)); \n" % ( "PS", i+1, 0, i+1, 0))
            self.fha.write("#endif \n\n")

        self.fha.write("  } \n")

