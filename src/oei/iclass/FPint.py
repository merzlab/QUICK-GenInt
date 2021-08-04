#!---------------------------------------------------------------------!
#! Written by Madu Manathunga on 07/16/2021                            !
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

# [f|p] class, subclass of OEint
class FPint(OEint):
    def gen_int(self):
        # write code paths for integrals. Note that we use C++ classes here.
        for m in range(0,self.max_m+1):
            if m == 0:
                self.fhc.write("\n/* FP true integral, m=%d */ \n" % (m)) 
                self.fhd.write("\n/* FP true integral, m=%d */ \n" % (m)) 
            else:
                self.fhc.write("\n/* FP auxilary integral, m=%d */ \n" % (m))
                self.fhd.write("\n/* FP auxilary integral, m=%d */ \n" % (m))              

            self.fhc.write("class FPint_%d{ \n" % (m))
            self.fhc.write("public: \n")

            # write class variables; convention being used is s=0, p=1-3, d=4-9, f=10-19, g=20-34
            for i in range(0,10):
                for j in range(0,3):
                    self.fhc.write("  QUICKDouble x_%d_%d; // %s, %s \n" % (i+10, j+1, self.f_lbl[i], self.p_lbl[j]))

            # write class functions
            self.fhc.write("  %s FPint_%d(QUICKDouble PAx, QUICKDouble PAy, QUICKDouble PAz,\n\
                QUICKDouble PBx, QUICKDouble PBy, QUICKDouble PBz, QUICKDouble PCx, QUICKDouble PCy, QUICKDouble PCz,\n\
                QUICKDouble Zeta, QUICKDouble* YVerticalTemp); \n" % (self.func_qualifier, m))          
            self.fhc.write("}; \n")

            # write function definitions
            self.fhd.write("%s FPint_%d::FPint_%d(QUICKDouble PAx, QUICKDouble PAy, QUICKDouble PAz,\n\
                QUICKDouble PBx, QUICKDouble PBy, QUICKDouble PBz, QUICKDouble PCx, QUICKDouble PCy, QUICKDouble PCz,\n\
                QUICKDouble Zeta, QUICKDouble* YVerticalTemp){ \n\n" % (self.func_qualifier, m, m))
            self.fhd.write("  DSint_%d ds_%d(PAx, PAy, PAz, PCx, PCy, PCz, Zeta, YVerticalTemp); // construct [d|s] for m=%d \n" % (m, m, m))            
            self.fhd.write("  DSint_%d ds_%d(PAx, PAy, PAz, PCx, PCy, PCz, Zeta, YVerticalTemp); // construct [d|s] for m=%d \n" % (m+1, m+1, m+1))
            self.fhd.write("  FSint_%d fs_%d(PAx, PAy, PAz, PCx, PCy, PCz, Zeta, YVerticalTemp); // construct [f|s] for m=%d \n" % (m, m, m))
            self.fhd.write("  FSint_%d fs_%d(PAx, PAy, PAz, PCx, PCy, PCz, Zeta, YVerticalTemp); // construct [f|s] for m=%d \n\n" % (m+1, m+1, m+1))

            for i in range(0,10):
                for j in range(0,3):
                    tmp_mcal=[params.Mcal[i+10][0], params.Mcal[i+10][1], params.Mcal[i+10][2]]
                    for k in range(0,3):
                        if params.Mcal[j+1][k] != 0:
                            self.fhd.write("  x_%d_%d = %s * fs_%d.x_%d_%d - %s * fs_%d.x_%d_%d; \n" % (i+10, j+1, self.PB[k], m, i+10, 0,\
                            self.PC[k], m+1, i+10, 0))


                            if tmp_mcal[k] != 0:
                                tmp_mcal[k] -= 1
                                tmp_i=params.trans[tmp_mcal[0]][tmp_mcal[1]][tmp_mcal[2]]
                                self.fhd.write("  x_%d_%d += 0.5/Zeta * %f * (ds_%d.x_%d_%d - ds_%d.x_%d_%d); \n" % (i+10, j+1, params.Mcal[i+10][k], m, tmp_i-1, 0, m+1, tmp_i-1, 0))

                            break
            self.fhd.write("\n } \n")


    # generate code to save computed [f|p] integral
    def save_int(self):
        self.fha.write("\n  /* FP integral, m=%d */ \n" % (0))
        self.fha.write("  if(I == 3 && J == 1){ \n")
        self.fha.write("    FPint_0 fp(PAx, PAy, PAz, PBx, PBy, PBz, PCx, PCy, PCz, Zeta, YVerticalTemp); \n")
        for i in range(0,10):
            for j in range(0,3):
                self.fha.write("    LOC2(store2, %d, %d, STOREDIM, STOREDIM) = fp.x_%d_%d;\n" % (i+10, j+1, i+10, j+1))

        # include print statements if debug option is on    
        if OEint.debug == 1:
            self.fha.write("\n#ifdef DEBUG_OEI \n")
            for i in range(0,10):
                for j in range(0,3):
                    self.fha.write("    printf(\"II %%d JJ %%d %s store2[%d,%d] = %%f \\n\", II, JJ, LOC2(store2, %d, %d, STOREDIM, STOREDIM)); \n" % ( "FP", i+10, j+1, i+10, j+1))
            self.fha.write("#endif \n\n")

        self.fha.write("  } \n")

