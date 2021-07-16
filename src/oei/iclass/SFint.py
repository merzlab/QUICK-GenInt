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

# [s|f] class, subclass of OEint
class SFint(OEint):
    def gen_int(self):
        # write code paths for integrals. Note that we use C++ classes here.
        for m in range(0,self.max_m+1):
            if m == 0:
                self.fhc.write("\n/* SF true integral, m=%d */ \n" % (m)) 
                self.fhd.write("\n/* SF true integral, m=%d */ \n" % (m)) 
            else:
                self.fhc.write("\n/* SF auxilary integral, m=%d */ \n" % (m))
                self.fhd.write("\n/* SF auxilary integral, m=%d */ \n" % (m))              

            self.fhc.write("class SFint_%d{ \n" % (m))
            self.fhc.write("public: \n")

            # write class variables; convention being used is s=0, p=1-3, d=4-9, f=10-19, g=20-34
            for i in range(0,10):
                self.fhc.write("  QUICKDouble x_%d_%d; // %s, %s \n" % (0, i+10, "S", self.f_lbl[i]))

            # write class functions
            self.fhc.write("  __device__ __inline__ SFint_%d(QUICKDouble PBx, QUICKDouble PBy, QUICKDouble PBz,\n\
                QUICKDouble PCx, QUICKDouble PCy, QUICKDouble PCz, QUICKDouble Zeta, QUICKDouble* YVerticalTemp); \n" % (m))          
            self.fhc.write("}; \n")

            # write function definitions
            self.fhd.write("__device__ __inline__ SFint_%d::SFint_%d(QUICKDouble PBx, QUICKDouble PBy, QUICKDouble PBz,\n\
                QUICKDouble PCx, QUICKDouble PCy, QUICKDouble PCz, QUICKDouble Zeta, QUICKDouble* YVerticalTemp){ \n\n" % (m, m))
            self.fhd.write("  SPint_%d sp_%d(PBx, PBy, PBz, PCx, PCy, PCz, YVerticalTemp); // construct [s|p] for m=%d \n" % (m, m, m))
            self.fhd.write("  SPint_%d sp_%d(PBx, PBy, PBz, PCx, PCy, PCz, YVerticalTemp); // construct [s|p] for m=%d \n" % (m+1, m+1, m+1))
            self.fhd.write("  SDint_%d sd_%d(PBx, PBy, PBz, PCx, PCy, PCz, YVerticalTemp); // construct [s|d] for m=%d \n" % (m, m, m))
            self.fhd.write("  SDint_%d sd_%d(PBx, PBy, PBz, PCx, PCy, PCz, YVerticalTemp); // construct [s|d] for m=%d \n\n" % (m+1, m+1, m+1))

            for i in range(0,10):
                tmp_mcal=[params.Mcal[i+10][0], params.Mcal[i+10][1], params.Mcal[i+10][2]]
                for j in range(0,3):
                    if params.Mcal[i+10][j] != 0:
                        tmp_mcal[j] -= 1
                        tmp_i=params.trans[tmp_mcal[0]][tmp_mcal[1]][tmp_mcal[2]]
                        self.fhd.write("  x_%d_%d = %s * sd_%d.x_%d_%d - %s * sd_%d.x_%d_%d; \n" % (0, i+10, self.PB[j], m, 0, tmp_i-1,\
                        self.PC[j], m+1, 0, tmp_i-1))

                        if params.Mcal[i+10][j] > 1:
                            tmp_mcal[j] = params.Mcal[i+10][j] - 2
                            tmp_i=params.trans[tmp_mcal[0]][tmp_mcal[1]][tmp_mcal[2]]
                            self.fhd.write("  x_%d_%d += 0.5/Zeta * %f * (sp_%d.x_%d_%d - sp_%d.x_%d_%d); \n" % (0, i+10, params.Mcal[i+10][j] - 1,\
                                m, 0, tmp_i-1, m+1, 0, tmp_i-1))

                        break
            self.fhd.write("\n } \n")

    # generate code to save computed [s|f] integral
    def save_int(self):
        self.fha.write("\n  /* SF integral, m=%d */ \n" % (0))
        self.fha.write("  if(I == 0 && J == 3){ \n")
        self.fha.write("    SFint_0 sf(PBx, PBy, PBz, PCx, PCy, PCz, Zeta, YVerticalTemp); \n")
        for i in range(0,10):
            self.fha.write("    LOC2(store, %d, %d, STOREDIM, STOREDIM) += sf.x_%d_%d;\n" % (0, i+10, 0, i+10))

        # include print statements if debug option is on    
        if OEint.debug == 1:
            self.fha.write("\n#ifdef DEBUG_OEI \n")
            for i in range(0,10):
                self.fha.write("    printf(\"II %%d JJ %%d %s store[%d,%d] = %%f \\n\", II, JJ, LOC2(store, %d, %d, STOREDIM, STOREDIM)); \n" % ( "SF", 0, i+10, 0, i+10))
            self.fha.write("#endif \n\n")

        self.fha.write("  } \n")

