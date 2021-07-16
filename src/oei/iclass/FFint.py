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

# [f|f] class, subclass of OEint
class FFint(OEint):
    def gen_int(self):
        # write code paths for integrals. Note that we use C++ classes here.
        for m in range(0,self.max_m+1):
            if m == 0:
                self.fhc.write("\n/* FF true integral, m=%d */ \n" % (m)) 
                self.fhd.write("\n/* FF true integral, m=%d */ \n" % (m)) 
            else:
                self.fhc.write("\n/* FF auxilary integral, m=%d */ \n" % (m))
                self.fhd.write("\n/* FF auxilary integral, m=%d */ \n" % (m))              

            self.fhc.write("class FFint_%d{ \n" % (m))
            self.fhc.write("public: \n")

            # write class variables; convention being used is s=0, p=1-3, d=4-9, f=10-19, g=20-34
            for i in range(0,10):
                for j in range(0,10):
                    self.fhc.write("  QUICKDouble x_%d_%d; // %s, %s \n" % (i+10, j+10, self.f_lbl[i], self.f_lbl[j]))

            # write class functions
            self.fhc.write("  __device__ __inline__ FFint_%d(QUICKDouble PAx, QUICKDouble PAy, QUICKDouble PAz,\n\
                QUICKDouble PBx, QUICKDouble PBy, QUICKDouble PBz, QUICKDouble PCx, QUICKDouble PCy, QUICKDouble PCz,\n\
                QUICKDouble Zeta, QUICKDouble* YVerticalTemp); \n" % (m))          
            self.fhc.write("}; \n")

            # write function definitions
            self.fhd.write("__device__ __inline__ FFint_%d::FFint_%d(QUICKDouble PAx, QUICKDouble PAy, QUICKDouble PAz,\n\
                QUICKDouble PBx, QUICKDouble PBy, QUICKDouble PBz, QUICKDouble PCx, QUICKDouble PCy, QUICKDouble PCz,\n\
                QUICKDouble Zeta, QUICKDouble* YVerticalTemp){ \n\n" % (m, m))
            self.fhd.write("  DDint_%d dd_%d(PAx, PAy, PAz, PBx, PBy, PBz, PCx, PCy, PCz, Zeta, YVerticalTemp); // construct [d|d] for m=%d \n" % (m, m, m))
            self.fhd.write("  FPint_%d fp_%d(PAx, PAy, PAz, PBx, PBy, PBz, PCx, PCy, PCz, Zeta, YVerticalTemp); // construct [f|p] for m=%d \n" % (m, m, m))            
            self.fhd.write("  FDint_%d fd_%d(PAx, PAy, PAz, PBx, PBy, PBz, PCx, PCy, PCz, Zeta, YVerticalTemp); // construct [f|d] for m=%d \n" % (m, m, m))
            self.fhd.write("  DDint_%d dd_%d(PAx, PAy, PAz, PBx, PBy, PBz, PCx, PCy, PCz, Zeta, YVerticalTemp); // construct [d|d] for m=%d \n" % (m+1, m+1, m+1))
            self.fhd.write("  FPint_%d fp_%d(PAx, PAy, PAz, PBx, PBy, PBz, PCx, PCy, PCz, Zeta, YVerticalTemp); // construct [f|p] for m=%d \n" % (m+1, m+1, m+1))             
            self.fhd.write("  FDint_%d fd_%d(PAx, PAy, PAz, PBx, PBy, PBz, PCx, PCy, PCz, Zeta, YVerticalTemp); // construct [f|d] for m=%d \n\n" % (m+1, m+1, m+1))

            for i in range(0,10):
                for j in range(0,10):
                    tmp_mcal1=[params.Mcal[i+10][0], params.Mcal[i+10][1], params.Mcal[i+10][2]]
                    tmp_mcal2=[params.Mcal[j+10][0], params.Mcal[j+10][1], params.Mcal[j+10][2]]

                    for k in range(0,3):
                        if params.Mcal[j+10][k] != 0:
                            tmp_mcal2[k] -= 1
                            tmp_j=params.trans[tmp_mcal2[0]][tmp_mcal2[1]][tmp_mcal2[2]]
                            iclass_obj="fd"
                            self.fhd.write("  x_%d_%d = %s * %s_%d.x_%d_%d - %s * %s_%d.x_%d_%d; \n" % (j+10, i+10, self.PB[k], iclass_obj, m, i+10, tmp_j-1,\
                            self.PC[k], iclass_obj, m+1, i+10, tmp_j-1))

                            if params.Mcal[j+10][k] > 1:
                                tmp_mcal2[k] = params.Mcal[j+10][k] - 2
                                tmp_j2=params.trans[tmp_mcal2[0]][tmp_mcal2[1]][tmp_mcal2[2]]
                                iclass_obj="fp"
                                self.fhd.write("  x_%d_%d += 0.5/Zeta * %f * (%s_%d.x_%d_%d - %s_%d.x_%d_%d); \n" % (j+10, i+10, params.Mcal[j+10][k]-1, iclass_obj, m, i+10, tmp_j2-1, \
                                iclass_obj, m+1, i+10, tmp_j2-1))

                            if tmp_mcal1[k] > 0:
                                tmp_mcal1[k] -= 1
                                tmp_i=params.trans[tmp_mcal1[0]][tmp_mcal1[1]][tmp_mcal1[2]]
                                iclass_obj="dd"
                                self.fhd.write("  x_%d_%d += 0.5/Zeta * %f * (%s_%d.x_%d_%d - %s_%d.x_%d_%d); \n" % (j+10, i+10, params.Mcal[i+10][k], iclass_obj, m, tmp_i-1, tmp_j-1,\
                                iclass_obj, m+1, tmp_i-1, tmp_j-1))
                            break
            self.fhd.write("\n } \n")

    # generate code to save computed [f|f] integral
    def save_int(self):
        self.fha.write("\n  /* FF integral, m=%d */ \n" % (0))
        self.fha.write("  if(I == 3 && J == 3){ \n")
        self.fha.write("    FFint_0 ff(PAx, PAy, PAz, PBx, PBy, PBz, PCx, PCy, PCz, Zeta, YVerticalTemp); \n")
        for i in range(0,10):
            for j in range(0,10):
                self.fha.write("    LOC2(store, %d, %d, STOREDIM, STOREDIM) += ff.x_%d_%d;\n" % (i+10, j+10, i+10, j+10))

        # include print statements if debug option is on    
        if OEint.debug == 1:
            self.fha.write("\n#ifdef DEBUG_OEI \n")
            for i in range(0,10):
                for j in range(0,10):
                    self.fha.write("    printf(\"II %%d JJ %%d %s store[%d,%d] = %%f \\n\", II, JJ, LOC2(store, %d, %d, STOREDIM, STOREDIM)); \n" % ( "FF", i+10, j+10, i+10, j+10))
            self.fha.write("#endif \n\n")

        self.fha.write("  } \n")

