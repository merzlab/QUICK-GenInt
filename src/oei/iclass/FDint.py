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

# [f|d] class, subclass of OEint
class FDint(OEint):
    def gen_int(self):
        # write code paths for integrals. Note that we use C++ classes here.
        for m in range(0,self.max_m+1):
            if m == 0:
                self.fhc.write("\n/* FD true integral, m=%d */ \n" % (m)) 
                self.fhd.write("\n/* FD true integral, m=%d */ \n" % (m)) 
            else:
                self.fhc.write("\n/* FD auxilary integral, m=%d */ \n" % (m))
                self.fhd.write("\n/* FD auxilary integral, m=%d */ \n" % (m))              

            self.fhc.write("class FDint_%d{ \n" % (m))
            self.fhc.write("public: \n")

            # write class variables; convention being used is s=0, p=1-3, d=4-9, f=10-19, g=20-34
            for i in range(0,10):
                for j in range(0,6):
                    self.fhc.write("  QUICKDouble x_%d_%d; // %s, %s \n" % (i+10, j+4, self.f_lbl[i], self.d_lbl[j]))

            # write class functions
            self.fhc.write("  __device__ __inline__ FDint_%d(QUICKDouble PAx, QUICKDouble PAy, QUICKDouble PAz,\n\
                QUICKDouble PBx, QUICKDouble PBy, QUICKDouble PBz, QUICKDouble PCx, QUICKDouble PCy, QUICKDouble PCz,\n\
                QUICKDouble Zeta, QUICKDouble* YVerticalTemp); \n" % (m))          
            self.fhc.write("}; \n")

            # write function definitions
            self.fhd.write("__device__ __inline__ FDint_%d::FDint_%d(QUICKDouble PAx, QUICKDouble PAy, QUICKDouble PAz,\n\
                QUICKDouble PBx, QUICKDouble PBy, QUICKDouble PBz, QUICKDouble PCx, QUICKDouble PCy, QUICKDouble PCz,\n\
                QUICKDouble Zeta, QUICKDouble* YVerticalTemp){ \n\n" % (m, m))
            self.fhd.write("  DPint_%d dp_%d(PAx, PAy, PAz, PBx, PBy, PBz, PCx, PCy, PCz, Zeta, YVerticalTemp); // construct [d|p] for m=%d \n" % (m, m, m))
            self.fhd.write("  FSint_%d fs_%d(PAx, PAy, PAz, PCx, PCy, PCz, Zeta, YVerticalTemp); // construct [f|s] for m=%d \n" % (m, m, m))
            self.fhd.write("  FPint_%d fp_%d(PAx, PAy, PAz, PBx, PBy, PBz, PCx, PCy, PCz, Zeta, YVerticalTemp); // construct [f|p] for m=%d \n" % (m, m, m))            
            self.fhd.write("  DPint_%d dp_%d(PAx, PAy, PAz, PBx, PBy, PBz, PCx, PCy, PCz, Zeta, YVerticalTemp); // construct [d|p] for m=%d \n" % (m+1, m+1, m+1))
            self.fhd.write("  FSint_%d fs_%d(PAx, PAy, PAz, PCx, PCy, PCz, Zeta, YVerticalTemp); // construct [f|s] for m=%d \n" % (m+1, m+1, m+1))
            self.fhd.write("  FPint_%d fp_%d(PAx, PAy, PAz, PBx, PBy, PBz, PCx, PCy, PCz, Zeta, YVerticalTemp); // construct [f|p] for m=%d \n\n" % (m+1, m+1, m+1))             

            for i in range(0,10):
                for j in range(0,6):
                    tmp_mcal1=[params.Mcal[i+10][0], params.Mcal[i+10][1], params.Mcal[i+10][2]]
                    tmp_mcal2=[params.Mcal[j+4][0], params.Mcal[j+4][1], params.Mcal[j+4][2]]

                    for k in range(0,3):
                        if params.Mcal[j+4][k] != 0:
                            tmp_mcal2[k] -= 1
                            tmp_j=params.trans[tmp_mcal2[0]][tmp_mcal2[1]][tmp_mcal2[2]]
                            iclass_obj="fp"
                            self.fhd.write("  x_%d_%d = %s * %s_%d.x_%d_%d - %s * %s_%d.x_%d_%d; \n" % (i+10, j+4, self.PB[k], iclass_obj, m, i+10, tmp_j-1,\
                            self.PC[k], iclass_obj, m+1, i+10, tmp_j-1))

                            if params.Mcal[j+4][k] > 1:
                                iclass_obj="fs"
                                self.fhd.write("  x_%d_%d += 0.5/Zeta * (%s_%d.x_%d_%d - %s_%d.x_%d_%d); \n" % (i+10, j+4, iclass_obj, m, i+10, 0, \
                                iclass_obj, m+1, i+10, 0))

                            if tmp_mcal1[k] > 0:
                                tmp_mcal1[k] -= 1
                                tmp_i=params.trans[tmp_mcal1[0]][tmp_mcal1[1]][tmp_mcal1[2]]
                                iclass_obj="dp"
                                self.fhd.write("  x_%d_%d += 0.5/Zeta * %f * (%s_%d.x_%d_%d - %s_%d.x_%d_%d); \n" % (i+10, j+4, params.Mcal[i+10][k], iclass_obj, m, tmp_i-1, tmp_j-1,\
                                iclass_obj, m+1, tmp_i-1, tmp_j-1))
                            break
            self.fhd.write("\n } \n")

    # generate code to save computed [f|d] integral
    def save_int(self):
        self.fha.write("\n  /* FD integral, m=%d */ \n" % (0))
        self.fha.write("  if(I == 3 && J == 2){ \n")
        self.fha.write("    FDint_0 fd(PAx, PAy, PAz, PBx, PBy, PBz, PCx, PCy, PCz, Zeta, YVerticalTemp); \n")
        for i in range(0,10):
            for j in range(0,6):
                self.fha.write("    LOC2(store, %d, %d, STOREDIM, STOREDIM) += fd.x_%d_%d;\n" % (i+10, j+4, i+10, j+4))

        # include print statements if debug option is on    
        if OEint.debug == 1:
            self.fha.write("\n#ifdef DEBUG_OEI \n")
            for i in range(0,10):
                for j in range(0,6):
                    self.fha.write("    printf(\"II %%d JJ %%d %s store[%d,%d] = %%f \\n\", II, JJ, LOC2(store, %d, %d, STOREDIM, STOREDIM)); \n" % ( "FD", i+10, j+4, i+10, j+4))
            self.fha.write("#endif \n\n")

        self.fha.write("  } \n")

