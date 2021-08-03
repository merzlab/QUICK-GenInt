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
                            self.fhd.write("  x_%d_%d = %s * %s_%d.x_%d_%d - %s * %s_%d.x_%d_%d; \n" % (i+4, j+4, self.PB[k], iclass_obj, m, i+4, tmp_j-1,\
                            self.PC[k], iclass_obj, m+1, i+4, tmp_j-1))

                            if params.Mcal[j+4][k] > 1:
                                iclass_obj="ds"
                                self.fhd.write("  x_%d_%d += 0.5/Zeta * %f * (%s_%d.x_%d_%d - %s_%d.x_%d_%d); \n" % (i+4, j+4, params.Mcal[j+4][k]-1, iclass_obj, m, i+4, 0, \
                                iclass_obj, m+1, i+4, 0))

                            if tmp_mcal1[k] > 0:
                                tmp_mcal1[k] -= 1
                                tmp_i=params.trans[tmp_mcal1[0]][tmp_mcal1[1]][tmp_mcal1[2]]
                                iclass_obj="pp"
                                self.fhd.write("  x_%d_%d += 0.5/Zeta * %f * (%s_%d.x_%d_%d - %s_%d.x_%d_%d); \n" % (i+4, j+4, params.Mcal[i+4][k], iclass_obj, m, tmp_i-1, tmp_j-1,\
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
                self.fha.write("    LOC2(store, %d, %d, STOREDIM, STOREDIM) += dd.x_%d_%d;\n" % (i+4, j+4, i+4, j+4))

        # include print statements if debug option is on    
        if OEint.debug == 1:
            self.fha.write("\n#ifdef DEBUG_OEI \n")
            for i in range(0,6):
                for j in range(0,6):
                    self.fha.write("    printf(\"II %%d JJ %%d %s store[%d,%d] = %%f \\n\", II, JJ, LOC2(store, %d, %d, STOREDIM, STOREDIM)); \n" % ( "DD", i+4, j+4, i+4, j+4))
            self.fha.write("#endif \n\n")

        self.fha.write("  } \n")

    # generate code to save [d|d] integral gradients
    def save_int_grad(self):
        self.fhga.write("\n  /* DD integral gradient, m=%d */ \n" % (0))
        self.fhga.write("  if(I == 2 && J == 2){ \n")
        self.fhga.write("    PDint_0 pd(PAx, PAy, PAz, PBx, PBy, PBz, PCx, PCy, PCz, Zeta, YVerticalTemp); \n")
        self.fhga.write("    DPint_0 dp(PAx, PAy, PAz, PBx, PBy, PBz, PCx, PCy, PCz, Zeta, YVerticalTemp); \n")
        self.fhga.write("    FDint_0 fd(PAx, PAy, PAz, PBx, PBy, PBz, PCx, PCy, PCz, Zeta, YVerticalTemp); \n")
        self.fhga.write("    DFint_0 df(PAx, PAy, PAz, PBx, PBy, PBz, PCx, PCy, PCz, Zeta, YVerticalTemp); \n\n")

        for i in range(0,6):
            for j in range(0,3):
                self.fhga.write("    LOC2(store2, %d, %d, STOREDIM, STOREDIM) = pd.x_%d_%d;\n" % (j+1, i+4, j+1, i+4))

        for i in range(0,6):
            for j in range(0,3):
                self.fhga.write("    LOC2(store2, %d, %d, STOREDIM, STOREDIM) = dp.x_%d_%d;\n" % (i+4, j+1, i+4, j+1))

        for i in range(0,10):
            for j in range(0,6):
                self.fhga.write("    LOC2(store2, %d, %d, STOREDIM, STOREDIM) = fd.x_%d_%d;\n" % (i+10, j+4, i+10, j+4))

        for i in range(0,10):
            for j in range(0,6):
                self.fhga.write("    LOC2(store2, %d, %d, STOREDIM, STOREDIM) = df.x_%d_%d;\n" % (j+4, i+10, j+4, i+10))

        if OEint.debug == 1:
            self.fhga.write("\n#ifdef DEBUG_OEI \n")

            for i in range(0,6):
                for j in range(0,3):
                    self.fhga.write("    printf(\"II %%d JJ %%d %s store2[%d,%d] = %%f \\n\", II, JJ, LOC2(store2, %d, %d, STOREDIM, STOREDIM)); \n" % ( "PD", j+1, i+4, j+1, i+4))

            for i in range(0,6):
                for j in range(0,3):
                    self.fhga.write("    printf(\"II %%d JJ %%d %s store2[%d,%d] = %%f \\n\", II, JJ, LOC2(store2, %d, %d, STOREDIM, STOREDIM)); \n" % ( "DP", i+4, j+1, i+4, j+1))

            for i in range(0,10):
                for j in range(0,6):
                    self.fhga.write("    printf(\"II %%d JJ %%d %s store2[%d,%d] = %%f \\n\", II, JJ, LOC2(store2, %d, %d, STOREDIM, STOREDIM)); \n" % ( "FD", i+10, j+4, i+10, j+4))

            for i in range(0,10):
                for j in range(0,6):
                    self.fhga.write("    printf(\"II %%d JJ %%d %s store[2%d,%d] = %%f \\n\", II, JJ, LOC2(store2, %d, %d, STOREDIM, STOREDIM)); \n" % ( "DF", j+4, i+10, j+4, i+10))

            self.fhga.write("#endif \n\n")

        self.fhga.write("  } \n")

