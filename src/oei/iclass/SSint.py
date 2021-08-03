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

# [s|s] class, a subclass of OEint
class SSint(OEint):

    # generate code to save computed [s|s] integral
    def save_int(self):
        self.fha.write("\n  /* SS integral, m=%d */ \n" % (0))
        self.fha.write("  if(I == 0 && J == 0){ \n")
        self.fha.write("    LOC2(store, 0, 0, STOREDIM, STOREDIM) += VY(0, 0, 0);\n")

        # include print statements if debug option is on 
        if OEint.debug == 1:
            self.fha.write("\n#ifdef DEBUG_OEI \n")
            self.fha.write("    printf(\"II %%d JJ %%d %s store[%d,%d] = %%f \\n\", II, JJ, LOC2(store, 0, 0, STOREDIM, STOREDIM)); \n" % ( "SS", 0, 0))
            self.fha.write("#endif \n\n")

        self.fha.write("  } \n")

    # generate code to save [s|s] integral gradients
    def save_int_grad(self):
        self.fhga.write("\n  /* SS integral gradient, m=%d */ \n" % (0))
        self.fhga.write("  if(I == 0 && J == 0){ \n")
        self.fhga.write("    SPint_0 sp(PBx, PBy, PBz, PCx, PCy, PCz, YVerticalTemp); \n")
        self.fhga.write("    PSint_0 ps(PAx, PAy, PAz, PCx, PCy, PCz, YVerticalTemp); \n\n")

        for i in range(0,3):                
            self.fhga.write("    LOC2(store2, %d, %d, STOREDIM, STOREDIM) = sp.x_%d_%d;\n" % (0, i+1, 0, i+1))

        for i in range(0,3):                
            self.fhga.write("    LOC2(store2, %d, %d, STOREDIM, STOREDIM) = ps.x_%d_%d;\n" % (i+1, 0, i+1, 0))

        if OEint.debug == 1:
            self.fhga.write("\n#ifdef DEBUG_OEI \n")
            for i in range(0,3):
                self.fhga.write("    printf(\"II %%d JJ %%d %s store2[%d,%d] = %%f \\n\", II, JJ, LOC2(store2, %d, %d, STOREDIM, STOREDIM)); \n" % ( "SP", 0, i+1, 0, i+1))
            for i in range(0,3):
                self.fhga.write("    printf(\"II %%d JJ %%d %s store2[%d,%d] = %%f \\n\", II, JJ, LOC2(store2, %d, %d, STOREDIM, STOREDIM)); \n" % ( "PS", i+1, 0, i+1, 0))            
            self.fhga.write("#endif \n\n")

        self.fhga.write("  } \n")


