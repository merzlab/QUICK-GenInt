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

# parent class for all one electron integrals
class OEint:
    # set file handlers
    fhc = 0 # file handler for class declarations
    fhd = 0 # file handler for function implementations
    fha= 0  # file handler for integral assembler
    fhga= 0 # file handler for integral gradient assembler
    debug=1 # include debug info in generated code, 0=no, 1=yes 

    # set function qualifiers
    func_qualifier='__device__ __inline__'

    # max_m ranges from 0 to a+b; where a and b are the angular momentum of i and j
    # of the integral being considered [i|j]. If max_m=0, the integral is a true integral
    # and m>0 results in auxliary integrals. 
    def __init__(self,max_m=0):
        self.max_m=max_m

        if self.fhc == 0 or self.fhd == 0:
            print("Warning: file handlers in OEint class are not set. \n")
        
        # set orbital labels for improving readability 
        self.s_lbl=("S")
        self.p_lbl=("Px", "Py", "Pz")
        self.d_lbl=("Dxy", "Dyz", "Dxz", "Dxx", "Dyy", "Dzz")
        self.f_lbl=("Fxyz", "Fxxy", "Fxyy", "Fxxz", "Fxzz", "Fyyz", "Fyzz", "Fxxx", "Fyyy", "Fzzz")
        self.g_lbl=("Gxxyy", "Gxxzz", "Gyyzz", "Gxxyz", "Gxyyz", "Gxyzz", "Gxxxz", "Gxzzz", "Gxxxy", "Gxyyy", "Gyyyz", "Gyzzz", "Gxxxx", "Gyyyy", "Gzzzz")

        # set labels for useful variables
        self.AA=("Ax", "Ay", "Az") # position of center A
        self.BB=("Bx", "By", "Bz") # position of center B
        self.CC=("Cx", "Cy", "Cz") # position of center C
        self.PP=("Px", "Py", "Pz") # quantity P in OS eqn 15
        self.PA=("PAx", "PAy", "PAz") # Px-Ax, Py-Ay, Pz-Az
        self.PB=("PBx", "PBy", "PBz") # Px-Bx, Py-By, Pz-Bz
        self.PC=("PCx", "PCy", "PCz") # Px-Cx, Py-Cy, Pz-Cz

    # generate source code for integrals
    def gen_int(self):
        pass 

    # save true integrals into store array
    def save_int(self):
        pass
