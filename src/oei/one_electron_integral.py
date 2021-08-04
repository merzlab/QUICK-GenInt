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
from src.oei.iclass.SSint import SSint
from src.oei.iclass.SPint import SPint
from src.oei.iclass.PSint import PSint
from src.oei.iclass.PPint import PPint
from src.oei.iclass.DSint import DSint
from src.oei.iclass.SDint import SDint
from src.oei.iclass.DPint import DPint
from src.oei.iclass.PDint import PDint
from src.oei.iclass.DDint import DDint
from src.oei.iclass.FSint import FSint
from src.oei.iclass.SFint import SFint
from src.oei.iclass.FPint import FPint
from src.oei.iclass.PFint import PFint
from src.oei.iclass.FDint import FDint
from src.oei.iclass.DFint import DFint
from src.oei.iclass.FFint import FFint

def write_oei(outdir, func_qualifier='__device__ __inline__'):

    # set files
    OEint.fhc = open(outdir+"/gpu_oei_classes.h",'w')
    OEint.fhd = open(outdir+"/gpu_oei_definitions.h",'w')
    OEint.fha= open(outdir+"/gpu_oei_assembler.h",'w')
    OEint.fhga= open(outdir+"/gpu_oei_grad_assembler.h",'w')

    # set function qualifiers
    OEint.func_qualifier=func_qualifier

    # write license info
    file_handler.write_license(OEint.fhc)
    file_handler.write_license(OEint.fhd)
    file_handler.write_license(OEint.fha)
    file_handler.write_license(OEint.fhga)

    # generate integral classes for systems containing only s, p, d and f functions.
    # generate [s|s], this is trivial and we will directly save the integral value from the driver. 
    ss=SSint()

    # generate [p|s]
    ps=PSint(5)
    ps.gen_int() 

    # generate [s|p]
    sp=SPint(5)
    sp.gen_int()

    # generate [p|p]
    pp=PPint(2)
    pp.gen_int()

    # generate [d|s]
    ds=DSint(4)
    ds.gen_int()

    # generate [s|d]
    sd=SDint(4)
    sd.gen_int() 

    # generate [d|p]
    dp=DPint(2)
    dp.gen_int()

    # generate [p|d]
    pd=PDint(2)
    pd.gen_int()

    # generate [d|d]
    dd=DDint(1)
    dd.gen_int()

    # generate [f|s]
    fs=FSint(3)
    fs.gen_int()

    # generate [s|f]
    sf=SFint(3)
    sf.gen_int()

    # generate [f|p]
    fp=FPint(2)
    fp.gen_int()

    # generate [p|f]
    pf=PFint(2)
    pf.gen_int()

    # generate [f|p]
    fd=FDint(1)
    fd.gen_int()    

    # generate [d|f]
    df=DFint(1)
    df.gen_int() 

    # generate [f|f]
    ff=FFint(0)
    ff.gen_int() 

    # Now we write the driver to save computed primitive integrals. The function parameters are as follows.
    # I, J - angular momentum of first and second shells (0, 1, 2, 3, and 4 for s, p, d, f and g respectively)
    # II, JJ - indices of shells
    # PAx, PAy, PAz - Difference between the first center (i.e. A) and weighted center (P) along x, y and z coordinates (i.e. Px-Ax, Py-Ay and Pz-Az)
    # PBx, PBy, PBz - Difference between the second center (i.e. B) and weighted center (P) along x, y and z coordinates (i.e. Px-Bx, Py-By and Pz-Bz)
    # PCx, PCy, PCz - Difference between the nuclei/point charge (i.e. C) and weighted center (P) along x, y and z coordinates (i.e. Px-Cx, Py-Cy and Pz-Cz)
    # Zeta - sum of the gaussian exponents on first and second centers

    # store is a 35x35 array where we would store primitive integral values. The anatomy is as follows: 
#      0     1       2       3       4        5         6        7        8        9        10        11        12        13        14        15        16        17        18        19        20         21         22         23         24          25         26        27         28         29         30         31         32         33         34  
# 0  [S|S], [S|Px], [S|Py], [S|Pz], [S|Dxy], [S|Dyz], [S|Dxz], [S|Dxx], [S|Dyy], [S|Dzz], [S|Fxyz], [S|Fxxy], [S|Fxyy], [S|Fxxz], [S|Fxzz], [S|Fyyz], [S|Fyzz], [S|Fxxx], [S|Fyyy], [S|Fzzz], [S|Gxxyy], [S|Gxxzz], [S|Gyyzz], [S|Gxxyz], [S|Gxyyz], [S|Gxyzz], [S|Gxxxz], [S|Gxzzz], [S|Gxxxy], [S|Gxyyy], [S|Gyyyz], [S|Gyzzz], [S|Gxxxx], [S|Gyyyy], [S|Gzzzz]
# 1  [Px|S], [Px|Px], [Px|Py], [Px|Pz], [Px|Dxy], [Px|Dyz], [Px|Dxz], [Px|Dxx], [Px|Dyy], [Px|Dzz], [Px|Fxyz], [Px|Fxxy], [Px|Fxyy], [Px|Fxxz], [Px|Fxzz], [Px|Fyyz], [Px|Fyzz], [Px|Fxxx], [Px|Fyyy], [Px|Fzzz], [Px|Gxxyy], [Px|Gxxzz], [Px|Gyyzz], [Px|Gxxyz], [Px|Gxyyz], [Px|Gxyzz], [Px|Gxxxz], [Px|Gxzzz], [Px|Gxxxy], [Px|Gxyyy], [Px|Gyyyz], [Px|Gyzzz], [Px|Gxxxx], [Px|Gyyyy], [Px|Gzzzz]
# 2  [Py|S], [Py|Px], [Py|Py], [Py|Pz], [Py|Dxy], [Py|Dyz], [Py|Dxz], [Py|Dxx], [Py|Dyy], [Py|Dzz], [Py|Fxyz], [Py|Fxxy], [Py|Fxyy], [Py|Fxxz], [Py|Fxzz], [Py|Fyyz], [Py|Fyzz], [Py|Fxxx], [Py|Fyyy], [Py|Fzzz], [Py|Gxxyy], [Py|Gxxzz], [Py|Gyyzz], [Py|Gxxyz], [Py|Gxyyz], [Py|Gxyzz], [Py|Gxxxz], [Py|Gxzzz], [Py|Gxxxy], [Py|Gxyyy], [Py|Gyyyz], [Py|Gyzzz], [Py|Gxxxx], [Py|Gyyyy], [Py|Gzzzz]
# 3  [Pz|S], [Pz|Px], [Pz|Py], [Pz|Pz], [Pz|Dxy], [Pz|Dyz], [Pz|Dxz], [Pz|Dxx], [Pz|Dyy], [Pz|Dzz], [Pz|Fxyz], [Pz|Fxxy], [Pz|Fxyy], [Pz|Fxxz], [Pz|Fxzz], [Pz|Fyyz], [Pz|Fyzz], [Pz|Fxxx], [Pz|Fyyy], [Pz|Fzzz], [Pz|Gxxyy], [Pz|Gxxzz], [Pz|Gyyzz], [Pz|Gxxyz], [Pz|Gxyyz], [Pz|Gxyzz], [Pz|Gxxxz], [Pz|Gxzzz], [Pz|Gxxxy], [Pz|Gxyyy], [Pz|Gyyyz], [Pz|Gyzzz], [Pz|Gxxxx], [Pz|Gyyyy], [Pz|Gzzzz]
# 4  [Dxy|S], [Dxy|Px], [Dxy|Py], [Dxy|Pz], [Dxy|Dxy], [Dxy|Dyz], [Dxy|Dxz], [Dxy|Dxx], [Dxy|Dyy], [Dxy|Dzz], [Dxy|Fxyz], [Dxy|Fxxy], [Dxy|Fxyy], [Dxy|Fxxz], [Dxy|Fxzz], [Dxy|Fyyz], [Dxy|Fyzz], [Dxy|Fxxx], [Dxy|Fyyy], [Dxy|Fzzz], [Dxy|Gxxyy], [Dxy|Gxxzz], [Dxy|Gyyzz], [Dxy|Gxxyz], [Dxy|Gxyyz], [Dxy|Gxyzz], [Dxy|Gxxxz], [Dxy|Gxzzz], [Dxy|Gxxxy], [Dxy|Gxyyy], [Dxy|Gyyyz], [Dxy|Gyzzz], [Dxy|Gxxxx], [Dxy|Gyyyy], [Dxy|Gzzzz]
# 5  [Dyz|S], [Dyz|Px], [Dyz|Py], [Dyz|Pz], [Dyz|Dxy], [Dyz|Dyz], [Dyz|Dxz], [Dyz|Dxx], [Dyz|Dyy], [Dyz|Dzz], [Dyz|Fxyz], [Dyz|Fxxy], [Dyz|Fxyy], [Dyz|Fxxz], [Dyz|Fxzz], [Dyz|Fyyz], [Dyz|Fyzz], [Dyz|Fxxx], [Dyz|Fyyy], [Dyz|Fzzz], [Dyz|Gxxyy], [Dyz|Gxxzz], [Dyz|Gyyzz], [Dyz|Gxxyz], [Dyz|Gxyyz], [Dyz|Gxyzz], [Dyz|Gxxxz], [Dyz|Gxzzz], [Dyz|Gxxxy], [Dyz|Gxyyy], [Dyz|Gyyyz], [Dyz|Gyzzz], [Dyz|Gxxxx], [Dyz|Gyyyy], [Dyz|Gzzzz]
# 6  [Dxz|S], [Dxz|Px], [Dxz|Py], [Dxz|Pz], [Dxz|Dxy], [Dxz|Dyz], [Dxz|Dxz], [Dxz|Dxx], [Dxz|Dyy], [Dxz|Dzz], [Dxz|Fxyz], [Dxz|Fxxy], [Dxz|Fxyy], [Dxz|Fxxz], [Dxz|Fxzz], [Dxz|Fyyz], [Dxz|Fyzz], [Dxz|Fxxx], [Dxz|Fyyy], [Dxz|Fzzz], [Dxz|Gxxyy], [Dxz|Gxxzz], [Dxz|Gyyzz], [Dxz|Gxxyz], [Dxz|Gxyyz], [Dxz|Gxyzz], [Dxz|Gxxxz], [Dxz|Gxzzz], [Dxz|Gxxxy], [Dxz|Gxyyy], [Dxz|Gyyyz], [Dxz|Gyzzz], [Dxz|Gxxxx], [Dxz|Gyyyy], [Dxz|Gzzzz]
# 7  [Dxx|S], [Dxx|Px], [Dxx|Py], [Dxx|Pz], [Dxx|Dxy], [Dxx|Dyz], [Dxx|Dxz], [Dxx|Dxx], [Dxx|Dyy], [Dxx|Dzz], [Dxx|Fxyz], [Dxx|Fxxy], [Dxx|Fxyy], [Dxx|Fxxz], [Dxx|Fxzz], [Dxx|Fyyz], [Dxx|Fyzz], [Dxx|Fxxx], [Dxx|Fyyy], [Dxx|Fzzz], [Dxx|Gxxyy], [Dxx|Gxxzz], [Dxx|Gyyzz], [Dxx|Gxxyz], [Dxx|Gxyyz], [Dxx|Gxyzz], [Dxx|Gxxxz], [Dxx|Gxzzz], [Dxx|Gxxxy], [Dxx|Gxyyy], [Dxx|Gyyyz], [Dxx|Gyzzz], [Dxx|Gxxxx], [Dxx|Gyyyy], [Dxx|Gzzzz]
# 8  [Dyy|S], [Dyy|Px], [Dyy|Py], [Dyy|Pz], [Dyy|Dxy], [Dyy|Dyz], [Dyy|Dxz], [Dyy|Dxx], [Dyy|Dyy], [Dyy|Dzz], [Dyy|Fxyz], [Dyy|Fxxy], [Dyy|Fxyy], [Dyy|Fxxz], [Dyy|Fxzz], [Dyy|Fyyz], [Dyy|Fyzz], [Dyy|Fxxx], [Dyy|Fyyy], [Dyy|Fzzz], [Dyy|Gxxyy], [Dyy|Gxxzz], [Dyy|Gyyzz], [Dyy|Gxxyz], [Dyy|Gxyyz], [Dyy|Gxyzz], [Dyy|Gxxxz], [Dyy|Gxzzz], [Dyy|Gxxxy], [Dyy|Gxyyy], [Dyy|Gyyyz], [Dyy|Gyzzz], [Dyy|Gxxxx], [Dyy|Gyyyy], [Dyy|Gzzzz]
# 9  [Dzz|S], [Dzz|Px], [Dzz|Py], [Dzz|Pz], [Dzz|Dxy], [Dzz|Dyz], [Dzz|Dxz], [Dzz|Dxx], [Dzz|Dyy], [Dzz|Dzz], [Dzz|Fxyz], [Dzz|Fxxy], [Dzz|Fxyy], [Dzz|Fxxz], [Dzz|Fxzz], [Dzz|Fyyz], [Dzz|Fyzz], [Dzz|Fxxx], [Dzz|Fyyy], [Dzz|Fzzz], [Dzz|Gxxyy], [Dzz|Gxxzz], [Dzz|Gyyzz], [Dzz|Gxxyz], [Dzz|Gxyyz], [Dzz|Gxyzz], [Dzz|Gxxxz], [Dzz|Gxzzz], [Dzz|Gxxxy], [Dzz|Gxyyy], [Dzz|Gyyyz], [Dzz|Gyzzz], [Dzz|Gxxxx], [Dzz|Gyyyy], [Dzz|Gzzzz]
# 10 [Fxyz|S], [Fxyz|Px], [Fxyz|Py], [Fxyz|Pz], [Fxyz|Dxy], [Fxyz|Dyz], [Fxyz|Dxz], [Fxyz|Dxx], [Fxyz|Dyy], [Fxyz|Dzz], [Fxyz|Fxyz], [Fxyz|Fxxy], [Fxyz|Fxyy], [Fxyz|Fxxz], [Fxyz|Fxzz], [Fxyz|Fyyz], [Fxyz|Fyzz], [Fxyz|Fxxx], [Fxyz|Fyyy], [Fxyz|Fzzz], [Fxyz|Gxxyy], [Fxyz|Gxxzz], [Fxyz|Gyyzz], [Fxyz|Gxxyz], [Fxyz|Gxyyz], [Fxyz|Gxyzz], [Fxyz|Gxxxz], [Fxyz|Gxzzz], [Fxyz|Gxxxy], [Fxyz|Gxyyy], [Fxyz|Gyyyz], [Fxyz|Gyzzz], [Fxyz|Gxxxx], [Fxyz|Gyyyy], [Fxyz|Gzzzz]
# 11 [Fxxy|S], [Fxxy|Px], [Fxxy|Py], [Fxxy|Pz], [Fxxy|Dxy], [Fxxy|Dyz], [Fxxy|Dxz], [Fxxy|Dxx], [Fxxy|Dyy], [Fxxy|Dzz], [Fxxy|Fxyz], [Fxxy|Fxxy], [Fxxy|Fxyy], [Fxxy|Fxxz], [Fxxy|Fxzz], [Fxxy|Fyyz], [Fxxy|Fyzz], [Fxxy|Fxxx], [Fxxy|Fyyy], [Fxxy|Fzzz], [Fxxy|Gxxyy], [Fxxy|Gxxzz], [Fxxy|Gyyzz], [Fxxy|Gxxyz], [Fxxy|Gxyyz], [Fxxy|Gxyzz], [Fxxy|Gxxxz], [Fxxy|Gxzzz], [Fxxy|Gxxxy], [Fxxy|Gxyyy], [Fxxy|Gyyyz], [Fxxy|Gyzzz], [Fxxy|Gxxxx], [Fxxy|Gyyyy], [Fxxy|Gzzzz]
# 12 [Fxyy|S], [Fxyy|Px], [Fxyy|Py], [Fxyy|Pz], [Fxyy|Dxy], [Fxyy|Dyz], [Fxyy|Dxz], [Fxyy|Dxx], [Fxyy|Dyy], [Fxyy|Dzz], [Fxyy|Fxyz], [Fxyy|Fxxy], [Fxyy|Fxyy], [Fxyy|Fxxz], [Fxyy|Fxzz], [Fxyy|Fyyz], [Fxyy|Fyzz], [Fxyy|Fxxx], [Fxyy|Fyyy], [Fxyy|Fzzz], [Fxyy|Gxxyy], [Fxyy|Gxxzz], [Fxyy|Gyyzz], [Fxyy|Gxxyz], [Fxyy|Gxyyz], [Fxyy|Gxyzz], [Fxyy|Gxxxz], [Fxyy|Gxzzz], [Fxyy|Gxxxy], [Fxyy|Gxyyy], [Fxyy|Gyyyz], [Fxyy|Gyzzz], [Fxyy|Gxxxx], [Fxyy|Gyyyy], [Fxyy|Gzzzz]
# 13 [Fxxz|S], [Fxxz|Px], [Fxxz|Py], [Fxxz|Pz], [Fxxz|Dxy], [Fxxz|Dyz], [Fxxz|Dxz], [Fxxz|Dxx], [Fxxz|Dyy], [Fxxz|Dzz], [Fxxz|Fxyz], [Fxxz|Fxxy], [Fxxz|Fxyy], [Fxxz|Fxxz], [Fxxz|Fxzz], [Fxxz|Fyyz], [Fxxz|Fyzz], [Fxxz|Fxxx], [Fxxz|Fyyy], [Fxxz|Fzzz], [Fxxz|Gxxyy], [Fxxz|Gxxzz], [Fxxz|Gyyzz], [Fxxz|Gxxyz], [Fxxz|Gxyyz], [Fxxz|Gxyzz], [Fxxz|Gxxxz], [Fxxz|Gxzzz], [Fxxz|Gxxxy], [Fxxz|Gxyyy], [Fxxz|Gyyyz], [Fxxz|Gyzzz], [Fxxz|Gxxxx], [Fxxz|Gyyyy], [Fxxz|Gzzzz]
# 14 [Fxzz|S], [Fxzz|Px], [Fxzz|Py], [Fxzz|Pz], [Fxzz|Dxy], [Fxzz|Dyz], [Fxzz|Dxz], [Fxzz|Dxx], [Fxzz|Dyy], [Fxzz|Dzz], [Fxzz|Fxyz], [Fxzz|Fxxy], [Fxzz|Fxyy], [Fxzz|Fxxz], [Fxzz|Fxzz], [Fxzz|Fyyz], [Fxzz|Fyzz], [Fxzz|Fxxx], [Fxzz|Fyyy], [Fxzz|Fzzz], [Fxzz|Gxxyy], [Fxzz|Gxxzz], [Fxzz|Gyyzz], [Fxzz|Gxxyz], [Fxzz|Gxyyz], [Fxzz|Gxyzz], [Fxzz|Gxxxz], [Fxzz|Gxzzz], [Fxzz|Gxxxy], [Fxzz|Gxyyy], [Fxzz|Gyyyz], [Fxzz|Gyzzz], [Fxzz|Gxxxx], [Fxzz|Gyyyy], [Fxzz|Gzzzz]
# 15 [Fyyz|S], [Fyyz|Px], [Fyyz|Py], [Fyyz|Pz], [Fyyz|Dxy], [Fyyz|Dyz], [Fyyz|Dxz], [Fyyz|Dxx], [Fyyz|Dyy], [Fyyz|Dzz], [Fyyz|Fxyz], [Fyyz|Fxxy], [Fyyz|Fxyy], [Fyyz|Fxxz], [Fyyz|Fxzz], [Fyyz|Fyyz], [Fyyz|Fyzz], [Fyyz|Fxxx], [Fyyz|Fyyy], [Fyyz|Fzzz], [Fyyz|Gxxyy], [Fyyz|Gxxzz], [Fyyz|Gyyzz], [Fyyz|Gxxyz], [Fyyz|Gxyyz], [Fyyz|Gxyzz], [Fyyz|Gxxxz], [Fyyz|Gxzzz], [Fyyz|Gxxxy], [Fyyz|Gxyyy], [Fyyz|Gyyyz], [Fyyz|Gyzzz], [Fyyz|Gxxxx], [Fyyz|Gyyyy], [Fyyz|Gzzzz]
# 16 [Fyzz|S], [Fyzz|Px], [Fyzz|Py], [Fyzz|Pz], [Fyzz|Dxy], [Fyzz|Dyz], [Fyzz|Dxz], [Fyzz|Dxx], [Fyzz|Dyy], [Fyzz|Dzz], [Fyzz|Fxyz], [Fyzz|Fxxy], [Fyzz|Fxyy], [Fyzz|Fxxz], [Fyzz|Fxzz], [Fyzz|Fyyz], [Fyzz|Fyzz], [Fyzz|Fxxx], [Fyzz|Fyyy], [Fyzz|Fzzz], [Fyzz|Gxxyy], [Fyzz|Gxxzz], [Fyzz|Gyyzz], [Fyzz|Gxxyz], [Fyzz|Gxyyz], [Fyzz|Gxyzz], [Fyzz|Gxxxz], [Fyzz|Gxzzz], [Fyzz|Gxxxy], [Fyzz|Gxyyy], [Fyzz|Gyyyz], [Fyzz|Gyzzz], [Fyzz|Gxxxx], [Fyzz|Gyyyy], [Fyzz|Gzzzz]
# 17 [Fxxx|S], [Fxxx|Px], [Fxxx|Py], [Fxxx|Pz], [Fxxx|Dxy], [Fxxx|Dyz], [Fxxx|Dxz], [Fxxx|Dxx], [Fxxx|Dyy], [Fxxx|Dzz], [Fxxx|Fxyz], [Fxxx|Fxxy], [Fxxx|Fxyy], [Fxxx|Fxxz], [Fxxx|Fxzz], [Fxxx|Fyyz], [Fxxx|Fyzz], [Fxxx|Fxxx], [Fxxx|Fyyy], [Fxxx|Fzzz], [Fxxx|Gxxyy], [Fxxx|Gxxzz], [Fxxx|Gyyzz], [Fxxx|Gxxyz], [Fxxx|Gxyyz], [Fxxx|Gxyzz], [Fxxx|Gxxxz], [Fxxx|Gxzzz], [Fxxx|Gxxxy], [Fxxx|Gxyyy], [Fxxx|Gyyyz], [Fxxx|Gyzzz], [Fxxx|Gxxxx], [Fxxx|Gyyyy], [Fxxx|Gzzzz]
# 18 [Fyyy|S], [Fyyy|Px], [Fyyy|Py], [Fyyy|Pz], [Fyyy|Dxy], [Fyyy|Dyz], [Fyyy|Dxz], [Fyyy|Dxx], [Fyyy|Dyy], [Fyyy|Dzz], [Fyyy|Fxyz], [Fyyy|Fxxy], [Fyyy|Fxyy], [Fyyy|Fxxz], [Fyyy|Fxzz], [Fyyy|Fyyz], [Fyyy|Fyzz], [Fyyy|Fxxx], [Fyyy|Fyyy], [Fyyy|Fzzz], [Fyyy|Gxxyy], [Fyyy|Gxxzz], [Fyyy|Gyyzz], [Fyyy|Gxxyz], [Fyyy|Gxyyz], [Fyyy|Gxyzz], [Fyyy|Gxxxz], [Fyyy|Gxzzz], [Fyyy|Gxxxy], [Fyyy|Gxyyy], [Fyyy|Gyyyz], [Fyyy|Gyzzz], [Fyyy|Gxxxx], [Fyyy|Gyyyy], [Fyyy|Gzzzz]
# 19 [Fzzz|S], [Fzzz|Px], [Fzzz|Py], [Fzzz|Pz], [Fzzz|Dxy], [Fzzz|Dyz], [Fzzz|Dxz], [Fzzz|Dxx], [Fzzz|Dyy], [Fzzz|Dzz], [Fzzz|Fxyz], [Fzzz|Fxxy], [Fzzz|Fxyy], [Fzzz|Fxxz], [Fzzz|Fxzz], [Fzzz|Fyyz], [Fzzz|Fyzz], [Fzzz|Fxxx], [Fzzz|Fyyy], [Fzzz|Fzzz], [Fzzz|Gxxyy], [Fzzz|Gxxzz], [Fzzz|Gyyzz], [Fzzz|Gxxyz], [Fzzz|Gxyyz], [Fzzz|Gxyzz], [Fzzz|Gxxxz], [Fzzz|Gxzzz], [Fzzz|Gxxxy], [Fzzz|Gxyyy], [Fzzz|Gyyyz], [Fzzz|Gyzzz], [Fzzz|Gxxxx], [Fzzz|Gyyyy], [Fzzz|Gzzzz]
# 20 [Gxxyy|S], [Gxxyy|Px], [Gxxyy|Py], [Gxxyy|Pz], [Gxxyy|Dxy], [Gxxyy|Dyz], [Gxxyy|Dxz], [Gxxyy|Dxx], [Gxxyy|Dyy], [Gxxyy|Dzz], [Gxxyy|Fxyz], [Gxxyy|Fxxy], [Gxxyy|Fxyy], [Gxxyy|Fxxz], [Gxxyy|Fxzz], [Gxxyy|Fyyz], [Gxxyy|Fyzz], [Gxxyy|Fxxx], [Gxxyy|Fyyy], [Gxxyy|Fzzz], [Gxxyy|Gxxyy], [Gxxyy|Gxxzz], [Gxxyy|Gyyzz], [Gxxyy|Gxxyz], [Gxxyy|Gxyyz], [Gxxyy|Gxyzz], [Gxxyy|Gxxxz], [Gxxyy|Gxzzz], [Gxxyy|Gxxxy], [Gxxyy|Gxyyy], [Gxxyy|Gyyyz], [Gxxyy|Gyzzz], [Gxxyy|Gxxxx], [Gxxyy|Gyyyy], [Gxxyy|Gzzzz]
# 21 [Gxxzz|S], [Gxxzz|Px], [Gxxzz|Py], [Gxxzz|Pz], [Gxxzz|Dxy], [Gxxzz|Dyz], [Gxxzz|Dxz], [Gxxzz|Dxx], [Gxxzz|Dyy], [Gxxzz|Dzz], [Gxxzz|Fxyz], [Gxxzz|Fxxy], [Gxxzz|Fxyy], [Gxxzz|Fxxz], [Gxxzz|Fxzz], [Gxxzz|Fyyz], [Gxxzz|Fyzz], [Gxxzz|Fxxx], [Gxxzz|Fyyy], [Gxxzz|Fzzz], [Gxxzz|Gxxyy], [Gxxzz|Gxxzz], [Gxxzz|Gyyzz], [Gxxzz|Gxxyz], [Gxxzz|Gxyyz], [Gxxzz|Gxyzz], [Gxxzz|Gxxxz], [Gxxzz|Gxzzz], [Gxxzz|Gxxxy], [Gxxzz|Gxyyy], [Gxxzz|Gyyyz], [Gxxzz|Gyzzz], [Gxxzz|Gxxxx], [Gxxzz|Gyyyy], [Gxxzz|Gzzzz]
# 22 [Gyyzz|S], [Gyyzz|Px], [Gyyzz|Py], [Gyyzz|Pz], [Gyyzz|Dxy], [Gyyzz|Dyz], [Gyyzz|Dxz], [Gyyzz|Dxx], [Gyyzz|Dyy], [Gyyzz|Dzz], [Gyyzz|Fxyz], [Gyyzz|Fxxy], [Gyyzz|Fxyy], [Gyyzz|Fxxz], [Gyyzz|Fxzz], [Gyyzz|Fyyz], [Gyyzz|Fyzz], [Gyyzz|Fxxx], [Gyyzz|Fyyy], [Gyyzz|Fzzz], [Gyyzz|Gxxyy], [Gyyzz|Gxxzz], [Gyyzz|Gyyzz], [Gyyzz|Gxxyz], [Gyyzz|Gxyyz], [Gyyzz|Gxyzz], [Gyyzz|Gxxxz], [Gyyzz|Gxzzz], [Gyyzz|Gxxxy], [Gyyzz|Gxyyy], [Gyyzz|Gyyyz], [Gyyzz|Gyzzz], [Gyyzz|Gxxxx], [Gyyzz|Gyyyy], [Gyyzz|Gzzzz]
# 23 [Gxxyz|S], [Gxxyz|Px], [Gxxyz|Py], [Gxxyz|Pz], [Gxxyz|Dxy], [Gxxyz|Dyz], [Gxxyz|Dxz], [Gxxyz|Dxx], [Gxxyz|Dyy], [Gxxyz|Dzz], [Gxxyz|Fxyz], [Gxxyz|Fxxy], [Gxxyz|Fxyy], [Gxxyz|Fxxz], [Gxxyz|Fxzz], [Gxxyz|Fyyz], [Gxxyz|Fyzz], [Gxxyz|Fxxx], [Gxxyz|Fyyy], [Gxxyz|Fzzz], [Gxxyz|Gxxyy], [Gxxyz|Gxxzz], [Gxxyz|Gyyzz], [Gxxyz|Gxxyz], [Gxxyz|Gxyyz], [Gxxyz|Gxyzz], [Gxxyz|Gxxxz], [Gxxyz|Gxzzz], [Gxxyz|Gxxxy], [Gxxyz|Gxyyy], [Gxxyz|Gyyyz], [Gxxyz|Gyzzz], [Gxxyz|Gxxxx], [Gxxyz|Gyyyy], [Gxxyz|Gzzzz]
# 24 [Gxyyz|S], [Gxyyz|Px], [Gxyyz|Py], [Gxyyz|Pz], [Gxyyz|Dxy], [Gxyyz|Dyz], [Gxyyz|Dxz], [Gxyyz|Dxx], [Gxyyz|Dyy], [Gxyyz|Dzz], [Gxyyz|Fxyz], [Gxyyz|Fxxy], [Gxyyz|Fxyy], [Gxyyz|Fxxz], [Gxyyz|Fxzz], [Gxyyz|Fyyz], [Gxyyz|Fyzz], [Gxyyz|Fxxx], [Gxyyz|Fyyy], [Gxyyz|Fzzz], [Gxyyz|Gxxyy], [Gxyyz|Gxxzz], [Gxyyz|Gyyzz], [Gxyyz|Gxxyz], [Gxyyz|Gxyyz], [Gxyyz|Gxyzz], [Gxyyz|Gxxxz], [Gxyyz|Gxzzz], [Gxyyz|Gxxxy], [Gxyyz|Gxyyy], [Gxyyz|Gyyyz], [Gxyyz|Gyzzz], [Gxyyz|Gxxxx], [Gxyyz|Gyyyy], [Gxyyz|Gzzzz]
# 25 [Gxyzz|S], [Gxyzz|Px], [Gxyzz|Py], [Gxyzz|Pz], [Gxyzz|Dxy], [Gxyzz|Dyz], [Gxyzz|Dxz], [Gxyzz|Dxx], [Gxyzz|Dyy], [Gxyzz|Dzz], [Gxyzz|Fxyz], [Gxyzz|Fxxy], [Gxyzz|Fxyy], [Gxyzz|Fxxz], [Gxyzz|Fxzz], [Gxyzz|Fyyz], [Gxyzz|Fyzz], [Gxyzz|Fxxx], [Gxyzz|Fyyy], [Gxyzz|Fzzz], [Gxyzz|Gxxyy], [Gxyzz|Gxxzz], [Gxyzz|Gyyzz], [Gxyzz|Gxxyz], [Gxyzz|Gxyyz], [Gxyzz|Gxyzz], [Gxyzz|Gxxxz], [Gxyzz|Gxzzz], [Gxyzz|Gxxxy], [Gxyzz|Gxyyy], [Gxyzz|Gyyyz], [Gxyzz|Gyzzz], [Gxyzz|Gxxxx], [Gxyzz|Gyyyy], [Gxyzz|Gzzzz]
# 26 [Gxxxz|S], [Gxxxz|Px], [Gxxxz|Py], [Gxxxz|Pz], [Gxxxz|Dxy], [Gxxxz|Dyz], [Gxxxz|Dxz], [Gxxxz|Dxx], [Gxxxz|Dyy], [Gxxxz|Dzz], [Gxxxz|Fxyz], [Gxxxz|Fxxy], [Gxxxz|Fxyy], [Gxxxz|Fxxz], [Gxxxz|Fxzz], [Gxxxz|Fyyz], [Gxxxz|Fyzz], [Gxxxz|Fxxx], [Gxxxz|Fyyy], [Gxxxz|Fzzz], [Gxxxz|Gxxyy], [Gxxxz|Gxxzz], [Gxxxz|Gyyzz], [Gxxxz|Gxxyz], [Gxxxz|Gxyyz], [Gxxxz|Gxyzz], [Gxxxz|Gxxxz], [Gxxxz|Gxzzz], [Gxxxz|Gxxxy], [Gxxxz|Gxyyy], [Gxxxz|Gyyyz], [Gxxxz|Gyzzz], [Gxxxz|Gxxxx], [Gxxxz|Gyyyy], [Gxxxz|Gzzzz]
# 27 [Gxzzz|S], [Gxzzz|Px], [Gxzzz|Py], [Gxzzz|Pz], [Gxzzz|Dxy], [Gxzzz|Dyz], [Gxzzz|Dxz], [Gxzzz|Dxx], [Gxzzz|Dyy], [Gxzzz|Dzz], [Gxzzz|Fxyz], [Gxzzz|Fxxy], [Gxzzz|Fxyy], [Gxzzz|Fxxz], [Gxzzz|Fxzz], [Gxzzz|Fyyz], [Gxzzz|Fyzz], [Gxzzz|Fxxx], [Gxzzz|Fyyy], [Gxzzz|Fzzz], [Gxzzz|Gxxyy], [Gxzzz|Gxxzz], [Gxzzz|Gyyzz], [Gxzzz|Gxxyz], [Gxzzz|Gxyyz], [Gxzzz|Gxyzz], [Gxzzz|Gxxxz], [Gxzzz|Gxzzz], [Gxzzz|Gxxxy], [Gxzzz|Gxyyy], [Gxzzz|Gyyyz], [Gxzzz|Gyzzz], [Gxzzz|Gxxxx], [Gxzzz|Gyyyy], [Gxzzz|Gzzzz]
# 28 [Gxxxy|S], [Gxxxy|Px], [Gxxxy|Py], [Gxxxy|Pz], [Gxxxy|Dxy], [Gxxxy|Dyz], [Gxxxy|Dxz], [Gxxxy|Dxx], [Gxxxy|Dyy], [Gxxxy|Dzz], [Gxxxy|Fxyz], [Gxxxy|Fxxy], [Gxxxy|Fxyy], [Gxxxy|Fxxz], [Gxxxy|Fxzz], [Gxxxy|Fyyz], [Gxxxy|Fyzz], [Gxxxy|Fxxx], [Gxxxy|Fyyy], [Gxxxy|Fzzz], [Gxxxy|Gxxyy], [Gxxxy|Gxxzz], [Gxxxy|Gyyzz], [Gxxxy|Gxxyz], [Gxxxy|Gxyyz], [Gxxxy|Gxyzz], [Gxxxy|Gxxxz], [Gxxxy|Gxzzz], [Gxxxy|Gxxxy], [Gxxxy|Gxyyy], [Gxxxy|Gyyyz], [Gxxxy|Gyzzz], [Gxxxy|Gxxxx], [Gxxxy|Gyyyy], [Gxxxy|Gzzzz]
# 29 [Gxyyy|S], [Gxyyy|Px], [Gxyyy|Py], [Gxyyy|Pz], [Gxyyy|Dxy], [Gxyyy|Dyz], [Gxyyy|Dxz], [Gxyyy|Dxx], [Gxyyy|Dyy], [Gxyyy|Dzz], [Gxyyy|Fxyz], [Gxyyy|Fxxy], [Gxyyy|Fxyy], [Gxyyy|Fxxz], [Gxyyy|Fxzz], [Gxyyy|Fyyz], [Gxyyy|Fyzz], [Gxyyy|Fxxx], [Gxyyy|Fyyy], [Gxyyy|Fzzz], [Gxyyy|Gxxyy], [Gxyyy|Gxxzz], [Gxyyy|Gyyzz], [Gxyyy|Gxxyz], [Gxyyy|Gxyyz], [Gxyyy|Gxyzz], [Gxyyy|Gxxxz], [Gxyyy|Gxzzz], [Gxyyy|Gxxxy], [Gxyyy|Gxyyy], [Gxyyy|Gyyyz], [Gxyyy|Gyzzz], [Gxyyy|Gxxxx], [Gxyyy|Gyyyy], [Gxyyy|Gzzzz]
# 30 [Gyyyz|S], [Gyyyz|Px], [Gyyyz|Py], [Gyyyz|Pz], [Gyyyz|Dxy], [Gyyyz|Dyz], [Gyyyz|Dxz], [Gyyyz|Dxx], [Gyyyz|Dyy], [Gyyyz|Dzz], [Gyyyz|Fxyz], [Gyyyz|Fxxy], [Gyyyz|Fxyy], [Gyyyz|Fxxz], [Gyyyz|Fxzz], [Gyyyz|Fyyz], [Gyyyz|Fyzz], [Gyyyz|Fxxx], [Gyyyz|Fyyy], [Gyyyz|Fzzz], [Gyyyz|Gxxyy], [Gyyyz|Gxxzz], [Gyyyz|Gyyzz], [Gyyyz|Gxxyz], [Gyyyz|Gxyyz], [Gyyyz|Gxyzz], [Gyyyz|Gxxxz], [Gyyyz|Gxzzz], [Gyyyz|Gxxxy], [Gyyyz|Gxyyy], [Gyyyz|Gyyyz], [Gyyyz|Gyzzz], [Gyyyz|Gxxxx], [Gyyyz|Gyyyy], [Gyyyz|Gzzzz]
# 31 [Gyzzz|S], [Gyzzz|Px], [Gyzzz|Py], [Gyzzz|Pz], [Gyzzz|Dxy], [Gyzzz|Dyz], [Gyzzz|Dxz], [Gyzzz|Dxx], [Gyzzz|Dyy], [Gyzzz|Dzz], [Gyzzz|Fxyz], [Gyzzz|Fxxy], [Gyzzz|Fxyy], [Gyzzz|Fxxz], [Gyzzz|Fxzz], [Gyzzz|Fyyz], [Gyzzz|Fyzz], [Gyzzz|Fxxx], [Gyzzz|Fyyy], [Gyzzz|Fzzz], [Gyzzz|Gxxyy], [Gyzzz|Gxxzz], [Gyzzz|Gyyzz], [Gyzzz|Gxxyz], [Gyzzz|Gxyyz], [Gyzzz|Gxyzz], [Gyzzz|Gxxxz], [Gyzzz|Gxzzz], [Gyzzz|Gxxxy], [Gyzzz|Gxyyy], [Gyzzz|Gyyyz], [Gyzzz|Gyzzz], [Gyzzz|Gxxxx], [Gyzzz|Gyyyy], [Gyzzz|Gzzzz]
# 32 [Gxxxx|S], [Gxxxx|Px], [Gxxxx|Py], [Gxxxx|Pz], [Gxxxx|Dxy], [Gxxxx|Dyz], [Gxxxx|Dxz], [Gxxxx|Dxx], [Gxxxx|Dyy], [Gxxxx|Dzz], [Gxxxx|Fxyz], [Gxxxx|Fxxy], [Gxxxx|Fxyy], [Gxxxx|Fxxz], [Gxxxx|Fxzz], [Gxxxx|Fyyz], [Gxxxx|Fyzz], [Gxxxx|Fxxx], [Gxxxx|Fyyy], [Gxxxx|Fzzz], [Gxxxx|Gxxyy], [Gxxxx|Gxxzz], [Gxxxx|Gyyzz], [Gxxxx|Gxxyz], [Gxxxx|Gxyyz], [Gxxxx|Gxyzz], [Gxxxx|Gxxxz], [Gxxxx|Gxzzz], [Gxxxx|Gxxxy], [Gxxxx|Gxyyy], [Gxxxx|Gyyyz], [Gxxxx|Gyzzz], [Gxxxx|Gxxxx], [Gxxxx|Gyyyy], [Gxxxx|Gzzzz]
# 33 [Gyyyy|S], [Gyyyy|Px], [Gyyyy|Py], [Gyyyy|Pz], [Gyyyy|Dxy], [Gyyyy|Dyz], [Gyyyy|Dxz], [Gyyyy|Dxx], [Gyyyy|Dyy], [Gyyyy|Dzz], [Gyyyy|Fxyz], [Gyyyy|Fxxy], [Gyyyy|Fxyy], [Gyyyy|Fxxz], [Gyyyy|Fxzz], [Gyyyy|Fyyz], [Gyyyy|Fyzz], [Gyyyy|Fxxx], [Gyyyy|Fyyy], [Gyyyy|Fzzz], [Gyyyy|Gxxyy], [Gyyyy|Gxxzz], [Gyyyy|Gyyzz], [Gyyyy|Gxxyz], [Gyyyy|Gxyyz], [Gyyyy|Gxyzz], [Gyyyy|Gxxxz], [Gyyyy|Gxzzz], [Gyyyy|Gxxxy], [Gyyyy|Gxyyy], [Gyyyy|Gyyyz], [Gyyyy|Gyzzz], [Gyyyy|Gxxxx], [Gyyyy|Gyyyy], [Gyyyy|Gzzzz]
# 34 [Gzzzz|S], [Gzzzz|Px], [Gzzzz|Py], [Gzzzz|Pz], [Gzzzz|Dxy], [Gzzzz|Dyz], [Gzzzz|Dxz], [Gzzzz|Dxx], [Gzzzz|Dyy], [Gzzzz|Dzz], [Gzzzz|Fxyz], [Gzzzz|Fxxy], [Gzzzz|Fxyy], [Gzzzz|Fxxz], [Gzzzz|Fxzz], [Gzzzz|Fyyz], [Gzzzz|Fyzz], [Gzzzz|Fxxx], [Gzzzz|Fyyy], [Gzzzz|Fzzz], [Gzzzz|Gxxyy], [Gzzzz|Gxxzz], [Gzzzz|Gyyzz], [Gzzzz|Gxxyz], [Gzzzz|Gxyyz], [Gzzzz|Gxyzz], [Gzzzz|Gxxxz], [Gzzzz|Gxzzz], [Gzzzz|Gxxxy], [Gzzzz|Gxyyy], [Gzzzz|Gyyyz], [Gzzzz|Gyzzz], [Gzzzz|Gxxxx], [Gzzzz|Gyyyy], [Gzzzz|Gzzzz] 

    # Yvertical - a 1D array holding boys function values,the size is Max_I+ Max_J+2.

    OEint.fha.write("%s void oei_vertical(int I, int J, int II, int JJ,QUICKDouble PAx, QUICKDouble PAy, QUICKDouble PAz,\n\
        QUICKDouble PBx, QUICKDouble PBy, QUICKDouble PBz, QUICKDouble PCx, QUICKDouble PCy, QUICKDouble PCz, QUICKDouble Zeta,\n\
        QUICKDouble* store, QUICKDouble* YVerticalTemp){ \n" % (func_qualifier))
    ss.save_int()
    ps.save_int()
    sp.save_int()
    pp.save_int()
    ds.save_int()
    sd.save_int()
    dp.save_int()
    pd.save_int()
    dd.save_int()
    fs.save_int()
    sf.save_int()
    fp.save_int()
    pf.save_int()
    fd.save_int()
    df.save_int()
    ff.save_int()

    OEint.fha.write("\n } \n")

    # Write the driver to save computed primitive integrals required for gradient calculation. The parameters are the same that we reported above. 
    OEint.fhga.write("%s void oei_grad_vertical(int I, int J, int II, int JJ,QUICKDouble PAx, QUICKDouble PAy, QUICKDouble PAz,\n\
        QUICKDouble PBx, QUICKDouble PBy, QUICKDouble PBz, QUICKDouble PCx, QUICKDouble PCy, QUICKDouble PCz, QUICKDouble Zeta,\n\
        QUICKDouble* store2, QUICKDouble* YVerticalTemp){ \n" % (func_qualifier))

    ss.save_int_grad()
    sp.save_int_grad()
    ps.save_int_grad()
    pp.save_int_grad()
    sd.save_int_grad()
    ds.save_int_grad()
    pd.save_int_grad()
    dp.save_int_grad()
    dd.save_int_grad()

    OEint.fhga.write("\n } \n") 

    OEint.fhc.close()
    OEint.fhd.close()
    OEint.fha.close()
    OEint.fhga.close()


