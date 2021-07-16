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

def write_oei(outdir):

    # set files
    OEint.fhc = open(outdir+"/gpu_oei_classes.h",'w')
    OEint.fhd = open(outdir+"/gpu_oei_definitions.h",'w')
    OEint.fha= open(outdir+"/gpu_oei_assembler.h",'w')

    # write license info
    file_handler.write_license(OEint.fhc)
    file_handler.write_license(OEint.fhd)
    file_handler.write_license(OEint.fha)

    # generate integral classes for systems containing only s, p and d functions.
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

    # write driver to use classes and save computed primitive integrals
    OEint.fha.write("__device__ __inline__ void OEint_vertical(int I, int J, int II, int JJ,QUICKDouble PAx, QUICKDouble PAy, QUICKDouble PAz,\n\
        QUICKDouble PBx, QUICKDouble PBy, QUICKDouble PBz, QUICKDouble PCx, QUICKDouble PCy, QUICKDouble PCz, QUICKDouble Zeta,\n\
        QUICKDouble* store, QUICKDouble* YVerticalTemp){ \n")
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

    OEint.fhc.close()
    OEint.fhd.close()
    OEint.fha.close()

