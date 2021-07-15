
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
#! This is the main source file of QUICK integral code generator.      !
#!---------------------------------------------------------------------!

import sys
import os

# Get the absolute path of oei directory and set module search path
GenIntHome=os.path.abspath(os.getcwd())
outdir_name="output"
outdir=os.path.join(GenIntHome,outdir_name)

try:
    os.mkdir(outdir)
except OSError as error:
    print(error)

#common_abs_path=os.path.abspath(os.getcwd())+"/src/common"
#oei_abs_path=os.path.abspath(os.getcwd())+"/src/oei"
#sys.path.append(common_abs_path)
#sys.path.append(oei_abs_path)

#print(sys.path)

#import src.common.params as params
#import src.oei.file_handler as file_handler
import src.oei.one_electron_integral as one_electron_integral

#from one_electron_integral import *
#from params import *

# generate one electron integral source code
one_electron_integral.write_oei(outdir)

