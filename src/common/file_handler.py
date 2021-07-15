#!---------------------------------------------------------------------!
#! Written by Madu Manathunga on 07/01/2021                            !
#!                                                                     !
#! Copyright (C) 2020-2021 Merz lab                                    !
#! Copyright (C) 2020-2021 Götz lab                                    !
#!                                                                     !
#! This Source Code Form is subject to the terms of the Mozilla Public !
#! License, v. 2.0. If a copy of the MPL was not distributed with this !
#! file, You can obtain one at http://mozilla.org/MPL/2.0/.            !
#!_____________________________________________________________________!

#!---------------------------------------------------------------------!
#! This source file contains classes for handling output files.        !
#!---------------------------------------------------------------------!

from datetime import date

class fhandler:
    def __init__(self, fname):
        self.fname=fname
        pass

def write_license(fh):
    today = date.today()
    dt = today.strftime("%d/%m/%Y")
    fh.write("/*\n")
    fh.write(" !---------------------------------------------------------------------!\n")
    fh.write(" ! Written by QUICK-GenInt code generator on %s                !\n" % (dt))
    fh.write(" !                                                                     !\n")
    fh.write(" ! Copyright (C) 2020-2021 Merz lab                                    !\n")
    fh.write(" ! Copyright (C) 2020-2021 Götz lab                                    !\n")
    fh.write(" !                                                                     !\n")
    fh.write(" ! This Source Code Form is subject to the terms of the Mozilla Public !\n")
    fh.write(" ! License, v. 2.0. If a copy of the MPL was not distributed with this !\n")
    fh.write(" ! file, You can obtain one at http://mozilla.org/MPL/2.0/.            !\n")
    fh.write(" !_____________________________________________________________________!\n")
    fh.write("*/\n\n")