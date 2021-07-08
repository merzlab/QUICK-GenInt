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
#! This source file contains parameters required for generating VRR    !
#! code paths.                                                         !
#!---------------------------------------------------------------------!

# The following arrays come from src/modules/quick_params_module. In the host OEI/ERI
# code they serve as maps to store/load integrals from VRR. 3D arrays are implemented
# as 3D tuples in python. 

# In the trans array, rows represent the x dimension, so below we have 8 rows, 0-7. The column of each row represent y
# dimension. For eg, 0th row has 8 columns and 7th row has 1 column. Elements encapsulated by brackets
# represent z dimension. For eg, there are 8 z components in row 0, column 0 and just 1 component in
# row 7, column 0.  
trans=(\
    ((1,4,10,20,35,56,84,120),(3,6,17,32,48,67,100),(9,16,23,42,73,106),(19,31,43,79,112),(34,49,74,113),(55,68,107),(83,101),(119,)),\
    ((2,7,15,28,50,69,102),(5,11,26,41,59,87),(13,25,36,60,88),(30,40,61,94),(52,58,89),(71,86),(104,)),\
    ((8,14,22,44,75,108),(12,24,37,62,90),(21,38,66,99),(46,64,98),(77,92),(110,)),\
    ((18,27,45,80,114),(29,39,63,95),(47,65,97),(81,96),(116,)),\
    ((33,51,76,115),(53,57,91),(78,93),(117,)),\
    ((54,70,109),(72,85),(111,)),\
    ((82,103),(105,)),\
    ((118,),)\
    )

# The following 3x120 2D array is used to indicate the cartesian components of basis functions. The columns of each row stand for
# X, Y and Z cartesian coodinates. Therefore, row 0 indicates S, row 1-3 indicate
# P ( Px=(1,0,0), Py=(0,1,0), Pz=(0,0,1)), row 4-9 indicate D (where Dxy=(1,1,0), Dyz=(0,1,1), Dxz=(1,0,1), Dxx=(2,0,0), Dyy=(0,2,0),
# Dzz=(0,0,2)), row 10-19 indicate F and row 20-34 indicate G.
   
Mcal=((0,0,0),\
    (1,0,0), (0,1,0), (0,0,1),\
    (1,1,0), (0,1,1), (1,0,1), (2,0,0), (0,2,0), (0,0,2),\
    (1,1,1), (2,1,0), (1,2,0), (2,0,1), (1,0,2), (0,2,1), (0,1,2), (3,0,0), (0,3,0), (0,0,3),\
    (2,2,0), (2,0,2), (0,2,2), (2,1,1), (1,2,1), (1,1,2), (3,0,1), (1,0,3), (3,1,0), (1,3,0), (0,3,1), (0,1,3), (4,0,0), (0,4,0), (0,0,4)\
    )

# compare stored values with fortran code, print in the same order
def print_trans():

    idx_x=(\
        0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,\
        1,1,1,1,1,2,2,2,2,2,2,3,3,3,4,1,2,2,3,1,\
        1,0,0,2,3,2,3,0,0,1,4,1,4,5,0,0,4,1,1,1,\
        1,2,3,2,3,2,0,0,1,5,1,5,0,0,2,4,2,4,0,3,\
        3,6,0,0,5,1,1,1,1,2,4,2,4,1,3,3,3,2,2,0,\
        0,1,6,1,6,0,0,2,5,2,5,0,0,3,4,3,4,7,0,0 \
        )

    idx_y=(\
        0,0,0,0,0,1,1,1,1,2,2,2,3,3,4,0,0,0,0,1,\
        1,1,2,2,3,0,0,0,1,1,2,0,0,1,0,2,1,2,1,3,\
        1,2,3,0,0,3,2,1,4,0,0,4,1,0,5,0,1,4,1,2,\
        3,1,1,3,2,2,1,5,0,0,5,1,2,4,0,0,4,2,3,0,\
        3,0,6,0,1,5,1,2,4,1,1,4,2,3,1,3,2,3,2,1,\
        6,0,0,6,1,2,5,0,0,5,2,3,4,0,0,4,3,0,7,0 \
        )

    idx_z=(\
        0,1,2,3,4,0,1,2,3,0,1,2,0,1,0,0,1,2,3,0,\
        1,2,0,1,0,0,1,2,0,1,0,0,1,0,0,2,2,1,1,1,\
        3,3,2,3,2,0,0,4,1,4,1,0,0,0,0,5,1,1,4,3,\
        2,3,2,1,1,2,5,1,5,1,0,0,4,2,4,2,0,0,3,3,\
        0,0,0,6,1,1,5,4,2,4,2,1,1,3,3,1,2,2,3,6,\
        1,6,1,0,0,5,2,5,2,0,0,4,3,4,3,0,0,0,0,7 \
        )
    for i in range(0,len(idx_x)):
        print("trans[%d][%d][%d] = %d" % (idx_x[i],idx_y[i],idx_z[i],trans[idx_x[i]][idx_y[i]][idx_z[i]]))

def print_mcal():
    for i in range(0, len(Mcal)):
        for j in range(0, len(Mcal[i])):
            print("Mcal[%d][%d] = %d" % (i, j, Mcal[i][j]))

#print_trans()
#print_mcal()
