from cps_fcns import *
from archive_fcns import *
from cntrd import *

list1=genfromtxt('workfile.csv',dtype=str,delimiter=',')
for row in list1:
    print 'Startin epoch %s...' % row[5]
    if 'cps' in row[0]:
        cps_ccf1,cps_ccf2,cps_ccf12,cps_images= \
                    cps_shift_estimate(row[0],row[14],row[15])
##        cpswave=cps_images['science']['w'][16]
##        cpsflux=cps_images['science']['f'][16]
        cps_R,cps_m,cps_fits,cps_vp,cps_vs,cps_cntr= \
                cps_todcor(cps_ccf1,cps_ccf2,cps_ccf12,row[6],row[7],cps_images)
    else:
        fpath=row[0].split('/')
        arc_ccf1,arc_ccf2,arc_ccf12,arc_images,orders=archive_shift_estimate(\
            fpath[1],fpath[-1][:18],row[14],row[15])
##        arcwave=arc_images['science']['w'][22]
##        arcflux=arc_images['science']['f'][22]
        arc_R,arc_m,arc_fits,arc_vp,arc_vs,arc_cntr=archive_todcor(arc_ccf1, \
            arc_ccf2,arc_ccf12,row[6],row[7],arc_images,orders)
        

