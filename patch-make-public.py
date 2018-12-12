from glob import glob
from os import path
import re
import copy
import doeopt
import iceutils

pwd = '/mnt/syno/shared/Designs'



ice = iceutils.ICESession(doeopt.ICE_POST, info=doeopt.ICE_INFO) # icetest
for d in ['SBCDE00001', 'SBCDE00002', 'SBCDE00003']:
    for l in open(path.join(pwd, d, d+'.oadi0')):
        m = l.split('\t')
        part = int(re.sub('SBC', '', m[0]))
        group_number = 2
        print 'Part', part, 'Adding permissions group', group_number
        ice.add_write_permission(part, group_number)
    
