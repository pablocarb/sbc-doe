'''Optional and external configuration'''
import sys
import os
sys.path.append('/mnt/SBC1/code/sbc-api')

ICE_INFO = os.path.join(os.getenv('HOME'), 'info/ice.json')
ICE_USER = 'Pablo Carbonell'
ICE_EMAIL = 'pablo.carbonell@manchester.ac.uk'
ICE_RETRIEVE = 'ice' # server for queries
ICE_POST = 'icetest' # server for registration
