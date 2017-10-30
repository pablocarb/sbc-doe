'''Optional and external configuration'''
import sys
import os
sys.path.append('/mnt/SBC1/code/ice-utils')

ICE_INFO = os.path.join(os.getenv('HOME'), 'info/ice2.json')
ICE_USER = 'Pablo Carbonell'
ICE_EMAIL = 'pablo.carbonell@manchester.ac.uk'
ICE_SERVER = "https://ice.synbiochem.co.uk"
ICE_SERVER_TEST =  "https://testice.synbiochem.co.uk:8443/rest"
ICE_RETRIEVE = 'ice' # server for queries
ICE_POST = 'icetest' # server for registration

