import requests
import json

def fetch(idtype='*'):
    headers = {'Content-Type': 'application/json', 'Accept': 'application/json'}
    data = {'type': idtype}
    response = requests.post("https://identity.synbiochem.co.uk/fetch/", data=json.dumps(data), headers=headers)
    return json.loads(response.text)

def register(idtype, idnum, who=None, description=None):
    headers = {'Content-Type': 'application/json', 'Accept': 'application/json'}
    data = {'type': idtype, 'id': idnum, 'who': who, 'description': description}
    response = requests.post("https://identity.synbiochem.co.uk/register/", data=json.dumps(data), headers=headers)
    return json.loads(response.text)

def reserve(idtype, blocksize, who=None, description=None):
    headers = {'Content-Type': 'application/json', 'Accept': 'application/json'}
    data = {'type': idtype, 'blocksize': blocksize, 'who': who, 'description': description}
    response = requests.post("https://identity.synbiochem.co.uk/reserve/", data=json.dumps(data), headers=headers)
    return json.loads(response.text)
