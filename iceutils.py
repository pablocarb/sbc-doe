import requests
import json

#BASE_URL = "https://ice.synbiochem.co.uk/rest"
#BASE_URL = "https://testice-env.eu-west-1.elasticbeanstalk.com:8443/rest"
#BASE_URL = "https://testice.synbiochem.co.uk:8443/rest"

class ICESession:
    def __init__(self, Server, info='/home/pablo/info/ice.json'):
        self.Server = Server
        self.info = info
        self.HEADERS = {'Accept' : 'application/json',
                        'Content-Type' : 'application/json'}
        self.SESSION_KEY = 'X-ICE-Authentication-SessionId'
        self.session_id = self.get_token()

    def credentials(self):
        x = json.load(open(self.info))[self.Server]
        self.BASE_URL = x['server']+'/rest'
        return {'email': x['email'], 'password': x['password']}

    def get_token(self):
        if self.SESSION_KEY in self.HEADERS:
            return self.HEADERS[SESSION_KEY]
        
        payload = self.credentials()
        response = requests.post(self.BASE_URL + '/accesstoken', data = json.dumps(payload),
                                 headers = self.HEADERS, verify=False)
        if response.status_code == 200:
            response_json = json.loads (response.text)
            session_id = response_json['sessionId']
            self.HEADERS[self.SESSION_KEY] = session_id
            return session_id
        return response.status_code


    def template_newpart(self, name, description,
                         creator, owner, investigator,
                         creatorEmail=None, ownerEmail=None, investigatorEmail=None):
        partData = {u'accessPermissions': [],
                    u'basePairCount': 0,
                    u'bioSafetyLevel': 1,
                    u'canEdit': True,
                    u'creator': creator,
                    u'creatorEmail': creatorEmail,
                    u'featureCount': 0,
                    u'fundingSource': u'BBSRC',
                    u'hasAttachment': False,
                    u'hasOriginalSequence': False,
                    u'hasSample': False,
                    u'hasSequence': False,
                    u'index': 0,
                    u'linkedParts': [],
                    u'links': [],
                    u'name': name,
                    u'owner': owner,
                    u'ownerEmail': ownerEmail,
                    u'parameters': [],
                    u'parents': [],
                    u'principalInvestigator': investigator,
                    u'principalInvestigatorEmail': investigatorEmail,
                    u'publicRead': False,
                    u'selectionMarkers': [],
                    u'shortDescription': description,
                    u'status': u'Complete',
                    u'type': u'PART',
                    u'visible': u'OK'
        }
        return partData


    #                u'principalInvestigatorId': 0,
    #                u'viewCount': 0,
    #                u'recordId': u'dc7c067f-44f4-4b32-a6cd-0d3a1a8262c9',
    #                u'partId': u'TEST_000002',
    #                u'creationTime': 1457021805938,
    #                u'modificationTime': 1457021805938,
    #                u'id': 2,
    #                u'ownerId': 1,
    #                u'creatorId': 0,

    def template_newplasmid(self, name, description,
                            creator, owner, investigator,
                            creatorEmail=None, ownerEmail=None, investigatorEmail=None):
        partData = {u'accessPermissions': [],
                    u'basePairCount': 0,
                    u'bioSafetyLevel': 1,
                    u'canEdit': True,
                    u'creator': creator,
                    u'creatorEmail': creatorEmail,
                    u'featureCount': 0,
                    u'fundingSource': u'BBSRC',
                    u'hasAttachment': False,
                    u'hasOriginalSequence': False,
                    u'hasSample': False,
                    u'hasSequence': False,
                    u'index': 0,
                    u'linkedParts': [],
                    u'links': [],
                    u'name': name,
                    u'owner': owner,
                    u'ownerEmail': ownerEmail,
                    u'parameters': [],
                    u'parents': [],
                    u'plasmidData': {u'circular': True},
                    u'principalInvestigator': investigator,
                    u'principalInvestigatorEmail': investigatorEmail,
                    u'publicRead': False,
                    u'selectionMarkers': [],
                    u'shortDescription': description,
                    u'status': u'Planned',
                    u'type': u'PLASMID',
                    u'visible': u'OK'
        }
        return partData
        #                u'creationTime': 1457088861013,
        #                u'recordId': u'3e61b464-8edf-4b45-b525-c68c05b0da74',
        #                u'creatorId': 0,
        #                u'id': 6,
        #                u'modificationTime': 1457088861013,
        #                u'viewCount': 0,
        #                u'partId': u'TEST_000006',
        #                u'ownerId': 1,
        #                u'principalInvestigatorId': 0,




    def get_part(self, part_number):
        url = self.BASE_URL + '/parts/' + str(part_number)
        response = requests.get(url, headers = self.HEADERS, verify=False)
        if response.status_code == 200:
            response_json = json.loads (response.text)
            return response_json
        return response.status_code

    def update_part(self, part_number, new_dict):
        url = self.BASE_URL + '/parts/' + str(part_number)
        response = requests.put(url, data=json.dumps(new_dict), headers = self.HEADERS, verify=False)
        if response.status_code == 200:
            response_json = json.loads (response.text)
            return response_json
        return response.status_code

    def rename_part(self, part_number):
        old_dict = get_part(part_number)
        if 'alias' in old_dict:
            old_name = old_dict['name']
            old_alias = old_dict['alias']
            old_references = ''
            if 'references' in old_dict:
                old_references = old_dict['references']

            new_name = old_alias
            new_alias = ''
            new_references = old_references + '\n\n' + old_name

            new_dict = old_dict.copy()
            new_dict['name'] = new_name
            new_dict['alias'] = new_alias
            new_dict['references'] = new_references

            return update_part(part_number, new_dict)


    def create_part(self, partData):
        url = self.BASE_URL + '/parts'
        response = requests.post(url, data=json.dumps(partData),
                                 headers = self.HEADERS, verify=False)
        if response.status_code == 200:
            response_json = json.loads (response.text)
            return response_json
        return response.status_code

    def get_sequence(self, part_number):
        url = self.BASE_URL +  '/parts/'+ str(part_number) + '/sequence'
        response = requests.get(url, headers = self.HEADERS, verify=False)
        if response.status_code == 200:
            response_json = json.loads (response.text)
            return response_json
        return response.status_code

    # Not working yet...
    def upload_sequence(self, part_number, rid, seqfile):
        HEADERS = {}
        for s in self.HEADERS:
            HEADERS[s] = self.HEADERS
        HEADERS['Accept'] = 'application/json'
        HEADERS['Content-Type'] = 'multipart/form-data'

    #    url = BASE_URL +  '/upload/'+ str(part_number) + '/sequence'
    #    url = BASE_URL +  '/file/sequence'

        url = self.BASE_URL +  '/uploads/'+ str(part_number) + '/sequence'
        response = requests.post(url, files={'file': open(seqfile, 'rb'), 'entryId': str(part_number)},
                             #                                         'entryRecordId': str(rid),
                             #                                         'entryType':'Plasmid'},
                                 headers = self.HEADERS, verify=False)
        if response.status_code == 200:
            response_json = json.loads (response.text)
            return response_json
        return response.status_code
    
