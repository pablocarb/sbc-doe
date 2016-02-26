# Some simple routines to create sbol files

import xml.etree.ElementTree as ET

# doe = [ component,...]

class sbol:

    def __init__(self):
        self.rdf = ET.Element('rdf:RDF')
        self.rdf.attrib['xmlns:rdf'] = "http://www.w3.org/1999/02/22-rdf-syntax-ns#"
        self.rdf.attrib['xmlns:dcterms'] = "http://purl.org/dc/terms/"
        self.rdf.attrib['xmlns:prov'] = "http://www.w3.org/ns/prov#"
        self.rdf.attrib['xmlns:sbol'] ="http://sbols.org/v2#"

    def collection(self, doeid, dlist):
        Collection = ET.SubElement(self.rdf, 'sbol:Collection')
        Collection.attrib['rdf:about:'] = doeid
        for m in dlist:
            x = ET.SubElement(Collection, 'sbol:member',
                              attrib={'rdf:resource':  m})

    def construct(self, cid, construct):
        construct = filter(lambda x: len(x) > 0, construct)
        construct = map(lambda x: 'http://ice.synbiochem.co.uk/'+x, construct)
        ComponentDefinition = ET.SubElement(self.rdf, 'sbol:ComponentDefinition',
                                            attrib={'rdf:about': cid})
        Type = ET.SubElement(ComponentDefinition, 'sbol:type',
                             attrib={'rdf:resource': "http://www.biopax.org/release/biopax-level3.owl#DnaRegion"})
        Component = ET.SubElement(ComponentDefinition, 'sbol:component')
        for c in construct:
            x = ET.SubElement(Component, 'sbol:Component', attrib={'rdf:about': c})
            y = ET.SubElement(x, 'sbol:access', attrib={'rdf:resource': "http://sbols.org/v2#public"})
            z = ET.SubElement(x, 'sbol:definition', attrib={'rdf:resource': c})
        SequenceConstraint = ET.SubElement(ComponentDefinition, 'sbol:sequenceConstraint')
        for i in range(1, len(construct)):
            s = ET.SubElement(SequenceConstraint, 'sbol:SequenceConstraint',
                              attrib={'rdf:about': cid+'/r'+str(i)})
            sr = ET.SubElement(s, 'sbol:restriction',
                               attrib={'rdf:resource': 'http://sbols.org/v2#precedes'})
            sr = ET.SubElement(s, 'sbol:subject',
                               attrib={'rdf:resource': construct[i-1]})
            sr = ET.SubElement(s, 'sbol:object',
                               attrib={'rdf:resource':  construct[i]})

    # Test to do all constructs together
    def construct2(self, cid, construct, defcomp=False):
        construct = filter(lambda x: len(x) > 0, construct)
        construct = map(lambda x: 'http://ice.synbiochem.co.uk/'+x, construct)
        ComponentDefinition = ET.SubElement(self.rdf, 'sbol:ComponentDefinition',
                                            attrib={'rdf:about': cid})
        Type = ET.SubElement(ComponentDefinition, 'sbol:type',
                             attrib={'rdf:resource': "http://www.biopax.org/release/biopax-level3.owl#DnaRegion"})
        SequenceConstraint = ET.SubElement(ComponentDefinition, 'sbol:sequenceConstraint')
        for i in range(1, len(construct)):
            s = ET.SubElement(SequenceConstraint, 'sbol:SequenceConstraint',
                              attrib={'rdf:about': cid+'/r'+str(i)})
            sr = ET.SubElement(s, 'sbol:restriction',
                               attrib={'rdf:resource': 'http://sbols.org/v2#precedes'})
            sr = ET.SubElement(s, 'sbol:subject',
                               attrib={'rdf:resource': construct[i-1]})
            sr = ET.SubElement(s, 'sbol:object',
                               attrib={'rdf:resource':  construct[i]})
        for c in construct:
            ComponentDefinition = ET.SubElement(self.rdf, 'sbol:ComponentDefinition',
                                                attrib={'rdf:about': c})
            ET.SubElement(ComponentDefinition, 'sbol:persistentIdentity',
                          attrib={'rdf:resource': c})
            Type = ET.SubElement(ComponentDefinition, 'sbol:type',
                                 attrib={'rdf:resource': "http://www.biopax.org/release/biopax-level3.owl#DnaRegion"})



    def serialize(self, ofile):
        s = '<?xml version="1.0" ?>\n'
        s += ET.tostring(self.rdf)
        of = open(ofile, 'w')
        of.write(s)
        of.close()
        
