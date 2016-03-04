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

    # This version instatiate all constructs together (see Example B.1.3, BBF RFC 108)
    # Each construct is a composite ComponentDefinition using SequenceConstraint
    # In principle, each subcomponent is defined using ComponentDefinition
    # We should perhaps encapsulate all of them within a Collection with the design id
    def construct2(self, cid, construct):
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


    # This version instatiate all constructs together.
    # It differs from Example B.1.3, BBF RFC 108, but seems to be the correct one
    # Each construct is a composite ComponentDefinition using SequenceConstraint
    # In principle, each subcomponent is defined using ComponentDefinition
    # We should perhaps encapsulate all of them within a Collection with the design id
    def construct3(self, cid, construct):
        construct = filter(lambda x: len(x) > 0, construct)
        construct = map(lambda x: 'http://ice.synbiochem.co.uk/'+x, construct)
        ComponentDefinition = ET.SubElement(self.rdf, 'sbol:ComponentDefinition',
                                            attrib={'rdf:about': cid})
        Type = ET.SubElement(ComponentDefinition, 'sbol:type',
                             attrib={'rdf:resource': "http://www.biopax.org/release/biopax-level3.owl#DnaRegion"})
        for i in range(1, len(construct)):
            SequenceConstraint = ET.SubElement(ComponentDefinition, 'sbol:sequenceConstraint')
            s = ET.SubElement(SequenceConstraint, 'sbol:SequenceConstraint',
                              attrib={'rdf:about': cid+'/r'+str(i)})
            sr = ET.SubElement(s, 'sbol:restriction',
                               attrib={'rdf:resource': 'http://sbols.org/v2#precedes'})
            sr = ET.SubElement(s, 'sbol:subject',
                               attrib={'rdf:resource': construct[i-1]})
            sr = ET.SubElement(s, 'sbol:object',
                               attrib={'rdf:resource':  construct[i]})
        for c in construct:
            component = ET.SubElement(ComponentDefinition, 'sbol:component')
            Comp = ET.SubElement(component, 'sbol:Component', attrib={'rdf:about': c})
            ET.SubElement(Comp, 'sbol:access', attrib={'rdf:resource': "http://sbols.org/v2#public"})
            ET.SubElement(Comp, 'sbol:definition', attrib={'rdf:resource': c})




    def serialize(self, ofile):
        import xml.dom.minidom as minidom
        s = ET.tostring(self.rdf, 'utf-8')
        reparsed = minidom.parseString(s)
        ns = reparsed.toprettyxml(indent='\t')
        of = open(ofile, 'w')
        of.write(s)
        of.close()
        
