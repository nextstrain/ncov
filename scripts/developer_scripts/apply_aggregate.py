
#Gisaid class is used to create aggregate that contain gisaid_epi_isl as root and 
#encapsulate/privatize the entity strain
class Gisaid:    
    def __init__(self, gisaid_epi_isl='?'):
        self.gisaid_epi_isl = gisaid_epi_isl
       #sself.genbank_accession = genbank_accession

    def addStrain(self, strainName,date):
        self.strain = strainName
        self.date = date

