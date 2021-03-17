
#Gisaid class is used to create aggregate that contain gisaid_epi_isl as root and 
#encapsulate/privatize the entity strain
class Gisaid:    
    def __init__(self, gisaid_epi_isl='?'):
        self.gisaid_epi_isl = gisaid_epi_isl
       #sself.genbank_accession = genbank_accession

    def addStrain(self, strainName,date):
        self.strain = strainName
        self.date = date


#######################################################
#Strain class is to store the strainname/id that can be later refernced to match with data in fasta file
# as part of aggregate strain should only be accesed via Gisaid
class Strain: 
    def __init__(self, strain):
        self.strain = strain

#######################################################
#Region class is independent 1 aggregate unit
# Region is aggregate root and once we have region it can have country,division,location
# but these entities dont exists without region 

class Region:
    def __init__(self, region):
        self.region = region
    
    def addData(self, region, country, division ='', location  = ''):
        self.region = region
        self.country = country 
        self.division = division
        self.location = location

 