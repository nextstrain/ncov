#!/bin/bash
#
# Call as: scripts/normalize_lab_metadata.sh data/metadata.tsv
# Creates: data/metadata.tsv.bak
#
# Normalize origin/submitter lab metadata by adjusting spacing,
# punctuation, and spelling variants.

METADATA_IN=$1

if [[ ! -r "$METADATA_IN" ]]
then
	echo "$0: input $METADATA_IN not found or not readable"
	exit 1
fi

METADATA_BAK=${METADATA_IN}.bak
METADATA_OUT=${METADATA_IN}

mv $METADATA_IN $METADATA_BAK

sed '
# Australia
##
# Centre for Infectious Diseases and Microbiology - Public Health (4)
# Centre for Infectious Diseases and Microbiology Public Health (30)
# Centre for Infectious Diseases and Microbiology- Public Health (1)
##
/^Australia/s/Microbiology- Public Health/Microbiology Public Health/;
/^Australia/s/Microbiology - Public Health/Microbiology Public Health/;


# Belgium
##
# KU Leuven, Clincal and Epidemiological Virology (1)
# KU Leuven, Clinical and Epidemiological Virology (42)
# Institute information  KU Leuven, Clinical and Epidemiological Virology (1)
##
#### Both origin and submitter, use /g
/^Belgium/s/KU Leuven, Clincal/KU Leuven, Clinical/g;
/^Belgium/s/Institute information  KU Leuven, Clinical and Epidemiological Virology/KU Leuven, Clinical and Epidemiological Virology/g;


# Brazil
## 
# Instituto Adolfo Lutz Interdisciplinary Procedures Center Strategic Laboratory (1)
# Instituto Adolfo Lutz, Interdiciplinary Procedures Center, Strategic Laboratory (12)
# Instituto Adolfo Lutz, Interdisciplinary Procedures Center, Strategic Laboratory (1)
# Instituto Adolfo Lutz, Interdisciplinary Procedures Center Strategic Laboratory (1)
##
/^Brazil/s/Adolfo Lutz Interdisciplinary/Adolfo Lutz, Interdisciplinary/;
/^Brazil/s/Adolfo Lutz, Interdisciplinary Procedures Center Strategic/Adolfo Lutz, Interdiciplinary Procedures Center, Strategic/;
/^Brazil/s/Interdiciplinary/Interdisciplinary/;
##
# Bioinformatics Laboratory - LNCC (9)
# Bioinformatics Laboratory / LNCC (2)
##
/^Brazil/s/Bioinformatics Laboratory \/ LNCC/Bioinformatics Laboratory - LNCC/;
##
# Universidade Federal do Rio de Janeiro (2)
# Universidade Federal do Rio de Janeiro - UFRJ (1)
##
/^Brazil/s/Universidade Federal do Rio de Janeiro - UFRJ/Universidade Federal do Rio de Janeiro/;


# Canada
## 
# Public Health Ontario Laboratories (55)
# Public Health Ontario Laboratory (2)
##
/^Canada/s/Public Health Ontario Laboratory/Public Health Ontario Laboratories/g;


# China
##
# National Institute for Viral Disease Control & Prevention, CCDC (6)
# National Institute for Viral Disease Control and Prevention, China CDC (3)
##
#### Not limited to top-level match because multiple: Sichuan, Wuhan, etc
s/National Institute for Viral Disease Control \& Prevention, CCDC/National Institute for Viral Disease Control and Prevention, China CDC/;
# Guangdong
##
# Guangdong Provincial Center for Diseases Control and Prevention; Guangdong Provincial Public Health (9)
# Guangdong Provincial Center for Diseases Control and Prevention; Guangdong Provinical Public Health (2)
# Guangdong Provincial Center for Diseases Control and Prevention;Guangdong Provincial Institute of Public Health (1)
# Guangdong Provincial Institution of Public Health, Guangdong Provinical Center for Disease Control and Prevention (41)
#submitter
# Guangdong Provincial Center for Disease Control and Prevention (1)
# Guangdong Provincial Center for Diseases Control and Prevention (5)
##
# Prefer:
# Guangdong Provincial Institution of Public Health, Guangdong Provincial Center for Disease Control and Prevention
#### Both origin and submitter, use /g, also multiple toplevel names so global
s/\tGuangdong Provincial Center for Disease Control and Prevention\t/\tGuangdong Provincial Institution of Public Health, Guangdong Provincial Center for Disease Control and Prevention\t/g;
s/\tGuangdong Provincial Center for Diseases Control and Prevention\t/\tGuangdong Provincial Institution of Public Health, Guangdong Provincial Center for Disease Control and Prevention\t/g;
s/\tGuangdong Provincial Institution of Public Health, Guangdong Provinical Center for Disease Control and Prevention\t/\tGuangdong Provincial Institution of Public Health, Guangdong Provincial Center for Disease Control and Prevention\t/g;
s/\tGuangdong Provincial Center for Diseases Control and Prevention;Guangdong Provincial Institute of Public Health\t/\tGuangdong Provincial Institution of Public Health, Guangdong Provincial Center for Disease Control and Prevention\t/g;
s/\tGuangdong Provincial Center for Diseases Control and Prevention; Guangdong Provinical Public Health\t/\tGuangdong Provincial Institution of Public Health, Guangdong Provincial Center for Disease Control and Prevention\t/g;
s/\tGuangdong Provincial Center for Diseases Control and Prevention; Guangdong Provincial Public Health\t/\tGuangdong Provincial Institution of Public Health, Guangdong Provincial Center for Disease Control and Prevention\t/g;
s/\tGuangdong Provincial Institution of Public Health\t/\tGuangdong Provincial Institution of Public Health, Guangdong Provincial Center for Disease Control and Prevention\t/g;


# Colombia
##
# Instituto Nacional de Salud Universidad Cooperativa de Colombia Instituto Alexander von Humboldt Imperial College-London London School of Hygiene & Tropical Medicine (1)
# Instituto Nacional de Salud, Universidad Cooperativa de Colombia, Instituto Alexander von Humboldt, Imperial College-London, London School of Hygiene & Tropical Medicine (1)
##
/^Colombia/s/Instituto Nacional de Salud Universidad Cooperativa de Colombia Instituto Alexander von Humboldt Imperial College-London London School of Hygiene \& Tropical Medicine/Instituto Nacional de Salud, Universidad Cooperativa de Colombia, Instituto Alexander von Humboldt, Imperial College-London, London School of Hygiene \& Tropical Medicine/;


# England
##
# Department of Infection, Immunity and Cardiovascular Disease, The Florey Institute,  The Medical School, University of Sheffield (1)
# Department of Infection, Immunity and Cardiovascular Disease, The Florey Institute, The Medical School, University of Sheffield (133)
##
/^England/s/The Florey Institute,  The Medical School/The Florey Institute, The Medical School/;


# France
##
# Institut des Agents Infectieux (IAI) Hospices Civils de Lyon (7)
# Institut des Agents Infectieux (IAI), Hospices Civils de Lyon (28)
##
/^France/s/Infectieux (IAI) Hospices/Infectieux (IAI), Hospices/;


# Hangzhou
##
# State Key Laboratory for Diagnosis and Treatment of Infectious Diseases, National Clinical Research Center for Infectious Diseases, First Affiliated Hospital, Zhejiang University School of Medicine, Hangzhou, China 310003 (8)
# State Key Laboratory for Diagnosis and Treatment of Infectious Diseases, National Clinical Research Center for Infectious Diseases, First Affiliated Hospital, Zhejiang University School of Medicine, Hangzhou, China. 310003 (3)
##
#### Both origin and submitter, use /g
/^Hangzhou/s/China\./China/g;
##
# Insepction Center of Hangzhou Center for Disease Control and Prevention (1)
# Inspection Center of Hanghzou Center for Disease Control and Prevention (1)
# Inspection Center of Hangzhou Center for Disease Control and Prevention (13)
##
/^Hangzhou/s/Insepction/Inspection/;
/^Hangzhou/s/Hanghzou/Hangzhou/;


# Hong Kong
##
# School of Public Health, The University of Hon g Kong (2)
# School of Public Health, The University of Hong Kong (7)
##
/^HongKong/s/University of Hon g Kong/University of Hong Kong/;
##
# Department of Clinical Pathology, Tuen Mun Hospital (1)
# Department of Clinical Pathology, Tuen Mun Hospital, 23 Tsing Chung Koon Road, Tuen Mun, N.T. (1)
##
/^HongKong/s/Department of Clinical Pathology, Tuen Mun Hospital, 23 Tsing Chung Koon Road, Tuen Mun, N\.T\./Department of Clinical Pathology, Tuen Mun Hospital/;
##
# Chinese University of Hong Kong, Hong Kong SAR, China (1)
# The Chinese University of Hong Kong, Hong Kong SAR, China (1)
##
/^HongKong/s/The Chinese University of Hong Kong/Chinese University of Hong Kong/;


# Hungary
##
# Bioinformatics Research Group, Szentagothai Research Centre (2)
# Bioinformatics Research Group, Szentagothai Research Centre, University of Pecs (1)
#
##
/^Hungary/s/Bioinformatics Research Group, Szentagothai Research Centre\t/Bioinformatics Research Group, Szentagothai Research Centre, University of Pecs\t/;
##
# Virological Research Group, Szentagothai Research Centre (2)
# Virological Research Group, Szentagothai Research Centre, University of Pecs (1)
##
/^Hungary/s/Virological Research Group, Szentagothai Research Centre\t/Virological Research Group, Szentagothai Research Centre, University of Pecs\t/;


# India
##
# Indian Council of Medical Research - National Institute of Virology (1)
# Indian Council of Medical Research-National Institute of Virology (1)
##
#### Both origin and submitter, use /g
/^India/s/Indian Council of Medical Research-National Institute of Virology/Indian Council of Medical Research - National Institute of Virology/g;


# Italy
##
# Department of Infectious Diseases, Istituto Superiore di Sanita, Roma , Italy (1)
# Department of Infectious Diseases, Istituto Superiore di Sanita, Rome, Italy (1)
##
/^Italy/s/Superiore di Sanita, Roma , Italy/Superiore di Sanita, Rome, Italy/;
##
# Laboratory of Molecular Virology International Center for Genetic Engineering and Biotechnology (ICGEB) (3)
# Laboratory of Molecular Virology International Center fro Genetic Engineering and Biotechnology (ICGEB) (1)
##
/^Italy/s/ fro Genetic Engineering/ for Genetic Engineering/;
## also sequences with many embedded double-quotes
#### just remove all the double-quotes
#### both origin/submitter so use /g
/^Italy/s/\"//g;
/^Italy/s/G\.Caporale/G\. Caporale/;
/^Italy/s/Abruzzo e Molise/Abruzzo e del Molise/;


# Japan
##
# Department of Virology III, National Institute of Infectious Diseases (1)
# Dept. of Virology III, National Institute of Infectious Diseases (4)
##
/^Japan/s/Dept\. of Virology III, National Institute of Infectious Diseases/Department of Virology III, National Institute of Infectious Diseases/;


# Luxembourg
##
# Laboratoire National de Sante, Microbiology, Virology (46)
# Laboratoire Nationale de Sante, Microbiology, Virology (9)
##
#### Both origin and submitter, use /g
/^Luxembourg/s/Laboratoire Nationale de Sante/Laboratoire National de Sante/g;


# Senegal
##
# Instirut Pasteur Dakar ()
# Institut Pasteur de Dakar (21)
# Institut pasteur Dakar (1)
##
#### Both origin and submitter, use /g
/^Senegal/s/Insti[rt]ut [pP]asteur Dakar/Institut Pasteur de Dakar/g;


# Spain
##
# Fundacion Jimenez Diaz (3)
# FUNDACION JIMENEZ DIAZ (2)
##
/^Spain/s/FUNDACION JIMENEZ DIAZ/Fundacion Jimenez Diaz/;
##
# HOSPITAL UNIVERSITARIO LA PAZ (1)
# Hospital Universitario La Paz (4)
##
/^Spain/s/HOSPITAL UNIVERSITARIO LA PAZ/Hospital Universitario La Paz/;
##
# HOSPITAL UNIVERSITARIO VIRGEN DE LAS NIEVES (2)
# Hospital Universitario Virgen de las Nieves (1)
##
/^Spain/s/HOSPITAL UNIVERSITARIO VIRGEN DE LAS NIEVES/Hospital Universitario Virgen de las Nieves/;
##
# Sequencing and Bioinformatics Service and Molecular Epidemiology Research Group. FISABIO-Public Health (44)
# Sequencing and Bioinformatics Service and Molecular Epidemiology Research Group. FISABIO-Public Health. (1)
##
/^Spain/s/FISABIO-Public Health\./FISABIO-Public Health/;
##
# Servicio Microbiologia, Hospital Clinico Universitario, Valencia (1)
# Servicio Microbiologia. Hospital Clinico Universitario. Valencia. (2)
# Servicio de Microbiologia. Hospital Clinico Universitario de Valencia (8)
##
/^Spain/s/Servicio Microbiologia\. Hospital Clinico Universitario\. Valencia\./Servicio de Microbiologia, Hospital Clinico Universitario de Valencia/;
/^Spain/s/Servicio Microbiologia, Hospital Clinico Universitario, Valencia/Servicio de Microbiologia, Hospital Clinico Universitario de Valencia/;
/^Spain/s/Servicio de Microbiologia. Hospital Clinico Universitario de Valencia/Servicio de Microbiologia, Hospital Clinico Universitario de Valencia/;


# Taiwan
##
# Department of Laboratory Medicine, Lin-Kou Chang Gung Memorial Hospital, Taoyuan, Taiwan (11)
# Department of Laboratory Medicine, Lin-Kou Chang Gung Memorial Hospital, Taoyuan, Taiwan. (1)
##
/^Taiwan/s/Taiwan\./Taiwan/;


# USA/TX
##
# Texas DSHS Lab Services (3)
# Texas Department of State Health Services (1)
# Texas Department of State Health Services Lab Services (1)
##
/^USA\/TX/s/Texas Department of State Health Services Lab Services/Texas DSHS Lab Services/;
/^USA\/TX/s/Texas Department of State Health Services/Texas DSHS Lab Services/;


# USA/VA
##
# Division of Consolidated Laboratories (5)
# Division of Consolidated Laboratories Services (1)
# Division of Consolidated Laboratory Services (19)
# VA DCLS (1)
# Virginia Division of Consolidated Laboratories (3)
#### Both origin and submitter, use /g
##
/^USA\/VA/s/VA DCLS\t/Virginia Division of Consolidated Laboratory Services\t/g;
/^USA\/VA/s/Virginia Division of Consolidated Laboratories/Virginia Division of Consolidated Laboratory Services/g;
/^USA\/VA/s/\tDivision of Consolidated Laboratories Services/\tVirginia Division of Consolidated Laboratory Services/g;
/^USA\/VA/s/\tDivision of Consolidated Laboratories/\tVirginia Division of Consolidated Laboratory Services/g;
/^USA\/VA/s/\tDivision of Consolidated Laboratory Services/\tVirginia Division of Consolidated Laboratory Services/g;


# USA/WI
##
# Gundersen Molecular Diagnostic Laboratory (2)
# Gundersen Molecular Diagnostics Laboratory (5)
##
/^USA\/WI/s/Gundersen Molecular Diagnostic Laboratory/Gundersen Molecular Diagnostics Laboratory/;
##
# University of Wisconsin - Madison AIDS Vaccine Research Laboratories (1)
# University of Wisconsin-Madison AIDS Vaccine Research Laboratories (15)
# University of Wisconsin-Madison AIDS Vaccine Research Laboratory (1)
# University of Wisconsin-Madison, AIDS Vaccine Research Laboratories (1)
# University of Wisconsin Madison, AIDS Vaccine Research Laboratories (1)
# AIDS Vaccine Research Laboratories (4)
#### Both origin and submitter, use /g
#### Also a Japan seq
##
/^USA\/WI/s/\tAIDS Vaccine Research Laboratories/\tUniversity of Wisconsin-Madison AIDS Vaccine Research Laboratories/g;
/^USA\/WI/s/University of Wisconsin - Madison AIDS Vaccine Research Laboratories/University of Wisconsin-Madison AIDS Vaccine Research Laboratories/g;
/^USA\/WI/s/Madison AIDS Vaccine Research Laboratory/Madison AIDS Vaccine Research Laboratories/g;
/^USA\/WI/s/Wisconsin-Madison, AIDS Vaccine Research Laboratories/Wisconsin-Madison AIDS Vaccine Research Laboratories/g;
/^USA\/WI/s/Wisconsin Madison, AIDS Vaccine Research Laboratories/Wisconsin-Madison AIDS Vaccine Research Laboratories/g;
/^Japan/s/University of Wisconsin Madison, AIDS Vaccine Research Laboratories/University of Wisconsin-Madison AIDS Vaccine Research Laboratories/;


# Vietnam
##
# National Influenza Center - National Institute of Hygiene and Epidemiology (NIHE) (1)
# National Influenza Center, National Institute of Hygiene and Epidemiology (NIHE) (5)
##
#### Both origin and submitter, use /g
/^Vietnam/s/National Influenza Center, National Institute of Hygiene/National Influenza Center - National Institute of Hygiene/g;

' $METADATA_BAK |

# extra sed command used to deal with single quote in Italy lab metadata
sed "
/^Italy/s/dell'Abruzzo/dellAbruzzo/g;
" - > $METADATA_OUT

exit 0
