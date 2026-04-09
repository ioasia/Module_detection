ALL_color_palette  <-  function () {
  library(magrittr)
  colorScheme <- list()
  
  colorScheme$colorPanel1 <- c(
    '1' = '#8dd3c7',
    '2' = '#b3b300',
    '3' = '#bebada',
    '4' = '#fb8072',
    '5' = '#80b1d3',
    '6' = '#fdb462',
    '7' = '#b3de69',
    '8' = '#fccde5'
  )
  
  colorScheme$colorPanel2 <- c(
    '1' = '#75b4c0',
    '2' = '#438644',
    '3' = '#8bc2a1',
    '4' = '#f0c0b5',
    '5' = '#ebe5e6',
    '6' = '#216100',
    '7' = '#8ab22e',
    '8' = '#bc8e46',
    '9' = '#e3d9ad',
    '10' = '#633e5a',
    '11' = '#cb95aa',
    '12' = '#fdae55',
    '13' = '#658cc8',
    '14' = '#007d92',
    '15' = '#62a399',
    '16' = '#fd7e52',
    '17' = '#f3166b',
    '18' = '#ab4252',
    '19' = '#fea617',
    '20' = '#f9c099'
  )
  
  colorScheme$blue <- 'deepskyblue3'
  
  colorScheme$white <- 'white'
  
  colorScheme$red <- 'firebrick3'
  
  colorScheme$samplesClusters$kmeans  <- colorScheme$colorPanel1
  colorScheme$samplesClusters$hierar    <-  colorScheme$colorPanel1
  
  colorScheme$featureClusters$hierar <- colorScheme$colorPanel2
  
  colorScheme$all.clc <- colorScheme$colorPanel2
  colorScheme$b.clc <- colorScheme$colorPanel2
  
  colorScheme$t.clc <- colorScheme$colorPanel2
  
  colorScheme$subtype <- c(
    'ABL1-ZMIZ1' = '#5C80BC',
    'B-Other' = '#9dcbf9',
    'BCL11B-TLX3' = '#71BC78',
    'BCR-ABL1' = '#FF9BAA',
    'CRLF2-IGH' = '#FF80ED',
    'EBV' = '#7c827c',
    'EP300-ZNF384' = '#8b4789',
    'ETV6-PDGFRB' = '#FD5E53',
    'ETV6-RUNX1' = '#FD7C6E',
    'Hypodiploid-LMO2' = '#1DACD6',
    'IGH-CRLF2' = '#1CD3A2',
    'IGH-MYC' = '#1CAC78',
    'KMT2A-AFF1' = '#FCD975',
    'KMT2A-EPS15' = '#7962f7',
    'KMT2A-FOXO4' = '#FCE883',
    'KMT2A-MLLT1' = '#FDFC74',
    'KMT2A-MLLT10' = '#FDFC74',
    'KMT2A-MLLT3' = '#89d14f',
    'KMT2A-SEPT2' = '#9b60db',
    'LCK-TRB' = '#685f4d',
    'LMO1-TCRD' = '#fce1b8',
    'LMO1-TRA' = '#c36d05',
    'LMO2-STAG2' = '#A2ADD0',
    'MEF2D-BCL9' = '#D68A59',
    'MEF2D-DAZAP1' = '#D68A59',
    'MEF2D-HNRNPUL1' = '#D68A59',
    'MYC-IGH' = '#d7cab8',
    'Near-haploid' = '#ECEABE',
    'NUP214-ABL1' = '#9D81BA',
    'P2RY8-CRLF2' = '#1CD3A2',
    'PAX5-ETV6' = '#FFA089',
    'POU2AF1-BCL6' = '#c6a88e',
    'PTMA-TMSB4X' = '#E0AFA0',
    'RAP1GDS1-NUP98' = '#a0c5e8',
    'SET-NUP214' = '#a3bced',
    'SFPQ-ABL' = '#f9f29d',
    'SFPQ-ABL1' = '#B5CA8D',
    'SIL-SCL' = '#EFCDB8',
    'STAG2-PTGER3' = '#A2ADD0',
    'STIL-TAL1' = '#fc9d2a',
    'T-Other' = '#bbbbfa',
    'TCF3-HLF' = '#F664AF',
    'TCF3-PBX1' = '#E6A8D7',
    'TRA-MYC' = '#FC6C85',
    'TRB-LCK' = '#1B264F',
    'TRB-NOTCH1' = '#b73133',
    'DUX4-IGH' = '#325d88',
    'Unknown' = 'gray',
    'NA' = 'gray'
  )
  
  
  
  colorScheme$grouped_subtype <- c(
    'ABL1r' = '#5C80BC',
    'B-Other' = '#9dcbf9',
    'CRLF2' = '#fce1b8',
    'DUX4-IGH' = '#325d88',
    'EBV' = '#7c827c',
    'ETP-ALL' = '#ECEABE',
    'ETV6-RUNX1' = '#FD7C6E',
    'ETV6r'  = '#FFA089',
    'KMT2Ar' = '#89d14f',
    'LMO' = '#A2ADD0',
    'MEF2Dr' = '#D68A59',
    'MYCr' = '#d7cab8',
    'PAX5-ETV6' = '#FFA089',
    'PAX5r' = '#e6ccb3',
    'STAG2' = '#A2ADD0',
    'T-Other' = '#bbbbfa',
    'TAL' = '#f9f29d',
    'TCF3-HLF' = '#F664AF',
    'TCF3-PBX1' = '#E6A8D7',
    'TCF3r' = '#808000',
    'TLX3' = '#E0AFA0',
    'TRB-NOTCH1' = '#b73133',
    'Unknown' = 'gray',
    'NA' = 'gray'
  )
  
  
  
  colorScheme$type <- c(
    'T-ALL' = '#fbb4ae',
    'B-ALL' = '#b3cde3',
    'preB' = '#ccebc5',
    'BCP-ALL' = '#ccebc5',
    'possibly_B' = 'gray',
    'EBV' = '#7c827c',
    'MPAL' = '#FDFC74',
    'LYMPHOMA' = '#984ea3',
    'Unknown' = 'gray',
    'NA' = 'gray'
  )
  
  colorScheme$lineage <- c(
    'T' = '#fbb4ae',
    'B' = '#ccebc5',
    'Unknown' = 'gray',
    'NA' = 'gray'
  )
  
  
  colorScheme$tissue <- c(
    'PB' = '#ff9999',
    'BM' = '#e6ccb3',
    'PE' = '#ff9999',
    'Unknown' = 'gray',
    'NA' = 'gray'
  )
  
  
  colorScheme$batch <- c(
    '1' = '#8b4789',
    '2' = '#63b8ff',
    '3' = '#33A02C',
    '4' = '#f58e2d',
    '5' = '#bbbb77',
    '6' = '#F664AF',
    '7' = 'tomato3',
    '8' = 'thistle4',
    'Unknown' = 'gray',
    'NA' = 'gray'
  )
  
  colorScheme$gender <- c('F' = '#ff99c2',
                          'M' = '#99ccff',
                          'Unknown' = 'gray',
                          'NA' = 'gray'
  )
  
  colorScheme$age <- c(
    'min' = '#99ebff',
    'max' = '#bbbb77'
  )
  
  
  
  
  colorScheme$subtype <- c(colorScheme$subtype,
                           colorScheme$subtype %>% setNames(gsub('-',
                                                                 '::',
                                                                 names(.))),
                           colorScheme$subtype %>% setNames(gsub('-',
                                                                 '.',
                                                                 names(.))))
  
  
  colorScheme$grouped_subtype <- c(colorScheme$grouped_subtype,
                                   colorScheme$grouped_subtype %>% setNames(gsub('-',
                                                                                 '::',
                                                                                 names(.))),
                                   colorScheme$grouped_subtype %>% setNames(gsub('-',
                                                                                 '.',
                                                                                 names(.))))
  
  
  colorScheme$Lineage <-  colorScheme$lineage
  colorScheme$Type <- colorScheme$type
  colorScheme$Subtype <- colorScheme$subtype
  colorScheme$Grouped_Subtype <- colorScheme$grouped_subtype
  colorScheme$Tissue <- colorScheme$tissue
  colorScheme$Gender <- colorScheme$gender
  colorScheme$sex <- colorScheme$gender
  colorScheme$Sex <- colorScheme$gender
  
  
  return(colorScheme)
}
