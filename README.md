# Tryps_Analysis
Analysis of Mass_spectrometry_data focusing on proteins involved in Antigenic variation in Trypanosomes.
##Analysis of Mass spec data

##import csv data

library(readxl)
library(readr)

library(tidyverse) 
library(ggpubr)
library(ggbeeswarm)
library(ggplot2)


Mass_spec_data1 <-read.csv('D:/MCGILL UNIVERSITY-IGOR PHD DETAILS/CESTARI_MASS_SPEC_DATA/Cestari_Lab_Tony_Final_Cleaned.csv')

view(Mass_spec_data1)

##p values vs number of peptides


Mass_spec_data1 %>% 
  ggplot(aes(x=pvalue, y=Number_of_Peptides))+
  facet_wrap(~Fractions)+
  geom_line() +
  geom_point()

Mass_spec_data1 %>% 
  ggplot(aes(x=pvalue, y=log10(Number_of_Peptides)))+
  facet_wrap(~Fractions)+
  geom_line() +
  geom_point()


##Log10SIn vs Pvalue

Mass_spec_data1 %>% 
  ggplot(aes(x=pvalue, y=Log2_SIn))+
  facet_wrap(~Fractions)+
  geom_point()

##Retrieve_data_reis_cestari_faria

Reis <-read_excel('D:/MCGILL UNIVERSITY-IGOR PHD DETAILS/Reis_Mass_Spectrometry_Data.xlsx')
Cestari_PIP5Pase_RAP1<-read_excel('D:/MCGILL UNIVERSITY-IGOR PHD DETAILS/Mass spec data_Cestari et al. Paper.xlsx',sheet = 'PIP5Pase & RAP1 MS')
Cestari_RAP1_10S_20S<-read_excel('D:/MCGILL UNIVERSITY-IGOR PHD DETAILS/Mass spec data_Cestari et al. Paper.xlsx',sheet = 'RAP1 10 & 20S MS')
Cestari_Proteins_Not_Present_in_Control<-read_excel('D:/MCGILL UNIVERSITY-IGOR PHD DETAILS/Mass spec data_Cestari et al. Paper.xlsx',sheet = 'Proteins not present in control')

Faria_Bloodstream<-read_excel('D:/MCGILL UNIVERSITY-IGOR PHD DETAILS/Cestari_Lab_Mass_Spectrometry_Faria_et_al.xlsx', sheet = 'Bloodstream_Proteomics')
na.omit(Faria_Bloodstream)
Faria_Procyclins<-read_excel('D:/MCGILL UNIVERSITY-IGOR PHD DETAILS/Cestari_Lab_Mass_Spectrometry_Faria_et_al.xlsx', sheet = 'Procyclins_Proteomics')
na.omit(Faria_Procyclins)
Fbs<-na.omit(Faria_Bloodstream)
Fpc<-na.omit(Faria_Procyclins)

View(Cestari_PIP5Pase_RAP1)

#cutout fold change>2

C_PIP5Pase_RAP1<-Cestari_PIP5Pase_RAP1 %>% filter(`Log2_FC`>3)
View(C_PIP5Pase_RAP1)


C_PIP5Pase_RAP1 %>% ggplot(aes(x=Protein_ID, y=Log2_FC)) +
  geom_bar(stat = 'identity', position = position_dodge())+
  theme(axis.text.x=element_text(angle=60,hjust = 1))+
  ggtitle('PIP5Pase&RAP1_Enriched FC>3')

C_PIP5Pase_RAP2<-Cestari_PIP5Pase_RAP1 %>% filter(`Log2_FC`>3&`p-value`<0.05)


C_PIP5Pase_RAP2 %>% ggplot(aes(x=Protein_ID, y=Log2_FC)) +
  geom_bar(stat = 'identity', position = position_dodge())+
  theme(axis.text.x=element_text(angle=60,hjust = 1))+
  ggtitle('PIP5Pase&RAP1_Enriched FC>3, pvalue<0.05')



C_RAP10_20S<-Cestari_RAP1_10S_20S %>% filter(`Log2_FC`>2)

C_PIP5Pase_RAP2<-Cestari_PIP5Pase_RAP1 %>% filter(`Log2_FC`>3&`p-value`<0.05)


C_RAP20 %>% ggplot(aes(x=Protein_ID, y=Log2_FC)) +
  geom_bar(stat = 'identity', position = position_dodge())+
  theme(axis.text.x=element_text(angle=60,hjust = 1))

C_RAP20A<-Cestari_PIP5Pase_RAP1 %>% filter(`Log2_FC`>3&`p-value`<0.05)


C_RAP20<-C_RAP10_20S %>% filter(Fractions>10)

view(C_RAP20)

C_RAP20S<-C_RAP20 %>% filter(`Log2_FC`>2&`p-value`<0.05)
view(C_RAP20S)



C_RAP20AS %>% ggplot(aes(x=Protein_ID, y=Log2_FC)) +
  geom_bar(stat = 'identity', position = position_dodge())+
  theme(axis.text.x=element_text(angle=60,hjust = 1))+
  ggtitle('RAP20S_Enriched FC>3, pvalue<0.05')

C_RAP10<-C_RAP10_20S %>% filter(Fractions<20)

C_RAP10S<-C_RAP10 %>% filter(`Log2_FC`>3& `p-value`<0.05)
view(C_RAP10S)

C_RAP10S %>% ggplot(aes(x=Protein_ID, y=Log2_FC)) +
  geom_bar(stat = 'identity', position = position_dodge())+
  theme(axis.text.x=element_text(angle=60,hjust = 1))+
  ggtitle('RAP10S_Enriched FC>3, pvalue<0.05')


view(C_RAP10)

##REIS

Reis_FC<-Reis %>% filter(`Log2_FC`>3)



Reis_FC %>% ggplot(aes(x=Protein_ID, y=Log2_FC)) +
  geom_bar(stat = 'identity', position = position_dodge())+
  theme(axis.text.x=element_text(angle=60,hjust = 1))+
  ggtitle('REIS_Enriched FC>3')

Reis_FC2<-Reis %>% filter(`Log2_FC`>3& `p-value`<0.05)

Reis_FC2 %>% ggplot(aes(x=Protein_ID, y=Log2_FC)) +
  geom_bar(stat = 'identity', position = position_dodge())+
  theme(axis.text.x=element_text(angle=60,hjust = 1))+
  ggtitle('REIS_Enriched FC>3, pvalue<0.05')


##FARIA
Faria_BS<-Fbs %>% filter(`Log2_FC`>0.5)

Faria_BS %>% ggplot(aes(x=Protein_ID, y=Log2_FC)) +
  geom_bar(stat = 'identity', position = position_dodge())+
  theme(axis.text.x=element_text(angle=60,hjust = 1))+
  ggtitle('Faria_Bloodstream_VEX1_Enriched FC>0.5')

Faria_PC<-Fpc %>% filter(`Log2_FC`>0.5)

Faria_PC %>% ggplot(aes(x=Protein_ID, y=Log2_FC)) +
  geom_bar(stat = 'identity', position = position_dodge())+
  theme(axis.text.x=element_text(angle=60,hjust = 1))+
  ggtitle('Faria_Procyclin_VEX1_Enriched FC>0.5')


##LogSIn vs Log2_FC
#Cestari

C_RAP20 %>% ggplot(aes(x=Log2_FC, y=`Log2 SIN_(Mean)`)) +
  geom_bar(stat = 'identity', position = position_dodge())+
  theme(axis.text.x=element_text(angle=60,hjust = 1))



###Total number of proteins identified in all fractions


Totals <-read_excel('D:/MCGILL UNIVERSITY-IGOR PHD DETAILS/Total proteins.xlsx')

view(Totals)

Totals %>% ggplot(aes(x=Fraction, y=Proteins_Identified)) +
  geom_bar(stat = 'identity', position = position_dodge())+
  theme(axis.text.x=element_text(angle=60,hjust = 1))+
  ggtitle('Total number of proteins identified from each fraction')

Totals2 <-read_excel('D:/MCGILL UNIVERSITY-IGOR PHD DETAILS/Total proteins.xlsx', sheet = 'Totals2')

Totals2 %>% ggplot(aes(x=Fraction, y=Proteins_Identified)) +
  geom_bar(stat = 'identity', position = position_dodge())+
  theme(axis.text.x=element_text(angle=60,hjust = 1))+
  ggtitle('Proteins following quantification')


Totals3 <-read_excel('D:/MCGILL UNIVERSITY-IGOR PHD DETAILS/Total proteins.xlsx', sheet = 'Common_proteins')




