#First we load the data
widedata = read.csv("C:/Users/838035/Desktop/MS.csv", header = T, sep = ";")

#Find those proteins which have 0 spectral count in all patients
spectral_count_sum_by_protein = apply(widedata[,-1], MARGIN = 1, sum)

#There are 881 proteins which 0 spectral counts
sum(spectral_count_sum_by_protein==0)

#We remove these 881 proteins because SPI cannot be calculated for these proteins
widedata = widedata[spectral_count_sum_by_protein>0,]

totalPatients = ncol(widedata) - 1
totalProteins = nrow(widedata)

#Now we convert wide format data to long format for analysis
longdata = reshape(widedata, varying = 2:13, timevar = 1, idvar = "PID", 
        direction="long", sep = "")

colnames(longdata) = c("Patient_ID", "Spectral_Count", "Protein_ID")

#Now we add a label for who is CF (TRUE) and who is not (FALSE)
longdata$CF = longdata$Patient_ID %in% seq(2,12,2)

#Order longdata by patient ID, then protein ID
longdata = longdata[order(longdata$Patient_ID, longdata$Protein_ID),]

#SPI function
spi_calculator = function(protein_data){
  mean_scf = mean(protein_data$Spectral_Count[protein_data$CF==TRUE])
  mean_snl = mean(protein_data$Spectral_Count[protein_data$CF==FALSE])
  
  nt_cf = sum(protein_data$CF == TRUE)
  nt_nl = sum(protein_data$CF == FALSE)
  
  nd_cf = sum(protein_data$Spectral_Count > 0 & protein_data$CF == TRUE)
  nd_nl = sum(protein_data$Spectral_Count > 0 & protein_data$CF == FALSE)
  
  left_part = (mean_scf / (mean_scf + mean_snl)) * (nd_cf/nt_cf)
  right_part = (mean_snl / (mean_scf + mean_snl)) * (nd_nl/nt_nl)
  
  spi = left_part - right_part
  return(spi)
}

spi_all_proteins = by(data = longdata, INDICES = longdata$Protein_ID, FUN = spi_calculator)

####################
#Now we do permutation analysis.
####################

#If we have N patients and K proteins, permutation over CF/NonCF should be done
# at patient level. That is entire patient should be permuted to be a CF/NonCF
#Thus total number of permutations are 2^N
#However 2 permutations where all patients are assigned CF and all patients
# are assigned NonCF are not possible as spi is not defined for them
#So total number of permutations are 2^N - 2

#Now here I generate 2^N - 2 combinations first
combinationGrid = as.matrix(expand.grid(rep(list(c(T,F)), totalPatients)))
combinationGrid = combinationGrid[-c(1, nrow(combinationGrid)),]

permuted_data_set_list = lapply(1:nrow(combinationGrid), function(i){
  newdata = longdata
  newdata$CF = rep(combinationGrid[i,], each=totalProteins)
  return(newdata)
})

#this takes hours to run
spi_all_proteins_all_permuted_datasets = sapply(permuted_data_set_list, function(permuted_data_set){
  by(data = permuted_data_set, INDICES = permuted_data_set$Protein_ID, FUN = spi_calculator)
})

#print the quantiles
round(quantile(as.numeric(spi_all_proteins_all_permuted_datasets), probs = seq(0,1, 0.005)),2)

save.image(file = file.choose())
