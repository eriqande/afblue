library(dplyr)
library(afblue)


#### New Monogamy:  separate sexes with monogamy; true replication of P; skip Fst and Ne ####
## This version adds Yank2 strategy
NLoci = 100
Ne = 500
S = 40
NGens = 10
NReps = 200
MaxSib = seq(1:NReps)
MaxFamily = 30
Familysize = seq(1:MaxFamily)
ProbFamily = 0.1
Sibcheck = 0
Mating = 3  ## 1 = random; 2 = monogamy; 3 = mixed
     if (Mating==3) {Familysize = seq(2:MaxFamily)+1}
SSP = rep(0,NLoci)
PBar = rep(0,NLoci)
RPP = rep(0,NLoci)
PBarYank = rep(0,NLoci)
SSPYank = rep(0,NLoci)
RPPYank = rep(0,NLoci)
PBar1 = rep(0,NLoci)
SSP1 = rep(0,NLoci)
RPP1 = rep(0,NLoci)
PBar2 = rep(0,NLoci)
SSP2 = rep(0,NLoci)
RPP2 = rep(0,NLoci)
PBar3 = rep(0,NLoci)
SSP3 = rep(0,NLoci)
RPP3 = rep(0,NLoci)
PBar4 = rep(0,NLoci)
SSP4 = rep(0,NLoci)
RPP4 = rep(0,NLoci)
PBar5 = rep(0,NLoci)
SSP5 = rep(0,NLoci)
RPP5 = rep(0,NLoci)
PBar6 = rep(0,NLoci)
SSP6 = rep(0,NLoci)
RPP6 = rep(0,NLoci)
SSYank = seq(1:NReps)
SS1 = seq(1:NReps)
SS2 = seq(1:NReps)
SS3 = seq(1:NReps)
SS4 = seq(1:NReps)
SS5 = seq(1:NReps)
SS6 = seq(1:NReps)
SS0 = rep(S,NReps)

PBarBLUE = rep(0,NLoci)
SSP_BLUE = rep(0,NLoci)
RPP_BLUE = rep(0,NLoci)
SS_BLUE = rep(S,NReps)

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}



#### functions for three mating models for first 10 generations #####

Monogamy <- function(parents)           { ## This is separate sexes with monogamy
  # N = number of inds (rows)
  N = dim(parents)[[1]]
  Half = N/2
  Males = seq(1:Half)
  Females = Males + Half
  nextgen = matrix(0,N,NLoci)
  for (i in 1:N)  {
     A = sample(Males,1) # pick a male parent
     B = A+Half # pick the female parent paired with that male
     for (j in 1:NLoci)  {  ## get offspring genotypes for each locus
        if (parents[A,j]==1) { ## heterozygote, so pick an allele
           nextgen[i,j] = nextgen[i,j] + rbinom(1,1,0.5)  } else {  ## homozygote, so take half of the genes
           nextgen[i,j] = nextgen[i,j] + 0.5*parents[A,j]  }  # end first parent
         if (parents[B,j]==1) { ## heterozygote, so pick an allele
           nextgen[i,j] = nextgen[i,j] + rbinom(1,1,0.5)  } else {  ## homozygote, so take half of the genes
           nextgen[i,j] = nextgen[i,j] + 0.5*parents[B,j]  }  # end first parent
           }  ## end for j
           }  ## end for i
  return(nextgen)               } # end function

Separate <- function(parents)           { ## This is separate sexes with random mating
  # N = number of inds (rows)
  N = dim(parents)[[1]]
  Half = N/2
  Males = seq(1:Half)
  Females = Males + Half
  nextgen = matrix(0,N,NLoci)
  for (i in 1:N)  {
     A = sample(Males,1) # pick a male parent
     B = sample(Females,1) # pick a female parent
     for (j in 1:NLoci)  {  ## get offspring genotypes for each locus
        if (parents[A,j]==1) { ## heterozygote, so pick an allele
           nextgen[i,j] = nextgen[i,j] + rbinom(1,1,0.5)  } else {  ## homozygote, so take half of the genes
           nextgen[i,j] = nextgen[i,j] + 0.5*parents[A,j]  }  # end first parent
         if (parents[B,j]==1) { ## heterozygote, so pick an allele
           nextgen[i,j] = nextgen[i,j] + rbinom(1,1,0.5)  } else {  ## homozygote, so take half of the genes
           nextgen[i,j] = nextgen[i,j] + 0.5*parents[B,j]  }  # end first parent
           }  ## end for j
           }  ## end for i
  return(nextgen)               } # end function

Mixed <- function(parents)           { ## This is separate sexes with mixed mating
    # N = number of inds (rows)
    N = dim(parents)[[1]]
    Half = N/2
    Males = seq(1:Half)
    Females = Males + Half
    Offspring = matrix(0,2*N,NLoci)

    OffNumber = 0
    while (OffNumber < N) {     ## produce at least N offspring
       A = sample(Males,1) # pick a male parent
       B = sample(Females,1) # pick a female parent
       C = rbinom(1,1,ProbFamily)
       if (C == 0) { Noff = 1  }  else {Noff = sample(Familysize,1)}
         for (i in 1:Noff)  {
         OffNumber = OffNumber + 1
       for (j in 1:NLoci)  {  ## get offspring genotypes for each locus
          if (parents[A,j]==1) { ## heterozygote, so pick an allele
             Offspring[OffNumber,j] = Offspring[OffNumber,j] + rbinom(1,1,0.5)  } else {  ## homozygote, so take half of the genes
             Offspring[OffNumber,j] = Offspring[OffNumber,j] + 0.5*parents[A,j]  }  # end first parent
           if (parents[B,j]==1) { ## heterozygote, so pick an allele
             Offspring[OffNumber,j] = Offspring[OffNumber,j] + rbinom(1,1,0.5)  } else {  ## homozygote, so take half of the genes
             Offspring[OffNumber,j] = Offspring[OffNumber,j] + 0.5*parents[B,j]  }  # end first parent
             }  ## end for j
             }  ## end for i
             } ## end while
      nextgen = Offspring[1:N,]  ## truncate offspring to desired N
 return(nextgen)            } # end function



#### functions for production of progeny, with random (MaxFamily = 1) or family-correlated samples (MaxFamily > 1) ####

  MonoFamily <- function(parents)           {  ## simulates family-correlated sampling with monogamy
    # N = number of inds (rows)
    N = dim(parents)[[1]]
    Half = N/2
    Males = seq(1:Half)
    Females = Males + Half
    S1 = matrix(0,2*S,2)
    Offspring = matrix(0,2*S,NLoci)

    OffNumber = 0
    while (OffNumber < S) {     ## produce at least S offspring
       A = sample(Males,1) # pick a male parent
       B = A+Half # pick the female parent paired with that male
       Noff = sample(Familysize,1)  # pick number of offspring in this family
         for (i in 1:Noff)  {
         OffNumber = OffNumber + 1
         S1[OffNumber,1] = A        ## record parents of each progeny
         S1[OffNumber,2] = B
       for (j in 1:NLoci)  {  ## get offspring genotypes for each locus
          if (parents[A,j]==1) { ## heterozygote, so pick an allele
             Offspring[OffNumber,j] = Offspring[OffNumber,j] + rbinom(1,1,0.5)  } else {  ## homozygote, so take half of the genes
             Offspring[OffNumber,j] = Offspring[OffNumber,j] + 0.5*parents[A,j]  }  # end first parent
           if (parents[B,j]==1) { ## heterozygote, so pick an allele
             Offspring[OffNumber,j] = Offspring[OffNumber,j] + rbinom(1,1,0.5)  } else {  ## homozygote, so take half of the genes
             Offspring[OffNumber,j] = Offspring[OffNumber,j] + 0.5*parents[B,j]  }  # end first parent
             }  ## end for j
             }  ## end for i
             } ## end while
    OffspringData = cbind(S1,Offspring)
    OffspringData = OffspringData[1:S,]  ## truncate offspring to desired sample size

  return(OffspringData)            } # end function

   RanFamily <- function(parents)           {  ## simulates family-correlated sampling with random mating
      # N = number of inds (rows)
      N = dim(parents)[[1]]
      Half = N/2
      Males = seq(1:Half)
      Females = Males + Half
      S1 = matrix(0,2*S,2)
      Offspring = matrix(0,2*S,NLoci)

      OffNumber = 0
      while (OffNumber < S) {     ## produce at least S offspring
         B = sample(Females,1) # pick a female parent
         Noff = sample(Familysize,1)  # pick number of maternal half sibs for this family
           for (i in 1:Noff)  {
           OffNumber = OffNumber + 1
           A = sample(Males,1)        # pick the male parent
           S1[OffNumber,1] = A        ## record parents of each progeny
           S1[OffNumber,2] = B
         for (j in 1:NLoci)  {  ## get offspring genotypes for each locus
            if (parents[A,j]==1) { ## heterozygote, so pick an allele
               Offspring[OffNumber,j] = Offspring[OffNumber,j] + rbinom(1,1,0.5)  } else {  ## homozygote, so take half of the genes
               Offspring[OffNumber,j] = Offspring[OffNumber,j] + 0.5*parents[A,j]  }  # end first parent
             if (parents[B,j]==1) { ## heterozygote, so pick an allele
               Offspring[OffNumber,j] = Offspring[OffNumber,j] + rbinom(1,1,0.5)  } else {  ## homozygote, so take half of the genes
               Offspring[OffNumber,j] = Offspring[OffNumber,j] + 0.5*parents[B,j]  }  # end first parent
               }  ## end for j
               }  ## end for i
               } ## end while
      OffspringData = cbind(S1,Offspring)
      OffspringData = OffspringData[1:S,]  ## truncate offspring to desired sample size

  return(OffspringData)            } # end function

   MixedFamily <- function(parents)           {  ## simulates same mixed mating model but allows samples larger than NS >= N
          # N = number of inds (rows)
          N = dim(parents)[[1]]
          Half = N/2
          Males = seq(1:Half)
          Females = Males + Half
          S1 = matrix(0,2*S,2)
          Offspring = matrix(0,2*S,NLoci)

          OffNumber = 0
          while (OffNumber < S) {     ## produce at least S offspring
             A = sample(Males,1) # pick a male parent
             B = sample(Females,1) # pick a female parent
             C = rbinom(1,1,ProbFamily)
             if (C == 0) { Noff = 1  }  else {Noff = sample(Familysize,1)}
               for (i in 1:Noff)  {
               OffNumber = OffNumber + 1
               S1[OffNumber,1] = A        ## record parents of each progeny
               S1[OffNumber,2] = B
             for (j in 1:NLoci)  {  ## get offspring genotypes for each locus
                if (parents[A,j]==1) { ## heterozygote, so pick an allele
                   Offspring[OffNumber,j] = Offspring[OffNumber,j] + rbinom(1,1,0.5)  } else {  ## homozygote, so take half of the genes
                   Offspring[OffNumber,j] = Offspring[OffNumber,j] + 0.5*parents[A,j]  }  # end first parent
                 if (parents[B,j]==1) { ## heterozygote, so pick an allele
                   Offspring[OffNumber,j] = Offspring[OffNumber,j] + rbinom(1,1,0.5)  } else {  ## homozygote, so take half of the genes
                   Offspring[OffNumber,j] = Offspring[OffNumber,j] + 0.5*parents[B,j]  }  # end first parent
                   }  ## end for j
                   }  ## end for i
                   } ## end while
      OffspringData = cbind(S1,Offspring)
      OffspringData = OffspringData[1:S,]  ## truncate offspring to desired sample size

    return(OffspringData)            } # end function




#### Initialize population as 100% heterozygotes so starting P = 0.5 at each locus ####
Parents0 = matrix(1,Ne,NLoci)

## We are simulating a Wright-Fisher ideal population with separate sexes

Parents1 = Parents0

##simulate NGens generations of random mating and genetic drift
for (k in 1:NGens)    {
###Get genotypes in generation k
  if (Mating == 1) {Parents2 = Separate(Parents1) } else if (Mating == 2) {
     Parents2 = Monogamy(Parents1) } else {
     Parents2 = Mixed(Parents1) }
Parents2 = Parents2[sample(nrow(Parents2), nrow(Parents2)),]  ##scramble the rows (individuals)
## offpsring become parents of next generation
Parents1 = Parents2
                      } # end for k

TrueP = colMeans(Parents1)/2  ## get allele frequencies in parents of sample
PP = TrueP*(1-TrueP)


#### Cycle over the Reps ####
for (jj in 1:NReps) {

###########################
## Produce a sample from the final generation
 if (Mating == 1) {Offspring = RanFamily(Parents1)  } else if (Mating == 2) {
     Offspring = MonoFamily(Parents1) } else {
     Offspring = MixedFamily(Parents1) }

Parentlist = Offspring[,1:2]
Genos = Offspring[,3:(NLoci+2)]
SampP = colMeans(Genos)/2  ## allele frequencies in total sample of progeny
RPP = RPP + (TrueP-SampP)^2
SSP = SSP + SampP^2
PBar = PBar + SampP

Sibs = matrix(0,S,S)
for (i in 1:(S-1))  {   ## Get matrix of sibships
for (j in (i+1):S)  {
   if (Parentlist[i,1]==Parentlist[j,1] | Parentlist[i,1]==Parentlist[j,2]) { Sibs[i,j] = Sibs[i,j]+1 }
   if (Parentlist[i,2]==Parentlist[j,1] | Parentlist[i,2]==Parentlist[j,2]) { Sibs[i,j] = Sibs[i,j]+1 }
   }  }  ## end for i and j

###############################################
##Implement Yank2 removal; remove individuals only if 2 full sibs are already in the sample

ReducedOffspring = Offspring[1,]
ilist = c(1)

for (i in 2:S)  {
sibyes = 0
for (j in 1:length(ilist)) {
  jk = ilist[j]
  if (Sibs[jk,i]>Sibcheck) { sibyes = sibyes + 1}
  }  # end for j
  if (sibyes < 2) {
    ilist = append(ilist,i)
    ReducedOffspring = rbind(ReducedOffspring,Offspring[i,])  }
  } # end for i

ReducedParentlistYank = ReducedOffspring[,1:2]
ReducedGenos = ReducedOffspring[,3:(NLoci+2)]
ReducedSampP = colMeans(ReducedGenos)/2  ## allele frequencies in sample of progeny that excludes siblings
RPPYank = RPPYank + (TrueP-ReducedSampP)^2
SSPYank = SSPYank + ReducedSampP^2
PBarYank = PBarYank + ReducedSampP

ReducedSYank = dim(ReducedGenos)[[1]]
SSYank[jj] = ReducedSYank


###############################################
##With probability Beta, exclude individuals who are siblings of another individual already in the sample

ReducedOffspring = Offspring[1,]
Beta = 0.25
ilist = c(1)

for (i in 2:S)  {
sibyes = 0
for (j in 1:length(ilist)) {
  jk = ilist[j]
  if (Sibs[jk,i]>Sibcheck) { sibyes = 1}
  }  # end for j
  if (sibyes == 0) {
    ilist = append(ilist,i)
    ReducedOffspring = rbind(ReducedOffspring,Offspring[i,])  } else {
  D = rbinom(1,1,Beta)
  if (D == 0) {ReducedOffspring = rbind(ReducedOffspring,Offspring[i,])  }  ## keep the sib with prob 1-Beta
  }  ## end else
  } # end for i

ReducedParentlist = ReducedOffspring[,1:2]
ReducedGenos = ReducedOffspring[,3:(NLoci+2)]
ReducedSampP = colMeans(ReducedGenos)/2  ## allele frequencies in sample of progeny that excludes siblings
RPP1 = RPP1 + (TrueP-ReducedSampP)^2
SSP1 = SSP1 + ReducedSampP^2
PBar1 = PBar1 + ReducedSampP

ReducedS = dim(ReducedGenos)[[1]]
SS1[jj] = ReducedS

######
ReducedOffspring = Offspring[1,]
Beta = 0.5
ilist = c(1)
for (i in 2:S)  {
sibyes = 0
for (j in 1:length(ilist)) {
  jk = ilist[j]
  if (Sibs[jk,i]>Sibcheck) { sibyes = 1}
  }  # end for j
  if (sibyes == 0) {
    ilist = append(ilist,i)
    ReducedOffspring = rbind(ReducedOffspring,Offspring[i,])  } else {
  D = rbinom(1,1,Beta)
  if (D == 0) {ReducedOffspring = rbind(ReducedOffspring,Offspring[i,])  }  ## keep the sib with prob 1-Beta
  }  ## end else
  } # end for i

ReducedParentlist = ReducedOffspring[,1:2]
ReducedGenos = ReducedOffspring[,3:(NLoci+2)]
ReducedSampP = colMeans(ReducedGenos)/2  ## allele frequencies in sample of progeny that excludes siblings

RPP2 = RPP2 + (TrueP-ReducedSampP)^2
SSP2 = SSP2 + ReducedSampP^2
PBar2 = PBar2 + ReducedSampP
ReducedS = dim(ReducedGenos)[[1]]
SS2[jj] = ReducedS

####
ReducedOffspring = Offspring[1,]
Beta = 0.75
ilist = c(1)
for (i in 2:S)  {
sibyes = 0
for (j in 1:length(ilist)) {
  jk = ilist[j]
  if (Sibs[jk,i]>Sibcheck) { sibyes = 1}
  }  # end for j
  if (sibyes == 0) {
    ilist = append(ilist,i)
    ReducedOffspring = rbind(ReducedOffspring,Offspring[i,])  } else {
  D = rbinom(1,1,Beta)
  if (D == 0) {ReducedOffspring = rbind(ReducedOffspring,Offspring[i,])  }  ## keep the sib with prob 1-Beta
  }  ## end else
  } # end for i

ReducedParentlist = ReducedOffspring[,1:2]
ReducedGenos = ReducedOffspring[,3:(NLoci+2)]
ReducedSampP = colMeans(ReducedGenos)/2  ## allele frequencies in sample of progeny that excludes siblings

RPP3 = RPP3 + (TrueP-ReducedSampP)^2
SSP3 = SSP3 + ReducedSampP^2
PBar3 = PBar3 + ReducedSampP
ReducedS = dim(ReducedGenos)[[1]]
SS3[jj] = ReducedS

####
ReducedOffspring = Offspring[1,]
Beta = 0.9
ilist = c(1)
for (i in 2:S)  {
sibyes = 0
for (j in 1:length(ilist)) {
  jk = ilist[j]
  if (Sibs[jk,i]>Sibcheck) { sibyes = 1}
  }  # end for j
  if (sibyes == 0) {
    ilist = append(ilist,i)
    ReducedOffspring = rbind(ReducedOffspring,Offspring[i,])  } else {
  D = rbinom(1,1,Beta)
  if (D == 0) {ReducedOffspring = rbind(ReducedOffspring,Offspring[i,])  }  ## keep the sib with prob 1-Beta
  }  ## end else
  } # end for i

ReducedParentlist = ReducedOffspring[,1:2]
ReducedGenos = ReducedOffspring[,3:(NLoci+2)]
ReducedSampP = colMeans(ReducedGenos)/2  ## allele frequencies in sample of progeny that excludes siblings

RPP4 = RPP4 + (TrueP-ReducedSampP)^2
SSP4 = SSP4 + ReducedSampP^2
PBar4 = PBar4 + ReducedSampP
ReducedS = dim(ReducedGenos)[[1]]
SS4[jj] = ReducedS

####
ReducedOffspring = Offspring[1,]
Beta = 0.95
ilist = c(1)
for (i in 2:S)  {
sibyes = 0
for (j in 1:length(ilist)) {
  jk = ilist[j]
  if (Sibs[jk,i]>Sibcheck) { sibyes = 1}
  }  # end for j
  if (sibyes == 0) {
    ilist = append(ilist,i)
    ReducedOffspring = rbind(ReducedOffspring,Offspring[i,])  } else {
  D = rbinom(1,1,Beta)
  if (D == 0) {ReducedOffspring = rbind(ReducedOffspring,Offspring[i,])  }  ## keep the sib with prob 1-Beta
  }  ## end else
  } # end for i

ReducedParentlist = ReducedOffspring[,1:2]
ReducedGenos = ReducedOffspring[,3:(NLoci+2)]
ReducedSampP = colMeans(ReducedGenos)/2  ## allele frequencies in sample of progeny that excludes siblings

RPP5 = RPP5 + (TrueP-ReducedSampP)^2
SSP5 = SSP5 + ReducedSampP^2
PBar5 = PBar5 + ReducedSampP
ReducedS = dim(ReducedGenos)[[1]]
SS5[jj] = ReducedS

####
ReducedOffspring = Offspring[1,]
Beta = 1
ilist = c(1)
for (i in 2:S)  {
sibyes = 0
for (j in 1:length(ilist)) {
  jk = ilist[j]
  if (Sibs[jk,i]>Sibcheck) { sibyes = 1}
  }  # end for j
  if (sibyes == 0) {
    ilist = append(ilist,i)
    ReducedOffspring = rbind(ReducedOffspring,Offspring[i,])  } else {
  D = rbinom(1,1,Beta)
  if (D == 0) {ReducedOffspring = rbind(ReducedOffspring,Offspring[i,])  }  ## keep the sib with prob 1-Beta
  }  ## end else
  } # end for i

ReducedParentlist = ReducedOffspring[,1:2]
ReducedGenos = ReducedOffspring[,3:(NLoci+2)]
ReducedSampP = colMeans(ReducedGenos)/2  ## allele frequencies in sample of progeny that excludes siblings
SSP6 = SSP6 + ReducedSampP^2
PBar6 = PBar6 + ReducedSampP
RPP6 = RPP6 + (TrueP-ReducedSampP)^2
ReducedS = dim(ReducedGenos)[[1]]
SS6[jj] = ReducedS

###############################################

##print(paste0("Replicate = ",jj))
##flush.console()

#### DO SOME BLUE STUFF

# first get the pedigree of the offspring
MyGenos <- Offspring[, -(1:2)]

ped <- as.data.frame(Offspring[,1:2]) %>%
  tbl_df() %>%
  mutate(id = paste("kid_", 1:n(), sep = ""),
         mom = paste("ma_", V1, sep = ""),
         dad = paste("pa_", V2, sep = ""),
         sex = "U") %>%
  select(id, mom, dad, sex)
rownames(MyGenos) <- ped$id

# then add the parents to the pedigree
full_ped <- bind_rows(ped,
          data_frame(id = unique(ped$mom), mom = NA, dad = NA, sex = "F"),
          data_frame(id = unique(ped$dad), mom = NA, dad = NA, sex = "M")
          )
# make the kinship matrix
Lmat <- matrix_L_from_pedigree(full_ped, ped$id)

# compute the weights for the kids
kid_wts <- weights_from_matrix_L(Lmat)

# then estimate P from that
ReducedSampP = colSums(MyGenos * kid_wts) / (2 * sum(kid_wts))  ## allele frequencies in sample of progeny weighted by the BLUE
SSP_BLUE = SSP_BLUE + ReducedSampP^2
PBarBLUE = PBarBLUE + ReducedSampP
RPP_BLUE = RPP_BLUE + (TrueP-ReducedSampP)^2
ReducedS = sum(kid_wts)  # this is actually the "effective sample size under BLUE"
SS_BLUE[jj] = ReducedS


}  ## end for jj


HMS = 1/(mean(1/SS0))
HMS1 = 1/(mean(1/SS1))
HMS2 = 1/(mean(1/SS2))
HMS3 = 1/(mean(1/SS3))
HMS4 = 1/(mean(1/SS4))
HMS5 = 1/(mean(1/SS5))
HMS6 = 1/(mean(1/SS6))
HMSYank = 1/(mean(1/SSYank))
HMSBlue <- 1/(mean(1/SS_BLUE))

AllS = c(HMS, HMSBlue, HMSYank, HMS1, HMS2, HMS3, HMS4, HMS5, HMS6)

MSE = RPP/NReps
RMSE = sqrt(MSE)
PBar = PBar/NReps
SSP = SSP/NReps
VarP = (SSP - PBar^2)*NReps/(NReps-1)
ES = PP*(NReps-1)/(2*MSE*NReps)
ES = 1/ES
ESS = 1/mean(ES,na.rm=TRUE)

MSEYank = RPPYank/NReps
RMSEYank = sqrt(MSEYank)
PBarYank = PBarYank/NReps
SSPYank = SSPYank/NReps
VarPYank = (SSPYank - PBarYank^2)*NReps/(NReps-1)
ESYank = PP*(NReps-1)/(2*MSEYank*NReps)
ESYank = 1/ESYank
ESSYank = 1/mean(ESYank,na.rm=TRUE)


MSEBlue = RPP_BLUE/NReps
RMSEBlue = sqrt(MSEBlue)
PBarBlue = PBarBLUE/NReps
SSPBlue = SSP_BLUE/NReps
VarPBlue = (SSPBlue - PBarBlue^2)*NReps/(NReps-1)
ESBlue = PP*(NReps-1)/(2*MSEBlue*NReps)
ESBlue = 1/ESBlue
ESSBlue = 1/mean(ESBlue,na.rm=TRUE)



MSE1 = RPP1/NReps
RMSE1 = sqrt(MSE1)
PBar1 = PBar1/NReps
SSP1 = SSP1/NReps
VarP1 = (SSP1 - PBar1^2)*NReps/(NReps-1)
ES1 = PP*(NReps-1)/(2*MSE1*NReps)
ES1 = 1/ES1
ESS1 = 1/mean(ES1,na.rm=TRUE)

MSE2 = RPP2/NReps
RMSE2 = sqrt(MSE2)
PBar2 = PBar2/NReps
SSP2 = SSP2/NReps
VarP2 = (SSP2 - PBar2^2)*NReps/(NReps-1)
ES2 = PP*(NReps-1)/(2*MSE2*NReps)
ES2 = 1/ES2
ESS2 = 1/mean(ES2,na.rm=TRUE)

MSE3 = RPP3/NReps
RMSE3 = sqrt(MSE3)
PBar3 = PBar3/NReps
SSP3 = SSP3/NReps
VarP3 = (SSP3 - PBar3^2)*NReps/(NReps-1)
ES3 = PP*(NReps-1)/(2*MSE3*NReps)
ES3 = 1/ES3
ESS3 = 1/mean(ES3,na.rm=TRUE)

MSE4 = RPP4/NReps
RMSE4 = sqrt(MSE4)
PBar4 = PBar4/NReps
SSP4 = SSP4/NReps
VarP4 = (SSP4 - PBar4^2)*NReps/(NReps-1)
ES4 = PP*(NReps-1)/(2*MSE4*NReps)
ES4 = 1/ES4
ESS4 = 1/mean(ES4,na.rm=TRUE)

MSE5 = RPP5/NReps
RMSE5 = sqrt(MSE5)
PBar5 = PBar5/NReps
SSP5 = SSP5/NReps
VarP5 = (SSP5 - PBar5^2)*NReps/(NReps-1)
ES5 = PP*(NReps-1)/(2*MSE5*NReps)
ES5 = 1/ES5
ESS5 = 1/mean(ES5,na.rm=TRUE)

MSE6 = RPP6/NReps
RMSE6 = sqrt(MSE6)
PBar6 = PBar6/NReps
SSP6 = SSP6/NReps
VarP6 = (SSP6 - PBar6^2)*NReps/(NReps-1)
ES6 = PP*(NReps-1)/(2*MSE6*NReps)
ES6 = 1/ES6
ESS6 = 1/mean(ES6,na.rm=TRUE)

AllESS = c(ESS, ESSBlue, ESSYank, ESS1, ESS2, ESS3, ESS4, ESS5, ESS6)

RRMSE = 1

RRMSEYank = RMSEYank/RMSE
RRMSEBlue = RMSEBlue/RMSE
RRMSE1 = RMSE1/RMSE
RRMSE2 = RMSE2/RMSE
RRMSE3 = RMSE3/RMSE
RRMSE4 = RMSE4/RMSE
RRMSE5 = RMSE5/RMSE
RRMSE6 = RMSE6/RMSE
YRRMSE = gm_mean(RRMSEYank)
BLUE_RRMSE = gm_mean(RRMSEBlue)
ARRMSE = gm_mean(RRMSE1)
BRRMSE = gm_mean(RRMSE2)
CRRMSE = gm_mean(RRMSE3)
DRRMSE = gm_mean(RRMSE4)
ERRMSE = gm_mean(RRMSE5)
FRRMSE = gm_mean(RRMSE6)

AllRRMSE = c(RRMSE, BLUE_RRMSE, YRRMSE, ARRMSE,BRRMSE,CRRMSE,DRRMSE,ERRMSE,FRRMSE)

All = cbind(AllRRMSE,AllS,AllESS)
rownames(All) = c("Full","BLUE", "Yank2","25","50","75","90","95","100")
All
cat(AllESS,sep="\n")
