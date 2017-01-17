### This code simulates all three mating models and also adds Yank2 strategy
### It calculates RRMSE for allele frequency, Fst, and Ne
### For calculation of effective sample size (ESS), see file "ESS.R"



NLoci = 100
Ne = 40
S = 100
NGens = 10
NReps = 10
MaxFamily = 9
ProbFamily = 0.5
Familysize = seq(1:MaxFamily)
Sibcheck = 1 ## 0 purges FS+HS; 1 purges only FS
Mating = 2  ## 1 = random; 2 = monogamy; 3 = mixed
     if (Mating==3) {Familysize = seq(2:MaxFamily)+1}
MaxSib = seq(1:NReps)
P2 = rep(0,NLoci)
SSP = rep(0,NLoci)
PBar = rep(0,NLoci)
PBarYank = rep(0,NLoci)
SSPYank = rep(0,NLoci)
RPPYank = rep(0,NLoci)
PBar6 = rep(0,NLoci)
SSP6 = rep(0,NLoci)
NeR = seq(1:NReps)
NeRYank = seq(1:NReps)
ReducedNeR = seq(1:NReps)
SamplePP = seq(1:NReps)
ReducedPP = seq(1:NReps)
PWOP = seq(1:NReps)
ReducedNeR1 = seq(1:NReps)
ReducedNeR2 = seq(1:NReps)
ReducedNeR3 = seq(1:NReps)
ReducedNeR4 = seq(1:NReps)
ReducedNeR5 = seq(1:NReps)
ReducedNeR6 = seq(1:NReps)
ReducedPPYank = seq(1:NReps)
ReducedPP1 = seq(1:NReps)
ReducedPP2 = seq(1:NReps)
ReducedPP3 = seq(1:NReps)
ReducedPP4 = seq(1:NReps)
ReducedPP5 = seq(1:NReps)
ReducedPP6 = seq(1:NReps)
PWOPYank = seq(1:NReps)
PWOP1 = seq(1:NReps)
PWOP2 = seq(1:NReps)
PWOP3 = seq(1:NReps)
PWOP4 = seq(1:NReps)
PWOP5 = seq(1:NReps)
PWOP6 = seq(1:NReps)
SS1 = seq(1:NReps)
SS2 = seq(1:NReps)
SS3 = seq(1:NReps)
SS4 = seq(1:NReps)
SS5 = seq(1:NReps)
SS6 = seq(1:NReps)
SS0 = rep(S,NReps)
SSYank = seq(1:NReps)
HS = seq(1:NReps)
FS = seq(1:NReps)
FSTP = seq(1:(NReps/2))
FSTO = seq(1:(NReps/2))
FSTYank = seq(1:(NReps/2))
FSTO1 = seq(1:(NReps/2))
FSTO2 = seq(1:(NReps/2))
FSTO3 = seq(1:(NReps/2))
FSTO4 = seq(1:(NReps/2))
FSTO5 = seq(1:(NReps/2))
FSTO6 = seq(1:(NReps/2))
FPP = seq(1:(NReps/2))
FPPYank = seq(1:(NReps/2))
FPP1 = seq(1:(NReps/2))
FPP2 = seq(1:(NReps/2))
FPP3 = seq(1:(NReps/2))
FPP4 = seq(1:(NReps/2))
FPP5 = seq(1:(NReps/2))
FPP6 = seq(1:(NReps/2))
QQ = seq(1:NLoci)
QQQ = seq(1:NLoci)

is.odd <- function(x) x %% 2 != 0

###################################################################
#### three mating models for first 10 generations #####

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


##############################
########## production of progeny, with random (MaxFamily = 1) or family-correlated samples (MaxFamily > 1)

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


#################################################################

###############################################################
GetLD <- function(Geno)     {
 UsedLoci = G
 rsq = matrix(NA,UsedLoci,UsedLoci)
 for (i in 1:(UsedLoci-1))     {
 for (j in (i+1):UsedLoci)  {
      if(var(Geno[,i])*var(Geno[,j])>0) {rsq[i,j] = cor(Geno[,i],Geno[,j])}   } } # end for i and j
 rsq = rsq^2
 rmean = mean(rsq,na.rm=TRUE)
return(rmean)      } # end function
######################################################################

###############################################################
GetFst <- function(yy)     {
OldP = yy[1,]
P2 = yy[2,]
SumP1 = 1 - (P2^2 + (1-P2)^2)
SumP2 = 1 - (OldP^2 + (1-OldP)^2)
Hexp = (SumP1+SumP2)/2
Pbar = (P2 + OldP)/2
Htot = 1 - Pbar^2 - (1-Pbar)^2
F = 1 - Hexp/Htot
Fst = mean(F,na.rm=TRUE)
return(Fst)  } # end function
###############################################################


##Initialize population as 100% heterozygotes so starting P = 0.5 at each locus
Parents0 = matrix(1,Ne,NLoci)

## Give initial population a distribution of allele frequences 0.2 - 0.5
x <- as.integer(0.6*Ne)
for (j in 1:NLoci)     {
   y = sample.int(x+1,1, replace = TRUE) - 1
   if (y > 0)         {
   Parents0[1:y,j] = 0   } # end if
                       } # end j

##scramble the rows (individuals) to randomize initial allele freqs for samples
Parents0 = Parents0[sample(nrow(Parents0), nrow(Parents0)),]

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

for (jj in 1:NReps) {

###########################
## Produce a sample from the final generation
 if (Mating == 1) {Offspring = RanFamily(Parents1)  } else if (Mating == 2) {
     Offspring = MonoFamily(Parents1) } else {
     Offspring = MixedFamily(Parents1) }

Parentlist = Offspring[,1:2]
Genos = Offspring[,3:(NLoci+2)]
SampP = colMeans(Genos)/2  ## allele frequencies in total sample of progeny
SPP = (TrueP-SampP)^2
P2 = P2+SPP
SSP = SSP + SampP^2
PBar = PBar + SampP
SamplePP[jj] = sqrt(mean(SPP))

zz = Genos[, SampP > 0.05 & SampP < 0.95, drop = F]
G = dim(zz)[[2]]
MeanRsq = GetLD(zz)
Rprime = MeanRsq - 1/(S-1)
if (S > 29) {NeR[jj] = (2/3 + sqrt(4/9-7.2*Rprime))/(2*Rprime)} else {NeR[jj] = (0.618 + sqrt(0.618^2-5.24*Rprime))/(2*Rprime)}

if (is.odd(jj)==TRUE) { ## then jj is odd, so skip Fst and record data for pop 1
    OldP = TrueP
    OldSampP = SampP  } else {  ## else = even, so compute Fst
  combo1 = rbind(TrueP, OldP)
  combo2 = rbind(SampP, OldSampP)
##eliminate loci monomorphic in both samples
for (i in 1:NLoci)  {
QQ[i] = combo1[1,i]+combo1[2,i]
QQQ[i] = combo2[1,i]+combo2[2,i]}
FP1 = combo1[,QQ<2 & QQ > 0]
FP2 = combo2[,QQQ<2 & QQQ > 0]
FSTP[jj/2] = GetFst(FP1)  ##True Fst in the parents
FSTO[jj/2] = GetFst(FP2) - 1/(2*S)  ##Fst in the sample
FPP[jj/2] = (FSTP[jj/2]-FSTO[jj/2])^2
} #end else

MaleSuccess = table(Parentlist[,1])
A = sum(MaleSuccess)
B = sum(MaleSuccess^2)
if (B>A) {MalePWOP = (A-1)/(B/A-1)} else {MalePWOP = 99999}
FemaleSuccess = table(Parentlist[,2])
C = sum(FemaleSuccess)
D = sum(FemaleSuccess^2)
if (D>C) {FemalePWOP = (C-1)/(D/C-1)} else {FemalePWOP = 99999}
PWOP[jj] = 4*MalePWOP*FemalePWOP/(MalePWOP + FemalePWOP) ## effective number of parents for sample

Sibs = matrix(0,S,S)
for (i in 1:(S-1))  {   ## Get matrix of sibships
for (j in (i+1):S)  {
   if (Parentlist[i,1]==Parentlist[j,1] | Parentlist[i,1]==Parentlist[j,2]) { Sibs[i,j] = Sibs[i,j]+1 }
   if (Parentlist[i,2]==Parentlist[j,1] | Parentlist[i,2]==Parentlist[j,2]) { Sibs[i,j] = Sibs[i,j]+1 }
   }  }  ## end for i and j

H = sum(Sibs == 1)
F = sum(Sibs == 2)
npairs = S*(S-1)/2
HS[jj] = H/npairs
FS[jj] = F/npairs
MaxSib[jj] = max(max(MaleSuccess),max(FemaleSuccess))

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
ReducedSampPYank = colMeans(ReducedGenos)/2  ## allele frequencies in sample of progeny that excludes siblings
ReducedSYank = dim(ReducedGenos)[[1]]
SSYank[jj] = ReducedSYank
RPPYank = (TrueP-ReducedSampPYank)^2
ReducedPPYank[jj] = sqrt(mean(RPPYank))

MaleSuccess = table(ReducedParentlistYank[,1])
A = sum(MaleSuccess)
B = sum(MaleSuccess^2)
if (B>A) {MalePWOP = (A-1)/(B/A-1)} else {MalePWOP = 99999}
FemaleSuccess = table(ReducedParentlistYank[,2])
C = sum(FemaleSuccess)
D = sum(FemaleSuccess^2)
if (D>C) {FemalePWOP = (C-1)/(D/C-1)} else {FemalePWOP = 99999}
PWOPYank[jj] = 4*MalePWOP*FemalePWOP/(MalePWOP + FemalePWOP) ## effective number of parents for sample

zzz = ReducedGenos[, ReducedSampPYank > 0.05 & ReducedSampPYank < 0.95, drop = F]
G = dim(zzz)[[2]]
MeanRsq2 = GetLD(zzz)
Rprime = MeanRsq2 - 1/(ReducedSYank-1)
if (ReducedSYank > 29) {NeRYank[jj] = (2/3 + sqrt(4/9-7.2*Rprime))/(2*Rprime)} else {NeRYank[jj] = (0.618 + sqrt(0.618^2-5.24*Rprime))/(2*Rprime)}
if (is.odd(jj)==TRUE) { ## then jj is odd, so skip Fst and record data for pop 1
    OldReducedSampPYank = ReducedSampPYank  } else {  ## else = even, so compute Fst
  combo2 = rbind(OldReducedSampPYank,ReducedSampPYank)
##eliminate loci monomorphic in both samples
for (i in 1:NLoci)  {
QQQ[i] = combo2[1,i]+combo2[2,i]}
FP2 = combo2[,QQQ<2 & QQQ > 0]
FSTYank[jj/2] = GetFst(FP2) - 1/(2*ReducedSYank)
FPPYank[jj/2] = (FSTP[jj/2]-FSTYank[jj/2])^2
} #end else


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

MaleSuccess = table(ReducedParentlist[,1])
A = sum(MaleSuccess)
B = sum(MaleSuccess^2)
if (B>A) {MalePWOP = (A-1)/(B/A-1)} else {MalePWOP = 99999}
FemaleSuccess = table(ReducedParentlist[,2])
C = sum(FemaleSuccess)
D = sum(FemaleSuccess^2)
if (D>C) {FemalePWOP = (C-1)/(D/C-1)} else {FemalePWOP = 99999}
PWOP1[jj] = 4*MalePWOP*FemalePWOP/(MalePWOP + FemalePWOP) ## effective number of parents for sample

RPP = (TrueP-ReducedSampP)^2
ReducedPP1[jj] = sqrt(mean(RPP))

ReducedS = dim(ReducedGenos)[[1]]
zzz = ReducedGenos[, ReducedSampP > 0.05 & ReducedSampP < 0.95, drop = F]
G = dim(zzz)[[2]]
MeanRsq2 = GetLD(zzz)
Rprime = MeanRsq2 - 1/(ReducedS-1)
if (ReducedS > 29) {ReducedNeR1[jj] = (2/3 + sqrt(4/9-7.2*Rprime))/(2*Rprime)} else {ReducedNeR1[jj] = (0.618 + sqrt(0.618^2-5.24*Rprime))/(2*Rprime)}
SS1[jj] = ReducedS

if (is.odd(jj)==TRUE) { ## then jj is odd, so skip Fst and record data for pop 1
    OldReducedSampP1 = ReducedSampP  } else {  ## else = even, so compute Fst
  combo2 = rbind(OldReducedSampP1,ReducedSampP)
##eliminate loci monomorphic in both samples
for (i in 1:NLoci)  {
QQQ[i] = combo2[1,i]+combo2[2,i]}
FP2 = combo2[,QQQ<2 & QQQ > 0]
FSTO1[jj/2] = GetFst(FP2) - 1/(2*ReducedS)
FPP1[jj/2] = (FSTP[jj/2]-FSTO1[jj/2])^2
} #end else

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

MaleSuccess = table(ReducedParentlist[,1])
A = sum(MaleSuccess)
B = sum(MaleSuccess^2)
if (B>A) {MalePWOP = (A-1)/(B/A-1)} else {MalePWOP = 99999}
FemaleSuccess = table(ReducedParentlist[,2])
C = sum(FemaleSuccess)
D = sum(FemaleSuccess^2)
if (D>C) {FemalePWOP = (C-1)/(D/C-1)} else {FemalePWOP = 99999}
PWOP2[jj] = 4*MalePWOP*FemalePWOP/(MalePWOP + FemalePWOP) ## effective number of parents for sample

RPP = (TrueP-ReducedSampP)^2
ReducedPP2[jj] = sqrt(mean(RPP))

ReducedS = dim(ReducedGenos)[[1]]
zzz = ReducedGenos[, ReducedSampP > 0.05 & ReducedSampP < 0.95, drop = F]
G = dim(zzz)[[2]]
MeanRsq2 = GetLD(zzz)
Rprime = MeanRsq2 - 1/(ReducedS-1)
if (ReducedS > 29) {ReducedNeR2[jj] = (2/3 + sqrt(4/9-7.2*Rprime))/(2*Rprime)} else {ReducedNeR2[jj] = (0.618 + sqrt(0.618^2-5.24*Rprime))/(2*Rprime)}
SS2[jj] = ReducedS

if (is.odd(jj)==TRUE) { ## then jj is odd, so skip Fst and record data for pop 1
    OldReducedSampP2 = ReducedSampP  } else {  ## else = even, so compute Fst
  combo2 = rbind(OldReducedSampP2,ReducedSampP)
##eliminate loci monomorphic in both samples
for (i in 1:NLoci)  {
QQQ[i] = combo2[1,i]+combo2[2,i]}
FP2 = combo2[,QQQ<2 & QQQ > 0]
FSTO2[jj/2] = GetFst(FP2) - 1/(2*ReducedS)
FPP2[jj/2] = (FSTP[jj/2]-FSTO2[jj/2])^2
} #end else

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

MaleSuccess = table(ReducedParentlist[,1])
A = sum(MaleSuccess)
B = sum(MaleSuccess^2)
if (B>A) {MalePWOP = (A-1)/(B/A-1)} else {MalePWOP = 99999}
FemaleSuccess = table(ReducedParentlist[,2])
C = sum(FemaleSuccess)
D = sum(FemaleSuccess^2)
if (D>C) {FemalePWOP = (C-1)/(D/C-1)} else {FemalePWOP = 99999}
PWOP3[jj] = 4*MalePWOP*FemalePWOP/(MalePWOP + FemalePWOP) ## effective number of parents for sample

RPP = (TrueP-ReducedSampP)^2
ReducedPP3[jj] = sqrt(mean(RPP))

ReducedS = dim(ReducedGenos)[[1]]
zzz = ReducedGenos[, ReducedSampP > 0.05 & ReducedSampP < 0.95, drop = F]
G = dim(zzz)[[2]]
MeanRsq2 = GetLD(zzz)
Rprime = MeanRsq2 - 1/(ReducedS-1)
if (ReducedS > 29) {ReducedNeR3[jj] = (2/3 + sqrt(4/9-7.2*Rprime))/(2*Rprime)} else {ReducedNeR3[jj] = (0.618 + sqrt(0.618^2-5.24*Rprime))/(2*Rprime)}
SS3[jj] = ReducedS

if (is.odd(jj)==TRUE) { ## then jj is odd, so skip Fst and record data for pop 1
    OldReducedSampP3 = ReducedSampP  } else {  ## else = even, so compute Fst
  combo2 = rbind(OldReducedSampP3,ReducedSampP)
##eliminate loci monomorphic in both samples
for (i in 1:NLoci)  {
QQQ[i] = combo2[1,i]+combo2[2,i]}
FP2 = combo2[,QQQ<2 & QQQ > 0]
FSTO3[jj/2] = GetFst(FP2) - 1/(2*ReducedS)
FPP3[jj/2] = (FSTP[jj/2]-FSTO3[jj/2])^2
} #end else

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

MaleSuccess = table(ReducedParentlist[,1])
A = sum(MaleSuccess)
B = sum(MaleSuccess^2)
if (B>A) {MalePWOP = (A-1)/(B/A-1)} else {MalePWOP = 99999}
FemaleSuccess = table(ReducedParentlist[,2])
C = sum(FemaleSuccess)
D = sum(FemaleSuccess^2)
if (D>C) {FemalePWOP = (C-1)/(D/C-1)} else {FemalePWOP = 99999}
PWOP4[jj] = 4*MalePWOP*FemalePWOP/(MalePWOP + FemalePWOP) ## effective number of parents for sample

RPP = (TrueP-ReducedSampP)^2
ReducedPP4[jj] = sqrt(mean(RPP))

ReducedS = dim(ReducedGenos)[[1]]
zzz = ReducedGenos[, ReducedSampP > 0.05 & ReducedSampP < 0.95, drop = F]
G = dim(zzz)[[2]]
MeanRsq2 = GetLD(zzz)
Rprime = MeanRsq2 - 1/(ReducedS-1)

if (ReducedS > 29) {ReducedNeR4[jj] = (2/3 + sqrt(4/9-7.2*Rprime))/(2*Rprime)} else {ReducedNeR4[jj] = (0.618 + sqrt(0.618^2-5.24*Rprime))/(2*Rprime)}
SS4[jj] = ReducedS

if (is.odd(jj)==TRUE) { ## then jj is odd, so skip Fst and record data for pop 1
    OldReducedSampP4 = ReducedSampP  } else {  ## else = even, so compute Fst
  combo2 = rbind(OldReducedSampP4,ReducedSampP)
##eliminate loci monomorphic in both samples
for (i in 1:NLoci)  {
QQQ[i] = combo2[1,i]+combo2[2,i]}
FP2 = combo2[,QQQ<2 & QQQ > 0]
FSTO4[jj/2] = GetFst(FP2) - 1/(2*ReducedS)
FPP4[jj/2] = (FSTP[jj/2]-FSTO4[jj/2])^2
} #end else

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

MaleSuccess = table(ReducedParentlist[,1])
A = sum(MaleSuccess)
B = sum(MaleSuccess^2)
if (B>A) {MalePWOP = (A-1)/(B/A-1)} else {MalePWOP = 99999}
FemaleSuccess = table(ReducedParentlist[,2])
C = sum(FemaleSuccess)
D = sum(FemaleSuccess^2)
if (D>C) {FemalePWOP = (C-1)/(D/C-1)} else {FemalePWOP = 99999}
PWOP5[jj] = 4*MalePWOP*FemalePWOP/(MalePWOP + FemalePWOP) ## effective number of parents for sample

RPP = (TrueP-ReducedSampP)^2
ReducedPP5[jj] = sqrt(mean(RPP))

ReducedS = dim(ReducedGenos)[[1]]
zzz = ReducedGenos[, ReducedSampP > 0.05 & ReducedSampP < 0.95, drop = F]
G = dim(zzz)[[2]]
MeanRsq2 = GetLD(zzz)
Rprime = MeanRsq2 - 1/(ReducedS-1)
if (ReducedS > 29) {ReducedNeR5[jj] = (2/3 + sqrt(4/9-7.2*Rprime))/(2*Rprime)} else {ReducedNeR5[jj] = (0.618 + sqrt(0.618^2-5.24*Rprime))/(2*Rprime)}
SS5[jj] = ReducedS

if (is.odd(jj)==TRUE) { ## then jj is odd, so skip Fst and record data for pop 1
    OldReducedSampP5 = ReducedSampP  } else {  ## else = even, so compute Fst
  combo2 = rbind(OldReducedSampP5,ReducedSampP)
##eliminate loci monomorphic in both samples
for (i in 1:NLoci)  {
QQQ[i] = combo2[1,i]+combo2[2,i]}
FP2 = combo2[,QQQ<2 & QQQ > 0]
FSTO5[jj/2] = GetFst(FP2) - 1/(2*ReducedS)
FPP5[jj/2] = (FSTP[jj/2]-FSTO5[jj/2])^2
} #end else

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

MaleSuccess = table(ReducedParentlist[,1])
A = sum(MaleSuccess)
B = sum(MaleSuccess^2)
if (B>A) {MalePWOP = (A-1)/(B/A-1)} else {MalePWOP = 99999}
FemaleSuccess = table(ReducedParentlist[,2])
C = sum(FemaleSuccess)
D = sum(FemaleSuccess^2)
if (D>C) {FemalePWOP = (C-1)/(D/C-1)} else {FemalePWOP = 99999}
PWOP6[jj] = 4*MalePWOP*FemalePWOP/(MalePWOP + FemalePWOP) ## effective number of parents for sample

RPP = (TrueP-ReducedSampP)^2
ReducedPP6[jj] = sqrt(mean(RPP))

ReducedS = dim(ReducedGenos)[[1]]
zzz = ReducedGenos[, ReducedSampP > 0.05 & ReducedSampP < 0.95, drop = F]
G = dim(zzz)[[2]]
MeanRsq2 = GetLD(zzz)
Rprime = MeanRsq2 - 1/(ReducedS-1)
if (ReducedS > 29) {ReducedNeR6[jj] = (2/3 + sqrt(4/9-7.2*Rprime))/(2*Rprime)} else {ReducedNeR6[jj] = (0.618 + sqrt(0.618^2-5.24*Rprime))/(2*Rprime)}
SS6[jj] = ReducedS

if (is.odd(jj)==TRUE) { ## then jj is odd, so skip Fst and record data for pop 1
    OldReducedSampP6 = ReducedSampP  } else {  ## else = even, so compute Fst
  combo2 = rbind(OldReducedSampP6,ReducedSampP)
##eliminate loci monomorphic in both samples
for (i in 1:NLoci)  {
QQQ[i] = combo2[1,i]+combo2[2,i]}
FP2 = combo2[,QQQ<2 & QQQ > 0]
FSTO6[jj/2] = GetFst(FP2) - 1/(2*ReducedS)
FPP6[jj/2] = (FSTP[jj/2]-FSTO6[jj/2])^2
} #end else

###############################################

print(paste0("Replicate = ",jj))
flush.console()

}  ## end for jj

HMNeR = 1/(mean(1/NeR,na.rm=TRUE))
HMPWOP = 1/(mean(1/PWOP,na.rm=TRUE))

HMNeRYank = 1/(mean(1/NeRYank,na.rm=TRUE))
HMPWOPYank = 1/(mean(1/PWOPYank,na.rm=TRUE))

HMReducedNeR1 = 1/(mean(1/ReducedNeR1,na.rm=TRUE))
HMPWOP1 = 1/(mean(1/PWOP1,na.rm=TRUE))

HMReducedNeR2 = 1/(mean(1/ReducedNeR2,na.rm=TRUE))
HMPWOP2 = 1/(mean(1/PWOP2,na.rm=TRUE))

HMReducedNeR3 = 1/(mean(1/ReducedNeR3,na.rm=TRUE))
HMPWOP3 = 1/(mean(1/PWOP3,na.rm=TRUE))

HMReducedNeR4 = 1/(mean(1/ReducedNeR4,na.rm=TRUE))
HMPWOP4 = 1/(mean(1/PWOP4,na.rm=TRUE))

HMReducedNeR5 = 1/(mean(1/ReducedNeR5,na.rm=TRUE))
HMPWOP5 = 1/(mean(1/PWOP5,na.rm=TRUE))

HMReducedNeR6 = 1/(mean(1/ReducedNeR6,na.rm=TRUE))
HMPWOP6 = 1/(mean(1/PWOP6,na.rm=TRUE))

TrueNe = Ne+0.5
if (Mating==3) {TrueNe = HMPWOP} ##PWOP tracks realized Ne for mixed mating model
DifNe = (1/NeR - 1/TrueNe)^2
RMSENe = sqrt(mean(DifNe,na.rm=TRUE))
DifNeYank = (1/NeRYank - 1/TrueNe)^2
RMSENeYank = sqrt(mean(DifNeYank,na.rm=TRUE))
DifNe1 = (1/ReducedNeR1 - 1/TrueNe)^2
RMSENe1 = sqrt(mean(DifNe1,na.rm=TRUE))
DifNe2 = (1/ReducedNeR2 - 1/TrueNe)^2
RMSENe2 = sqrt(mean(DifNe2,na.rm=TRUE))
DifNe3 = (1/ReducedNeR3 - 1/TrueNe)^2
RMSENe3 = sqrt(mean(DifNe3,na.rm=TRUE))
DifNe4 = (1/ReducedNeR4 - 1/TrueNe)^2
RMSENe4 = sqrt(mean(DifNe4,na.rm=TRUE))
DifNe5 = (1/ReducedNeR5 - 1/TrueNe)^2
RMSENe5 = sqrt(mean(DifNe5,na.rm=TRUE))
DifNe6 = (1/ReducedNeR6 - 1/TrueNe)^2
RMSENe6 = sqrt(mean(DifNe6,na.rm=TRUE))

HMS = 1/(mean(1/SS0))
HMSYank = 1/(mean(1/SSYank))
HMS1 = 1/(mean(1/SS1))
HMS2 = 1/(mean(1/SS2))
HMS3 = 1/(mean(1/SS3))
HMS4 = 1/(mean(1/SS4))
HMS5 = 1/(mean(1/SS5))
HMS6 = 1/(mean(1/SS6))

SibNe = 4/(mean(HS) + 2*mean(FS))

AllPWOP = c(HMPWOP, HMPWOPYank, HMPWOP1, HMPWOP2, HMPWOP3, HMPWOP4, HMPWOP5, HMPWOP6)
AllS = c(HMS, HMSYank, HMS1, HMS2, HMS3, HMS4, HMS5, HMS6)
AllNe2 = c(RMSENe/RMSENe, RMSENeYank/RMSENe, RMSENe1/RMSENe, RMSENe2/RMSENe, RMSENe3/RMSENe, RMSENe4/RMSENe, RMSENe5/RMSENe, RMSENe6/RMSENe)
AllNe = c(HMNeR, HMNeRYank, HMReducedNeR1,HMReducedNeR2,HMReducedNeR3,HMReducedNeR4,HMReducedNeR5,HMReducedNeR6)
AllPP = c(1, mean(ReducedPPYank)/mean(SamplePP), mean(ReducedPP1)/mean(SamplePP), mean(ReducedPP2)/mean(SamplePP),mean(ReducedPP3)/mean(SamplePP),mean(ReducedPP4)/mean(SamplePP),mean(ReducedPP5)/mean(SamplePP),mean(ReducedPP6)/mean(SamplePP))
AllFPP = c(1,sqrt(mean(FPPYank))/sqrt(mean(FPP)),sqrt(mean(FPP1))/sqrt(mean(FPP)),sqrt(mean(FPP2))/sqrt(mean(FPP)),sqrt(mean(FPP3))/sqrt(mean(FPP)),sqrt(mean(FPP4))/sqrt(mean(FPP)),sqrt(mean(FPP5))/sqrt(mean(FPP)),sqrt(mean(FPP6))/sqrt(mean(FPP)))
NeRatio = AllNe/TrueNe
All = cbind(AllS,AllPP,AllFPP,AllPWOP,AllNe,NeRatio,AllNe2)
rownames(All) = c("Full","Yank2","25","50","75","90","95","100")
colnames(All) = c("Mean S","RRMSE P","RRMSE Fst","PWOP","HM Ne^","Ne^/Ne","RRMSE 1/Ne")
All
mean(FS)
mean(HS)
mean(MaxSib)
SibNe
