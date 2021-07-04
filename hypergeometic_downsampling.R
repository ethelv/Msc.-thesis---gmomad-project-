setwd("~/Documents/Research/all chrom")

# Read all chrom singelton data ----
df <- read.csv("all_chr_singelton.csv")

# Choose only canonical transcripts ----
stringsAsFactors= FALSE
canonical_transcript_list<- read.csv("canonical transcript list.csv", header = F)
colnames(canonical_transcript_list)<- c("Chrom", "Start", "End", "Number", "Transcript", "Gene" )

m<- unlist(strsplit(as.character(canonical_transcript_list$Transcript),".",fixed = TRUE))
canonical_transcript_list$Transcript<- m[seq(1, length(m), 2)]

### Choose only autosomal chromosomes ---- 
df_aut<- df[which(!(df$CHROM %in% c("X","Y"))),]

#take only synonymous variants ----
df_aut_syn<- df_aut[which(df_aut$Consequence=="synonymous_variant"),]
df_aut_syn<- df_aut_syn[which(df_aut_syn$Gene %in% canonical_transcript_list$Transcript),]

sum(df_aut_syn$sum_AC_female)
#277507
sum(df_aut_syn$sum_AC_male)
#349322
sum(df_aut_syn$sum_AC)
#626829

# Down-sampling syn

y <- c(10000,25000,50000,75000,100000,150000,200000,300000,400000,500000,600000)
x <- y/ sum(df_aut_syn$sum_AC)
r_sq_syn<- c()
for (val in x) {
  #print(val)
  syn_new_AC<- trunc(df_aut_syn$sum_AC* val)
  syn_new_AC_female<- c()
  for (i in c(1:length(syn_new_AC))){
    syn_new_AC_female[i]<- rhyper(1,df_aut_syn$sum_AC_female[i], df_aut_syn$sum_AC_male[i],syn_new_AC[i])
  }
  syn_new_AC_male  = syn_new_AC - syn_new_AC_female
  r_sq_syn[match(val,x)]<- summary(lm(syn_new_AC_male ~ syn_new_AC_female))$r.squared
}
plot(y, r_sq_syn, main = "synonymous autosome", xlab = "sample size", ylab = "r.squared")

#take only missense variants ----
df_aut_mis<- df_aut[which(df_aut$Consequence=="missense_variant"),]
df_aut_mis<- df_aut_mis[which(df_aut_mis$Gene %in% canonical_transcript_list$Transcript),]

sum(df_aut_mis$sum_AC_female)
#629022
sum(df_aut_mis$sum_AC_male)
#790638
sum(df_aut_mis$sum_AC)
#1419660

y <-  c(10000,25000,50000,75000,100000,150000,200000,300000,400000,500000,600000)
x <- y/ sum(df_aut_mis$sum_AC)
r_sq_mis<- c()
for (val in x) {
  print(val)
  miss_new_AC<- trunc(df_aut_mis$sum_AC* val)
  miss_new_AC_female<- c()
  for (i in c(1:length(miss_new_AC))){
    miss_new_AC_female[i]<- rhyper(1,df_aut_mis$sum_AC_female[i], df_aut_mis$sum_AC_male[i],miss_new_AC[i])
  }
  miss_new_AC_male  = miss_new_AC - miss_new_AC_female
  r_sq_mis[match(val,x)]<- summary(lm(miss_new_AC_male ~ miss_new_AC_female))$r.squared
}
plot(y, r_sq_mis, main = "missense autosomes", xlab = "sample size", ylab = "r.squared")


plot(y, r_sq_mis-r_sq_syn, main = "r.squared missense - r.squared synonymous", ylim = c(-0.06, 0.09),
     ylab= "r-squared difference", xlab = "sample size")
abline(0,0)


#take only lof variants ----
df_aut_lof<- df_aut[which(df_aut$Consequence%in%
                                         c("stop_gained","splice_donor_variant",
                                           "splice_acceptor_variant","frameshift_variant")),]
df_aut_lof<- df_aut_lof[which(df_aut_lof$Gene %in% canonical_transcript_list$Transcript),]
library(plyr)
df_aut_lof<- ddply(df_aut_lof,"Gene",numcolwise(sum))

sum(df_aut_lof$sum_AC_female)
#68792
sum(df_aut_lof$sum_AC_male)
#85409
sum(df_aut_lof$sum_AC)
#154201

y_lof <-  c(5000, 10000,25000,50000,65000,75000,90000,100000,125000, 150000)
x <- y_lof/ sum(df_aut_lof$sum_AC)
r_sq_lof<- c()
for (val in x) {
  print(val)
  lof_new_AC<- trunc(df_aut_lof$sum_AC* val)
  lof_new_AC_female<- c()
  for (i in c(1:length(lof_new_AC))){
    lof_new_AC_female[i]<- rhyper(1,df_aut_lof$sum_AC_female[i], df_aut_lof$sum_AC_male[i],lof_new_AC[i])
  }
  lof_new_AC_male  = lof_new_AC - lof_new_AC_female
  r_sq_lof[match(val,x)]<- summary(lm(lof_new_AC_male ~ lof_new_AC_female))$r.squared
}

plot(y_lof, r_sq_lof, main = "LoF singeltones", xlab = "sample size", ylab = "r.squared")

#plot all types of variants----
plot( y, r_sq_syn,ylim = c(0, 1), xlab = "sample size", ylab = expression(paste("R"^"2")),
      type="b", col="red" ,main = bquote(R^2~ "Between the Sexes Autosomes"))
lines(y, r_sq_mis, type="b", col="green")
lines( y_lof, r_sq_lof, type="b", col="blue")
legend("bottomright",  legend = c("Synonymous", "Missense", "LoF"), lty=c(1,1),col = c("red", "green", "blue"))

# X chromosome ---- 
df_X<- df[which(df$CHROM == "X"),]

#take only synonymous variants ----
df_X_syn<- df_X[which(df_X$Consequence=="synonymous_variant"),]
df_X_syn<- df_X_syn[which(df_X_syn$Gene %in% canonical_transcript_list$Transcript),]

sum(df_X_syn$sum_AC_female)
#10206
sum(df_X_syn$sum_AC_male)
#6178
sum(df_X_syn$sum_AC)
#16384


# downsampling syn

y_x_syn <- c(1000,2000,2500,3500,5000,6500,8000,9500,11000,13000,14500,16000)
x <- y_x_syn/ sum(df_X_syn$sum_AC)
r_sq_syn_x<- c()
for (val in x) {
  print(val)
  syn_new_AC_x<- trunc(df_X_syn$sum_AC* val)
  syn_new_AC_female_x<- c()
  for (i in c(1:length(syn_new_AC_x))){
    syn_new_AC_female_x[i]<- rhyper(1,df_X_syn$sum_AC_female[i], df_X_syn$sum_AC_male[i],syn_new_AC_x[i])
  }
  syn_new_AC_male_x  = syn_new_AC_x - syn_new_AC_female_x
  r_sq_syn_x[match(val,x)]<- summary(lm(syn_new_AC_male_x ~ syn_new_AC_female_x))$r.squared
}
plot(y_x_syn, r_sq_syn_x, main = "synonymous singeltones", xlab = "sample size", ylab = "r.squared")

#take only missense variants ----
df_x_mis<- df_X[which(df_X$Consequence=="missense_variant"),]
df_x_mis<- df_x_mis[which(df_x_mis$Gene %in% canonical_transcript_list$Transcript),]

sum(df_x_mis$sum_AC_female)
#21675
sum(df_x_mis$sum_AC_male)
#13022
sum(df_x_mis$sum_AC)
#34697



# downsampling miss

y_x_mis <- c(1000,2000,2500,3500,5000,6500,8000,9500,11000,13000,14500,16000)
x <- y_x_mis/ sum(df_x_mis$sum_AC)
r_sq_mis_x<- c()
for (val in x) {
  print(val)
  miss_new_AC_x<- trunc(df_x_mis$sum_AC* val)
  miss_new_AC_female_x<- c()
  for (i in c(1:length(miss_new_AC_x))){
    miss_new_AC_female_x[i]<- rhyper(1,df_x_mis$sum_AC_female[i], df_x_mis$sum_AC_male[i],miss_new_AC_x[i])
  }
  miss_new_AC_male_x  = miss_new_AC_x - miss_new_AC_female_x
  r_sq_mis_x[match(val,x)]<- summary(lm(miss_new_AC_male_x ~ miss_new_AC_female_x))$r.squared
}
plot(y_x_mis, r_sq_mis_x, main = "missense singeltones", xlab = "sample size", ylab = "r.squared")


#take only lof variants ----
df_x_lof<- df_X[which(df_X$Consequence %in%
                                     c("stop_gained","splice_donor_variant",
                                       "splice_acceptor_variant","frameshift_variant")),]
df_x_lof<- df_x_lof[which(df_x_lof$Gene %in% canonical_transcript_list$Transcript),]
library(plyr)
df_x_lof<- ddply(df_x_lof,"Gene",numcolwise(sum))

sum(df_x_lof$sum_AC_female)
#1347
sum(df_x_lof$sum_AC_male)
#798
sum(df_x_lof$sum_AC)
#2145

y_lof_x <-  c(100,500,1000,1250,1500,1800,2140)
x <- y_lof_x/ sum(df_x_lof$sum_AC)
r_sq_lof_x<- c()
for (val in x) {
  print(val)
  lof_new_AC_x<- trunc(df_x_lof$sum_AC* val)
  lof_new_AC_female_x<- c()
  for (i in c(1:length(lof_new_AC))){
    lof_new_AC_female_x[i]<- rhyper(1,df_x_lof$sum_AC_female[i], df_x_lof$sum_AC_male[i],lof_new_AC_x[i])
  }
  lof_new_AC_male_x  = lof_new_AC_x - lof_new_AC_female_x
  r_sq_lof_x[match(val,x)]<- summary(lm(lof_new_AC_male_x ~ lof_new_AC_female_x))$r.squared
}
plot(y_lof_x, r_sq_lof_x, main = "LoF singeltones", xlab = "sample size", ylab = "r.squared")

#plot all types of variants----
plot( y_x_syn, r_sq_syn_x,ylim = c(0, 1), xlab = "sample size", ylab = expression(paste("R"^"2")),
      type="b", col="red" , main = bquote(R^2~ "Between the Sexes X chromosome"))
lines(y_x_syn, r_sq_mis_x, type="b", col="green")
lines( y_lof_x, r_sq_lof_x, type="b", col="blue")
legend("bottomright",  legend = c("Synonymous", "Missense", "LoF"), lty=c(1,1),
       col = c("red", "green", "blue"))

