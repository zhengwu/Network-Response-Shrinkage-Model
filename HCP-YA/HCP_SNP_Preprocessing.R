## the genetic pre-process part

./plink --noweb --bfile /home/ysm/zhao_yize/yz875/brainnetwork/hcp_plink_format/HCP_GeneticData/matrix/MEGA_Chip --extract /home/ysm/zhao_yize/yz875/brainnetwork/hcp_plink_format/HCP_GeneticData/matrix/SNPname.txt --maf 0.05 --mind 0.1 --geno 0.1 --hwe 0.000001 --make-founders --indep-pairwise 100 10 0.2 --out prune

./plink --noweb --bfile /home/ysm/zhao_yize/yz875/brainnetwork/hcp_plink_format/HCP_GeneticData/matrix/MEGA_Chip  --maf 0.05 --mind 0.1 --geno 0.1 --hwe 0.000001 --extract prune.prune.in --recodeA --out fSNP

./plink --noweb --bfile /home/ysm/zhao_yize/yz875/brainnetwork/hcp_plink_format/HCP_GeneticData/matrix/MEGA_Chip  --maf 0.05 --mind 0.1 --geno 0.1 --hwe 0.000001 --extract prune.prune.in --pca 10

## find SNP-set

find.blocks=function(geno.matrix.scaled,r2=.05,look.ahead=100,min.prop=0.5){
 m=dim(geno.matrix.scaled)[2]
 n=dim(geno.matrix.scaled)[1]
 cluster.begin=1
 cluster.end=NULL

 last=0

 while (last<(m-look.ahead)){
 begin=last+1
 end=begin+look.ahead-1
 prop.vector=rep(0,look.ahead)
 x=cor(geno.matrix.scaled[,begin:end])^2 > r2
 x[is.na(x)]=0
 for (i in (1:dim(x)[1])){
  prop.vector[i]=mean(x[1:i,1:i])
 }
 plot(prop.vector)
 diff=c(0,prop.vector[2:look.ahead]-prop.vector[1:(look.ahead-1)])
 ind=max(2,which.max(cumsum(((diff<=0)&(prop.vector>min.prop)))))
 last=last+ind-1
 cluster.end=c(cluster.end,last)
 cluster.begin=c(cluster.begin,last+1)
 #print(last)
 #print(cluster.begin)
 #print(cluster.end)
 }
 cluster.end=c(cluster.end,m)
 return(list(cluster.begin=cluster.begin,cluster.end=cluster.end))
 }

scaleSNP=scale(fSNP)
blocksize=find.blocks(scaleSNP,r2=.02,look.ahead=1000,min.prop=0.5)
save(blocksize,file='Bsize_new.rdata')