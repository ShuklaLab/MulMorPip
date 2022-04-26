rm(list=ls())
df=ukbtools::ukb_df("ukb46777")
library("icd")
dia=grep("f41270_0", names(df), value = T)
v=df[c("eid",dia)]
#write.table(v, file="all.tsv", row.names=F, quote=F)

#v=read.table("ra.tsv", header=T)

v_cmb=comorbid_charlson(v,return_df=T)
#v_cmb=comorbid_ahrq(v,return_df=T)
#v_cmb=comorbid_elix(v,return_df=T)
#v_cmb=comorbid_quan_deyo(v,return_df=T)
#v_cmb=comorbid_quan_elix(v,return_df=T)
#v_cmb=comorbid_pccc_pcs(v,return_df=T)
#v_cmb=comorbid_pccc_dx(v,return_df=T)
#v_cmb=comorbid_ccs(v,return_df=T)
#v_cmb=comorbid_hcc(v,return_df=T)

######## Data process
data_cmb1=v_cmb[-c(1)]
row.names(data_cmb1)=v_cmb$eid
data_cmb1$nod=rowSums(data_cmb1)
data_cmb=data_cmb1[data_cmb1$nod>1,]
lst=read.table("w48433_20210809.csv", header=F)
data_cmb=data_cmb[!(row.names(data_cmb) %in% lst$V1),]
table(data_cmb$nod)
data_cmb$nod<-NULL

###Desity plot wrt a
#plot(res.mca$ind$coord[,1],res.mca$ind$coord[,2], col=rgb(0,0,0,(a+15)/1000))

library("FactoMineR")
res.mca=MCA(data_cmb, graph=F)
#print(res.mca)

plot(res.mca$ind$coord[,1],res.mca$ind$coord[,2])

library("factoextra")
eig.val <- get_eigenvalue(res.mca)
fviz_screeplot(res.mca, addlabels = T, ylim = c(0, 45))
fviz_mca_biplot(res.mca, 
               repel = TRUE, # Avoid text overlapping (slow if many point)
               ggtheme = theme_minimal())

## PCA plot
pc=data.frame(res.mca$ind$coord[,1:2])
colnames(pc)=c("pc1","pc2")
pc$map=pc$pc1+1.8*pc$pc2
#hist(pc$map)
pc$map1=1.8*pc$pc1-pc$pc2
plot(pc$map,pc$map1)
abline(v=0)
abline(v=.7)
abline(v=1.4)
abline(v=2.1)
pc$clust=1
pc$clust[pc$map>0]=2
pc$clust[pc$map>.7]=3
pc$clust[pc$map>1.4]=4
pc$clust[pc$map>2.1]=5
plot(pc$pc1,pc$pc2, col=pc$clust)

table(pc$clust)

##### Validation
library(caret)
set.seed(200)
idx=createDataPartition(pc$clust,p=.8,list=F)
nrow(pc[-idx,])
table(pc[-idx,]$clust)
library(rpart)
library(rpart.plot)
fit=rpart(pc[idx,]$clust~.,data=data_cmb[idx,], method="class")
rpart.plot(fit,extra=104)
y_test=predict(fit,data_cmb[-idx,],type="class")
length(y_test)
table(y_test)
cm=table(pc[-idx,]$clust,y_test)
cm
acc=(cm[1,1]+cm[2,2]+cm[3,3]+cm[4,4]+cm[5,5])/length(y_test)
acc
#library(jaccard)
#jaccard(pc[-idx,]$clust,y_test)
#jaccard.test(pc[-idx,]$clust,y_test)

# check for claster 3 missclass
a=data_cmb[-idx,]
a$clust=pc[-idx,]$clust
a$y_test=y_test
b=a[(a$clust==3)&(a$y_test==2),]
table(b$Dementia)

## MCA plot
mca.test=MCA(data_cmb[-idx,],graph=F)
plot(mca.test$ind$coord[,1],mca.test$ind$coord[,2],col=y_test)
plot(mca.test$ind$coord[,1],mca.test$ind$coord[,2],col=pc[-idx,]$clust)


node=data.frame()
for (i in 1:16){
  for (j in (i+1):17){
    res=chisq.test(table(data_cmb[[colnames(data_cmb)[i]]], data_cmb[[colnames(data_cmb)[j]]]),correct=F)
    print(paste(colnames(data_cmb)[i],colnames(data_cmb)[j],res$p.value))
    row=c(i,j,res$p.value)
    node=rbind(node,row)
}}

data_cmb1=data_cmb[pc$clust==1,]
node=data.frame()
for (i in 1:16){
  for (j in (i+1):17){
    tab=table(data_cmb1[[colnames(data_cmb1)[i]]], data_cmb1[[colnames(data_cmb1)[j]]])
    if(length(tab)==4){
      res=chisq.test(tab,correct=F)
      print(paste(colnames(data_cmb1)[i],colnames(data_cmb1)[j],res$p.value))
      row=c(i,j,res$p.value)
      node=rbind(node,row)
}}}
write.table(node, file="node1.tsv", row.names=F, quote=F)

data_cmb1=data_cmb[pc$clust==2,]
node=data.frame()
for (i in 1:16){
  for (j in (i+1):17){
    tab=table(data_cmb1[[colnames(data_cmb1)[i]]], data_cmb1[[colnames(data_cmb1)[j]]])
    if(length(tab)==4){
      res=chisq.test(tab,correct=F)
      print(paste(colnames(data_cmb1)[i],colnames(data_cmb1)[j],res$p.value))
      row=c(i,j,res$p.value)
      node=rbind(node,row)
}}}
write.table(node, file="node2.tsv", row.names=F, quote=F)

data_cmb1=data_cmb[pc$clust==3,]
node=data.frame()
for (i in 1:16){
  for (j in (i+1):17){
    tab=table(data_cmb1[[colnames(data_cmb1)[i]]], data_cmb1[[colnames(data_cmb1)[j]]])
    if(length(tab)==4){
      res=chisq.test(tab,correct=F)
      print(paste(colnames(data_cmb1)[i],colnames(data_cmb1)[j],res$p.value))
      row=c(i,j,res$p.value)
      node=rbind(node,row)
}}}
write.table(node, file="node3.tsv", row.names=F, quote=F)


data_cmb1=data_cmb[pc$clust==4,]
node=data.frame()
for (i in 1:16){
  for (j in (i+1):17){
    tab=table(data_cmb1[[colnames(data_cmb1)[i]]], data_cmb1[[colnames(data_cmb1)[j]]])
    if(length(tab)==4){
      res=chisq.test(tab,correct=F)
      print(paste(colnames(data_cmb1)[i],colnames(data_cmb1)[j],res$p.value))
      row=c(i,j,res$p.value)
      node=rbind(node,row)
}}}
write.table(node, file="node4.tsv", row.names=F, quote=F)

data_cmb1=data_cmb[pc$clust==5,]
node=data.frame()
for (i in 1:16){
  for (j in (i+1):17){
    tab=table(data_cmb1[[colnames(data_cmb1)[i]]], data_cmb1[[colnames(data_cmb1)[j]]])
    if(length(tab)==4){
      res=chisq.test(tab,correct=F)
      print(paste(colnames(data_cmb1)[i],colnames(data_cmb1)[j],res$p.value))
      row=c(i,j,res$p.value)
      node=rbind(node,row)
}}}
write.table(node, file="node5.tsv", row.names=F, quote=F)

data_cmb1=data_cmb[data_cmb$MI==T,]
node=data.frame()
for (i in 1:16){
  for (j in (i+1):17){
    tab=table(data_cmb1[[colnames(data_cmb1)[i]]], data_cmb1[[colnames(data_cmb1)[j]]])
    if(length(tab)==4){
      res=chisq.test(tab,correct=F)
      print(paste(colnames(data_cmb1)[i],colnames(data_cmb1)[j],res$p.value))
      row=c(i,j,res$p.value)
      node=rbind(node,row)
}}}
write.table(node, file="mi.tsv", row.names=F, quote=F)

des=MCA(data_cmb1,graph=F)
fviz_mca_biplot(des)

data_cmb1=data_cmb[data_cmb$CHF==T,]
node=data.frame()
for (i in 1:16){
  for (j in (i+1):17){
    tab=table(data_cmb1[[colnames(data_cmb1)[i]]], data_cmb1[[colnames(data_cmb1)[j]]])
    if(length(tab)==4){
      res=chisq.test(tab,correct=F)
      print(paste(colnames(data_cmb1)[i],colnames(data_cmb1)[j],res$p.value))
      row=c(i,j,res$p.value)
      node=rbind(node,row)
}}}
write.table(node, file="chf.tsv", row.names=F, quote=F)

des=MCA(data_cmb1,graph=F)
fviz_mca_biplot(des)

data_cmb1=data_cmb[data_cmb$PVD==T,]
node=data.frame()
for (i in 1:16){
  for (j in (i+1):17){
    tab=table(data_cmb1[[colnames(data_cmb1)[i]]], data_cmb1[[colnames(data_cmb1)[j]]])
    if(length(tab)==4){
      res=chisq.test(tab,correct=F)
      print(paste(colnames(data_cmb1)[i],colnames(data_cmb1)[j],res$p.value))
      row=c(i,j,res$p.value)
      node=rbind(node,row)
}}}
write.table(node, file="pvd.tsv", row.names=F, quote=F)

des=MCA(data_cmb1,graph=F)
fviz_mca_biplot(des)


data_cmb1=data_cmb[data_cmb$Stroke==T,]
node=data.frame()
for (i in 1:16){
  for (j in (i+1):17){
    tab=table(data_cmb1[[colnames(data_cmb1)[i]]], data_cmb1[[colnames(data_cmb1)[j]]])
    if(length(tab)==4){
      res=chisq.test(tab,correct=F)
      print(paste(colnames(data_cmb1)[i],colnames(data_cmb1)[j],res$p.value))
      row=c(i,j,res$p.value)
      node=rbind(node,row)
}}}
write.table(node, file="str.tsv", row.names=F, quote=F)

des=MCA(data_cmb1,graph=F)
fviz_mca_biplot(des)

data_cmb1=data_cmb[data_cmb$Dementia==T,]
node=data.frame()
for (i in 1:16){
  for (j in (i+1):17){
    tab=table(data_cmb1[[colnames(data_cmb1)[i]]], data_cmb1[[colnames(data_cmb1)[j]]])
    if(length(tab)==4){
      res=chisq.test(tab,correct=F)
      print(paste(colnames(data_cmb1)[i],colnames(data_cmb1)[j],res$p.value))
      row=c(i,j,res$p.value)
      node=rbind(node,row)
}}}
write.table(node, file="dem.tsv", row.names=F, quote=F)

des=MCA(data_cmb1,graph=F)
fviz_mca_biplot(des)


data_cmb1=data_cmb[data_cmb$Pulmonary==T,]
node=data.frame()
for (i in 1:16){
  for (j in (i+1):17){
    tab=table(data_cmb1[[colnames(data_cmb1)[i]]], data_cmb1[[colnames(data_cmb1)[j]]])
    if(length(tab)==4){
      res=chisq.test(tab,correct=F)
      print(paste(colnames(data_cmb1)[i],colnames(data_cmb1)[j],res$p.value))
      row=c(i,j,res$p.value)
      node=rbind(node,row)
}}}
write.table(node, file="pul.tsv", row.names=F, quote=F)

des=MCA(data_cmb1,graph=F)
fviz_mca_biplot(des)

data_cmb1=data_cmb[data_cmb$Rheumatic==T,]
node=data.frame()
for (i in 1:16){
  for (j in (i+1):17){
    tab=table(data_cmb1[[colnames(data_cmb1)[i]]], data_cmb1[[colnames(data_cmb1)[j]]])
    if(length(tab)==4){
      res=chisq.test(tab,correct=F)
      print(paste(colnames(data_cmb1)[i],colnames(data_cmb1)[j],res$p.value))
      row=c(i,j,res$p.value)
      node=rbind(node,row)
}}}
write.table(node, file="rhe.tsv", row.names=F, quote=F)

des=MCA(data_cmb1,graph=F)
fviz_mca_biplot(des)

data_cmb1=data_cmb[data_cmb$PUD==T,]
node=data.frame()
for (i in 1:16){
  for (j in (i+1):17){
    tab=table(data_cmb1[[colnames(data_cmb1)[i]]], data_cmb1[[colnames(data_cmb1)[j]]])
    if(length(tab)==4){
      res=chisq.test(tab,correct=F)
      print(paste(colnames(data_cmb1)[i],colnames(data_cmb1)[j],res$p.value))
      row=c(i,j,res$p.value)
      node=rbind(node,row)
}}}
write.table(node, file="pud.tsv", row.names=F, quote=F)

des=MCA(data_cmb1,graph=F)
fviz_mca_biplot(des)

data_cmb1=data_cmb[data_cmb$LiverMild==T,]
node=data.frame()
for (i in 1:16){
  for (j in (i+1):17){
    tab=table(data_cmb1[[colnames(data_cmb1)[i]]], data_cmb1[[colnames(data_cmb1)[j]]])
    if(length(tab)==4){
      res=chisq.test(tab,correct=F)
      print(paste(colnames(data_cmb1)[i],colnames(data_cmb1)[j],res$p.value))
      row=c(i,j,res$p.value)
      node=rbind(node,row)
}}}
write.table(node, file="lm.tsv", row.names=F, quote=F)

des=MCA(data_cmb1,graph=F)
fviz_mca_biplot(des)


data_cmb1=data_cmb[data_cmb$DM==T,]
node=data.frame()
for (i in 1:16){
  for (j in (i+1):17){
    tab=table(data_cmb1[[colnames(data_cmb1)[i]]], data_cmb1[[colnames(data_cmb1)[j]]])
    if(length(tab)==4){
      res=chisq.test(tab,correct=F)
      print(paste(colnames(data_cmb1)[i],colnames(data_cmb1)[j],res$p.value))
      row=c(i,j,res$p.value)
      node=rbind(node,row)
}}}
write.table(node, file="dm.tsv", row.names=F, quote=F)

des=MCA(data_cmb1,graph=F)
fviz_mca_biplot(des)


data_cmb1=data_cmb[data_cmb$DMcx==T,]
node=data.frame()
for (i in 1:16){
  for (j in (i+1):17){
    tab=table(data_cmb1[[colnames(data_cmb1)[i]]], data_cmb1[[colnames(data_cmb1)[j]]])
    if(length(tab)==4){
      res=chisq.test(tab,correct=F)
      print(paste(colnames(data_cmb1)[i],colnames(data_cmb1)[j],res$p.value))
      row=c(i,j,res$p.value)
      node=rbind(node,row)
}}}
write.table(node, file="dmc.tsv", row.names=F, quote=F)

des=MCA(data_cmb1,graph=F)
fviz_mca_biplot(des)


data_cmb1=data_cmb[data_cmb$Paralysis==T,]
node=data.frame()
for (i in 1:16){
  for (j in (i+1):17){
    tab=table(data_cmb1[[colnames(data_cmb1)[i]]], data_cmb1[[colnames(data_cmb1)[j]]])
    if(length(tab)==4){
      res=chisq.test(tab,correct=F)
      print(paste(colnames(data_cmb1)[i],colnames(data_cmb1)[j],res$p.value))
      row=c(i,j,res$p.value)
      node=rbind(node,row)
}}}
write.table(node, file="par.tsv", row.names=F, quote=F)

des=MCA(data_cmb1,graph=F)
fviz_mca_biplot(des)


data_cmb1=data_cmb[data_cmb$Renal==T,]
node=data.frame()
for (i in 1:16){
  for (j in (i+1):17){
    tab=table(data_cmb1[[colnames(data_cmb1)[i]]], data_cmb1[[colnames(data_cmb1)[j]]])
    if(length(tab)==4){
      res=chisq.test(tab,correct=F)
      print(paste(colnames(data_cmb1)[i],colnames(data_cmb1)[j],res$p.value))
      row=c(i,j,res$p.value)
      node=rbind(node,row)
}}}
write.table(node, file="ren.tsv", row.names=F, quote=F)

des=MCA(data_cmb1,graph=F)
fviz_mca_biplot(des)


data_cmb1=data_cmb[data_cmb$Cancer==T,]
node=data.frame()
for (i in 1:16){
  for (j in (i+1):17){
    tab=table(data_cmb1[[colnames(data_cmb1)[i]]], data_cmb1[[colnames(data_cmb1)[j]]])
    if(length(tab)==4){
      res=chisq.test(tab,correct=F)
      print(paste(colnames(data_cmb1)[i],colnames(data_cmb1)[j],res$p.value))
      row=c(i,j,res$p.value)
      node=rbind(node,row)
}}}
write.table(node, file="can.tsv", row.names=F, quote=F)

des=MCA(data_cmb1,graph=F)
fviz_mca_biplot(des)


data_cmb1=data_cmb[data_cmb$LiverSevere==T,]
node=data.frame()
for (i in 1:16){
  for (j in (i+1):17){
    tab=table(data_cmb1[[colnames(data_cmb1)[i]]], data_cmb1[[colnames(data_cmb1)[j]]])
    if(length(tab)==4){
      res=chisq.test(tab,correct=F)
      print(paste(colnames(data_cmb1)[i],colnames(data_cmb1)[j],res$p.value))
      row=c(i,j,res$p.value)
      node=rbind(node,row)
}}}
write.table(node, file="ls.tsv", row.names=F, quote=F)

des=MCA(data_cmb1,graph=F)
fviz_mca_biplot(des)


data_cmb1=data_cmb[data_cmb$Mets==T,]
node=data.frame()
for (i in 1:16){
  for (j in (i+1):17){
    tab=table(data_cmb1[[colnames(data_cmb1)[i]]], data_cmb1[[colnames(data_cmb1)[j]]])
    if(length(tab)==4){
      res=chisq.test(tab,correct=F)
      print(paste(colnames(data_cmb1)[i],colnames(data_cmb1)[j],res$p.value))
      row=c(i,j,res$p.value)
      node=rbind(node,row)
}}}
write.table(node, file="met.tsv", row.names=F, quote=F)

des=MCA(data_cmb1,graph=F)
fviz_mca_biplot(des)


data_cmb1=data_cmb[data_cmb$HIV==T,]
node=data.frame()
for (i in 1:16){
  for (j in (i+1):17){
    tab=table(data_cmb1[[colnames(data_cmb1)[i]]], data_cmb1[[colnames(data_cmb1)[j]]])
    if(length(tab)==4){
      res=chisq.test(tab,correct=F)
      print(paste(colnames(data_cmb1)[i],colnames(data_cmb1)[j],res$p.value))
      row=c(i,j,res$p.value)
      node=rbind(node,row)
}}}
write.table(node, file="hiv.tsv", row.names=F, quote=F)

des=MCA(data_cmb1,graph=F)
fviz_mca_biplot(des)

####### For comorbidity matrix plot
library(plot.matrix)
tt=matrix(0,17,17)
pp=matrix(0,17,17)
rownames(tt)=colnames(data_cmb)
colnames(tt)=colnames(data_cmb)
rownames(pp)=colnames(data_cmb)
colnames(pp)=colnames(data_cmb)
data_cmb1=data_cmb[pc$clust==5,]
for (i in 1:17){
  for (j in 1:17){
    res1=table(data_cmb1[[colnames(data_cmb1)[i]]], data_cmb1[[colnames(data_cmb1)[j]]])
    a=try(res1["TRUE","TRUE"],silent=T)
    if(class(a)=="integer"){
       tt[j,i]=a
       pp[j,i]=a/nrow(data_cmb1)
}}}
#par(mfrow=c(2,1),mar=c(5.1,4.1,4.1,4.1))
#plot(tt, digit=0, gray=T)
png("c.png", width=1400, height=1400)
plot(pp, col=rev(gray.colors(20)))
dev.off()
############### End of plot


para=v[data_cmb$Paralysis==T & pc$clust==1,]
write.table(para, file="par1.tsv", row.names=F, quote=F)
#icd10_map_charlson[12]
para=v[data_cmb$Stroke==T & pc$clust==1,]
write.table(para, file="str1.tsv", row.names=F, quote=F)
#icd10_map_charlson[4]
para=v[data_cmb$Dementia==T & pc$clust==1,]
write.table(para, file="dem1.tsv", row.names=F, quote=F)
#icd10_map_charlson[5]

para=v[data_cmb$Paralysis==T & pc$clust==2,]
write.table(para, file="par2.tsv", row.names=F, quote=F)
para=v[data_cmb$Stroke==T & pc$clust==2,]
write.table(para, file="str2.tsv", row.names=F, quote=F)
para=v[data_cmb$Dementia==T & pc$clust==2,]
write.table(para, file="dem2.tsv", row.names=F, quote=F)

para=v[data_cmb$Paralysis==T & pc$clust==3,]
write.table(para, file="par3.tsv", row.names=F, quote=F)
para=v[data_cmb$Stroke==T & pc$clust==3,]
write.table(para, file="str3.tsv", row.names=F, quote=F)
para=v[data_cmb$Dementia==T & pc$clust==3,]
write.table(para, file="dem3.tsv", row.names=F, quote=F)

para=v[data_cmb$Paralysis==T & pc$clust==4,]
write.table(para, file="par4.tsv", row.names=F, quote=F)
para=v[data_cmb$Stroke==T & pc$clust==4,]
write.table(para, file="str4.tsv", row.names=F, quote=F)
para=v[data_cmb$Dementia==T & pc$clust==4,]
write.table(para, file="dem4.tsv", row.names=F, quote=F)

para=v[data_cmb$Paralysis==T & pc$clust==5,]
write.table(para, file="par5.tsv", row.names=F, quote=F)
para=v[data_cmb$Stroke==T & pc$clust==5,]
write.table(para, file="str5.tsv", row.names=F, quote=F)
para=v[data_cmb$Dementia==T & pc$clust==5,]
write.table(para, file="dem5.tsv", row.names=F, quote=F)


write.table(node, file="node5.tsv", row.names=F, quote=F)

write.table(node, file="node.tsv", row.names=F, quote=F)
write.table(node, file="node_ra.tsv", row.names=F, quote=F)


###########Barplot
par(mfrow=c(5,1))
barplot(100*colSums(data.matrix(data_cmb[pc$clust==1,]))/nrow(data_cmb[pc$clust==1,]), ylim=c(0,100))
barplot(100*colSums(data.matrix(data_cmb[pc$clust==2,]))/nrow(data_cmb[pc$clust==2,]), ylim=c(0,100))
barplot(100*colSums(data.matrix(data_cmb[pc$clust==3,]))/nrow(data_cmb[pc$clust==3,]), ylim=c(0,100))
barplot(100*colSums(data.matrix(data_cmb[pc$clust==4,]))/nrow(data_cmb[pc$clust==4,]), ylim=c(0,100))
barplot(100*colSums(data.matrix(data_cmb[pc$clust==5,]))/nrow(data_cmb[pc$clust==5,]), ylim=c(0,100))

c1=100*colSums(data.matrix(data_cmb[pc$clust==1,]))/nrow(data_cmb[pc$clust==1,])
c2=100*colSums(data.matrix(data_cmb[pc$clust==2,]))/nrow(data_cmb[pc$clust==2,])
c3=100*colSums(data.matrix(data_cmb[pc$clust==3,]))/nrow(data_cmb[pc$clust==3,])
c4=100*colSums(data.matrix(data_cmb[pc$clust==4,]))/nrow(data_cmb[pc$clust==4,])
c5=100*colSums(data.matrix(data_cmb[pc$clust==5,]))/nrow(data_cmb[pc$clust==5,])

png("a.png", width=1400)
barplot(rbind(c1,c2,c3,c4,c5), beside=T, legend=c("Cluster 1","Cluster 2","Cluster 3","Cluster 4","Cluster 5"))
dev.off()

pie(colSums(data_cmb[pc$clust==1,]))
####### End of bar plot


####### Boxplot
wt=c(1,1,1,1,1,1,1,1,1,1,2,2,2,2,3,6,6)
mat=data.matrix(data_cmb)
ci=mat%*%wt

#ci=charlson(v)

library(ggplot2)
ggplot(data=pc, aes(x=pc1,y=pc2,color=ci))+geom_point(alpha=0.1) +scale_color_gradient(low="white",high="black")

aar=grep("f21022_0", names(df), value = T)
aad=grep("f40007_0", names(df), value = T)
sex=grep("f31_0", names(df), value=T)
mde=grep("f26410_0", names(df), value=T)
demo=df[c("eid",aar,aad,sex,mde)]
age=demo[demo$eid %in% row.names(data_cmb),]
colnames(age)=c("eid","aar","aad","sex","mde")
age$clust=pc$clust
age$ci=ci
boxplot(aar~clust, data=age)
age.aov=aov(aar~clust, data=age)
summary(age.aov)
table(age$sex,age$clust)
round(prop.table(table(age$sex,age$clust),2),2)
plot(prop.table(table(age$sex,age$clust),2)[1,], type='l')
boxplot(ci~clust, data=age)
mean(age[age$sex=="Male",]$mde, na.rm=T)
### End of boxplot

fviz_mca_ind(res.mca, label="none", habillage=c("Paralysis", "Stroke", "Dementia"), palette=rep(c("red","blue"),3))

###### Subgroup
a=v[v$eid %in% rownames(data_cmb),]
a$G041=(rowSums(sapply(a, function(i) grepl("G041",i)))>0)
a$G114=(rowSums(sapply(a, function(i) grepl("G114",i)))>0)
a$G80=(rowSums(sapply(a, function(i) grepl("G80",i)))>0)
a$G81=(rowSums(sapply(a, function(i) grepl("G81",i)))>0)
a$G82=(rowSums(sapply(a, function(i) grepl("G82",i)))>0)
a$G83=(rowSums(sapply(a, function(i) grepl("G83",i)))>0)
a$G45=(rowSums(sapply(a, function(i) grepl("G45",i)))>0)
a$G46=(rowSums(sapply(a, function(i) grepl("G46",i)))>0)
a$H34=(rowSums(sapply(a, function(i) grepl("H34",i)))>0)
a$I60=(rowSums(sapply(a, function(i) grepl("I60",i)))>0)
a$I62=(rowSums(sapply(a, function(i) grepl("I62",i)))>0)
a$I63=(rowSums(sapply(a, function(i) grepl("I63",i)))>0)
a$I64=(rowSums(sapply(a, function(i) grepl("I64",i)))>0)
a$I65=(rowSums(sapply(a, function(i) grepl("I65",i)))>0)
a$I66=(rowSums(sapply(a, function(i) grepl("I66",i)))>0)
a$I67=(rowSums(sapply(a, function(i) grepl("I67",i)))>0)
a$I68=(rowSums(sapply(a, function(i) grepl("I68",i)))>0)
a$I69=(rowSums(sapply(a, function(i) grepl("I69",i)))>0)
a$F00=(rowSums(sapply(a, function(i) grepl("F00",i)))>0)
a$F01=(rowSums(sapply(a, function(i) grepl("F01",i)))>0)
a$F02=(rowSums(sapply(a, function(i) grepl("F02",i)))>0)
a$F03=(rowSums(sapply(a, function(i) grepl("F03",i)))>0)
a$F05=(rowSums(sapply(a, function(i) grepl("F05",i)))>0)
a$G30=(rowSums(sapply(a, function(i) grepl("G30",i)))>0)
a$G31=(rowSums(sapply(a, function(i) grepl("G31",i)))>0)
data_cmb2=a[,225:249]
c1=100*colSums(data.matrix(data_cmb2[pc$clust==1,]))/nrow(data_cmb[pc$clust==1,])
c2=100*colSums(data.matrix(data_cmb2[pc$clust==2,]))/nrow(data_cmb[pc$clust==2,])
c3=100*colSums(data.matrix(data_cmb2[pc$clust==3,]))/nrow(data_cmb[pc$clust==3,])
c4=100*colSums(data.matrix(data_cmb2[pc$clust==4,]))/nrow(data_cmb[pc$clust==4,])
c5=100*colSums(data.matrix(data_cmb2[pc$clust==5,]))/nrow(data_cmb[pc$clust==5,])
png("s.png",width=1400)
barplot(rbind(c1,c2,c3,c4,c5), beside=T, legend=c("Cluster 1","Cluster 2","Cluster 3","Cluster 4","Cluster 5"))
dev.off()

library(plot.matrix)
tt=matrix(0,25,25)
pp=matrix(0,25,25)
rownames(tt)=colnames(data_cmb2)
colnames(tt)=colnames(data_cmb2)
rownames(pp)=colnames(data_cmb2)
colnames(pp)=colnames(data_cmb2)
data_cmb3=data_cmb2[pc$clust==1,]
for (i in 1:25){
  for (j in 1:25){
    res1=table(data_cmb3[[colnames(data_cmb3)[i]]], data_cmb3[[colnames(data_cmb3)[j]]])
    a=try(res1["TRUE","TRUE"],silent=T)
    if(class(a)=="integer"){
       tt[j,i]=a
       pp[j,i]=a/nrow(data_cmb3)
}}}
png("b.png", width=1400, height=1400)
plot(pp, col=rev(gray.colors(20)))
dev.off()

######End of subdisease


#####n/w
options(stringsAsFactors=F)
node=data.frame()
tt=matrix(0,17,17)
rownames(tt)=colnames(data_cmb)
colnames(tt)=colnames(data_cmb)
data_cmb1=data_cmb[pc$clust==1,]
for (i in 1:16){
  for (j in (i+1):17){
    res1=table(data_cmb1[[colnames(data_cmb1)[i]]], data_cmb1[[colnames(data_cmb1)[j]]])
    a=try(res1["TRUE","TRUE"],silent=T)
    if(class(a)=="integer"){
      tt[j,i]=a
      b=a/nrow(data_cmb1)
      row=c(colnames(data_cmb)[i],colnames(data_cmb)[j],b)
      node=rbind(node,row)

}}}
write.table(node, file="mat.tsv", row.names=F, quote=F)

c1=colSums(data_cmb[pc$clust==1,])/nrow(data_cmb[pc$clust==1,])
c2=colSums(data_cmb[pc$clust==2,])/nrow(data_cmb[pc$clust==2,])
c3=colSums(data_cmb[pc$clust==3,])/nrow(data_cmb[pc$clust==3,])
c4=colSums(data_cmb[pc$clust==4,])/nrow(data_cmb[pc$clust==4,])
c5=colSums(data_cmb[pc$clust==5,])/nrow(data_cmb[pc$clust==5,])
c=rbind(c1,c2,c3,c4,c5)
write.table(c, file="num.tsv", row.names=F, quote=F)

