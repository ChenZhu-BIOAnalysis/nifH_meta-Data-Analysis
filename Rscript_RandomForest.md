## Analysis for RandomForest (The code for the Figure2)
```{r cv, echo=TRUE}
setwd("E:/R/meta_nifh/random_forest_2/")
design = read.table("env2.txt", header=T, row.names=1, sep="\t")
design$group = design$groupID
otutab = read.table("otutab.txt", header=T, row.names=1, sep="\t")
otutab[1:10,1:10]

#otutab <- otutab[which(rowSums(otutab) >= 0.0005), ]

#Randomforest classification apply to training set
# install.packages("randomForest")
library(randomForest)
library(rfPermute)
library(rfUtilities)

otutab_t = as.data.frame(t(otutab))
otutab_t[1:10,1:10]

train <- sample(nrow(otutab_t), nrow(otutab_t)*0.7)
otu_train <- otutab_t[train, ]
otu_train[1:10,1:10]
otu_test <- otutab_t[-train, ]
otu_test[1:10,1:10]


### otu_train_modeling####
otutab_t=otu_train
# Set classification info.
otutab_t$group = factor(design[rownames(otutab_t),]$group)

# set random seed for reproducible
set.seed(315)

# RandomForest Classification
otutab_t.rf= randomForest(group ~ ., data=otutab_t, importance=TRUE, proximity=TRUE)
print(otutab_t.rf)

##
varImpPlot(otutab_t.rf, main = "Top feature importance", n.var = 23)
write.table(otutab_t.rf$confusion, file = "genus_confusion.txt", sep = "\t", quote = F, row.names = T, col.names = T)
imp = as.data.frame(round(importance(otutab_t.rf), 2))
imp=imp[order(imp$MeanDecreaseAccuracy,decreasing = F),]
write.table(imp, file = "genus_imp.txt", sep = "\t", quote = F, row.names = T, col.names = NA)

##
###ntree = 500
set.seed(315)
otutab_t.rfP <- rfPermute(group ~ ., data = otutab_t, ntree = 500, na.action = na.omit, nrep = 50, num.cores = 1)
frimp.scaled <- rp.importance(otutab_t.rfP, scale = TRUE)
frimp.scaled
write.csv(frimp.scaled,"RandomForest_individual_Pvalue.csv")


otutab_t <- na.omit(otutab_t)
set.seed(315)
frrf.mdl <- randomForest(x=otutab_t[,1:(ncol(otutab_t)-1)], y=otutab_t[,ncol(otutab_t)], ntree = 500)
frrf.perm <- rf.significance(frrf.mdl, otutab_t[,1:(ncol(otutab_t)-1)], nperm=99, ntree=500)


#Cross validation to choose porper number of features

# great time consumption, need 5-30 minutes
n = ncol(otutab_t)-1
myotutab_t= otutab_t[1:n]
set.seed(315)
result= rfcv(myotutab_t, otutab_t$group, cv.fold=5, scale = "log", step = 0.9)
# result$n.var
# length(result$n.var)
with(result, plot(n.var, error.cv, log="x", type="o", lwd=2))

result1 = result
error.cv = data.frame(num = result$n.var, error.1 =  result$error.cv)
error1.cv=error.cv
write.csv(error1.cv,"error1_cv.csv")
for (i in 316:(315+4)){
  print(i)
  set.seed(i)
  result= rfcv(myotutab_t, otutab_t$group, cv.fold=5, scale = "log", step = 0.9)
  error.cv = cbind(error.cv, result$error.cv)
}
error2.cv=error.cv
write.csv(error2.cv,"error2_cv.csv")

n.var = error.cv$num
error.cv = error.cv[,2:6]
colnames(error.cv) = paste('err',1:5,sep='.')
err.mean = apply(error.cv,1,mean)
allerr = data.frame(num=n.var,err.mean=err.mean,error.cv)
write.table(allerr, file = "genus_rfcv.txt", sep = "\t", quote = F, row.names = T, col.names = NA)

### number of features selected,number determined by allerr
optimal = 32
library(ggplot2)
p = ggplot() + 
  geom_line(aes(x = allerr$num, y = allerr$err.1), colour = 'grey') + 
  geom_line(aes(x = allerr$num, y = allerr$err.2), colour = 'grey') + 
  geom_line(aes(x = allerr$num, y = allerr$err.3), colour = 'grey') + 
  geom_line(aes(x = allerr$num, y = allerr$err.4), colour = 'grey') + 
  geom_line(aes(x = allerr$num, y = allerr$err.5), colour = 'grey') + 
  geom_line(aes(x = allerr$num, y = allerr$err.mean), colour = 'black') + 
  geom_vline(xintercept = optimal, colour='black', lwd=0.36, linetype="dashed") + 
  #  geom_hline(yintercept = min(allerr$err.mean), colour='black', lwd=0.36, linetype="dashed") + 
  coord_trans(x = "log2") +
  scale_x_continuous(breaks = c(1, 2, 5, 10, 20, 30, 50, 100, 200)) + # , max(allerr$num)
  labs(title=paste('Training set (n = ', dim(otutab_t)[1],')', sep = ''), 
       x='Number of genus', y='Cross-validation error rate') + 
  annotate("text", x = optimal, y = max(allerr$err.mean), label=paste("optimal = ", optimal, sep="")) + 
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.line.x=element_line(size=.5, colour="black"),
        axis.line.y=element_line(size=.5, colour="black"),
        axis.ticks=element_line(color="black"),
        axis.text=element_text(color="black", size=7),
        legend.position="right",
        legend.background=element_blank(),
        legend.key=element_blank(),
        legend.text= element_text(size=7),
        text=element_text( size=7))
p
ggsave(p, file = "genus_rfcv.pdf", width = 89, height = 50, unit = 'mm')

##
taxonomy = read.table("abun.txt", header=T, row.names=1, sep="\t")
head(taxonomy)
taxonomy=taxonomy[,c(4,1)]
##taxonomy=taxonomy[!duplicated(taxonomy$Species), ]
imp$OTU<-rownames(imp)
a=merge(imp,taxonomy,by.x="OTU",all=TRUE)
write.table(a, file = "genus_imp_phylum.txt", sep = "\t", quote = F, row.names = T, col.names = NA)

optimal = 32

imp = read.table("genus_imp_phylum.txt", header=T, row.names= 1, sep="\t")
imp1=imp[order(imp[,9]),]
imp2 = tail(imp1, n = optimal)
imp2$Genus = factor(imp2$OTU, levels = (imp2$OTU))

p = ggplot(imp2, aes(x = Genus, y = MeanDecreaseAccuracy, fill = Class)) + 
  geom_bar(stat = "identity") + 
  coord_flip() +
  theme(panel.background=element_blank(),
        panel.grid=element_blank(),
        axis.line.x=element_line(size=.5, colour="black"),
        axis.line.y=element_line(size=.5, colour="black"),
        axis.ticks=element_line(color="black"),
        axis.text=element_text(color="black", size=7),
        legend.position="right",
        legend.background=element_blank(),
        legend.key=element_blank(),
        legend.text= element_text(size=7),
        text=element_text(size=7))
p

 ##c. otu_test Test 验证####
 #  Select by manual set group
 otutab =otu_test
 otutab[1:10,1:10]
 
  
 idx = rownames(design) %in% rownames(otutab)
 sub_design = design[idx,]
 head(sub_design)
 sub_otutab = otutab
 head(sub_otutab)[,1:10]
 
 otutab_t = as.data.frame(sub_otutab)
 head(otutab_t)[,1:10]
 otutab_t$group = factor(design[rownames(otutab_t),]$group)
 head(otutab_t)[,1:10]

 
 set.seed(315)
 library(randomForest)
 otutab.pred = predict(otutab_t.rf, otutab_t )  
 pre_tab = table(observed=otutab_t[,"group"], predicted=otutab.pred) 
 pre_tab
 # save prediction result
 predict = data.frame(group = otutab_t[,"group"], predicted=otutab.pred)
 write.table(predict, file = "RF_prediction.txt", quote = F, row.names = T, col.names = NA, sep = "\t")
 head(predict)
 
 predict$result = ifelse(predict$group == predict$predicted, 1, 0)
 predict$predict = ifelse(predict$predicted == "crust" | predict$predicted == "Grassland"| predict$predicted == "Farmland"| predict$predicted == "Forest" | predict$predicted == "Tundra"| predict$predicted == "Mine"| predict$predicted == "riparian",1,2)
 
 predict$predict = ifelse(predict$predicted == "crust",1,ifelse(predict$predicted == "Grassland",2,ifelse(predict$predicted == "Farmland",3,ifelse(predict$predicted == "Forest",4,ifelse(predict$predicted == "Tundra",5,ifelse(predict$predicted == "Mine",6,7))))))
 
 

 predict$predict1 = ifelse(predict$result == "1", predict$predict,0)
 head(predict)
 
##pheatmap drawing for predict
column = 26
row = round(length(predict$predict1)/column + 0.5)
i = column * row - length(predict$predict1)
predict1 = c(predict$predict1, rep(NA, i))
matrix = matrix(predict1, nrow = row, ncol = column, byrow = T)

pheatmap(matrix, color = c("#040000","#F9766E", "#00BFC4","#02759E") , cluster_rows = F, cluster_cols = F, cellwidth = 15, cellheight = 12, 
          filename = "Genus_test.pdf")
 
pheatmap(matrix, color = c("#040000","#0077BB","#33BBEE","#EE7733","#009988","#CC3311","#EE3377","#BBBBBB") , cluster_rows = F, cluster_cols = F, cellwidth = 15, cellheight = 12, filename = "Genus_test111.pdf")



