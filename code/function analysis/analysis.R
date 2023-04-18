<<<<<<< HEAD
## R version for analysis is 4.0.1 ###
=======
## 使用R的版本为4.0.1,建议使用相同版本，以免出现包不兼容等问题 ###
## 如提示缺少某某函数，请先安装对应的包 ##########################
>>>>>>> 5b969249d9ccc5d24eeb90abd7b72afdc772d8f9

############################ 分期验证 ############################
library(Hmisc)
library(ggpubr)
library(ggbeeswarm)

<<<<<<< HEAD
### validation set ####
=======
### 训练集 ####
>>>>>>> 5b969249d9ccc5d24eeb90abd7b72afdc772d8f9
m <-read.csv("E:\\LUAD_model\\anova.csv",head=T,sep=',')
m <-read.csv("E:\\LUAD_model\\anova1.csv",head=T,sep=',')
m <-read.csv("E:\\LUAD_model\\anova2.csv",head=T,sep=',')
m <-read.csv("E:\\LUAD_model\\anova3.csv",head=T,sep=',')
ggplot(m, aes(x=stage, y=score),color=stage) +
  geom_quasirandom(aes(colour=factor(stage)))+stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="crossbar", width=0.5)+stat_compare_means(method = "anova")

<<<<<<< HEAD
### validation set1 ###
=======
### 验证集1 ###
>>>>>>> 5b969249d9ccc5d24eeb90abd7b72afdc772d8f9
m <-read.csv("E:\\LUAD_model\\otherdata\\GSE31210\\anova.csv",head=T,sep=',')
m <-read.csv("E:\\LUAD_model\\otherdata\\GSE31210\\anova1.csv",head=T,sep=',')
m <-read.csv("E:\\LUAD_model\\otherdata\\GSE31210\\anova2.csv",head=T,sep=',')
m <-read.csv("E:\\LUAD_model\\otherdata\\GSE31210\\anova3.csv",head=T,sep=',')
ggplot(m, aes(x=stage, y=score),color=stage) +
  geom_quasirandom(aes(colour=factor(stage)))+stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="crossbar", width=0.5)+stat_compare_means(method = "anova")

<<<<<<< HEAD
### validation set2 ###
=======
### 验证集2 ###
>>>>>>> 5b969249d9ccc5d24eeb90abd7b72afdc772d8f9
m <-read.csv("E:\\LUAD_model\\otherdata\\validation\\anova.csv",head=T,sep=',')
m <-read.csv("E:\\LUAD_model\\otherdata\\validation\\anova1.csv",head=T,sep=',')
m <-read.csv("E:\\LUAD_model\\otherdata\\validation\\anova2.csv",head=T,sep=',')
m <-read.csv("E:\\LUAD_model\\otherdata\\validation\\anova3.csv",head=T,sep=',')
m$stage <- factor(m$stage,levels=c("WELL","Moderate","POOR"))
ggplot(m, aes(x=stage, y=score),color=stage) +
  geom_quasirandom(aes(colour=factor(stage)))+stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), geom="crossbar", width=0.5)+stat_compare_means(method = "anova")

<<<<<<< HEAD
############################# survival ###################################
=======
############################# 生存分析 ###################################
>>>>>>> 5b969249d9ccc5d24eeb90abd7b72afdc772d8f9
library(survival)
library(survminer)
survial <- read.csv("E:\\LUAD_model\\survialdata.csv",head=T,sep=',',stringsAsFactors = FALSE) #TCGA-LUAD
survial <- read.csv("E:\\LUAD_model\\otherdata\\validation\\survialdata.csv",head=T,sep=',',stringsAsFactors = FALSE) #GSE68465
survial <- read.csv("E:\\LUAD_model\\otherdata\\GSE31210\\survialdata.csv",head=T,sep=',',stringsAsFactors = FALSE) #GSE31210
sfit <- survfit(Surv(time, status)~cluster, data=survial)
ggsurvplot(sfit, risk.table=TRUE, conf.int=FALSE, pval=TRUE, xlab="Time", ggtheme = theme_light(),ncensor.plot= TRUE,palette = c("#FF0000","#00FF00","#00002C","#FF1AB9","#FFD300","#005800","#8484FF","#9E4F46","#00FFC1")) #conf.int为TRUE则绘制置信区间，risk.table绘制风险表，


<<<<<<< HEAD
########################## enrichR ########################################
=======
########################## enrichR富集分析########################################
>>>>>>> 5b969249d9ccc5d24eeb90abd7b72afdc772d8f9
library(enrichR)
symbol <- read.csv("E:\\LUAD_model\\symbol_ID2.csv",head=F,sep=',')

enrichr <- enrichr(symbol[,1],databases = c('GO_Biological_Process_2018',
                                            'GO_Molecular_Function_2018',
                                            'GO_Cellular_Component_2018',
                                            'Reactome_2016'))
enrichr$GO_Biological_Process_2018 <- dplyr::filter(enrichr$GO_Biological_Process_2018,Adjusted.P.value <0.05)
enrichr$GO_Molecular_Function_2018 <- dplyr::filter(enrichr$GO_Molecular_Function_2018,Adjusted.P.value <0.05)
enrichr$GO_Cellular_Component_2018 <- dplyr::filter(enrichr$GO_Cellular_Component_2018,Adjusted.P.value <0.05)
enrichr$Reactome_2016 <- dplyr::filter(enrichr$Reactome_2016,Adjusted.P.value <0.05)
for (i in 1:50) {
  enrichr$GO_Biological_Process_2018$Term[i] <- substr(enrichr$GO_Biological_Process_2018$Term[i],1,regexpr("\\(",enrichr$GO_Biological_Process_2018$Term[i])-2)
}
for (i in 1:50) {
  enrichr$GO_Molecular_Function_2018$Term[i] <- substr(enrichr$GO_Molecular_Function_2018$Term[i],1,regexpr("\\(",enrichr$GO_Molecular_Function_2018$Term[i])-2)
}
for (i in 1:50) {
  enrichr$GO_Cellular_Component_2018$Term[i] <- substr(enrichr$GO_Cellular_Component_2018$Term[i],1,regexpr("\\(",enrichr$GO_Cellular_Component_2018$Term[i])-2)
}
for (i in 1:50) {
  enrichr$Reactome_2016$Term[i] <- substr(enrichr$Reactome_2016$Term[i],1,regexpr("Homo sapiens",enrichr$Reactome_2016$Term[i])-2)
}

# write.csv(enrichr$Reactome_2016,"E:\\LUAD_model\\enrich\\Reactome2016.csv")
plotEnrich(enrichr[[1]], showTerms = 20, numChar = 70, y = "Count",orderBy = "P.value") + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
plotEnrich(enrichr[[2]], showTerms = 20, numChar = 70, y = "Count",orderBy = "P.value") + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
plotEnrich(enrichr[[3]], showTerms = 20, numChar = 70, y = "Count",orderBy = "P.value") + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
plotEnrich(enrichr[[4]], showTerms = 20, numChar = 70, y = "Count",orderBy = "P.value") + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())

############################ GRN for 314 genes #######################################
fea <- read.csv("E:\\LUAD_model\\GRN\\cell_cycle_genes2.csv",head=F,sep=',')
colname <- read.csv("E:\\LUAD_model\\GRN\\colname.csv",head=F,sep=',') 
colname <- as.matrix(colname)
symbol_ID1 <- read.csv("E:\\LUAD_model\\gene_id\\symbol_ID.csv",head=F,sep=',')
table(fea[,1] %in% symbol_ID1[,1])
data <- read.csv("E:\\LUAD_model\\gene_id\\gene_expression_all1.csv",head=T,row.names = 1,sep=',')
sample_order <- read.csv("E:\\LUAD_model\\GRN\\sample_order.csv",head=F,sep=',')
sample_order <- as.matrix(sample_order)


data353 <- data[match(fea[,1],symbol_ID1[,1]),]
data353 <- na.omit(data353)
data353 <- data353[,-c(6,24,87,143,170,183,212,262,413,485,516,529)]
data353 <- data353[,sample_order]
data353 <- log2(data353) #使用以2为底log函数
data353[data353==-Inf] <- 0
fea <- as.matrix(fea)
dimnames(data353) <- list(fea,colname)
#write.csv(data353,"E:\\LUAD_model\\GRN\\GRN_input_data_new2.csv")

inference <- read.table("E:\\LUAD_model\\GRN\\GRNVBEM\\GRNVBEM-master\\GRN_input_data_new2__AR1MA1_GRN_inference.txt",head=T,sep=' ') 
inference <- inference[which(inference[,4]>0.3 | inference[,4] < -0.3),]
#write.table(inference, "E:\\LUAD_model\\GRN\\GRNVBEM\\GRNVBEM-master\\GRN_input_data_new2__AR1MA1_GRN_inference(1).txt")

<<<<<<< HEAD
############################# BUB1B expression with the change of pseudo time #########################
=======
############################# BUB1B表达随伪时间值变化 #########################
>>>>>>> 5b969249d9ccc5d24eeb90abd7b72afdc772d8f9
library(Hmisc)
library(ggpubr)
library(ggbeeswarm)
data <- read.csv("E:\\LUAD_model\\gene_id\\gene_expression.csv",head=F,sep=',')
symbol_ID1 <- read.csv("E:\\LUAD_model\\gene_id\\symbol_ID.csv",head=F,sep=',')
sample_order <- read.csv("E:\\LUAD_model\\GRN\\sample_order.csv",head=F,sep=',')
sample_order <- as.matrix(sample_order)
data353 <- data[match('BUB1B',symbol_ID1[,1]),] #BUB1B  FAM13A
data353 <- na.omit(data353)
data353 <- data353[,-c(6,24,87,143,170,183,212,262,413,485,516,529)]
data353 <- data353[,sample_order]
y <- as.matrix(data353[1,])
x <- read.csv("E:\\LUAD_model\\fscore_521.csv",head=F,sep=',')
x <- (x-x[1,])/(x[521,]-x[1,])
x <- as.matrix(x)
my.data <- data.frame(x,y)
b <- ggplot(my.data, aes(x = x, y = y)) + xlab("Pseudotime value") + ylab("Expression value")+theme_bw()
b <- b+ geom_point()  + geom_smooth(method = "lm",formula=y~I(x*x), color = "blue", fill = "yellow",size = 1)+stat_cor(method = "pearson")
b <- b + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
b

############################# BUB1B survial #####################################
library(survival)
library(survminer)
symbol <- read.csv("E:\\LUAD_model\\gene_id\\symbol_ID.csv",head=F,sep=',')
survial <- read.csv("E:\\LUAD_model\\GRN\\survial.csv",head=T,sep=',')
data <- read.csv("E:\\LUAD_model\\gene_id\\gene_expression.csv",head=F,sep=',')
data <- as.matrix(data)
data1 <- data[10038,] #BUB1B
survial[which(data1 > median(data1)),3] <- 1
survial[which(data1 <= median(data1)),3] <- 0
names(survial) <- c("time","status","cluster")
survival <- na.omit(survial)
survival[which(survival[,1] > 1825),1] <- 1825
survival[which(survival[,1] > 1825),2] <- 1
survival[,1] <- survival[,1]/365

sfit <- survfit(Surv(time, status)~cluster, data=survival)
ggsurvplot(sfit, risk.table=FALSE, conf.int=FALSE, pval=TRUE, xlab="Time", ggtheme = theme_light(),ncensor.plot= FALSE) #conf.int为TRUE则绘制置信区间，risk.table绘制风险表，

<<<<<<< HEAD
################################## gene expression correlation ####################################
=======
################################## 基因表达相关性 ####################################
>>>>>>> 5b969249d9ccc5d24eeb90abd7b72afdc772d8f9
symbol_ID1 <- read.csv("E:\\LUAD_model\\gene_id\\symbol_ID.csv",head=F,sep=',')
data <- read.csv("E:\\LUAD_model\\gene_id\\gene_expression.csv",head=F,sep=',')

data1 <- data[match('BUB1B',symbol_ID1[,1]),] 
data1 <- data1[,-c(6,24,87,143,170,183,212,262,413,485,516,529)]

data2 <- data[match('BUB1',symbol_ID1[,1]),] 
data2<- data2[,-c(6,24,87,143,170,183,212,262,413,485,516,529)]

res <- rbind(data1,data2)
res <- t(res)
res <- as.data.frame(res)

library(Hmisc)
library(ggpubr)
library(ggbeeswarm)
b <- ggplot(res, aes(x = res[,1], y = res[,2])) + xlab("BUB1B Expression value") + ylab("BUB1 Expression value")+theme_bw()
#b <- b+ geom_point()  + geom_smooth(method = "lm",formula=y~I(x*x), color = "blue", fill = "yellow",size = 1)+stat_cor(method = "pearson")
b <- b+ geom_point()  + geom_smooth(method = "lm",formula=y~x, color = "blue", fill = "yellow",size = 1)+stat_cor(method = "pearson")
b <- b + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
b

##############
data1 <- data[match('BUB1B',symbol_ID1[,1]),] #BUB1B  FAM13A
data1 <- data1[,-c(6,24,87,143,170,183,212,262,413,485,516,529)]

data2 <- data[match('BUB3',symbol_ID1[,1]),] #BUB1B  FAM13A
data2<- data2[,-c(6,24,87,143,170,183,212,262,413,485,516,529)]

res <- rbind(data1,data2)
res <- t(res)
res <- as.data.frame(res)
b <- ggplot(res, aes(x = res[,1], y = res[,2])) + xlab("BUB1B Expression value") + ylab("BUB3 Expression value")+theme_bw()
#b <- b+ geom_point()  + geom_smooth(method = "lm",formula=y~I(x*x), color = "blue", fill = "yellow",size = 1)+stat_cor(method = "pearson")
b <- b+ geom_point()  + geom_smooth(method = "lm",formula=y~x, color = "blue", fill = "yellow",size = 1)+stat_cor(method = "pearson")
b <- b + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
b

##############
data1 <- data[match('BUB1',symbol_ID1[,1]),] #BUB1B  FAM13A
data1 <- data1[,-c(6,24,87,143,170,183,212,262,413,485,516,529)]

data2 <- data[match('STAT3',symbol_ID1[,1]),] #BUB1B  FAM13A
data2<- data2[,-c(6,24,87,143,170,183,212,262,413,485,516,529)]

res <- rbind(data1,data2)
res <- t(res)
res <- as.data.frame(res)
b <- ggplot(res, aes(x = res[,1], y = res[,2])) + xlab("BUB1 Expression value") + ylab("STAT3 Expression value")+theme_bw()
#b <- b+ geom_point()  + geom_smooth(method = "lm",formula=y~I(x*x), color = "blue", fill = "yellow",size = 1)+stat_cor(method = "pearson")
b <- b+ geom_point()  + geom_smooth(method = "lm",formula=y~x, color = "blue", fill = "yellow",size = 1)+stat_cor(method = "pearson")
b <- b + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
b

##############
data1 <- data[match('STAT3',symbol_ID1[,1]),] #BUB1B  FAM13A
data1 <- data1[,-c(6,24,87,143,170,183,212,262,413,485,516,529)]

data2 <- data[match('BUB3',symbol_ID1[,1]),] #BUB1B  FAM13A
data2<- data2[,-c(6,24,87,143,170,183,212,262,413,485,516,529)]

res <- rbind(data1,data2)
res <- t(res)
res <- as.data.frame(res)
b <- ggplot(res, aes(x = res[,1], y = res[,2])) + xlab("STAT3 Expression value") + ylab("BUB3 Expression value")+theme_bw()
#b <- b+ geom_point()  + geom_smooth(method = "lm",formula=y~I(x*x), color = "blue", fill = "yellow",size = 1)+stat_cor(method = "pearson")
b <- b+ geom_point()  + geom_smooth(method = "lm",formula=y~x, color = "blue", fill = "yellow",size = 1)+stat_cor(method = "pearson")
b <- b + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
b

<<<<<<< HEAD
############################### variant num with pseudo time change ###########################
=======
############################### 突变数量随伪时间值的变化 ###########################
>>>>>>> 5b969249d9ccc5d24eeb90abd7b72afdc772d8f9
library(Hmisc)
library(ggpubr)
library(ggbeeswarm)
number <- read.csv("E:\\LUAD_model\\mut_number.csv",head=F,sep=',')
score<- read.csv("E:\\LUAD_model\\mut_score.csv",head=F,sep=',')
number <- as.matrix(number)
score <- as.matrix(score)
score <- (score-min(score))/(max(score)-min(score))

plot(score,number)
df <- cbind(score,number)
df <- as.data.frame(df)
colnames(df) <- c("score","number")
b <- ggplot(df, aes(x = score, y = number)) + xlab("score") + ylab("mutation number")
b <- b+ geom_point(size=2,alpha=0.7)  + geom_smooth(method = "lm",formula=y~I(x*x), color = "blue", fill = "yellow",size = 1.5)+stat_cor(method = "pearson")
b <- b + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank(),panel.background = element_blank(),axis.line = element_line(colour = "black"))
b

<<<<<<< HEAD
######################## construct subclone ######################################
=======
######################## 构建亚克隆 ######################################
>>>>>>> 5b969249d9ccc5d24eeb90abd7b72afdc772d8f9
v = read.table("E:\\LUAD_model\\subClone\\data\\vcf.txt",header=T)
c = read.table("E:\\LUAD_model\\subClone\\data\\copyNumber.txt",header=T)
sample = read.table("E:\\LUAD_model\\subClone\\data\\sample_del.txt",header=F)
outPath <- "E:\\LUAD_model\\subClone"
out_clusters <- paste("clusters", ".txt", sep='')
out_pdf <- paste("result", ".pdf", sep='')
out_summary <- paste("summary", ".txt", sep='')
for(i in 1:501){
  sample_i = sample[i,]
  vaf <- v[v$Sample == sample_i,c(2,3,4,5,6)]
  cn <- c[c$sample == sample_i,c(2,3,4,5)]
  sc = sciClone(vafs=list(vaf),
                sampleNames= sample_i,
                copyNumberCalls=list(cn),
                cnCallsAreLog2=TRUE,
                annotation = list(v[v$Sample == sample_i,c(2,3,7)]),
                useSexChrs = TRUE,
                doClustering = TRUE,
                clusterMethod = 'bmm',
                copyNumberMargins = 0.25,
                minimumDepth = 30)
  out_fileName1 <- paste(i,out_clusters,sep='_')
  out_fileName2 <- paste(i,out_pdf,sep='_')
  out_fileName3 <- paste(i,out_summary,sep='_')
  out_filePath1 <- paste(outPath,out_fileName1,sep='\\')
  out_filePath2 <- paste(outPath,out_fileName2,sep='\\')
  out_filePath3 <- paste(outPath,out_fileName3,sep='\\')
  writeClusterTable(sc,out_filePath1)
  sc.plot1d(sc,out_filePath2)
  writeClusterSummaryTable(sc,out_filePath3)
}
<<<<<<< HEAD
######################## subclone num ###############################
=======
######################## 统计亚克隆数量 ###############################
>>>>>>> 5b969249d9ccc5d24eeb90abd7b72afdc772d8f9
read_clusters <- paste("clusters", ".txt", sep='')
read_fileName <- matrix(0,501,1)
outPath <- "E:\\LUAD_model\\subClone"
for (i in 1:501) {
  fileNames <- paste(i,read_clusters,sep='_')
  read_fileName[i,] <- paste(outPath,fileNames,sep='\\')
}
data <- lapply(read_fileName, function(x){
  read.table(x, header=T)})
count <- matrix(0,501,1)
for(i in 1:501){
  gene1 <- subset(data[i][[1]],select = c("cluster"))
  count[i,1] <- length(unique(gene1[,1]))-1
}
count[count==0] <- 1
#write.csv(count,"E:\\LUAD_model\\subClone\\result\\subclone_number.csv")
count <- read.csv("E:\\LUAD_model\\subClone\\result\\subclone_number.csv")
p1<-ggplot(data=count, aes(x=disp, y=mpg)) +geom_point(color="#d7191c") +geom_smooth(method="lm",color="#1a9641") +geom_text(aes(x=400, y=32,label=paste("R","=",signif(r,3),seq="")),color="#fdae61")

<<<<<<< HEAD
##################### num of subclone with the change of pseudo time ###################################
=======
##################### 亚克隆数量随伪时间的变化 ###################################
>>>>>>> 5b969249d9ccc5d24eeb90abd7b72afdc772d8f9
library(Hmisc)
library(ggpubr)
library(ggbeeswarm)
m <-read.csv("E:\\LUAD_model\\subClone\\result\\subclone_number.csv",
             head=T,sep=',')
m$cluster <- as.factor(m$number)
m <- m[-151,]
m <- m[-239,]
g <- ggplot2::ggplot(m,ggplot2::aes(x=order,y=cluster,fill=cluster))
#g <- ggplot2::ggplot(m,ggplot2::aes(x=cluster,y=order,fill=cluster))
g <- g + ggplot2::geom_boxplot()
g <- g + ggplot2::labs(x = 'Pseudotime order',y='Subclone number',fill='Subclone\n Number')+theme_bw()
g <- g + theme(panel.grid.major=element_blank(),panel.grid.minor=element_blank())
g

<<<<<<< HEAD
##################### relationship between subclone and sample num ################################
=======
##################### 亚克隆数量与样本数量的关系 ################################
>>>>>>> 5b969249d9ccc5d24eeb90abd7b72afdc772d8f9
n <- matrix(0,9,2)
for (i in 1:9) {
  n[i,1] <- i
  n[i,2] <- length(which(m[,2]==i))
}
colnames(n) <- c("number","snumber")
n <- as.data.frame(n)
n$number <- as.factor(n$number)
ggplot(data=n,mapping=aes(x=number,y=snumber,fill=number,group=factor(1)))+
  geom_bar(stat="identity",width=0.5)+
  labs(x = 'Subclone number',y='Sample number',fill='Subclone\n Number') + theme_bw()

################################# path1 clone ###############################
read_clusters <- paste("clusters", ".txt", sep='')
path <- read.table("E:\\LUAD_model\\subClone\\data\\path1.txt")
read_fileName <- matrix(0,167,1)
outPath <- "E:\\LUAD_model\\subClone"
for (i in 1:167) {
  fileNames <- paste(path[i,1],read_clusters,sep='_')
  read_fileName[i,] <- paste(outPath,fileNames,sep='\\')
}
data <- lapply(read_fileName, function(x){
  read.table(x, header=T)})
luad <- read.table("E:\\LUAD_model\\subClone\\data\\CGC_genes.txt",head=F)
result <- matrix(0,1,2)
for (i in 1:153) {
  luad_i <- luad[i,1]
  for(j in 1:167){
    gene1 <- subset(data[j][[1]],select = c("cluster","Hugo_Symbol"))
    result1 <- matrix(0,1,2)
    if(luad_i %in% gene1[,2] & (gene1[match(luad_i,gene1[,2]),1] %in% 1)){
      result1[1,1] <- "clone"
      result1[1,2] <- luad_i
    } else if(!(gene1[match(luad_i,gene1[,2]),1] %in% NA)){
      result1[1,1] <- "subclone"
      result1[1,2] <- luad_i
    }else{
      next;
    }
    result <- rbind(result,result1)
  }
}
result <- result[-1,]
result <- as.data.frame(result)
names(result) <- c("clone","gene")
ggplot(result,aes(x=reorder(gene, as.numeric(clone=='clone')), fill = clone)) + #x轴的分类为clarity，填充颜色为color（J和H）
  theme_classic() + #设置主题
  geom_bar(position = position_fill()) +
  scale_fill_manual(values = c('#5B72E0', '#F18565')) +
  theme(axis.text.x = element_text(size=10),axis.text.y = element_text(size=6))+
  labs(y = 'Percent',x = 'Genes') + #设置y轴名为‘Percent’
  coord_flip() #旋转坐标轴


################################# path2 clone #################################
read_clusters <- paste("clusters", ".txt", sep='')
path <- read.table("E:\\LUAD_model\\subClone\\data\\path2.txt")
read_fileName <- matrix(0,135,1)
outPath <- "E:\\LUAD_model\\subClone"
for (i in 1:135) {
  fileNames <- paste(path[i,1],read_clusters,sep='_')
  read_fileName[i,] <- paste(outPath,fileNames,sep='\\')
}
data <- lapply(read_fileName, function(x){
  read.table(x, header=T)})
luad <- read.table("E:\\LUAD_model\\subClone\\data\\CGC_genes.txt",head=F)
result <- matrix(0,1,2)
for (i in 1:153) {
  luad_i <- luad[i,1]
  for(j in 1:135){
    gene1 <- subset(data[j][[1]],select = c("cluster","Hugo_Symbol"))
    result1 <- matrix(0,1,2)
    if(luad_i %in% gene1[,2] & (gene1[match(luad_i,gene1[,2]),1] %in% 1)){
      result1[1,1] <- "clone"
      result1[1,2] <- luad_i
    } else if(!(gene1[match(luad_i,gene1[,2]),1] %in% NA)){
      result1[1,1] <- "subclone"
      result1[1,2] <- luad_i
    }else{
      next;
    }
    result <- rbind(result,result1)
  }
}
result <- result[-1,]
result <- as.data.frame(result)
names(result) <- c("clone","gene")
result[which(result[,1]=="clone"),3] <- "subclone"
result[which(result[,1]=="subclone"),3] <- "clone"
names(result) <- c("clone","gene","subclone")
ggplot(result,aes(x=reorder(gene, as.numeric(subclone=='subclone')), fill = clone)) + #x轴的分类为clarity，填充颜色为color（J和H）
  geom_bar(position = position_fill()) +
  theme_classic() + #设置主题
  scale_fill_manual(values = c('#5B72E0', '#F18565')) +
  theme(axis.text.x = element_text(size=8),axis.text.y = element_text(size=6))+
  labs(y = 'Percent',x = 'Genes') + #设置y轴名为‘Percent’
  coord_flip() #旋转坐标轴
##################################### path3 clone ##########################
read_clusters <- paste("clusters", ".txt", sep='')
path <- read.table("E:\\LUAD_model\\subClone\\data\\path3.txt")
read_fileName <- matrix(0,265,1)
outPath <- "E:\\LUAD_model\\subClone"
for (i in 1:265) {
  fileNames <- paste(path[i,1],read_clusters,sep='_')
  read_fileName[i,] <- paste(outPath,fileNames,sep='\\')
}
data <- lapply(read_fileName, function(x){
  read.table(x, header=T)})
luad <- read.table("E:\\LUAD_model\\subClone\\data\\CGC_genes.txt",head=F)
result <- matrix(0,1,2)
for (i in 1:153) {
  luad_i <- luad[i,1]
  for(j in 1:265){
    gene1 <- subset(data[j][[1]],select = c("cluster","Hugo_Symbol"))
    result1 <- matrix(0,1,2)
    if(luad_i %in% gene1[,2] & (gene1[match(luad_i,gene1[,2]),1] %in% 1)){
      result1[1,1] <- "clone"
      result1[1,2] <- luad_i
    } else if(!(gene1[match(luad_i,gene1[,2]),1] %in% NA)){
      result1[1,1] <- "subclone"
      result1[1,2] <- luad_i
    }else{
      next;
    }
    result <- rbind(result,result1)
  }
}
result <- result[-1,]
result <- as.data.frame(result)
names(result) <- c("clone","gene")
ggplot(result,aes(x=reorder(gene, as.numeric(clone=='clone')), fill = clone)) + #x轴的分类为clarity，填充颜色为color（J和H）
  geom_bar(position = position_fill()) +
  theme_classic() + #设置主题
  scale_fill_manual(values = c('#5B72E0', '#F18565')) +
  theme(axis.text.x = element_text(size=8),axis.text.y = element_text(size=6))+
  labs(y = 'Percent',x = 'Genes') + #设置y轴名为‘Percent’
  coord_flip() #旋转坐标轴

