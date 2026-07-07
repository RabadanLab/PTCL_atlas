library(ggplot2)

original<-read.csv("output/subtypes_1.0_1.csv",row.names=1)
SEQ<-c(0.05,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.1,1.2,1.3,1.4,1.5)
str_SEQ<-c("0.05","0.1","0.2","0.3","0.4","0.5","0.6","0.7","0.8","0.9","1.0","1.1","1.2","1.3","1.4","1.5")
resolution<-sapply(str_SEQ,function(x){paste0("res",x)})

wtd_50<-NULL
wtd_60<-NULL
wtd_75<-NULL
wtd_90<-NULL

pdf("plot_jaccard.pdf")

for(RES in resolution){

	print(RES)
	array_index<-NULL
	array_value<-NULL

	for(ORDER in seq(1,100,1)){

		if(!file.exists(paste0("output/subtypes_0.9_",ORDER,".csv"))){
			next
		}

		downsample<-read.csv(paste0("output/subtypes_0.9_",ORDER,".csv"),row.names=1)
		primary<-original[rownames(original) %in% rownames(downsample),]

		count_table<-table(data.frame(primary[RES],downsample[RES]))

		jaccard_table<-count_table
		for(i in 1:nrow(count_table)){
			for(j in 1:ncol(count_table)){
				jaccard_table[i,j]=count_table[i,j]/(sum(count_table[i,])+sum(count_table[,j])-count_table[i,j]+0.000001)
			}
		}

		maximal_jaccard<-apply(jaccard_table,1,max)
		array_index<-c(array_index,names(maximal_jaccard))
		array_value<-c(array_value,as.vector(maximal_jaccard))
	}


        df<-data.frame(array_index,array_value)
        p<-ggplot(df,aes(x=as.numeric(array_index),y=array_value))+
	        theme_bw()+
        	#geom_boxplot(alpha=1,size=0.75,width=0.25)+
                #annotate("rect", xmin = 1, xmax = 1.5, ymin = 0, ymax = 1, alpha = 1, fill = 'white') +
                #annotate("rect", xmin = 2, xmax = 2.5, ymin = 0, ymax = 1, alpha = 1, fill = 'white') +
                #annotate("rect", xmin = 3, xmax = 3.5, ymin = 0, ymax = 1, alpha = 1, fill = 'white') +
                #annotate("rect", xmin = 4, xmax = 4.5, ymin = 0, ymax = 1, alpha = 1, fill = 'white') +
                #annotate("rect", xmin = 5, xmax = 5.5, ymin = 0, ymax = 1, alpha = 1, fill = 'white') +
                #annotate("rect", xmin = 6, xmax = 6.5, ymin = 0, ymax = 1, alpha = 1, fill = 'white') +
                #annotate("rect", xmin = 7, xmax = 7.5, ymin = 0, ymax = 1, alpha = 1, fill = 'white') +
                #annotate("rect", xmin = 8, xmax = 8.5, ymin = 0, ymax = 1, alpha = 1, fill = 'white') +
                #annotate("rect", xmin = 9, xmax = 9.5, ymin = 0, ymax = 1, alpha = 1, fill = 'white') +
                geom_point(aes(x=as.numeric(array_index), y=array_value),
                       alpha = 0.1, position = position_jitter(width = 0.1))+
                geom_hline(yintercept=0.5, linetype="dashed", color = "dark green")+
                geom_hline(yintercept=0.6, linetype="dashed", color = "blue")+
                geom_hline(yintercept=0.75, linetype="dashed", color = "dark red")+
		ggtitle(paste0("res",RES))

       	#pdf(paste0("./plots_jaccard/plot_jaccard_",RES,".pdf"))
        print(p)
        #dev.off()


	w<-NULL
	for(item in unique(array_index)){
		w<-c(w,sum(original[RES]==item)/nrow(original))
	}

	over_50<-NULL
	over_60<-NULL
	over_75<-NULL
	over_90<-NULL

	for(item in unique(array_index)){
		v<-array_value[array_index==item]
		over_50<-c(over_50,sum(v>0.5)/length(v))
		over_60<-c(over_60,sum(v>0.6)/length(v))
		over_75<-c(over_75,sum(v>0.75)/length(v))
		over_90<-c(over_90,sum(v>0.9)/length(v))
	}

	wtd_50<-c(wtd_50,sum(over_50*w))
	wtd_60<-c(wtd_60,sum(over_60*w))
	wtd_75<-c(wtd_75,sum(over_75*w))
	wtd_90<-c(wtd_90,sum(over_90*w))
}

dev.off()

pdf("plot_perc.pdf")
plot(SEQ,wtd_50,type="l",lwd=2,col="dark green",xlab="Cluster number",ylab="Percantage of jaccard score > 0.50")
points(SEQ,wtd_50,pch=8,col="dark green")
plot(SEQ,wtd_60,type="l",lwd=2,col="blue",xlab="Cluster number",ylab="Percantage of jaccard score > 0.60")
points(SEQ,wtd_60,pch=8,col="blue")
plot(SEQ,wtd_75,type="l",lwd=2,col="dark red",xlab="Cluster number",ylab="Percantage of jaccard score > 0.75")
points(SEQ,wtd_75,pch=8,col="dark red")
plot(SEQ,wtd_90,type="l",lwd=2,col="purple",xlab="Cluster number",ylab="Percantage of jaccard score > 0.90")
points(SEQ,wtd_90,pch=8,col="purple")
dev.off()

