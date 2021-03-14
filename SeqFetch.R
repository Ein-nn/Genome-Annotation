# 
args <- commandArgs(T)
if(length(args)==0)
{
  print("No arguments supplied.")
}else
{
  p <- as.numeric(args[1])
  t <- as.numeric(args[2])
  gff3 <- args[3]

  seqFetch<-function(gff3File,pace,times){  
    # fna文件读入
    fileName<-"./S288C.fna"
    # 读入序列，每个元素存入一行
    seq <- readLines(fileName)
    # 去除空行
    seq <- seq[seq != ""] 
    # 正则匹配（regular expression）注释行,是注释行为1，否则为-1
    is.anno <- regexpr("^>", seq, perl=T) 
    # 注释内容
    seq.anno <- seq[ which(is.anno == 1) ] 
    # 序列内容
    seq.content <- seq[ which(is.anno == -1) ] 
    ##--计算每条序列内容所占的行数，便于后来拼接--##
    # 注释行行号
    start <- which(is.anno == 1) 
    if(length(start)==1 && start == 1){
      end <- length(seq) 
    }else{
      # 第二条记录注释行到最后一条记录注释行行号减一，即为每条记录结束行号，这里会统计少一行——最后一行的结束未统计
      end <- start[ 2:length(start) ]-1 
      # 末尾添加一行：所有序列结束行  
      end <- c(end, length(seq) )  
    }
    # 每条记录所占行号
    distance <- end - start 
    # 生成一个一到记录总个数的向量
    index <- 1:length(start) 
    index <- rep(index, distance) # 分组标签
    index <- as.factor(index)
    seqs <- tapply(seq.content, index, paste, collapse="")  # 拼接每条序列内容，返回一个列表，列表每个元素为一条序列的内容
    #seqs[[1]]
    seq.content<-as.character( seqs ) # 将列表转换为向量，向量每个元素为一条序列的内容
    seq.len <- nchar(seq.content) # 获得序列长度
    seq.ID <- substring(seq.anno,2,12)
    seq.ID[4] <- "NC_001136.10"
    result <- data.frame( seq.ID, seq.content,stringsAsFactors=FALSE ) 
    #View(result) 
    #write.table(result,"fna_table.txt",sep = "\t",quote = FALSE,row.names = FALSE)
    
    #gffFile<-read.table("./test.gff3",header = FALSE,fill = TRUE,stringsAsFactors = FALSE,comment.char = "#")
    gffFile<-read.table(gff3File,header = FALSE,fill = TRUE,stringsAsFactors = FALSE,comment.char = "#")
    gene<-gffFile[which(gffFile$V3=="gene"),]
    #pace<-500
    k=1
    fa<-data.frame()
    for (i in 1:dim(gene)[1]) {
      chr = strsplit(gene[i,1], ";")
      for (j in 1:dim(result)[1]) {
        if(chr[[1]][1]==result[j,1]){
          if(gene[i,4]>pace & ( nchar(result[j,2])-gene[i,5])>pace){#延长范围在染色体内（不超过染色个体起始终止位置）
            fa[k,1]<-gene[i,1]
            fa[k,2]<-gene[i,4]
            fa[k,3]<-gene[i,5]
            fa[k,4]<-gene[i,7]
            fa[k,5]<-substring(result[j,2],gene[i,4]-pace,gene[i,5]+pace) 
            k=k+1
          }
          if(gene[i,4]<pace & ( nchar(result[j,2])-gene[i,5])>pace){#延长范围在超过染色体(start)
            fa[k,1]<-gene[i,1]
            fa[k,2]<-gene[i,4]
            fa[k,3]<-gene[i,5]
            fa[k,4]<-gene[i,7]
            fa[k,5]<-substring(result[j,2],1,gene[i,5]+pace) 
            k=k+1
          }
          if(gene[i,4]>pace & ( nchar(result[j,2])-gene[i,5])<pace){#延长范围在超过染色(end)
            fa[k,1]<-gene[i,1]
            fa[k,2]<-gene[i,4]
            fa[k,3]<-gene[i,5]
            fa[k,4]<-gene[i,7]
            fa[k,5]<-substring(result[j,2],gene[i,4]-pace,nchar(result[j,2])) 
            k=k+1
          }
        }
      }
    }
    output<-data.frame()
    l=1
    for(i in 1:dim(fa)[1]){
      output[l,1]<-paste(">",fa[i,1],";",sep="")
      output[l+1,1]<-fa[i,5]
      l=l+2
    }
    #View(output)
    #return(fa)
    outputFileName<-paste("extend_",times*pace,"bp",".fasta",sep = "")
    write.table(output,outputFileName,sep = "",quote = FALSE,
                col.names = FALSE,row.names = FALSE,na = "")
  }
  fasta<-seqFetch(gff3,p,t)
}
