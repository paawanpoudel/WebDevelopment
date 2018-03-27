#load('/usr/share/blueprint/R_commands/gene_and_tx_summary.rda')

#/usr/share/blueprint/R_commands/gene_and_tx_summary.rda

create_heatmap_gene <- function() {

	#load('/home/pp376/RServe/datadir/DE_genes_08012014.rda')
       
	cells.order <- c("HSC", "MPP", "CLP", "CMP", "GMP", "MEP", "EB", "MK")

	# matching needs to be done
	index <- (de.genes$ensembl_gene_id %in% gene.sel | de.genes$hgnc_symbol %in% gene.sel)
	full.data.sel=de.genes[index,]

	data2plot = exp(as.matrix(full.data.sel[,grep("^mu",colnames(full.data.sel))]))
	colnames(data2plot) <- sub("mu_([[:alpha:]]+)_([[:digit:]])","\\1", colnames(data2plot))
	data2plot = data2plot[, order(match(colnames(data2plot), cells.order))]

	data2plot = (data2plot - apply(data2plot, 1, mean)) /apply(data2plot, 1, sd)
	heatmap.2(data2plot,
    	labCol=colnames(data2plot),
	col=colorRampPalette(brewer.pal(9,"RdYlBu"),bias=2)(100),mar=c(6,12),
    	Colv=F,
    	Rowv=T,
    	trace="none",
    	density.info="none",
    	dendrogram="none",
    	ColSideColors=rep(brewer.pal(8,"Paired"),c(3,3,3,3,3,4,3,3)),
    	distfun=function(x) as.dist(1-cor(t(x),method="pear")),
    	hclustfun=function(x) 
   	hclust(x,method="ward"), labRow=full.data.sel$external_gene_id, lhei=c(2,8),colsep=c(seq(3,15,3),19,22),keysize=1, cexRow=1.5)

	dev.off()


}
create_heatmap_transcript <- function() {

	#load('/home/pp376/RServe/datadir/DE_genes_08012014.rda')

	cells.order <- c("HSC", "MPP", "CLP", "CMP", "GMP", "MEP", "EB", "MK")

	# matching needs to be done
	index <- (de.tx$ensembl_transcript_id %in% gene.sel)
	full.data.sel=de.tx[index,]

	data2plot = exp(as.matrix(full.data.sel[,grep("^mu",colnames(full.data.sel))]))
	colnames(data2plot) <- sub("mu_([[:alpha:]]+)_([[:digit:]])","\\1", colnames(data2plot))
	data2plot = data2plot[, order(match(colnames(data2plot), cells.order))]

	data2plot = (data2plot - apply(data2plot, 1, mean)) /apply(data2plot, 1, sd)
	heatmap.2(data2plot,
    	labCol=colnames(data2plot),
	col=colorRampPalette(brewer.pal(9,"RdYlBu"),bias=2)(100),mar=c(6,15),
    	Colv=F,
    	Rowv=T,
    	trace="none",
    	density.info="none",
    	dendrogram="none",
    	ColSideColors=rep(brewer.pal(8,"Paired"),c(3,3,3,3,3,4,3,3)),
    	distfun=function(x) as.dist(1-cor(t(x),method="pear")),
    	hclustfun=function(x) 
   	hclust(x,method="ward"), labRow=full.data.sel$ensembl_transcript_id, lhei=c(2,8),colsep=c(seq(3,15,3),19,22),keysize=1, cexRow=1.5)

	dev.off()


}


riverplot_gene <-function() {

options(stringsAsFactors=F)

require(ggplot2)
require(reshape)
require(RColorBrewer)

createData2Plot <- function(genes, multiplier = 10, white.line.points = c(), white.space = 0.1)
{
    cells.order <- c("HSC", "MPP", "CLP", "CMP", "GMP", "MEP", "EB", "MK")
    
    cell.type.expression <- genes
    cell.type.expression <- merge(cell.type.expression, de.genes, by.x="gene.name", by.y="hgnc_symbol")
    
    cell.type.expression <- unique(cell.type.expression[,grep("gene.name|col|ensembl_gene_id|^mean",colnames(cell.type.expression))])
    colnames(cell.type.expression) <- sub("mean_([[:alpha:]]+)","\\1", colnames(cell.type.expression))
    
    cell.type.expression <- melt(cell.type.expression, variable_name = "cell.type", id=c("gene.name", "ensembl_gene_id", "col"))
    
## This is where I normalise the expression per proportion

cell.type.expression$value <- log(exp(as.numeric(cell.type.expression$value)) + 1, 2)
    
cell.type.expression$cell.type <- factor(cell.type.expression$cell.type, levels=cells.order)
print(cell.type.expression)
 
cell.type.expression$gene.name <- factor(cell.type.expression$gene.name, levels = genes$gene.name)
cell.type.expression <- cell.type.expression[order(cell.type.expression$gene.name), ]
    
cum.sum <- data.frame()

for (ct in levels(cell.type.expression$cell.type))
{
max.fpkm.sum <- cumsum(rev(cell.type.expression$value[cell.type.expression$cell.type==ct]))
min.fpkm.sum <- c(0, max.fpkm.sum[1:(length(max.fpkm.sum)-1)])
        
### If I want to have space between all of them
if (length(white.line.points) == 0)
{
white.line <- c(0, cumsum(rep(white.space, (length(max.fpkm.sum)-1))))
max.fpkm.sum <- max.fpkm.sum + white.line
min.fpkm.sum <- min.fpkm.sum + white.line
} else
{
for (ind in cumsum(rev(white.line.points[-1])))
{
max.fpkm.sum[(ind + 1):length(max.fpkm.sum)] <- max.fpkm.sum[(ind + 1):length(max.fpkm.sum)] + white.space
min.fpkm.sum[(ind + 1):length(min.fpkm.sum)] <- min.fpkm.sum[(ind + 1):length(min.fpkm.sum)] + white.space
}
}
        if(length(levels(cell.type.expression$gene.name))>1)
        {
            cum.sum <- rbind(cum.sum, cbind(gene.name = rev(levels(cell.type.expression$gene.name)), cell.type = rep(ct, length(levels(cell.type.expression$gene.name))), cum.fpkm = max.fpkm.sum, min.fpkm = min.fpkm.sum))
        } else
        {
            cum.sum <- rbind(cum.sum, cbind(gene.name = rev(levels(cell.type.expression$gene.name)), cell.type = rep(ct, length(levels(cell.type.expression$gene.name))), cum.fpkm = max.fpkm.sum, min.fpkm = 0))
        }
}
    
cell.type.expression <- merge(cell.type.expression, cum.sum, by = c("cell.type", "gene.name"))
cell.type.expression <- cell.type.expression[order(cell.type.expression$gene.name), ]
cell.type.expression$cum.fpkm <- multiplier*(as.numeric(cell.type.expression$cum.fpkm))
cell.type.expression$min.fpkm <- multiplier*(as.numeric(cell.type.expression$min.fpkm))
    
return(cell.type.expression)
}

plotGeneListsinTreeStruture <- function(cell.type.expression)
{
cell.type.coordinates <- data.frame(cbind(cell.type = c("HSC", "MPP", "MPP", "CLP", "CMP", "CMP", "GMP", "MEP", "MEP", "EB", "MK"), branch = c("MEGA", "MEGA", "LYMPHO", "LYMPHO", "MEGA", "GRAN-MONO", "GRAN-MONO", "ERY", "MEGA", "ERY", "MEGA"), x = c(0, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4), y = c(0, 0, 0, 40, -40, -40, -10, -70, -70, -50, -120)))
    
    
text <- data.frame(cbind(labels = c("HSC", "MPP", "CLP", "CMP", "GMP", "MEP", "EB", "MK"), x = c(0, 0.9, 2, 1.9, 3, 2.9, 4, 4), y= c(-10, -10, 27, -50, 10, -80, -5, -130)))
text$x <- as.numeric(text$x)
text$y <- as.numeric(text$y)
for (ct in unique(cell.type.expression$cell.type))
{
mean.norm <- (max(cell.type.expression$cum.fpkm[cell.type.expression$cell.type==ct])/2)
print(paste("For cell type ", ct, "substract: ", mean.norm, sep=""))
cell.type.expression$cum.fpkm[cell.type.expression$cell.type==ct] <- cell.type.expression$cum.fpkm[cell.type.expression$cell.type==ct] - mean.norm
cell.type.expression$min.fpkm[cell.type.expression$cell.type==ct] <- cell.type.expression$min.fpkm[cell.type.expression$cell.type==ct] - mean.norm
text$y[text$labels==ct] <- text$y[text$labels==ct] - mean.norm
}
    
data2plot <- merge(cell.type.coordinates, cell.type.expression, by = "cell.type")
for(clmn in c("x", "y", "value", "cum.fpkm")) data2plot[, clmn] <- as.numeric(data2plot[, clmn])
data2plot <- data2plot[order(data2plot$gene.name), ]
    
q <- ggplot(data=data2plot, aes(x=x, group=gene.name))
for (i in unique(data2plot$branch))
{
q <- q + geom_ribbon(data = data2plot[data2plot$branch==i, ], aes(x = x, ymin= y + min.fpkm, ymax= y + cum.fpkm, fill = gene.name))
}
q <- q + theme_bw() + labs(x="", y="") + scale_x_continuous(breaks=NULL, labels=NULL) + scale_y_continuous(breaks=NULL, labels=NULL) + scale_fill_manual(name="", limits = data2plot$gene.name, values = data2plot$col)
    
q <- q + guides(fill = guide_legend(nrow = 2))
q <- q + annotate('text', x = text$x, y = text$y, label =  text$labels,  col = 'black', size=7) + theme(legend.position='none', panel.border = element_rect(colour = "white")) + theme(legend.position='bottom', legend.text = element_text(size = 14, face = "bold"))
    
return(q)
}


### need to input this from the play

#genes <- (de.genes$ensembl_gene_id %in% gene.sel | de.genes$hgnc_symbol %in% gene.sel)

#genes=genes.sel
index <- (de.genes$ensembl_gene_id %in% gene.sel | de.genes$hgnc_symbol %in% gene.sel)
genes=de.genes$hgnc_symbol[index]
#genes <- c("ITGA2B", "RUNX1", "MYB", "SPI1", "GATA1", "MPEG1", "TAL1", "MPO", "LYZ")

trans.colours <- rev(brewer.pal(ifelse(length(genes)>3, length(genes), 3), "Paired"))
trans.colours <- trans.colours[1:length(genes)]

trans <- data.frame(gene.name=genes, col=trans.colours)

cell.type.expression <- createData2Plot(trans, multiplier = 2, white.space = 0)

plot <- plotGeneListsinTreeStruture(cell.type.expression)
 
print(plot)
dev.off()

}


riverplot_transcript <-function() {

options(stringsAsFactors=F)

require(ggplot2)
require(reshape)
require(RColorBrewer)

createData2Plot_t <- function(genes, multiplier = 10, white.line.points = c(), white.space = 0.1)
    {
        cells.order <- c("HSC", "MPP", "CLP", "CMP", "GMP", "MEP", "EB", "MK")
        
        cell.type.expression <- genes
        cell.type.expression <- merge(cell.type.expression, de.tx, by.x = "gene.name", by.y = "ensembl_transcript_id")
        
        cell.type.expression <- unique(cell.type.expression[,grep("gene.name|col|ensembl_gene_id|hgnc_symbol|^mean",colnames(cell.type.expression))])
        colnames(cell.type.expression) <- sub("mean_([[:alpha:]]+)","\\1", colnames(cell.type.expression))
        
        cell.type.expression <- melt(cell.type.expression, variable_name = "cell.type", id=c("gene.name", "ensembl_gene_id", "hgnc_symbol", "col"))
        
        ## This is where I normalise the expression per proportion
        
        cell.type.expression$value <- log(exp(as.numeric(cell.type.expression$value)) + 1, 2)
        
        cell.type.expression$cell.type <- factor(cell.type.expression$cell.type, levels=cells.order)
        print(cell.type.expression)
        
        cell.type.expression$gene.name <- factor(cell.type.expression$gene.name, levels = genes$gene.name)
        cell.type.expression <- cell.type.expression[order(cell.type.expression$gene.name), ]
        
        cum.sum <- data.frame()
        
        for (ct in levels(cell.type.expression$cell.type))
        {
            max.fpkm.sum <- cumsum(rev(cell.type.expression$value[cell.type.expression$cell.type==ct]))
            min.fpkm.sum <- c(0, max.fpkm.sum[1:(length(max.fpkm.sum)-1)])
            
            ### If I want to have space between all of them
            if (length(white.line.points) == 0)
            {
                white.line <- c(0, cumsum(rep(white.space, (length(max.fpkm.sum)-1))))
                max.fpkm.sum <- max.fpkm.sum + white.line
                min.fpkm.sum <- min.fpkm.sum + white.line
            } else
            {
                for (ind in cumsum(rev(white.line.points[-1])))
                {
                    max.fpkm.sum[(ind + 1):length(max.fpkm.sum)] <- max.fpkm.sum[(ind + 1):length(max.fpkm.sum)] + white.space
                    min.fpkm.sum[(ind + 1):length(min.fpkm.sum)] <- min.fpkm.sum[(ind + 1):length(min.fpkm.sum)] + white.space
                }
            }
            if(length(levels(cell.type.expression$gene.name))>1)
            {
                cum.sum <- rbind(cum.sum, cbind(gene.name = rev(levels(cell.type.expression$gene.name)), cell.type = rep(ct, length(levels(cell.type.expression$gene.name))), cum.fpkm = max.fpkm.sum, min.fpkm = min.fpkm.sum))
            } else
            {
                cum.sum <- rbind(cum.sum, cbind(gene.name = rev(levels(cell.type.expression$gene.name)), cell.type = rep(ct, length(levels(cell.type.expression$gene.name))), cum.fpkm = max.fpkm.sum, min.fpkm = 0))
            }
        }
        
        cell.type.expression <- merge(cell.type.expression, cum.sum, by = c("cell.type", "gene.name"))
        cell.type.expression <- cell.type.expression[order(cell.type.expression$gene.name), ]
        cell.type.expression$cum.fpkm <- multiplier*(as.numeric(cell.type.expression$cum.fpkm))
        cell.type.expression$min.fpkm <- multiplier*(as.numeric(cell.type.expression$min.fpkm))
        
        return(cell.type.expression)
    }
    
plotGeneListsinTreeStruture_t <- function(cell.type.expression)
    {
        cell.type.coordinates <- data.frame(cbind(cell.type = c("HSC", "MPP", "MPP", "CLP", "CMP", "CMP", "GMP", "MEP", "MEP", "EB", "MK"), branch = c("MEGA", "MEGA", "LYMPHO", "LYMPHO", "MEGA", "GRAN-MONO", "GRAN-MONO", "ERY", "MEGA", "ERY", "MEGA"), x = c(0, 1, 1, 2, 2, 2, 3, 3, 3, 4, 4), y = c(0, 0, 0, 40, -40, -40, -10, -70, -70, -50, -120)))
        
        
        text <- data.frame(cbind(labels = c("HSC", "MPP", "CLP", "CMP", "GMP", "MEP", "EB", "MK"), x = c(0, 0.9, 2, 1.9, 3, 2.9, 4, 4), y= c(-10, -10, 27, -50, 10, -80, -5, -130)))
        text$x <- as.numeric(text$x)
        text$y <- as.numeric(text$y)
        for (ct in unique(cell.type.expression$cell.type))
        {
            mean.norm <- (max(cell.type.expression$cum.fpkm[cell.type.expression$cell.type==ct])/2)
            print(paste("For cell type ", ct, "substract: ", mean.norm, sep=""))
            cell.type.expression$cum.fpkm[cell.type.expression$cell.type==ct] <- cell.type.expression$cum.fpkm[cell.type.expression$cell.type==ct] - mean.norm
            cell.type.expression$min.fpkm[cell.type.expression$cell.type==ct] <- cell.type.expression$min.fpkm[cell.type.expression$cell.type==ct] - mean.norm
            text$y[text$labels==ct] <- text$y[text$labels==ct] - mean.norm
        }
        
        data2plot <- merge(cell.type.coordinates, cell.type.expression, by = "cell.type")
        for(clmn in c("x", "y", "value", "cum.fpkm")) data2plot[, clmn] <- as.numeric(data2plot[, clmn])
        data2plot <- data2plot[order(data2plot$gene.name), ]
        
        q <- ggplot(data=data2plot, aes(x=x, group=gene.name))
        for (i in unique(data2plot$branch))
        {
            q <- q + geom_ribbon(data = data2plot[data2plot$branch==i, ], aes(x = x, ymin= y + min.fpkm, ymax= y + cum.fpkm, fill = gene.name))
        }
        q <- q + theme_bw() + labs(x="", y="") + scale_x_continuous(breaks=NULL, labels=NULL) + scale_y_continuous(breaks=NULL, labels=NULL) + scale_fill_manual(name="", limits = data2plot$gene.name, values = data2plot$col)
        
        q <- q + guides(fill = guide_legend(nrow = ceiling(length(levels(data2plot$gene.name))/4)))
        q <- q + annotate('text', x = text$x, y = text$y, label =  text$labels,  col = 'black', size=7) + theme(legend.position='none', panel.border = element_rect(colour = "white")) + theme(legend.position='bottom', legend.text = element_text(size = 10, face = "bold"))
        
        return(q)
    }

### need to input this from the play

#genes <- (de.genes$ensembl_gene_id %in% gene.sel | de.genes$hgnc_symbol %in% gene.sel)

#genes=genes.sel
index <- (de.tx$ensembl_transcript_id %in% gene.sel)
genes=de.tx$ensembl_transcript_id[index]
#genes <- c("ITGA2B", "RUNX1", "MYB", "SPI1", "GATA1", "MPEG1", "TAL1", "MPO", "LYZ")

trans.colours <- rev(brewer.pal(ifelse(length(genes)>3, length(genes), 3), "Paired"))
trans.colours <- trans.colours[1:length(genes)]

trans <- data.frame(gene.name=genes, col=trans.colours)

cell.type.expression <- createData2Plot_t(trans, multiplier = 2, white.space = 0)

plot <- plotGeneListsinTreeStruture_t(cell.type.expression)
 
print(plot)
dev.off()

}






