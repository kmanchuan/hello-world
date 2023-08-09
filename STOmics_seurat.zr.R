########
#
#   Stereo_seurat.R
#   
#   1. load stereo-seq format matrix and convert to seurat object
#   2. add code to do seurat process,
#
########

### Get the parameters
parser = argparse::ArgumentParser(description = 'Script for converting Stereo-seq matrix to seurat format')
parser$add_argument('-i', '--input', dest = 'input', help = 'input tsv filename')
parser$add_argument('-b', '--binsize', dest = 'binsize', default = 1, type = 'integer', help = 'bin size to binning, your input should be in binSize1 if you set this')
parser$add_argument('-s', '--sample', dest = 'sample', help = 'sample ID, will be used as output prefix and seurat object ident')
parser$add_argument('-t', '--tissue', dest = 'tissue', help = 'csv format file listed cell ID which can be used to lasso, should in x_y format')
parser$add_argument('-o', '--out', dest = 'outdir', help = 'directory where to save the output files, all output files will be indexed by sample ID')

parser$add_argument('-r', '--rds', dest = 'rds', help = 'save the pre-seurat matrix in rds format')
parser$add_argument('-5', '--h5ad', dest = 'h5ad', help = 'save the pre-seurat matrix in h5ad format')

parser$add_argument('--minCount', dest = 'minCount', default = 0, type = 'integer', help = 'minimum UMI number')
parser$add_argument('--maxCount', dest = 'maxCount', type = 'integer', help = 'maximum UMI number')
parser$add_argument('--minFeature', dest = 'minFeature', default = 0, type = 'integer', help = 'minimum Feature number')
parser$add_argument('--maxFeature', dest = 'maxFeature', type = 'integer', help = 'maximum Feature number')
parser$add_argument('--vg', dest = 'vg', default = 3000, type = 'integer', help = 'number of variable genes, default 3000')
parser$add_argument('--pc', dest = 'pc', default = 30, type = 'integer', help = 'number of PC to use, default 30')
parser$add_argument('--resolution', dest = 'resolution', default = 0.5, help = 'cluster resolution, default 0.5')

parser$add_argument('--pointSize', dest = 'pointSize', default = 0.2, help = 'point size of spatial plot, default 0.2')
parser$add_argument('--colors', dest = 'colors', default = 70, type = 'integer', help = 'colors palette, one of c(25, 70), default 70')
opts = parser$parse_args()

library(Seurat)
library(SeuratObject)
library(SeuratDisk)
library(ggplot2)
library(patchwork)
library(dplyr)
library(data.table)
library(Matrix)
library(rjson)
library(RColorBrewer)

#opts$sample <- paste0('sample_', gsub('-', '_', opts$sample))
opts$sample <- paste0('sample_', opts$sample)
#paste0连接函数，默认sep=""(空）
#gsub(aa,bb,xx)替换字符串函数，将xx中的bb替换成aa

opts$pointSize <- as.numeric(opts$pointSize)
opts$resolution <- as.numeric(opts$resolution)
#as.numeric数值型转换函数

dir.create(opts$outdir, recursive=TRUE)

data <- fread(file = opts$input)
#fread函数，读取行列统一的大文件

#' group counts into bins
data$x <- trunc(data$x / opts$binsize) * opts$binsize
data$y <- trunc(data$y / opts$binsize) * opts$binsize
#trunc取整函数，只取整数部分
#将横纵坐标根据binsize取整，比如，如果binszie=10,x=2401,y=4522,
#那么处理后，x=2400,y=4520
gc()
if ('MIDCounts' %in% colnames(data)) {
# aa %in% bb，判断aa是否存在于bb中
    data <- data[, .(counts=sum(MIDCounts)), by = .(geneID, x, y)]
} else if ('MIDCount' %in% colnames(data)) {
    data <- data[, .(counts=sum(MIDCount)), by = .(geneID, x, y)]
} else {
    data <- data[, .(counts=sum(UMICount)), by = .(geneID, x, y)]
}

#' create sparse matrix from stereo
gc()
data$cell <- paste0(opts$sample, ':', data$x, '_', data$y)
data$geneIdx <- match(data$geneID, unique(data$geneID))
data$cellIdx <- match(data$cell, unique(data$cell))
#match(aa,bb)返回bb在aa中的位置

gc()
if (! is.null(opts$binsize)){
#is.null(aa) 判断aa是否为空，如果空，返回Ture；
#! is.null(aa)，相反地，判断aa是否非空，如果非空，返回True；
    #write.table(data, file = paste0(opts$outdir, '/', opts$sample, '_bin', opts$binsize, '.tsv'), 
    #            quote = FALSE, sep = '\t', row.names = FALSE )
}
#write.table(aa,file=xx,sep ="\t",quote=FALSE,row.names =TRUE, col.names =TRUE)
#把内容aa输出到文件xx中
#
print ("this sparseMatrix")
gc()
mat <- sparseMatrix(i = data$geneIdx, j = data$cellIdx, x = data$counts, 
                    dimnames = list(unique(data$geneID), unique(data$cell)))
#sparseMatrix稀疏矩阵函数
print ("this cell_coords")
gc()
cell_coords <- unique(data[, c('cell', 'x', 'y')])
#unique去除重复函数，删除cell,x和y都一样的行
print ("this rownames")
gc()
rownames(cell_coords) <- cell_coords$cell
#rownames 重命名行名

#cell_coords$cell <- NULL
print ("this seurat_spatialObj")
gc()
seurat_spatialObj <- CreateSeuratObject(counts = mat, project = 'Stereo', assay = 'Spatial', 
                                        names.delim = ':', meta.data = cell_coords)

#' create pseudo image
cell_coords$x <- cell_coords$x - min(cell_coords$x) + 1
cell_coords$y <- cell_coords$y - min(cell_coords$y) + 1

tissue_lowres_image <- matrix(1, max(cell_coords$y), max(cell_coords$x))
#matrix(aa,x,y)以aa为输入向量，创建一个x行y列的矩阵
#构造一个seruat image

tissue_positions_list <- data.frame(row.names = cell_coords$cell,
                                    tissue = 1,
                                    row = cell_coords$y, col = cell_coords$x,
                                    imagerow = cell_coords$y, imagecol = cell_coords$x)


scalefactors_json <- toJSON(list(fiducial_diameter_fullres = opts$binsize,
                                 tissue_hires_scalef = 1,
                                 tissue_lowres_scalef = 1))
#toJSON: 把json格式 转换成 list格式


#' function to create image object
generate_spatialObj <- function(image, scale.factors, tissue.positions, filter.matrix = TRUE){
    if (filter.matrix) {
        tissue.positions <- tissue.positions[which(tissue.positions$tissue == 1), , drop = FALSE]
    }

    unnormalized.radius <- scale.factors$fiducial_diameter_fullres * scale.factors$tissue_lowres_scalef

    spot.radius <- unnormalized.radius / max(dim(x = image))

    return(new(Class = 'VisiumV1', 
               image = image, 
               scale.factors = scalefactors(spot = scale.factors$tissue_hires_scalef, 
                                            fiducial = scale.factors$fiducial_diameter_fullres, 
                                            hires = scale.factors$tissue_hires_scalef, 
                                            lowres = scale.factors$tissue_lowres_scalef), 
               coordinates = tissue.positions, 
               spot.radius = spot.radius))
}

spatialObj <- generate_spatialObj(image = tissue_lowres_image, 
                                  scale.factors = fromJSON(scalefactors_json), 
                                  tissue.positions = tissue_positions_list)
#可以理解为构建一个spatial背景

#' import image into seurat object
spatialObj <- spatialObj[Cells(x = seurat_spatialObj)]
DefaultAssay(spatialObj) <- 'Spatial'

seurat_spatialObj[['slice1']] <- spatialObj

#' filter out empty cell
seurat_spatialObj <- subset(seurat_spatialObj, subset = nCount_Spatial > 0)
#subset(aa, aa$xx>=2)从aa中获取xx列大于2的行

############# here we finished the creatation of seurat spatial object using pseudo image

#' save seurat spatial object
if (! is.null(opts$hd5)) {
    genes <- as.data.frame(rownames(seurat_spatialObj), row.names = rownames(seurat_spatialObj))
    names(genes) <- 'gene'
    
    cells <- as.data.frame(colnames(seurat_spatialObj), row.names = colnames(seurat_spatialObj))
    names(cells) <- 'cell'

    row <- seurat_spatialObj@images$slice1@coordinates$row
    col <- seurat_spatialObj@images$slice1@coordinates$col
    coordinates <- list(matrix(c(row, col), ncol=2)); names(coordinates) <- 'spatial'
    
    adata <- anndata::AnnData(X = seurat_spatialObj@assays$Spatial@counts, obs = genes, var = cells, varm = coordinates)
    adata <- adata$T
    
    adata$write_h5ad(filename = opts$h5ad)
}

if (! is.null(opts$rds)) {
    saveRDS(seurat_spatialObj, file = opts$rds)
}

######## below prepare for seurat processing
#' function to label cells based on lasso list, need further test
Lasso <- function(object, filename, sample){
    df <- read.csv(filename, header = FALSE, col.names = c('barcode', 'sample'), skip = 1, quote = '')
    df$barcode <- paste(sample, df$barcode, sep = '_')
    object$location <- ifelse(colnames(object) %in% df$barcode, 'tissue', 'nontissue')
    return(object)
}

#' create heatmap and cluster color palette for figure plot function
heatmap_Palette <- colorRampPalette(rev(brewer.pal(11, 'Spectral')))
#if (opts$colors == 25){
#    cluster_Palette <- c('dodgerblue2', '#E31A1C', 'green4', '#6A3D9A', '#FF7F00', 'black', 'gold1', 
#                         'skyblue2', '#FB9A99', 'palegreen2', '#CAB2D6', '#FDBF6F', 'gray70', 'khaki2', 
#                         'maroon', 'orchid1', 'deeppink1', 'blue1', 'steelblue4', 'darkturquoise', 
#                         'green1', 'yellow4', 'yellow3','darkorange4', 'brow')
#} else if (opts$colors == 70){
#    qual_col_pals <- brewer.pal.info[brewer.pal.info$category == 'qual',]
#    cluster_Palette <- unique(unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals))))
#}
#cluster_Palette <- c('#33a02c','#ff69b4','#1f78b4','#ff4500','#b2df8a','#a6cee3','#a020f0','#ff7f00','#66cdaa','#db7093','#1e90ff','#fdbf6f','#6a3d9a','#ffd700','#b15928','#8dd3c7','#ffffb3','#bebada','#ffc1c1','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5','#ffed6f','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf','#999999','#fbb4ae','#b3cde3','#ccebc5','#decbe4','#fed9a6','#ffffcc','#e5d8bd','#fddaec','#f2f2f2','#8c510a','#bf812d','#dfc27d','#f6e8c3','#f5f5f5','#c7eae5','#80cdc1','#35978f','#01665e')

cluster_Palette <- c('#33a02c','#ff69b4','#1f78b4','#ff4500','#b2df8a','#a6cee3','#a020f0','#ff7f00','#66cdaa','#db7093','#1e90ff','#fdbf6f','#6a3d9a','#ffd700','#b15928','#8dd3c7','#ffffb3','#bebada','#ffc1c1','#80b1d3','#fdb462','#b3de69','#fccde5','#d9d9d9','#bc80bd','#ccebc5','#ffed6f','#4daf4a','#984ea3','#ff7f00','#ffff33','#a65628','#f781bf','#999999','#fbb4ae','#b3cde3','#ccebc5','#decbe4','#fed9a6','#ffffcc','#e5d8bd','#fddaec','#f2f2f2','#8c510a','#bf812d','#dfc27d','#f6e8c3','#f5f5f5','#c7eae5','#80cdc1','#35978f','#01665e','#B84D64', '#864A68', '#EE7072', '#E32D32', '#998B95', '#5E549A', '#8952A0', '#4552A0', '#384B97', '#2B3B72', '#911310', '#384C99', '#9B8E8C', '#7CA878', '#35A132', '#6B70B0', '#3D6AAA', '#394D9B', '#75ACC3', '#20ACBD', '#38509F', '#959897','#FF0000','#00FF00','#0000FF','#FFFF00','#FF00FF','#00FFFF','#800080','#008080','#808000','#FFC0CB','#FFA500','#FF1493','#00FF7F','#7FFFD4','#B8860B','#8B008B','#CD5C5C','#9ACD32','#7B68EE','#FFA07A','#F08080','#FFFFE0','#8B4513','#2E8B57','#FF8C00','#9370DB')


#' custom spatial plot on specific feature
iPlot <- function(object, features, pt.size = 0.2){
    plot <- ggplot(object@meta.data, aes_string(x = 'x', y = 'y', color = features)) +
            geom_point(shape = 19, size = pt.size) + theme_void() +
            theme(axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), 
#axis.text: 坐标轴文本； axis.ticks: 坐标轴刻度线；  panel.grid: 网格线 
                  axis.title = element_blank(), axis.line = element_blank(), plot.margin=margin(t=1,b=1,r=1,l=1))  
    if (features %in% c('nCount_Spatial', 'nFeature_Spatial')){
        plot <- plot + scale_color_gradientn(colours = heatmap_Palette(100)) + 
                guides(colour = guide_legend(override.aes = list(size=3), nrow = 10))
    }else if(features %in% c('seurat_clusters')){
        plot <- plot + scale_color_manual(values = cluster_Palette) +
                guides(colour = guide_legend(override.aes = list(size=3), nrow = 10))
    }
    plot <- plot + theme_void() + coord_fixed()
    return(plot)
}

#' basic statistics plot function
Preview <- function(object, feature){
    p1 <- VlnPlot(object, features = feature, pt.size = 0.2) + NoLegend() + theme(axis.text.x=element_blank(), axis.title.x=element_blank()) 
    #p2 <- SpatialFeaturePlot(obj, pt.size.factor = opts$pointSize, features = feature, stroke = 0) + theme(legend.position = 'right')
    p2 <- iPlot(object, features = feature, pt.size = opts$pointSize) 
    patch <- p1 | p2
    return(patch)
}

#' seurat clustering workflow
Clustering <- function(object, dims = 30){
    object <- RunPCA(object)
    object <- FindNeighbors(object, dims = 1:dims)
    object <- FindClusters(object, verbose = FALSE, resolution = opts$resolution)
    object <- RunUMAP(object, dims = 1:dims)
    return(object)
}

############ start seurat processing

obj <- seurat_spatialObj 

#' filter object based on nCount and nFeature
if (is.null(opts$maxCount)){
    opts$maxCount <- max(obj$nCount_Spatial)
}
if (is.null(opts$maxFeature)){
    opts$maxFeature <- max(obj$nFeature_Spatial)
}
obj <- subset(obj, subset = nCount_Spatial >= opts$minCount & nCount_Spatial <= opts$maxCount & 
              nFeature_Spatial >= opts$minFeature & nFeature_Spatial <= opts$maxFeature)

#' lasso object based on tissue parameter, need further test
if (!is.null(opts$tissue)){
    obj <- Lasso(obj, opts$tissue, opts$sample)
}

#' basic statistics plot on nCount and nFeature
plot1 <- Preview(obj, 'nCount_Spatial')
plot2 <- Preview(obj, 'nFeature_Spatial')

#' do normalization and clustering
gc()
print ("this SCTransform")
obj <- SCTransform(obj, assay = 'Spatial', variable.features.n = as.numeric(opts$vg), 
                   return.only.var.genes = FALSE, n_genes=NULL, min_cells=5, method='qpoisson')
gc()
print ("this Clustering")
obj <- Clustering(obj, dims = opts$pc)


#' find all markers of each seurat_cluster and save to file
markers <- FindAllMarkers(obj, min.pct = 0.1, logfc.threshold = 0.25)
markers <- markers[order(markers$cluster,markers$avg_log2FC,decreasing=TRUE),]
write.table(format(markers,digits=3), paste0(opts$outdir, '/', opts$sample, '_bin', opts$binsize, '_AllMarkers.xls'), sep = '\t', quote = FALSE, row.names =FALSE)

#' if lasso argment is specified, find markers between lasso cells and others, need further test
if ('location' %in% attributes(obj)$name){
    markers <- FindMarkers(obj, ident.1 = 'tissue', ident.2 = 'nontissue', group.by = 'location')
    write.table(markers, paste0(opts$outdir, '/', opts$sample, '_bin', opts$binsize, '_TissueMarkers.xls'), sep='\t', quote = FALSE, row.names =FALSE)
}

#obj$x<-sapply(rownames(obj@meta.data),function(x)as.numeric(unlist(strsplit(x,"_"))[2]))
#obj$y<-sapply(rownames(obj@meta.data),function(x)as.numeric(unlist(strsplit(x,"_"))[3]))

h <- (max(obj@meta.data$y)-min(obj@meta.data$y))%/%opts$binsize
w <- (max(obj@meta.data$x)-min(obj@meta.data$x))%/%opts$binsize
s <- w/h
ps <- 3000/(max(h,w)^1.6)
ps_dim <- 5000/(max(h,w)^2)

print('---point_size using---')
print(ps)

#' plot cluster results on umap and spatial
plot3 <- DimPlot(obj, reduction="umap", cols = cluster_Palette, pt.size=ps_dim, label = TRUE, label.size=7) + theme(plot.title=element_blank(),axis.title.x=element_text(size=15),axis.title.y=element_text(size=15),axis.text.x=element_text(vjust=1,size=15),axis.text.y=element_text(vjust=1,size=15), plot.margin=unit(c(1,1,1,1),"cm")) + guides(colour = guide_legend(override.aes = list(size=5), nrow = 15))
plot4 <- iPlot(obj, feature = 'seurat_clusters', pt.size = ps)

pdf(paste0(opts$out, '/', opts$sample,'_bin', opts$binsize, '.pdf'), width = 8, height = 9)
plot1 / plot2 / (plot3 + plot4)
#patchwork::wrap_plots(plot1, plot2, plot3|plot4,ncol=1,widths=c(3,3),heights=c(3,3,4) )
dev.off()
png(paste0(opts$out, '/', opts$sample,'_bin', opts$binsize, '_nstat.png'), width = 700, height = 600)
plot1 / plot2 
dev.off()
pdf(paste0(opts$out, '/', opts$sample,'_bin', opts$binsize, '_nstat.pdf'), width = 9, height = 8)
plot1 / plot2 
dev.off()
png(paste0(opts$out, '/', opts$sample,'_bin', opts$binsize, '_cluster.png'), width = 2*900*s, height=1200*s)
plot3|plot4
dev.off()
#ggsave(file = paste0(opts$out, '/', opts$sample,'_bin', opts$binsize, '_cluster.pdf'), plot = plot3|plot4, width = 2*6*s, height=9)
pdf(paste0(opts$out, '/', opts$sample,'_bin', opts$binsize, '_cluster.pdf'), width = 2*9*s, height=12)
plot3|plot4
dev.off()

#' plot each cluster on spatial
allclusters <- sort(unique(obj$seurat_clusters))
tmplist <- list()
for ( i in 1:length(allclusters)){
        each <- allclusters[[i]]
            plottmp <- obj@meta.data %>% mutate(highlight_flag=ifelse(seurat_clusters==each,i-1,"others")) %>% ggplot(aes(x=x,y=y)) +
            geom_point(shape = 19, size =0.2, aes(color = highlight_flag )) + 
            scale_color_manual(values = c('red','skyblue2')) +
            theme(panel.background=element_blank(), axis.text = element_blank(), axis.ticks = element_blank(), panel.grid = element_blank(), axis.title = element_blank(), axis.line = element_blank(),legend.key = element_blank(), legend.position = 'right', legend.title = element_blank(), plot.margin=unit(c(3,0,3,0),"pt"))  +  coord_fixed() +
            guides(colour = guide_legend(override.aes = list(size=3), nrow = 2))
            tmplist[[i]] <- plottmp
}
ggsave(paste0(opts$out, '/', opts$sample,'_bin', opts$binsize, "_each_cluster.pdf"),patchwork::wrap_plots(tmplist,ncol = 4),width = 12, height = 12, scale =1)
ggsave(paste0(opts$out, '/', opts$sample,'_bin', opts$binsize, "_each_cluster.png"),patchwork::wrap_plots(tmplist,ncol = 4),width = 12, height = 12, scale =1)

#ggpubr::ggarrange(plotlist = tmplist, ncol = 3, byrow = TRUE)
#pdf(paste0(opts$out, '/', opts$sample,'_bin', opts$binsize, "_each_cluster.pdf"),width=12)
#patchwork::wrap_plots(tmplist,ncol = 3)
#dev.off()
#png(paste0(opts$out, '/', opts$sample,'_bin', opts$binsize, "_each_cluster.png"), width = 900)
#patchwork::wrap_plots(tmplist,ncol = 3)
#dev.off()

#' do heatmap for top cluster markers
#topmkg <- markers %>% group_by(cluster) %>% slice(1:3)
#%>% 是data.table包中的运算符,类似于管道符：markers 中 cluster这列，每个group，前3个值
top3mk <- format(markers,digits=3) %>% group_by(cluster) %>% slice(1:3)
write.table(top3mk, paste0(opts$outdir, '/', opts$sample, '_bin', opts$binsize, '_top3_Markers.xls'), sep = '\t', quote = FALSE, row.names =FALSE)
top1mk <- format(markers,digtis=3) %>% group_by(cluster) %>% slice(1:1)
write.table(top1mk, paste0(opts$outdir, '/', opts$sample, '_bin', opts$binsize, '_top1_Markers.xls'), sep = '\t', quote = FALSE, row.names =FALSE)
plot5 <- DoHeatmap(obj, features = unique(top3mk$gene),size = 2) + NoLegend()
pdf(paste0(opts$out, '/', opts$sample, '_bin', opts$binsize,'_heatmapClusters.pdf'),width = 9, height = 9)
plot5
dev.off()
png(paste0(opts$out, '/', opts$sample, '_bin', opts$binsize,'_heatmapClusters.png'), width = 900, height = 700)
plot5
dev.off()


#' do vlnplot for top cluster markers
plot6 <- VlnPlot(obj, idents=sort(unique(obj@meta.data$seurat_clusters)), features=unique(top3mk$gene), stack=TRUE, pt.size=0) + guides(fill="none")
pdf(paste0(opts$out, '/', opts$sample, '_bin', opts$binsize,'_vlnplotClusters.pdf'),width = 12, height = 9)
plot6
dev.off()
png(paste0(opts$out, '/', opts$sample, '_bin', opts$binsize,'_vlnplotClusters.png'),width = 900, height = 700)
plot6
dev.off()

#' do dotplot for top cluster markers
#plot7 <- DotPlot(obj, idents=sort(unique(obj@meta.data$seurat_clusters)), features=unique(top3mk$gene))+ RotatedAxis() + scale_x_discrete("") + scale_y_discrete("")
plot7 <- DotPlot(obj, idents=sort(unique(obj@meta.data$seurat_clusters)), features=unique(top3mk$gene))+ scale_x_discrete("") + scale_y_discrete("") + theme(axis.text.x = element_text(angle = 90,hjust = -0.1,vjust = 0.8))
pdf(paste0(opts$out, '/', opts$sample, '_bin', opts$binsize,'_dotplotClusters.pdf'),width = 12, height = 9)
plot7
dev.off()
png(paste0(opts$out, '/', opts$sample, '_bin', opts$binsize,'_dotplotClusters.png'), width = 900, height = 700)
plot7
dev.off()

#' save processed seurat object to rds file
#saveRDS(obj, file = paste0(opts$outdir, '/', opts$sample, '_bin', opts$binsize, '_seurat.rds'))
#' save processed seurat object to loom file
#obj.loom <- as.loom(obj, filename = paste0(opts$outdir, '/', opts$sample, '_bin', opts$binsize, '_seurat.loom'), 
#        overwrite = TRUE)
#obj.loom$close_all()
print ("task finished")
