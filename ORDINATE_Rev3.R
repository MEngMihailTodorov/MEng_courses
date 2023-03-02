#ACHTUNG: From the RStudio menu, click on "Session" and then "Set Working Directory" to "To Source File Location"
#Script for calculation of beta diversity (Principle Coordinate Analysis/Nonmetric Distance Scaling)

#Sometimes bray-curtis distance doesn't work with the newer version of phyloseq
#install it using:
#library("devtools")
#install_github("joey711/phyloseq")


library(phyloseq)
library(ggrepel)
library(vegan)
library(ggplot2)
library(ape)
library(phangorn)
library(stringr)
library(grid)
library(biomformat)
library(pals)
library(yarrr)

#PARAMETERS ###########################
which_level<-"Otus" #Phylum Class Order Family Genus Otus
kind <- "se" #sd se (sd is for drawing ellipse based on sd, se is for drawing ellipse based on standard errors)
which_method<-"PCOA" #NMDS or PCOA
which_distance<-"bray" #bray unifrac wunifrac userprovided
user_provided_distance<-NULL #If you want to use your own distance matrix then the value should be the file that you'll read from
#user_provided_distance<-read.csv("../../Data/ko_metagenome.dist",sep="\t",header=T,row.names=1,check.names=FALSE)
library(biomformat);b_<-read_biom("../../Data/collated_feature_w_tax.biom");physeq<-merge_phyloseq(otu_table(as(biom_data(b_),"matrix"),taxa_are_rows=TRUE),tax_table(as(observation_metadata(b_),"matrix")))
meta_table<-read.csv("../../Data/meta_data_rev2.csv",header=T,row.names=1)
#Load the tree using ape package
OTU_tree <- read.tree("../../Data/tree.nwk")
exclude_legends=FALSE
draw_mean_values_text=FALSE
point_size=10
point_opacity=1
draw_glow=FALSE
point_glow_opacity=0.1
point_glow_differential=2
draw_confidence_intervals=TRUE
draw_ellipses_and_not_polygons=FALSE
pairwise_connections_and_not_longitudinal_connections=FALSE#FALSE
should_connections_end_in_arrows=TRUE
should_we_reverse_direction_of_arrows=TRUE
opacity_ellipses_polygons=0.45
linesize_ellipses_polygons=1
linetype_ellipses_polygons="dashed" #blank solid dashed dotted dotdash longdash twodash
linking_samples_line_size=1
linking_samples_line_opacity=1.0
linking_samples_linetype="solid" #blank solid dashed dotted dotdash longdash twodash
legend_text_size=22
legend_title_size=24
axis_title_size=34
text_size=22
axis_text_size=30
height_image=40 #15
width_image=60
use_provided_colors=TRUE
#/PARAMETERS ###########################

abund_table<-otu_table(physeq)
abund_table<-t(abund_table)
#Uncomment if you'd like to get rid of samples below a certain library size
abund_table<-abund_table[rowSums(abund_table)>=5000,]


OTU_taxonomy<-data.frame(as(tax_table(physeq),"matrix"))
colnames(OTU_taxonomy)<-c("Kingdom","Phylum","Class","Order","Family","Genus","Otus")

#Ensure that all columns of OTU_taxonomy are character and not factors
OTU_taxonomy[] <- lapply(OTU_taxonomy, function(x) as.character(x))
OTU_taxonomy[is.na(OTU_taxonomy)]<-""
OTU_taxonomy$Otus<-gsub("D_6__|s__","",OTU_taxonomy$Otus)
OTU_taxonomy$Genus<-gsub("D_5__|g__","",OTU_taxonomy$Genus)
OTU_taxonomy$Family<-gsub("D_4__|f__","",OTU_taxonomy$Family)
OTU_taxonomy$Order<-gsub("D_3__|o__","",OTU_taxonomy$Order)
OTU_taxonomy$Class<-gsub("D_2__|c__","",OTU_taxonomy$Class)
OTU_taxonomy$Phylum<-gsub("D_1__|p__","",OTU_taxonomy$Phylum)
OTU_taxonomy$Kingdom<-gsub("D_0__|k__","",OTU_taxonomy$Kingdom)

#Remove singletons and adjust OTU_taxonomy
abund_table<-abund_table[,colSums(abund_table)>1]
OTU_taxonomy<-OTU_taxonomy[colnames(abund_table),]

#get rid of contaminants with "Unassigned", "Chloroplast" and "Mitochondria" assignment", and "non classified" at Phylum level
abund_table<-abund_table[,!(OTU_taxonomy$Kingdom %in% c("Unassigned") | OTU_taxonomy$Phylum=="" | OTU_taxonomy$Order %in% c("Chloroplast") | OTU_taxonomy$Family %in% c("Mitochondria"))]
OTU_taxonomy<-OTU_taxonomy[colnames(abund_table),]

#extract subset of abund_table for which samples also exists in meta_table
abund_table<-abund_table[rownames(abund_table) %in% rownames(meta_table),]
#when reducing the abund_table, there is a high likelihood that an OTU was only present in a sample that is removed, so we shrink
#the abund_table to get rid of empty columns
abund_table<-abund_table[,colSums(abund_table)>0]
#make your meta_table smaller by only considering samples that appear in abund_table

meta_table<-meta_table[rownames(abund_table),]
#make OTU_taxonomy smaller by only considering OTUs that appear in abund_table
OTU_taxonomy<-OTU_taxonomy[colnames(abund_table),]
#At this point we have abund_table, meta_table, and OTU_taxonomy are ready and their dimensions should match
#/DATA IMPORT############################################################


#PARAMETERS CHANGE THE GROUPING COLUMN AS YOU DESIRE############################
#In the hypothesis space, all you need is to select the rows in meta_table you are interested in
#and then allocate a column to meta_table$Groups that you want to use, meta_table$Connections to connect them
#additionally if you provide a second meta_table$Subconnections, it will connect the group averages
#You can use meta_table$Type to assign shape and PERMANOVA_variables to give variables

label="Hypothesis1"
meta_table<-meta_table[meta_table$Material %in% c(
  "Microplastic_fibers",                           
  "PE",                                              
  "PP",                                              
  "PET",                                             
  "Pellets",                                         
  "PS",                                              
  "Plastic_marine_debris_PMD",                       
  "nylon_wood",                                      
  "nylon_PS",                                        
  "nylon_PE",                                        
  "water",                                           
  "wood",                                            
  "additivated_PE_with_prooxidant_OXO",              
  "poly_3_hydroxybutyrate_co_3_hydroxyvalerate_PHBV",
  "sediment",                                        
  "LDPE",                                            
  "PVC",                                             
  "glass",                                           
  "rubber",                                          
  "metal",                                           
  "HDPE",                                            
  "Polylactic acid_PLA",                             
  "Polyurethane_PU",                                 
  "Polymethyl_methacrylate_PMMA",                    
  "PAHs",                                            
  "PCBs",                                            
  "stone"
),]
#First provide grouping column
meta_table$Groups<-as.character(meta_table$Material)
#The colours in the the next instruction match the factors for meta_table$Groups
meta_table$Groups<-factor(meta_table$Groups,c(
  "Microplastic_fibers",                           
  "PE",                                              
  "PP",                                              
  "PET",                                             
  "Pellets",                                         
  "PS",                                              
  "Plastic_marine_debris_PMD",                       
  "nylon_wood",                                      
  "nylon_PS",                                        
  "nylon_PE",                                        
  "water",                                           
  "wood",                                            
  "additivated_PE_with_prooxidant_OXO",              
  "poly_3_hydroxybutyrate_co_3_hydroxyvalerate_PHBV",
  "sediment",                                        
  "LDPE",                                            
  "PVC",                                             
  "glass",                                           
  "rubber",                                          
  "metal",                                           
  "HDPE",                                            
  "Polylactic acid_PLA",                             
  "Polyurethane_PU",                                 
  "Polymethyl_methacrylate_PMMA",                    
  "PAHs",                                            
  "PCBs",                                            
  "stone"
))
colours <- c(
  "#f86c6c", 
  "#fccece", 
  "#994b00", 
  "#cc6500", 
  "#ff7f00", 
  "#ffb266", 
  "#ccaf00", 
  "#ffdb00", 
  "#ffe966", 
  "#fef7cc", 
  "#00660b", 
  "#009a17", 
  "#00cc17", 
  "#ccfed1", 
  "#002766", 
  "#004ecc", 
  "#3380ff", 
  "#ccdffe", 
  "#570099", 
  "#9200ff", 
  "#bd66ff", 
  "#e9ccfe", 
  "#9BFF00",
  "#00A0FF",
  "#0024FF",
  "#FF00BB",
  "#E5FF00",
  
  #Next colors are for lines mainly used in the PCoA script
  "#000080","#4876FF","#CAE1FF","#9FB6CD","#1E90FF","#00F5FF","#00C957",grey.colors(1000));
#meta_table$Type is for shapes
meta_table$Type<-"Salinity_range_2"
meta_table$Connections<-NULL
meta_table$Subconnections<-NULL
PERMANOVA_variables<-c(
  "Study_number", "Material", "Temperature_range_2",	
  "Temperature_range",	"Temperature",	"Salinity_range_2",	"Salinity_range", "Salinity", 
  "pH", "geo_loc_name_country"
)

#PARAMETERS CHANGE THE GROUPING COLUMN AS YOU DESIRE############################

#Adjust abund_table to contain only those rows that got selected in the Hypothesis space
abund_table<-abund_table[rownames(meta_table),]

#After adjustment, get rid of OTUs that are all empty
abund_table<-abund_table[,colSums(abund_table)>0]
#Adjust OTU taxonomy
OTU_taxonomy<-OTU_taxonomy[colnames(abund_table),]


#COLLATE OTUS AT A PARTICULAR LEVEL#######################################
new_abund_table<-NULL
if(which_level=="Otus"){
  new_abund_table<-abund_table
} else {
  list<-unique(OTU_taxonomy[,which_level])
  new_abund_table<-NULL
  for(i in list){
    tmp<-data.frame(rowSums(abund_table[,rownames(OTU_taxonomy)[OTU_taxonomy[,which_level]==i],drop=FALSE]))
    if(i==""){colnames(tmp)<-c("__Unknowns__")} else {
      colnames(tmp)<-paste("",i,sep="")
      colnames(tmp)<-gsub(";+$","",paste(sapply(OTU_taxonomy[OTU_taxonomy[,which_level]==i,][1,1:which(colnames(OTU_taxonomy)==which_level)],as.character),collapse=";"))
    }
    if(is.null(new_abund_table)){new_abund_table<-tmp} else {new_abund_table<-cbind(tmp,new_abund_table)}
  }
}

new_abund_table<-as.data.frame(as(new_abund_table,"matrix"))
abund_table<-new_abund_table
#/COLLATE OTUS AT A PARTICULAR LEVEL#######################################


#Convert the data to phyloseq format
OTU = otu_table(as.matrix(abund_table), taxa_are_rows = FALSE)
TAX = tax_table(as.matrix(OTU_taxonomy))
SAM = sample_data(meta_table)

physeq<-NULL
if(which_level=="Otus"){
  physeq<-merge_phyloseq(phyloseq(OTU, TAX),SAM,OTU_tree)
} else {
  physeq<-merge_phyloseq(phyloseq(OTU),SAM)
}



#Reference: http://stackoverflow.com/questions/13794419/plotting-ordiellipse-function-from-vegan-package-onto-nmds-plot-created-in-ggplo
#Data frame df_ell contains values to show ellipses. It is calculated with function veganCovEllipse which is hidden in vegan package. This function is applied to each level of NMDS (group) and it uses also function cov.wt to calculate covariance matrix.
veganCovEllipse<-function (cov, center = c(0, 0), scale = 1, npoints = 100) 
{
  theta <- (0:npoints) * 2 * pi/npoints
  Circle <- cbind(cos(theta), sin(theta))
  t(center + scale * t(Circle %*% chol(cov)))
}

#coloring function  
gg_color_hue<-function(n){
  hues=seq(15,375,length=n+1)
  hcl(h=hues,l=65,c=100)[1:n]
}

sol<-NULL
#BUG FIX when you get NA in distance matrices
if(which_distance=="bray"){
  a<-phyloseq::distance(physeq,"bray")
  if(sum(is.na(a))>0){
    x = as.matrix(a)
    x = x[rowSums(is.na(x)) == 0, colSums(is.na(x)) == 0, drop = FALSE]
    meta_table<-meta_table[rownames(x),]
    a<-as.dist(x) 
  }
  if(which_method=="NMDS"){
    sol<-metaMDS(a, k = 2, trymax = 50)    
  } else if (which_method=="PCOA"){
    sol<-cmdscale(a,eig=T)
  }
  
} else if(which_distance=="wunifrac" & which_level=="Otus"){
  a<-phyloseq::distance(physeq,"wunifrac")
  if(sum(is.na(a))>0){
    x = as.matrix(a)
    x = x[rowSums(is.na(x)) == 0, colSums(is.na(x)) == 0, drop = FALSE]
    meta_table<-meta_table[rownames(x),]
    a<-as.dist(x) 
  }
  if(which_method=="NMDS"){
    sol<-metaMDS(a, k = 2, trymax = 50)    
  } else if (which_method=="PCOA"){
    sol<-cmdscale(a,eig=T)
  }
  
} else if(which_distance=="unifrac" & which_level=="Otus"){
  a<-phyloseq::distance(physeq,"unifrac")
  if(sum(is.na(a))>0){
    x = as.matrix(a)
    x = x[rowSums(is.na(x)) == 0, colSums(is.na(x)) == 0, drop = FALSE]
    meta_table<-meta_table[rownames(x),]
    a<-as.dist(x) 
  }
  if(which_method=="NMDS"){
    sol<-metaMDS(a, k = 2, trymax = 50)
  } else if (which_method=="PCOA"){
    sol<-cmdscale(a,eig=T)
  }
} else if(which_distance=="userprovided"){
  a<-as.dist(user_provided_distance[rownames(meta_table),rownames(meta_table)])
  if(sum(is.na(a))>0){
    x = as.matrix(a)
    x = x[rowSums(is.na(x)) == 0, colSums(is.na(x)) == 0, drop = FALSE]
    meta_table<-meta_table[rownames(x),]
    a<-as.dist(x) 
  }
  if(which_method=="NMDS"){
    sol<-metaMDS(a, k = 2, trymax = 50)
  } else if (which_method=="PCOA"){
    sol<-cmdscale(a,eig=T)
  }
}


if(!is.null(sol)){
  ORDINATION=data.frame(x=sol$points[,1],y=sol$points[,2],meta_table)
  
  plot.new()
  ord<-ordiellipse(sol, meta_table$Groups,display = "sites", kind = kind, conf = 0.95, label = T)
  dev.off()
  
  
  #Generate ellipse points
  df_ell <- data.frame()
  for(g in levels(ORDINATION$Groups)){
    if(g!="" && (g %in% names(ord))){
      
      tryCatch(df_ell <- rbind(df_ell, cbind(as.data.frame(with(ORDINATION[ORDINATION$Groups==g,],
                                                                veganCovEllipse(ord[[g]]$cov,ord[[g]]$center,ord[[g]]$scale)))
                                             ,Groups=g)),error=function(e) NULL)
    }
  }
  
  if (sum(dim(df_ell))>0){
    colnames(df_ell)<-c("x","y","Groups")
  }
  
  #Enforce df_ell to follow the same grouping even if some of the levels are dropped
  df_ell$Groups<-factor(as.character(df_ell$Groups),levels=levels(meta_table$Groups))
  
  
  #Generate mean values from ORDINATION plot grouped on
  ORDINATION.mean=aggregate(ORDINATION[,1:2],list(group=ORDINATION$Groups),mean)
  
  #Connecting samples based on lines with meta_table$Connections and meta_table$Subconnections ###########
  ORDINATION_lines<-NULL
  
  #Check if meta_table$Subconnections exists
  if(!is.null(meta_table$Subconnections)){
    #Step 1, populate Connections_IDs and get rid of singletons using ORDINATION$Connections
    Connections_IDs<-as.character(ORDINATION$Connections)
    Connections_IDs<-names(table(Connections_IDs)[table(Connections_IDs)>1])
    Subconnections_IDs_mask<-as.character(sapply(Connections_IDs,function(x){if(length(unique(ORDINATION[ORDINATION$Connections==x,"Subconnections"]))<2) x else "___APPROVED__"}))
    #Step 2, loop through each Connection_IDs and see we can find multiple ORDINATION$Subconnections and then filter Connections_IDs further
    Connections_IDs<-Connections_IDs[!Connections_IDs %in% Subconnections_IDs_mask]
    Connections_IDs<-sort(Connections_IDs)
    if(length(Connections_IDs)>0){
      for(p in Connections_IDs){
        Subconnections_IDs<-unique(ORDINATION[ORDINATION$Connections==p,"Subconnections"])
        Subconnections_IDs<-sort(Subconnections_IDs)
        S<-NULL
        if(pairwise_connections_and_not_longitudinal_connections){
          #Get pair-wise combinations from Subconnections_IDs
          S<-combn(Subconnections_IDs,2)
        } else{
          #Get longitudinal connections from Subconnections_IDs
          S<-t(cbind(Subconnections_IDs[-length(Subconnections_IDs)],Subconnections_IDs[-1]))
        }
        for(ii in 1:ncol(S)){
          tmp<-data.frame(t(colMeans(ORDINATION[ORDINATION$Connections==p & ORDINATION$Subconnections==S[1,ii],c("x","y"),drop=F])),t(colMeans(ORDINATION[ORDINATION$Connections==p & ORDINATION$Subconnections==S[2,ii],c("x","y"),drop=F])),p)
          colnames(tmp)<-c("xfrom","yfrom","xto","yto","ID")
          if(is.null(ORDINATION_lines)){ORDINATION_lines<-tmp} else {ORDINATION_lines<-rbind(ORDINATION_lines,tmp)}
        }
      }
    }
  } else {
    #To connect lines between samples we need to extract connections
    Connections_IDs<-as.character(ORDINATION$Connections)
    #Next we filter out Connections_IDs that are singletons and also uniquify them
    Connections_IDs<-names(table(Connections_IDs)[table(Connections_IDs)>1])
    Connections_IDs<-sort(Connections_IDs)
    if(length(Connections_IDs)>0){
      #We iterate through the IDs one at a time
      for(p in Connections_IDs){
        rownames_list<-rownames(ORDINATION[ORDINATION$Connections %in% p,,drop=F])
        S<-NULL
        if(pairwise_connections_and_not_longitudinal_connections){
          #Get pair-wise combinations from rownames_list
          S<-combn(rownames_list,2)
        } else{
          #Get longitudinal connections from rownames_list
          S<-t(cbind(rownames_list[-length(rownames_list)],rownames_list[-1]))
        }
        for(ii in 1:ncol(S)){
          tmp<-cbind(ORDINATION[S[1,ii],c("x","y")],ORDINATION[S[2,ii],c("x","y")],p)
          colnames(tmp)<-c("xfrom","yfrom","xto","yto","ID")
          if(is.null(ORDINATION_lines)){ORDINATION_lines<-tmp} else {ORDINATION_lines<-rbind(ORDINATION_lines,tmp)}
        }
      }
    }
  }
  #/#Connecting samples based on lines with meta_table$Connections and meta_table$Subconnections ###########
  
  
  cols=gg_color_hue(length(unique(ORDINATION$Groups))) #ORDINATION$Groups
  
  p<-ggplot(data=ORDINATION,aes(x,y,colour=Groups)) #colour=Groups
  
  p<-p + geom_point(aes(ORDINATION$x,ORDINATION$y,colour=ORDINATION$Groups),inherit.aes=F,alpha=point_opacity,size=point_size)
  if(!is.null(meta_table$Type)){
    p<-p + geom_point(aes(ORDINATION$x,ORDINATION$y,shape=Type),inherit.aes=F,size=point_size-6,colour="black") #colour=ORDINATION$Groups
  }
  if(draw_glow){
    p<-p + geom_point(alpha=point_glow_opacity,size = point_size+point_glow_differential,show.legend=FALSE)
  }
  
  p<-p+theme_bw()
  if(draw_mean_values_text){
    p<-p+geom_label_repel(size=12,inherit.aes = FALSE,aes(x=x,y= y,label=group), data=ORDINATION.mean,
                          box.padding   = 0.35, 
                          point.padding = 0.5,
                          segment.color = 'white',
                          max.overlaps =50,
                          #colour=colours[as.numeric(ORDINATION.mean$group)],
                          colour="white",
                          fill="black",
                          vjust=0.3
                          )
    }
  if (sum(dim(df_ell))>0){
    if(draw_confidence_intervals){
      if(draw_ellipses_and_not_polygons){
        p<-p+ geom_path(data=df_ell, aes(x=x, y=y), size=linesize_ellipses_polygons, linetype=linetype_ellipses_polygons,alpha=opacity_ellipses_polygons, show.legend = FALSE)
      } else {
        p<-p+ geom_polygon(data=df_ell, aes(x=x, y=y,fill=Groups), size=linesize_ellipses_polygons, linetype=linetype_ellipses_polygons,alpha=opacity_ellipses_polygons, show.legend = FALSE) #fill=Groups
      }
    }
  }
  if(exclude_legends){
    p<-p+guides(colour=FALSE)
  }
  
  
  
  if(which_method=="NMDS"){
    p<-p+xlab("NMDS1")+ylab("NMDS2")+ggtitle(paste("Stress=",sprintf("%.4g",sol$stress),sep=""))  
  } else if (which_method=="PCOA"){
    p<-p+xlab(paste("Dim1 (",sprintf("%.4g",(sol$eig[1]/sum(sol$eig))*100),"%)",sep=""))+ylab(paste("Dim2 (",sprintf("%.4g",(sol$eig[2]/sum(sol$eig))*100),"%)",sep=""))
  }
  p<-p+theme(legend.title=element_text(size=legend_title_size),
             legend.text=element_text(size=legend_text_size),
             text = element_text(size=text_size),
             axis.text=element_text(size=axis_text_size),
             axis.title=element_text(size=axis_title_size))
  
  
  if(use_provided_colors){
    p<-p+scale_color_manual("Groups",values=colours) #"Groups"
    p<-p+scale_fill_manual("Groups",values=colours[unique(as.numeric(df_ell$Groups))]) #"Groups"
  }
  if(!is.null(meta_table$Type)){
    # This sets shapes as symbols the , c is when you want to add gap/skip certain symbols
    #p<-p+scale_shape_manual("Type",values=c(c(17:25),c(33:127)))
    
    #p<-p+scale_shape_manual("Type",values=c(c(6),c(8)
    #))
    p<-p+scale_shape_manual("Type",values=c(c(17:20),
                                            c(22:25),
                                            c(0:14),
                                           c(35),
                                           c(42:43)
   ))
    # This changes shapes to alphabet letters
    #p<-p+scale_shape_manual("Type",values=c(65:122))
  }
  
  #only draw lines connecting dots if the lines are available
  if(!is.null(ORDINATION_lines)){
    arrow<-NULL
    if(should_connections_end_in_arrows){
      arrow=arrow(length=unit(0.2,"inches"))
    }
    cl=nlevels(meta_table$Groups)
    for(i in unique(meta_table$Connections)){
      cl<-cl+1
      if(should_we_reverse_direction_of_arrows){
        p<-p+geom_segment(data=ORDINATION_lines[as.character(ORDINATION_lines$ID)==i,],inherit.aes=FALSE,aes(x=xto,y=yto,xend=xfrom,yend=yfrom),colour=colours[cl],size=linking_samples_line_size,alpha=linking_samples_line_opacity,linetype=linking_samples_linetype, show.legend = NA,arrow=arrow)
      } else {
        p<-p+geom_segment(data=ORDINATION_lines[as.character(ORDINATION_lines$ID)==i,],inherit.aes=FALSE,aes(x=xfrom,y=yfrom,xend=xto,yend=yto),colour=colours[cl],size=linking_samples_line_size,alpha=linking_samples_line_opacity,linetype=linking_samples_linetype, show.legend = NA,arrow=arrow)
      }
      p<-p + guides(colour = guide_legend(override.aes = list(shape = 15)))
    }
  }
  
  p<-p+theme(legend.key = element_rect(fill = "grey"))
  pdf(paste(which_method,"_",which_distance,"_",which_level,"_",label,".pdf",sep=""),width=width_image,height=height_image)
  print(p)
  dev.off()
  
  dist<-NULL
  if(which_distance=="userprovided"){
    dist<-as.dist(user_provided_distance[rownames(meta_table),rownames(meta_table)])
  } else {
    dist<-phyloseq::distance(physeq,which_distance)
  }
  #BUGFIX for NA
  if(sum(is.na(dist))>0){
    x = as.matrix(dist)
    x = x[rowSums(is.na(x)) == 0, colSums(is.na(x)) == 0, drop = FALSE]
    dist<-as.dist(x) 
  }
  
  capture.output(adonis2(as.formula(paste("dist ~",paste(PERMANOVA_variables,collapse="+"))), data=meta_table[rownames(as.matrix(dist)),]),file=paste("ADONIS_",which_distance,"_",which_level,"_",label,".txt",sep=""))
  
}
save.image(paste("ORDINATE_",which_distance,"_",which_level,"_",label,".RData",sep=""))