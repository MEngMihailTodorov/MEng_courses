#ACHTUNG: From the RStudio menu, click on "Session" and then "Set Working Directory" to "To Source File Location"
#Script for calculation of diversity indices, and also supports paired data
#v1.2 (All metrics are now being saved)


#Note: For meta-analysis work, please be mindful of how you load the biom file
# In the scripts where you use phyloseq, you have two options to 
# load the biom file (both will generate the physeq object). The
# reason why we have these two options is mainly because write_biom()
# from biomformat package saves biom file in a format that is not readable
# using import_biom() function. So we are introducing an extra step:
#
# physeq<-import_biom("collated_feature_w_tax.biom")
# library(biomformat);b_<-read_biom("collated_feature_w_tax.biom");physeq<-merge_phyloseq(otu_table(as(biom_data(b_),"matrix"),taxa_are_rows=TRUE),tax_table(as(observation_metadata(b_),"matrix")))

library(phyloseq) #Bioconductor 
library(stringr)
library(data.table)
library(vegan)
library(ggplot2)
library(grid) #We need grid to draw the arrows
library(ggbeeswarm)

#PARAMETERS ###########################
which_level="Otus" #Otus Genus Family Order Class Phylum
library(biomformat);b_<-read_biom("../../Data/collated_feature_w_tax.biom");physeq<-merge_phyloseq(otu_table(as(biom_data(b_),"matrix"),taxa_are_rows=TRUE),tax_table(as(observation_metadata(b_),"matrix")))
meta_table<-read.csv("../../Data/meta_data.csv",header=T,row.names=1)
text_size=16
axis_text_size=14
strip_text_size=18
increment_divider=3
exclude_pvalues_text_from_drawing=FALSE
legends_position_bottom=FALSE
exclude_legends=TRUE #FALSE
pairwise_text_size=10
number_of_rows=1
legend_text_size=20
legend_title_size=22
axis_title_size=20
height_image=36
width_image=60
smoothing_alpha=0.3
use_provided_colors=TRUE
turn_smoothing_on=FALSE
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
OTU_taxonomy$Kingdom<-gsub("D_0__|d__","",OTU_taxonomy$Kingdom)

#Remove singletons and adjust OTU_taxonomy
abund_table<-abund_table[,colSums(abund_table)>1]
OTU_taxonomy<-OTU_taxonomy[colnames(abund_table),]

#get rid of contaminants with "Unassigned", "Chloroplast" and "Mitochondria" assignment", and "non classified" at Phylum level
abund_table<-abund_table[,!(OTU_taxonomy$Kingdom %in% c("Unassigned") | OTU_taxonomy$Phylum=="" | OTU_taxonomy$Order %in% c("Chloroplast") | OTU_taxonomy$Family %in% c("Mitochondria"))]


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
#and then allocate a column to meta_table$Groups that you want to use.

label="Hypothesis1"
meta_table<-meta_table[meta_table$Material_updated %in% c(
  "Microplastic_fibers",                             
  "Marine_plastic_litter_MPL_PE",                    
  "Polypropylene_PP",
  "Single_use_polyethylene_terephthalate_PET",
  "Pellets",
  "Polystyrene_PS",
  "Polyethylene_PE",
  "Plastic_marine_debris_PMD",
  "nylon_wood",
  "nylon_PS",
  "nylon_PE",
  "additivated_PE_with_prooxidant_OXO",
  "poly_3_hydroxybutyrate_co_3_hydroxyvalerate_PHBV",
  "Polyethylene_terephthalate_PET",
  "LDPE",
  "Polyvinyl_chloride_PVC",
  "Polyvinyl_chloride_PVC_DEHP",
  "HDPE",
  "Polyvinyl_chloride_PVC_DINP",
  "Macroplastics",
  "Polylactic_acid_PLA",
  "Polyurethane_foam_PUF",
  "Polyurethane_PU",
  "Polymethyl_methacrylate_PMMA",
  "PAHs",
  "PCBs",
  "PHBV_Polyhydroxyalkanoate_PHA"
  ),]
#First provide grouping column
meta_table$Groups<-as.character(meta_table$Material_updated)
#The colours in the the next instruction match the factors for meta_table$Groups
meta_table$Groups<-factor(meta_table$Groups,c(
  "Microplastic_fibers",                             
  "Marine_plastic_litter_MPL_PE",                    
  "Polypropylene_PP",
  "Single_use_polyethylene_terephthalate_PET",
  "Pellets",
  "Polystyrene_PS",
  "Polyethylene_PE",
  "Plastic_marine_debris_PMD",
  "nylon_wood",
  "nylon_PS",
  "nylon_PE",
  "additivated_PE_with_prooxidant_OXO",
  "poly_3_hydroxybutyrate_co_3_hydroxyvalerate_PHBV",
  "Polyethylene_terephthalate_PET",
  "LDPE",
  "HDPE",
  "Polyvinyl_chloride_PVC",
  "Polyvinyl_chloride_PVC_DEHP",
  "Polyvinyl_chloride_PVC_DINP",
  "Macroplastics",
  "Polylactic_acid_PLA",
  "Polyurethane_foam_PUF",
  "Polyurethane_PU",
  "Polymethyl_methacrylate_PMMA",
  "PAHs",
  "PCBs",
  "PHBV_Polyhydroxyalkanoate_PHA"
))
colours <- c(
  "#920606", 
  "#f40a0a", 
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
  "#686666", 
  "#c0bfbf",
   #Next colors are for lines mainly used in the PCoA script
   "#000080","#4876FF","#CAE1FF","#9FB6CD","#1E90FF","#00F5FF","#00C957",grey.colors(1000));
#meta_table$Type is for shapes
meta_table$Type<-NULL
provide_your_own_pvalue_combinations<-TRUE
provided_combination<-cbind(
  #Cross-sectional comparisons
  combn(c(
    "Microplastic_fibers",                             
    "Marine_plastic_litter_MPL_PE",                    
    "Polypropylene_PP",
    "Single_use_polyethylene_terephthalate_PET",
    "Pellets",
    "Polystyrene_PS",
    "Polyethylene_PE",
    "Plastic_marine_debris_PMD",
    "nylon_wood",
    "nylon_PS",
    "nylon_PE",
    "additivated_PE_with_prooxidant_OXO",
    "poly_3_hydroxybutyrate_co_3_hydroxyvalerate_PHBV",
    "Polyethylene_terephthalate_PET",
    "LDPE",
    "HDPE",
    "Polyvinyl_chloride_PVC",
    "Polyvinyl_chloride_PVC_DEHP",
    "Polyvinyl_chloride_PVC_DINP",
    "Macroplastics",
    "Polylactic_acid_PLA",
    "Polyurethane_foam_PUF",
    "Polyurethane_PU",
    "Polymethyl_methacrylate_PMMA",
    "PAHs",
    "PCBs",
    "PHBV_Polyhydroxyalkanoate_PHA"
  ),2)
)
meta_table$Connections<-NULL


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
      #colnames(tmp)<-paste("",i,sep="")
      colnames(tmp)<-gsub(";+$","",paste(sapply(OTU_taxonomy[OTU_taxonomy[,which_level]==i,][1,1:which(colnames(OTU_taxonomy)==which_level)],as.character),collapse=";"))
    }
    if(is.null(new_abund_table)){new_abund_table<-tmp} else {new_abund_table<-cbind(tmp,new_abund_table)}
  }
}

new_abund_table<-as.data.frame(as(new_abund_table,"matrix"))
abund_table<-new_abund_table
#/COLLATE OTUS AT A PARTICULAR LEVEL#######################################

grouping_column<-"Groups"

#Calculate Richness
R<-vegan::rarefy(abund_table,min(rowSums(abund_table)))
df_R<-data.frame(sample=names(R),value=R,measure=rep("Richness",length(R)))

#Calculate Shannon entropy
H<-vegan::diversity(abund_table)
df_H<-data.frame(sample=names(H),value=H,measure=rep("Shannon",length(H)))

#Calculate Simpson diversity index
simp <- vegan::diversity(abund_table, "simpson")
df_simp<-data.frame(sample=names(simp),value=simp,measure=rep("Simpson",length(simp)))

#Calculate Fisher alpha
alpha <- vegan::fisher.alpha(abund_table)
df_alpha<-data.frame(sample=names(alpha),value=alpha,measure=rep("Fisher alpha",length(alpha)))

#Calculate Pielou's evenness
S <- vegan::specnumber(abund_table)
J <- H/log(S)
df_J<-data.frame(sample=names(J),value=J,measure=rep("Pielou's evenness",length(J)))

#Uncomment to retain everything
df<-rbind(df_R,df_H,df_simp,df_alpha,df_J)

rownames(df)<-NULL

#Incorporate categorical data in df
df<-data.frame(df,meta_table[as.character(df$sample),])

#To do anova, we will convert our data.frame to data.table

#Since we can't pass a formula to data.table, I am creating
#a dummy column .group. so that I don't change names in the formula
dt<-data.table(data.frame(df,.group.=df[,grouping_column]))

#I am also specifying a p-value cutoff for the ggplot2 strips
pValueCutoff<-0.05
pval<-dt[, list(pvalue = sprintf("%.2g", 
                                 tryCatch(summary(aov(value ~ .group.))[[1]][["Pr(>F)"]][1],error=function(e) NULL))), 
         by=list(measure)]

#Filter out pvals that we don't want
pval<-pval[!pval$pvalue=="",]
pval<-pval[as.numeric(pval$pvalue)<=pValueCutoff,]

#I am using sapply to generate significances for pval$pvalue using the cut function.
pval$pvalue<-sapply(as.numeric(pval$pvalue),function(x){as.character(cut(x,breaks=c(-Inf, 0.001, 0.01, 0.05, Inf),label=c("***", "**", "*", "")))})

#Update df$measure to change the measure names if the grouping_column has more than three classes
if(length(unique(as.character(meta_table[,grouping_column])))>2){
  df$measure<-as.character(df$measure)
  if(dim(pval)[1]>0){
    for(i in seq(1:dim(pval)[1])){
      df[df$measure==as.character(pval[i,measure]),"measure"]=paste(as.character(pval[i,measure]),as.character(pval[i,pvalue]))
    }
  }
  df$measure<-as.factor(df$measure)
}

#Get all possible combination of values in the grouping_column
s<-combn(unique(as.character(df[,grouping_column])),2)

#df_pw will store the pair-wise p-values
df_pw<-NULL
for(k in unique(as.character(df$measure))){
  #We need to calculate the coordinate to draw pair-wise significance lines
  #for this we calculate bas as the maximum value
  bas<-max(df[(df$measure==k),"value"])
  
  #Calculate increments as % of the maximum values
  inc<-0.05*(bas-min(df[(df$measure==k),"value"]))
  
  #Give an initial increment
  bas<-bas+inc
  for(l in 1:dim(s)[2]){
    
    tmp<-NULL
    #If it is paired-data we are interested in
    if(is.null(meta_table$Connections)){
      #Do a normal anova
      cat("Do a normal anova\n")
      tmp<-c(k,s[1,l],s[2,l],bas,paste(sprintf("%.2g",tryCatch(summary(aov(as.formula(paste("value ~",grouping_column)),data=df[(df$measure==k) & (df[,grouping_column]==s[1,l] | df[,grouping_column]==s[2,l]),] ))[[1]][["Pr(>F)"]][1],error=function(e) NULL)),"",sep=""))
    } else {
      #Do a paired anova with Error(Connections/Groups)
      cat("Do a paired anova with Error(Connections/Groups)\n")
      tmp<-c(k,s[1,l],s[2,l],bas,paste(sprintf("%.2g",tryCatch(summary(aov(as.formula(paste("value ~",grouping_column,"+",paste("Error(","Connections/",grouping_column,")",sep=""))),data=df[(df$measure==k) & (df[,grouping_column]==s[1,l] | df[,grouping_column]==s[2,l]),] ))[[1]][[1]][["Pr(>F)"]][1],error=function(e) NULL)),"",sep=""))
    }
    
    
    #Ignore if anova fails
    if(!is.na(as.numeric(tmp[length(tmp)]))){
      
      #Only retain those pairs where the p-values are significant
      if(as.numeric(tmp[length(tmp)])<0.05){
        if(is.null(df_pw)){df_pw<-tmp}else{df_pw<-rbind(df_pw,tmp)}
        
        #Generate the next position
        bas<-bas+inc
      }
    }
  }  
}

if(!is.null(df_pw)){
  if(sum(class(df_pw) %in% c("character"))>0){
    df_pw<-t(as.matrix(df_pw))
  }
  df_pw<-data.frame(row.names=NULL,df_pw)
  names(df_pw)<-c("measure","from","to","y","p")
}


#Draw the boxplots
p<-ggplot(aes_string(x=grouping_column,y="value",color=grouping_column,fill=grouping_column,color=grouping_column,group=grouping_column),data=df)
#p<-p+geom_boxplot(outlier.size=0)+geom_jitter(position = position_jitter(height = 0, width=0))
#p<-p+geom_violin(alpha=0.5)
p<-p+geom_boxplot(alpha=0.3)
#p<-p+geom_beeswarm(size=3)
p<-p+theme_bw()
p<-p+geom_point(size=3,alpha=0.2)
p<-p+facet_wrap(~measure,scales="free_y",nrow=number_of_rows)+ylab("Observed Values")+xlab("Samples")

if(!is.null(df_pw)){
  #This loop will generate the lines and signficances
  for(i in 1:dim(df_pw)[1]){
    p<-p+geom_path(inherit.aes=F,aes(x,y),data = data.frame(x = c(which(levels(df[,grouping_column])==as.character(df_pw[i,"from"])),which(levels(df[,grouping_column])==as.character(df_pw[i,"to"]))), y = c(as.numeric(as.character(df_pw[i,"y"])),as.numeric(as.character(df_pw[i,"y"]))), measure=c(as.character(df_pw[i,"measure"]),as.character(df_pw[i,"measure"]))), color="black",lineend = "butt",arrow = arrow(angle = 90, ends = "both", length = unit(0.1, "inches")))
    p<-p+geom_text(inherit.aes=F,aes(x=x,y=y,label=label),data=data.frame(x=(which(levels(df[,grouping_column])==as.character(df_pw[i,"from"]))+which(levels(df[,grouping_column])==as.character(df_pw[i,"to"])))/2,y=as.numeric(as.character(df_pw[i,"y"])),measure=as.character(df_pw[i,"measure"]),label=as.character(cut(as.numeric(as.character(df_pw[i,"p"])),breaks=c(-Inf, 0.001, 0.01, 0.05, Inf),label=c("***", "**", "*", "")))),size=pairwise_text_size)
    if(exclude_pvalues_text_from_drawing){
      p<-p+geom_text(inherit.aes=F,aes(x=x,y=y,label=label),data=data.frame(x=(which(levels(df[,grouping_column])==as.character(df_pw[i,"from"]))+which(levels(df[,grouping_column])==as.character(df_pw[i,"to"])))/2,y=as.numeric(as.character(df_pw[i,"y"])),measure=as.character(df_pw[i,"measure"]),label=paste("p=",as.character(as.numeric(as.character(df_pw[i,"p"])))),sep=""),size=pairwise_text_size,vjust=-1)
    }
  }
}



if(use_provided_colors){
  p<-p+scale_color_manual(grouping_column,values=colours)
  p<-p+scale_fill_manual(grouping_column,values=colours)
}

#If connections exist, then connect them with a line
if(!is.null(meta_table$Connections)){
  p<-p+geom_line(aes(x=Groups,y=value,group=Connections),colour="black",alpha=0.9,linetype="dotted",data=df)
}

#if crashes at panel.margin change it to panel.spacing, if crashes at panel.spacing, change it to panel.margin
p<-p+theme(strip.background = element_rect(fill="white"))+theme(panel.spacing = unit(2, "lines"),
                                                                strip.text = element_text(size=strip_text_size),
                                                                legend.text=element_text(size=legend_text_size),
                                                                text = element_text(size=text_size),
                                                                axis.text=element_text(size=axis_text_size),
                                                                axis.title=element_text(size=axis_title_size),
                                                                axis.text.x = element_text(angle = 90, hjust = 1))
if(legends_position_bottom){
  p<-p+theme(legend.key = element_blank(),  #removes the box around each legend item
             legend.position = "bottom", #legend at the bottom
             legend.direction = "horizontal",
             legend.box = "horizontal",
             legend.box.just = "centre")
}
if(exclude_legends){
  p<-p+guides(colour="none")
  p<-p+guides(fill="none")
}

pdf(paste("ANOVA_diversity_",which_level,"_",label,".pdf",sep=""),height=height_image,width=width_image)
print(p)
dev.off()

