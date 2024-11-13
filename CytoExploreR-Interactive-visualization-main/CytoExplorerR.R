#############################################
#Title: "Interactive Analysis of Cytometry Data in R"
#Author: "Pritam Kumar Panda"
#Institution: "German Cancer Research Center DKFZ"
#Country: "Germany"
#############################################


#Set your working directory
setwd("/Users/pritam/facs_data_analysis_R/FACS/CytoExploreR")

#initiate renv if you have migrated to another computer
#renv::init()

#increase memory and use cache of larger files
knitr::opts_chunk$set(cache = TRUE, warning = FALSE, 
                      message = FALSE, cache.lazy = FALSE)

#Install openblas and lapack for graphics rendering
#brew install openblas
#brew install lapack
#export LDFLAGS="-L/usr/local/opt/openblas/lib"                            
#export CPPFLAGS="-I/usr/local/opt/openblas/include"

#Install all the packages
install.packages("tidyverse")
install.packages("knitr") 
install.packages("ggplot2")
install.packages("remotes") 
install.packages("BiocManager")
install.packages("devtools") 
install.packages("kableExtra")
install.packages("devtools")

#Bioconductor packages if (!require("BiocManager", quietly = TRUE))
install.packages("BiocManager") 
BiocManager::install("flowCore")
BiocManager::install("flowWorkspace") 
BiocManager::install("openCyto")
BiocManager::install("flowAI") 
BiocManager::install("ggcyto")
BiocManager::install("CytoML")

# CytoExploreR

devtools::install_github("DillonHammill/CytoExploreR")
devtools::install_github("DillonHammill/CytoExploreRData")

#Load Required Packages
#CytoExploreR is closely integrated with the RGLab's suite of cytometry packages, including flowCore, flowWorkspace
#and openCyto. This means that these packages will be automatically
#loaded when CytoExploreR is loaded and any functions from these packages
# are completely compatible with CytoExploreR. \# To demonstrate the
#features of CytoExploreR we will also need to load CytoExploreRData
#which contains the Compensation and Activation datasets.


# Load required packages
library(flowCore)
library(CytoML)
library(flowAI)
library(flowWorkspace)
library(ggcyto)
library(tidyverse)
library(openCyto)
library(knitr)
library(kableExtra)
library(dplyr)
library(CytoExploreR)
library(CytoExploreRData)

# #Save FCS Files
# Users should create a new folder in the current working directory and add
#their FCS files into this folder. It is recommended that you store your
#compensation controls and samples in separate folders as these files will be analysed separately. 
# The following snippet outlines the steps required to download FCS files for the Compensation and
#Activation datasets and save them to appropriately named folders in your current working directory.

# Compensation FCS Files
cyto_save(Compensation, 
          save_as = "Compensation-Samples")

# Activation FCS Files
cyto_save(Activation,
          save_as = "Activation-Samples")

#Compensation of Fluorescent Spillover
##Prepare compensation controls
# Setup compensation controls
gs <- cyto_setup("Compensation-Samples",
                 gatingTemplate = "Compensation-gatingTemplate.csv")

#Abbreviations used
# FSC: Forward Scatter
# FSC.A: Forward Scatter Area
# FSC.H: Forward Scatter Height
# SSC: Side Scatter
# SSC.A: Side Scatter Area
# SSC.H: Side Scatter Height
# FITC: Fluorescein Isothiocyanate
# APC: Allophycocyanin
# APC.A: Allophycocyanin Area
# APC.R700.A: APC R700 Area
# APC.H7.A: APC H7 Area
# APC is another fluorescent dye used in flow cytometry
# PE: Phycoerythrin
# PE.CF594.A: Phycoerythrin CF594 Area
# PE.Cy5.A: Phycoerythrin Cy5 Area
# PE.Cy7.A: Phycoerythrin Cy7 Area
# PE.A: Phycoerythrin Area
# BV: Brilliant Violet


# Transform fluorescent channels - default logicle transformations
gs <- cyto_transform(gs, type="logicle")

# Gate Cells
cyto_gate_draw(gs,
               parent = "root",
               alias = "Cells",
               channels = c("FSC-A", "SSC-A"))

#Many type of gates available. You can specify by using type=boundary"
#type ="polygon"," ellipsoid", interval, threshold, quadrant, 
#also specify multiple gates : type =c("rectangle, "ellipse")

# Gate Single Cells
cyto_gate_draw(gs,
               parent = "Cells",
               alias = "Single Cells",
               channels = c("FSC-A", "FSC-H"))

# Points to consider
# Alexa Fluor 488-A: This is a green-fluorescent dye often used to label antibodies, peptides, or proteins. The "A" at the end usually denotes the form associated with an antibody conjugate.
# PE-A: Stands for Phycoerythrin-A, a bright orange-fluorescent dye conjugated to an antibody. It's used as a fluorescent marker in cell biology.
# PE-Texas Red-A: A tandem dye that combines Phycoerythrin (PE) with Texas Red, enhancing the emission spectrum into the red wavelengths when conjugated with an antibody.
# 7-AAD-A: 7-Aminoactinomycin D is a fluorescent dye that binds selectively to GC-rich regions of DNA, often used in viability assays and cell cycle studies.
# PE-Cy7-A: Another tandem dye where Phycoerythrin is coupled with a cyanine dye (Cy7), shifting its fluorescence emission to the far-red region, conjugated with antibodies for use in applications like flow cytometry.
# Alexa Fluor 405-A: A violet-fluorescent dye that's conjugated with antibodies and other proteins for use in a variety of assays.
# Alexa Fluor 430-A: A dye with fluorescence emission in the green spectrum when conjugated to a target molecule such as an antibody.
# Qdot 605-A: Quantum dot 605 is a nanoparticle with distinct optical and electronic properties, emitting at a wavelength of 605 nm, and is often used in biological labeling and imaging.
# Alexa Fluor 647-A: A red-fluorescent dye that's used widely for labeling and is useful in multi-color applications due to its far-red emission.
# Alexa Fluor 700-A: A near-infrared fluorescent dye, which allows for detection of targets in the near-infrared range.
# APC-Cy7-A: Allophycocyanin (APC) conjugated with cyanine dye 7 (Cy7), emitting in the far-red to near-infrared range, often used for flow cytometry.

#Compute Spillover Matrix
spill <- cyto_spillover_compute(gs,
                                parent = "Single Cells",
                                spillover = "Spillover-Matrix.csv")

# channel description
# Compensation-7AAD.fcs: 7-AAD is typically excited by a red laser (633 or 640 nm) and its emission is detected in the far-red area of the spectrum, so you should select the channel that detects far-red.
# Compensation-AF700.fcs: Alexa Fluor 700 is excited by a red laser and emits in the near-infrared region, so choose the channel that detects near-infrared fluorescence.
# Compensation-APC-Cy7.fcs: APC-Cy7 is a tandem dye that is excited by a red laser and emits in the far-red to near-infrared range. Select the channel that is appropriate for far-red to near-infrared emission.
# Compensation-APC.fcs: Allophycocyanin (APC) is also excited by a red laser and emits in the far-red, so choose the channel that detects far-red fluorescence.
# Compensation-FITC.fcs: FITC is typically excited by a blue laser (488 nm) and emits in the green region of the spectrum. Select the channel that corresponds to green fluorescence.
# Compensation-PE.fcs: Phycoerythrin (PE) is best excited by a blue or green laser and emits in the yellow to orange region. Choose the channel detecting the orange part of the spectrum.
# Compensation-Unstained.fcs: For the unstained control, you would typically look at all the channels to ensure that there is no fluorescence being detected, as this sample should not have any fluorescent labels.

#selecting channels
# Compensation-7AAD.fcs: Select the channel 7-AAD-A because 7-AAD is detected in this channel.
# Compensation-AF700.fcs: Choose Alexa Fluor 700-A as this channel is designed for the detection of Alexa Fluor 700 dye.
# Compensation-APC-Cy7.fcs: You should use the APC-Cy7-A channel, which corresponds to the APC-Cy7 tandem dye.
# Compensation-APC.fcs: There isn't a channel explicitly listed for APC in the second image. APC is typically detected in a channel similar to Alexa Fluor 647 due to their close fluorescence emission spectra, so you would normally choose Alexa Fluor 647-A. However, if your flow cytometer setup has a dedicated APC channel not listed, you should use that one.
# Compensation-FITC.fcs: Select Alexa Fluor 488-A because FITC and Alexa Fluor 488 have similar excitation and emission properties and are often detected in the same channel.
# Compensation-PE.fcs: Choose PE-A as this channel is intended for the detection of the PE dye.
# Compensation-Unstained.fcs: For the unstained control, you would check all channels to ensure that there is no fluorescence being detected.

spill

#Interactively Edit Spillover Matrices
# Open CytoExploreR spillover matrix editor
spill <- cyto_spillover_edit(gs,
                             parent = "Single Cells",
                             spillover = "Spillover-Matrix.csv")

#Visualise Compensation
cyto_plot_compensation(gs,
                       parent = "Single Cells")


# Visualise compensated data
cyto_plot_compensation(gs,
                       parent = "Single Cells",
                       spillover = "Spillover-Matrix.csv",
                       compensate = TRUE)

#Samples
###########################################################################
#Now analyze the samples (treated samples)
# Load and annotate samples
gs <- cyto_setup("Activation-Samples",
                 gatingTemplate = "Activation-gatingTemplate.csv")

#If you have messed up , then you can clean the matrix
#gs <- cyto_setup("Activation-Samples",
#gatingTemplate = "Activation-gatingTemplate.csv", clean = TRUE)

# Apply compensation
gs <- cyto_compensate(gs,
                      spillover = "Spillover-Matrix.csv")

# Transform fluorescent channels - default logicle transformations
gs <- cyto_transform(gs)

# Gate Cells in the samples
cyto_gate_draw(gs,
               parent = "root",
               alias = "Cells",
               channels = c("FSC-A","SSC-A"))

# Gate Single Cells
cyto_gate_draw(gs,
               parent = "Cells",
               alias = "Single Cells",
               channels = c("FSC-A","FSC-H"))

# Extract unstained control 
NIL <- cyto_extract(gs, "Single Cells")[[33]]

# Gate Live Cells
cyto_gate_draw(gs,
               parent = "Single Cells",
               alias = c("Dead Cells", "Live Cells"),
               channels = c("Hoechst-405", "Hoechst-430"),
               type = "rectangle",
               negate = TRUE,
               overlay = NIL)

#If you have done something wrong, then you can edit the template.
cyto_gatingTemplate_edit(gs,
                         gatingTemplate ="Activation-gatingTemplate.csv")

# Gate T Cells and Dedritic Cells
cyto_gate_draw(gs,
               parent = "Live Cells",
               alias = c("T Cells", "Dendritic Cells"),
               channels = c("CD11c", "Va2"),
               type = c("ellipse", "rectangle"),
               overlay = NIL)

# Gate CD4 & CD8 T Cells
cyto_gate_draw(gs,
               parent = "T Cells",
               alias = c("CD4 T Cells", "CD8 T Cells"),
               channels = c("CD4", "CD8"),
               type = "r")

# Extract CD4 T Cells
CD4 <- cyto_extract(gs, "CD4 T Cells")

# Extract naive CD4 T Cells
CD4_naive <- cyto_select(CD4, OVAConc = 0)

# Gate CD69+ CD4 T Cells
cyto_gate_draw(gs,
               parent = "CD4 T Cells",
               alias = "CD69+ CD4 T Cells",
               channels = c("CD44", "CD69"),
               type = "rectangle",
               overlay = CD4_naive)


# Extract CD4 T Cells
CD4 <- cyto_extract(gs, "CD4 T Cells")

# Extract naive CD4 T Cells
CD4_naive <- cyto_select(CD4, OVAConc = 0)

# Gate CD69+ CD4 T Cells
cyto_gate_draw(gs,
               parent = "CD4 T Cells",
               alias = "CD69+ CD4 T Cells",
               channels = c("CD44", "CD69"),
               type = "rectangle",
               overlay = CD4_naive)
# Gate CD69+ CD8 T Cells
cyto_gate_draw(gs,
               parent = "CD8 T Cells",
               alias = "CD69+ CD8 T Cells",
               channels = c("CD44", "CD69"),
               type = "rectangle",
               contour_lines = 15)

# Gating Tree
cyto_plot_gating_tree(gs[[32]],
                      stat = "freq")

# Gating scheme
cyto_plot_gating_scheme(gs[32],
                        back_gate = TRUE,
                        gate_track = TRUE)


#Export Statistics
# Compute medFI - exclude unstained control
cyto_stats_compute(gs[1:32],
                   alias = c("CD69+ CD4 T Cells",
                             "CD69+ CD8 T Cells"),
                   stat = "median",
                   channels = c("CD44", "CD69"))





