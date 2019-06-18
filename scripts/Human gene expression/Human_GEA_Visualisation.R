
###############################################################################################################################
###############################################################################################################################

### Project: PhenoExpression ##################################################################################################
### Script: Human_GEA_Visualisation.R ###################################################################################################
### Purpose: To create venn diagrams to show statisticallt significant genes taht may overlap betwene thresholds ###############################################################################
### Author: Friha Zafar ####################################################################################################
### Date: 18/06/2019 ##########################################################################################################

###############################################################################################################################
###############################################################################################################################

## load packages #############################################################################################################

library( ggplot2 );library( scales );library(RColorBrewer);library(cowplot)

###############################################################################################################################################

# A plotting R script produced by the REVIGO server at http://revigo.irb.hr/
# If you found REVIGO useful in your work, please cite the following reference:
# Supek F et al. "REVIGO summarizes and visualizes long lists of Gene Ontology
# terms" PLoS ONE 2011. doi:10.1371/journal.pone.0021800

########################### TPM >1 MF : ######################################


# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","log10_p_value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0001540","beta-amyloid binding", 0.191, 3.062,-4.495, 1.531,-2.1423,0.897,0.000),
                     c("GO:0005215","transporter activity", 7.829,-1.845, 7.379, 3.132,-5.8409,0.905,0.000),
                     c("GO:0022835","transmitter-gated channel activity", 0.283,-6.640, 0.522, 1.699,-9.9192,0.203,0.000),
                     c("GO:0060089","molecular transducer activity", 9.545, 4.750, 0.153, 3.218,-1.9497,0.907,0.000),
                     c("GO:0016595","glutamate binding", 0.058, 2.949, 4.758, 1.041,-1.4861,0.888,0.004),
                     c("GO:0035240","dopamine binding", 0.075,-1.669,-7.232, 1.146,-1.4185,0.888,0.090),
                     c("GO:0015081","sodium ion transmembrane transporter activity", 0.884,-5.686,-0.428, 2.188,-1.4861,0.347,0.496),
                     c("GO:0015108","chloride transmembrane transporter activity", 0.532,-6.490,-0.805, 1.968,-1.4185,0.369,0.564),
                     c("GO:0022839","ion gated channel activity", 0.260,-5.964, 0.792, 1.663,-4.5321,0.318,0.658));

one.data.MF.1<- data.frame(revigo.data);
names(one.data.MF.1) <- revigo.names;
one.data.MF.1 <- one.data.MF.1 [(one.data.MF.1$plot_X != "null" & one.data.MF.1$plot_Y != "null"), ];
one.data.MF.1$plot_X <- as.numeric( as.character(one.data.MF.1$plot_X) );
one.data.MF.1$plot_Y <- as.numeric( as.character(one.data.MF.1$plot_Y) );
one.data.MF.1$plot_size <- as.numeric( as.character(one.data.MF.1$plot_size) );
one.data.MF.1$log10_p_value <- as.numeric( as.character(one.data.MF.1$log10_p_value) );
one.data.MF.1$frequency <- as.numeric( as.character(one.data.MF.1$frequency) );
one.data.MF.1$uniqueness <- as.numeric( as.character(one.data.MF.1$uniqueness) );
one.data.MF.1$dispensability <- as.numeric( as.character(one.data.MF.1$dispensability) );

# --------------------------------------------------------------------------
# Names of the axes, sizes of the numbers and letters, names of the columns,
# etc. can be changed below

P_1_MF <- ggplot( data = one.data.MF.1 );
pla3.cc = rev(brewer.pal(9, "OrRd"))[2:5]

P_1_MF = P_1_MF + geom_point( aes( plot_X, plot_Y, colour = log10_p_value, size = plot_size), alpha = I(0.6) ) 
#P_1_MF = P_1_MF + scale_colour_distiller(palette="Greys")
#P_1_MF = P_1_MF + scale_color_viridis(option = "F")
P_1_MF = P_1_MF + scale_color_gradientn(colors = pla3.cc)
P_1_MF = P_1_MF + scale_size(limits=c(0,5),range=c(1, 15),guide="none")
P_1_MF =  P_1_MF + xlim(-10,10)
P_1_MF = P_1_MF + ylim(-10,10)
P_1_MF = P_1_MF  + theme_bw()
P_1_MF = P_1_MF + labs (y = "semantic space y", x = "semantic space x")
P_1_MF_withlegend = P_1_MF + theme(legend.position=c(0.89,0.2),legend.text=element_text(size=8),legend.box="horizontal")


## LABEL THE CIRCLES
terms.P_1_MF= one.data.MF.1 [ one.data.MF.1$description %in% c("ion gated channel activity",
                                                               "transmitter-gated channel activity"), ]

terms.P_1_MFb= one.data.MF.1 [ one.data.MF.1$description %in% c("transporter activity","beta-amyloid binding",
                                                                "molecular transducer activity",
                                                                "glutamate binding",
                                                                "dopamine binding",
                                                                "sodium ion transmembrane transporter activity"), ]



## CHANGE FONT OF THE LABELS

P_1_MF = P_1_MF + geom_text( data = terms.P_1_MF, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 4)
P_1_MF = P_1_MF + geom_text( data = terms.P_1_MFb, aes(plot_X, plot_Y, label = description), 
                             colour = I(alpha("black", 0.85)), 
                             size = 4,fontface = "bold") 

## ADD TITLE ETC TO PLOT

P_1_MF = P_1_MF + ggtitle("Gene Ontology Molecular Function Threshold >/= 1 TPM")
P_1_MF = P_1_MF + theme(plot.title = element_text(hjust=0.5,vjust=1,size=12,face="bold"));P_1_MF

P_1_MF;

## SAVE PLOT

save_plot("./Plots/Human/GENE EXPRESSION ANALYSIS/P_1_MF.jpeg",
          P_1_MF ,base_height= 5.5 ,base_aspect_ratio = 1.5) 


################################################################################


########################### TPM >1 CC : ######################################

library( ggplot2 );
library( scales );
library(RColorBrewer);
# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","log10_p_value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0030054","cell junction", 6.490,-0.663,-6.120, 3.080,-5.0408,0.934,0.000),
                     c("GO:0030424","axon", 2.098, 5.625,-1.020, 2.590,-3.1271,0.634,0.000),
                     c("GO:0033010","paranodal junction", 0.027,-3.371,-5.059, 0.778,-1.3804,0.930,0.000),
                     c("GO:0045202","synapse", 4.392, 2.587,-6.346, 2.910,-2.6962,0.933,0.000),
                     c("GO:0098982","GABA-ergic synapse", 0.005, 4.268, 5.999, 0.301,-4.8484,0.688,0.000),
                     c("GO:1902495","transmembrane transporter complex", 1.736,-5.085, 3.861, 2.508,-8.8705,0.637,0.000),
                     c("GO:0009986","cell surface", 4.132,-5.832,-3.220, 2.884,-1.3804,0.930,0.003),
                     c("GO:1902711","GABA-A receptor complex", 0.092,-6.151, 1.414, 1.255,-3.9083,0.806,0.173),
                     c("GO:0044459","plasma membrane part",14.164,-3.164, 5.596, 3.418,-5.1175,0.762,0.211),
                     c("GO:0045211","postsynaptic membrane", 1.190, 0.902, 4.718, 2.344,-4.0456,0.410,0.419),
                     c("GO:0032809","neuronal cell body membrane", 0.103, 2.165, 0.204, 1.301,-1.7453,0.583,0.458),
                     c("GO:0044224","juxtaparanode region of axon", 0.054, 5.670,-1.670, 1.041,-2.4967,0.653,0.495),
                     c("GO:0044305","calyx of Held", 0.022, 4.509, 2.435, 0.699,-1.6159,0.521,0.496),
                     c("GO:0099061","integral component of postsynaptic density membrane", 0.016, 1.145, 5.100, 0.602,-3.8615,0.495,0.621),
                     c("GO:0098978","glutamatergic synapse", 0.994, 3.253, 6.170, 2.265,-3.3317,0.593,0.661),
                     c("GO:0098984","neuron to neuron synapse", 1.082, 3.530, 5.837, 2.303,-2.9139,0.591,0.667),
                     c("GO:0099056","integral component of presynaptic membrane", 0.045, 2.443, 2.855, 0.954,-2.9335,0.424,0.673),
                     c("GO:0098802","plasma membrane receptor complex", 0.919,-3.210, 2.760, 2.233,-1.9070,0.641,0.680),
                     c("GO:0036477","somatodendritic compartment", 3.548, 5.132,-1.359, 2.818,-2.9335,0.642,0.684),
                     c("GO:0017146","NMDA selective glutamate receptor complex", 0.054,-3.728, 2.912, 1.041,-2.8108,0.571,0.684));

## CHANGE FONT OF THE LABELS

P_1_CC = P_1_CC + geom_text( data = terms.P_1_CC, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 4)
P_1_CC = P_1_CC + geom_text( data = terms.P_1_CCb, aes(plot_X, plot_Y, label = description), 
                             colour = I(alpha("black", 0.85)), 
                             size = 4,fontface = "bold") 


## ADD TITLE ETC TO PLOT

P_1_CC = P_1_CC + ggtitle("Gene Ontology Cellular Components Threshold >/= 1 TPM")
P_1_CC = P_1_CC + theme(plot.title = element_text(hjust=0.5,vjust=1,size=12,face="bold"));P_1_CC


P_1_CC;

## SAVE PLOT

save_plot("./Plots/Human/GENE EXPRESSION ANALYSIS/P_1_CC.jpeg",
          P_1_CC ,base_height= 5.5 ,base_aspect_ratio = 1.5) 

#error message:
# While parsing your data, warning(s) were encountered:
# Go term 99240 was not found in the current version of the GeneOntology, dated 22:12:2016 16:59. GO term will be skipped.
# Go term 98691 was not found in the current version of the GeneOntology, dated 22:12:2016 16:59. GO term will be skipped.

##########################################################################################################

########################### TPM >1 BP : ##################################################################

revigo.names <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","log10_p_value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0098660","inorganic ion transmembrane transport", 4.409,-5.885,-2.085, 2.884,-2.0718,0.677,0.000),
                     c("GO:0099505","regulation of presynaptic membrane potential", 0.219,-6.323, 4.547, 1.591,-1.8734,0.811,0.000),
                     c("GO:0099537","trans-synaptic signaling", 3.485,-2.336, 5.183, 2.782,-6.6624,0.702,0.000),
                     c("GO:0071420","cellular response to histamine", 0.046, 2.637, 0.763, 0.954,-2.3017,0.905,0.002),
                     c("GO:0030534","adult behavior", 0.819, 1.782,-5.937, 2.155,-1.8734,0.721,0.013),
                     c("GO:0019233","sensory perception of pain", 0.537, 5.928,-1.660, 1.973,-1.7194,0.809,0.104),
                     c("GO:0022010","central nervous system myelination", 0.087, 4.965, 3.621, 1.204,-1.8734,0.837,0.115),
                     c("GO:0035235","ionotropic glutamate receptor signaling pathway", 0.150, 0.579, 6.923, 1.431,-4.1264,0.795,0.130),
                     c("GO:0015844","monoamine transport", 0.427,-3.734,-2.328, 1.875,-1.3473,0.803,0.203),
                     c("GO:0097480","establishment of synaptic vesicle localization", 0.767,-5.573,-3.918, 2.127,-1.8734,0.807,0.211),
                     c("GO:0099003","vesicle-mediated transport in synapse", 0.767,-4.385,-4.606, 2.127,-1.8734,0.817,0.218),
                     c("GO:1904862","inhibitory synapse assembly", 0.017, 4.002, 5.099, 0.602,-1.5010,0.849,0.228),
                     c("GO:0045163","clustering of voltage-gated potassium channels", 0.017, 4.945, 4.417, 0.602,-1.3351,0.845,0.288),
                     c("GO:0060134","prepulse inhibition", 0.075, 6.104,-1.052, 1.146,-1.5010,0.792,0.387),
                     c("GO:0031644","regulation of neurological system process", 0.387, 5.335,-2.471, 1.833,-1.3902,0.786,0.446),
                     c("GO:0098976","excitatory chemical synaptic transmission", 0.006,-2.570, 5.738, 0.301,-1.8734,0.778,0.494),
                     c("GO:0042756","drinking behavior", 0.058, 1.130,-6.494, 1.041,-1.3351,0.796,0.540),
                     c("GO:0042749","regulation of circadian sleep/wake cycle", 0.127, 2.473,-5.519, 1.362,-1.6193,0.737,0.578),
                     c("GO:0022898","regulation of transmembrane transporter activity", 1.166,-6.323,-2.536, 2.307,-1.5010,0.717,0.609),
                     c("GO:0060080","inhibitory postsynaptic potential", 0.081,-3.405, 4.967, 1.176,-1.8734,0.662,0.618),
                     c("GO:0023061","signal release", 2.441,-4.061, 1.431, 2.627,-1.5978,0.639,0.619),
                     c("GO:0072511","divalent inorganic cation transport", 2.499,-6.188,-1.526, 2.637,-1.5010,0.704,0.638),
                     c("GO:0097553","calcium ion transmembrane import into cytosol", 0.629,-5.638,-0.448, 2.041,-1.5010,0.609,0.645),
                     c("GO:0048169","regulation of long-term neuronal synaptic plasticity", 0.133,-3.791, 4.550, 1.380,-1.5010,0.709,0.649),
                     c("GO:0007628","adult walking behavior", 0.162, 1.589,-6.345, 1.462,-1.8175,0.740,0.653));

one.data.BP_1 = data.frame(revigo.data);
names(one.data.BP_1) = revigo.names;
one.data.BP_1$plot_X = as.numeric( as.character(one.data.BP_1$plot_X) );
one.data.BP_1$plot_Y = as.numeric( as.character(one.data.BP_1$plot_Y) );
one.data.BP_1$plot_size = as.numeric( as.character(one.data.BP_1$plot_size) );
one.data.BP_1$log10_p_value = as.numeric( as.character(one.data.BP_1$log10_p_value));
one.data.BP_1$frequency = as.numeric( as.character(one.data.BP_1$frequency) );
one.data.BP_1$uniqueness = as.numeric( as.character(one.data.BP_1$uniqueness) );
one.data.BP_1$dispensability = as.numeric( as.character(one.data.BP_1$dispensability) );

####################################################################################################################

p_1_BP = ggplot( data = one.data.BP_1);

pla3 = rev(brewer.pal(9, "OrRd"))[2:5]

p_1_BP = p_1_BP + geom_point( aes( plot_X, plot_Y, colour = log10_p_value, size = plot_size), alpha = I(0.6) ) 
#p_1_BP = p_1_BP + scale_colour_distiller(palette="Greys")
#p_1_BP = p_1_BP + scale_color_viridis(option = "F")
p_1_BP = p_1_BP + scale_color_gradientn(colors = pla3)
p_1_BP = p_1_BP + scale_size(limits=c(0,5),range=c(1, 15),guide="none")
p_1_BP =  p_1_BP + xlim(-10,10)
p_1_BP = p_1_BP + ylim(-10,10)
p_1_BP = p_1_BP  + theme_bw()
p_1_BP = p_1_BP + labs (y = "semantic space y", x = "semantic space x")
p_1_BP_withlegend = p_1_BP + theme(legend.position=c(0.89,0.2),legend.text=element_text(size=8),legend.box="horizontal")

## LABEL THE CIRCLES

terms.p_1_BP= one.data.BP_1 [ one.data.BP_1$description %in% c("trans-synaptic signaling",
                                                               "clustering of voltage-gated potassium channels",
                                                               "vesicle-mediated transport in synapse",
                                                               "ionotropic glutamate receptor signaling pathway"), ]



terms.p_1_BPb= one.data.BP_1 [ one.data.BP_1$description %in% c("signal release","adult behavior","adult walking behavior",
                                                                "central nervous system myelination" , 
                                                                "sensory perception of pain",
                                                                "monoamine transport",
                                                                "inorganic ion transmembrane transport",
                                                                "inhibitory synapse assembly",
                                                                "regulation of circadian sleep/wake cycle"), ]
## CHANGE FONT OF THE LABELS

p_1_BP = p_1_BP + geom_text( data = terms.p_1_BP, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 4)
p_1_BP = p_1_BP + geom_text( data = terms.p_1_BPb, aes(plot_X, plot_Y, label = description), 
                             colour = I(alpha("black", 0.85)), 
                             size = 4,fontface = "bold") 

## ADD TITLE ETC TO PLOT

p_1_BP = p_1_BP + ggtitle("Gene Ontology Biological Process Threshold >/= 1 TPM")
p_1_BP = p_1_BP + theme(plot.title = element_text(hjust=0.5,vjust=1,size=12,face="bold"));p_1_BP

## SAVE PLOT

save_plot("./Plots/Human/GENE EXPRESSION ANALYSIS/p_1_BP.jpeg",
          p_1_BP ,base_height= 5.5 ,base_aspect_ratio = 1.5) 



########################### TPM >1 MF : ######################################


# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","log10_p_value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0001540","beta-amyloid binding", 0.191,-2.973,-4.642, 1.531,-3.4862,0.897,0.000),
                     c("GO:0005215","transporter activity", 7.829,-2.902, 4.612, 3.132,-7.5998,0.905,0.000),
                     c("GO:0022835","transmitter-gated channel activity", 0.283, 6.706, 0.419, 1.699,-12.0761,0.203,0.000),
                     c("GO:0060089","molecular transducer activity", 9.545, 1.771,-7.358, 3.218,-3.2614,0.907,0.000),
                     c("GO:0016595","glutamate binding", 0.058,-4.681,-0.001, 1.041,-2.7378,0.888,0.004),
                     c("GO:0035240","dopamine binding", 0.075, 1.881, 7.254, 1.146,-2.5976,0.888,0.090),
                     c("GO:0015081","sodium ion transmembrane transporter activity", 0.884, 5.756,-0.537, 2.188,-2.7135,0.347,0.496),
                     c("GO:0015108","chloride transmembrane transporter activity", 0.532, 6.563,-0.909, 1.968,-2.6075,0.369,0.564),
                     c("GO:0022839","ion gated channel activity", 0.260, 6.029, 0.685, 1.663,-6.2118,0.318,0.658));

one.data.MF.1<- data.frame(revigo.data);
names(one.data.MF.1) <- revigo.names;
one.data.MF.1 <- one.data.MF.1 [(one.data.MF.1$plot_X != "null" & one.data.MF.1$plot_Y != "null"), ];
one.data.MF.1$plot_X <- as.numeric( as.character(one.data.MF.1$plot_X) );
one.data.MF.1$plot_Y <- as.numeric( as.character(one.data.MF.1$plot_Y) );
one.data.MF.1$plot_size <- as.numeric( as.character(one.data.MF.1$plot_size) );
one.data.MF.1$log10_p_value <- as.numeric( as.character(one.data.MF.1$log10_p_value) );
one.data.MF.1$frequency <- as.numeric( as.character(one.data.MF.1$frequency) );
one.data.MF.1$uniqueness <- as.numeric( as.character(one.data.MF.1$uniqueness) );
one.data.MF.1$dispensability <- as.numeric( as.character(one.data.MF.1$dispensability) );

# --------------------------------------------------------------------------
# Names of the axes, sizes of the numbers and letters, names of the columns,
# etc. can be changed below

P_1_MF <- ggplot( data = one.data.MF.1 );
pla3.cc = rev(brewer.pal(9, "OrRd"))[2:5]

P_1_MF = P_1_MF + geom_point( aes( plot_X, plot_Y, colour = log10_p_value, size = plot_size), alpha = I(0.6) ) 
#P_1_MF = P_1_MF + scale_colour_distiller(palette="Greys")
#P_1_MF = P_1_MF + scale_color_viridis(option = "F")
P_1_MF = P_1_MF + scale_color_gradientn(colors = pla3.cc)
P_1_MF = P_1_MF + scale_size(limits=c(0,5),range=c(1, 15),guide="none")
P_1_MF =  P_1_MF + xlim(-8.5,9.5)
P_1_MF = P_1_MF + ylim(-7,7)
P_1_MF = P_1_MF  + theme_bw()
P_1_MF = P_1_MF + labs (y = "semantic space y", x = "semantic space x")
P_1_MF_withlegend = P_1_MF + theme(legend.position=c(0.89,0.2),legend.text=element_text(size=8),legend.box="horizontal")


## LABEL THE CIRCLES
terms.P_1_MF= one.data.MF.1 [ one.data.MF.1$description %in% c("ion gated channel activity",
                                                               "transmitter-gated channel activity"), ]

terms.P_1_MFb= one.data.MF.1 [ one.data.MF.1$description %in% c("transporter activity","beta-amyloid binding",
                                                                "molecular transducer activity",
                                                                "glutamate binding",
                                                                "dopamine binding",
                                                                "sodium ion transmembrane transporter activity"), ]



## CHANGE FONT OF THE LABELS

P_1_MF = P_1_MF + geom_text( data = terms.P_1_MF, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 4)
P_1_MF = P_1_MF + geom_text( data = terms.P_1_MFb, aes(plot_X, plot_Y, label = description), 
                             colour = I(alpha("black", 0.85)), 
                             size = 4,fontface = "bold") 

## ADD TITLE ETC TO PLOT

P_1_MF = P_1_MF + ggtitle("Gene Ontology Molecular Function Threshold >/= 1 TPM")
P_1_MF = P_1_MF + theme(plot.title = element_text(hjust=0.5,vjust=1,size=12,face="bold"));P_1_MF

P_1_MF;

## SAVE PLOT

save_plot("./Plots/Human/GENE EXPRESSION ANALYSIS/P_1_MF.jpeg",
          P_1_MF ,base_height= 5.5 ,base_aspect_ratio = 1.5) 


################################################################################



########################### TPM >0.1 MF : ######################################


# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","log10_p_value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0001227","transcriptional repressor activity, RNA polymerase II transcription regulatory region sequence-specific binding", 1.057,-0.900,-6.819, 2.265,-1.4383,0.909,0.000),
                     c("GO:0005215","transporter activity", 7.829,-2.273, 5.849, 3.132,-2.1218,0.949,0.000),
                     c("GO:0008310","single-stranded DNA 3'-5' exodeoxyribonuclease activity", 0.029,-2.519, 0.330, 0.778,-1.3425,0.935,0.000),
                     c("GO:0022835","transmitter-gated channel activity", 0.283, 5.707, 0.169, 1.699,-8.2762,0.345,0.000),
                     c("GO:0042166","acetylcholine binding", 0.098,-5.966, 1.956, 1.255,-1.6003,0.931,0.000),
                     c("GO:0060089","molecular transducer activity", 9.545, 3.279,-7.323, 3.218,-2.1932,0.950,0.000),
                     c("GO:0001515","opioid peptide activity", 0.023,-4.923,-4.293, 0.699,-1.3425,0.917,0.004),
                     c("GO:0001162","RNA polymerase II intronic transcription regulatory region sequence-specific DNA binding", 0.046, 4.572, 6.769, 0.954,-1.3425,0.944,0.004),
                     c("GO:0019894","kinesin binding", 0.220, 0.309, 6.895, 1.591,-1.3425,0.941,0.023),
                     c("GO:0016594","glycine binding", 0.087,-4.795, 3.785, 1.204,-1.5362,0.896,0.152),
                     c("GO:0003918","DNA topoisomerase type II (ATP-hydrolyzing) activity", 0.017,-1.749,-1.354, 0.602,-1.3425,0.935,0.190),
                     c("GO:0031771","type 1 hypocretin receptor binding", 0.006,-5.647,-3.140, 0.301,-1.3425,0.877,0.245),
                     c("GO:0008506","sucrose:proton symporter activity", 0.017, 6.268,-2.615, 0.602,-1.3425,0.590,0.362),
                     c("GO:0008066","glutamate receptor activity", 0.156, 4.565, 3.491, 1.447,-1.3425,0.728,0.390),
                     c("GO:0004993","G-protein coupled serotonin receptor activity", 0.162, 3.587, 3.367, 1.462,-1.3425,0.724,0.391),
                     c("GO:0015108","chloride transmembrane transporter activity", 0.532, 6.748,-1.830, 1.968,-1.4418,0.546,0.474),
                     c("GO:0022839","ion gated channel activity", 0.260, 6.539,-1.269, 1.663,-4.2238,0.516,0.658),
                     c("GO:0003700","transcription factor activity, sequence-specific DNA binding", 7.095,-0.525,-6.856, 3.090,-1.3425,0.912,0.683),
                     c("GO:0015157","oligosaccharide transmembrane transporter activity", 0.017, 5.868,-2.536, 0.602,-1.3425,0.617,0.691));

one.data.MF_0.1 <- data.frame(revigo.data);
names(one.data.MF_0.1) <- revigo.names;
one.data.MF_0.1 <- one.data.MF_0.1 [(one.data.MF_0.1$plot_X != "null" & one.data.MF_0.1$plot_Y != "null"), ];
one.data.MF_0.1$plot_X <- as.numeric( as.character(one.data.MF_0.1$plot_X) );
one.data.MF_0.1$plot_Y <- as.numeric( as.character(one.data.MF_0.1$plot_Y) );
one.data.MF_0.1$plot_size <- as.numeric( as.character(one.data.MF_0.1$plot_size) );
one.data.MF_0.1$log10_p_value <- as.numeric( as.character(one.data.MF_0.1$log10_p_value) );
one.data.MF_0.1$frequency <- as.numeric( as.character(one.data.MF_0.1$frequency) );
one.data.MF_0.1$uniqueness <- as.numeric( as.character(one.data.MF_0.1$uniqueness) );
one.data.MF_0.1$dispensability <- as.numeric( as.character(one.data.MF_0.1$dispensability) );
#head(one.data);


# --------------------------------------------------------------------------
# Names of the axes, sizes of the numbers and letters, names of the columns,
# etc. can be changed below

P_0.1_MF <- ggplot( data = one.data.MF_0.1 );
pla3.cc = rev(brewer.pal(9, "OrRd"))[2:5]

P_0.1_MF = P_0.1_MF + geom_point( aes( plot_X, plot_Y, colour = log10_p_value, size = plot_size), alpha = I(0.6) ) 
#P_0.1_MF = P_0.1_MF + scale_colour_distiller(palette="Greys")
#P_0.1_MF = P_0.1_MF + scale_color_viridis(option = "F")
P_0.1_MF = P_0.1_MF + scale_color_gradientn(colors = pla3.cc)
P_0.1_MF = P_0.1_MF + scale_size(limits=c(0,5),range=c(1, 15),guide="none")
P_0.1_MF =  P_0.1_MF + xlim(-10,10)
P_0.1_MF = P_0.1_MF + ylim(-10,10)
P_0.1_MF = P_0.1_MF  + theme_bw()
P_0.1_MF = P_0.1_MF + labs (y = "semantic space y", x = "semantic space x")
P_0.1_MF_withlegend = P_0.1_MF + theme(legend.position=c(0.89,0.2),legend.text=element_text(size=8),legend.box="horizontal")



## LABEL THE CIRCLES
terms.P_0.1_MF= one.data.MF_0.1 [ one.data.MF_0.1$description %in% c("transmitter-gated channel activity",
                                                                     "ion gated channel activity",
                                                                     "molecular transducer activity",
                                                                     "transporter activity","glutamate receptor activity",
                                                                     "single-stranded DNA 3'-5' exodeoxyribonuclease activity"), ]

terms.P_0.1_MFb= one.data.MF_0.1[ one.data.MF_0.1$description %in% c("transcription factor activity, sequence-specific DNA binding",
                                                                     "acetylcholine binding",
                                                                     "opioid peptide activity"), ]


## CHANGE FONT OF THE LABELS

P_0.1_MF = P_0.1_MF + geom_text( data = terms.P_0.1_MF, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 4)
P_0.1_MF = P_0.1_MF + geom_text( data = terms.P_0.1_MFb, aes(plot_X, plot_Y, label = description), 
                                 colour = I(alpha("black", 0.85)), 
                                 size = 4,fontface = "bold") 

## ADD TITLE ETC TO PLOT

P_0.1_MF = P_0.1_MF + ggtitle("Gene Ontology Molecular Function Threshold >/= 0.1 TPM")
P_0.1_MF = P_0.1_MF + theme(plot.title = element_text(hjust=0.5,vjust=1,size=12,face="bold"));P_0.1_MF

P_0.1_MF;

## SAVE PLOT
library(cowplot)

save_plot("./Plots/Human/GENE EXPRESSION ANALYSIS/P_0.1_MF.jpeg",
          P_0.1_MF ,base_height= 5.5 ,base_aspect_ratio = 1.5) 




################################################################################



########################### TPM >0.1 CC : ######################################


library( ggplot2 );
library( scales );

# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","log10_p_value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0030054","cell junction", 6.490,-2.142,-4.491, 3.080,-2.5900,0.893,0.000),
                     c("GO:0036477","somatodendritic compartment", 3.548, 3.205,-3.831, 2.818,-4.1229,0.839,0.000),
                     c("GO:1902495","transmembrane transporter complex", 1.736,-4.838, 5.273, 2.508,-4.3592,0.538,0.000),
                     c("GO:0099061","integral component of postsynaptic density membrane", 0.016, 2.368, 5.397, 0.602,-4.0278,0.335,0.133),
                     c("GO:1902711","GABA-A receptor complex", 0.092,-5.815, 1.806, 1.255,-3.5113,0.667,0.173),
                     c("GO:0098590","plasma membrane region", 5.165,-1.071, 6.800, 2.980,-1.6843,0.555,0.287),
                     c("GO:0098982","GABA-ergic synapse", 0.005, 5.788, 3.054, 0.301,-2.5900,0.592,0.322),
                     c("GO:0045211","postsynaptic membrane", 1.190, 2.161, 4.751, 2.344,-2.5900,0.302,0.621),
                     c("GO:0098978","glutamatergic synapse", 0.994, 5.391, 4.653, 2.265,-1.8450,0.493,0.661),
                     c("GO:0098984","neuron to neuron synapse", 1.082, 5.336, 4.086, 2.303,-1.8101,0.491,0.667),
                     c("GO:0099056","integral component of presynaptic membrane", 0.045, 2.413, 3.151, 0.954,-1.6843,0.350,0.673),
                     c("GO:0099634","postsynaptic specialization membrane", 0.049, 2.968, 4.928, 1.000,-1.3375,0.385,0.678),
                     c("GO:0098802","plasma membrane receptor complex", 0.919,-2.409, 4.162, 2.233,-1.7774,0.456,0.680),
                     c("GO:0017146","NMDA selective glutamate receptor complex", 0.054,-2.998, 4.484, 1.041,-1.9811,0.416,0.684));

one.data.CC_0.1 <- data.frame(revigo.data);
names(one.data.CC_0.1) <- revigo.names;
one.data.CC_0.1 <- one.data.CC_0.1 [(one.data.CC_0.1$plot_X != "null" & one.data.CC_0.1$plot_Y != "null"), ];
one.data.CC_0.1$plot_X <- as.numeric( as.character(one.data.CC_0.1$plot_X) );
one.data.CC_0.1$plot_Y <- as.numeric( as.character(one.data.CC_0.1$plot_Y) );
one.data.CC_0.1$plot_size <- as.numeric( as.character(one.data.CC_0.1$plot_size) );
one.data.CC_0.1$log10_p_value <- as.numeric( as.character(one.data.CC_0.1$log10_p_value) );
one.data.CC_0.1$frequency <- as.numeric( as.character(one.data.CC_0.1$frequency) );
one.data.CC_0.1$uniqueness <- as.numeric( as.character(one.data.CC_0.1$uniqueness) );
one.data.CC_0.1$dispensability <- as.numeric( as.character(one.data.CC_0.1$dispensability) );
#head(one.data.CC_0.1);


# --------------------------------------------------------------------------
# Names of the axes, sizes of the numbers and letters, names of the columns,
# etc. can be changed below

P_0.1_CC <- ggplot( data = one.data.CC_0.1 );
pla3.cc = rev(brewer.pal(9, "OrRd"))[2:5]

P_0.1_CC = P_0.1_CC + geom_point( aes( plot_X, plot_Y, colour = log10_p_value, size = plot_size), alpha = I(0.6) ) 
#P_0.1_CC = P_0.1_CC + scale_colour_distiller(palette="Greys")
#P_0.1_CC = P_0.1_CC + scale_color_viridis(option = "F")
P_0.1_CC = P_0.1_CC + scale_color_gradientn(colors = pla3.cc)
P_0.1_CC = P_0.1_CC + scale_size(limits=c(0,5),range=c(1, 15),guide="none")
P_0.1_CC =  P_0.1_CC + coord_cartesian(xlim=c(-10,10))
P_0.1_CC = P_0.1_CC + ylim=c(-10,10)
P_0.1_CC = P_0.1_CC  + theme_bw()
P_0.1_CC = P_0.1_CC + labs (y = "semantic space y", x = "semantic space x")
P_0.1_CC_withlegend = P_0.1_CC + theme(legend.position=c(0.89,0.2),legend.text=element_text(size=8),legend.box="horizontal")



## LABEL THE CIRCLES
terms.P_0.1_CC= one.data.CC_0.1 [ one.data.CC_0.1$description %in% c("transmembrane transporter complex",
                                                                     "cell junction",
                                                                     "somatodendritic compartment"), ]




terms.P_0.1_CCb= one.data.CC_0.1[ one.data.CC_0.1$description %in% c("plasma membrane region",
                                                                     "neuron to neuron synapse",
                                                                     "postsynaptic membrane",
                                                                     "plasma membrane receptor complex"), ]


## CHANGE FONT OF THE LABELS

P_0.1_CC = P_0.1_CC + geom_text( data = terms.P_0.1_CC, aes(plot_X, plot_Y, label = description), colour = I(alpha("black", 0.85)), size = 4)
P_0.1_CC = P_0.1_CC + geom_text( data = terms.P_0.1_CCb, aes(plot_X, plot_Y, label = description), 
                                 colour = I(alpha("black", 0.85)), 
                                 size = 4,fontface = "bold") 

## ADD TITLE ETC TO PLOT

P_0.1_CC = P_0.1_CC + ggtitle("Gene Ontology Cellular Component Threshold >/= 0.1 TPM")
P_0.1_CC = P_0.1_CC + theme(plot.title = element_text(hjust=0.5,vjust=1,size=12,face="bold"));P_0.1_CC

P_0.1_CC;

## SAVE PLOT
library(cowplot)

save_plot("./Plots/Human/GENE EXPRESSION ANALYSIS/P_0.1_CC.jpeg",
          P_0.1_CC ,base_height= 5.5 ,base_aspect_ratio = 1.5) 




################################################################################



########################### TPM >0.1 BP : ######################################


library( ggplot2 );
library( scales );


# --------------------------------------------------------------------------
# Here is your data from REVIGO. Scroll down for plot configuration options.

revigo.names <- c("term_ID","description","frequency_%","plot_X","plot_Y","plot_size","log10_p_value","uniqueness","dispensability");
revigo.data <- rbind(c("GO:0071420","cellular response to histamine", 0.046,-5.332,-1.426, 0.954,-2.1446,0.379,0.000),
                     c("GO:0099537","trans-synaptic signaling", 3.485, 3.934,-3.852, 2.782,-2.1446,0.285,0.002),
                     c("GO:0035235","ionotropic glutamate receptor signaling pathway", 0.150, 3.752, 3.760, 1.431,-1.6794,0.384,0.130),
                     c("GO:0035095","behavioral response to nicotine", 0.040,-5.387, 0.896, 0.903,-1.3351,0.394,0.376));

one.data.BP_0.1 <- data.frame(revigo.data);
names(one.data.BP_0.1) <- revigo.names;
one.data.BP_0.1 <- one.data.BP_0.1 [(one.data.BP_0.1$plot_X != "null" & one.data.BP_0.1$plot_Y != "null"), ];
one.data.BP_0.1$plot_X <- as.numeric( as.character(one.data.BP_0.1$plot_X) );
one.data.BP_0.1$plot_Y <- as.numeric( as.character(one.data.BP_0.1$plot_Y) );
one.data.BP_0.1$plot_size <- as.numeric( as.character(one.data.BP_0.1$plot_size) );
one.data.BP_0.1$log10_p_value <- as.numeric( as.character(one.data.BP_0.1$log10_p_value) );
one.data.BP_0.1$frequency <- as.numeric( as.character(one.data.BP_0.1$frequency) );
one.data.BP_0.1$uniqueness <- as.numeric( as.character(one.data.BP_0.1$uniqueness) );
one.data.BP_0.1$dispensability <- as.numeric( as.character(one.data.BP_0.1$dispensability) );
#head(one.data.BP_0.1);


# --------------------------------------------------------------------------
# Names of the axes, sizes of the numbers and letters, names of the columns,
# etc. can be changed below

P_0.1_BP<- ggplot( data = one.data.BP_0.1 );
pla3.cc = rev(brewer.pal(9, "OrRd"))[2:5]

P_0.1_BP = P_0.1_BP + geom_point( aes( plot_X, plot_Y, colour = log10_p_value, size = plot_size), alpha = I(0.6) ) 
#P_0.1_BP = P_0.1_BP + scale_colour_distiller(palette="Greys")
#P_0.1_BP = P_0.1_BP + scale_color_viridis(option = "F")
P_0.1_BP = P_0.1_BP + scale_color_gradientn(colors = pla3.cc)
P_0.1_BP = P_0.1_BP + scale_size(limits=c(0,5),range=c(1, 15),guide="none")
P_0.1_BP =  P_0.1_BP + xlim(-10,10)
P_0.1_BP = P_0.1_BP + ylim(-10,10)
P_0.1_BP = P_0.1_BP  + theme_bw()
P_0.1_BP = P_0.1_BP + labs (y = "semantic space y", x = "semantic space x")
P_0.1_BP_withlegend = P_0.1_BP + theme(legend.position=c(0.89,0.2),legend.text=element_text(size=8),legend.box="horizontal")



## LABEL THE CIRCLES


terms.P_0.1_BPb= one.data.BP_0.1[ one.data.BP_0.1$description %in% c("ionotropic glutamate receptor signaling pathway",
                                                                     "behavioral response to nicotine",
                                                                     "trans-synaptic signaling",
                                                                     "cellular response to histamine"), ]


## CHANGE FONT OF THE LABELS

P_0.1_BP = P_0.1_BP + geom_text( data = terms.P_0.1_BPb, aes(plot_X, plot_Y, label = description), 
                                 colour = I(alpha("black", 0.85)), 
                                 size = 4,fontface = "bold") 

## ADD TITLE ETC TO PLOT

P_0.1_BP = P_0.1_BP + ggtitle("Gene Ontology Cellular Component Threshold >/= 0.1 TPM")
P_0.1_BP = P_0.1_BP + theme(plot.title = element_text(hjust=0.5,vjust=1,size=12,face="bold"));P_0.1_BP

P_0.1_BP;

## SAVE PLOT
library(cowplot)

save_plot("./Plots/Human/GENE EXPRESSION ANALYSIS/P_0.1_BP.jpeg",
          P_0.1_BP ,base_height= 5.5 ,base_aspect_ratio = 1.5) 




################################################################################










