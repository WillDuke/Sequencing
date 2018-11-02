lapply(c("tidyverse", "magrittr", "ComplexHeatmap", "devtools"),
       require, character.only = TRUE)
install.packages("data.table")
#load in data, file = alterations_across_samples.tsv
#file: alterations_across_samples.tsv
alt_data <- read.delim(file.choose(),
                       header = TRUE, stringsAsFactors = FALSE)

#save dimnames
alt_rownames <- alt_data$Sample.ID
alt_colnames <- colnames(alt_data)[-c(1:4)]
alt_colnames[35:36] <- c("mut_num", "driver")

#remove extra columns, set Sample ID column to rownames,
#remove mutation details, convert to matrix,
#order by # concomittant mutations and driver type, transpose,
#convert to matrix, carry along rownames (should just not use them eek),
#and clean up rownames
alt_data <- alt_data %>%
  select(-Study.ID, -Patient.ID, -Altered) %>% column_to_rownames("Sample.ID") %>%
  as.matrix() %>%
  str_replace("MUT:[\\s\\S]*$", "MUT") %>%
  matrix(nrow = 1002, ncol = 34) %>%
  as.data.frame(row.names = alt_rownames) %>%
  mutate(mut_num = apply(alt_data, 1, function(x){sum(str_count(x, "MUT"))}),
         driver = ifelse(V1 == "MUT", "JAK2",
                  ifelse(V2 == "MUT", "CALR",
                  ifelse(V3 == "MUT", "MPL", "")))) %>%
  mutate(driver = factor(driver, levels = c("JAK2", "CALR", "MPL"), ordered = TRUE)) %>%
  set_colnames(alt_colnames) %>%
  set_rownames(alt_rownames) %>%
  rownames_to_column("Sample.IDs") %>%
  filter(driver != "") %>% arrange(mut_num, driver) %>%
  select(-mut_num, -driver) %>%
  column_to_rownames("Sample.IDs") %>%
  as.matrix() %>%
  t() %>%
  set_rownames(gsub("\\..*","",rownames(.)))

#create a list with column order and rendering details for oncoprint,
#and other info about the data (number of samples and genes, number of patients)
supp <- list(
  alter_fun = list(
    background = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
    },
    MUT = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "blue", col = NA))
    }
  ),
  col = c("MUT" = "blue"),
  order = list(column_order = colnames(alt_data), row_order = rownames(alt_data)),
  dims = dim(alt_data),
  patient_n = length(unique(gsub("(.*)S.*","\\1", colnames(alt_data)))))

#create pdf of oncoprint
pdf(file = "/Users/Will/Desktop/plot3.pdf", bg = "transparent")
ht <- oncoPrint(alt_data,
                alter_fun = supp$alter_fun, col = supp$col,
                remove_empty_columns = TRUE,
                column_order = supp$order$column_order,
                row_order = supp$order$row_order,
                column_title = "cBioPortal PTs with mutant JAK2, CALR, or MPL",
                heatmap_legend_param = list(title = "Alterations", at = c("MUT"),
                                            labels = c("Mutation")))
draw(ht, heatmap_legend_side = "bottom")
dev.off()

#####
#CALR first implementation

#doing it again with subset of CALR and separating CALR types

CALR_dat <- read.delim(file.choose(),
                       header = TRUE, stringsAsFactors = FALSE)
CALR_dat2 <- CALR_dat %>% select(-Study.ID, -Patient.ID, -Altered) %>%
  filter(CALR..MUT.FUSION. != "") %>% select(Sample.ID, CALR..MUT.FUSION.) %>%
  rename(CALR = CALR..MUT.FUSION.)
CALR_sub <- as.data.frame(alt_data) %>%
  rownames_to_column("Genes") %>%
  t()
colnames(CALR_sub) <- CALR_sub[1,]
CALR_sub <- CALR_sub[-1,]
rest_of_data <- CALR_sub %>% as.data.frame() %>% rownames_to_column("Sample.ID") %>% filter(CALR == "MUT") %>%
  select(-JAK2, - MPL, - CALR)

CALR_detail <- merge(CALR_dat2, rest_of_data, by = "Sample.ID") %>%
  mutate(CALR_type = ifelse(CALR == "MUT: K385fs*;", "Type_II", "Type_I")) %>%
  select(-CALR) %>% rename(CALR = CALR_type) %>% select(Sample.ID, CALR, 2:33) %>% t()
colnames(CALR_detail) <- as.character(CALR_detail[1,])
CALR_detail <- CALR_detail[-1,]
CALR_detail <- as.matrix(CALR_detail)

CALR_ordered <- CALR_detail %>%
  as.data.frame() %>%
  rownames_to_column("IDs")%>%
  t() %>%
  as.data.frame() %>%
  set_colnames(rownames(CALR_detail)) %>%
  rownames_to_column("Sample.ID") %>%
  slice(-1) %>%
  mutate(mut_num = apply(., MARGIN = 1, FUN =
                           function(x){sum(str_count(x, "(MUT|Type_I|Type_II)"))})) %>%
  arrange(CALR, mut_num) %>% t() %>% as.data.frame() %>%
  set_colnames(colnames(CALR_detail)) %>% rownames_to_column("IDs") %>% slice(c(-1, -34)) %>%
  column_to_rownames("IDs")


#found type classifications from COSMIC
write.csv(levels(as.factor(CALR_detail$CALR)),
          file = "/Users/Will/R_Files/Projects/Sequencing/CALR_mutations.csv")
rownames(CALR_ordered)

CALR_supp <- list(
  alter_fun = list(
    Type_I = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "blue", col = NA))
    },
    Type_II = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "#008000", col = NA))
    },
    MUT = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "red", col = NA))
    },
    background = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
    }
  ),
  col = c("Type_II" = "#008000", "Type_I" = "blue", "MUT" = "red"),
  order = list(column_order = colnames(CALR_ordered), row_order = rownames(CALR_ordered)),
  row_reorder <- c("CALR", "TET2", "DNMT3A", "IDH1", "IDH2", 
                   "ASXL1", "EZH2", "BCOR", "BCORL1", "SRSF2", "SF3B1", 
                   "U2AF1", "ZRSR2", "SH2B3", "CBL", "NRAS", " KRAS", "PTPN11", 
                   "GNB1", "GNAS", "FLT3", "RIT1", "KIT", "BRAF", "TP53", 
                   "RUNX1", "ETV6", "GATA2", "SMC1", "STAG2", "RAD21", "WT1"),
  dims = dim(CALR_ordered),
  patient_n = length(unique(gsub("(.*)S.*","\\1", colnames(CALR_ordered)))))

save(CALR_detail, file = "/Users/Will/R_Files/Projects/Sequencing/CALR_detail.rda")
save(CALR_ordered, file = "/Users/Will/R_Files/Projects/Sequencing/CALR_ordered.rda")

#create pdf of oncoprint
pdf(file = "/Users/Will/Desktop/plot6.pdf", bg = "transparent")
ht <- oncoPrint(CALR_ordered,
                alter_fun = CALR_supp$alter_fun, col = CALR_supp$col,
                remove_empty_columns = FALSE,
                column_order = CALR_supp$order$column_order,
                row_order = CALR_supp$order$row_reorder,
                column_title = "cBioPortal PTs with mutant CALR Types I and II",
                heatmap_legend_param = list(title = "Alterations", at = c("Type_I", "Type_II", "MUT"),
                            labels = c("CALR Type I", "CALR Type II", "Concomitant Mutation"),
                            nrow = 1, title_position = "leftcenter"))
draw(ht, heatmap_legend_side = "bottom")
dev.off()

#without explicit ordering
ht_test <- oncoPrint(CALR_ordered, alter_fun = alter_fun, col = CALR_supp$col)
draw(ht_test, heatmap_legend_side = "bottom")
#####

#recreate with different colors for JAK2, CALR, and MPL

color_diff <- alt_data %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate(JAK2_col = str_replace(alt_data[1,], "MUT", "JAK2"),
         CALR_col = str_replace(alt_data[2,], "MUT", "CALR"),
         MPL_col = str_replace(alt_data[3,], "MUT", "MPL")) %>% 
  select(-JAK2, -CALR, -MPL) %>% 
  rename(JAK2 = JAK2_col, CALR = CALR_col, MPL = MPL_col) %>% 
  select(JAK2, CALR, MPL, 1:31) %>% 
  as.matrix() %>% 
  t() %>% 
  as.data.frame() %>% 
  set_colnames(colnames(alt_data)) 

#create supp

index <- c(1:3,15:26,14,4:13, 34, 30:32, 27:29, 33)

col_supp <- list(
  alter_fun = list(
    JAK2 = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#E41A1C", col = NA))
    },
    CALR = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#377EB8", col = NA))
    },
    MPL = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#4DAF4A", col = NA))
    },
    MUT = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "#984EA3", col = NA))
    },
    background = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#F7F7F7", col = NA))
    }
  ),
  col = c("JAK2" = "#E41A1C", "CALR" = "#377EB8", "MPL" = "#4DAF4A","MUT" = "#984EA3"),
  order = list(column_order = colnames(color_diff), row_order = rownames(color_diff)[index]),
  row_reorder <- rownames(color_diff)[index],
  dims = dim(color_diff),
  patient_n = length(unique(gsub("(.*)S.*","\\1", colnames(color_diff)))))

#create pdf of oncoprint
pdf(file = "./Figures/plot2.pdf", bg = "transparent")
ht <- oncoPrint(color_diff,
                alter_fun = col_supp$alter_fun, col = col_supp$col,
                remove_empty_columns = TRUE,
                column_order = col_supp$order$column_order,
                row_order = col_supp$order$row_order,
                column_title = "cBioPortal PTs with mutant JAK2, CALR, or MPL",
                heatmap_legend_param = list(title = "Alterations", at = c("MUT"),
                                            labels = c("Mutation")))
draw(ht, show_heatmap_legend = FALSE)
dev.off()

library("RColorBrewer")
display.brewer.all()
brewer.pal(n = 4, name = "Set1")



#####

#recreate with CALR colored by group

CALR_grouped <- CALR_ordered %>% 
  t() %>% 
  as.data.frame() %>% 
  mutate_at(.vars = vars(CBL,NRAS,KRAS,PTPN11,GNB1,GNAS,FLT3,RIT1,BRAF,SH2B3), 
            .funs = funs(str_replace(., "MUT", "SIGNAL"))) %>% 
  mutate_at(.vars = vars(TET2,DNMT3A,IDH1,IDH2,ASXL1,EZH2,BCOR,BCORL1), 
            .funs = funs(str_replace(., "MUT", "EPI"))) %>% 
  mutate_at(.vars = vars(SRSF2,SF3B1,U2AF1,ZRSR2),
            .funs = funs(str_replace(., "MUT", "SPLICE"))) %>% 
  mutate_at(.vars = vars(SMC1A,STAG2,RAD21),
             .funs = funs(str_replace(., "MUT", "COHESIN"))) %>% 
  mutate_at(.vars = vars(RUNX1,ETV6,GATA2),
            .funs = funs(str_replace(., "MUT", "TF"))) %>% 
  t() %>% 
  as.data.frame() %>% 
  set_colnames(colnames(CALR_ordered))

save(CALR_grouped, file = "/Users/Will/R_Files/Projects/Sequencing/R_Data/CALR_grouped.rda")
save(color_diff, file = "/Users/Will/R_Files/Projects/Sequencing/R_Data/color_diff.rda")
index2 <- c(1,13:24,12,2:11, 32, 28:30, 25:27, 31)
Group_supp <- list(
  alter_fun = list(
    Type_I = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "#E41A1C", col = NA))
    },
    Type_II = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "#377EB8", col = NA))
    },
    SIGNAL = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "#4DAF4A", col = NA))
    },
    EPI = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "#984EA3", col = NA))
    },
    SPLICE = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "#FF7F00", col = NA))
    },
    COHESIN = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "#FFFF33", col = NA))
    },
    TF = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "#A65628", col = NA))
    },
    MUT = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h*0.33, gp = gpar(fill = "#F781BF", col = NA))
    },
    background = function(x, y, w, h) {
      grid.rect(x, y, w-unit(0.5, "mm"), h-unit(0.5, "mm"), gp = gpar(fill = "#CCCCCC", col = NA))
    }
  ),
  col = c("Type_I" = "#E41A1C", "Type_II" = "#377EB8", 
          "SIGNAL" = "#4DAF4A", "EPI" = "#984EA3", 
          "SPLICE" = "#FF7F00", "COHESIN" = "#FFFF33", 
          "TF" = "#A65628", "MUT" = "#F781BF"),
  order = list(column_order = colnames(CALR_grouped), row_order = rownames(CALR_grouped)[index2]),
  row_reorder <- c("CALR", "TET2", "DNMT3A", "IDH1", "IDH2", 
                   "ASXL1", "EZH2", "BCOR", "BCORL1", "SRSF2", "SF3B1", 
                   "U2AF1", "ZRSR2", "SH2B3", "CBL", "NRAS", " KRAS", "PTPN11", 
                   "GNB1", "GNAS", "FLT3", "RIT1", "KIT", "BRAF", "TP53", 
                   "RUNX1", "ETV6", "GATA2", "SMC1", "STAG2", "RAD21", "WT1"),
  dims = dim(CALR_grouped),
  patient_n = length(unique(gsub("(.*)S.*","\\1", colnames(CALR_grouped)))))

#create pdf of oncoprint
pdf(file = "/Users/Will/Desktop/plot7.pdf", bg = "white")
tiff(file = "/Users/Will/R_Files/Projects/Sequencing/Figures/test.tiff", 
     bg = "white", res = 300, width = 2200, height = 2500)
ht <- oncoPrint(CALR_grouped,
                alter_fun = Group_supp$alter_fun, col = Group_supp$col,
                remove_empty_columns = FALSE,
                column_order = Group_supp$order$column_order,
                row_order = Group_supp$order$row_order,
                column_title = "cBioPortal PTs with mutant CALR Types I and II",
                heatmap_legend_param = 
                  list(title = "Alterations", at = c("Type_I", "Type_II", "EPI","SPLICE",
                                                     "SIGNAL", "COHESIN", "MUT", "TF"),
                                            labels = c("CALR Type I", "CALR Type II", 
                                                       "Epigenetic","Signaling",  "Splicing", "Cohesin", 
                                                       "Other", "Transcription Factor"),
                                            nrow = 1, title_position = "topcenter"))
draw(ht, heatmap_legend_side = "bottom")
dev.off()
?tiff
         