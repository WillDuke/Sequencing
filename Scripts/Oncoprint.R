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

#create a list with column order and endering details for oncoprint,
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

#####################################

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
  select(-CALR) %>% rename(CALR = CALR_type) %>% select(Sample.ID, CALR, 3:33) %>% t()
colnames(CALR_detail) <- as.character(CALR_detail[1,])
CALR_detail <- CALR_detail[-1,]
CALR_detail <- as.matrix(CALR_detail)
dim(CALR_detail)

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
  set_colnames(colnames(CALR_detail)) %>% rownames_to_column("IDs") %>% slice(c(-1, -33)) %>%
  column_to_rownames("IDs")


#found type classifications from COSMIC
write.csv(levels(as.factor(CALR_detail$CALR)),
          file = "/Users/Will/R_Files/Projects/Sequencing/CALR_mutations.csv")


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
  dims = dim(CALR_ordered),
  patient_n = length(unique(gsub("(.*)S.*","\\1", colnames(CALR_ordered)))))

save(CALR_detail, file = "/Users/Will/R_Files/Projects/Sequencing/CALR_detail.rda")
save(CALR_ordered, file = "/Users/Will/R_Files/Projects/Sequencing/CALR_ordered.rda")


#create pdf of oncoprint
pdf(file = "/Users/Will/Desktop/plot6.pdf", bg = "transparent")
ht <- oncoPrint(CALR_ordered,
                alter_fun = CALR_supp$alter_fun, col = CALR_supp$col,
                remove_empty_columns = TRUE,
                column_order = CALR_supp$order$column_order,
                row_order = CALR_supp$order$row_order,
                column_title = "cBioPortal PTs with mutant CALR Types I and II",
                heatmap_legend_param = list(title = "Alterations", at = c("Type_I", "Type_II", "MUT"),
                            labels = c("CALR Type I", "CALR Type II", "Concomitant Mutation"),
                            nrow = 1, title_position = "leftcenter"))
draw(ht, heatmap_legend_side = "bottom")
dev.off()

#without explicit ordering
ht_test <- oncoPrint(CALR_ordered, alter_fun = alter_fun, col = CALR_supp$col)
draw(ht_test, heatmap_legend_side = "bottom")



