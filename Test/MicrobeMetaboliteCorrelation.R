setwd("~/Documents/R_WorkSpace/m2interact")
source('R/Heatmap.R')

sample_data <- load.meta.data('Data/IHMP/hmp2_metadata.csv', tax_column = 2)
# sample_data <- microbe_sample_data[, c(75, 34, 40)]
microbe_abundance_table <- load.abundance.data('Data/iHMP/taxonomic_profiles.tsv_AbundanceTable_2019-07-09.csv')
microbe_abundance_table <- microbe_abundance_table[, -c(1)]

# read ihmp metabolite data
# iHMP_metabolomics_HILIC-neg_060517
ihmp <- read.csv(
  
  'Data/iHMP/iHMP_metabolomics_HILIC-neg_060517.csv',
  
  header = TRUE,
  
  fill = TRUE,
  
  comment.char = "" ,
  
  check.names = TRUE,
  
  stringsAsFactors = FALSE
  
)

ihmp_col_names <- ihmp[4,-c(1:7)]

# filter out for only hmdb id
ihmp  <- ihmp[(ihmp$X.5) != "",]

# get column_names
ihmp_row_names <- ihmp[-1,6]

ihmp <- ihmp[-1,-c(1:7)]

ihmp <- apply(as.matrix(ihmp), 2, as.numeric)

# clean up ihmp data
colnames(ihmp) <- ihmp_col_names
rownames(ihmp) <- ihmp_row_names

ihmp[is.na(ihmp)] <- 0

# Match up samples

for(i in 1:ncol(microbe_abundance_table)) {
  a <- colnames(microbe_abundance_table)
  for(j in 1:nrow(sample_data)) {
    if(colnames(microbe_abundance_table)[i] == rownames(sample_data)[j])
      a[i] <- sample_data[j,3]
  }
  colnames(microbe_abundance_table) <- a
}

a <- microbe_abundance_table[, colnames(microbe_abundance_table) %in% colnames(ihmp)]
q <- ihmp[, colnames(ihmp) %in% colnames(microbe_abundance_table)]

a <- a[, !duplicated(colnames(a))]
q <- q[, !duplicated(colnames(q))]

a <- a[, colnames(a) %in% colnames(q)]
q <- q[, colnames(q) %in% colnames(a)]

z <- cor(x = t(a), y = t(q), method = 'spearman')

z[is.na(z)] <- 0

s_meta <- metabolite_data[rownames(metabolite_data) %in% colnames(z), ]
f_meta <- microbe_data[rownames(microbe_data) %in% row.names(z), ]

create.heatmap(data=z, sample_meta = metabolite_data[,c(91,92)], feature_meta = microbe_data[,c(7, 25, 24)], show = TRUE, omit_na = FALSE)
# Get single microbe
z.r <- cor(t(a[1,]), t(q))
z.r.x <- list()

a.b <- a[,colnames(q)]

for(i in 1:nrow(a.b)) {
  a.l <- a.b[i, , drop = FALSE]
  a.l <- cor(t(a.l), t(q), method = 'spearman')
  z.r.x <- rbindlist(list(z.r.x, as.data.frame(a.l)), use.names = TRUE)
}

z.r.x[is.na(z.r.x)] <- 0
z.r.x <- as.data.table(lapply(z.r.x, as.numeric))
z.r.x <- as.matrix(as.data.frame.data.frame(z.r.x))
row.names(z.r.x) <- rownames(a.b)
create.heatmap(data = z.r.x,  sample_meta = metabolite_data[,c(91,92)], feature_meta = microbe_data[,c(7, 25, 24)], show = TRUE, omit_na = FALSE)


a.l <- cor(t(a.b), t(q), method = 'spearman')
create.heatmap(data = a.l,  sample_meta = metabolite_data[,c(91,92)], feature_meta = microbe_data[,c(7, 25, 24)], show = TRUE, omit_na = FALSE)
# Determine if correlation with metabolites

# Look at only metabolites that microbe known to use

# Look at only metabolites that microbe known to produce