
# Fetch p-value tables from official MINE website
# Turn them into internal data

# This script recquires an internet connection to work
# ------------------------------------------------------------------------------------

library(RCurl);
library(stringr);
library(readr);
library(devtools);

# Get all the lines of the source code of the page containing the links to the pval-tables
# ------------------------------------------------------------------------------------
website <- "http://www.exploredata.net"
internet_directory <- file.path(website, "Downloads/P-Value-Tables");
source_code <- readLines(internet_directory);

# ------------------------------------------------------------------------------------
# Get each link to a pval-table
table_match <- str_extract(source_code, '/ftp/Pvalues/.*[.]csv');
table_relative_link <- table_match[-which(is.na(table_match))];
downloadable_link <- paste(website, table_relative_link, sep = "");
n_table <- length(downloadable_link);

# ------------------------------------------------------------------------------------
# The Key of each table is the sample size
sample_size <- str_extract(downloadable_link, '[0-9]{2,4}');

# ------------------------------------------------------------------------------------
# Read and store each pval-table
pval_tables <- list();
for(i in 1:n_table) {
  curr_file <- read.csv(downloadable_link[i], sep = ",");
  end <- length(curr_file$X) - 1;
  table <- data.frame(MIC = curr_file$X[11:end],
                      pval = curr_file$X.1[11:end],
                      conf = curr_file$X.2[11:end]);
  pval_tables[[i]] <- apply(table, 2, as.numeric)
}
names(pval_tables) <- sample_size;

# ------------------------------------------------------------------------------------
# Save object in mage/R/sysdata.rda
use_data(sample_size, pval_tables, internal = TRUE, overwrite = TRUE);
