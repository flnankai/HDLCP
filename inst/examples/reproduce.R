# Run manually after installing HDLCP and DMS. This is not executed on installation.
library(HDLCP)
results <- reproduce_paper(output_dir="hdlcp-results",workers=25L)
# Application (download UCI files first):
# reproduce_secom("secom.data","secom_labels.data",workers=25L)
