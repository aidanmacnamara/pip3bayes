require(stringr)

# correct input data so that headers are:
# 'readout'_'experiment'_'run'

data = read.csv("../share/data/H2O2_data_all_wide.csv")
# remove first column
data = data[,-1]
headers = names(data)[-1]
readouts = gsub("^.*\\.([[:alnum:]_]+)$", "\\1", headers)
readouts = str_replace(readouts, "_[[:alnum:]]+$", "")
expers = gsub("^Exp([[:digit:]]+).*$", "\\1", headers)
runs = gsub("^.*ROI([[:digit:]]+)\\..*$", "\\1", headers)
names(data)[-1] = paste(readouts,expers,runs,sep="_")
write.csv(data, file="../share/data/H2O2_data_all_wide_corrected.csv")

data = read.csv("../share/data/SingleCell.csv")
data = data[,-1]
headers = names(data)[-1]
headers = str_replace(headers, 'fk506\\.', 'fk506_1\\.')
readouts = gsub("^.*\\.([[:alnum:]]+)\\..*$", "\\1", headers)
expers = gsub("^cd1_fk506_([[:digit:]]+).*$", "\\1", headers)
runs = gsub("^.*?([[:digit:]]+)$", "\\1", headers)
names(data)[-1] = paste(readouts,expers,runs,sep="_")
write.csv(data, file="../share/data/SingleCell_corrected.csv")
