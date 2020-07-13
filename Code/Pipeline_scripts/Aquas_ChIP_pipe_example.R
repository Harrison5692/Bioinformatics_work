source("/storage/cylin/grail/projects/Lagor/chipseq_aquas_pipeline_fixed.R")

##############
project.dir <- "/storage/cylin/grail/projects/Lagor/aquas_pipe_hg38_CHIP_good"
peaks.dirname <- "macsEnriched"


########
dt.name <- "hg38_liver_CHIP_data_table.txt"
datainfo.location <- file.path(project.dir, "data_tables", dt.name)

##Step1:Alignment and pre peak call QC
#pushAquasAlignQC(project.dir, datainfo.location)

##Step2:Call peaks
#collectAquasQC(project.dir)
#mvBams(project.dir, datainfo.location)
#callAquasPeaks(project.dir, datainfo.location)

##Step3: Get QC afterthe peak call
mvPeakFiles(project.dir, datainfo.location)
#createSignalTrack(project.dir, datainfo.location)
df <- getCustomQC(project.dir, datainfo.location, peaks.dirname)
output.filename <- sprintf("post_peak_call_qc_%s", dt.name)
output.location <- file.path(project.dir, "qc_summary", output.filename)
write.table(df, file = output.location, col.names = T, 
	row.names = F, quote = F, sep = "\t")
