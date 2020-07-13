#==============Remove the rows of a data table=============

#-------remove if 'na'--------

df<- df[-which(is.na(df$start_pc)), ]


#-------remove if blank-------

 df[!(is.na(df$start_pc) | df$start_pc==""), ]
