################# Return table if at least of of a group as a certain value ####################
################################################################################################
####################### for filtering out lowly expressed genes ################################

# example: This will keep the old table but remove rows if the rows in either of the 3 columns has a value less than 10


new_table <- subset(old_table, old_table$column_name & old_table$column_name2 & old_table$column_name3 > 10)