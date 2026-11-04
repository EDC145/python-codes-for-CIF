






files = list.files(path = "D:/research v3.0/project associate/Complex Interface packing/7. Bulk Packing attempt - 1/Only_Atomic_Coords/", 
                   pattern = ".txt")

df = data.frame()

for (i in files){
  print (paste("Reading ",i))
  a <- read.csv(paste0("D:/research v3.0/project associate/Complex Interface packing/7. Bulk Packing attempt - 1/Only_Atomic_Coords/",i), sep = " ", header = FALSE)
  xmax = max(a[,2])
  xmin = min(a[,2])
  ymax = max(a[,3])
  ymin = min(a[,3])
  zmax = max(a[,4])
  zmin = min(a[,4])
  t <- data.frame(xmin, xmax, ymin, ymax,zmin, zmax, i,stringsAsFactors = FALSE)
  df = rbind(df,t)
}


write.table(df, "Input_for_voro.txt", row.names = FALSE, quote = FALSE)









