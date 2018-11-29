roiNames <- read.table('../roi_names.txt')
names(roiNames) <- "ROI"

roiNames$voxelNumber <- rep(0, nrow(roiNames))
roiNames$voxelSizeMM <- rep(0, nrow(roiNames))

for (i in 1:nrow(roiNames)){
#       print(roiNames$ROI[i])
        message <- paste0("fslstats ../Harvard-Oxford_ROIs/", roiNames$ROI[i], " -V > temp.txt")
#       print(message)
        system(message)
        voxNum <- read.table("temp.txt")
#       print(voxNum)
        roiNames$voxelNumber[i] <- voxNum[1,1]
        roiNames$voxelSizeMM[i] <- voxNum[1,2]
}

write.csv(roiNames, file="roi_voxel_size.csv", col.names=TRUE)




