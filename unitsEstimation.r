# Copyright 2017 Martha Dais Ferreira, Rodrigo Fernandes de Mello
# 
# This file is part of ImageFalseNearest.
# 
# ImageFalseNearest is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# ImageFalseNearest is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with ImageFalseNearest.  If not, see <http://www.gnu.org/licenses/>.
# 
# See the file "COPYING" for the text of the license.
# 
# Contact: 
# 	Martha D. Ferreira: daismf@icmc.usp.br
# 	Rodrigo Fernandes de Mello: mello@icmc.usp.br

require(png)
require(jpeg)
require(pixmap)
require(ImageFalseNearest)

#
# Choosing best number of units considering all color channels
#  Input:
#     folder_initial: folder with the FNN results
#            maskRow: row of the mask size
#            maskCol: column of the mask size
#
#  Output:
#	   ret$maskRow: row of the mask size
#	   ret$maskCol: column of the mask size
# 	     ret$posMR: k best number of units
#
analise_NK_all <- function(folder_initial,maskRow,maskCol){
	
	#Reading FNN result
	cat(folder_initial,"\n")
	normMR = as.matrix(read.table(paste(folder_initial,"meanResultsAllColors.dat",sep="")))

	#Sorting values and taking the position of the first zero
	posMR = sort.list(normMR,dec=F)[1]
	
	#output
	ret = list()
	ret$maskRow = maskRow
	ret$maskCol = maskCol
	ret$posMR = posMR
	
	ret
}

#
# FNN to Units Evaluation
# Input:
#      dataset: the name of the dataset
#     listfile: a file containing paths to dataset images (PNG)
#       maxRow: list of possible numbers of row mask size
#       maxCol: list of possible numbers of column mask size
#            k: number of the best mask sizes to be returned
#    threshold: threshold to stop the images analysis
#     nKernels: maximum number to evaluate the units
#         endT: number of neighbors to be considered in FNN
#
# Output:
#    ret
#            ret$datset: dataset name
#    ret$folder_initial: folder create to output files
#         ret$bestUnits: best number of units found
#
NK_evaluation <- function(dataset, listfile, maskRow, maskCol, k=5, threshold=1e-4, nKernels=50, endT=1){
	
	cat("-------Compute FNN for number of units estimation.-------\n")
	
	#Creating folders to output
	folder_initial = paste(dataset,"/fnnNKernelEstimate/",sep="")
	system(paste("mkdir ",dataset,sep=""))
	system(paste("mkdir ",folder_initial,sep=""))

	#Reading files
	Files = read.table(listfile)
	sizeFiles = length(Files[[1]])
	cat("Analysing ",sizeFiles," files\n")
	
	#Taking images dimension
	channels=1
	# img = readPNG(as.character(Files[[1]][1])) #png
	# img = as.matrix(readJPEG(as.character(Files[[1]][1]))) #jpg
	img = read.pnm(file = as.character(Files[[1]][1]),cellres=1)@grey #pgm
	if(length(dim(img))>2){
		channels=3
		img = img[,,1]
	}

	#Preparing variables
	result_color = list()
	result_all = list()
	
	#randomize examples
	ids = sample(length(Files[[1]]))
	files = Files[[1]][ids]
	
	#Computing for possible numbers of mask sizes
	for(currMask in 1:length(maskRow)){
		
		mean_normMR = rep(0, nKernels)
		mean_normFH = rep(0, nKernels)
	
		#Computing for three color channels
		for(color in 1:channels){
			#1-red 2-green 3-blue
			if(color == 1) color_name = "RED"
			if(color == 2) color_name = "GREEN"
			if(color == 3) color_name = "BLUE"
			if(channels == 1) color_name = "GRAY"
			cat("\nAnalysing color ",color_name," for mask ",maskRow[currMask],"x",maskCol[currMask],"\n")

			#folder for mask
			folder_mask = paste(folder_initial,"Mask-",maskRow[currMask],"-",maskCol[currMask],"/",sep="")
			system(paste("mkdir ",folder_mask,sep=""))
			
			#folder for color
			folder = paste(folder_mask,color_name,"/",sep="")
			system(paste("mkdir ",folder,sep=""))
			
			#Preparing variables
			sumResultskernels = rep(0, nKernels)
			histogramAll = rep(0, nKernels)
			histogramAll_before = histogramAll

			eucNorm = c()
			curr_norm = 1
			i=1
			while (((curr_norm > threshold) || (curr_norm == 0)) && (i < sizeFiles)){
				if(channels==3){
					# img = readPNG(as.character(files[i]))[,,color] #png
					# img = as.matrix(readJPEG(as.character(files[i]))) #jpg
					img = read.pnm(file = as.character(files[i]),cellres=1)@grey #pgm
				} else {
					# img = readPNG(as.character(files[i])) #png
					# img = as.matrix(readJPEG(as.character(files[i]))) #jpg
					img = read.pnm(file = as.character(files[i]),cellres=1)@grey #pgm
				}
				
				#Taking file name for output
				fileFirstname = unlist(strsplit(as.character(files[i]), "/",fixed = TRUE))
				nameId = length(fileFirstname)
				fileFirstname = unlist(strsplit(fileFirstname[nameId], ".",fixed = TRUE))[1]
				cat("\tfile[",i,"]=",fileFirstname,"\n")
				
				#Appling FNN for units evaluation in a single image considering the current color
				cat("\t\tCalculating number of kernels...\n")
				fractionMap = rep(0, nKernels)
				resultsFNN = assess.number.of.kernels(img,mask.nrows=maskRow[currMask], mask.ncols=maskCol[currMask], end.T=endT, nkernels = 1:nKernels, plot.opt=F)$heatmap
				for (k in 1:nrow(resultsFNN)) {
						fractionMap[k] = mean(resultsFNN[k,])
				}
				
				#Saving results for each image considering the current color
				#nameFile = paste(folder,fileFirstname,".dat",sep="")
				#write.table(fractionMap,nameFile,row.names=F,col.names=F)
				
				#Summing results
				sumResultskernels = sumResultskernels + fractionMap
				
				#Saving results
				Resultskernels = sumResultskernels/i
				write.table(Resultskernels,paste(folder,"outputMeanResults.dat",sep=""),row.names=F,col.names=F)
				
				#Computing Euclidean 
				histNorm = Resultskernels/max(Resultskernels)
				A = histogramAll_before - histNorm 
				curr_norm = sum(A^2)
				if(is.nan(curr_norm)) curr_norm=0
				eucNorm = c(eucNorm, curr_norm)
				cat("\t\tEuclidean Norm: ",curr_norm,"\n")
				
				histogramAll_before = histNorm
				i=i+1
			}
			
			#Average FNN results considering the color under analysis
			normMR = Resultskernels/max(Resultskernels)
			#Summing results
			mean_normMR = mean_normMR + normMR
		}
		
		#Computing the mean of all		
		mean_normMR = mean_normMR/channels
		
		#Saving files
		cat("Saving Files...")
		nameFile = paste(folder_mask,"meanResultsAllColors.dat",sep="")
		write.table(mean_normMR,nameFile,row.names=F,col.names=F)
		
		#compute possible number of kernels values
		result_all[[currMask]] = analise_NK_all(folder_mask,maskRow[currMask],maskCol[currMask])
	}
	
	#Printing and Saving
	infoPrint = c()
	cat("\nBest Number of Units:\n")
	for(i in 1:length(result_all)){
		cat("\tMask:",result_all[[i]]$maskRow,"x",result_all[[i]]$maskCol,"\t")
		cat("Number of Units: ",result_all[[i]]$posMR,"\n")
		infoPrint = rbind(infoPrint,c(result_all[[i]]$maskRow,result_all[[i]]$maskCol,result_all[[i]]$posMR))
	}
	nameFile = paste(folder_initial,"bestUnits.dat",sep="")
	write.table(infoPrint,nameFile,row.names=F,col.names=F)
	
	#output
	ret = list()
	ret$dataset = dataset
	ret$folder_initial = folder_initial
	ret$bestUnits = result_all
	
	ret 
}

