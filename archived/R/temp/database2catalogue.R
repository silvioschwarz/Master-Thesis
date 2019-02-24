# Set working directory

path = getwd()
setwd(path)

rm(list = ls())

#needed for distance calculations
library(fields)


DMBI04_2_Catalogue = function(filename_in, eqlist, filename_out){
                
                data = read.csv(filename_in, dec=",", as.is = TRUE)
                
                eq = read.csv(eqlist, sep = "/t", header = TRUE)
                
                data$IsOr = as.numeric(data$IsOr)
                data = data[!is.na(data$IsOr),]
                data = data[!is.na(data$Gi),]
                
                Date = character(length = nrow(data))
                R = numeric(length = nrow(data))
                
                
                for(i in 1:nrow(data)){
                        Date[i] = format(as.Date(paste(data$Gi[i],data$Me[i], data$An[i]),format("%d %m %Y")), format="%d/%m/%Y")
                        
                        R[i] = diag(rdist.earth(matrix(c(data$LonOr[i], data$LatOr[i]), ncol = 2), 
                                                matrix(c(data$LonEp[i], data$LatEp[i]), ncol = 2), miles = FALSE, R = 6371))
                } 
                
                daten = data.frame(ID = 1, Datum = Date, Ix = data$Ix, I0 = data$Io, Mw = data$Mw, Mw_err = data$Daw,
                                   Is = data$IsOr/10,  R = R, LatEp = data$LatEp, LonEp = data$LonEp, LatOr = data$LatOr, LonOr = data$LonOr )
                
                test = daten[1,]
                id = 1
                
                for(i in 1:nrow(eq)){
                        
                        
                        test1 = daten[intersect(which(Date == as.character(eq$Date[i])),which(data$Io == eq$I0[i])),]
                        
                        # check on the number of observations
                        if(nrow(test1) < 10){
                                test = test
                        }
                        else {test1$ID = id;
                              id = id +1;
                              test = rbind(test,test1)}
                        
                        
                }
                test = test[-1,]
                
                x = strsplit(filename_in, split = "/")[[1]]
                x = x[1:length(x)-1]
                file_out = paste(paste(x,collapse = "/"), filename_out, sep = "/")
                
                write.table(test, file = file_out, row.names = FALSE)
                
                return(test)
                
}



######
eqlist= "./Masterarbeit/Masterarbeitdata/Rotondi 2004/D47.txt"
filename_in = "./Masterarbeit/Masterarbeitdata/Rotondi 2004/dbmi04.csv"
filename_out = "D47.dat"
catalogued47  = DMBI04_2_Catalogue(filename_in, eqlist, filename_out)

eqlist2= "./Masterarbeit/Masterarbeitdata/Rotondi 2004/DZ-47.txt"
filename_in2 = "./Masterarbeit/Masterarbeitdata/Rotondi 2004/dbmi04.csv"
filename_out2 = "DZ-47.dat"
catalogueDZ47  = DMBI04_2_Catalogue(filename_in2, eqlist2, filename_out2)
