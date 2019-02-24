pb = winProgressBar(title="Example progress bar", 
                     label="0% done", 
                     min=0, 
                     max=100, 
                     initial=0)
for(i in 1:100) {
  Sys.sleep(0.1) # slow down the code for illustration purposes
  info = sprintf("%d%% done", round((i/100)*100)) 
  setWinProgressBar(pb, i/(100)*100, label=info) 
} 
close(pb)
