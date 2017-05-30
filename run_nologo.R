weightYPD<-function(x){
  weightY<-c(A=0.05493370, C=0.01259810, D=0.05851865, E=0.06538412, F=0.04417271, G=0.04974184, H=0.02167693, I=0.06556848, K=0.07348676, L=0.09500875, M=0.02080609, N=0.06162645, P=0.04377389, Q=0.03954864, R=0.04435365, S=0.08984322, T=0.05917161, V=0.05559829, W=0.01039398, Y=0.03379412) 
  weightY[x]
}


read_params2 <- function (file) {
  d<-read.delim(file,header=FALSE);
  f<-as.vector(d$V2);
  names(f)<-d$V1;
  return (f)
  
}

read_spacerparams<-function(file){
  d<-read.table(file, header=TRUE, sep="\t", row.names=1, as.is=F)
  f<-as.matrix(d)
  return(f)
}

pr_given_nologo <- function (x, fhyd, fss, fspa) { #depends on the weightYPD function defined above.
  L<-0; 
  for (s1 in 1:3) {
    prs1<-prod(as.numeric(fspa[x[(1+1):(1+s1)]]));
    ts1<-fhyd[[x[1]]] * fss[s1,1]* fhyd[[x[1+s1+1]]] * prs1
    
    for (s2 in 1:3) {
      prs2<-prod(as.numeric(fspa[x[(1+s1+1+1):(1+s1+1+s2)]]))        
      ts2<-ts1 * fss[s2,2] * fhyd[[x[1+s1+1+s2+1]]] * prs2
      
      for (s3 in 1:3) {
        prs3<-prod(as.numeric(fspa[x[(1+s1+1+s2+1+1):(1+s1+1+s2+1+s3)]]))
        prstail<-1
        if ((1+s1+1+s2+1+s3+1)<length(x)) {
          prstail<-prod(weightYPD(x[(1+s1+1+s2+1+s3+1+1):length(x)]))
        }
        term <- ts2 * fss[s3,3]* fhyd[[x[1+s1+1+s2+1+s3+1]]] * prs3 * prstail 
        
        L<-L+term
      }
    }
  }
  return (L);	
}


fss<-read_spacerparams('spacer_length35.txt');
fspa<-read_params('spacer_model35.txt');
fhyd<-read_params('hydrophobic_model35.txt');


cutoff<-0;


matchnum<-0; matchdata<-NULL; 

for (n in 1:length(sequencelist))
{ 
  for (i in 1:(length(sequencelist[[n]])-12))
  {
    pepbeg <- i;
    pepend <- i+12;
    sequencelist[[n]][pepbeg:pepend]->peptidewindow;
    Lbg<-prod(weightYPD(peptidewindow));
    Lnologo<-pr_given_nologo(peptidewindow,fhyd,fss,fspa);
    pepscore<-log(Lnologo/Lbg);
    #mllist2Alan[[n]][[i]]<-pepscore;
    
    if (pepscore>=cutoff) { #store the peptide windows that pass the cutoff.
      seq<-paste(peptidewindow,collapse="");
      matchdata[[matchnum+1]]<-c(names(sequencelist[n]),pepbeg,pepend,pepscore,seq);
      matchnum<-matchnum+1;
    }
  }
}



do.call("rbind", matchdata)->matchdata_df
write.table(matchdata_df, "Nologo35NESscutoff0.txt")

