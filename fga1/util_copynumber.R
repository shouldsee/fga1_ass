
find_seg <- function(probe,segEND){
  which((probe <= segEND[-1]) & (probe > segEND[-length(segEND)]))
}

assign_copynumber <- function(probidx, p_vct, L_vct, L_chromosome = tail(probidx,1), maxIter = 5, guess = 1.5){
  
  L_avg = p_vct %*% L_vct 
  seg_count = (L_chromosome / L_avg) 
  seg_states = c()
  segEND = c(0)
  # segEND_lst = vector(maxIter,'list')
  for (i in 1:maxIter){
    # sample_segment( N = seg_count * 1.5 %f% floor, p_vct, L_vct, s_vct = c(1,2,3) )
    N = (seg_count * guess) %f% floor
    s_vct =  c(1,2,3)
    seg_states_new = sample( s_vct, N ,replace = T,prob = p_vct)
    seg_states = c(seg_states, seg_states_new)
    segL_vct =  rpois( 1 : (seg_states%f%length), L_vct[seg_states])
    # segL_vct = sample(c(1,2,3), N ,replace = T,prob = p_vct)
    segEND_new = tail(segEND,1) + cumsum(segL_vct)
    segEND = c(segEND,segEND_new)
    allEND = tail(segEND,1) 
    if  ( L_chromosome  - allEND < 0 ){
      ### Length achieved 
      break
    }else{
      if (i==maxIter){
        #### This should be rare
        stop("[ERROR]: Failed to generate a chromosom of wanted length, try increase 'maxIter' ")
      }else{
        ###  update estimated chromosome
        seg_count = (L_chromosome - allEND) / L_avg
      }
    }
  }
  # segEND[-1]
  # segEND
  
  # probe_states=sapply(probes, FUN = function(x){find_seg(x,segEND)})
  probe_states=sapply(probeidx, FUN = function(x){find_seg(x,segEND)})
  
  
  # %f%as.vector
  cpnum = seg_states[probe_states]
}

log_convert <- function(vct){
  log2(vct/2)
}

noisy <- function(vct, sigma = .5, relative = F){
  if (relative){
    sigma =  abs(vct) * sigma
  }
  vct + rnorm(vct, mean = 0, sd =sigma)
}


######  Transform a per-segment representation into a per-probe representation.
flatten <- function(idx, ENDs,states){
  if (length(ENDs) == length(states)){
    ENDs = c(0,ENDs)
  }
  probe_states=lapply( idx, FUN = function(x){find_seg(x,ENDs)}) %f%unlist
}


######  Extract fitted segments and flatten into per-probe representation
extract_fitted <- function(CNAoutput, idxs = 1:max(CNAoutput$segRows$endRow)){
  states = CNAoutput$output$seg.mean
  fitted = states[
    flatten(
      idx = idxs,
      ENDs = CNAoutput$segRows$endRow,
      states = states
    )
    
    ]
  
}


crossref<- function( fit,model, logical=T){
  TP = model  & fit
  FP = !model & fit
  TN = !model & !fit
  FN = model  & !fit
  if (!logical){
    TP = sum(TP)
    FP = sum(FP)
    TN = sum(TN)
    FN = sum(FN)
  }
  return(list(TP = TP,
              FP = FP,
              TN = TN,
              FN = FN))
}

thres2PR <- function(thres, fitted,model){
  fitted_binary = fitted < thres
  model_binary = model == -1
  stats = crossref(fitted_binary, model_binary, logical = F)
  # print((stats$TP + stats$FP) == 0)
  if ((stats$TP + stats$FP) == 0){
    return(NULL)
  }else{
    P = stats$TP /(stats$TP + stats$FP)
    R = stats$TP /(stats$TP + stats$FN)
    return(list(P=P,R=R))
  }
  
}


calc_mad <- function(vct){
  median(abs(vct - median(vct) ))
}
seg_thres <- function(signal, MIN_LEN=5, Z_MIN=4, seg = 1, SUPER_THRES  = -0.5 , DEBUG = F){  
  if (length(signal) < 20){
    sprintf("[WARNING] Not enough observation (N = %d < 20) to perform clustering reiliably", length(signal))%f%warning
  }
  vct = sort(signal)
  # start = 0
  # vcti = vct[-(1:start)]
  vcti = vct
  x = c()
  
  if (DEBUG){
    xs = c()
    ys = c()
    zs = c()
  }
  outdf = data.frame(
    avg = NULL,
    stdev=NULL,
    N = NULL,
    t.pval=NULL,
    include = NULL
  )
  
  # superthres = -0.5
  thres_curr = NA
  for (ii in  1:length(vcti)){
    i  = vcti[ii]
    if (length(x) >= MIN_LEN){
      
      mad = calc_mad(x)
      z = (i-mean(x)) / mad
      p = pnorm(z)
      
      if (!is.na(z) & seg){
        if (z > Z_MIN | ii == length(vcti)){
          # ?t.test
          # pval = 
          mx = mean(x)
          ##### Simple thresholding the average 
          # include = (mx < superthres)
          
          ##### A modified one-tail t-test using MAD as an estimator for sample standard deviation.
          ##### This makes p-value more robust w.r.t outlier.
          pval = pt( (mx - SUPER_THRES)/ (mad/sqrt(length(x))),df = length(x) -1, lower = T)
          include = pval < 1E-3
          
          
          ##### By default, R uses sd() to estimate sample stdev.
          # stat = t.test(x, mu = superthres, alternative = 'smaller')
          # pval = stat$p.value
          # paste0(round(mean(x),5),' ',round(mad,5),' ',length(x))%f%print
          # print(pval) 
          # include = 
          outdf = rbind(outdf,list(avg = mean(x) 
                                   ,stdev = mad
                                   ,N = length(x)
                                   ,include = include))
          if (include){
            thres_curr = max(thres_curr, x, na.rm=T)
          }
          #### Reset the sample vector "x"
          x = c()
        }
      }
    }else{
      z = 0
    }
    x  = c(x,i)
    if (DEBUG){
      ys = c(ys,z)
      xs = c(xs,i)
      zs = c(zs,mad)
    }
    # zs = c(zs,sd(x))
    # print(p)
  }
  # if (sum(outdf$include)==0){
  #   thres_curr = NA
  # }
  if (!DEBUG){
    return(thres_curr)
  }else{
    print(outdf)
    print(thres_curr)
    try({
      plot(xs, ylab = 'logR',xlab = 'Segment Index', ylim = c(-1.5,1.5)
           , main = tail( df$data %f% names,1))
      abline(h = -0.5, col= 2,lty = 2)
      abline(h = thres_curr, col= 3,lty = 2)
      par(new = T)
      # plot(ys,axes = F)
      grid()
      # par(new = T)
      # # zs = exp(zs)
      # plot(zs, ylab = '',axes = F,type = 'l')
    })
    return(list(thres_curr,outdf))
    
  }
  
}


