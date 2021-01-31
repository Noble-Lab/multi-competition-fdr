library(multicomp)

parameter_file = commandArgs(trailingOnly = T)
arg_table = read.table(parameter_file, sep = "=")

fasta = arg_table[arg_table$V1 == "fasta_file",]$V2
ms2 = arg_table[arg_table$V1 == "spectra_file",]$V2
scoring = arg_table[arg_table$V1 == "scoring",]$V2
n_decoys = strtoi(arg_table[arg_table$V1 == "n_decoys",]$V2)
fdr_threshold = arg_table[arg_table$V1 == "fdr_threshold",]$V2
output_file = arg_table[arg_table$V1 == "output_file",]$V2

CRUX="`echo ./crux-*/bin/crux | head -n1 | awk '{print $1}'`"
val = 0

if("name" %in% arg_table$V1)
{
   dataSet = arg_table[arg_table$V1 == "name",]$V2
}else{
   dataSet = "dataSet"
}

if("mz_bin_width" %in% arg_table$V1)
{
   mz_bin_width = arg_table[arg_table$V1 == "mz_bin_width",]$V2
}else{
   mz_bin_width = 0.02
}

if("precursor" %in% arg_table$V1)
{
   precursor = arg_table[arg_table$V1 == "precursor",]$V2
}else{
   precursor = 50
}

if("keep_psms" %in% arg_table$V1)
{
  keep_psms = arg_table[arg_table$V1 == "keep_psms",]$V2
}else{
   keep_psms = "T"
}

if("seed" %in% arg_table$V1)
{
  seed = arg_table[arg_table$V1 == "seed",]$V2
}else{
  seed = "NULL"
}
# Preprocessing

system(paste("mkdir -p ", dataSet, sep = ""))

for(decoy in 0:(n_decoys - 1))
{
  index = paste(dataSet, "/dcy", decoy, ".index", sep = "")
  if(!file.exists(index))
  {
    system(paste(CRUX, " tide-index --allow-dups T --seed ", decoy," --overwrite T --output-dir ", index, " ", fasta, " ", index, sep = ""))
  }
  
  if(scoring %in% c("xcorr", "pvalue"))
  {
    root = paste(dataSet, ".", scoring, sep = "")
    if(scoring == "xcorr")
    {
      params = paste("--auto-mz-bin-width warn --mz-bin-width ", mz_bin_width, sep = "") 
    }
    else
    {
      params="--exact-p-value T --mz-bin-width 1.0005079"
    }
    tideResults = paste(root, "/dcy", decoy, ".tide-search.decoy.txt", sep = "")
    if(!file.exists(tideResults))
    {
      system(paste(CRUX, " tide-search --output-dir ", root, " --fileroot dcy", decoy, " --overwrite T --concat F --top-match 1 --precursor-window " , precursor, " --precursor-window-type ppm --auto-precursor-window warn ", params, " ", ms2, " ", index, sep = ""))
    }
    
    if(decoy == 0)
    {
      tideResults = paste(root, "/dcy", decoy, ".tide-search.target.txt", sep = "")
      for(method in c("tdc", "mix-max"))
      {
        confidenceResults = paste(root, "/", method, ".dcy", decoy, ".assign-confidence.target.txt", sep = "")
        if(!file.exists(confidenceResults))
        {
          system(paste(CRUX, " extract-rows --comparison eq --column-type int ", tideResults, " \"xcorr rank\" 1 > target.tmp", sep = ""))
          system(paste(CRUX, " extract-rows --comparison eq --column-type int \`echo ", tideResults, " | sed \'s/target/decoy/\'\` \"xcorr rank\" 1 > decoy.tmp", sep = ""))
          system(paste(CRUX, " assign-confidence --estimation-method ", method, " --output-dir ", root, " --fileroot ", method, ".dcy", decoy, " target.tmp", sep = ""))
          system("rm target.tmp decoy.tmp")
          
        }
      }
    }
  }else{
	  print("ERROR: Invalid scoring method, please use \"xcorr\" or \"pvalue\"")
	}
}
	
# Remove duplicated targets
if(n_decoys > 1)
{
	for(decoy in 1:(n_decoys - 1))
	{
  		system(paste("rm -f */dcy", decoy, ".tide-search.target.txt", sep = ""))
	}

	for(decoy in 0:(n_decoys - 1))
	{
  		system(paste("rm -f */dcy", decoy, ".tide-search.log.txt", sep = ""))
 		 system(paste("rm -f */dcy", decoy, ".tide-search.params.txt", sep = ""))
	}
}

if(scoring == "xcorr")
{
  is_xcorr = T
} else if(scoring == "pvalue")
{
  is_xcorr = F
} else
{
  print("Invalid scoring type. Defaulting to xcorr score")
  is_xcorr = T
}


get_peptide_scores = function(dat_list)
{
  tar = dat_list[[1]]
  total_peps = tar$sequence
  for(i in 2:length(dat_list))
  {
    total_peps = c(total_peps, dat_list[[i]]$original.target.sequence)
  }
  uniq_total_peps = unique(total_peps)
  
  tar_scores_sequence = uniq_total_peps
  tar_scores_scores = c()
  dec_score_matrix = matrix(-100, nrow = length(uniq_total_peps), ncol = length(dat_list) - 1)
  
  for(i in 1:length(uniq_total_peps))
  {
    if(uniq_total_peps[i] %in% tar$sequence)
    {
      tar_scores_scores[i] = max(tar[which(tar$sequence == uniq_total_peps[i]),]$xcorr.score) 
    }
    else{
      tar_scores_scores[i] = -100
    }
    
    for(j in 2:length(dat_list))
    {
      if(uniq_total_peps[i] %in% dat_list[[j]]$original.target.sequence)
      {
        dec_score_matrix[i, j-1] = max(dat_list[[j]][which(dat_list[[j]]$original.target.sequence == uniq_total_peps[i]),]$xcorr.score) 
      }
      else
      {
        dec_score_matrix[i, j-1] = -100
      }
    }
  }
  
  tar_scores = data.frame(sequence = tar_scores_sequence, xcorr.score = tar_scores_scores)
  return(list(tar_scores, dec_score_matrix))
}


read_all = function(min, max, data_name, use_xcorr)
{
  if(use_xcorr)
  {
    df.list = list()
    df.list[[1]] = read.csv(paste("./", data_name, ".xcorr/dcy0.tide-search.target.txt", sep = ""), sep = "\t")
    for(i in min:max)
    {
      df.list[[i+2]] = read.csv(paste("./", data_name, ".xcorr/dcy", toString(i),".tide-search.decoy.txt", sep = ""), sep = "\t")
    }
    return(df.list) 
  }
  else{
    df.list = list()
    df.list[[1]] = read.csv(paste("./", data_name, ".pvalue/dcy0.tide-search.target.txt", sep = ""), sep = "\t")
    for(i in min:max)
    {
      df.list[[i+2]] = read.csv(paste("./", data_name, ".pvalue/dcy", toString(i),".tide-search.decoy.txt", sep = ""), sep = "\t")
    }
    return(df.list) 
  }
}

dat_list = read_all(0, n_decoys - 1 , dataSet, is_xcorr) 

if(!is_xcorr)
{
    for(i in 1:length(dat_list))
    {
		old_name = names(dat_list[[i]])
		old_name[which(old_name == "exact.p.value")] = "xcorr.score"
		dat_list[[i]] = setNames(dat_list[[i]], old_name)
		dat_list[[i]]$xcorr.score = -1* dat_list[[i]]$xcorr.score
    }
}

if(seed != "NULL")
{
  set.seed(as.integer(seed))
}

pep_scores = get_peptide_scores(dat_list)

tar_scores = matrix(pep_scores[[1]]$xcorr.score, nrow = 1)
dec_scores = t(pep_scores[[2]])

alpha_range = as.double(fdr_threshold)
res_index = multidecoy_comp(tar_scores, dec_scores, alpha_range)[[1]]
res_peps = pep_scores[[1]][res_index,]$sequence

write(res_peps, output_file)

if(keep_psms == "F")
{
   system(paste("rm -r ", dataSet, ".xcorr", sep = ""))
   system(paste("rm -r ", dataSet, ".pvalue", sep = ""))
   system(paste("rm -r ", dataSet, sep = ""))
}

