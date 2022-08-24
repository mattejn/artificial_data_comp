library('FNN')
trans = function(x)
  asinh(x / 10)
norm = function(x)
  (x - min(x)) / (max(x) - min(x))
recur_unpack = function(pars, df, todel) {
  par_rows = df[df$pop %in% pars, ]
  new_pars = c()
  for (r in 1:nrow(par_rows)) {
    if (par_rows[r, 'pop'] %in% df[, 'parent_pop']) {
      par_row = par_rows[r, ]
      children_rows = df[df$parent_pop == par_row$pop, ]
      children_rows$relative_cnt = children_rows$relative_cnt * par_row$relative_cnt
      for (i in 1:nrow(children_rows)) {
        for (k in 1:ncol(children_rows)) {
          #there is definitely a way to do this better
          if (is.na(children_rows[i, k])) {
            children_rows[i, k] = par_row[1, k]
          }
        }
      }
      df[rownames(children_rows), ] = children_rows
      children = children_rows$pop
      for (child in children) {
        if (nrow(df[df$parent_pop %in% child, ]) != 0) {
          new_pars = c(new_pars, child)
        }
        todel = c(todel, as.integer(rownames(par_row)))
      }
    }
    if (r == nrow(par_rows) &&
        length(new_pars) == 0) {
      return(df[-unique(todel), ])
    }
  }
  recur_unpack(new_pars, df, todel)
}

gate <- function(data,pair, neighbors,cutoff_perc) {
  df=data[,pair]
  nn=get.knn(df, k=neighbors, algorithm=c("kd_tree", "cover_tree", "CR", "brute"))$nn.dist
  data=cbind(data,rowSums(nn)/ncol(nn))
  n=ncol(data)
  new=data[ data[,n] < quantile(data[,n] , cutoff_perc ) , ]
  new=new[,-n]
  return(new)
}

generate_artificial_cytometry_data = function(spectra, pheno, n, stronk, uncertainty=1,dead_fraction=0.175) {
  spectra = spectra / sqrt(rowSums(spectra ^ 2)) 
  d = ncol(spectra) 
  n_pheno = nrow(pheno)
  pheno[is.na(pheno)] = 0
  rownames(pheno) = paste0(pheno$parent_pop, pheno$pop)
  bkp = rownames(pheno)
  pheno$relative_cnt = pheno$relative_cnt * (1 / sum(pheno$relative_cnt))
  counts_pheno = pheno$relative_cnt
  pheno = pheno[, -c(1:3)]
  dead = pheno #duplicate
  rownames(dead) = paste0(rep('dead_', n_pheno), rownames(pheno))
  pheno[, 'LIVE DEAD Blue'] = rep(0, n_pheno) #LD expr for live cells is ~0
  pheno = pheno[, rownames(spectra)]
  
  exp_prob_mtx_normal = pheno[rep(1:n_pheno , round(counts_pheno * n *
                                                      (1-dead_fraction))), ]
  if(dead_fraction==0){
    full_prob_mtx=exp_prob_mtx_normal
  }
  else{
  cnts_dead = counts_pheno * runif(n_pheno, 0.1, 0.25) 
  cnts_dead = cnts_dead / sum(cnts_dead)
  exp_prob_mtx_dead = dead[rep(1:n_pheno , round(cnts_dead * n * dead_fraction)), ]
  dead_exp=runif(nrow(exp_prob_mtx_dead), 0.3, 0.6)
  exp_prob_mtx_dead = exp_prob_mtx_dead * dead_exp 
  dead_exp = (1-dead_exp) + runif(nrow(exp_prob_mtx_dead), 0.3, 0.5)
  dead_exp[dead_exp>1]=1
  exp_prob_mtx_dead[, 'LIVE DEAD Blue']=dead_exp
  exp_prob_mtx_dead = exp_prob_mtx_dead[, rownames(spectra)]
  full_prob_mtx = rbind(exp_prob_mtx_normal, exp_prob_mtx_dead)
  }
  nr = nrow(full_prob_mtx) 
  weak = (1:ncol(full_prob_mtx))[!(1:ncol(full_prob_mtx) %in% stronk)]
  k_stronk = length(stronk)
  k_weak = nrow(spectra) - k_stronk
  
  exprs_stronk =
     10 ^ ((4-uncertainty) + uncertainty * (matrix(
      rbinom(nr * k_stronk, 1, as.matrix(full_prob_mtx[, stronk])), nr, k_stronk
    ) + rnorm(nr * k_stronk, sd = 0.2)))
  exprs_weak =
    10 ^ (2 + 2 * (matrix(
      rbinom(nr * k_weak, 1, as.matrix(full_prob_mtx[, weak])), nr, k_weak
    ) + rnorm(nr * k_weak, sd = 0.2)))
  exprs = matrix(,
                 nrow = nrow(full_prob_mtx),
                 ncol = ncol(full_prob_mtx))
  exprs[, stronk] = exprs_stronk
  exprs[, weak] = exprs_weak
  rownames(exprs) = rownames(full_prob_mtx)
  
  gt = exprs 
 for (row in 1:nrow(exprs)) {
    exprs[row, ] = exprs[row, ] + rnorm(ncol(exprs), sd = sum(exprs[row, ]) *
                                          0.001)}
  colnames(exprs) = colnames(full_prob_mtx) 
  rm(full_prob_mtx)
  
  emitted = exprs %*% spectra
  
  received = emitted + rnorm(length(emitted), sd = 0.0005 * sqrt(rowSums(emitted ^
                                                                           2)))
  colnames(gt) = rownames(spectra)
  return(list(received, gt, nr))
  
}

generate_artificial_cytometry_data_scaled = function(spectra, pheno, n, stronk, uncertainty=1,dead_fraction=0.175) {
  spectra = spectra / sqrt(rowSums(spectra ^ 2)) 
  d = ncol(spectra) 
  n_pheno = nrow(pheno)
  pheno[is.na(pheno)] = 0
  rownames(pheno) = paste0(pheno$parent_pop, pheno$pop)
  bkp = rownames(pheno)
  pheno$relative_cnt = pheno$relative_cnt * (1 / sum(pheno$relative_cnt))
  counts_pheno = pheno$relative_cnt
  pheno = pheno[, -c(1:3)]
  dead = pheno #duplicate
  rownames(dead) = paste0(rep('dead_', n_pheno), rownames(pheno))
  pheno[, 'LIVE DEAD Blue'] = rep(0, n_pheno) #LD expr for live cells is ~0
  pheno = pheno[, rownames(spectra)]
  
  exp_prob_mtx_normal = pheno[rep(1:n_pheno , round(counts_pheno * n *
                                                      (1-dead_fraction))), ]
  if(dead_fraction==0){
    full_prob_mtx=exp_prob_mtx_normal
  }
  else{
    cnts_dead = counts_pheno * runif(n_pheno, 0.1, 0.25) 
    cnts_dead = cnts_dead / sum(cnts_dead)
    exp_prob_mtx_dead = dead[rep(1:n_pheno , round(cnts_dead * n * dead_fraction)), ]
    dead_exp=runif(nrow(exp_prob_mtx_dead), 0.3, 0.6)
    exp_prob_mtx_dead = exp_prob_mtx_dead * dead_exp 
    dead_exp = (1-dead_exp) + runif(nrow(exp_prob_mtx_dead), 0.3, 0.5)
    dead_exp[dead_exp>1]=1
    exp_prob_mtx_dead[, 'LIVE DEAD Blue']=dead_exp
    exp_prob_mtx_dead = exp_prob_mtx_dead[, rownames(spectra)]
    full_prob_mtx = rbind(exp_prob_mtx_normal, exp_prob_mtx_dead)
  }
  nr = nrow(full_prob_mtx) 
  weak = (1:ncol(full_prob_mtx))[!(1:ncol(full_prob_mtx) %in% stronk)]
  k_stronk = length(stronk)
  k_weak = nrow(spectra) - k_stronk
  
  exprs_stronk =
    2.2* 10 ^ ((4-uncertainty) + uncertainty * (matrix(
      rbinom(nr * k_stronk, 1, as.matrix(full_prob_mtx[, stronk])), nr, k_stronk
    ) + rnorm(nr * k_stronk, sd = 0.15)))
  exprs_stronk =
    2.5*10 ^ ((4-uncertainty) + uncertainty * (matrix(
      rbinom(nr * k_stronk, 1, as.matrix(full_prob_mtx[, stronk])), nr, k_stronk
    ) + rnorm(nr * k_stronk, sd = 0.2)))
  exprs_weak =
    2.5*10 ^ (2 + 2 * (matrix(
      rbinom(nr * k_weak, 1, as.matrix(full_prob_mtx[, weak])), nr, k_weak
    ) + rnorm(nr * k_weak, sd = 0.2)))
  exprs = matrix(,
                 nrow = nrow(full_prob_mtx),
                 ncol = ncol(full_prob_mtx))
  exprs[, stronk] = exprs_stronk
  exprs[, weak] = exprs_weak
  rownames(exprs) = rownames(full_prob_mtx)
  
  gt = exprs 
  for (row in 1:nrow(exprs)) {
    exprs[row, ] = exprs[row, ] + rnorm(ncol(exprs), sd = sum(exprs[row, ]) *
                                          0.002)}
  colnames(exprs) = colnames(full_prob_mtx) 
  rm(full_prob_mtx)
  
  emitted = exprs %*% spectra
  
  received = emitted + rnorm(length(emitted), sd = 0.001 * sqrt(rowSums(emitted ^
                                                                          2)))
  colnames(gt) = rownames(spectra)
  return(list(received, gt, nr))
  
}

get_data_struct = function(spectra, pheno, n, stronk) {
  data = generate_artificial_cytometry_data(spectra, pheno, n, stronk)
  res = list()
  spec_nms = row.names(spectra)
  nspectra = nrow(spectra)
  res[['received']] = data[[1]]
  res[['gt_trans']] = as.data.frame(trans(data[[2]]))
  colnames(res[['gt_trans']]) = spec_nms
  res[['n']] = data[[3]]
  rm(data)
  res[['cell_brightness']] = rowSums(res[['received']])
  res[['ngd']] = list()
  res[['ngd']][['unmix_trans']] = as.data.frame(trans(
    nougadmt(
      res[['received']],
      spectra = spectra,
      snw = 11, #18
      spw = 105, #276
      nw = 358, #1996
      start = 66, #210
      iters = 500, 
      alpha = 0.015, #0.001
      accel=0.085, #0.710
      nthreads = detectCores()
    )$unmixed
  ))
  res[['ols']] = list()
  res[['ols']][['unmix_trans']] = as.data.frame(trans(t(lm(
    t(res[['received']]) ~ t(spectra) + 0
  )$coefficients)))
  res[['wols']] = list()
  res[['wols']][['unmix_trans']] = res[['ols']][['unmix_trans']]
  for (i in 1:res[['n']]) {
    ws = norm(pmax(0, res[['received']][i, ] / sum(res[['received']][i, ])))
    res[['wols']][['unmix_trans']][i, ] = trans(t(lm(t(res[['received']][i, , drop = F]) ~ t(spectra) + 0, weights =
                                                       ws)$coefficients))
  }
  for (toclamp in c('ngd', 'wols', 'ols'))
  {
    clp = paste0('zeroclamp_', toclamp)
    res[[clp]] = list()
    res[[clp]][['unmix_trans']] = as.data.frame(res[[toclamp]][['unmix_trans']])
    res[[clp]][['unmix_trans']][res[[clp]][['unmix_trans']] < 0] = 0
  }
  for (method in c('ngd',
                   'ols',
                   'wols',
                   'zeroclamp_ngd',
                   'zeroclamp_ols',
                   'zeroclamp_wols')) {
    res[[method]][['sq_err']] = as.data.frame((res[[method]][['unmix_trans']] - res[['gt_trans']]) ^
                                                2)
    colnames(res[[method]][['sq_err']]) = spec_nms
    res[[method]][['mean_cell_err']] = rowSums(res[[method]][['sq_err']]) /
      nspectra
    res[[method]][['mean_spec_err']] = colSums(res[[method]][['sq_err']]) /
      res[['n']]
    res[[method]][['mse']] = sum(res[[method]][['sq_err']]) / (res[['n']] *
                                                                 nspectra)
    colnames(res[[method]][['unmix_trans']]) = spec_nms
  }
  return(res)
}