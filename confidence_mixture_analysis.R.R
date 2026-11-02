# Confidence mixture commentary: reproducible analyses
# Author: Kym Michelle McCormick
# Version: 3.2 | 2026-11-02 | R with recommended/base packages only
#
#
# MODEL: independent observations / MAX of raw evidence; independent faces.
# Innocent evidence N(0,1); culprit N(mu,sigma^2). EV sigma=1; UV sigma free.
# c is an absolute noise-referenced threshold, NOT centred SDT response bias.
# Theoretical comparison: SIX faces; equal lineup base rates; fixed c=2 versus
# c chosen to achieve exactly 95% population PPV. No Monte Carlo simulation.
# Empirical fits: EIGHT faces; five free ordered criteria in each of six groups;
# c5 governs confidence >80. UV mu is a noise-standardised mean, not da.
# Empirical Table 2 PPV uses ACTUAL group-specific lineup allocation, NOT 50:50.
# Innocent-suspect rates use CA_any/8, assuming exchangeable innocent faces.
# CIs are observed-Hessian Wald/delta intervals; no bootstrap/equivalence test.
# G2 p-values are asymptotic; inspect sparse-cell and convergence diagnostics.
#
# OUTPUTS: table1_theoretical.csv, table2_observed.csv, table3_parameters_fit.csv,
# table4_criterion_differences.csv, observed_counts.csv, cell_diagnostics.csv,
# fit_checkpoint.rds, analysis_settings.txt, session_info.txt, input_md5.txt.
# Tables contain unrounded values; multiply probability columns by 100 for %.
# Figures, simulation-first-crossing searches, and obsolete two-group variants
# are omitted because they are not needed to reproduce the manuscript tables.
# Existing output files with these names will be replaced on rerunning.
#
# References:
# McCormick, K. M. (in press). Qualitative constraints on models of eyewitness
# identification. Psychological Review. Appendix A.
# Wixted, J. T., Vul, E., Mickes, L., & Wilson, B. M. (2018).
# Models of lineup memory. Cognitive Psychology, 105, 81-114.
# https://doi.org/10.1016/j.cogpsych.2018.06.001
# Pennekamp, P., Moore, K. N., Mansour, J. K., & Lampinen, J. M. (2026).
# When confident eyewitnesses are mistaken: A mixture model account for the
# confidence-accuracy relationship. Memory, 34(10), 1265-1275.
# https://doi.org/10.1080/09658211.2026.2729472
# This is an alternative-policy comparison, not an implementation of the mixture.
#
# MIT License
# Copyright (c) 2026 Kym Michelle McCormick
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

din_2 <- read.csv("din_2.csv", stringsAsFactors = FALSE)

analysis_cells <- c(paste0("CP_culprit_", 1:5), paste0("CP_filler_", 1:5),
                    "CP_reject", paste0("CA_any_", 1:5), "CA_reject")

participant_counts <- function(data) {
  required <- c("memory", "expectation", "target", "Test_T1_suspectIdentified",
                "confidence_rating")
  dat <- as.data.frame(data)
  if (!all(required %in% names(dat))) stop("Missing columns: ",
      paste(setdiff(required, names(dat)), collapse=", "))
  for (nm in required[1:4]) dat[[nm]] <- as.character(dat[[nm]])
  if (anyNA(dat[required[1:4]]) || any(vapply(dat[required[1:4]],
      function(x) any(!nzchar(trimws(x))), logical(1)))) stop("Missing condition/response codes.")
  if (!all(dat$memory %in% c("S","W","VW")) ||
      !all(dat$expectation %in% c("H","L")) || !all(dat$target %in% c("P","A")))
    stop("Unexpected condition codes. Supply only the final S/W/VW x H/L x P/A sample.")
  choose <- dat$Test_T1_suspectIdentified != "Silhouette"
  confidence <- suppressWarnings(as.numeric(as.character(dat$confidence_rating)))
  if (any(!is.finite(confidence[choose])) ||
      any(confidence[choose] < 0 | confidence[choose] > 100))
    stop("Selected faces require finite numeric confidence from 0 to 100.")
  bin <- cut(confidence, c(-Inf,20,40,60,80,Inf), labels=FALSE)
  cell <- rep(NA_character_, nrow(dat))
  cell[!choose & dat$target=="P"] <- "CP_reject"
  cell[!choose & dat$target=="A"] <- "CA_reject"
  ix <- choose & dat$target=="P" & dat$Test_T1_suspectIdentified=="F68"
  cell[ix] <- paste0("CP_culprit_", bin[ix])
  ix <- choose & dat$target=="P" & dat$Test_T1_suspectIdentified!="F68"
  cell[ix] <- paste0("CP_filler_", bin[ix])
  ix <- choose & dat$target=="A"
  cell[ix] <- paste0("CA_any_", bin[ix])
  counts <- as.data.frame(table(
    memory=factor(dat$memory, levels=c("S","W","VW")),
    expectation=factor(dat$expectation, levels=c("H","L")),
    cell=factor(cell, levels=analysis_cells)), stringsAsFactors=FALSE)
  names(counts)[4] <- "observed"
  stopifnot(sum(counts$observed)==nrow(dat))
  counts
}

validate_counts <- function(counts) {
  req <- c("memory","expectation","cell","observed")
  if (!all(req %in% names(counts))) stop("Count input requires: ",paste(req,collapse=", "))
  counts <- counts[req]
  for (nm in req[1:3]) counts[[nm]] <- as.character(counts[[nm]])
  if (anyNA(counts[req[1:3]])) stop("Missing count labels.")
  if (!is.numeric(counts$observed) || any(!is.finite(counts$observed)) ||
      any(counts$observed<0 | counts$observed!=floor(counts$observed)))
    stop("Counts must be finite nonnegative integers.")
  expected <- expand.grid(memory=c("S","W","VW"), expectation=c("H","L"),
                          cell=analysis_cells, stringsAsFactors=FALSE)
  key <- function(x) paste(x$memory,x$expectation,x$cell,sep="|")
  if (anyDuplicated(key(counts)) || nrow(counts)!=nrow(expected) ||
      !setequal(key(counts),key(expected)))
    stop("Supply exactly one row for each of 17 cells in all six groups, including zeros.")
  counts[match(key(expected),key(counts)), ]
}

# Stable numerical integration for tiny tail probabilities at weak memory.
exact_rates <- function(dprime, criterion, m=6L) {
  anchor <- max(criterion,dprime)
  scaled <- integrate(function(x) exp(dnorm(x,dprime,1,log=TRUE)-
    dnorm(anchor,dprime,1,log=TRUE)+(m-1)*pnorm(x,log.p=TRUE)),
    criterion,Inf,rel.tol=1e-10,abs.tol=1e-12,subdivisions=1000L)$value
  logH <- log(scaled)+dnorm(anchor,dprime,1,log=TRUE)
  logPhi <- pnorm(criterion,log.p=TRUE)
  logF <- log(-expm1(m*logPhi))-log(m)
  H <- exp(logH); F <- exp(logF)
  c(H=H,F=F,PPV_equal=plogis(logH-logF),suspect_frequency_equal=.5*(H+F))
}

theoretical_comparison <- function(out_dir) {
  ds <- c(0,seq(.4,2,.2))
  rows <- lapply(ds,function(d) {
    fixed <- exact_rates(d,2)
    criterion <- if(d==0) NA_real_ else uniroot(function(c)
      exact_rates(d,c)["PPV_equal"]-.95,c(-8,12),tol=1e-9)$root
    adjusted <- if(is.na(criterion)) rep(NA_real_,4) else exact_rates(d,criterion)
    data.frame(dprime=d,lineup_size=6L,fixed_c=2,
      H_fixed=unname(fixed[1]),F_fixed=unname(fixed[2]),PPV_fixed=unname(fixed[3]),
      suspect_frequency_fixed=unname(fixed[4]),c_PPV95=criterion,
      H_PPV95=unname(adjusted[1]),F_PPV95=unname(adjusted[2]),
      PPV_PPV95=unname(adjusted[3]),suspect_frequency_PPV95=unname(adjusted[4]))
  })
  result <- do.call(rbind,rows)
  z <- exact_rates(0,2)
  stopifnot(abs(z["H"]-z["F"])<1e-10,abs(z["PPV_equal"]-.5)<1e-10,
    diff(range(result$F_fixed))<1e-12,
    all(abs(result$PPV_PPV95[-1]-.95)<1e-7),
    all(diff(result$H_fixed)>0),all(diff(result$PPV_fixed)>0))
  write.csv(result,file.path(out_dir,"table1_theoretical.csv"),row.names=FALSE)
  result
}

empirical_sdt <- function(counts, out_dir,
                          memories=c("S","W","VW"),expectations=c("H","L")) {
  counts <- validate_counts(counts)
  start_time <- proc.time()[["elapsed"]]
  unpack <- function(p, model) list(c = c(p[1], p[1] + cumsum(exp(p[2:5]))),
                                  mu = p[6], sigma = if (model == "EVSDT") 1 else exp(p[7]))
  cells <- c(paste0("CP_culprit_",1:5), paste0("CP_filler_",1:5), "CP_reject",
             paste0("CA_any_",1:5), "CA_reject")
  probabilities <- function(p, model) {
    u <- unpack(p, model)
    if (any(!is.finite(c(u$c,u$mu,u$sigma))) || u$sigma < .05 || u$sigma > 10 ||
        any(abs(u$c) > 15) || abs(u$mu) > 15) stop("Extreme parameter proposal")
    bounds <- c(u$c, Inf)
    h <- f <- a <- numeric(5)
    for (j in 1:5) {
      h[j] <- integrate(function(x) dnorm(x,u$mu,u$sigma)*pnorm(x)^7,
                         bounds[j],bounds[j+1],rel.tol=1e-7,abs.tol=1e-10)$value
      f[j] <- integrate(function(x) 7*dnorm(x)*pnorm(x)^6*pnorm(x,u$mu,u$sigma),
                         bounds[j],bounds[j+1],rel.tol=1e-7,abs.tol=1e-10)$value
      a[j] <- exp(8*pnorm(bounds[j+1],log.p=TRUE)) - exp(8*pnorm(bounds[j],log.p=TRUE))
    }
    pr <- c(h,f,pnorm(u$c[1],u$mu,u$sigma)*pnorm(u$c[1])^7,a,pnorm(u$c[1])^8)
    if (any(!is.finite(pr)) || any(pr < 0) ||
        max(abs(c(sum(pr[1:11]),sum(pr[12:17]))-1)) > 1e-5) stop("Invalid probabilities")
    pr
  }
  fits <- list(); estimates <- diagnostics <- residuals <- counts_out <- observed_rows <- list()
  for (ex in expectations) for (mem in memories) {
    g <- counts[counts$memory==mem & counts$expectation==ex, ]
    y <- g$observed[match(cells,g$cell)]
    ncp <- sum(y[1:11]); nca <- sum(y[12:17]); ncell <- c(rep(ncp,11),rep(nca,6))
    counts_out[[length(counts_out)+1L]] <- data.frame(memory=mem,expectation=ex,cell=cells,observed=y)
    if(ncp<=0 || nca<=0) stop("Both lineup types require positive totals: ",mem,"_",ex)
    observed_H <- y[5]/ncp; observed_F <- y[16]/nca/8
    denominator <- y[5]+y[16]/8
    observed_rows[[length(observed_rows)+1L]] <- data.frame(memory=mem,expectation=ex,
      n_CP=ncp,n_CA=nca,culprit_high_count=y[5],CA_any_high_count=y[16],
      estimated_innocent_high_count=y[16]/8,H_high=observed_H,F_high=observed_F,
      PPV_experimental=if(denominator>0) y[5]/denominator else NA_real_)
    ev_start <- NULL
    for (model in c("EVSDT","UVSDT")) {
      key <- paste(model,mem,ex,sep="_")
      objective <- function(p) {
        pr <- tryCatch(probabilities(p,model),error=function(e) NULL)
        if (is.null(pr)) return(1e50)
        -sum(y*log(pmax(pr,1e-300)))
      }
      base <- c(.3,rep(log(.5),4),1.5)
      starts <- if (model == "EVSDT") lapply(c(.7,1.5,2.5), function(d) {p<-base;p[6]<-d;p}) else
        lapply(c(.7,1,1.4), function(s) c(ev_start,log(s)))
      runs <- list()
      for (i in seq_along(starts)) {
        cat(sprintf("%s: start %d/%d; elapsed %.1f s\n",key,i,length(starts),proc.time()[["elapsed"]]-start_time)); flush.console()
        runs[[i]] <- tryCatch(optim(starts[[i]],objective,method="BFGS",
          control=list(maxit=600,reltol=1e-9)),error=function(e) NULL)
      }
      valid <- Filter(function(x) !is.null(x) && is.finite(x$value) && x$value < 1e49,runs)
      converged <- Filter(function(x) x$convergence == 0,valid)
      if (!length(converged)) stop("No converged fit for ",key,". Completed earlier fits were checkpointed.")
      best <- converged[[which.min(vapply(converged,function(x)x$value,numeric(1)))]]
      lower_unconverged <- any(vapply(valid,function(x)x$value < best$value-1e-4,logical(1)))
      if (lower_unconverged) warning(key,": an unconverged run has a lower objective; inspect fit.")
      if (model == "EVSDT") ev_start <- best$par
      u <- unpack(best$par,model); pr <- probabilities(best$par,model); expected <- ncell*pr
      hs <- tryCatch(optimHess(best$par,objective),error=function(e) NULL)
      covariance <- NULL; condition <- NA_real_
      if (!is.null(hs) && all(is.finite(hs))) {
        eig <- eigen((hs+t(hs))/2,symmetric=TRUE,only.values=TRUE)$values
        if (min(eig)>0) {
          condition <- max(eig)/min(eig)
          if (condition < 1e10) covariance <- tryCatch(solve(hs),error=function(e) NULL)
        }
      }
      gradient_c5 <- c(1,exp(best$par[2:5]),0,if(model=="UVSDT") 0)
      c5se <- if (is.null(covariance)) NA_real_ else sqrt(drop(t(gradient_c5)%*%covariance%*%gradient_c5))
      muse <- if (is.null(covariance)) NA_real_ else sqrt(covariance[6,6])
      sigse <- if (model=="EVSDT") NA_real_ else if(is.null(covariance)) NA_real_ else u$sigma*sqrt(covariance[7,7])
      H <- pr[5]; F <- pr[16]/8
      estimates[[length(estimates)+1L]] <- data.frame(model,memory=mem,expectation=ex,
        mu=u$mu,mu_se=muse,sigma=u$sigma,sigma_se=sigse,
        c1=u$c[1],c2=u$c[2],c3=u$c[3],c4=u$c[4],c5=u$c[5],c5_se=c5se,
        c5_lower=u$c[5]-1.96*c5se,c5_upper=u$c[5]+1.96*c5se,
        H_high=H,F_high=F,PPV_equal_base_rates=H/(H+F),PPV_experimental=(ncp*H)/(ncp*H+nca*F),
        observed_H_high=y[5]/ncp,observed_F_high=y[16]/nca/8)
      positive <- y>0
      G2 <- 2*sum(y[positive]*log(y[positive]/pmax(expected[positive],1e-300)))
      df <- 15-length(best$par)
      diagnostics[[length(diagnostics)+1L]] <- data.frame(model,memory=mem,expectation=ex,
        n_CP=ncp,n_CA=nca,nll=best$value,k=length(best$par),AIC=2*best$value+2*length(best$par),
        G2,df,p_asymptotic=pchisq(G2,df,lower.tail=FALSE),
        cells_expected_below5=sum(expected<5),minimum_expected=min(expected),hessian_condition=condition,
        intervals_available=!is.null(covariance),converged_starts=length(converged),
        lower_unconverged_run=lower_unconverged)
      residuals[[length(residuals)+1L]] <- data.frame(model,memory=mem,expectation=ex,
        cell=cells,observed=y,expected,pearson=(y-expected)/sqrt(pmax(expected,1e-300)))
      fits[[key]] <- list(fit=best,parameters=u,covariance=covariance,counts=y,probabilities=pr,runs=runs)
      saveRDS(list(fits=fits,estimates=estimates,diagnostics=diagnostics),file.path(out_dir,"fit_checkpoint.rds"))
    }
  }
  est <- do.call(rbind,estimates); diag <- do.call(rbind,diagnostics)
  diffs <- list()
  for (model in c("EVSDT","UVSDT")) for (ex in expectations) {
    e <- est[est$model==model & est$expectation==ex,]
    # Compare each included memory condition with the first (default S).
    ref <- e[e$memory==memories[1],]
    for (mem in memories[-1]) {
      q <- e[e$memory==mem,]; difference <- q$c5-ref$c5
      se <- sqrt(q$c5_se^2+ref$c5_se^2)
      diffs[[length(diffs)+1L]] <- data.frame(model,expectation=ex,
        contrast=paste(mem,"minus",memories[1]),c5_difference=difference,se,
        lower=difference-1.96*se,upper=difference+1.96*se)
    }
  }
  difference_table <- do.call(rbind,diffs)
  observed_table <- do.call(rbind,observed_rows)
  combined <- merge(est,diag,by=c("model","memory","expectation"),sort=FALSE)
  order_key <- paste(rep(expectations,each=length(memories)),rep(memories,times=length(expectations)))
  combined <- combined[order(match(combined$model,c("EVSDT","UVSDT")),
      match(paste(combined$expectation,combined$memory),order_key)), ]
  write.csv(observed_table,file.path(out_dir,"table2_observed.csv"),row.names=FALSE)
  write.csv(combined,file.path(out_dir,"table3_parameters_fit.csv"),row.names=FALSE)
  write.csv(diag,file.path(out_dir,"absolute_fit.csv"),row.names=FALSE)
  write.csv(difference_table,file.path(out_dir,"table4_criterion_differences.csv"),row.names=FALSE)
  write.csv(do.call(rbind,residuals),file.path(out_dir,"cell_diagnostics.csv"),row.names=FALSE)
  write.csv(do.call(rbind,counts_out),file.path(out_dir,"observed_counts.csv"),row.names=FALSE)
  if(any(diag$cells_expected_below5>0))
    warning("Some expected cells are below five; inspect absolute_fit.csv and cell_diagnostics.csv. G2 p-values are asymptotic.")
  if(any(!diag$intervals_available)) warning("Some Hessian intervals unavailable; inspect diagnostics.")
  message("Empirical fits completed. Wald intervals assume independent groups; nonsignificance is not equivalence.")
  invisible(list(observed=observed_table,estimates=est,fit=diag,
                 differences=difference_table,fits=fits))
}

run_analysis <- function(input=NULL,out_dir="confidence_mixture_results",data=NULL) {
  dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
  if(!dir.exists(out_dir)) stop("Cannot create output directory.")
  settings <- c("Theoretical lineup size: 6; empirical lineup size: 8.",
    "Theory: independent EVSDT; c=2 or exact population PPV=.95; equal base rates.",
    "Empirical: six separate groups; EVSDT and UVSDT; confidence bins 0-20/21-40/41-60/61-80/>80 (continuous ratings use open-left upper bins).",
    "Observed Table 2 PPV: actual condition-specific lineup allocation; false IDs estimated as CA selections / 8.",
    "No random sampling is used; no RNG seed is needed.",
    "All reported probabilities in CSVs are proportions, not percentages.",
    "Wald/delta CIs; asymptotic G2 tests; no bootstrap or equality/equivalence test.")
  writeLines(settings,file.path(out_dir,"analysis_settings.txt"))
  on.exit(capture.output(sessionInfo(),file=file.path(out_dir,"session_info.txt")))
  theory <- theoretical_comparison(out_dir)
  if(is.null(input) && is.null(data)) {
    if(file.exists("observed_counts.csv")) input <- "observed_counts.csv"
    else if(file.exists("data.csv")) input <- "data.csv"
  }
  if(!is.null(input)) {
    if(!is.null(data)) stop("Supply input OR data, not both.")
    if(!file.exists(input)) stop("Input file not found: ",input)
    data <- read.csv(input,stringsAsFactors=FALSE,check.names=FALSE)
    capture.output(tools::md5sum(input),file=file.path(out_dir,"input_md5.txt"))
  }
  if(is.null(data)) {
    message("Theory completed. Empirical analysis SKIPPED: supply observed_counts.csv, a participant CSV, or din_2.")
    return(invisible(list(theory=theory,empirical=NULL)))
  }
  counts <- if(all(c("cell","observed") %in% names(data))) validate_counts(data) else participant_counts(data)
  empirical <- empirical_sdt(counts,out_dir)
  message("Outputs: ",normalizePath(out_dir))
  invisible(list(theory=theory,empirical=empirical))
}

# Set options(confidence_mixture.autorun=FALSE) before source() to load functions
# without running. Then call run_analysis(input="observed_counts.csv") explicitly.
if(isTRUE(getOption("confidence_mixture.autorun",TRUE))) {
  args <- if(sys.nframe()==0L) commandArgs(trailingOnly=TRUE) else character()
  if(length(args)>2L) stop("Usage: Rscript confidence_mixture_analysis.R [input.csv] [output_directory]")
  input_path <- if(length(args)>=1L) args[1] else NULL
  output_path <- if(length(args)>=2L) args[2] else "confidence_mixture_results"
  # Explicit files override din_2. Sourcing with din_2 preserves the original workflow.
  existing_data <- if(is.null(input_path)) get0("din_2",inherits=TRUE) else NULL
  confidence_mixture_results <- run_analysis(input=input_path,out_dir=output_path,data=existing_data)
}

# Figure 1: fixed criterion versus criterion maintaining 95% PPV
# Six-member lineups; equal target-present/target-absent base rates.

make_theoretical_figure <- function(
    out_dir = "confidence_mixture_results") {
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  # No finite 95%-PPV criterion exists at d' = 0.
  # Evaluate that policy from d' = 0.4, matching Table 1.
  fixed_d <- seq(0, 2, length.out = 201)
  adjusted_d <- seq(0.4, 2, length.out = 161)
  
  fixed <- do.call(rbind, lapply(fixed_d, function(d) {
    rates <- exact_rates(d, criterion = 2, m = 6L)
    data.frame(
      dprime = d,
      criterion = 2,
      H = unname(rates["H"]),
      F = unname(rates["F"]),
      PPV = unname(rates["PPV_equal"])
    )
  }))
  
  adjusted <- do.call(rbind, lapply(adjusted_d, function(d) {
    criterion <- uniroot(
      function(c) {
        exact_rates(d, criterion = c, m = 6L)["PPV_equal"] - 0.95
      },
      interval = c(-8, 12),
      tol = 1e-9
    )$root
    
    rates <- exact_rates(d, criterion = criterion, m = 6L)
    
    data.frame(
      dprime = d,
      criterion = criterion,
      H = unname(rates["H"]),
      F = unname(rates["F"]),
      PPV = unname(rates["PPV_equal"])
    )
  }))
  
  stopifnot(
    all(abs(adjusted$PPV - 0.95) < 1e-7),
    diff(range(fixed$F)) < 1e-12
  )
  
  draw_figure <- function() {
    old_par <- par(no.readonly = TRUE)
    on.exit(par(old_par))
    
    # Shared legend above four panels.
    layout(
      matrix(c(1, 1, 2, 3, 4, 5), nrow = 3, byrow = TRUE),
      heights = c(0.45, 3, 3)
    )
    
    colours <- c("#0072B2", "#D55E00")
    line_types <- c(1, 2)
    
    par(mar = c(0, 0, 0, 0))
    plot.new()
    legend(
      "center",
      legend = c(
        "Fixed criterion (c = 2)",
        "Criterion maintaining 95% PPV"
      ),
      col = colours,
      lty = line_types,
      lwd = 2.2,
      horiz = TRUE,
      bty = "n",
      cex = 0.95
    )
    
    par(
      mar = c(4.3, 4.8, 2.5, 1),
      mgp = c(2.9, 0.8, 0),
      las = 1,
      family = "sans",
      cex = 0.95
    )
    
    panel <- function(field, title, ylabel, ylim,
                      multiplier = 1, reference = NULL) {
      
      plot(
        NA,
        xlim = c(0, 2),
        ylim = ylim,
        xaxs = "i",
        yaxs = "i",
        xaxt = "n",
        xlab = expression("Memory strength (" * d * prime * ")"),
        ylab = ylabel,
        bty = "l"
      )
      
      axis(1, at = seq(0, 2, 0.4))
      
      if (!is.null(reference)) {
        abline(h = reference, col = "grey70", lty = 3)
      }
      
      lines(
        fixed$dprime,
        multiplier * fixed[[field]],
        col = colours[1],
        lty = line_types[1],
        lwd = 2.2
      )
      
      lines(
        adjusted$dprime,
        multiplier * adjusted[[field]],
        col = colours[2],
        lty = line_types[2],
        lwd = 2.2
      )
      
      title(main = title, adj = 0, cex.main = 1)
    }
    
    panel(
      "criterion",
      "A  High-confidence criterion",
      "Criterion (innocent-evidence SD units)",
      ylim = c(0, 8)
    )
    
    panel(
      "H",
      "B  Culprit identification",
      "High-confidence culprit ID (%)",
      ylim = c(0, 60),
      multiplier = 100
    )
    
    panel(
      "F",
      "C  Innocent-suspect identification",
      "High-confidence innocent-suspect ID (%)",
      ylim = c(0, 3),
      multiplier = 100
    )
    
    panel(
      "PPV",
      "D  Positive predictive value",
      "High-confidence suspect PPV (%)",
      ylim = c(45, 100),
      multiplier = 100,
      reference = 95
    )
  }
  
  save_plot <- function(path, format) {
    if (format == "pdf") {
      pdf(
        path,
        width = 9,
        height = 7,
        family = "Helvetica",
        useDingbats = FALSE
      )
    } else {
      tiff(
        path,
        width = 9,
        height = 7,
        units = "in",
        res = 600,
        compression = "lzw"
      )
    }
    
    on.exit(dev.off())
    draw_figure()
  }
  
  save_plot(
    file.path(out_dir, "figure1_theoretical_comparison.pdf"),
    "pdf"
  )
  
  save_plot(
    file.path(out_dir, "figure1_theoretical_comparison.tiff"),
    "tiff"
  )
  
  message("Figure saved in: ", normalizePath(out_dir))
  invisible(list(fixed = fixed, adjusted = adjusted))
}

figure1_data <- make_theoretical_figure()
