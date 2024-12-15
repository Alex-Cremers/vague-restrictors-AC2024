# Load packages and define functions ----

library(tidyverse)
library(XML) # to read degree labels
library(lme4) # for simple exploratory models and prior smoothing
library(mvtnorm) # to generate multivariate normal samples for graphs
library(mclust) # to get density of multivariate (faster than mvtnorm)
library(parallel) # to speed things up
library(xtable) # latex table export
library(loo) # leave-one-out model comparison
library(cubature) # faster integrals
library(cmdstanr) # Stan interface
library(glue) # always need some glue

options(mc.cores = parallel::detectCores()-2L)
options(tibble.width=Inf)
options(dplyr.summarise.inform = FALSE)

# Functions for permutation test:
source('functions/get_all_clusters.R')
source('functions/get_H0_dist.R')
source('functions/get_signif_clusters.R')

# Functions for RSA modelling:
log1mexp <- copula::log1mexp
log1pexp <- copula::log1pexp
logSumExp <- matrixStats::logSumExp
source('functions/mean_log_L0.R')
source('functions/mean_L0.R')
source('functions/compute_utilities.R')
# These two are not actually used in the Stan model, but allow us to 
# compute the model predictions in R if desired:
source('functions/S1.R')
source('functions/L1.R')

# Standard error function for graphs
se <- function(x){sd(x,na.rm=T)/sqrt(length(x[!is.na(x)]))}

# Palettes for graphs
TwoColorPalette = c(rgb(.5,0,.7),rgb(.9,.65,0))
FourColorPalette = c("#F5793A","#A95AA1","#85C0F9","#0F2080")
FiveColorPalette = c("#F5793A","#A95AA1","#85C0F9","#0F2080", "#3EE4A7")

# Function to decode XML labels
decode <- function(x) {
  if(x==""){return("")} else {
    xmlValue(getNodeSet(htmlParse(x, asText = TRUE,encoding="UTF8"), "//p|//span")[[1]])
  }
}

# Prepare a folder for graphs
suppressWarnings(dir.create("graphs/"))

# Read and process the data ----

grad_restr_data <- read_csv("vague_restrictors_anon_results.csv", show_col_types = FALSE)

# Load complete details about each context:
context_info <- read_csv("vague_restrictors_context_info.csv", show_col_types = FALSE)

# Format numeric degrees to merge with main data:
numeric_degrees <- context_info %>%
  select(ContextName,A:H) %>%
  pivot_longer(
    A:H,
    names_to = "Probe",
    values_to = "Degree"
  ) %>%
  rename(Adjective=ContextName)

probes <- toupper(letters[1:8])
names(probes) <- paste0("Deg",1:8)

degree_labels <- context_info %>%
  select(ContextName,starts_with("Deg")) %>%
  pivot_longer(
    starts_with("Deg"),
    names_to = "Probe",
    values_to = "DegreeLabel"
  ) %>%
  rename(Adjective=ContextName) %>%
  mutate(DegreeLabel=if_else(is.na(DegreeLabel),"",DegreeLabel),
         DecodedLabel = as.character(sapply(DegreeLabel,decode)),
         Probe = probes[Probe]
  )
rm(probes)

# Merge information about degrees and normalize everything:
grad_restr_data <- grad_restr_data %>%
  left_join(numeric_degrees, by = c("Adjective", "Probe")) %>%
  mutate(
    Adjective = factor(Adjective,
                       levels=c("tall","powerful","hot","expensive","large","young","late","spicy","profitable","complete","safe","full")),
    AdjClass = if_else(Adjective %in% c("tall","powerful","young","hot","expensive","large"), "Relative", "Absolute"),
    AdjMonotonicity = if_else(Adjective%in%c("young","safe","full"),"Negative","Positive"),
    AdjRating=AdjRating/100,
    Prior=Prior/100,
    Posterior=Posterior/100,
    pol_degree = if_else(AdjMonotonicity=="Positive",Degree,-Degree)
  ) %>%
  mutate(AdjClass = factor(AdjClass, levels = c("Relative", "Absolute")))

# Define errors on extreme elements of the scale
grad_restr_data <- grad_restr_data %>%
  mutate(
    Error = case_when(
      Probe=="B" & (AdjRating >= 0.5)  & !(Adjective%in%c("young","full","safe")) ~ T,
      Probe=="B" & (AdjRating <= 0.5)  &  (Adjective%in%c("young","full","safe")) ~ T,
      Probe=="H" & (AdjRating <= 0.5)  & !(Adjective%in%c("young","full","safe")) ~ T,
      Probe=="H" & (AdjRating >= 0.5)  &  (Adjective%in%c("young","full","safe")) ~ T,
      T ~ F
    )
  )

# Show error rate on each adjective:
grad_restr_data %>%
  filter(Probe %in% c("B", "H")) %>%
  group_by(Adjective) %>%
  summarize(ErrorRate = mean(Error)) %>%
  arrange(-ErrorRate)

# Check what happens with each adjective
grad_restr_data %>%
  na.omit() %>%
  # filter(Probe %in% c("B", "H")) %>%
  ggplot(aes(x = Probe, y = AdjRating, group = Subject)) +
  facet_wrap(.~Adjective) +
  geom_line(alpha = .5) +
  geom_hline(yintercept = 0.5, col = 'red') +
  theme_bw()
# Conclusion: many participant consider even the lower temperatures 'hot'
# On adjectives with inverted scale, we see occasional confusion about monotonicity, which should definitely be excluded.


# Loosen our condition for relative adjectives: as long as rating increase along probes, accept.
final_error <- grad_restr_data %>%
  group_by(Subject, AssignmentId) %>%
  summarize(Error = any(
    (AdjClass[1] == "Absolute" & any(Error)),
    AdjClass[1] == "Relative" & AdjMonotonicity == "Positive" & any(Error) & !(AdjRating[Probe == "B"] <= AdjRating[!Probe %in% c("B", "H")] & AdjRating[!Probe %in% c("B", "H")] <= AdjRating[Probe == "H"]),
    AdjClass[1] == "Relative" & AdjMonotonicity == "Negative" & any(Error) & !(AdjRating[Probe == "B"] >= AdjRating[!Probe %in% c("B", "H")] & AdjRating[!Probe %in% c("B", "H")] >= AdjRating[Probe == "H"])
  ))

# Track participants who left all sliders untouched (didn't actually happen in this experiment)
SliderUntouched <- grad_restr_data %>%
  group_by(Subject, AssignmentId) %>%
  summarize(SliderError = all(AdjRating == 0.5))


# Check number of errors on adjective rating for extreme items:
ParticipantInfo <- grad_restr_data %>%
  group_by(AssignmentId,Subject) %>%
  summarize(AdjError=sum(Error),WorkTimeInSeconds=first(WorkTimeInSeconds), Quantifier = first(Quantifier), Adjective = first(Adjective)) %>%
  select(Subject,AssignmentId,Quantifier,Adjective,AdjError,WorkTimeInSeconds) %>%
  left_join(final_error)

# Proportion of HITs excluded on the basis of adjective rating errors:
mean(ParticipantInfo$AdjError>0)
sum(ParticipantInfo$AdjError>0)

# After new criteria following reviewer's remarks:
mean(ParticipantInfo$Error)
sum(ParticipantInfo$Error)


## Remove excluded HITs ----

grad_restr_data <- grad_restr_data %>%
  filter(AssignmentId%in%ParticipantInfo$AssignmentId[ParticipantInfo$Error==FALSE] & !is.na(Probe))


# First graphs and analyses ----

## Check correlation between Adjective rating and Priors ----

stat_grad_restr_data <- grad_restr_data %>%
  mutate(
    AdjRating=as.numeric(scale(AdjRating)),
    Prior=as.numeric(scale(Prior)),
    Posterior = as.numeric(scale(Posterior))
  )

# A lot of correlations, in various directions.
PriorRatingModel <- lmer(Prior~1+scale(AdjRating):Adjective+(1+AdjRating|Subject)+(1|Adjective),data=stat_grad_restr_data)
summary(PriorRatingModel)

PriorRatingModel0 <- lmer(Prior~1+(1+AdjRating|Subject)+(1|Adjective),data=stat_grad_restr_data)
anova(PriorRatingModel0, PriorRatingModel)

# We cannot assume that the prior is independent from adjective rating, which will complicate modeling a bit.

rm(stat_grad_restr_data)

## Plot the adjective acceptability (Q1) ----

# A lot of setup to make a nice graph with proper x-axis label 
# reflecting the degree labels actually seen by the participants

plot_data <- grad_restr_data %>%
  left_join(select(degree_labels,-DegreeLabel)) %>%
  mutate(Adjective = factor(Adjective,
                            levels=c("tall","powerful","hot","expensive","large","young","late","spicy","profitable","complete","safe","full")),
         DecodedLabel = gsub(" ","",DecodedLabel)
  ) %>% 
  group_by(Adjective) %>%
  mutate(ProbeRank = as.integer(factor(Probe)))


breaks_fun <- function(x) {
  count_breaks <<- count_breaks + 0.5
  switch(
    floor(count_breaks),
    sort(unique(plot_data$ProbeRank[plot_data$Adjective=="tall"])),
    sort(unique(plot_data$ProbeRank[plot_data$Adjective=="powerful"])),
    sort(unique(plot_data$ProbeRank[plot_data$Adjective=="hot"])),
    sort(unique(plot_data$ProbeRank[plot_data$Adjective=="expensive"])),
    sort(unique(plot_data$ProbeRank[plot_data$Adjective=="large"])),
    sort(unique(plot_data$ProbeRank[plot_data$Adjective=="young"])),
    sort(unique(plot_data$ProbeRank[plot_data$Adjective=="late"])),
    sort(unique(plot_data$ProbeRank[plot_data$Adjective=="spicy"])),
    sort(unique(plot_data$ProbeRank[plot_data$Adjective=="profitable"])),
    sort(unique(plot_data$ProbeRank[plot_data$Adjective=="complete"])),
    sort(unique(plot_data$ProbeRank[plot_data$Adjective=="safe"])),
    sort(unique(plot_data$ProbeRank[plot_data$Adjective=="full"]))
  )
}

profitable_labels <- unique(plot_data$DecodedLabel[plot_data$Adjective=="profitable"][order(plot_data$ProbeRank[plot_data$Adjective=="profitable"])])
full_labels <- unique(plot_data$DecodedLabel[plot_data$Adjective=="full"][order(plot_data$ProbeRank[plot_data$Adjective=="full"])])

spicy_scale <- c("","\uD83C\uDF36","\uD83C\uDF36\uD83C\uDF36","\uD83C\uDF36\uD83C\uDF36\uD83C\uDF36\uD83C\uDF36")

labels_fun <- function(x) {
  count_label <<- count_label + 1L
  switch(
    count_label,
    unique(plot_data$DecodedLabel[plot_data$Adjective=="tall"][order(plot_data$ProbeRank[plot_data$Adjective=="tall"])]),
    unique(plot_data$DecodedLabel[plot_data$Adjective=="powerful"][order(plot_data$ProbeRank[plot_data$Adjective=="powerful"])]),
    unique(plot_data$DecodedLabel[plot_data$Adjective=="hot"][order(plot_data$ProbeRank[plot_data$Adjective=="hot"])]),
    unique(plot_data$DecodedLabel[plot_data$Adjective=="expensive"][order(plot_data$ProbeRank[plot_data$Adjective=="expensive"])]),
    unique(plot_data$DecodedLabel[plot_data$Adjective=="large"][order(plot_data$ProbeRank[plot_data$Adjective=="large"])]),
    unique(plot_data$DecodedLabel[plot_data$Adjective=="young"][order(plot_data$ProbeRank[plot_data$Adjective=="young"])]),
    unique(plot_data$DecodedLabel[plot_data$Adjective=="late"][order(plot_data$ProbeRank[plot_data$Adjective=="late"])]),
    spicy_scale,
    profitable_labels,
    unique(plot_data$DecodedLabel[plot_data$Adjective=="complete"][order(plot_data$ProbeRank[plot_data$Adjective=="complete"])]),
    unique(plot_data$DecodedLabel[plot_data$Adjective=="safe"][order(plot_data$ProbeRank[plot_data$Adjective=="safe"])]),
    full_labels
  )
}


zero_data <- tibble(
  Adjective = c("tall","powerful","hot","expensive","large","late","spicy","profitable","complete","young","safe","full"),
  zero = c(NaN, NaN, NaN, NaN, NaN, 1.8, 1, 1.9, 4, NaN, 0, 1)) %>%
  mutate(Adjective = factor(Adjective, levels = levels(plot_data$Adjective)))

# First results graph: (uncomment lines to export to file)
count_breaks <- 1
count_label <- 0
figure1 <- plot_data %>%
  ggplot(aes(y=AdjRating,x=ProbeRank,group=Probe)) +
  facet_wrap(~Adjective,scales = "free_x",dir = "h",ncol=6) +
  geom_boxplot(outlier.shape = NA,width=.5,col=TwoColorPalette[1]) +
  geom_vline(data=na.omit(zero_data),aes(xintercept = zero),linetype=2,color=TwoColorPalette[2])+
  theme_bw()+
  scale_x_continuous(name="Measure",breaks = breaks_fun,minor_breaks = NULL,labels=labels_fun, limits = c(NA, NA))+
  scale_y_continuous(name="Adjective acceptability",labels=scales::percent)+
  # theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5,size=6))+
  theme(axis.text.x = element_text(angle=60,hjust=1,vjust=1,size=6))+
  theme(text = element_text(family = "Arial Unicode MS"))
# quartz(width=10,height=5)
figure1
# ps = gridSVG::grid.export("graphs/Figure1_wide.svg", addClass=T,progress = T, exportJS="inline")


# Compare prior and posterior across everything:
grad_restr_data %>%
  pivot_longer(cols=c(Posterior,Prior),names_to = "Type",values_to = "PredRating") %>%
  mutate(Type = case_when(
    Type == "Prior" ~ "Prior",
    Type == "Posterior" & Quantifier == "Every" ~ "Posterior (every)",
    Type == "Posterior" & Quantifier == "No" ~ "Posterior (no)"
  )) %>%
  ggplot(aes(x=AdjRating,y=PredRating,color=Type,fill=Type))+
  # facet_wrap(.~Adjective,labeller=labeller(Adjective = label_value, Quantifier = label_value,.multi_line = F)) + 
  geom_point() + 
  #geom_abline(intercept = 0,slope=1,linetype=2) + 
  geom_smooth(span=1.5,se=T)+
  theme_bw() + 
  scale_color_manual(values=FourColorPalette)+
  scale_fill_manual(values=FourColorPalette)+
  scale_y_continuous(name="Predicate rating",breaks=seq(0,1,by=.5),minor_breaks = NULL,labels = scales::percent)+
  scale_x_continuous(name="Adjective rating",breaks=seq(0,1,by=.5),minor_breaks = NULL,labels = scales::percent)+
  coord_cartesian(ylim=c(0,1))



## Theory-neutral analysis ----

# Bootstrap method adapted from Maris & Oostenveld
# smoothing added on x-axis as each series is measured at different x-values

# Compute deviation between prior and posterior
dev_data <- grad_restr_data %>%
  arrange(Quantifier, AdjClass, AdjRating) %>%
  mutate(Deviation = Posterior-Prior) %>%
  select(Quantifier, AdjClass, Subject, AdjRating, Deviation, Posterior, Prior)

# Find significant clusters:
tb_clusters <- dev_data %>%
  group_by(Quantifier, AdjClass) %>%
  reframe(get_signif_clusters(AdjRating, Deviation, id = Subject, cores = 10L, N = 1000L))

# Convert back from rank to actual adjective ratings:
AdjRating_indices <- dev_data %>%
  group_by(Quantifier, AdjClass) %>%
  filter(!duplicated(AdjRating)) %>%
  mutate(Index = seq_along(AdjRating)) %>%
  select(Quantifier, AdjClass, Index, AdjRating)

tb_clusters <- tb_clusters %>%
  left_join(AdjRating_indices, by = join_by(Quantifier == Quantifier, AdjClass == AdjClass, start == Index)) %>%
  rename(xmin = AdjRating) %>%
  left_join(AdjRating_indices, by = join_by(Quantifier == Quantifier, AdjClass == AdjClass, end == Index)) %>%
  rename(xmax = AdjRating) %>%
  mutate(
    ymin = ifelse(stat<0, -Inf, 0),
    ymax = ifelse(stat>0, Inf, 0),
  )

# Visualize clusters (uncomment top and bottom line to export pdf)
# pdf("graphs/clusters.pdf", width = 8, height = 6)
grad_restr_data %>%
  pivot_longer(cols=c(Posterior,Prior),names_to = "Type",values_to = "PredRating") %>%
  ggplot(aes(x=AdjRating,y=PredRating,color=Type))+
  facet_grid(AdjClass~Quantifier,labeller=labeller(Quantifier = label_value)) +
  # facet_wrap(.~Adjective,labeller=labeller(Adjective = label_value, Quantifier = label_value,.multi_line = F)) + 
  geom_rect(data = tb_clusters %>% filter(p_val < .05), aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=Inf, fill = factor(sign(stat))), inherit.aes = FALSE) +
  geom_point(alpha = 0.7) + 
  #geom_abline(intercept = 0,slope=1,linetype=2) + 
  geom_smooth(span=0.5,se=T, fill = 'grey')+
  theme_bw() + 
  scale_color_manual(values=FourColorPalette)+
  scale_fill_manual(values = c(rgb(1,0,0,.25), rgb(0,1,0,0.25)), guide = "none") +
  # scale_fill_manual(values=FourColorPalette)+
  scale_y_continuous(name="Predicate rating",breaks=seq(0,1,by=.5),minor_breaks = NULL,labels = scales::percent)+
  scale_x_continuous(name="Adjective rating",breaks=seq(0,1,by=.5),minor_breaks = NULL,labels = scales::percent)+
  coord_cartesian(ylim=c(0,1))
# dev.off()

# Generate a LaTeX table:
tb_clusters %>%
  filter(p_val < .05) %>%
  mutate(
    cluster_range = glue("[{format(xmin, digits = 2)}, {format(xmax, digits = 2)}]"),
    p = ifelse(p_val == 0, "<.001",format(p_val, digits = 3))
  ) %>%
  select(Quantifier, AdjClass, cluster_range, stat, p) %>%
  xtable(align = "cllcrr", digits = 0) %>%
  print(include.rownames=FALSE, booktabs = TRUE)



# Check interaction between Adj Class and Posterior-Prior.
int_data <- grad_restr_data %>%
  mutate(Deviation = ifelse(Quantifier == "Every", Posterior-Prior, Prior-Posterior)) %>%
  group_by(Quantifier, Subject, AdjRating) %>%
  mutate(Id = seq(n())) %>%
  ungroup() %>%
  pivot_wider(id_cols = c(Quantifier, Subject, AdjRating, Id), names_from = AdjClass, values_from = Deviation) %>%
  arrange(Quantifier, AdjRating)

tb_int_clusters <- int_data %>%
  group_by(Quantifier) %>%
  reframe(get_signif_clusters(AdjRating, y=Relative, y2=Absolute, id=Subject, cores = 10L, N = 1000L))

AdjRating_int_indices <- dev_data %>%
  group_by(Quantifier) %>%
  filter(!duplicated(AdjRating)) %>%
  mutate(Index = seq_along(AdjRating)) %>%
  select(Quantifier, Index, AdjRating)

tb_int_clusters <- tb_int_clusters %>%
  left_join(AdjRating_int_indices, by = join_by(start == Index, Quantifier == Quantifier)) %>%
  rename(xmin = AdjRating) %>%
  left_join(AdjRating_int_indices, by = join_by(end == Index, Quantifier == Quantifier)) %>%
  rename(xmax = AdjRating) %>%
  mutate(
    ymin = ifelse(stat<0, -Inf, 0),
    ymax = ifelse(stat>0, Inf, 0),
  )


# Plot the cluster(s) of interaction with Adjective class:
grad_restr_data %>%
  mutate(Deviation = ifelse(Quantifier == "Every", Posterior-Prior, Prior-Posterior)) %>%
  ggplot(aes(x=AdjRating,y=Deviation,color=AdjClass))+
  facet_grid(Quantifier~.) +
  # facet_wrap(.~Adjective,labeller=labeller(Adjective = label_value, Quantifier = label_value,.multi_line = F)) + 
  geom_rect(data = tb_int_clusters %>% filter(p_val < .05), aes(xmin=xmin, xmax=xmax, ymin=-Inf, ymax=Inf, fill = factor(sign(stat))), inherit.aes = FALSE) +
  geom_point(alpha = 0.7) + 
  #geom_abline(intercept = 0,slope=1,linetype=2) + 
  geom_smooth(span=0.5,se=T, fill = 'grey')+
  theme_bw() + 
  scale_color_manual(values=FourColorPalette)+
  scale_fill_manual(values = c(rgb(1,0,0,.25), rgb(0,1,0,0.25)), guide = "none") +
  # scale_fill_manual(values=FourColorPalette)+
  scale_y_continuous(name="Posterior - Prior (flipped for 'no')",breaks=seq(-1,1,by=.5),minor_breaks = NULL,labels = scales::percent)+
  scale_x_continuous(name="Adjective rating",breaks=seq(0,1,by=.5),minor_breaks = NULL,labels = scales::percent)+
  coord_cartesian(ylim=c(-1,1))



# Modelling ----

## Step 1: Fit 1st and 2nd order vagueness ----

# apply lemon squeezer, compute various useful quantities
stat_grad_restr_data <- grad_restr_data %>%
  mutate(
    AdjRating=.005+AdjRating*.99,
    Prior=.005+Prior*.99,
    Posterior=.005+Posterior*.99,
    pol_degree = if_else(AdjMonotonicity=="Positive",Degree,-Degree)
  ) %>%
  group_by(Adjective) %>%
  mutate(
    mean_pol_degree=mean(pol_degree),
    sd_degree=sd(pol_degree),
    pol_norm_degree=(pol_degree-mean_pol_degree)/sd_degree,
    NormDegree=if_else(AdjMonotonicity=="Positive",pol_norm_degree,-pol_norm_degree)) %>%
  ungroup()

# Don't recompute if we already have a saved file, as this takes quite a few minutes
if (file.exists("precomputed_data/Theta_params_gaussian.csv")) {
  cat("Loading existing files\n")
  Theta_param <- read_csv("precomputed_data/Theta_params_gaussian.csv", show_col_types = FALSE) %>% column_to_rownames("Adjective")
  tb_subj_Theta <- read_csv("precomputed_data/Theta_subjects_gaussian.csv", show_col_types = FALSE)
} else {
  cat("Recomputing Theta parameters...\n")
  stan_Theta_model <- cmdstan_model("models/Theta_model.stan")
  
  Theta_param <- tibble(Adjective=as.character(unique(stat_grad_restr_data$Adjective)),m_mu=NA,m_sig=NA,s_mu=NA,s_sig=NA,rho=NA)
  Subj_Theta <- list()
  for (i in 1:nrow(Theta_param)){
    adj=Theta_param$Adjective[i]
    cat(glue("Computing Theta and Omega for {adj}...\n"))
    stan_adj_data <- stat_grad_restr_data %>%
      filter(Adjective==adj) %>%
      mutate(Subject=as.numeric(factor(Subject)))
    data = list(
      N = nrow(stan_adj_data),
      S = n_distinct(stan_adj_data$Subject),
      y = stan_adj_data$AdjRating,
      degree = as.numeric(stan_adj_data$pol_norm_degree),
      subject = stan_adj_data$Subject
    )
    fit_adj <- stan_Theta_model$sample(
      data = data,
      chains = 8,
      parallel_chains = 14,
      iter_warmup = 1000,
      iter_sampling = 2000,
      adapt_delta = 0.95,
      show_messages = F
    )
    # Extract posterior mean of the main parameters:
    param <- fit_adj$summary(c("m_mu","m_sig","s_mu","s_sig","L_u[2,1]")) %>% pull(mean)
    Subj_Theta[[adj]] <- fit_adj$summary(c("mu", "sigma")) %>% 
      select(variable, mean) %>% 
      mutate(subject = sub("(mu|sigma)\\[(\\d+)\\]", "\\2", variable),
             variable = if_else(grepl("mu", variable), "mu", "sigma"))
    Theta_param[i,2:6] <- t(param)
  }
  
  # Get mapping of subject number by adjective
  Subj_num <- stat_grad_restr_data %>%
    group_by(Adjective) %>%
    mutate(Subject_num=as.numeric(factor(Subject))) %>%
    select(Subject, Adjective, Subject_num) %>%
    distinct() %>%
    arrange(Adjective, Subject_num)
  
  tb_subj_Theta <- Theta_param$Adjective %>%
    lapply(function(adj){Subj_Theta[[adj]] %>% mutate(Adjective = adj)}) %>%
    do.call(rbind, .) %>%
    mutate(Subject_num = as.numeric(subject)) %>%
    left_join(Subj_num) %>%
    select(-c(Subject_num, subject))
  
  # Run SVD on the correlation matrix (used later to facilitate integration)
  # We only need to save the eigen values and the first row of the rotation matrix
  decompose <- function(s_mu,s_sig,rho){
    out = tibble(l1=numeric(),l2=numeric(),V11=numeric(),V12=numeric())
    for (i in 1:length(s_mu)){
      Decomposition <- eigen(matrix(c(s_mu[i]^2,rho[i]*s_mu[i]*s_sig[i],rho[i]*s_mu[i]*s_sig[i],s_sig[i]^2),nrow=2),symmetric=T)
      out[i,] = t(matrix(c(Decomposition$values,as.vector(t(Decomposition$vectors))))[1:4])
    }
    return(out)
  }
  Theta_param <- Theta_param %>%
    mutate((decompose)(s_mu,s_sig,rho)) %>%
    column_to_rownames("Adjective")
  
  # Save hyperparameters fitted to each adjective:
  write_csv(Theta_param %>% rownames_to_column(var = "Adjective"), "precomputed_data/Theta_params_gaussian.csv")
  
  # Save individual values fitted to each subject:
  tb_subj_Theta <- tb_subj_Theta %>%
    pivot_wider(names_from = variable, values_from = mean)
  write_csv(tb_subj_Theta, "precomputed_data/Theta_subjects_gaussian.csv")
  
}


### Visualize inferred threshold distribution ----

# Simulate the distribution to visualize second order vagueness

# Adjust number of samples (100 gives a good trade-off between representativeness and speed)
N_sim=100

simulated_theta_dist <- tibble(Adjective=character(),mu=numeric(),sigma=numeric())
for(i in 1:nrow(Theta_param)){
  samples <- with(Theta_param[i,],rmvnorm(n=N_sim,c(m_mu,m_sig),matrix(c(s_mu^2,rho*s_mu*s_sig,rho*s_mu*s_sig,s_sig^2),nrow=2)))
  simulated_theta_dist <- bind_rows(simulated_theta_dist,tibble(Adjective=row.names(Theta_param)[i],mu=samples[,1],sigma=exp(samples[,2])))
}

simulated_theta_dist <- simulated_theta_dist %>%
  mutate(
    Adjective = factor(Adjective,
                       levels=c("tall","powerful","hot","expensive","large","late","spicy","profitable","complete","young","safe","full"))
  )

pol_pnorm <- function(x,mean,sd,polarity){
  pnorm(polarity*x,mean,sd)
}


plot_data <- stat_grad_restr_data %>%
  left_join(select(degree_labels,-DegreeLabel), by = c("Adjective", "Probe")) %>%
  # mutate(NormDegree=if_else(Adjective %in% c("young","full","safe"),-NormDegree,NormDegree)) %>%
  mutate(Adjective = factor(Adjective,
                            levels=c("tall","powerful","hot","expensive","large","young","late","spicy","profitable","complete","safe","full")),
         DecodedLabel = gsub(" ","",DecodedLabel)
  ) %>% 
  group_by(Adjective) %>%
  mutate(ProbeRank = as.integer(factor(Probe))) %>%
  ungroup()


# We use normalized degrees as X axis, since we can't do the categorical trick we used for custom x-axis labels earlier.
# pdf("graphs/fitted_vague_denotations.pdf", width=7,height=6)
plot_data %>%
  ggplot(aes(y=AdjRating,x=NormDegree,group=Probe)) +
  facet_wrap(~Adjective,scales = "free_x",dir = "h") +
  geom_boxplot(outlier.shape = NA,width=.2,col=TwoColorPalette[1]) +
  mapply(function(mean,sd,adj) {
    range=1.15*range(plot_data$NormDegree[plot_data$Adjective==adj])
    stat_function(data=tibble(Adjective=adj,NormDegree=range,AdjRating=c(-1,2),Probe="I"),
                  fun = pol_pnorm,
                  args = list(mean = mean, sd = sd,polarity=if_else(adj%in%c("full","young","safe"),-1,1)),
                  colour = TwoColorPalette[2],
                  alpha=0.1,
                  n=25)}, 
    mean = simulated_theta_dist$mu, sd = simulated_theta_dist$sigma, adj=simulated_theta_dist$Adjective)+
  geom_boxplot(data=plot_data,outlier.shape = NA,width=.2,col=TwoColorPalette[1]) +
  scale_x_continuous(name="Normalized degree", limits = c(NA, NA))+
  scale_y_continuous(labels=scales::percent,name="") +
  theme_bw()+
  theme(axis.text.x = element_text(angle=90,hjust=1,vjust=0.5,size=6))
# dev.off()


## Step 2: Define and precompute Utilities ----

# joint PDF of (mu,sigma) with correlation (a bit faster than dmvnorm):
joint_pdf <- function(mu,sigma,m_mu,m_sig,s_mu,s_sig,rho){
  exp(
    -(((mu-m_mu)/s_mu)^2 - 2*rho*(mu-m_mu)*(log(sigma)-m_sig)/(s_mu*s_sig) + ((log(sigma)-m_sig)/s_sig)^2)/(2*(1-rho^2))
  )/(
    2*pi*s_mu*s_sig*sigma*sqrt(1-rho^2)
  )
}


# Define a tibble of extended data that includes the items that were not probed:
# (we need all 8 items on the scales, not just the 3 that were tested)
extended_data <- expand_grid(AssignmentId=unique(stat_grad_restr_data$AssignmentId),Probe=LETTERS[1:8]) %>%
  left_join(select(stat_grad_restr_data,-c(Degree,NormDegree,Order,WorkTimeInSeconds))) %>%
  mutate(Adjective = as.character(Adjective)) %>%
  group_by(AssignmentId) %>%
  mutate(
    Subject=max(Subject,na.rm=T),
    Adjective=max(Adjective,na.rm=T),
    Quantifier=max(Quantifier,na.rm=T),
    mean_pol_degree=max(mean_pol_degree,na.rm=T),
    sd_degree=max(sd_degree,na.rm=T),
    AdjMonotonicity=(max(AdjMonotonicity,na.rm=T)),
    AdjClass=(max(as.character(AdjClass),na.rm=T))
  ) %>%
  ungroup() %>%
  left_join(numeric_degrees) %>%
  mutate(
    pol_degree=if_else(AdjMonotonicity=="Positive", Degree,-Degree),
    pol_norm_degree = (pol_degree-mean_pol_degree)/sd_degree,
    NormDegree=if_else(AdjMonotonicity=="Positive", pol_norm_degree,-pol_norm_degree),
    ZeroDegree = -mean_pol_degree/sd_degree,
    Adjective = `contrasts<-`(factor(Adjective),11 , contr.sum(12))
  )


# Prior smoothing using a simple mixed model.
prior_smoothing <- function(degrees,priors,subject){
  if(sd(priors,na.rm=T)==0){return(rep(max(priors,na.rm=T),length(priors)))}
  X <- data.frame(prior = priors[!is.na(priors)], degree = degrees[!is.na(priors)], subject = factor(subject[!is.na(priors)]))
  smoothing_model <- lmer(qlogis(prior) ~ degree + (1+degree||subject), data = X)
  return(as.numeric(plogis(predict(smoothing_model,data.frame(degree = degrees,subject = subject)))))
}

# Apply the smoothing (we get a singularity with 'powerful')
extended_data <- extended_data %>%
  group_by(Adjective) %>%
  mutate(smoothed_prior = prior_smoothing(NormDegree,Prior, Subject)) %>%
  ungroup()

# Visualize the results: (clear correlations in some cases, not much in others)
# For 'powerful', the singularity is visible (no random intercepts so all curves cross at a single point)
extended_data %>%
  ggplot(aes(x = Degree, y = smoothed_prior, col = factor(AssignmentId), group = AssignmentId)) +
  facet_wrap(.~Adjective, scales = "free_x") +
  geom_line( alpha = .7) + 
  geom_point(data = na.omit(extended_data), aes(y = Prior), alpha = .8)+
  scale_color_discrete(guide = 'none')

# Inspect quality of the fit. Not great for 'young' and 'late'
extended_data %>%
  filter(!is.na(Prior)) %>%
  ggplot(aes(x = Prior, y = smoothed_prior)) +
  facet_wrap(.~Adjective) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, col = 'red') +
  theme_bw()

# Prepare a tibble with all the possible worlds given 8 individuals plus some useful things for the model computations:
# (highest element satisfying and not satisfying the predicate)
for (i in 1:8){
  assign(paste0("X",i),c(0,1))
}
world_template <- as_tibble(expand.grid(mget(paste0("X",1:8)))) %>%
  mutate(
    last_not_pred = max.col(across(starts_with("X"))==0, 'last'), 
    last_pred = max.col(across(starts_with("X"))==1, 'last'),
    last_not_pred = if_else(rowSums(across(starts_with("X")))==8,as.integer(NA),last_not_pred),
    last_pred = if_else(rowSums(across(starts_with("X")))==0,as.integer(NA),last_pred)
  )


## Step 3: Pre-compute bare utilities ----

# We need to pre-compute the utilities for each trial and each message in all 256 worlds.
# Best solution will be to store this in a 3d array
# At S1, we consider all messages
# At L1, we only keep the target and its relevant interpretations (4 if including absolute adjs, 2 if working on relative only)

# list of message-parse pairs:
messages = c("quant","quant_adj_lit","quant_adj_exh","quant_adj_min_lit","quant_adj_min_exh","list_WE","list_SE","null") # with absolute
messages = c("quant","quant_adj_lit","quant_adj_exh","list_WE","list_SE","null") # relative only

extended_data2 <- extended_data %>%
  filter(AdjClass == "Relative")

# We need to reassign a number to each subject and assignment so that there are no gaps after subsetting:
extended_data2 <- extended_data2 %>%
  mutate(AssignmentId2 = as.numeric(factor(AssignmentId)),
         # Convert probes to numeric and switch the order of probes for adjectives with negative polarity:
         probe_num = match(Probe,LETTERS[1:8]),
         probe_num = if_else(AdjMonotonicity=="Negative", 9L-probe_num,probe_num)
  ) %>%
  arrange(AssignmentId2,probe_num)

id_list <- unique(extended_data2$AssignmentId2)

get_utilities <- function(id){
  tmp <- extended_data2 %>%
    filter(AssignmentId2==id)
  quant_sum=ifelse(as.character(first(tmp$Quantifier))=="No",0,8)
  arr <- with(tmp,
              compute_utilities(
                pol_norm_degree,
                smoothed_prior,
                as.character(Adjective)[1],
                as.character(Quantifier)[1],
                strict_threshold=NA
              )) %>%
    mutate(
      U_null = log_prior,
      U_list_WE = WE_info,
      U_quant = if_else(rowSums(across(starts_with("X")))==quant_sum,0,-Inf),
      U_targ_lit = mean_logL0_lit,
      U_targ_exh = mean_logL0_exh,
      U_targ_strict_lit = logL0_strict_lit,
      U_targ_strict_exh = logL0_strict_exh
    ) %>%
    unite("world",(starts_with("X")),sep="") %>%
    select(world,starts_with("U_")) %>%
    column_to_rownames("world")
  return(do.call(cbind,arr))
}


# Saving precomputed values, even though it's not that long to recompute.
if (file.exists("precomputed_data/precomputed_utilities_rel.rds")) {
  all_utilities <- readRDS("precomputed_data/precomputed_utilities_rel.rds")
} else {
  # Takes about 1min (after fixing integrands)
  t <- Sys.time()
  all_utilities <- mclapply(id_list,get_utilities, mc.cores = 14L, mc.preschedule = FALSE)
  Sys.time()-t
  saveRDS(all_utilities,file="precomputed_data/precomputed_utilities_rel.rds")
}


# Sanity check: bare utilities should not excceed 0
max(sapply(all_utilities, max)) <= 0


## Define and fit RSA-SvI model ----
# (relative only)

# Get information by HIT to avoid redundancies:
pre_stan_data <- extended_data2 %>%
  group_by(AssignmentId2) %>%
  summarize(
    Adjective = as.character(first(Adjective)),
    Quantifier = as.character(first(Quantifier))
  ) %>%
  mutate(
    Adjective_num = as.numeric(factor(Adjective)),
    Quantifier_num = if_else(Quantifier=="Every",1,0)
  )

# Number of positive answers in each world (to compute cost of list answers)
n_list <- world_template %>% mutate(n_list=rowSums(across(starts_with("X")))) %>% pull(n_list)

# Quickly retrieve which worlds to sum to get the probability for a given probe:
worlds_index <- list(NULL)
for (i in 1:8){
  worlds_index[[i]] = which(as.numeric(unlist(world_template[,i]))==1)
}
worlds_indices <- do.call(rbind,worlds_index)

# Generate the data list for Stan
stan_model_data <- list(
  N = nrow(extended_data2),
  N_Id = n_distinct(extended_data2$AssignmentId2),
  N_adj = n_distinct(extended_data2$Adjective),
  N_measured = sum(!is.na(extended_data2$Posterior)),
  measured_indices = which(!is.na(extended_data2$Posterior)),
  Id = extended_data2$AssignmentId2, #[!is.na(extended_data2$Posterior)],
  adjective = pre_stan_data$Adjective_num,
  quantifier = pre_stan_data$Quantifier_num,
  all_utilities = lapply(all_utilities, function(x){x[, 1:5]}),
  n_list = n_list,
  y = extended_data2$Posterior[!is.na(extended_data2$Posterior)],
  degree = extended_data2$pol_norm_degree, #[!is.na(extended_data2$Posterior)],
  probe = extended_data2$probe_num, #[!is.na(extended_data2$Posterior)],
  worlds_indices = worlds_indices
)

# Define an initialization
RSA_SvI_init_function <- function()list(
  lambda = runif(1,2.5,4),
  cost_null=runif(1,0,10),
  cost_only=runif(1,0,2),
  cost_atom=runif(1,0.2,0.8),
  cost_every=runif(1,0,1),
  cost_no=runif(1,0,1),
  cost_adj=runif(stan_model_data[["N_adj"]],0,2),
  eps=runif(1,.1,.4)
)

RSA_SvI_model <- cmdstan_model("models/RSA_SvI_relative_model.stan")

suppressWarnings(dir.create("L1_samples"))

t <- Sys.time() # couple hours on my computer with only relatives
RSA_SvI_fit <- RSA_SvI_model$sample(
  data = stan_model_data,
  init = RSA_SvI_init_function,
  chains = 12, 
  parallel_chains = 12,
  iter_sampling = 1000,
  iter_warmup=1000,
  output_dir = "L1_samples"
)
Sys.time()-t

# Extract posterior distribution of relevant parameters:
RSA_SvI_posteriors <- RSA_SvI_fit$summary(variables = c("cost_null", "cost_atom", "cost_only", "cost_every", "cost_no", "cost_adj","lambda"))

# reverse mapping numerical adj to actual adjective
tb_adj_num <- pre_stan_data %>%
  select(Adjective, Adjective_num) %>%
  distinct()

# Make a LaTeX table
as_tibble(RSA_SvI_posteriors) %>%
  rename(parameter = variable) %>%
  mutate(Adjective_num = as.numeric(str_extract(parameter, "\\d"))) %>%
  left_join(tb_adj_num, by = "Adjective_num") %>%
  mutate(parameter = if_else(is.na(Adjective), parameter, paste0("cost_", Adjective))) %>%
  select(parameter, mean, sd, q5, q95) %>%
  xtable() %>%
  print(include.rownames = FALSE, booktabs = TRUE)

# Posterior on EXH parse:
pEXH_posterior <- RSA_SvI_fit$summary(variables = "pExh")
plot(sort(pEXH_posterior$mean))
plot(density(pEXH_posterior$mean))
summary(pEXH_posterior$mean)

loo_L1 <- RSA_SvI_fit$loo()
plot(loo_L1)


## Define and fit literal models for comparison ----

L0_posteriors <- function(
    degrees, # expected length 8
    priors, # expected length 8
    adj, # to retrieve Theta parameters, expected length 1
    quant, # expected length 1
    strict_threshold # expected length 1, used for the min-std interpretations
){
  # Make sure everything is properly sorted first:
  ord <- order(degrees)
  degrees=degrees[ord]
  priors=priors[ord]
  log_pred_prior = log(priors)
  clog_pred_prior = log(1-priors)
  worlds <- world_template %>%
    mutate(
      log_prior = (X1*log_pred_prior[1] + (1-X1)*clog_pred_prior[1]+
                     X2*log_pred_prior[2] + (1-X2)*clog_pred_prior[2]+
                     X3*log_pred_prior[3] + (1-X3)*clog_pred_prior[3]+
                     X4*log_pred_prior[4] + (1-X4)*clog_pred_prior[4]+
                     X5*log_pred_prior[5] + (1-X5)*clog_pred_prior[5]+
                     X6*log_pred_prior[6] + (1-X6)*clog_pred_prior[6]+
                     X7*log_pred_prior[7] + (1-X7)*clog_pred_prior[7]+
                     X8*log_pred_prior[8] + (1-X8)*clog_pred_prior[8]
      ),
      last_falsifier = if(quant=="No") last_pred else last_not_pred,
      last_false_degree = degrees[last_falsifier],
      last_false_degree = if_else(is.na(last_false_degree),-Inf,last_false_degree)
    )
  Th <- Theta_param[as.character(adj),]
  dmax <- max(degrees)
  # Use the thing computed above to compute mean_log_L0:
  worlds <- worlds %>%
    mutate(mean_L0(last_falsifier,log_prior,last_false_degree,degrees,Th))
  return(worlds)
}

get_L0_posteriors <- function(id){
  tmp <- extended_data2 %>%
    filter(AssignmentId2==id)
  quant_sum=ifelse(as.character(first(tmp$Quantifier))=="No",0,8)
  arr <- with(tmp,
              L0_posteriors(
                pol_norm_degree,
                smoothed_prior,
                as.character(Adjective)[1],
                as.character(Quantifier)[1],
                strict_threshold=NA
              )) %>%
    unite("world",(starts_with("X")),sep="") %>%
    select(world,log_prior, meanL0_lit, meanL0_exh) %>%
    column_to_rownames("world")
  return(do.call(cbind,arr))
}

# Precompute mean L0 values for each participant
if (file.exists("precomputed_data/precomputed_L0_post_rel.rds")) {
  all_L0_posteriors <- readRDS("precomputed_data/precomputed_L0_post_rel.rds")
} else {
  # Takes 3 min
  t <- Sys.time()
  all_L0_posteriors <- mclapply(id_list,get_L0_posteriors, mc.cores = 12L, mc.preschedule = FALSE)
  Sys.time()-t
  saveRDS(all_L0_posteriors,file="precomputed_data/precomputed_L0_post_rel.rds")
}

# Define the data list for Stan
stan_model_data_L0 <- list(
  N = nrow(extended_data2),
  N_Id = n_distinct(extended_data2$AssignmentId2),
  N_adj = n_distinct(extended_data2$Adjective),
  Id = extended_data2$AssignmentId2, #[!is.na(extended_data2$Posterior)],
  adjective = pre_stan_data$Adjective_num,
  quantifier = pre_stan_data$Quantifier_num,
  all_L0_posteriors = all_L0_posteriors,
  y = extended_data2$Posterior[!is.na(extended_data2$Posterior)],
  measured_indices = which(!is.na(extended_data2$Posterior)),
  N_measured = sum(!is.na(extended_data2$Posterior)),
  degree = extended_data2$pol_norm_degree, #[!is.na(extended_data2$Posterior)],
  probe = extended_data2$probe_num, #[!is.na(extended_data2$Posterior)],
  worlds_indices = worlds_indices
)


# Fit pure lit and pure exh L0 models:
L0_lit_model <- cmdstan_model("models/relative_L0lit_model.stan")

suppressWarnings(dir.create("L0_lit_samples"))

t <- Sys.time()
L0_lit_fit <- L0_lit_model$sample(
  data = stan_model_data_L0,
  chains = 12, 
  parallel_chains = 12,
  iter_sampling = 1000,
  iter_warmup=1000,
  output_dir = "L0_lit_samples"
)
Sys.time()-t

loo_L0_lit <- L0_lit_fit$loo(cores = 12L)
plot(loo_L0_lit)



L0_exh_model <- cmdstan_model("models/relative_L0exh_model.stan")

suppressWarnings(dir.create("L0_exh_samples"))

t <- Sys.time()
L0_exh_fit <- L0_exh_model$sample(
  data = stan_model_data_L0,
  chains = 12, 
  parallel_chains = 12,
  iter_sampling = 1000,
  iter_warmup = 1000,
  output_dir = "L0_exh_samples"
)
Sys.time()-t

loo_L0_exh <- L0_exh_fit$loo(cores = 12L)
plot(loo_L0_exh)


loo_compare(loo_L0_lit, loo_L0_exh)


## Compare the results ----


# Make a pretty LaTeX table:
loo::loo_compare(loo_L1, loo_L0_lit, loo_L0_exh)  %>%
  as_tibble() %>%
  mutate(Model = c("L0exh", "L0lit", "L1")) %>%
  select(Model, elpd_loo, elpd_diff, se_diff, p_loo) %>%
  xtable() %>%
  print(include.rownames = FALSE, booktabs = TRUE)

# Check the PSIS assumption: all k values should be below 0.7 ideally (we're far below here)
max(c(pareto_k_influence_values(loo_L1), pareto_k_influence_values(loo_L0_lit), pareto_k_influence_values(loo_L0_exh)))


# Make pretty graph ----

# Lengthy code, but worth it!

# Extract predictions:
stan_L1_predictions <- RSA_SvI_fit$summary(variables = "prediction") %>% as_tibble() %>% select(mean,sd)
stan_L0lit_predictions <- L0_lit_fit$summary(variables = "prediction") %>% as_tibble() %>% select(mean,sd)
stan_L0exh_predictions <- L0_exh_fit$summary(variables = "prediction") %>% as_tibble() %>% select(mean,sd)

extended_data2 <- extended_data2 %>% mutate(
  prediction_L1 = stan_L1_predictions$mean, 
  prediction_L0lit = stan_L0lit_predictions$mean, 
  prediction_L0exh = stan_L0exh_predictions$mean)


extended_data2 %>%
  filter(!is.na(Posterior)) %>%
  pivot_longer(c(prediction_L0exh, prediction_L1), values_to = "Prediction", names_to = "Model") %>%
  ggplot(aes(x = Posterior, y = Prediction, col = Model)) + 
  facet_grid(Quantifier ~ Adjective) + 
  geom_point(alpha = .7) +
  geom_smooth(method = "lm") +
  geom_abline(slope = 1, intercept = 0) +
  scale_color_manual(values = TwoColorPalette) +
  theme_bw()

plot_data2 <- extended_data2 %>%
  # filter(!is.na(Posterior)) %>%
  select(AssignmentId, Probe, Quantifier, Adjective, Degree, smoothed_prior, Posterior, prediction_L1, prediction_L0exh, prediction_L0lit) %>%
  pivot_longer(c(smoothed_prior, Posterior, prediction_L1, prediction_L0exh, prediction_L0lit), names_to = "Model", values_to = "Prediction") %>%
  mutate(Model = case_when(
    Model == "prediction_L0exh" ~ "L0exh",
    Model == "prediction_L0lit" ~ "L0lit",
    Model == "prediction_L1" ~ "L1",
    Model == "Posterior" ~ "measured",
    Model == "smoothed_prior" ~ "prior"
  )) %>%
  left_join(select(degree_labels,-DegreeLabel)) %>%
  mutate(Adjective = factor(Adjective,
                            levels=c("tall","powerful","hot","expensive","large","young")),
         DecodedLabel = gsub(" ","",DecodedLabel)
  ) %>% 
  group_by(Adjective) %>%
  mutate(ProbeRank = as.integer(factor(Probe))) %>%
  ungroup()

ribbon_data1 <- plot_data2 %>% 
  filter(!Model %in% c("measured", "prior")) %>%
  group_by(Adjective, Quantifier, Model, Degree, ProbeRank) %>%
  summarise(y = median(Prediction), ymin = quantile(Prediction, 0.25), ymax = quantile(Prediction, 0.75))
ribbon_data2 <- plot_data2 %>% 
  filter(Model == "prior") %>%
  group_by(Adjective, Model, Degree, ProbeRank) %>%
  summarise(y = median(Prediction), ymin = quantile(Prediction, 0.25), ymax = quantile(Prediction, 0.75))


breaks_fun <- function(x) {
  count_breaks <<- count_breaks + 0.5
  switch(
    floor(count_breaks),
    sort(unique(plot_data2$ProbeRank[plot_data2$Adjective=="tall"])),
    sort(unique(plot_data2$ProbeRank[plot_data2$Adjective=="powerful"])),
    sort(unique(plot_data2$ProbeRank[plot_data2$Adjective=="hot"])),
    sort(unique(plot_data2$ProbeRank[plot_data2$Adjective=="expensive"])),
    sort(unique(plot_data2$ProbeRank[plot_data2$Adjective=="large"])),
    sort(unique(plot_data2$ProbeRank[plot_data2$Adjective=="young"]))
  )
}

labels_fun <- function(x) {
  count_label <<- count_label + 1L
  switch(
    count_label,
    unique(plot_data2$DecodedLabel[plot_data2$Adjective=="tall"][order(plot_data2$ProbeRank[plot_data2$Adjective=="tall"])]),
    unique(plot_data2$DecodedLabel[plot_data2$Adjective=="powerful"][order(plot_data2$ProbeRank[plot_data2$Adjective=="powerful"])]),
    unique(plot_data2$DecodedLabel[plot_data2$Adjective=="hot"][order(plot_data2$ProbeRank[plot_data2$Adjective=="hot"])]),
    unique(plot_data2$DecodedLabel[plot_data2$Adjective=="expensive"][order(plot_data2$ProbeRank[plot_data2$Adjective=="expensive"])]),
    unique(plot_data2$DecodedLabel[plot_data2$Adjective=="large"][order(plot_data2$ProbeRank[plot_data2$Adjective=="large"])]),
    unique(plot_data2$DecodedLabel[plot_data2$Adjective=="young"][order(plot_data2$ProbeRank[plot_data2$Adjective=="young"])])
  )
}

p <- plot_data2 %>%
  ggplot(aes(x=ProbeRank,y=Prediction, col = Model, fill = Model, linetype = (Model == "prior"))) +
  facet_grid(Quantifier~Adjective, scales = 'free_x')+
  coord_cartesian(ylim = c(0,1)) +
  scale_color_manual(values = FiveColorPalette[c(1, 2, 5, 4, 4)]) +
  scale_fill_manual(values = c(FiveColorPalette[c(1, 2, 5, 4)], "transparent", FiveColorPalette[4])) +
  scale_linetype_manual(values = c(1, 2), guide = "none") +
  scale_y_continuous(labels = scales::percent) +
  theme_bw()


for (adj in unique(extended_data2$Adjective)) {
  p <- p + geom_boxplot(
    data = plot_data2 %>% filter(Model == "measured" & Adjective == adj & !is.na(Prediction)), 
    aes(group = Degree), 
    alpha = .3, width = 1)
}

p <- p+geom_ribbon(data = ribbon_data1, aes(y = y, ymin = ymin, ymax = ymax, colour = after_scale(alpha(colour, 0.7))), alpha = .3)+
  geom_line(data = ribbon_data1, aes(y = y))+
  geom_line(data = ribbon_data2, aes(y = y))


# Don't forget to reset counters each time you run this
count_breaks <- 1
count_label <- 0
p <- p +
  scale_x_continuous(name="Measure",breaks = breaks_fun, minor_breaks = NULL,labels=labels_fun, limits = c(NA, NA))+
  theme(axis.text.x = element_text(angle=60,hjust=1,vjust=1,size=6))+
  theme(text = element_text(family = "Arial Unicode MS"))

# quartz(width = 16, height = 7)
print(p)
# ps = gridSVG::grid.export("graphs/prediction_graph.svg", addClass=T,progress = T, exportJS="inline")


