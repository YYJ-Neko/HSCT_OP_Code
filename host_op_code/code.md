# #######Chi-square & Competing Risk Model
library(cmprsk)
library(ggplot2)
library(readxl)

df1 <- read.csv("data.csv")

df <- df1[,c(2,15,16)]
df <- na.omit(df)
time <- df[,3]
event <- df[,2]
group <- df[,1]

cr_model <- cuminc(time,event,group)

plot_data <- data.frame()
for(i in seq_along(cr_model)){
  temp <- cr_model[[i]]
  if(is.list(temp)&&all(c("time","est")%in%names(temp))){
    temp_data <- data.frame(time=temp$time,cumulative_incidence=temp$est,group=names(cr_model)[i])
    plot_data <- rbind(plot_data,temp_data)
  }
}

ggplot(plot_data,aes(x=time,y=cumulative_incidence,color=group))+
geom_line()+
labs(x="Time",y="Cumulative Incidence",title="Cumulative Incidence Curve",color="Group")+
theme_minimal()

df <- df1[,c(1,23,24,9,50)]
time <- df[,3]
event <- df[,2]
group <- df[,1]
cov1 <- df[,4]
cov2 <- df[,5]
df <- data.frame(time,event,group,cov1,cov2)
df$cov1 <- as.factor(df$cov1)
df$cov2 <- as.factor(df$cov2)
df_clean <- na.omit(df)

X <- model.matrix(~group+cov1+cov2,data=df_clean)[,-1]
df_clean$fstatus <- as.numeric(factor(df_clean$event))

fit_multi <- crr(ftime=df_clean$time,fstatus=df_clean$fstatus,cov1=X)
summary(fit_multi)
# #######Alpha Diversity
library(vegan)
library(picante)
library(ggpubr)
library(ggplot2)
library(dplyr)

otu <- read.csv("otu.csv", row.names = 1)
group <- read.csv("group.csv", row.names = 1)

Shannon <- diversity(otu, index = "shannon", MARGIN = 2)
Simpson <- diversity(otu, index = "simpson", MARGIN = 2)
Richness <- specnumber(otu, MARGIN = 2)
obs_chao_ace <- t(estimateR(t(otu)))

index <- data.frame(
  Shannon = Shannon,
  Simpson = Simpson,
  Richness = Richness,
  Chao = obs_chao_ace[,2],
  ACE = obs_chao_ace[,4],
  Sobs = obs_chao_ace[,1]
)
index$Pielou <- Shannon/log(Richness,2)
index$Goods_coverage <- 1-colSums(otu==1)/colSums(otu)
index$sample <- rownames(index)

df <- merge(index, group, by = "sample")
write.csv(df, "alpha_diversity.csv", row.names = F)
# #######Propensity Score Matching
library(tidyverse)
library(MatchIt)

df <- read.csv("psm_raw.csv", fileEncoding = "GBK")
otu <- read.csv("otu.csv", fileEncoding = "GBK")

catVars <- c("V1","V2","V3","V4","V5","V6","V7","V8","V9","V10")
df[catVars] <- lapply(df[catVars], as.factor)
df$group <- as.factor(df$group)

m.out <- matchit(group ~ V1 + V2 + V3 + V4 + V5 + V6 + V7,
                 data = df,
                 method = "nearest",
                 ratio = 1)
summary(m.out)
plot(m.out, type = "jitter")

matched_data <- match.data(m.out)
write.csv(matched_data, "psm_matched.csv", row.names = FALSE)

rownames(matched_data) <- matched_data[[1]]
rownames(otu) <- otu[[1]]
otu_matched <- otu[rownames(matched_data),]
write.csv(otu_matched, "otu_psm.csv")
# #######Microbial Network Analysis
library(igraph)
library(foreach)
library(doParallel)
library(abind)

network_speed_up <- function(otu_table,threads){
  otu_table2 <- apply(otu_table,2,rank)
  r <- function(rx,ry){
    n <- length(rx)
    lxy <- sum((rx-mean(rx))*(ry-mean(ry)))
    lxx <- sum((rx-mean(rx))^2)
    lyy <- sum((ry-mean(ry))^2)
    r <- lxy/sqrt(lxx*lyy)
    t <- (r * sqrt(n - 2))/sqrt(1 - r^2)
    p <- -2 * expm1(pt(abs(t), (n - 2), log.p = TRUE))
    return(c(r,p))
  }
  arraybind <- function(...) abind(...,along = 3,force.array=TRUE)
  nc <- ncol(otu_table)
  registerDoParallel(cores = threads)
  corr <- foreach(i = 1:nc,.combine = "arraybind") %dopar% {
    corr1 <- matrix(0,2,nc)
    for(j in 1:nc) if(j>i) corr1[,j] <- r(otu_table2[,i],otu_table2[,j])
    corr1
  }
  rr <- corr[1,,]; rr <- rr+t(rr); diag(rr) <- 1
  pp <- corr[2,,]; lp <- lower.tri(pp); pp[lp] <- p.adjust(pp[lp],"fdr"); pp <- pp+t(pp)
  dimnames(rr) <- dimnames(pp) <- list(colnames(otu_table),colnames(otu_table))
  return(list(r=rr,p=pp))
}

otu <- read.csv("otu.csv",row.names=1)
otu1 <- t(otu)
abundance_threshold <- 300000
sample_threshold <- ceiling(0.01*nrow(otu1))
filtered_otu <- otu1; filtered_otu[filtered_otu<abundance_threshold] <- 0; filtered_otu[filtered_otu>=abundance_threshold] <- 1
otu2 <- otu1[,colSums(filtered_otu)>=sample_threshold]
otu <- t(otu2)

res <- network_speed_up(otu,detectCores())
s <- res$r * ifelse(res$p<=0.05,1,0)
s[is.na(s)] <- 0

m <- graph.adjacency(s,weighted=T,mode="undirected",diag=F)
m <- simplify(m)
m <- delete.vertices(m,names(degree(m)[degree(m)==0]))
E(m)$correlation <- E(m)$weight
E(m)$weight <- abs(E(m)$weight)

tax <- read.csv("tax.csv",row.names=1)
V(m)$tax <- tax[V(m)$name,1]

edge <- data.frame(as_edgelist(m))
edge <- data.frame(source=edge[[1]],target=edge[[2]],weight=E(m)$weight,cor=E(m)$correlation)
write.csv(edge,"network.edge.csv",row.names=F)
node <- data.frame(id=V(m)$name,tax=V(m)$tax)
write.csv(node,"network.node.csv",row.names=F)

clu <- spinglass.community(m, spins=200, update.rule="simple", implementation="neg")
V(m)$module <- membership(clu)
write.csv(data.frame(id=V(m)$name,module=V(m)$module),"module.csv",row.names=F)

calcTC <- function(g) sum(abs(E(g)$weight))
initTC <- calcTC(m)
keyTaxa <- data.frame(id=V(m)$name,deltaTC=NA)
for(n in V(m)$name){
  gtemp <- delete_vertices(m,n)
  keyTaxa[keyTaxa$id==n,"deltaTC"] <- initTC - calcTC(gtemp)
}
keyTaxa <- keyTaxa[order(-keyTaxa$deltaTC),]
# #######Survival Analysis
library(survival)
library(survminer)

dat <- read.csv("surv_raw.csv", fileEncoding = "GBK")
fit <- survfit(Surv(time, status) ~ Group, data = dat)
ggsurvplot(fit, data = dat, pval = TRUE) |> ggsave("km.pdf", width = 4, height = 6)

cox_dat <- read.csv("cox_raw.csv", fileEncoding = "GBK")
univ <- do.call(rbind, lapply(colnames(cox_dat)[-c(1,2)], function(x){
  s <- summary(coxph(as.formula(sprintf("Surv(OS, Status) ~ %s", x)), data = cox_dat))
  data.frame(var=x, HR=s$coefficients[2], L95=s$conf.int[3], H95=s$conf.int[4], P=s$coefficients[5])
}))
write.csv(univ, "univ_cox.csv", row.names=F)

mul <- summary(coxph(Surv(OS, Status) ~ ., data = cox_dat))
mul_res <- data.frame(HR=mul$coefficients[,2], L95=mul$conf.int[,3], H95=mul$conf.int[,4], P=mul$coefficients[,5])
write.csv(mul_res, "mult_cox.csv", row.names=F)
# #######Lefse
library(microeco)

otu <- read.csv("otu.csv",row.names=1)
sam <- read.csv("sample.csv",row.names=1,fileEncoding="GBK")
tax <- read.csv("tax.csv",row.names=1)
sam$Group <- factor(sam$Group)

d <- microtable$new(otu_table=otu,sample_table=sam,tax_table=tax)
lf <- trans_diff$new(d,method="lefse",group="Group",taxa_level="Species",alpha=0.05,p_adjust_method="non")
write.csv(lf$res_diff,"lefse_res.csv",row.names=F)

p <- lf$plot_diff_bar(color_values=c("#669aba","#be1420"),group_order=c("G1","G2")) +
  labs(y = "LDA Score (= 0 )")
ggsave("lefse.pdf",p,width=8,height=6)
# #######Mantel's Test
library(linkET)
library(FD)

env <- read.csv("env.csv",row.names=1)
spe <- read.csv("spe.csv",row.names=1)
spe <- spe[rownames(env),]
spe[] <- lapply(spe, as.numeric)

mantel <- mantel_test(spe, env, mantel_fun = 'mantel', spec_dist = 'bray', env_dist = 'euclidean',
                      spec_select = list(M1 = 1:12,M2 = 13:46,M3 = 47:82)) %>%
  mutate(rd = cut(r, breaks = c(-Inf,0.2,0.4,Inf), labels = c("< 0.2","0.2 - 0.4",">= 0.4")),
         pd = cut(p, breaks = c(-Inf,0.01,0.05,Inf), labels = c("< 0.01","0.01 - 0.05",">= 0.05")))
write.csv(mantel, "mantel_result.csv", row.names=F)
# #######ReportScore
library(ReporterScore)
library(dplyr)

otu <- read.csv("ko_otu.csv",row.names=1)
group <- read.csv("ko_group.csv",row.names=1)

ko_p <- ko.test(otu,"ShannonGroup",group)
ko_z <- pvalue2zs(ko_p,mode="directed")
rs_res <- get_reporter_score(ko_z)
res_all <- reporter_score(otu,"ShannonGroup",group,mode="directed")

write.csv(rs_res,"reporter_score.csv",row.names=F)
# #######Univariate Analysis
library(readxl,dplyr,broom,forestplot,grid)
df <- read_excel("dat.xlsx") %>% select(-V1,-V2) %>% mutate(Y=as.factor(Y))
cv <- setdiff(names(df),c("Y","V3","V4","V5","V6"))
df <- df %>% mutate(across(all_of(cv),as.factor))
vars <- setdiff(names(df),"Y")
uni_res <- lapply(vars,function(v){
  m <- try(glm(as.formula(paste("Y~",v)),data=df,family=binomial),silent=T)
  if(inherits(m,"try-error")) return(NULL)
  tidy(m,exponentiate=T,conf.int=T) %>% filter(term!="(Intercept)") %>% mutate(variable=v)
}) %>% bind_rows()
write.csv(uni_res,"uni_all.csv",row.names=F)
uni_res_p <- uni_res %>% filter(p.value<0.05) %>% mutate(across(where(is.numeric),round,3))
write.csv(uni_res_p,"uni_sig.csv",row.names=F)
data <- read.csv("uni_sig.csv",stringsAsFactors=F)
flip_or_row <- function(d,pat,lab=NULL){
  idx <- grep(pat,d$Variable,ignore.case=T)
  if(!length(idx)) return(d)
  for(i in idx){
    d$OR[i] <- 1/d$OR[i];d$Lower[i] <- 1/d$Upper[i];d$Upper[i] <- 1/d$Lower[i]
    if(!is.null(lab)) d$Variable[i] <- lab
  }
  d
}
data <- flip_or_row(data,"V3","V3")
data <- flip_or_row(data,"V4","V4 abnormal")
data$Significance <- ifelse(data$Pvalue<0.001,"***",ifelse(data$Pvalue<0.01,"**",ifelse(data$Pvalue<0.05,"*","")))
tb1 <- rbind(c("Variable","OR","95% CI","P Value"),cbind(data$Variable,sprintf("%.3f",data$OR),paste0("(",sprintf("%.3f",data$Lower),"-",sprintf("%.3f",data$Upper),")"),data$Significance))
forestplot(labeltext=tb1,mean=c(NA,data$OR),lower=c(NA,data$Lower),upper=c(NA,data$Upper),is.summary=c(T,rep(F,nrow(data))),xlog=F,col=fpColors(box="#2c7bb6",lines="black"),txt_gp=fpTxtGp(label=gpar(fontsize=12),ticks=gpar(fontsize=11),title=gpar(fontsize=14,fontface="bold")),lineheight=unit(1.2,"cm"),zero=1,xticks=seq(0,6,1),clip=c(0.5,3),boxsize=0.15,graph.pos=4,line.margin=0.2,align="l")
grid.text("Univariate Logistic Regression Analysis",x=0.5,y=0.75,gp=gpar(fontsize=15,fontface="bold"))
data$Pvalue_formatted <- ifelse(data$Pvalue<0.001,"<0.001",sprintf("%.3f",data$Pvalue))
tb2 <- rbind(c("Variable","OR","95% CI","P Value"),cbind(data$Variable,sprintf("%.3f",data$OR),paste0("(",sprintf("%.3f",data$Lower),"-",sprintf("%.3f",data$Upper),")"),data$Pvalue_formatted))
forestplot(labeltext=tb2,mean=c(NA,data$OR),lower=c(NA,data$Lower),upper=c(NA,data$Upper),is.summary=c(T,rep(F,nrow(data))),xlog=F,col=fpColors(box="#2c7bb6",lines="black"),txt_gp=fpTxtGp(label=gpar(fontsize=12),ticks=gpar(fontsize=11),title=gpar(fontsize=14,fontface="bold")),lineheight=unit(1.2,"cm"),zero=1,xticks=seq(0,6,1),clip=c(0.5,3),boxsize=0.15,graph.pos=4,line.margin=0.2,align="l")
grid.text("Univariate Logistic Regression Analysis",x=0.5,y=0.75,gp=gpar(fontsize=15,fontface="bold"))
# #######Machine Learning
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler, PolynomialFeatures
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score, roc_curve, accuracy_score, precision_score, recall_score, f1_score
from tqdm import tqdm
from scipy.stats import loguniform
import time
import shap
import matplotlib as mpl
import matplotlib.colors as mcolors

mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
mpl.rcParams['svg.fonttype'] = 'none'
mpl.rcParams['text.usetex'] = False
mpl.rcParams['font.family'] = 'sans-serif'
mpl.rcParams['font.sans-serif'] = ['Arial', 'DejaVu Sans']
mpl.rcParams['axes.unicode_minus'] = False

GRADIENT_START_COLOR = "#bdbdbd"
GRADIENT_END_COLOR = "#eb4a40"
GRADIENT_MID_COLOR = None
GRADIENT_TOP_USES_END = True
DATA_PATH = "input.csv"
TEST_DATA_PATH = "input_test.csv"
OUTPUT_DIR = "output"
RANDOM_STATE = 42
TOP_N_FEATURES = 10
N_ESTIMATORS = 5

def _make_cmap(start_hex, end_hex, mid_hex=None):
    if mid_hex is None:
        colors = [start_hex, end_hex]
    else:
        colors = [start_hex, mid_hex, end_hex]
    return mcolors.LinearSegmentedColormap.from_list("cmap", colors, N=256)

def _ranks_to_colors(n, cmap, top_uses_end=True):
    if n <= 1:
        return [cmap(1.0 if top_uses_end else 0.0)]
    pos = np.linspace(0.0, 1.0, n)
    if top_uses_end:
        pos = pos[::-1]
    return [cmap(p) for p in pos]

def create_output_dir(directory):
    os.makedirs(directory, exist_ok=True)

def timer(func):
    def wrapper(*args, **kwargs):
        s = time.time()
        res = func(*args, **kwargs)
        e = time.time()
        return res
    return wrapper

def plot_roc_curves(fpr_list, tpr_list, mean_fpr, mean_tpr, auc_scores, output_dir):
    plt.figure(figsize=(10,8))
    for i,(fpr,tpr) in enumerate(zip(fpr_list,tpr_list)):
        plt.plot(fpr,tpr,lw=1,alpha=0.3)
    plt.plot(mean_fpr,mean_tpr,color='b',lw=2,alpha=0.8)
    plt.plot([0,1],[0,1],linestyle='--',lw=2,color='r',alpha=0.8)
    plt.xlim([-0.05,1.05])
    plt.ylim([-0.05,1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.tight_layout()
    plt.savefig(f"{output_dir}/roc_cv.png", dpi=300)
    plt.close()

def plot_test_roc_curve(fpr, tpr, auc_score, output_dir):
    plt.figure(figsize=(10,8))
    plt.plot(fpr,tpr,color='b',lw=2,alpha=0.8)
    plt.plot([0,1],[0,1],linestyle='--',lw=2,color='r',alpha=0.8)
    plt.xlim([-0.05,1.05])
    plt.ylim([-0.05,1.05])
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate')
    plt.tight_layout()
    plt.savefig(f"{output_dir}/roc_test.png", dpi=300)
    plt.close()

def plot_feature_importance(names, scores, output_dir, top_n):
    plt.figure(figsize=(10,8))
    cmap = _make_cmap(GRADIENT_START_COLOR, GRADIENT_END_COLOR, GRADIENT_MID_COLOR)
    colors = _ranks_to_colors(len(names), cmap, GRADIENT_TOP_USES_END)
    plt.barh(range(len(names)), scores, color=colors)
    plt.yticks(range(len(names)), names)
    plt.xlabel('Value')
    plt.gca().invert_yaxis()
    plt.grid(axis='x', alpha=0.3)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/feature_importance_{top_n}.png", dpi=300)
    plt.close()

@timer
def load_data(path):
    df = pd.read_csv(path)
    target = df.columns[-1]
    X = df.drop(target, axis=1).iloc[:,1:]
    y = df[target].values
    numerical_cols = X.select_dtypes(exclude=['object','category']).columns
    X_num = X[numerical_cols].copy()
    feature_names = numerical_cols.tolist()
    return X_num, y, feature_names

@timer
def select_features(X, y, feature_names):
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    poly = PolynomialFeatures(degree=2, include_bias=False, interaction_only=True)
    X_poly = poly.fit_transform(X_scaled)
    poly_names = poly.get_feature_names_out(feature_names).tolist()
    l1 = LogisticRegression(C=0.1, penalty='l1', solver='liblinear', max_iter=2000, class_weight='balanced', random_state=RANDOM_STATE)
    l1.fit(X_poly, y)
    imp = np.abs(l1.coef_[0])
    sorted_pairs = sorted(zip(poly_names, imp), key=lambda x:x[1], reverse=True)
    top_features = [x[0] for x in sorted_pairs[:TOP_N_FEATURES]]
    top_imp = [x[1] for x in sorted_pairs[:TOP_N_FEATURES]]
    plot_feature_importance(top_features, top_imp, OUTPUT_DIR, TOP_N_FEATURES)
    pd.DataFrame({'Rank':range(1,TOP_N_FEATURES+1),'Feature':top_features,'Importance':top_imp}).to_csv(f"{OUTPUT_DIR}/feature_selection.csv", index=False)
    with open(f"{OUTPUT_DIR}/top_features.txt","w") as f:
        for ft in top_features:
            f.write(ft+"\n")
    return top_features, poly

def preprocess_fold(X_train, X_test, orig_names, fixed_features, poly):
    scaler = StandardScaler()
    Xtr_scaled = scaler.fit_transform(X_train)
    Xte_scaled = scaler.transform(X_test)
    Xtr_poly = poly.transform(Xtr_scaled)
    Xte_poly = poly.transform(Xte_scaled)
    poly_names = poly.get_feature_names_out(orig_names)
    idx = [i for i,n in enumerate(poly_names) if n in fixed_features]
    Xtr_sel = Xtr_poly[:,idx]
    Xte_sel = Xte_poly[:,idx]
    sel_names = [n for n in fixed_features if n in poly_names]
    return Xtr_sel, Xte_sel, sel_names, scaler

@timer
def train_ensemble(X, y, orig_names, fixed_features, poly):
    param_sets = [{'C':loguniform(1e-4,1e-2).rvs(random_state=i),'penalty':'l1'} for i in range(N_ESTIMATORS//2)] + [{'C':loguniform(1e-2,1e2).rvs(random_state=i+N_ESTIMATORS//2),'penalty':'l2'} for i in range(N_ESTIMATORS - N_ESTIMATORS//2)]
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)
    X_poly = poly.transform(X_scaled)
    poly_names = poly.get_feature_names_out(orig_names)
    idx = [i for i,n in enumerate(poly_names) if n in fixed_features]
    X_sel = X_poly[:,idx]
    models = []
    for i,p in enumerate(param_sets):
        m = LogisticRegression(C=p['C'], penalty=p['penalty'], solver='liblinear', class_weight='balanced', max_iter=2000, random_state=RANDOM_STATE+i)
        m.fit(X_sel,y)
        models.append(m)
    return models, scaler, [n for n in fixed_features if n in poly_names], param_sets

@timer
def evaluate_cv(X, y, orig_names, fixed_features, poly, out_dir):
    kf = StratifiedKFold(n_splits=5, shuffle=True, random_state=RANDOM_STATE)
    fold_aucs = []
    fpr_list, tpr_list = [], []
    param_sets = [{'C':loguniform(1e-4,1e-2).rvs(random_state=i),'penalty':'l1'} for i in range(N_ESTIMATORS//2)] + [{'C':loguniform(1e-2,1e2).rvs(random_state=i+N_ESTIMATORS//2),'penalty':'l2'} for i in range(N_ESTIMATORS - N_ESTIMATORS//2)]
    for train_idx, test_idx in tqdm(kf.split(X,y), total=5):
        Xtr_raw, Xte_raw = X.iloc[train_idx], X.iloc[test_idx]
        ytr, yte = y[train_idx], y[test_idx]
        Xtr_sel, Xte_sel, _, _ = preprocess_fold(Xtr_raw, Xte_raw, orig_names, fixed_features, poly)
        fold_models = []
        for i,p in enumerate(param_sets):
            m = LogisticRegression(C=p['C'], penalty=p['penalty'], solver='liblinear', class_weight='balanced', max_iter=2000, random_state=RANDOM_STATE+i)
            m.fit(Xtr_sel,ytr)
            fold_models.append(m)
        probs = np.mean([m.predict_proba(Xte_sel)[:,1] for m in fold_models], axis=0)
        auc = roc_auc_score(yte, probs)
        fold_aucs.append(auc)
        fpr,tpr,_ = roc_curve(yte, probs)
        fpr_list.append(fpr)
        tpr_list.append(tpr)
    mean_fpr = np.linspace(0,1,100)
    mean_tpr = np.zeros_like(mean_fpr)
    for fpr,tpr in zip(fpr_list,tpr_list):
        mean_tpr += np.interp(mean_fpr, fpr, tpr)
    mean_tpr /= len(fpr_list)
    plot_roc_curves(fpr_list, tpr_list, mean_fpr, mean_tpr, fold_aucs, out_dir)
    final_ensemble, train_scaler, final_ft, final_params = train_ensemble(X,y,orig_names,fixed_features,poly)
    scaler_shap = StandardScaler().fit(X)
    X_scaled_shap = scaler_shap.transform(X)
    X_poly_shap = poly.transform(X_scaled_shap)
    poly_names_shap = poly.get_feature_names_out(orig_names)
    idx_shap = [i for i,n in enumerate(poly_names_shap) if n in fixed_features]
    X_shap_sel = X_poly_shap[:,idx_shap]
    ref_model = LogisticRegression(C=0.1, penalty='l1', solver='liblinear', max_iter=2000, class_weight='balanced', random_state=RANDOM_STATE)
    ref_model.fit(X_shap_sel,y)
    return {'auc':np.mean(fold_aucs), 'fold_aucs':fold_aucs, 'ref_model':ref_model, 'ensemble':final_ensemble, 'params':final_params, 'features':final_ft, 'X_shap':X_shap_sel, 'scaler':train_scaler, 'poly':poly}

@timer
def evaluate_test(X_test, y_test, orig_names, fixed_features, poly, train_scaler, ensemble, params, out_dir):
    Xte_scaled = train_scaler.transform(X_test)
    Xte_poly = poly.transform(Xte_scaled)
    poly_names = poly.get_feature_names_out(orig_names)
    idx = [i for i,n in enumerate(poly_names) if n in fixed_features]
    Xte_sel = Xte_poly[:,idx]
    probs = np.mean([m.predict_proba(Xte_sel)[:,1] for m in ensemble], axis=0)
    pred = (probs >= 0.5).astype(int)
    auc = roc_auc_score(y_test, probs)
    acc = accuracy_score(y_test, pred)
    prec = precision_score(y_test, pred, zero_division=0)
    rec = recall_score(y_test, pred, zero_division=0)
    f1 = f1_score(y_test, pred, zero_division=0)
    fpr,tpr,_ = roc_curve(y_test, probs)
    plot_test_roc_curve(fpr,tpr,auc,out_dir)
    pd.DataFrame({'Metric':['AUC','Accuracy','Precision','Recall','F1'],'Value':[auc,acc,prec,rec,f1]}).to_csv(f"{out_dir}/test_metrics.csv", index=False)
    pred_df = pd.DataFrame({'y_true':y_test,'prob':probs,'pred':pred})
    for i,p in enumerate([m.predict_proba(Xte_sel)[:,1] for m in ensemble]):
        pred_df[f'model_{i}_prob'] = p
    pred_df.to_csv(f"{out_dir}/test_predictions.csv", index=False)
    pd.DataFrame(params).to_csv(f"{out_dir}/ensemble_params.csv", index=False)
    return {'auc':auc,'acc':acc,'prec':prec,'rec':rec,'f1':f1}

@timer
def shap_analysis(model, X, feature_names, out_dir):
    masker = shap.maskers.Independent(data=X)
    explainer = shap.LinearExplainer(model, masker=masker)
    shap_vals = explainer.shap_values(X)
    plt.figure(figsize=(12,8))
    shap.summary_plot(shap_vals, X, feature_names=feature_names, show=False)
    plt.tight_layout()
    plt.savefig(f"{out_dir}/shap_beeswarm.png", dpi=300)
    plt.close()
    shap_imp = np.mean(np.abs(shap_vals), axis=0)
    sorted_idx = np.argsort(shap_imp)[::-1][:TOP_N_FEATURES]
    sorted_ft = [feature_names[i] for i in sorted_idx]
    sorted_imp = shap_imp[sorted_idx]
    plot_feature_importance(sorted_ft, sorted_imp, out_dir, 'shap')
    pd.DataFrame({'Feature':feature_names,'SHAP':shap_imp}).sort_values('SHAP',ascending=False).to_csv(f"{out_dir}/shap_importance.csv", index=False)
    return shap_vals

def main():
    create_output_dir(OUTPUT_DIR)
    X, y, orig_ft = load_data(DATA_PATH)
    fixed_ft, poly = select_features(X,y,orig_ft)
    res = evaluate_cv(X,y,orig_ft,fixed_ft,poly,OUTPUT_DIR)
    pd.DataFrame({'MeanAUC':[res['auc']],'StdAUC':[np.std(res['fold_aucs'])]}).to_csv(f"{OUTPUT_DIR}/cv_summary.csv", index=False)
    shap_analysis(res['ref_model'], res['X_shap'], res['features'], OUTPUT_DIR)
    if os.path.exists(TEST_DATA_PATH):
        Xte, yte, test_ft = load_data(TEST_DATA_PATH)
        if test_ft == orig_ft:
            test_res = evaluate_test(Xte,yte,orig_ft,fixed_ft,poly,res['scaler'],res['ensemble'],res['params'],OUTPUT_DIR)
            pd.DataFrame({'Dataset':['CV','Test'],'AUC':[res['auc'],test_res['auc']]}).to_csv(f"{OUTPUT_DIR}/final_summary.csv", index=False)

if __name__ == "__main__":
    main()





