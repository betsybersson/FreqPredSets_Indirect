#################################
## Replication codes for eBird data application in Bersson and Hoff 2023b
#################################
#################################
library(tidyverse)
library(parallel)
library(reshape2)
library(maps)
library(ggmap)
library(xtable)
library(forcats)
library(ggthemes)

source("./pred_functions.R")
source("./helper_functions.R")
#################################

#################################
## prediction error
alpha = .05
#################################


#################################
## Preliminaries
#################################

## read in data
df = readRDS("./data/ebird_NC_May2023.RDS")

# get some parameters
county.names = rownames(df)
n.counties = length(county.names)
bird.names = colnames(df)
n.birds = length(bird.names)

## set up base map
states <- map_data("state")
state_df <- subset(states, region %in% c("north carolina"))
counties <- map_data("county")
county_df <- subset(counties, region == "north carolina")

# get county info
countydf = read_csv("./data/NC_fips_and_locs.txt") %>% 
  mutate(county = paste0("US-NC-",COUNTYFP)) %>% 
  select(-STATEFP,-COUNTYFP,-STNAME,-POPULATION) %>%
  filter(county %in% county.names) %>% 
  arrange(factor(county,levels=county.names)) %>% 
  mutate(subregion = tolower(COUNAME))

# get county code joined to poly files
county_df = left_join(county_df, select(countydf,subregion,county),
                      by="subregion")
#################################


#################################
## EDA
#################################

# analyis within species
quantile(colSums(df),c(.5))
# analysis within counties
range(rowSums(df))

## plot sample size map
# lets get fab matrix into tidy form
N = rowSums(df)
N.df = tibble(N) %>% 
  mutate(county = names(N), N = log(N)) %>% 
  data.frame()

county_df_bird = left_join(county_df,N.df,
                           by = c("county" = "county"))

ss_plot = ggplot(data = state_df,
            mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(1.3) + 
  geom_polygon(color = "grey",lwd = 0.25) + 
  geom_polygon(data = county_df_bird,
               color = "grey",lwd = 0.15,
               aes(fill = N)) +
  scale_fill_gradient2(midpoint = median(log(N)),
                       low ="red",mid = "white", high = "#0000ff",
                       name = "log(N)") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title= element_blank(), 
        legend.position="right",
        plot.title = element_text(hjust = 0.5,size=15,face="bold",family="Times"),
        legend.text = element_text(size=15,family="Times"),
        plot.margin=grid::unit(c(0,0,0,0), "mm")) 
ss_plot
#################################

#################################
## Obtain direct and indirect prediction sets
#################################
## obtain direct prediction sets for each county
direct.set = matrix(0,nrow = n.counties, ncol = n.birds)
colnames(direct.set) = bird.names
rownames(direct.set) = county.names
for (j in 1:n.counties){
  direct.set.temp = MN_indir_pred(df[j,],rep(0,n.birds),
                                  alpha,bird.names)$pred
  
  direct.set[j,which(bird.names %in% direct.set.temp)] = 1
}

## obtain indirect prediction sets for each county

## hyperprior parameter estimation
## first obtain 5 nearest neighbors of each county
Adj.ind = get_nn(countydf$LATITUDE,countydf$LONGITUDE,
                 n.nn = 5,
                 countydf$county)

## run EM algorithm in parallel across counties
gamma0.spatial = mclapply(1:n.counties,
                          function(j)polyaMLE(df[Adj.ind[j,],],
                                              method = "separate")) %>%
  unlist() %>%
  matrix(nrow = n.counties,byrow=T)

## get indirect sets now
indirect.set = matrix(0,nrow = n.counties, ncol = n.birds)
colnames(indirect.set) = bird.names
rownames(indirect.set) = county.names
for (j in 1:n.counties){
  gamma.j = gamma0.spatial[j,]
  
  set.temp = MN_indir_pred(df[j,],gamma.j,
                           alpha,bird.names)$pred
  
  indirect.set[j,which(bird.names %in% set.temp)] = 1
}
#################################

#################################
## Cardinality comparisons
#################################
## plot cardinality ratio
N.cat = ifelse(N>quantile(N,.25),"",as.character(round(N)))
card.mat = cbind(rowSums(direct.set),rowSums(indirect.set))
colnames(card.mat) = c("direct","indirect")
card.mat = data.frame(as.data.frame(card.mat),"N" = N.cat)

# lets get cardinality matrix into tidy form
card.df = card.mat %>% 
  mutate(Ratio = indirect/direct) %>% 
  select(Ratio,N) %>% 
  mutate(county = rownames(card.mat)) 

centroid <- aggregate(cbind(long,lat) ~ county, data=county_df, FUN=mean) 
names(centroid)[2:3] = c("cent.lon","cent.lat")

county_df_bird = left_join(county_df,card.df,
                           by = c("county" = "county")) %>% 
  left_join(centroid)

ratio_plot = ggplot(data = state_df,
            mapping = aes(x = long, y = lat, group = group)) + 
  coord_fixed(1.3) + 
  geom_polygon(color = "grey",lwd = 0.25) + 
  geom_polygon(data = county_df_bird,
               color = "grey",lwd = 0.15,
               aes(fill = Ratio)) +
  scale_fill_gradient2(midpoint = 1,limits = c(0,1.5),
                       low ="red",mid = "white", high = "#0000ff") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_blank(),
        axis.text = element_blank(), axis.ticks = element_blank(),
        axis.title= element_blank(), 
        legend.position="right",
        plot.title = element_text(hjust = 0.5,size=15,face="bold",family="Times"),
        legend.text = element_text(size=15,family="Times"),
        plot.margin=grid::unit(c(0,0,0,0), "mm")) +
  geom_text(data = county_df_bird,size = 3,family="Times",
            aes(x=cent.lon,y=cent.lat,
                group = county, label = N))
ratio_plot

## analysis of direct vs indirect cardinality
nrow(filter(card.df,Ratio>1)) # number of counties where indirect bigger
nrow(filter(card.df,Ratio<1)) # number of counties where indirect smaller
nrow(filter(card.df,Ratio==1)) # number of counties where indirect and default equal

filter(card.mat,indirect==n.birds|direct==n.birds) %>% arrange(as.numeric(N)) # trivial sets and sample size
80/n.birds*100 # perc of first non trivial fab set length

## get number of indirect and direct sets that totally agree
len.intersect = c()
for ( j in 1:n.counties){
  len.intersect = c(len.intersect, sum((direct.set[j,]==1) &(indirect.set[j,]==1)))
}
cbind(select(card.mat,-N),N,len.intersect) %>% 
  mutate(Ratio = direct/indirect) %>% 
  filter(Ratio == 1) %>% 
  filter((len.intersect == direct) & (len.intersect == indirect)) %>% nrow()


#################################

#################################
## In-depth county analysis
#################################
## select county
county = "US-NC-155"  # "US-NC-087"  #
county.name = "Robeson" # "Haywood" # 
#################################

ind = which(rownames(df)==county)
county.compare = data.frame("X" = df[ind,],
                            "gamma" = gamma0.spatial[ind,],
                            "Species" = colnames(df),
                            "dir.incl" = direct.set[ind,],
                            "fab.incl" = indirect.set[ind,])
rownames(county.compare) = NULL

# summary
card.mat[ind,]

county.compare = county.compare %>% 
  mutate(XpA = (X + gamma)/(N[ind] + sum(gamma0.spatial[ind,]) ),
         X = X/N[ind] ) %>% # standardize
  arrange((XpA)) %>%
  select(Species,X,XpA,dir.incl,fab.incl) %>% 
  filter(!((dir.incl==0) & (fab.incl==0))) %>% 
  mutate(Species = fct_reorder(Species, (XpA))) %>%  # sort in order of emp distribution
  mutate(cols = case_when((dir.incl==1) & (fab.incl==1) ~ "purple",
                          (fab.incl==1) == 1 ~ "red",
                          (dir.incl==1) == 1 ~ "blue",
                          TRUE ~ "gray"))

ggplot(melt(select(county.compare,Species,X,XpA),variable.name = "Method",
            value.name = "Value"), 
       aes(x = Species,y = (Value),fill= Method)) +   # Fill column
  geom_col(width = .6,position = "dodge") +   # draw the bars
  theme_tufte() +
  theme(plot.title = element_text(hjust = .5),
        axis.ticks = element_blank(),
        axis.text.x = element_text(color = (county.compare$cols),
                                   angle = 45, hjust = 1,vjust = 1.13,
                                   family="Times"),
        legend.title = element_blank(),
        plot.margin=grid::unit(c(.2,0,0,8.4), "mm"),
        axis.title.x = element_blank(),
        axis.title.y = element_text(family="Times"),
        ) + 
  ylab("Empirical Probability Mass") +
  scale_fill_manual(values=c("blue", "red"),
                    labels=c("MLE","Post.Pred."))

BIRDs = c("Eastern Kingbird","Pine Warbler","Double-crested Cormorant","Chipping Sparrow")
bird.perc = row_standardize(df)[c(ind,Adj.ind[ind,]),which(bird.names%in%BIRDs)]*100 %>% round(2)
rownames(bird.perc) = rownames(df[c(ind,Adj.ind[ind,]),])
bird.perc2 = apply(bird.perc,2,function(j)paste0(round(j,2),"%"))
rownames(bird.perc2) = rownames(bird.perc)

xtable(rbind(bird.perc2,
       round(gamma0.spatial[ind,which(bird.names%in%BIRDs)],2))[,c(1,4,3,2)]
)
#################################
