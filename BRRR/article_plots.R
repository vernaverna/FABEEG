# collection of plots for teh article.

setwd("/projects/FABEEG/BRRR")
library("penalizedLDA")
library("ggplot2")
library("dplyr")
library("cvms")
library("RColorBrewer")
library("corrplot")
library("viridis")

# For viridis color scheme, see e.g. 
#   https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html

# # # # # # # # # # # # # # # #
#   load visualization data   #
# # # # # # # # # # # # # # # #

load("results/full/NEW_alldata_2N2_BRRR_K12.RData")
Y <- res$data$phenotypes
X <- res$data$genotypes
subj <- row.names(X)
ages = read.csv('data/new_age_df.csv')
ages <- ages[,-1]
ages[ages==" "] <- NA #replace empty strings with NA-values
ages = ages[ages$File %in% unique(subj)==TRUE, ] #to get rid of some extra subjects that should not be there
ages['Sex'] <- as.factor(ages$Sex)

age_df <- na.omit(ages)

inv_G <- res$scaling #inv(average(Gamma))
lat_space <- Y%*%inv_G #obs x latent components
lat_map = as.data.frame(rbind(lat_space))

# ========================================================
# FIGURE 6
# make a stacked histogram with subject info
age_df <- ages %>% mutate(age_group = cut(Age, breaks=14))
age_df <- na.omit(age_df)
prop_df <- age_df %>% group_by(age_group, Sex) %>% summarise(n = n()) %>% mutate(prop=n/sum(n) )
prop_df['prop'] <- round(prop_df$prop, 2)
prop_df$prop[prop_df$Sex=='M'] <- NA
prop_df <- prop_df %>% mutate(total=sum(n) )
age_df <- merge(x=age_df,y=prop_df, by=c('age_group', 'Sex'), all=TRUE)

pg <- ggplot(age_df, aes(x=age_group, fill=Sex)) + 
  geom_bar(position='stack') + theme_bw() + ylab("Frequency") + xlab('Age') +
  geom_text(aes(y=total, label=prop), vjust=-0.4, color='gray26', size=5) + 
  scale_fill_manual(values=c("#DCE319FF", "#33638DFF")) + # scale_fill_hue(l=70, c=90) + 
  scale_x_discrete(labels=c("0","","","4","","","8","","","12","","","16","")) +
  theme(legend.key = element_rect(linewidth = 16), 
        legend.text = element_text(size = 13),
        legend.title = element_text(face = "bold", size = 18),
        axis.text.x = element_text(size=13), axis.text.y = element_text(size=13),
        axis.title.x = element_text(size=18),
        axis.title.y = element_text(size=18),
        panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())

ggsave(file='demographics.pdf', plot=pg, width=10, height=7) 
# NONNIH.

# =================================================================================

# FIGURE 5
# tSNE projection of the subspace colored by age and sex.

library("Rtsne")
D <- normalize_input(lat_space)
tsne <- Rtsne(D, perplexity=200, theta=0.0, check_duplicates = FALSE, max_iter=2000) 

nsubj <- length(unique(subj))
reps = 2
# adding together N2 mappings  plus some covariates
lat_map = as.data.frame(rbind(lat_space))
lat_map['X1'] = tsne$Y[,1]
lat_map['X2'] = tsne$Y[,2]
lat_map['spectra'] = c(rep('N1B', nsubj), rep('N1A', nsubj))
lat_map['age'] = rep(ages$Age, reps) #get age data
lat_map['group'] = rep(round(ages$Age, 0), reps) #get age data
lat_map['sex'] = rep(ages$Sex, reps)
lat_map['cap'] = rep(ages$Cap, reps)
lat_map['subject'] = subj[1:nsubj]



p <- ggplot(data=lat_map, aes(x=V2, y=V3, shape=spectra, colour=age)) + 
            geom_point(alpha=0.3, size=2) +
            ggtitle("Subjects in latent mapping ") + 
            scale_color_viridis() +
            xlab("t-SNE #1") + ylab("t-SNE #2") + #xlim(-42,48) + ylim(-20,25) +
            theme(legend.position="none") + theme_bw() + 
            theme(panel.border = element_blank(), 
                  panel.grid.major = element_blank(),
                  panel.grid.minor = element_blank(), 
                  axis.line = element_line(colour = "black"))
p
ggsave("figures/NEW_latmap_all_2N2_Age.pdf", width=7.4, height=5.2)

q <- ggplot(data=lat_map, aes(x=V2, y=V3, shape=spectra, colour=sex)) + 
  geom_point(alpha=0.5, size=3) +
  ggtitle("Subjects in latent mapping ") + 
  scale_color_viridis(discrete=TRUE) +
  xlab("t-SNE #1") + ylab("t-SNE #2") + #xlim(-42,48) + ylim(-20,25) +
  theme(legend.position="none") + theme_bw() + 
  theme(panel.border = element_blank(), 
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black"))
q
ggsave("figures/NEW_latmap_all_2N2_Sex.pdf", width=7.4, height=5.2)









