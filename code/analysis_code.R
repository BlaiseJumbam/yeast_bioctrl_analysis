source("e:/professdev/r_projects/r_packages.R")

# Set preferred packages for the functions that conflict ------------------
conflicts_prefer(ggpubr::mutate)
conflicts_prefer(dplyr::rename)
conflicts_prefer(dplyr::filter) 
conflicts_prefer(dplyr::select)
conflicts_prefer(dplyr::summarise)
conflicts_prefer(dplyr::arrange)
conflicts_prefer(dplyr::desc)

# check sheet names in excel spreadsheet ----------------------------------
excel_sheets('data/all_data.xlsx')

# Confrontation Assay Analysis --------------------------------------------

# import, clean, and export data ------------------------------------------
read_excel("data/all_data.xlsx", sheet=1) |>
  mutate(treatmnt=case_when(isolate_id=='C'~'CTRL',
                           isolate_id=='357'~'JB543', # Pichia
                           isolate_id=='461'~'JB524', # Pichia
                           isolate_id=='506'~'JB624', # Zygoascus
                           isolate_id=='707'~'JB767', # Pichia
                           isolate_id=='957'~'JB951'),# Aureobasidium
         diam_mm=(radius_mm*2)-5) |>  # subtract initial 5mm disc used
  select(-c(inoc_date,col_date,isolate_id,radius_mm)) |> 
  relocate(treatmnt, .before=pathogen) |> #changes position of isolate
  group_by(pathogen, experiment) |> 
  mutate(max_diam=max(diam_mm),
         perc_in=100 * (max_diam-diam_mm)/max_diam,
         log_fdiam=log10(diam_mm),
         sqrt_fdiam=sqrt(diam_mm)) |>
  ungroup() |> 
  rename(iz=inhibition_zone,ydiam_mm=yeast_diameter) |>
  write.xlsx('output/bccf_antagon_clean_exp1-3.xlsx')

# import cleaned data -----------------------------------------------------

# Botrytis cinerea --------------------------------------------------------

# subset data by experiment and explore
bc_antagon_clean1 <- read_xlsx('output/bccf_antagon_clean_exp1-3.xlsx') |> 
  filter(pathogen=='Bc') |>
  filter(experiment==1) |> 
  select(-c(pathogen,dpi,iz,ydiam_mm,diam_mm,max_diam))

bc_antagon_clean2 <- read_xlsx('output/bccf_antagon_clean_exp1-3.xlsx') |> 
  filter(pathogen=='Bc') |>
  filter(experiment==2) |> 
  select(-c(pathogen,dpi,iz,ydiam_mm,diam_mm,max_diam))

bc_antagon_clean3 <- read_xlsx('output/bccf_antagon_clean_exp1-3.xlsx') |> 
  filter(pathogen=='Bc') |>
  filter(experiment==3) |> 
  select(-c(pathogen,dpi,iz,ydiam_mm,diam_mm,max_diam))

# check variables and make appropriate conversions ------------------------
str(bc_antagon_clean3)
bc_antagon_clean3$treatmnt <- as.factor(bc_antagon_clean3$treatmnt)
glimpse(bc_antagon_clean3)

# perform eda -------------------------------------------------------------
# summary statistics
bc_antagon_clean1 |> 
  group_by(treatmnt, experiment) |> 
  get_summary_stats(perc_in, type='mean_sd')

# boxplot
bc_antagon_clean1 %>%
  group_by(treatmnt,experiment) %>% 
  ggboxplot(., x='treatmnt', y='perc_in') + 
  facet_wrap(~experiment)

# histogram
bc_antagon_clean1 |> 
  group_by(treatmnt,experiment)%>%
  gghistogram(., x='perc_in') +
  facet_wrap(~experiment)

# identify outliers
bc_antagon_clean1 |> 
  group_by(treatmnt, experiment) |> 
  identify_outliers(perc_in)  # two value in CTRL identified as extreme values

# check for normality -----------------------------------------------------

bc_antagon_clean1_aov <- aov(perc_in ~ treatmnt, data=bc_antagon_clean1)
summary(bc_antagon_clean1_aov)

# qqplot
ggqqplot(residuals(bc_antagon_clean1_aov)) # one way
plot(bc_antagon_clean1_aov, 2) # second way

# shapiro test
shapiro_test(residuals(bc_antagon_clean1_aov)) # p is not signif

bc_antagon_clean1_aov_residuals <- residuals(object=bc_antagon_clean1_aov)
shapiro.test(x=bc_antagon_clean1_aov_residuals)

# check shapiro for each treatment level
bc_antagon_clean1 |>
  group_by(treatmnt) |> 
  shapiro_test(perc_in) # p here is also signif for all levels

# qqplot for each group level 
ggqqplot(bc_antagon_clean1, 'perc_in', facet.by='treatmnt')

# homogeneity of variance
plot(bc_antagon_clean1_aov, 1)
leveneTest(perc_in ~ treatmnt, data = bc_antagon_clean1)

# check for skewness and kurtosis
describe(bc_antagon_clean1$perc_in, na.rm = T, type = 2)

# kruskal wallis for non-normally distributed data
# edited for all experiments
kruskal.test(perc_in ~ treatmnt, data = bc_antagon_clean1) #signif

# post hoc test
# edited code for all experiments
bc_antagon_pwt1 <- pairwise.wilcox.test(bc_antagon_clean1$perc_in, 
                                        bc_antagon_clean1$treatmnt,
                                        p.adjust.method = "BH")

bc_antagon_pwt_table1 <- fullPTable(bc_antagon_pwt1$p.value) # table of p-values
bc_antagon_letters1 <- multcompLetters(bc_antagon_pwt_table1) # signif letters

# summarize the data to plot ----------------------------------------------
# edited code for all experiments
bc_antagon_treatmnt_summry1 <- bc_antagon_clean1 |>  
  group_by(treatmnt) |>  
  summarise(mean_perc_in = mean(perc_in),
            sd_perc_in=sd(perc_in), 
            n_perc_in = n(),
            se_perc_in = sd(perc_in)/sqrt(n()))
view(bc_antagon_treatmnt_summry1)

# convert cld to list and merge with summarized data ----------------------
# edited for all experiments
bc_antagon_cld1 <- as.data.frame.list(bc_antagon_letters1) # converts to list
bc_antagon_treatmnt_summry1$wilcox_cl <- bc_antagon_cld1$Letters # add lets-sumary

# save the summarized and merged data -------------------------------------
# edited this for all experiments
bc_antagon_treatmnt_summry1 |> 
  arrange(desc(mean_perc_in)) |>  
  write_csv("output/bc_antagon_treatmnt_summry1.csv")

# plotting percent confrontation of Bc by yeasts --------------------------
# edited code for all experiments
bc_antagon_plot1 <- bc_antagon_treatmnt_summry1 %>%
  arrange(desc(mean_perc_in)) %>%
  ggplot(., aes(x=factor(treatmnt, 
                         levels=c("CTRL","JB951","JB543","JB767","JB524")),
                y=mean_perc_in)) +
  geom_bar(stat="identity", position="dodge", alpha=0.7,
           color="#000", fill="#f8f8ff", width=0.7) +
  geom_errorbar(aes(ymin = mean_perc_in - se_perc_in, 
                    ymax = mean_perc_in + se_perc_in),
                position = position_dodge(0.9), width=0.25,
                show.legend=FALSE) +
  geom_text(aes(label=wilcox_cl), position=position_dodge(0.9), size=4, 
            vjust=-2.3, hjust=0.3, color="#000") +
  scale_y_continuous(limits = c(0, 101), expand = c(0,0)) +
  theme_classic() +
  theme(
    axis.title=element_text(size=12, colour="#000", family='sans'),
    axis.title.x=element_blank(),
    axis.text=element_text(size=11, colour="#000", family='sans'),
    axis.ticks.x=element_blank(),
    axis.ticks.y=element_line(color="#000")
  ) +
  labs(y="MIR (%)", x = "")

ggsave("output/bc_antagon_plot1.tif",bc_antagon_plot1,height=4.5,width=4,dpi=300)
ggsave("output/bc_antagon_plot1.pdf",bc_antagon_plot1,height=4.5,width=4,dpi=300)

# pooled confrontation assay analysis for Bc ------------------------------

bc_antagon_clean <- read_xlsx('output/bccf_antagon_clean_exp1-3.xlsx') |> 
  filter(pathogen=='Bc') |>
  select(-c(pathogen,dpi,iz,ydiam_mm,diam_mm,max_diam,experiment))

# check variables and make appropriate conversions ------------------------
str(bc_antagon_clean)
bc_antagon_clean$treatmnt <- as.factor(bc_antagon_clean$treatmnt)
glimpse(bc_antagon_clean)

# perform eda -------------------------------------------------------------

# boxplot
bc_antagon_clean %>%
  group_by(treatmnt) %>% 
  ggboxplot(., x='treatmnt', y='perc_in') # we see 2 outliers in the CTRL

# histogram
bc_antagon_clean |> 
  group_by(treatmnt)%>%
  gghistogram(., x='perc_in')

# identify outliers and extreme values
bc_antagon_clean |> 
  group_by(treatmnt) |> 
  identify_outliers(perc_in)  # two identified as extreme values in the CTRL

# check for normality -----------------------------------------------------

bc_antagon_clean_aov <- aov(perc_in ~ treatmnt, data=bc_antagon_clean)
summary(bc_antagon_clean_aov)

# qqplot
ggqqplot(residuals(bc_antagon_clean_aov)) # one way
plot(bc_antagon_clean_aov, 2) # second way

# shapiro test
shapiro_test(residuals(bc_antagon_clean_aov)) # p is not signif

bc_antagon_clean_aov_residuals <- residuals(object=bc_antagon_clean_aov)
shapiro.test(x=bc_antagon_clean_aov_residuals)

# check shapiro for each treatment level
bc_antagon_clean |>
  group_by(treatmnt) |> 
  shapiro_test(perc_in) # p here is also signif for all levels

# qqplot for each group level 
ggqqplot(bc_antagon_clean, 'perc_in', facet.by='treatmnt')

#homogeneity of variance
plot(bc_antagon_clean_aov, 1)
leveneTest(perc_in ~ treatmnt, data = bc_antagon_clean)

#check for skewness and kurtosis
describe(bc_antagon_clean$perc_in, na.rm = T, type = 2)

#kruskal wallis for non-normally distributed data
kruskal.test(perc_in ~ treatmnt, data = bc_antagon_clean) #signif

#post hoc test
bc_antagon_pwt <- pairwise.wilcox.test(bc_antagon_clean$perc_in, 
                                        bc_antagon_clean$treatmnt,
                                        p.adjust.method = "BH")

bc_antagon_pwt_table <- fullPTable(bc_antagon_pwt$p.value) # table of p-values
bc_antagon_letters <- multcompLetters(bc_antagon_pwt_table) # signif letters

# summarize the data to plot ----------------------------------------------
bc_antagon_treatmnt_summry <- bc_antagon_clean |>  
  group_by(treatmnt) |>  
  summarise(mean_perc_in = mean(perc_in),
            sd_perc_in=sd(perc_in), 
            n_perc_in = n(),
            se_perc_in = sd(perc_in)/sqrt(n()))
view(bc_antagon_treatmnt_summry)

# convert cld to list and merge with summarized data ----------------------
bc_antagon_cld <- as.data.frame.list(bc_antagon_letters) # converts to list
bc_antagon_treatmnt_summry$wilcox_cl <- bc_antagon_cld$Letters # add lets-sumary

# save the summarized and merged data -------------------------------------
bc_antagon_treatmnt_summry |> 
  arrange(desc(mean_perc_in)) |>  
  write_csv("output/bc_antagon_treatmnt_summry.csv")

# plotting pooled percent confrontation of Bc by yeasts -------------------

bc_antagon_plot <- bc_antagon_treatmnt_summry %>%
  arrange(desc(mean_perc_in)) %>%
  ggplot(., aes(x=factor(treatmnt, 
                         levels=c("CTRL","JB951","JB543","JB767","JB524")),
                y=mean_perc_in)) +
  geom_bar(stat="identity", position="dodge", alpha=0.7,
           color="#000", fill="#f8f8ff", width=0.7) +
  geom_errorbar(aes(ymin = mean_perc_in - se_perc_in, 
                    ymax = mean_perc_in + se_perc_in),
                position = position_dodge(0.9), width=0.25,
                show.legend=FALSE) +
  geom_text(aes(label=wilcox_cl), position=position_dodge(0.9), size=4, 
            vjust=-2.3, hjust=0.3, color="#000") +
  scale_y_continuous(limits = c(0, 101), expand = c(0,0)) +
  theme_classic() +
  theme(
    axis.title=element_text(size=12, colour="#000", family='sans'),
    axis.title.x=element_blank(),
    axis.text=element_text(size=11, colour="#000", family='sans'),
    axis.ticks.x=element_blank(),
    axis.ticks.y=element_line(color="#000")
  ) +
  labs(y="MIR (%)", x = "")

ggsave("output/bc_antagon_plot.tif",bc_antagon_plot,height=4.5,width=4,dpi=300)
ggsave("output/bc_antagon_plot.pdf",bc_antagon_plot,height=4.5,width=4,dpi=300)


# Colletotrichum fioriniae ------------------------------------------------

# subset data by experiment and explore
cf_antagon_clean1 <- read_xlsx('output/bccf_antagon_clean_exp1-3.xlsx') |> 
  filter(pathogen=='Cf') |>
  filter(experiment==1) |> 
  select(-c(pathogen,dpi,iz,ydiam_mm,diam_mm,max_diam))
# exp2
cf_antagon_clean2 <- read_xlsx('output/bccf_antagon_clean_exp1-3.xlsx') |> 
  filter(pathogen=='Cf') |>
  filter(experiment==2) |> 
  select(-c(pathogen,dpi,iz,ydiam_mm,diam_mm,max_diam))
# exp3
cf_antagon_clean3 <- read_xlsx('output/bccf_antagon_clean_exp1-3.xlsx') |> 
  filter(pathogen=='Cf') |>
  filter(experiment==3) |> 
  select(-c(pathogen,dpi,iz,ydiam_mm,diam_mm,max_diam))

# check variables and make appropriate conversions ------------------------
# edited for all experiments
str(cf_antagon_clean1)
cf_antagon_clean1$treatmnt <- as.factor(cf_antagon_clean1$treatmnt)
glimpse(cf_antagon_clean1)

# perform eda -------------------------------------------------------------
# summary statistics
cf_antagon_clean1 |> 
  group_by(treatmnt, experiment) |> 
  get_summary_stats(perc_in, type='mean_sd')

# boxplot
cf_antagon_clean1 %>%
  group_by(treatmnt,experiment) %>% 
  ggboxplot(., x='treatmnt', y='perc_in') + 
  facet_wrap(~experiment)

# histogram
cf_antagon_clean1 |> 
  group_by(treatmnt,experiment)%>%
  gghistogram(., x='perc_in') +
  facet_wrap(~experiment)

# identify outliers
cf_antagon_clean1 |> 
  group_by(treatmnt, experiment) |> 
  identify_outliers(perc_in)  # no extreme values

# check for normality -----------------------------------------------------
# edited for all experiments
cf_antagon_clean1_aov <- aov(perc_in ~ treatmnt, data=cf_antagon_clean1)
summary(cf_antagon_clean1_aov)

# qqplot
ggqqplot(residuals(cf_antagon_clean1_aov)) # one way
plot(cf_antagon_clean1_aov, 2) # second way

# shapiro test
shapiro_test(residuals(cf_antagon_clean1_aov)) # p is not signif

cf_antagon_clean1_aov_residuals <- residuals(object=cf_antagon_clean1_aov)
shapiro.test(x=cf_antagon_clean1_aov_residuals)

# check shapiro for each treatment level
cf_antagon_clean1 |>
  group_by(treatmnt) |> 
  shapiro_test(perc_in) # p here is also signif for all levels

# qqplot for each group level 
ggqqplot(cf_antagon_clean1, 'perc_in', facet.by='treatmnt')

# homogeneity of variance
plot(cf_antagon_clean1_aov, 1)
leveneTest(perc_in ~ treatmnt, data = cf_antagon_clean1)

# check for skewness and kurtosis
describe(cf_antagon_clean1$perc_in, na.rm = T, type = 2)

# kruskal wallis for non-normally distributed data
# edited for all experiments
kruskal.test(perc_in ~ treatmnt, data = cf_antagon_clean1) #signif

# post hoc test
# edited for all experiments
cf_antagon_pwt1 <- pairwise.wilcox.test(cf_antagon_clean1$perc_in, 
                                       cf_antagon_clean1$treatmnt,
                                       p.adjust.method = "BH")

cf_antagon_pwt_table1 <- fullPTable(cf_antagon_pwt1$p.value) # table of p-values
cf_antagon_letters1 <- multcompLetters(cf_antagon_pwt_table1) # signif letters

# summarize the data to plot ----------------------------------------------
# edited for all experiments
cf_antagon_treatmnt_summry1 <- cf_antagon_clean1 |>  
  group_by(treatmnt) |>  
  summarise(mean_perc_in = mean(perc_in),
            sd_perc_in=sd(perc_in), 
            n_perc_in = n(),
            se_perc_in = sd(perc_in)/sqrt(n()))
view(cf_antagon_treatmnt_summry1)

# convert cld to list and merge with summarized data ----------------------
# edited for all experiments
cf_antagon_cld1 <- as.data.frame.list(cf_antagon_letters1) # converts to list
cf_antagon_treatmnt_summry1$wilcox_cl <- cf_antagon_cld1$Letters 

# save the summarized and merged data -------------------------------------
# edited for all experiments
cf_antagon_treatmnt_summry1 |> 
  arrange(desc(mean_perc_in)) |>  
  write_csv("output/cf_antagon_treatmnt_summry1.csv")

# plotting percent confrontation of Cf by yeasts --------------------------
# edited code for all experiments
cf_antagon_plot1 <- cf_antagon_treatmnt_summry1 %>%
  arrange(desc(mean_perc_in)) %>%
  ggplot(., aes(x=factor(treatmnt, 
                         levels=c("CTRL","JB951","JB543","JB767","JB624")),
                y=mean_perc_in)) +
  geom_bar(stat="identity", position="dodge", alpha=0.7,
           color="#000", fill="#f8f8ff", width=0.7) +
  geom_errorbar(aes(ymin = mean_perc_in - se_perc_in, 
                    ymax = mean_perc_in + se_perc_in),
                position = position_dodge(0.9), width=0.25,
                show.legend=FALSE) +
  geom_text(aes(label=wilcox_cl), position=position_dodge(0.9), size=4, 
            vjust=-2.3, hjust=0.3, color="#000") +
  scale_y_continuous(breaks=c(0,20,40,60,80,100),
                     limits=c(0,101),expand = c(0,0)) +
  theme_classic() +
  theme(
    axis.title=element_text(size=12, colour="#000", family='sans'),
    axis.title.x=element_blank(),
    axis.text=element_text(size=11, colour="#000", family='sans'),
    axis.ticks.x=element_blank(),
    axis.ticks.y=element_line(color="#000")
  ) +
  labs(y="MIR (%)", x = "")

ggsave("output/cf_antagon_plot1.tif",cf_antagon_plot1,height=4.5,width=4,dpi=300)
ggsave("output/cf_antagon_plot1.pdf",cf_antagon_plot1,height=4.5,width=4,dpi=300)

# pooled confrontation assay analysis for Cf ------------------------------

cf_antagon_clean <- read_xlsx('output/bccf_antagon_clean_exp1-3.xlsx') |> 
  filter(pathogen=='Cf') |>
  select(-c(pathogen,dpi,iz,ydiam_mm,diam_mm,max_diam,experiment))

# perform eda -------------------------------------------------------------

# boxplot
cf_antagon_clean %>%
  group_by(treatmnt) %>% 
  ggboxplot(., x='treatmnt', y='perc_in') # we see 2 outliers in the CTRL

# histogram
cf_antagon_clean |> 
  group_by(treatmnt)%>%
  gghistogram(., x='perc_in')

# identify outliers and extreme values
cf_antagon_clean |> 
  group_by(treatmnt) |> 
  identify_outliers(perc_in) # no extreme values

# check for normality -----------------------------------------------------

cf_antagon_clean_aov <- aov(perc_in ~ treatmnt, data=cf_antagon_clean)
summary(cf_antagon_clean_aov)

# qqplot
ggqqplot(residuals(cf_antagon_clean_aov)) # one way
plot(cf_antagon_clean_aov, 2) # second way

# shapiro test
shapiro_test(residuals(cf_antagon_clean_aov)) # p is not signif

cf_antagon_clean_aov_residuals <- residuals(object=cf_antagon_clean_aov)
shapiro.test(x=cf_antagon_clean_aov_residuals)

# check shapiro for each treatment level
cf_antagon_clean |>
  group_by(treatmnt) |> 
  shapiro_test(perc_in) # p here is also signif for all levels

# qqplot for each group level 
ggqqplot(cf_antagon_clean, 'perc_in', facet.by='treatmnt')

#homogeneity of variance
plot(cf_antagon_clean_aov, 1)
leveneTest(perc_in ~ treatmnt, data = cf_antagon_clean)

#check for skewness and kurtosis
describe(cf_antagon_clean$perc_in, na.rm = T, type = 2)

#kruskal wallis for non-normally distributed data
kruskal.test(perc_in ~ treatmnt, data = cf_antagon_clean) #signif

#post hoc test
cf_antagon_pwt <- pairwise.wilcox.test(cf_antagon_clean$perc_in, 
                                       cf_antagon_clean$treatmnt,
                                       p.adjust.method = "BH")

cf_antagon_pwt_table <- fullPTable(cf_antagon_pwt$p.value) # table of p-values
cf_antagon_letters <- multcompLetters(cf_antagon_pwt_table) # signif letters

# summarize the data to plot ----------------------------------------------
cf_antagon_treatmnt_summry <- cf_antagon_clean |>  
  group_by(treatmnt) |>  
  summarise(mean_perc_in = mean(perc_in),
            sd_perc_in=sd(perc_in), 
            n_perc_in = n(),
            se_perc_in = sd(perc_in)/sqrt(n()))
view(cf_antagon_treatmnt_summry)

# convert cld to list and merge with summarized data ----------------------
cf_antagon_cld <- as.data.frame.list(cf_antagon_letters) # converts to list
cf_antagon_treatmnt_summry$wilcox_cl <- cf_antagon_cld$Letters # add lets-sumary

# save the summarized and merged data -------------------------------------
cf_antagon_treatmnt_summry |> 
  arrange(desc(mean_perc_in)) |>  
  write_csv("output/cf_antagon_treatmnt_summry.csv")


# plotting pooled percent confrontation of Cf by yeasts -------------------

cf_antagon_plot <- cf_antagon_treatmnt_summry %>%
  arrange(desc(mean_perc_in)) %>%
  ggplot(., aes(x=factor(treatmnt, 
                         levels=c("CTRL","JB951","JB543","JB767","JB624")),
                y=mean_perc_in)) +
  geom_bar(stat="identity", position="dodge", alpha=0.7,
           color="#000", fill="#f8f8ff", width=0.7) +
  geom_errorbar(aes(ymin = mean_perc_in - se_perc_in, 
                    ymax = mean_perc_in + se_perc_in),
                position = position_dodge(0.9), width=0.25,
                show.legend=FALSE) +
  geom_text(aes(label=wilcox_cl), position=position_dodge(0.9), size=4, 
            vjust=-2.3, hjust=0.3, color="#000") +
  scale_y_continuous(limits = c(0, 101), expand = c(0,0)) +
  theme_classic() +
  theme(
    axis.title=element_text(size=12, colour="#000", family='sans'),
    axis.title.x=element_blank(),
    axis.text=element_text(size=11, colour="#000", family='sans'),
    axis.ticks.x=element_blank(),
    axis.ticks.y=element_line(color="#000")
  ) +
  labs(y="MIR (%)", x = "")

ggsave("output/cf_antagon_plot.tif",cf_antagon_plot,height=4.5,width=4,dpi=300)
ggsave("output/cf_antagon_plot.pdf",cf_antagon_plot,height=4.5,width=4,dpi=300)


# Pathogen on Yeast Assay Analysis ----------------------------------------

# import, clean, and export data ------------------------------------------
bccf_poy_clean <- read_excel("data/all_data.xlsx", sheet=2) |> 
  mutate(treatmnt=case_when(isolate_id=='C'~'CTRL',
                           isolate_id=='357'~'JB543',
                           isolate_id=='461'~'JB524',
                           isolate_id=='506'~'JB624',
                           isolate_id=='707'~'JB767',
                           isolate_id=='957'~'JB951'),
         diam_mm=diam_mm-5) |> # subtract initial 5mm disc used
  select(-c(inoc_date,col_date,isolate_id)) |> 
  relocate(treatmnt, .before=pathogen) #|> #changes position of isolate
  #write.xlsx('output/bccf_poy_clean_exp1-2.xlsx')

# Botrytis
# subset data for two experiments
bc_dir_clean1 <- bccf_poy_clean |> 
  filter(pathogen=='Bc') |> #summarise(max_diam=max(diam_mm))
  filter(experiment==1) |> 
  mutate(max_diam=max(diam_mm)) |> 
  group_by(treatmnt) |> 
  mutate(perc_in=100*(max_diam-diam_mm)/max_diam,
         log_diam=log10(diam_mm + 1),
         sqrt_diam=sqrt(diam_mm))

bc_dir_clean2 <- bccf_poy_clean |> 
  filter(pathogen=='Bc') |> #summarise(max_diam=max(diam_mm))
  filter(experiment==2) |> 
  mutate(max_diam=max(diam_mm)) |> 
  group_by(treatmnt) |> 
  mutate(perc_in=100*(max_diam-diam_mm)/max_diam,
         log_diam=log10(diam_mm + 1),
         sqrt_diam=sqrt(diam_mm))
  
# check variables and make appropriate conversions ------------------------
# edit for both experiments
str(bc_dir_clean1)
bc_dir_clean1$treatmnt <- as.factor(bc_dir_clean1$treatmnt)
glimpse(bc_dir_clean1)  

# perform eda -------------------------------------------------------------

# edit for both experiments

# summary statistics
bc_dir_clean1 |> 
  group_by(treatmnt) |> 
  get_summary_stats(perc_in, type='mean_sd')

# boxplot
bc_dir_clean1 |> 
  ggboxplot(x='treatmnt', y='perc_in') # we see 2 outliers in the CTRL

# histogram
bc_dir_clean1 |> 
  gghistogram(x='perc_in')

# identify outliers
bc_dir_clean1 |> 
  group_by(treatmnt) |> 
  identify_outliers(perc_in)

# check for normality -----------------------------------------------------

bc_dir_clean1_aov <- aov(perc_in ~ treatmnt, data=bc_dir_clean1)
summary(bc_dir_clean1_aov)

# qqplot
ggqqplot(residuals(bc_dir_clean1_aov)) # one way
plot(bc_dir_clean1_aov, 2) # second way

# shapiro test
shapiro_test(residuals(bc_dir_clean1_aov)) # p is signif

bc_dir_clean1_aov_residuals <- residuals(object=bc_dir_clean1_aov)
shapiro.test(x = bc_dir_clean1_aov_residuals)

# qqplot for each group level 
ggqqplot(bc_dir_clean1, 'perc_in', facet.by='treatmnt')

# homogeneity of variance
plot(bc_dir_clean1_aov, 1)
leveneTest(perc_in ~ treatmnt, data = bc_dir_clean1)

# check for skewness and kurtosis
describe(bc_dir_clean1$perc_in, na.rm = T, type = 2)

# kruskal wallis for non-normally distributed data

kruskal.test(perc_in ~ treatmnt, data = bc_dir_clean1) #signif

# post hoc test
# edit for both experiments
bc_dir_pwt1 <- pairwise.wilcox.test(bc_dir_clean1$perc_in, 
                                       bc_dir_clean1$treatmnt,
                                       p.adjust.method = "BH")
bc_dir_pwt_table1 <- fullPTable(bc_dir_pwt1$p.value) # table of p-values
bc_dir_letters1 <- multcompLetters(bc_dir_pwt_table1) # sig. dif. letters

# summarize the data to plot ----------------------------------------------
# edit for both experiments
bc_dir_treatmnt_summry1 <- bc_dir_clean1 |>  
  group_by(treatmnt) |>  
  summarise(mean_perc_in = mean(perc_in),
            sd_perc_in=sd(perc_in), 
            n_perc_in = n(),
            se_perc_in = sd(perc_in)/sqrt(n()))

view(bc_dir_treatmnt_summry1)

# convert cld to list and merge with summarized data ----------------------
# edit for both experiments
bc_dir_cld1 <- as.data.frame.list(bc_dir_letters1) # converts to list
bc_dir_treatmnt_summry1$wilcox_cl <- bc_dir_cld1$Letters

# save the summarized and merged data -------------------------------------
# edit for both experiments
bc_dir_treatmnt_summry1 |> 
  arrange(desc(mean_perc_in)) |>  
  write_csv("output/bc_poy_treatmnt_summry1.csv")


# plotting percent pathogen on yeast of Bc by yeasts ----------------------

# edit for both experiments
bc_dir_plot1 <- bc_dir_treatmnt_summry1 |> 
  arrange(desc(mean_perc_in)) |> 
  ggplot(aes(factor(x=treatmnt, 
                    levels=c('CTRL','JB543','JB767','JB524','JB951')), 
             y=mean_perc_in)) + 
  geom_errorbar(aes(ymin = mean_perc_in - se_perc_in, 
                    ymax = mean_perc_in + se_perc_in),
                position = position_dodge(0.9), width=0.25,
                show.legend=FALSE) +
  geom_line(aes(group=1)) +
  geom_point() +
  geom_text(aes(label=wilcox_cl), position=position_dodge(0.9), size=4, 
            vjust=-2.3, hjust=0.3, color="#000") +
  scale_y_continuous(breaks=c(0,20,40,60,80,100),
                     limits=c(0,115),expand = c(0,0)) +
  theme_classic() +
  theme(
    axis.title=element_text(size=12, colour="#000", family='sans'),
    axis.title.x=element_blank(),
    axis.text=element_text(size=11, colour="#000", family='sans'),
    axis.ticks=element_line(color="#000"),
    panel.grid.major.x=element_line(),
    panel.grid.major.y=element_line()
  ) +
  labs(y='MRGI (%)', x='') 

ggsave("output/bc_poy_plot1.tif",bc_dir_plot1,height=4,width=4,dpi=300)
ggsave("output/bc_poy_plot1.pdf",bc_dir_plot1,height=4,width=4,dpi=300)

# pooled exp1 and 2 data for Botrytis -------------------------------------

bc_dir_clean <- bccf_poy_clean |> 
  filter(pathogen=='Bc') |> #summarise(max_diam=max(diam_mm))
  mutate(max_diam=max(diam_mm)) |> 
  group_by(treatmnt) |> 
  mutate(perc_in=100*(max_diam-diam_mm)/max_diam,
         log_diam=log10(diam_mm + 1),
         sqrt_diam=sqrt(diam_mm)) |> 
  ungroup()

# check variables and make appropriate conversions ------------------------
str(bc_dir_clean)
bc_dir_clean$treatmnt <- as.factor(bc_dir_clean$treatmnt)
glimpse(bc_dir_clean)  

# perform eda -------------------------------------------------------------
# summary statistics
bc_dir_clean |> 
  group_by(treatmnt) |> 
  get_summary_stats(perc_in, type='mean_sd')

# boxplot
bc_dir_clean |> 
  ggboxplot(x='treatmnt', y='perc_in') # we see 2 outliers in the CTRL

# histogram
bc_dir_clean |> 
  gghistogram(x='perc_in')

# identify outliers
bc_dir_clean |> 
  group_by(treatmnt) |> 
  identify_outliers(perc_in)

# check for normality -----------------------------------------------------

bc_dir_clean_aov <- aov(perc_in ~ treatmnt, data=bc_dir_clean)
summary(bc_dir_clean_aov)

# qqplot
ggqqplot(residuals(bc_dir_clean_aov)) # one way
plot(bc_dir_clean_aov, 2) # second way

# shapiro test
shapiro_test(residuals(bc_dir_clean_aov)) # p is signif

bc_dir_clean_aov_residuals <- residuals(object=bc_dir_clean_aov)
shapiro.test(x = bc_dir_clean_aov_residuals)

# check shapiro for each treatment level
bc_dir_clean |>
  group_by(treatmnt) |> 
  shapiro_test(perc_in) 

# qqplot for each group level 
ggqqplot(bc_dir_clean, 'perc_in', facet.by='treatmnt')

# homogeneity of variance
plot(bc_dir_clean_aov, 1)
leveneTest(perc_in ~ treatmnt, data = bc_dir_clean)

# check for skewness and kurtosis
describe(bc_dir_clean$perc_in, na.rm = T, type = 2)

# kruskal wallis for non-normally distributed data

kruskal.test(perc_in ~ treatmnt, data = bc_dir_clean) #signif

# post hoc test
 
bc_dir_pwt <- dunnTest(bc_dir_clean$perc_in, 
                        bc_dir_clean$treatmnt, method = "bh") 
bc_dir_cld <- cldList(P.adj ~ Comparison, data=bc_dir_pwt$res)

# summarize the data to plot ----------------------------------------------
bc_dir_treatmnt_summry <- bc_dir_clean |>  
  group_by(treatmnt) |>  
  summarise(mean_perc_in = mean(perc_in),
            sd_perc_in=sd(perc_in), 
            n_perc_in = n(),
            se_perc_in = sd(perc_in)/sqrt(n()))

# convert cld to list and merge with summarized data ----------------------
bc_dir_treatmnt_summry$dunn_cl <- bc_dir_cld$Letter 

# save the summarized and merged data -------------------------------------
bc_dir_treatmnt_summry |> 
  arrange(desc(mean_perc_in)) |>  
  write_csv("output/bc_poy_treatmnt_summry.csv")

# plotting pooled percent pathogen on yeast of Bc by yeasts ---------------

bc_poy_plot <- bc_dir_treatmnt_summry |> 
  arrange(desc(mean_perc_in)) |> 
  ggplot(aes(factor(x=treatmnt, 
                    levels=c('CTRL','JB951','JB543','JB767','JB524')), 
             y=mean_perc_in)) + 
  geom_errorbar(aes(ymin = mean_perc_in - se_perc_in, 
                    ymax = mean_perc_in + se_perc_in),
                position = position_dodge(0.9), width=0.25,
                show.legend=FALSE) +
  geom_line(aes(group=1)) +
  geom_point() +
  geom_text(aes(label=dunn_cl), position=position_dodge(0.9), size=4, 
            vjust=-2.3, hjust=0.3, color="#000") +
  scale_y_continuous(breaks=c(0,20,40,60,80,100),
                     limits=c(0,115),expand = c(0,0)) +
  theme_classic() +
  theme(
    axis.title=element_text(size=12, colour="#000", family='sans'),
    axis.title.x=element_blank(),
    axis.text=element_text(size=11, colour="#000", family='sans'),
    axis.ticks=element_line(color="#000"),
    panel.grid.major.x=element_line(),
    panel.grid.major.y=element_line()
  ) +
  labs(y='MRGI (%)', x='') 

ggsave("output/bc_poy_plot_pooled.tif",bc_poy_plot,height=4,width=4,dpi=300)
ggsave("output/bc_poy_plot_pooled.pdf",bc_poy_plot,height=4,width=4,dpi=300)


# Colletotrichum

cf_dir_clean1 <- bccf_poy_clean |> 
  filter(pathogen=='Cf') |> #summarise(max_diam=max(diam_mm))
  filter(experiment==1) |> 
  mutate(max_diam=max(diam_mm)) |> 
  group_by(treatmnt) |> 
  mutate(perc_in=100*(max_diam-diam_mm)/max_diam,
         log_diam=log10(diam_mm + 1),
         sqrt_diam=sqrt(diam_mm))

cf_dir_clean2 <- bccf_poy_clean |> 
  filter(pathogen=='Cf') |> #summarise(max_diam=max(diam_mm))
  filter(experiment==2) |> 
  mutate(max_diam=max(diam_mm)) |> 
  group_by(treatmnt) |> 
  mutate(perc_in=100*(max_diam-diam_mm)/max_diam,
         log_diam=log10(diam_mm + 1),
         sqrt_diam=sqrt(diam_mm))

# check variables and make appropriate conversions ------------------------
# edit for both experiments
str(cf_dir_clean1)
cf_dir_clean1$treatmnt <- as.factor(cf_dir_clean1$treatmnt)
glimpse(cf_dir_clean1)  

# perform eda -------------------------------------------------------------

# edit for both exp

# summary statistics
cf_dir_clean1 |> 
  group_by(treatmnt) |> 
  get_summary_stats(perc_in, type='mean_sd')

# boxplot
cf_dir_clean1 |> 
  ggboxplot(x='treatmnt', y='perc_in') # we see 2 outliers in the CTRL

# histogram
cf_dir_clean1 |> 
  gghistogram(x='perc_in')

# identify outliers
cf_dir_clean1 |> 
  group_by(treatmnt) |> 
  identify_outliers(perc_in)  # one value in JV543 is identified as ev

# check for normality -----------------------------------------------------

# edit for both exp
cf_dir_clean1_aov <- aov(perc_in ~ treatmnt, data=cf_dir_clean1)
summary(cf_dir_clean1_aov)

# qqplot
ggqqplot(residuals(cf_dir_clean1_aov)) # one way
plot(cf_dir_clean1_aov, 2) # second way

# shapiro test
shapiro_test(residuals(cf_dir_clean1_aov)) # p is signif

cf_dir_clean1_aov_residuals <- residuals(object=cf_dir_clean1_aov)
shapiro.test(x=cf_dir_clean1_aov_residuals)

# check shapiro for each treatment level
cf_dir_clean1 |>
  group_by(treatmnt) |> 
  shapiro_test(perc_in) # p here is also signif for all levels

# qqplot for each group level 
ggqqplot(cf_dir_clean1, 'perc_in', facet.by='treatmnt')

# homogeneity of variance
plot(cf_dir_clean1_aov, 1)
leveneTest(perc_in ~ treatmnt, data = cf_dir_clean1)

# check for skewness and kurtosis
describe(cf_dir_clean1$perc_in, na.rm = T, type = 2)

# kruskal wallis for non-normally distributed data

kruskal.test(perc_in ~ treatmnt, data = cf_dir_clean1) #signif

# post hoc test
# edit for both experiments
cf_dir_pwt1 <- pairwise.wilcox.test(cf_dir_clean1$perc_in, 
                                    cf_dir_clean1$treatmnt,
                                    p.adjust.method = "BH")
cf_dir_pwt_table1 <- fullPTable(cf_dir_pwt1$p.value) # table of p-values

cf_dir_letters1 <- multcompLetters(cf_dir_pwt_table1) # sig. dif. letters

# summarize the data to plot ----------------------------------------------
# edit for both experiments
cf_dir_treatmnt_summry1 <- cf_dir_clean1 |>  
  group_by(treatmnt) |>  
  summarise(mean_perc_in = mean(perc_in),
            sd_perc_in=sd(perc_in), 
            n_perc_in = n(),
            se_perc_in = sd(perc_in)/sqrt(n()))

# convert cld to list and merge with summarized data ----------------------
# edit for both exp
cf_dir_cld1 <- as.data.frame.list(cf_dir_letters1) # converts to list
cf_dir_treatmnt_summry1$wilcox_cl <- cf_dir_cld1$Letters

# save the summarized and merged data -------------------------------------
# edit for both exp
cf_dir_treatmnt_summry1 |> 
  arrange(desc(mean_perc_in)) |>  
  write_csv("output/cf_poy_treatmnt_summry1.csv")

# plotting percent pathogen on yeast of Bc by yeasts ----------------------

# edit for both experiments
cf_dir_plot1 <- cf_dir_treatmnt_summry1 |> 
  arrange(desc(mean_perc_in)) |> 
  ggplot(aes(factor(x=treatmnt, 
                    levels=c('CTRL','JB951','JB524','JB543','JB767')), 
             y=mean_perc_in)) + 
  geom_errorbar(aes(ymin = mean_perc_in - se_perc_in, 
                    ymax = mean_perc_in + se_perc_in),
                position = position_dodge(0.9), width=0.25,
                show.legend=FALSE) +
  geom_line(aes(group=1)) +
  geom_point() +
  geom_text(aes(label=wilcox_cl), position=position_dodge(0.9), size=4, 
            vjust=-2.3, hjust=0.3, color="#000") +
  scale_y_continuous(breaks=c(0,20,40,60,80,100),
                     limits=c(0,115),expand = c(0,0)) +
  theme_classic() +
  theme(
    axis.title=element_text(size=12, colour="#000", family='sans'),
    axis.title.x=element_blank(),
    axis.text=element_text(size=11, colour="#000", family='sans'),
    axis.ticks=element_line(color="#000"),
    panel.grid.major.x=element_line(),
    panel.grid.major.y=element_line()
  ) +
  labs(y='MRGI (%)', x='') 

ggsave("output/cf_poy_plot1.tif",cf_dir_plot1,height=3,width=3,dpi=300)
ggsave("output/cf_poy_plot1.pdf",cf_dir_plot1,height=3,width=3,dpi=300)


# pooled exp1 and 2 data for Colletotrichum -------------------------------

cf_dir_clean <- bccf_poy_clean |> 
  filter(pathogen=='Cf') |> #summarise(max_diam=max(diam_mm))
  mutate(max_diam=max(diam_mm)) |> 
  group_by(treatmnt) |> 
  mutate(perc_in=100*(max_diam-diam_mm)/max_diam,
         log_diam=log10(diam_mm + 1),
         sqrt_diam=sqrt(diam_mm)) |> 
  ungroup()

# check variables and make appropriate conversions ------------------------
str(cf_dir_clean)
cf_dir_clean$treatmnt <- as.factor(cf_dir_clean$treatmnt)
glimpse(cf_dir_clean)  

# perform eda -------------------------------------------------------------
# summary statistics
cf_dir_clean |> 
  group_by(treatmnt) |> 
  get_summary_stats(perc_in, type='mean_sd')

# boxplot
cf_dir_clean |> 
  ggboxplot(x='treatmnt', y='perc_in') # we see 2 outliers in the CTRL

# histogram
cf_dir_clean |> 
  gghistogram(x='perc_in')

# identify outliers
cf_dir_clean |> 
  group_by(treatmnt) |> 
  identify_outliers(perc_in)

# check for normality -----------------------------------------------------

cf_dir_clean_aov <- aov(perc_in ~ treatmnt, data=cf_dir_clean)
summary(cf_dir_clean_aov)

# qqplot
ggqqplot(residuals(cf_dir_clean_aov)) # one way
plot(cf_dir_clean_aov, 2) # second way

# shapiro test
shapiro_test(residuals(cf_dir_clean_aov)) # p is signif

cf_dir_clean_aov_residuals <- residuals(object=cf_dir_clean_aov)
shapiro.test(x = cf_dir_clean_aov_residuals)

# check shapiro for each treatment level
cf_dir_clean |>
  group_by(treatmnt) |> 
  shapiro_test(perc_in) 

# qqplot for each group level 
ggqqplot(cf_dir_clean, 'perc_in', facet.by='treatmnt')

#homogeneity of variance
plot(cf_dir_clean_aov, 1)
leveneTest(perc_in ~ treatmnt, data = cf_dir_clean)

#check for skewness and kurtosis
describe(cf_dir_clean$perc_in, na.rm = T, type = 2)

#kruskal wallis for non-normally distributed data

kruskal.test(perc_in ~ treatmnt, data = cf_dir_clean) #signif


#post hoc test
cf_dir_pwt <- dunnTest(cf_dir_clean$perc_in, 
                       cf_dir_clean$treatmnt, method = "bh")
cf_dir_cld <- cldList(P.adj ~ Comparison, data=cf_dir_pwt$res)

# summarize the data to plot ----------------------------------------------
# experiment 1
cf_dir_treatmnt_summry <- cf_dir_clean |>  
  group_by(treatmnt) |>  
  summarise(mean_perc_in = mean(perc_in),
            sd_perc_in=sd(perc_in), 
            n_perc_in = n(),
            se_perc_in = sd(perc_in)/sqrt(n()))

# convert cld to list and merge with summarized data ----------------------
cf_dir_treatmnt_summry$dunn_cl <- cf_dir_cld$Letter 

# save the summarized and merged data -------------------------------------
cf_dir_treatmnt_summry |> 
  arrange(desc(mean_perc_in)) |>  
  write_csv("output/cf_poy_treatmnt_summry.csv")

# plotting pooled percent pathogen on yeast of Cf by yeasts ---------------

cf_poy_plot <- cf_dir_treatmnt_summry |> 
  arrange(desc(mean_perc_in)) |> 
  ggplot(aes(factor(x=treatmnt, 
                    levels=c('CTRL','JB951','JB543','JB767','JB624')), 
             y=mean_perc_in)) + 
  geom_errorbar(aes(ymin = mean_perc_in - se_perc_in, 
                    ymax = mean_perc_in + se_perc_in),
                position = position_dodge(0.9), width=0.25,
                show.legend=FALSE) +
  geom_line(aes(group=1)) +
  geom_point() +
  geom_text(aes(label=dunn_cl), position=position_dodge(0.9), size=4, 
            vjust=-2.3, hjust=0.3, color="#000") +
  scale_y_continuous(breaks=c(0,20,40,60,80,100),
                     limits=c(0,115),expand = c(0,0)) +
  theme_classic() +
  theme(
    axis.title=element_text(size=12, colour="#000", family='sans'),
    axis.title.x=element_blank(),
    axis.text=element_text(size=11, colour="#000", family='sans'),
    axis.ticks=element_line(color="#000"),
    panel.grid.major.x=element_line(),
    panel.grid.major.y=element_line()
  ) +
  labs(y='MRGI (%)', x='') 

ggsave("output/cf_poy_plot_pooled.tif",cf_poy_plot,height=4,width=4,dpi=300)
ggsave("output/cf_poy_plot_pooled.pdf",cf_poy_plot,height=4,width=4,dpi=300)


# Volatile organic compounds (VOCs) ---------------------------------------

# import, clean, and export data ------------------------------------------
bccf_voc_data_clean <- read_excel("data/all_data.xlsx", sheet=3) |> 
  mutate(treatmnt=case_when(isolate_id=='C'~'CTRL',
                            isolate_id=='357'~'JB543',
                            isolate_id=='461'~'JB524',
                            isolate_id=='506'~'JB624',
                            isolate_id=='707'~'JB767',
                            isolate_id=='957'~'JB951'),
         diam_mm=diam_mm-5) |> # subtract initial 5mm disc used
  select(-c(date_inoculation,date_data_col,isolate_id)) |> 
  relocate(treatmnt, .before=pathogen) #|> #changes position of isolate
  #write.xlsx('output/bccf_voc_clean_exp1-2.xlsx')


################################# Botrytis #####################################

# Subset Bc data for analysis ---------------------------------------------
bc_voc_clean1 <- bccf_voc_data_clean |> 
  filter(pathogen=='Bc') |>
  filter(experiment==1) |> 
  mutate(max_diam=max(diam_mm)) |> 
  group_by(treatmnt) |> 
  mutate(perc_in=100*(max_diam-diam_mm)/max_diam,
         log_diam=log10(diam_mm + 1),
         sqrt_diam=sqrt(diam_mm)) |>
    ungroup()

bc_voc_clean2 <- bccf_voc_data_clean |> 
  filter(pathogen=='Bc') |> 
  filter(experiment==2) |> 
  mutate(max_diam=max(diam_mm)) |>
  group_by(treatmnt) |> 
  mutate(perc_in=100*(max_diam-diam_mm)/max_diam,
         log_diam=log10(diam_mm + 1),
         sqrt_diam=sqrt(diam_mm)) |>
  ungroup()

# check variables and make appropriate conversions ------------------------
# edit for both exp
str(bc_voc_clean1)
bc_voc_clean1$treatmnt <- as.factor(bc_voc_clean1$treatmnt)
glimpse(bc_voc_clean1)  

# perform eda -------------------------------------------------------------
# edit for both exp

# summary statistics
bc_voc_clean1 |> 
  group_by(treatmnt) |> 
  get_summary_stats(perc_in, type='mean_sd')

# boxplot
bc_voc_clean1 %>% 
  ggboxplot(., x='treatmnt', y='perc_in') 

# histogram
bc_voc_clean1 %>%
  gghistogram(., x='perc_in')

# identify outliers
bc_voc_clean1 |> 
  group_by(treatmnt) |> 
  identify_outliers(perc_in)

# check for normality -----------------------------------------------------

bc_voc_clean1_aov <- aov(perc_in ~ treatmnt, data=bc_voc_clean1)
summary(bc_voc_clean1_aov)

# qqplot
ggqqplot(residuals(bc_voc_clean1_aov)) # one way
plot(bc_voc_clean1_aov, 2) # second way

# shapiro test
shapiro_test(residuals(bc_voc_clean1_aov)) # p is signif

bc_voc_clean1_aov_residuals <- residuals(object=bc_voc_clean1_aov)
shapiro.test(x = bc_voc_clean1_aov_residuals)

# check shapiro for each treatment level
bc_voc_clean1 |>
  group_by(treatmnt) |> 
  shapiro_test(perc_in) # p here is also signif for all levels

# qqplot for each group level 
ggqqplot(bc_voc_clean1, 'perc_in', facet.by='treatmnt')

# homogeneity of variance
plot(bc_voc_clean1_aov, 1)
leveneTest(perc_in ~ treatmnt, data = bc_voc_clean1)

# check for skewness and kurtosis
describe(bc_voc_clean1$perc_in, na.rm = T, type = 2)

# kruskal wallis for non-normally distributed data
kruskal.test(perc_in ~ treatmnt, data = bc_voc_clean1) #signif

# post hoc test
# edit for both exp
bc_voc_pwt1 <- dunnTest(bc_voc_clean1$perc_in, 
                        bc_voc_clean1$treatmnt, method = "bh")

bc_voc_cld1 <- cldList(P.adj ~ Comparison, data=bc_voc_pwt1$res)

# summarize the data to plot ----------------------------------------------
# edit for both exp
bc_voc_treatmnt_summry1 <- bc_voc_clean1 |>  
  group_by(treatmnt) |>  
  summarise(mean_perc_in = mean(perc_in),
            sd_perc_in=sd(perc_in), 
            n_perc_in = n(),
            se_perc_in = sd(perc_in)/sqrt(n()))

# convert cld to list and merge with summarized data ----------------------
# edit for both exp
bc_voc_treatmnt_summry1$dunn_cl <- bc_voc_cld1$Letter

# save the summarized and merged data -------------------------------------
# edit for both exp
bc_voc_treatmnt_summry1 |> 
  arrange(desc(mean_perc_in)) |> 
  write_csv("output/bc_voc_treatmnt_summry1.csv")

# Plotting percent growth inhibition of Bc --------------------------------
# edit for both exp

bc_voc_plot1 <- bc_voc_treatmnt_summry1 |> 
  arrange(desc(mean_perc_in)) |> 
  ggplot(aes(factor(x=treatmnt,levels=c('CTRL','JB524','JB767','JB543','JB951')), 
             y=mean_perc_in,
             ymin = mean_perc_in - se_perc_in, 
             ymax = mean_perc_in + se_perc_in)) + 
  geom_linerange() +
  geom_line(aes(group=1)) +
  geom_point(size=1) +
  geom_text(aes(label=dunn_cl), position=position_dodge(0.9), size=4, 
            vjust=-2.3, hjust=0.3, color="#000") +
  scale_y_continuous(breaks=c(0,25,50,75,100),
                     limits=c(0,101),expand = c(0,0)) +
  theme_bw() +
  theme(
    axis.title=element_text(size=12, colour="#000", family='sans'),
    axis.title.x=element_blank(),
    axis.text=element_text(size=11, colour="#000", family='sans'),
    axis.ticks.x=element_blank(),
    axis.ticks.y=element_line(color="#000")
  ) +
  labs(y='MGI (%)', x='') 

ggsave("output/bc_voc_plot1.tif",bc_voc_plot1,height=4,width=4,dpi=300)
ggsave("output/bc_voc_plot1.pdf",bc_voc_plot1,height=4,width=4,dpi=300)


# pooled exp1 and 2 data for Botrytis -------------------------------------

bc_voc_clean <- bccf_voc_data_clean |> 
  filter(pathogen=='Bc') |> 
  mutate(max_diam=max(diam_mm)) |>
  group_by(treatmnt) |> 
  mutate(perc_in=100*(max_diam-diam_mm)/max_diam,
         log_diam=log10(diam_mm + 1),
         sqrt_diam=sqrt(diam_mm)) |>
  ungroup()

# check variables and make appropriate conversions ------------------------
str(bc_voc_clean)
bc_voc_clean$treatmnt <- as.factor(bc_voc_clean$treatmnt)
glimpse(bc_voc_clean)  

# perform eda -------------------------------------------------------------
# summary statistics
bc_voc_clean |> 
  group_by(treatmnt) |> 
  get_summary_stats(perc_in, type='mean_sd')

# boxplot
bc_voc_clean %>% 
  ggboxplot(., x='treatmnt', y='perc_in') 

# histogram
bc_voc_clean %>%
  gghistogram(., x='perc_in')

# identify outliers
bc_voc_clean |> 
  group_by(treatmnt) |> 
  identify_outliers(perc_in)

# check for normality -----------------------------------------------------

bc_voc_clean_aov <- aov(perc_in ~ treatmnt, data=bc_voc_clean)
summary(bc_voc_clean_aov)

# qqplot
ggqqplot(residuals(bc_voc_clean_aov)) # one way
plot(bc_voc_clean_aov, 2) # second way

# shapiro test
shapiro_test(residuals(bc_voc_clean_aov)) # p is signif

bc_voc_clean_aov_residuals <- residuals(object=bc_voc_clean_aov)
shapiro.test(x = bc_voc_clean_aov_residuals)

# check shapiro for each treatment level
bc_voc_clean |>
  group_by(treatmnt) |> 
  shapiro_test(perc_in) # p here is also signif for all levels exp JB951

# qqplot for each group level 
ggqqplot(bc_voc_clean, 'perc_in', facet.by='treatmnt')

#homogeneity of variance
plot(bc_voc_clean_aov, 1)
leveneTest(perc_in ~ treatmnt, data = bc_voc_clean)

#check for skewness and kurtosis
describe(bc_voc_clean$perc_in, na.rm = T, type = 2)

#kruskal wallis for non-normally distributed data
kruskal.test(perc_in ~ treatmnt, data = bc_voc_clean) #signif

#post hoc test
bc_voc_pwt <- dunnTest(bc_voc_clean$perc_in, 
                        bc_voc_clean$treatmnt, method = "bh")

bc_voc_cld <- cldList(P.adj ~ Comparison, data=bc_voc_pwt$res)

# summarize the data to plot ----------------------------------------------
bc_voc_treatmnt_summry <- bc_voc_clean |>  
  group_by(treatmnt) |>  
  summarise(mean_perc_in = mean(perc_in),
            sd_perc_in=sd(perc_in), 
            n_perc_in = n(),
            se_perc_in = sd(perc_in)/sqrt(n()))

# convert cld to list and merge with summarized data ----------------------
bc_voc_treatmnt_summry$dunn_cl <- bc_voc_cld$Letter

# save the summarized and merged data -------------------------------------
bc_voc_treatmnt_summry |> 
  arrange(desc(mean_perc_in)) |> 
  write.xlsx("output/bc_voc_treatmnt_summry.xlsx")


# plotting percent voc growth inhibition of Bc ----------------------------

bc_voc_plot <- bc_voc_treatmnt_summry |> 
  arrange(desc(mean_perc_in)) |> 
  ggplot(aes(factor(x=treatmnt,levels=c('CTRL','JB524','JB767','JB543','JB951')), 
             y=mean_perc_in,
             ymin = mean_perc_in - se_perc_in, 
             ymax = mean_perc_in + se_perc_in)) + 
  geom_linerange() +
  geom_line(aes(group=1)) +
  geom_point(size=1) +
  geom_text(aes(label=dunn_cl), position=position_dodge(0.9), size=4, 
            vjust=-2.3, hjust=0.3, color="#000") +
  scale_y_continuous(breaks=c(0,25,50,75,100),
                     limits=c(0,101),expand = c(0,0)) +
  theme_bw() +
  theme(
    axis.title=element_text(size=12, colour="#000", family='sans'),
    axis.title.x=element_blank(),
    axis.text=element_text(size=11, colour="#000", family='sans'),
    axis.ticks.x=element_blank(),
    axis.ticks.y=element_line(color="#000")
  ) +
  labs(y='MGI (%)', x='') 

ggsave("output/bc_voc_plot_pooled.tif",bc_voc_plot,height=4,width=4,dpi=300)
ggsave("output/bc_voc_plot_pooled.pdf",bc_voc_plot,height=4,width=4,dpi=300)


########################### Colletotrichum ################################


# Subset Cf data for analysis ---------------------------------------------
cf_voc_clean1 <- bccf_voc_data_clean |> 
  filter(pathogen=='Cf') |>
  filter(experiment==1) |> 
  mutate(max_diam=max(diam_mm)) |> 
  group_by(treatmnt) |> 
  mutate(perc_in=100*(max_diam-diam_mm)/max_diam,
         log_diam=log10(diam_mm + 1),
         sqrt_diam=sqrt(diam_mm)) |>
  ungroup()

cf_voc_clean2 <- bccf_voc_data_clean |> 
  filter(pathogen=='Cf') |> 
  filter(experiment==2) |> 
  mutate(max_diam=max(diam_mm)) |>
  group_by(treatmnt) |> 
  mutate(perc_in=100*(max_diam-diam_mm)/max_diam,
         log_diam=log10(diam_mm + 1),
         sqrt_diam=sqrt(diam_mm)) |>
  ungroup()

# check variables and make appropriate conversions ------------------------
# edit for both exp
str(cf_voc_clean1)
cf_voc_clean1$treatmnt <- as.factor(cf_voc_clean1$treatmnt)
glimpse(cf_voc_clean1)  

# perform eda -------------------------------------------------------------
# edit for both exp

# summary statistics
cf_voc_clean1 |> 
  group_by(treatmnt) |> 
  get_summary_stats(perc_in, type='mean_sd')

# boxplot
cf_voc_clean1 %>% 
  ggboxplot(., x='treatmnt', y='perc_in') 

# histogram
cf_voc_clean1 %>%
  gghistogram(., x='perc_in')

# identify outliers
cf_voc_clean1 |> 
  group_by(treatmnt) |> 
  identify_outliers(perc_in)

# check for normality -----------------------------------------------------
cf_voc_clean1_aov <- aov(perc_in ~ treatmnt, data=cf_voc_clean1)
summary(cf_voc_clean1_aov)

# qqplot
ggqqplot(residuals(cf_voc_clean1_aov)) # one way
plot(cf_voc_clean1_aov, 2) # second way

# shapiro test
shapiro_test(residuals(cf_voc_clean1_aov)) # p is signif

cf_voc_clean1_aov_residuals <- residuals(object=cf_voc_clean1_aov)
shapiro.test(x = cf_voc_clean1_aov_residuals)

# check shapiro for each treatment level
cf_voc_clean1 |>
  group_by(treatmnt) |> 
  shapiro_test(perc_in) # p here is also signif for all levels

# qqplot for each group level 
ggqqplot(cf_voc_clean1, 'perc_in', facet.by='treatmnt')

# homogeneity of variance
plot(cf_voc_clean1_aov, 1)
leveneTest(perc_in ~ treatmnt, data = cf_voc_clean1)

# check for skewness and kurtosis
describe(cf_voc_clean1$perc_in, na.rm = T, type = 2)

# kruskal wallis for non-normally distributed data
kruskal.test(perc_in ~ treatmnt, data = cf_voc_clean1) #signif

# post hoc test
# edit for both exp
cf_voc_pwt1 <- dunnTest(cf_voc_clean1$perc_in, 
                        cf_voc_clean1$treatmnt, method = "bh")

cf_voc_cld1 <- cldList(P.adj ~ Comparison, data=cf_voc_pwt1$res)

# summarize the data to plot ----------------------------------------------
# edit for both exp
cf_voc_treatmnt_summry1 <- cf_voc_clean1 |>  
  group_by(treatmnt) |>  
  summarise(mean_perc_in = mean(perc_in),
            sd_perc_in=sd(perc_in), 
            n_perc_in = n(),
            se_perc_in = sd(perc_in)/sqrt(n()))

# convert cld to list and merge with summarized data ----------------------
# edit for both exp
cf_voc_treatmnt_summry1$dunn_cl <- cf_voc_cld1$Letter

# save the summarized and merged data -------------------------------------
# edit for both exp
cf_voc_treatmnt_summry1 |> 
  arrange(desc(mean_perc_in)) |> 
  write_csv("output/cf_voc_treatmnt_summry1.csv")

# plotting percent voc growth inhibition of Cf ----------------------------

# edit for both exp
cf_voc_plot1 <- cf_voc_treatmnt_summry1 |> 
  arrange(desc(mean_perc_in)) |> 
  ggplot(aes(factor(x=treatmnt,levels=c('CTRL','JB543','JB767','JB624','JB951')), 
             y=mean_perc_in,
             ymin = mean_perc_in - se_perc_in, 
             ymax = mean_perc_in + se_perc_in)) + 
  geom_linerange() +
  geom_line(aes(group=1)) +
  geom_point(size=1) +
  geom_text(aes(label=dunn_cl), position=position_dodge(0.9), size=4, 
            vjust=-2.3, hjust=0.3, color="#000") +
  scale_y_continuous(breaks=c(0,25,50,75,100),
                     limits=c(0,101),expand = c(0,0)) +
  theme_bw() +
  theme(
    axis.title=element_text(size=12, colour="#000", family='sans'),
    axis.title.x=element_blank(),
    axis.text=element_text(size=11, colour="#000", family='sans'),
    axis.ticks.x=element_blank(),
    axis.ticks.y=element_line(color="#000")
  ) +
  labs(y='MGI (%)', x='') 

ggsave("output/cf_voc_plot1.tif",cf_voc_plot1,height=4,width=4,dpi=300)
ggsave("output/cf_voc_plot1.pdf",cf_voc_plot1,height=4,width=4,dpi=300)

# pooled exp1 and 2 data for Colletotrichum -------------------------------

cf_voc_clean <- bccf_voc_data_clean |> 
  filter(pathogen=='Cf') |> 
  mutate(max_diam=max(diam_mm)) |>
  group_by(treatmnt) |> 
  mutate(perc_in=100*(max_diam-diam_mm)/max_diam,
         log_diam=log10(diam_mm + 1),
         sqrt_diam=sqrt(diam_mm)) |>
  ungroup()

# check variables and make appropriate conversions ------------------------
str(cf_voc_clean)
cf_voc_clean$treatmnt <- as.factor(cf_voc_clean$treatmnt)
glimpse(cf_voc_clean)  

# perform eda -------------------------------------------------------------
# summary statistics
cf_voc_clean |> 
  group_by(treatmnt) |> 
  get_summary_stats(perc_in, type='mean_sd')

# boxplot
cf_voc_clean %>% 
  ggboxplot(., x='treatmnt', y='perc_in') 

# histogram
cf_voc_clean %>%
  gghistogram(., x='perc_in')

# identify outliers
cf_voc_clean |> 
  group_by(treatmnt) |> 
  identify_outliers(perc_in) # no extreme values identified

# check for normality -----------------------------------------------------

cf_voc_clean_aov <- aov(perc_in ~ treatmnt, data=cf_voc_clean)
summary(cf_voc_clean_aov)

# qqplot
ggqqplot(residuals(cf_voc_clean_aov)) # one way
plot(cf_voc_clean_aov, 2) # second way

# shapiro test
shapiro_test(residuals(cf_voc_clean_aov)) # p is signif

cf_voc_clean_aov_residuals <- residuals(object=cf_voc_clean_aov)
shapiro.test(x = cf_voc_clean_aov_residuals)

# check shapiro for each treatment level
cf_voc_clean |>
  group_by(treatmnt) |> 
  shapiro_test(perc_in) # p here is also signif for all levels exp JB951

# qqplot for each group level 
ggqqplot(cf_voc_clean, 'perc_in', facet.by='treatmnt')

#homogeneity of variance
plot(cf_voc_clean_aov, 1)
leveneTest(perc_in ~ treatmnt, data = cf_voc_clean)

#check for skewness and kurtosis
describe(cf_voc_clean$perc_in, na.rm = T, type = 2)

#kruskal wallis for non-normally distributed data
kruskal.test(perc_in ~ treatmnt, data = cf_voc_clean) #signif

#post hoc test
cf_voc_pwt <- dunnTest(cf_voc_clean$perc_in, 
                       cf_voc_clean$treatmnt, method = "bh") #'bonferroni'

cf_voc_cld <- cldList(P.adj ~ Comparison, data=cf_voc_pwt$res)

# summarize the data to plot ----------------------------------------------
cf_voc_treatmnt_summry <- cf_voc_clean |>  
  group_by(treatmnt) |>  
  summarise(mean_perc_in = mean(perc_in),
            sd_perc_in=sd(perc_in), 
            n_perc_in = n(),
            se_perc_in = sd(perc_in)/sqrt(n()))

# convert cld to list and merge with summarized data ----------------------
cf_voc_treatmnt_summry$dunn_cl <- cf_voc_cld$Letter

# save the summarized and merged data -------------------------------------
cf_voc_treatmnt_summry |> 
  arrange(desc(mean_perc_in)) |> 
  write.xlsx("output/cf_voc_treatmnt_summry.xlsx")


# Plotting pooled percent voc growth inhibition of Cf ---------------------

cf_voc_plot <- cf_voc_treatmnt_summry |> 
  arrange(desc(mean_perc_in)) |> 
  ggplot(aes(factor(x=treatmnt,levels=c('CTRL','JB543','JB767','JB624','JB951')), 
             y=mean_perc_in,
             ymin = mean_perc_in - se_perc_in, 
             ymax = mean_perc_in + se_perc_in)) + 
  geom_linerange() +
  geom_line(aes(group=1)) +
  geom_point(size=1) +
  geom_text(aes(label=dunn_cl), position=position_dodge(0.9), size=4, 
            vjust=-2.3, hjust=0.3, color="#000") +
  scale_y_continuous(breaks=c(0,25,50,75,100),
                     limits=c(0,101),expand = c(0,0)) +
  theme_bw() +
  theme(
    axis.title=element_text(size=12, colour="#000", family='sans'),
    axis.title.x=element_blank(),
    axis.text=element_text(size=11, colour="#000", family='sans'),
    axis.ticks.x=element_blank(),
    axis.ticks.y=element_line(color="#000")
  ) +
  labs(y='MGI (%)', x='') 

ggsave("output/cf_voc_plot_pooled.tif",cf_voc_plot,height=4,width=4,dpi=300)
ggsave("output/cf_voc_plot_pooled.pdf",cf_voc_plot,height=4,width=4,dpi=300)


# Correlation Test --------------------------------------------------------

library(ggpubr)

# Inhibition zone vs percent inhibition -----------------------------------
cor_iz_pi_plot <- read_xlsx('output/bccf_antagon_clean_exp1-3.xlsx') |> 
  select(pathogen,ydiam_mm,iz,perc_in) |> 
  ggscatter(x="iz", y="perc_in", 
            add="reg.line", conf.int=TRUE, 
            cor.coef=TRUE, cor.method="pearson",
            xlab="Inhibition zone (mm)", ylab="Percent inhibition (%)") +
  theme_bw() +
  theme(
    axis.title=element_text(size=12, colour="#000", family='sans'),
    axis.text=element_text(size=11, colour="#000", family='sans'),
    axis.ticks=element_line(color="#000")
  )

ggsave("output/cor_iz_pi_sorted.tif",cor_iz_pi_plot,height=4,width=4,dpi=300)
ggsave("output/cor_iz_pi_sorted.pdf",cor_iz_pi_plot,height=5,width=6,dpi=300)

# Inhibition zone vs yeast diameter ---------------------------------------
cor_iz_yd_plot <- read_xlsx('output/bccf_antagon_clean_exp1-3.xlsx') |> 
  select(pathogen,ydiam_mm,iz,perc_in) |>  
  filter(iz>0) |>  
  ggscatter(x="iz", y="ydiam_mm", 
            add="reg.line", conf.int=TRUE, 
            cor.coef=TRUE, cor.method="pearson",
            xlab="Inhibition zone (mm)", ylab="Yeast diameter (mm)") +
  theme_bw() +
  theme(
    axis.title=element_text(size=12, colour="#000", family='sans'),
    axis.text=element_text(size=11, colour="#000", family='sans'),
    axis.ticks=element_line(color="#000")
  )

ggsave("output/cor_iz_yd_sorted.tif",cor_iz_yd_plot,height=4,width=4,dpi=300)
ggsave("output/cor_iz_yd_sorted.pdf",cor_iz_yd_plot,height=4,width=4,dpi=300)


# Spore Germination Data Analysis -----------------------------------------

# import and clean data for analysis --------------------------------------

spore_germ <- read_xlsx('data/all_data.xlsx', sheet = 4) |> 
  clean_names() |> 
  select(-c(reps,7:8)) |> 
  pivot_longer(-c(1:2,6), names_to = 'treatmnt', values_to = 'value') |> 
  arrange(treatmnt) #|> view()
  
spore_germ$experiment <- as.factor(spore_germ$experiment)
#spore_germ$treatmnt <- as.factor(spore_germ$treatmnt)
glimpse(spore_germ)  
  
# Botrytis

# edit for both exp
spore_germ_bc1 <- spore_germ |> 
  filter(pathogen=='Bc') |> 
  filter(experiment=='1') |> 
  group_by(yeast_id,treatmnt) |> 
  summarise(N=n(),
         perc_germ = (100*(sum(value)/N))+.5, .groups = 'drop') |> 
  select(-N)

# plotting
bc_sg_plot1 <- spore_germ_bc1 %>% 
  mutate(treatmnt=fct_reorder(treatmnt,desc(perc_germ))) %>% 
  group_by(yeast_id, treatmnt) %>%
  ggbarplot(., x = "yeast_id", y = "perc_germ", add = "mean_se",
            fill = "treatmnt",  
            position = position_dodge(0.8)) +
  stat_compare_means(label.y = 82.5, size = 3) +
  stat_compare_means(aes(label = after_stat(p.signif)),
                     method = "wilcox.test", 
                     ref.group = "sdw",
                     label.y = 20.5,
                     size = 4) +
  scale_y_continuous(limits = c(0,85),
                     expand = c(0,0)) +
  scale_fill_manual(name=NULL,
                    breaks=c("sdw","wcs","cfs"),
                    labels=c("SDW","WCS","CFS"),
                    values=c("#483D8B", "#1E90FF","#F0F0F0")) +
  labs(x="",
       y="Spore germination (%)")

ggsave("output/bc_sg_plot1.tif",bc_sg_plot1,dpi=300,width=4,height=4)

#post hoc test
bc_sg_pwt1 <- pairwise.wilcox.test(spore_germ_bc1$perc_germ, 
                                           spore_germ_bc1$treatmnt,
                                        p.adjust.method = "BH")

bc_sg_tab1 <- fullPTable(bc_sg_pwt1$p.value) # tab of p-values
bc_sg_let1 <- multcompLetters(bc_sg_tab1)
bc_sg_cld1 <- as.data.frame.list(bc_sg_let1)

# pooled spore germination data for B. cinerea ----------------------------

# for 767
sg_bc_767_pooled <- spore_germ |> 
  filter(pathogen=='Bc') |>
  filter(yeast_id=='JB767') |> 
  select(-c(experiment,pathogen)) 

#kruskal wallis for non-normally distributed data
kruskal.test(value ~ treatmnt, data = sg_bc_767_pooled) #signif

# post hoc test
# pooled data
sg_bc_767_pwt <- dunnTest(sg_bc_767_pooled$value, 
                          sg_bc_767_pooled$treatmnt, method = "bh")

sg_bc_767_cld <- cldList(P.adj ~ Comparison, data=sg_bc_767_pwt$res)

# summarize the data to plot ----------------------------------------------
sg_bc_767_summry <- sg_bc_767_pooled |>  
  group_by(treatmnt) |>  
  summarise(mean_sg = mean(value),
            sd_sg=sd(value), 
            n_sg = n(),
            se_sg = 1 + sd(value)/sqrt(n()),
            perc_sg = (100*(sum(value)/n_sg))+ 1, .groups = 'drop')

# convert cld to list and merge with summarized data ----------------------
sg_bc_767_summry$dunn_cl <- sg_bc_767_cld$Letter

# save the summarized and merged data -------------------------------------
sg_bc_767_summry |> 
  arrange(desc(perc_sg)) |> 
  write.xlsx("output/sg_bc_767_summry.xlsx")

# Plotting percent germination rate of bc ---------------------------------
sg_bc_767_pooled <- sg_bc_767_summry |> 
  arrange(desc(perc_sg)) |> 
  ggplot(aes(factor(x=treatmnt,levels=c('sdw','wcs','cfs')), 
             y=perc_sg, fill=treatmnt)) + 
  geom_errorbar(aes(ymin = perc_sg - se_sg, 
                    ymax = perc_sg + se_sg),
                position = position_dodge(0.9), width=0.25,
                show.legend=FALSE) +
  geom_bar(stat="identity", position="dodge", alpha=.9,
           color="#000", width=0.5, show.legend=F) +
  geom_text(aes(label=dunn_cl), position=position_dodge(0.9), size=4, 
            vjust=-2.3, hjust=0.3, color="#000") +
  scale_y_continuous(breaks=c(0,20,40,60,80,100),
                     limits=c(0,100),expand = c(0,0)) +
  scale_fill_manual(values=c("#483D8B","#6CA6CD","#F0F0F0","#1E90FF")) +
  scale_x_discrete(name=NULL,
                   breaks=c("sdw","wcs","cfs"),
                   labels=c("SDW","WCS","CFS")) +
  theme_bw() +
  theme(
    axis.title=element_text(size=12, colour="#000", family='sans'),
    axis.title.x=element_blank(),
    axis.text=element_text(size=11, colour="#000", family='sans'),
    axis.ticks.x=element_blank(),
    axis.ticks.y=element_line(color="#000")
  ) +
  labs(y='SGR (%)', x='') 

ggsave("output/sg_bc_767_pooled.tif",sg_bc_767_pooled,height=4,width=4,dpi=300)
ggsave("output/sg_bc_767_pooled.pdf",sg_bc_767_pooled,height=4,width=4,dpi=300)

# for 951
sg_bc_951_pooled <- spore_germ |> 
  filter(pathogen=='Bc') |>
  filter(yeast_id=='JB951') |> 
  select(-c(experiment,pathogen)) 

#kruskal wallis for non-normally distributed data
kruskal.test(value ~ treatmnt, data = sg_bc_951_pooled) #signif

# post hoc test
sg_bc_951_pwt <- dunnTest(sg_bc_951_pooled$value, 
                          sg_bc_951_pooled$treatmnt, method = "bh")

sg_bc_951_cld <- cldList(P.adj ~ Comparison, data=sg_bc_951_pwt$res)

# summarize the data to plot ----------------------------------------------
sg_bc_951_summry <- sg_bc_951_pooled |>  
  group_by(treatmnt) |>  
  summarise(mean_sg = mean(value),
            sd_sg=sd(value), 
            n_sg = n(),
            se_sg = 1 + sd(value)/sqrt(n()),
            perc_sg = (100*(sum(value)/n_sg))+ 1, .groups = 'drop')

# convert cld to list and merge with summarized data ----------------------
sg_bc_951_summry$dunn_cl <- sg_bc_951_cld$Letter

# save the summarized and merged data -------------------------------------
sg_bc_951_summry |> 
  arrange(desc(perc_sg)) |> 
  write.xlsx("output/sg_bc_951_summry.xlsx")

# Plotting percent germination rate of bc ---------------------------------
sg_bc_951_pooled <- sg_bc_951_summry |> 
  arrange(desc(perc_sg)) |> 
  ggplot(aes(factor(x=treatmnt,levels=c('sdw','wcs','cfs')), 
             y=perc_sg, fill=treatmnt)) + 
  geom_errorbar(aes(ymin = perc_sg - se_sg, 
                    ymax = perc_sg + se_sg),
                position = position_dodge(0.9), width=0.25,
                show.legend=FALSE) +
  geom_bar(stat="identity", position="dodge", alpha=.9,
           color="#000", width=0.5, show.legend=F) +
  geom_text(aes(label=dunn_cl), position=position_dodge(0.9), size=4, 
            vjust=-2.3, hjust=0.3, color="#000") +
  scale_y_continuous(breaks=c(0,20,40,60,80,100),
                     limits=c(0,100),expand = c(0,0)) +
  scale_fill_manual(values=c("#483D8B","#6CA6CD","#F0F0F0","#1E90FF")) +
  scale_x_discrete(name=NULL,
                   breaks=c("sdw","wcs","cfs"),
                   labels=c("SDW","WCS","CFS")) +
  theme_bw() +
  theme(
    axis.title=element_text(size=12, colour="#000", family='sans'),
    axis.title.x=element_blank(),
    axis.text=element_text(size=11, colour="#000", family='sans'),
    axis.ticks.x=element_blank(),
    axis.ticks.y=element_line(color="#000")
  ) +
  labs(y='SGR (%)', x='') 

ggsave("output/sg_bc_951_pooled.tif",sg_bc_951_pooled,height=4,width=4,dpi=300)
ggsave("output/sg_bc_951_pooled.pdf",sg_bc_951_pooled,height=4,width=4,dpi=300)


# Colletotrichum

# edit for both exp
spore_germ_cf1 <- spore_germ |> 
  filter(pathogen=='Cf') |> 
  filter(experiment=='1') |> 
  group_by(yeast_id,treatmnt) |> 
  summarise(N=n(),
            perc_germ = (100*(sum(value)/N))+.5, .groups = 'drop') |> 
  select(-N)

# plotting
cf_sg_plot1 <- spore_germ_cf1 %>% 
  mutate(treatmnt=fct_reorder(treatmnt,desc(perc_germ))) %>% 
  group_by(yeast_id, treatmnt) %>%
  ggbarplot(., x = "yeast_id", y = "perc_germ", add = "mean_se",
            fill = "treatmnt",  
            position = position_dodge(0.8)) +
  stat_compare_means(label.y = 82.5, size = 3) +
  stat_compare_means(aes(label = after_stat(p.signif)),
                     method = "wilcox.test", 
                     ref.group = "sdw",
                     label.y = 20.5,
                     size = 4) +
  scale_y_continuous(limits = c(0,90),
                     expand = c(0,0)) +
  scale_fill_manual(name=NULL,
                    breaks=c("sdw","cfs","wcs"),
                    labels=c("SDW","CFS","WCS"),
                    values=c("#483D8B", "#1E90FF","#F0F0F0")) +
  labs(x="",
       y="Spore germination (%)")

ggsave("output/cf_sg_plot1.tif",cf_sg_plot1,dpi=300,width=4,height=4)

#post hoc test
cf_sg_pwt1 <- pairwise.wilcox.test(spore_germ_cf1$perc_germ, 
                                           spore_germ_cf1$treatmnt,
                                           p.adjust.method = "BH")

cf_sg_tab1 <- fullPTable(cf_sg_pwt1$p.value) # tab of p-values
cf_sg_let1 <- multcompLetters(cf_sg_tab1)
cf_sg_cld1 <- as.data.frame.list(cf_sg_let1)

# pooled spore germination data for C. fioriniae --------------------------

# for 767
sg_cf_767_pooled <- spore_germ |> 
  filter(pathogen=='Cf') |>
  filter(yeast_id=='JB767') |> 
  select(-c(experiment,pathogen)) 

#kruskal wallis for non-normally distributed data
kruskal.test(value ~ treatmnt, data = sg_cf_767_pooled) #signif

# post hoc test
sg_cf_767_pwt <- dunnTest(sg_cf_767_pooled$value, 
                          sg_cf_767_pooled$treatmnt, method = "bh")

sg_cf_767_cld <- cldList(P.adj ~ Comparison, data=sg_cf_767_pwt$res)

# summarize the data to plot ----------------------------------------------
sg_cf_767_summry <- sg_cf_767_pooled |>  
  group_by(treatmnt) |>  
  summarise(mean_sg = mean(value),
            sd_sg=sd(value), 
            n_sg = n(),
            se_sg = 1 + sd(value)/sqrt(n()),
            perc_sg = (100*(sum(value)/n_sg))+ 1, .groups = 'drop')

# convert cld to list and merge with summarized data ----------------------
sg_cf_767_summry$dunn_cl <- sg_cf_767_cld$Letter

# save the summarized and merged data -------------------------------------
sg_cf_767_summry |> 
  arrange(desc(perc_sg)) |> 
  write.xlsx("output/sg_cf_767_summry.xlsx")

# Plotting percent germination rate of Cf ---------------------------------
sg_cf_767_pooled <- sg_cf_767_summry |> 
  arrange(desc(perc_sg)) |> 
  ggplot(aes(factor(x=treatmnt,levels=c('sdw','cfs','wcs')), 
             y=perc_sg, fill=treatmnt)) + 
  geom_errorbar(aes(ymin = perc_sg - se_sg, 
                    ymax = perc_sg + se_sg),
                position = position_dodge(0.9), width=0.25,
                show.legend=FALSE) +
  geom_bar(stat="identity", position="dodge", alpha=.9,
           color="#000", width=0.5, show.legend=F) +
  geom_text(aes(label=dunn_cl), position=position_dodge(0.9), size=4, 
            vjust=-2.3, hjust=0.3, color="#000") +
  scale_y_continuous(breaks=c(0,20,40,60,80,100),
                     limits=c(0,100),expand = c(0,0)) +
  scale_fill_manual(values=c("#483D8B","#6CA6CD","#F0F0F0","#1E90FF")) +
  scale_x_discrete(name=NULL,
                   breaks=c("sdw","cfs","wcs"),
                   labels=c("SDW","CFS","WCS")) +
  theme_bw() +
  theme(
    axis.title=element_text(size=12, colour="#000", family='sans'),
    axis.title.x=element_blank(),
    axis.text=element_text(size=11, colour="#000", family='sans'),
    axis.ticks.x=element_blank(),
    axis.ticks.y=element_line(color="#000")
  ) +
  labs(y='SGR (%)', x='') 

ggsave("output/sg_cf_767_pooled.tif",sg_cf_767_pooled,height=4,width=4,dpi=300)
ggsave("output/sg_cf_767_pooled.pdf",sg_cf_767_pooled,height=4,width=4,dpi=300)

# for 951
sg_cf_951_pooled <- spore_germ |> 
  filter(pathogen=='Cf') |>
  filter(yeast_id=='JB951') |> 
  select(-c(experiment,pathogen)) 

#kruskal wallis for non-normally distributed data
kruskal.test(value ~ treatmnt, data = sg_cf_951_pooled) #signif

# post hoc test
sg_cf_951_pwt <- dunnTest(sg_cf_951_pooled$value, 
                          sg_cf_951_pooled$treatmnt, method = "bh")

sg_cf_951_cld <- cldList(P.adj ~ Comparison, data=sg_cf_951_pwt$res)

# summarize the data to plot ----------------------------------------------
sg_cf_951_summry <- sg_cf_951_pooled |>  
  group_by(treatmnt) |>  
  summarise(mean_sg = mean(value),
            sd_sg=sd(value), 
            n_sg = n(),
            se_sg = 1 + sd(value)/sqrt(n()),
            perc_sg = (100*(sum(value)/n_sg))+ 1, .groups = 'drop')

# convert cld to list and merge with summarized data ----------------------
sg_cf_951_summry$dunn_cl <- sg_cf_951_cld$Letter

# save the summarized and merged data -------------------------------------
sg_cf_951_summry |> 
  arrange(desc(perc_sg)) |> 
  write.xlsx("output/sg_cf_951_summry.xlsx")

# Plotting percent germination rate of Cf ---------------------------------
sg_cf_951_pooled <- sg_cf_951_summry |> 
  arrange(desc(perc_sg)) |> 
  ggplot(aes(factor(x=treatmnt,levels=c('sdw','cfs','wcs')), 
             y=perc_sg, fill=treatmnt)) + 
  geom_errorbar(aes(ymin = perc_sg - se_sg, 
                    ymax = perc_sg + se_sg),
                position = position_dodge(0.9), width=0.25,
                show.legend=FALSE) +
  geom_bar(stat="identity", position="dodge", alpha=.9,
           color="#000", width=0.5, show.legend=F) +
  geom_text(aes(label=dunn_cl), position=position_dodge(0.9), size=4, 
            vjust=-2.3, hjust=0.3, color="#000") +
  scale_y_continuous(breaks=c(0,20,40,60,80,100),
                     limits=c(0,100),expand = c(0,0)) +
  scale_fill_manual(values=c("#483D8B","#6CA6CD","#F0F0F0","#1E90FF")) +
  scale_x_discrete(name=NULL,
                   breaks=c("sdw","cfs","wcs"),
                   labels=c("SDW","CFS","WCS")) +
  theme_bw() +
  theme(
    axis.title=element_text(size=12, colour="#000", family='sans'),
    axis.title.x=element_blank(),
    axis.text=element_text(size=11, colour="#000", family='sans'),
    axis.ticks.x=element_blank(),
    axis.ticks.y=element_line(color="#000")
  ) +
  labs(y='SGR (%)', x='') 

ggsave("output/sg_cf_951_pooled.tif",sg_cf_951_pooled,height=4,width=4,dpi=300)
ggsave("output/sg_cf_951_pooled.pdf",sg_cf_951_pooled,height=4,width=4,dpi=300)



# Table 1. Total Number of Yeasts Screened --------------------------------

read_xlsx('data/all_data.xlsx', sheet=5) |> 
  clean_names() |> 
  filter(complete.cases(isolate)) |> 
  filter(organism !='bacteria') |> 
  filter(name !='Control') |> 
  distinct(isolate, .keep_all=T) |> 
  group_by(name) |> 
  summarise(n=n(), .groups='drop') |> 
  write.xlsx('output/unique_yeast_screened_table1.xlsx')


# Supplementary Tabe1 Summary ---------------------------------------------

# pre-screened yeasts on Botrytis -----------------------------------------
read_xlsx('data/all_data.xlsx', sheet=6) |> 
  clean_names() |>
  mutate(diam_mm = radius_mm * 2) |> #distinct(organism)
  filter(organism !='bacteria') |> 
  select(-c(organism,radius_mm,ydiam,pattern,iz,dpi,ino_date,rec_date)) |> 
  filter(pathogen=='Bc') |> 
  group_by(sample,isolate,name,pathogen) |> 
  summarise(mean_diam = mean(diam_mm),
            sd_diam=sd(diam_mm), 
            N = n(),
            se_diam = sd(diam_mm)/sqrt(n()), .groups='drop') |> 
  arrange(name) |> 
  write.xlsx('output/bc_prescreened_yeast_sup_tab1.xlsx')

# pre-screened yeasts on Colletotrichum -------------------------------------
read_xlsx('data/all_data.xlsx', sheet=6) |> 
  clean_names() |>
  mutate(diam_mm = radius_mm * 2) |> #distinct(organism)
  filter(organism !='bacteria') |> 
  select(-c(organism,radius_mm,ydiam,pattern,iz,dpi,ino_date,rec_date)) |> 
  filter(pathogen=='Cf') |> 
  group_by(sample,isolate,name,pathogen) |> 
  summarise(mean_diam = mean(diam_mm),
            sd_diam=sd(diam_mm), 
            N = n(),
            se_diam = sd(diam_mm)/sqrt(n()), .groups='drop') |> 
  arrange(name) |> 
  write.xlsx('output/cf_prescreened_yeast_sup_tab1.xlsx')



# Detached Fruit Assay Analysis -------------------------------------------

# lesion diameter for Botrytis and Colletotrichum -------------------------

ld_bot_col <- read_xlsx("data/all_data.xlsx", sheet=7) %>% 
  select(treatment, av_ld, pathogen, dpi) |>
  filter(treatment %in% c(357,506,707,957,461,'SDI','Cf','Bc','switch'))

glimpse(ld_bot_col)

ld_bot_col$treatment <- as.factor(ld_bot_col$treatment)

ld_bot_col$dpi <- as.factor(ld_bot_col$dpi)

ld_bot_col$pathogen <- as.factor(ld_bot_col$pathogen)

ld_bot_col_aov <- aov(av_ld ~ treatment * pathogen * dpi, data = ld_bot_col)
summary(ld_bot_col_aov)

# Homogeneity of variance -------------------------------------------------
plot(ld_bot_col_aov, 1)
leveneTest(av_ld ~ treatment, data = ld_bot_col)

# Checking for normality --------------------------------------------------
plot(ld_bot_col_aov, 2)
ld_bot_col_aov_residuals <- residuals(object = ld_bot_col_aov)
shapiro.test(x = ld_bot_col_aov_residuals)

# Check for skewness and kurtosis -----------------------------------------
describe(ld_bot_col$av_ld, na.rm = T, type = 2)

# Kruskal-Wallis for non-normally distributed data ------------------------
# Overall
kruskal.test(av_ld ~ treatment, data = ld_bot_col) #signif
kruskal.test(av_ld ~ pathogen, data = ld_bot_col) #signif
kruskal.test(av_ld ~ dpi, data = ld_bot_col) #signif

# Kruskal-Wallis for Colletotrichum ---------------------------------------
ld_col <- ld_bot_col |> 
  filter(pathogen=='Colletotrichum')

kruskal.test(av_ld ~ treatment, data = ld_col) #signif
kruskal.test(av_ld ~ dpi, data = ld_col) #signif

# Kruskal-Wallis for Botrytis ---------------------------------------------
ld_bot <- ld_bot_col |> 
  filter(pathogen=='Botrytis')

kruskal.test(av_ld ~ treatment, data = ld_bot) #signif
kruskal.test(av_ld ~ dpi, data = ld_bot) #signif

# Post hoc test -----------------------------------------------------------
pairwise.wilcox.test(ld_bot_col$av_ld, ld_bot_col$treatment,
                     p.adjust.method = "BH")

# a. Plotting effect of different yeasts on Colletotrichum rot ------------
ld_col_plot <- ld_bot_col |>   
  filter(pathogen!="Botrytis") |> 
  mutate(treatment=case_when(treatment=='957'~'JB951',treatment=='707'~'JB767',
                             treatment=='506'~'JB624',treatment=='357'~'JB543',
                             treatment=='Cf'~'Cf only',treatment=='SDI'~'SDW',
                             treatment=='switch'~'Switch')) |> 
  mutate(treatment=fct_reorder(treatment,desc(av_ld))) |>  
  group_by(treatment, dpi) |> 
  ggbarplot(x = "treatment", y = "av_ld", add = "mean_se",
            fill = "dpi",  
            position = position_dodge(0.8)) +
  stat_compare_means(label.y = 19.5, size = 3) +
  stat_compare_means(aes(label = after_stat(p.signif)),
                     method = "wilcox.test", 
                     ref.group = "Cf only",
                     label.y = 20.5,
                     size = 4) +
  scale_y_continuous(limits = c(0,23),
                     expand = c(0,0)) +
  scale_fill_manual(name=NULL,
                    breaks=c("4", "7"),
                    labels=c("4 dpi", "7 dpi"),
                    values=c("#000","#ffffff")) +
  labs(x="",
       y="Average Lesion Diameter (mm)") +
  theme_classic() +
  theme(
    legend.text=element_text(face="plain",family="sans",size=9,color="#000"),
    legend.spacing.y=unit(1, "cm"),
    legend.position=c(.9,.8),
    axis.text=element_text(face="plain",family="sans",size=10,color="#000"),
    axis.ticks.x=element_blank(),
    axis.title =element_text(face="plain",family="sans",size=11,color="#000")
  ) +
  guides(fill = guide_legend(byrow = TRUE))

ggsave("output/ld_by_Cf.tif", ld_col_plot, dpi=300, width=5, height=5)
ggsave("output/ld_by_Cf.pdf", ld_col_plot, dpi=300, width=5, height=5)

# b. Plotting effect of different yeasts on Botrytis rot ------------------
ld_bot_plot <- ld_bot_col |>   
  filter(pathogen!="Colletotrichum") |> #distinct(treatment)
  mutate(treatment=case_when(treatment=='957'~'JB951',treatment=='707'~'JB767',
                             treatment=='461'~'JB524',treatment=='357'~'JB543',
                             treatment=='Bc'~'Bc only',treatment=='SDI'~'SDW',
                             treatment=='switch'~'Switch')) |> 
  mutate(treatment=fct_reorder(treatment,desc(av_ld))) |>  
  group_by(treatment, dpi) |> 
  ggbarplot(x = "treatment", y = "av_ld", add = "mean_se",
            fill = "dpi",  
            position = position_dodge(0.8)) +
  stat_compare_means(label.y = 23.5, size = 3) +
  stat_compare_means(aes(label = after_stat(p.signif)),
                     method = "wilcox.test", 
                     ref.group = "Bc only",
                     label.y = 25,
                     size = 4) +
  scale_y_continuous(limits = c(0,26),
                     expand = c(0,0)) +
  scale_fill_manual(name=NULL,
                    breaks=c("4", "7"),
                    labels=c("4 dpi", "7 dpi"),
                    values=c("#000","#ffffff")) +
  labs(x="",
       y="Average Lesion Diameter (mm)") +
  theme_classic() +
  guides(fill = guide_legend(byrow = TRUE)) +
  theme(
    legend.text=element_text(face="plain",family="sans",size=9,color="#000"),
    legend.spacing.y=unit(1.5, "lines"),
    #legend.spacing.x=unit(10,'cm'),
    legend.position=c(.9,.9),
    axis.text=element_text(face="plain",family="sans",size=10,color="#000"),
    axis.ticks.x=element_blank(),
    axis.title =element_text(face="plain",family="sans",size=11,color="#000")
  ) 
ggsave("output/ld_by_bc.tif", ld_bot_plot, dpi=300, width=5, height=5)
ggsave("output/ld_by_bc.pdf", ld_bot_plot, dpi=300, width=5, height=5)

#################INCIDENCE CALCULATION#########################################

bot_col_incidence <- ld_bot_col |> 
  filter(dpi==7) |> 
  mutate(pres_abs=if_else(av_ld < 1, 0, 1)) %>%
  group_by(treatment, pathogen) %>% 
  summarise(N=n(),
            incidence=100*sum(pres_abs)/length(pres_abs), .groups = "drop")


glimpse(bot_col_incidence)
bot_col_incidence_aov <- aov(incidence ~ treatment*pathogen, 
                             data = bot_col_incidence)
summary(bot_col_incidence_aov)

# Homogeneity of variance -------------------------------------------------
plot(bot_col_incidence_aov, 1)
leveneTest(incidence ~ treatment, data = bot_col_incidence)

# Checking for normality --------------------------------------------------
plot(bot_col_incidence_aov, 2)

bot_col_incidence_aov_residuals <- residuals(object = bot_col_incidence_aov)
shapiro.test(x = bot_col_incidence_aov_residuals)

# Kruskal-Wallis for non-normally distributed data ------------------------
kruskal.test(incidence ~ treatment, data = bot_col_incidence) #signif
kruskal.test(incidence ~ pathogen, data = bot_col_incidence) #signif

# Post hoc test -----------------------------------------------------------
pairwise.wilcox.test(bot_col_incidence$incidence, 
                     bot_col_incidence$treatment,
                     p.adjust.method = "BH")

# c. Plotting Incidence after 7days for Colletotrichum --------------------
col_inciden_plot <- bot_col_incidence |> 
  select(-N) |>  
  filter(pathogen!="Botrytis") |> 
  mutate(treatment=case_when(treatment=='957'~'JB951',treatment=='707'~'JB767',
                             treatment=='506'~'JB624',treatment=='357'~'JB543',
                             treatment=='Cf'~'Cf alone',treatment=='SDI'~'SDW',
                             treatment=='switch'~'Switch')) |>
  mutate(treatment=fct_reorder(treatment,desc(incidence))) |>  
  ggbarplot(x = "treatment", y = "incidence", add = "mean_se",
            fill = "treatment",  
            position = position_dodge(0.8)) +
  stat_compare_means(label.y = 105, size = 4) +
  scale_y_continuous(limits = c(0,110),
                     expand = c(0,0)) +
  scale_fill_manual(name=NULL, 
                    values=c('#f0f0f0','#f0f0f0','#f0f0f0','#f0f0f0',
                             '#f0f0f0','#f0f0f0','#f0f0f0','#f0f0f0')) +
  labs(x="",
       y="Incidence (%)") +
  theme_classic() +
  theme(
    legend.position="none",
    axis.text=element_text(face="plain",family="sans",size=10,color="black"),
    axis.ticks.x=element_blank(),
    axis.title =element_text(face="plain",family="sans",size=11,color="black")
  )
ggsave("output/incidence_Cf_rot.tif",col_inciden_plot,dpi=300,width=5,height=5)
ggsave("output/incidence_Cf_rot.pdf",col_inciden_plot,dpi=600,width=5,height=5)

# d. Plotting Incidence after 7days for Botrytis --------------------------
bot_inciden_plot <- bot_col_incidence |> 
  select(-N) |>  
  filter(pathogen!="Colletotrichum") |> 
  mutate(treatment=case_when(treatment=='957'~'JB951',treatment=='707'~'JB767',
                             treatment=='461'~'JB524',treatment=='357'~'JB543',
                             treatment=='Bc'~'Bc alone',treatment=='SDI'~'SDW',
                             treatment=='switch'~'Switch')) |> 
  mutate(treatment=fct_reorder(treatment,desc(incidence))) %>% 
  ggbarplot(., x = "treatment", y = "incidence", add = "mean_se",
            fill = "treatment",  
            position = position_dodge(0.8)) +
  stat_compare_means(label.y = 105, size = 4) +
  scale_y_continuous(limits = c(0,110),
                     expand = c(0,0)) +
  scale_fill_manual(name=NULL, 
                    values=c('#f0f0f0','#f0f0f0','#f0f0f0','#f0f0f0',
                             '#f0f0f0','#f0f0f0','#f0f0f0','#f0f0f0')) +
  labs(x="",
       y="Incidence (%)") +
  theme_classic() +
  theme(
    legend.position="none",
    axis.text=element_text(face="plain",family="sans",size=10,color="black"),
    axis.ticks.x=element_blank(),
    axis.title =element_text(face="plain",family="sans",size=11,color="black")
  )
ggsave("output/incidence_bc_rot.tif",bot_inciden_plot,dpi=300,width=5,height=5)
ggsave("output/incidence_bc_rot.pdf",bot_inciden_plot,dpi=600,width=5,height=5)

################# SEVERITY CALCULATIONA#########################################

bot_col_severity <- read_xlsx(
  "data/all_data.xlsx", sheet=8) |> #distinct(treatment)
  filter(dpi==7) |> 
  group_by(treatment,dpi,pathogen) |> 
  mutate(num_ratings=rating_scale*total_rating,
         max_rating=max(rating_scale),
         dis_sev=sum(num_ratings)/(sum(total_rating)*max_rating)*100) %>% 
  ungroup()

glimpse(bot_col_severity)

bot_col_severity$treatment <- as.factor(bot_col_severity$treatment)
bot_col_severity$pathogen <- as.factor(bot_col_severity$pathogen)
bot_col_severity_aov <- aov(dis_sev ~ treatment*pathogen, 
                            data = bot_col_severity)
summary(bot_col_severity_aov)

# Homogeneity of variance -------------------------------------------------
plot(bot_col_severity_aov, 1)

# Checking for normality --------------------------------------------------
plot(bot_col_severity_aov, 2)

bot_col_severity_aov_residuals <- residuals(object = bot_col_severity_aov)
shapiro.test(x = bot_col_severity_aov_residuals)

# Check for skewness and kurtosis -----------------------------------------
describe(bot_col_severity$dis_sev, na.rm = T, type = 2)

# Kruskal-Wallis for non-normally distributed data ------------------------
kruskal.test(dis_sev ~ treatment, data = bot_col_severity) #signif
kruskal.test(dis_sev ~ pathogen, data = bot_col_severity) #signif

# Post hoc test -----------------------------------------------------------
pairwise.wilcox.test(bot_col_severity$dis_sev, 
                     bot_col_severity$treatment,
                     p.adjust.method = "BH")

my_colors <- c("#0D0D0D","#696969","#A9A9A9","#C9C9C9",
               "#E5E5E5","#F0F0F0","#FFFFFF")

# e. Severity after 7 days for Botrytis -----------------------------------
bot_severity_plot <- bot_col_severity |>  
  select(treatment, pathogen, dis_sev) |> 
  filter(pathogen=="Botrytis") |> #distinct(treatment)
  mutate(treatment=case_when(treatment=='957'~'JB951',treatment=='707'~'JB767',
                             treatment=='461'~'JB524',treatment=='357'~'JB543',
                             treatment=='Bc'~'Bc alone',treatment=='SDI'~'SDW',
                             treatment=='Switch'~'Switch')) |>
  mutate(treatment=fct_reorder(treatment,desc(dis_sev))) |>  
  ggbarplot(x = "treatment", y = "dis_sev", add = "mean_se",
            fill = "treatment",  
            position = position_dodge(0.8)) +
  stat_compare_means(label.y = 103, size = 4) +
  stat_compare_means(aes(label = after_stat(p.signif)),
                     method = "wilcox.test", 
                     ref.group = "Bc alone",
                     label.y = 108,
                     size = 4) +
  scale_y_continuous(breaks=c(0,20,40,60,80,100),
                     limits = c(0,112),
                     expand = c(0,0)) +
  scale_fill_manual(name=NULL, values=my_colors) +
  labs(x="",
       y="Disease Severity Index (%)") +
  theme_classic() +
  theme(
    legend.position="none",
    axis.text=element_text(face="plain",family="sans",size=10,color="#000"),
    axis.ticks.x=element_blank(),
    axis.title =element_text(face="plain",family="sans",size=11,color="#000")
  )
ggsave("output/bot_severity_plot.tif",bot_severity_plot,dpi=300,width=5,height=5)
ggsave("output/bot_severity_plot.pdf",bot_severity_plot,dpi=300,width=5,height=5)

# f. Severity after 7 days for Colletotrichum -----------------------------
col_severity_plot <- bot_col_severity |>  
  select(treatment, pathogen, dis_sev) |> 
  filter(pathogen=="Colletotrichum") |> #distinct(treatment)
  mutate(treatment=case_when(treatment=='957'~'JB951',treatment=='707'~'JB767',
                             treatment=='506'~'JB624',treatment=='357'~'JB543',
                             treatment=='Cf'~'Cf alone',treatment=='SDI'~'SDW',
                             treatment=='Switch'~'Switch')) |> 
  mutate(treatment=fct_reorder(treatment,desc(dis_sev))) |>  
  ggbarplot(x = "treatment", y = "dis_sev", add = "mean_se",
            fill = "treatment",  
            position = position_dodge(0.8)) +
  stat_compare_means(label.y = 90, size = 4) +
  stat_compare_means(aes(label = after_stat(p.signif)),
                     method = "wilcox.test", 
                     ref.group = "Cf alone",
                     label.y = 100,
                     size = 4) +
  scale_y_continuous(breaks=c(0,20,40,60,80,100),
                     limits = c(0,112),
                     expand = c(0,0)) +
  scale_fill_manual(name=NULL,  values=my_colors) +
  labs(x="",
       y="Disease Severity Index (%)") +
  theme_classic() +
  theme(
    legend.position="none",
    axis.text=element_text(face="plain",family="sans",size=10,color="#000"),
    axis.ticks.x=element_blank(),
    axis.title =element_text(face="plain",family="sans",size=11,color="#000")
  )
ggsave("output/col_severity_plot.tif",col_severity_plot,dpi=300,width=5,height=5)
ggsave("output/col_severity_plot.pdf",col_severity_plot,dpi=300,width=5,height=5)


# Grape Juice Fermentation Data Analysis ----------------------------------

# fermentation plot1
ferm_data_plot <- read_excel('data/data.xlsx', sheet=9) |>
  mutate('1'=day1-day1,'2'=day1-day2,'3'=day1-day3,'4'=day1-day4,'5'=day1-day5,
         '6'=day1-day6,'7'=day1-day7,'8'=day1-day8,'9'=day1-day9,'10'=day1-day10,
         '11'=day1-day11,'12'=day1-day12,'13'=day1-day13,'14'=day1-day14,
         '15'=day1-day15,'16'=day1-day16,'17'=day1-day17,'18'=day1-day18,
         '19'=day1-day19,'20'=day1-day20) |>
  select(1,22:41) |> 
  pivot_longer(-strain, names_to='days', values_to='weight') %>% 
  ggplot(., aes(x=factor(days,level=c("1","2","3","4","5","6","7","8","9","10",
                                      "11","12","13","14","15","16","17","18",
                                      "19","20")), 
                y=weight, group=strain, color=strain)) +
  geom_line() +
  geom_point(size=1) +
  scale_color_manual(name='Strain',
                     breaks=c('CTRL','EC-1118','JB524','JB540','JB624','JB767'),
                     labels=c('CTRL','EC-1118','JB524','JB540','JB624','JB767'),
                     values=c("#B452CD","#0000ff","#ff0000","#000000","#00ced1",
                              '#228b22',"#6495ed","#8b4513","#bdb76b","#a9a9a9",
                              "#e9967a")) +
  labs(x="Days after incubation",
       y="Weight loss (g)") +
  theme_classic() +
  theme(
    legend.background=element_blank(),
    legend.text=element_text(size=10, color="#000", family='sans'),
    legend.title=element_text(size=12, color="#000", family='sans'),
    legend.position=c(.9,.65),
    axis.text=element_text(size=11, color="#000", family='sans'),
    axis.title=element_text(size=12, color="#000", family='sans')
  )
ggsave("output/ferm_plot.tif",ferm_data_plot,width=5,height=5,dpi=300)
ggsave("output/ferm_plot.pdf",ferm_data_plot,width=5,height=5,dpi=600)
